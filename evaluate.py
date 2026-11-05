"""Measure the detectors honestly.

For every detector and threshold we report
  * false alarms per hour on CLEAN data (alarm events, not samples)
  * detection rate: fraction of injected faults with at least one alarm inside the fault window
  * detection delay: minutes from fault start to the first alarm (median over detected faults)

Then we calibrate each detector to the same false-alarm budget and compare delay at that budget.
Comparing detectors at different false-alarm rates would be meaningless: any detector can be made
"fast" by letting it cry wolf.

Usage:  python evaluate.py [--runs 100] [--target-fa 0.5] [--seed 0] [--outdir results]
"""
import argparse
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from detectors import DETECTORS, estimate_baseline
from simulate import (CHANNELS, FAULT_KINDS, N_SAMPLES, SAMPLE_SECONDS, WARMUP,
                      simulate_run)

MIN_PER_SAMPLE = SAMPLE_SECONDS / 60.0


# ---------------------------------------------------------------- data + scores
def build_scores(n_runs, seed0):
    """Simulate everything and compute detector scores once.

    Returns scores[detector][(group, channel)] -> 2-D array (runs, samples) and
    meta[(group, channel)] -> (starts, ends) arrays for faulted groups.
    group is "clean" or a fault kind.
    """
    scores = {d: {} for d in DETECTORS}
    meta = {}
    runs_by_key = {}
    counter = seed0
    for ch_name, ch in CHANNELS.items():
        for group in ("clean",) + FAULT_KINDS:
            runs = []
            for _ in range(n_runs):
                kind = None if group == "clean" else group
                runs.append(simulate_run(ch, seed=counter, fault_kind=kind))
                counter += 1
            key = (group, ch_name)
            runs_by_key[key] = runs
            if group != "clean":
                meta[key] = (np.array([r.fault.start for r in runs]),
                             np.array([r.fault.end for r in runs]))
            baselines = [estimate_baseline(r.x, r.warmup) for r in runs]
            for dname, (fn, _) in DETECTORS.items():
                scores[dname][key] = np.vstack(
                    [fn(r.x, mu, sg) for r, (mu, sg) in zip(runs, baselines)])
    return scores, meta, runs_by_key


# ---------------------------------------------------------------- metrics
def false_alarms_per_hour(score_2d, thr):
    """Alarm EVENTS (rising edges) per hour of clean monitoring after warm-up."""
    a = score_2d > thr                      # NaN > thr is False
    seg = a[:, WARMUP - 1:]
    edges = seg[:, 1:] & ~seg[:, :-1]
    hours = a.shape[0] * (a.shape[1] - WARMUP) * SAMPLE_SECONDS / 3600.0
    return edges.sum() / hours


def detection_stats(score_2d, starts, ends, thr):
    """(detection rate, median delay in minutes) for one group of faulted runs."""
    a = score_2d > thr
    idx = np.arange(a.shape[1])
    inwin = (idx[None, :] >= starts[:, None]) & (idx[None, :] < ends[:, None])
    # Credit a detection only for an alarm that STARTS inside the fault window
    # (a rising edge). An alarm already ringing before the fault does not count.
    prev = np.concatenate([np.zeros((a.shape[0], 1), dtype=bool), a[:, :-1]], axis=1)
    hit = (a & ~prev) & inwin
    detected = hit.any(axis=1)
    first = hit.argmax(axis=1)
    delays = (first[detected] - starts[detected]) * MIN_PER_SAMPLE
    return detected.mean(), (float(np.median(delays)) if detected.any() else np.nan)


def evaluate_threshold(scores, meta, dname, thr, channels):
    """Pooled metrics over the given channels."""
    clean = np.vstack([scores[dname][("clean", c)] for c in channels])
    out = {"fa_per_hr": false_alarms_per_hour(clean, thr)}
    for kind in FAULT_KINDS:
        sc = np.vstack([scores[dname][(kind, c)] for c in channels])
        st = np.concatenate([meta[(kind, c)][0] for c in channels])
        en = np.concatenate([meta[(kind, c)][1] for c in channels])
        rate, delay = detection_stats(sc, st, en, thr)
        out[f"{kind}_rate"] = rate
        out[f"{kind}_delay_min"] = delay
    return out


def sweep(scores, meta, channels, label):
    rows = []
    for dname, (_, grid) in DETECTORS.items():
        for thr in grid:
            r = evaluate_threshold(scores, meta, dname, thr, channels)
            rows.append({"channel": label, "detector": dname, "threshold": float(thr), **r})
    return pd.DataFrame(rows)


def calibrate(df, target_fa):
    """Most sensitive threshold whose false-alarm rate is <= target (per detector, per channel)."""
    picks = []
    for (ch, det), g in df.groupby(["channel", "detector"]):
        ok = g[g["fa_per_hr"] <= target_fa].sort_values("threshold")
        if len(ok):
            row = ok.iloc[0].copy()
            row["meets_budget"] = True
        else:
            row = g.sort_values("threshold").iloc[-1].copy()
            row["meets_budget"] = False
        picks.append(row)
    return pd.DataFrame(picks).reset_index(drop=True)


# ---------------------------------------------------------------- reporting
def markdown_summary(cal, target_fa, n_runs):
    lines = [f"# Detector comparison at a false-alarm budget of {target_fa} per hour",
             "",
             f"{n_runs} runs per channel per fault type; each run is 4 h sampled every "
             f"{SAMPLE_SECONDS} s (1 h warm-up, 1 h fault window). Delay = median minutes from "
             "fault start to first alarm, over detected faults only. 'n/a' = never detected.",
             ""]
    header = "| channel | detector | thr | FA/hr | " + " | ".join(
        f"{k} detected / delay (min)" for k in FAULT_KINDS) + " |"
    lines += [header, "|" + "---|" * (4 + len(FAULT_KINDS))]
    for _, r in cal.iterrows():
        cells = []
        for k in FAULT_KINDS:
            d = r[f"{k}_delay_min"]
            cells.append(f"{r[f'{k}_rate']*100:.0f}% / " + ("n/a" if np.isnan(d) else f"{d:.1f}"))
        flag = "" if r["meets_budget"] else " (!)"
        lines.append(f"| {r['channel']} | {r['detector']} | {r['threshold']:.2f}{flag} | "
                     f"{r['fa_per_hr']:.2f} | " + " | ".join(cells) + " |")
    lines += ["", "(!) = no threshold in the grid met the false-alarm budget; the least sensitive "
              "one is shown."]
    return "\n".join(lines)


def plot_tradeoff(pooled, path):
    fig, axes = plt.subplots(2, len(FAULT_KINDS), figsize=(13, 7), sharex=True)
    for j, kind in enumerate(FAULT_KINDS):
        for det, g in pooled.groupby("detector"):
            g = g.sort_values("fa_per_hr")
            axes[0, j].plot(g["fa_per_hr"], g[f"{kind}_rate"] * 100, marker="o", ms=3, label=det)
            axes[1, j].plot(g["fa_per_hr"], g[f"{kind}_delay_min"], marker="o", ms=3, label=det)
        axes[0, j].set_title(f"{kind} fault")
        axes[0, j].set_ylim(-3, 103)
        axes[1, j].set_xlabel("false alarms per hour (clean data)")
        axes[0, j].set_xscale("symlog", linthresh=0.01)
        axes[1, j].set_xscale("symlog", linthresh=0.01)
    axes[0, 0].set_ylabel("faults detected (%)")
    axes[1, 0].set_ylabel("median detection delay (min)")
    axes[0, 0].legend()
    fig.suptitle("Catching faults fast vs crying wolf (all channels pooled)")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_examples(cal, runs_by_key, path):
    """One example run per channel/fault, with alarms from each calibrated detector."""
    combos = [("pH", "drift"), ("temperature", "step"), ("weight", "stuck")]
    fig, axes = plt.subplots(len(combos), 1, figsize=(12, 9))
    colors = {"zscore": "tab:red", "ewma": "tab:green", "cusum": "tab:purple", "flatline": "tab:orange"}
    for ax, (ch_name, kind) in zip(axes, combos):
        run = runs_by_key[(kind, ch_name)][0]
        mu, sg = estimate_baseline(run.x, run.warmup)
        m = run.minutes
        ax.plot(m, run.x, lw=0.8, color="0.3")
        ax.axvspan(run.fault.start * MIN_PER_SAMPLE, run.fault.end * MIN_PER_SAMPLE,
                   color="gold", alpha=0.25, label=f"{kind} fault")
        for k, det in enumerate(DETECTORS):
            row = cal[(cal.channel == ch_name) & (cal.detector == det)].iloc[0]
            sc = DETECTORS[det][0](run.x, mu, sg)
            alarm = np.where(sc > row["threshold"])[0]
            alarm = alarm[alarm >= run.warmup]
            if len(alarm):
                ax.scatter(m[alarm], np.full(len(alarm), run.x.min() - (k + 1) * 0.35 * sg),
                           s=6, color=colors[det], label=det)
        ax.set_title(f"{ch_name}: {kind}")
        ax.set_ylabel(run.channel.unit)
        ax.legend(loc="upper left", ncol=5, fontsize=8)
    axes[-1].set_xlabel("minutes")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", type=int, default=100, help="runs per channel per group")
    ap.add_argument("--target-fa", type=float, default=0.5, help="false alarms per hour budget")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--outdir", default="results")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    t0 = time.perf_counter()
    scores, meta, runs_by_key = build_scores(args.runs, args.seed)
    print(f"Simulated and scored in {time.perf_counter() - t0:.1f}s")

    per_channel = pd.concat([sweep(scores, meta, [c], c) for c in CHANNELS], ignore_index=True)
    pooled = sweep(scores, meta, list(CHANNELS), "all")
    full = pd.concat([per_channel, pooled], ignore_index=True)
    full.to_csv(os.path.join(args.outdir, "tradeoff.csv"), index=False)

    cal = calibrate(per_channel, args.target_fa)
    cal.to_csv(os.path.join(args.outdir, "calibrated.csv"), index=False)
    md = markdown_summary(cal, args.target_fa, args.runs)
    with open(os.path.join(args.outdir, "summary.md"), "w") as f:
        f.write(md + "\n")
    print(md)

    plot_tradeoff(pooled, os.path.join(args.outdir, "tradeoff.png"))
    plot_examples(cal, runs_by_key, os.path.join(args.outdir, "example_runs.png"))
    print(f"\nWrote tradeoff.csv, calibrated.csv, summary.md, tradeoff.png, example_runs.png "
          f"to {args.outdir}/")


if __name__ == "__main__":
    main()
