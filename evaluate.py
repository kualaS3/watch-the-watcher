"""Measure the detectors honestly.

For every detector and threshold we report
  * false alarms per hour on CLEAN data (alarm events, not samples)
  * detection rate: fraction of injected faults with an alarm that starts inside the fault window
  * detection delay: minutes from fault start to the first alarm (median over detected faults)
  * lead over the routine 4-hourly manual check: how many minutes earlier the detector alarmed
  * for drifts: how many minutes of warning before the reading actually left the standard range

Detectors are compared at the same false-alarm budget. Comparing at different false-alarm rates
would be meaningless: any detector can be made "fast" by letting it cry wolf.

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
from simulate import (CHANNELS, DRIFT_RAMP, FAULT_KINDS, SAMPLE_SECONDS, WARMUP, breach_index,
                      check_indices, manual_check_index, simulate_run)

MIN_PER_SAMPLE = SAMPLE_SECONDS / 60.0
PERSIST_WINDOW = 720     # offset / stuck must be caught within 2 h of starting


def detection_window_end(run, breach):
    """Where the 'did the detector catch it?' window closes.

    spike  : shortly after the spike ends
    drift  : when the reading actually leaves the range (after that, ordinary limit checks
             would catch it, so a detector only earns credit for EARLY warning)
    offset / stuck : 2 h after the fault starts
    Bounding the window matters: a long window lets random false alarms look like detections.
    """
    f, n = run.fault, len(run.x)
    if f.kind == "spike":
        return f.end
    if f.kind == "drift":
        return breach if breach >= 0 else min(n, f.start + DRIFT_RAMP)
    return min(n, f.start + PERSIST_WINDOW)


# ---------------------------------------------------------------- data + scores
def build_scores(n_runs, seed0):
    """Simulate everything and compute detector scores once.

    Returns
      scores[detector][(group, channel)] -> 2-D array (runs, samples)
      meta[(kind, channel)] -> (starts, ends, breach_idx, manual_idx) arrays (-1 = never)
      runs_by_key[(group, channel)] -> list of Run
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
                breach = np.array([breach_index(r) for r in runs])
                meta[key] = (np.array([r.fault.start for r in runs]),
                             np.array([detection_window_end(r, b) for r, b in zip(runs, breach)]),
                             breach,
                             np.array([manual_check_index(r) for r in runs]))
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


def _median(v):
    return float(np.median(v)) if len(v) else np.nan


def _first_rising_in_window(a, starts, ends):
    """Boolean (detected) and index of the first alarm that STARTS inside each run's window."""
    idx = np.arange(a.shape[1])
    inwin = (idx[None, :] >= starts[:, None]) & (idx[None, :] < ends[:, None])
    # An alarm already ringing before the window does not count: only rising edges do.
    prev = np.concatenate([np.zeros((a.shape[0], 1), dtype=bool), a[:, :-1]], axis=1)
    hit = (a & ~prev) & inwin
    return hit.any(axis=1), hit.argmax(axis=1)


def detection_stats(score_2d, meta, thr, clean_2d):
    """All per-fault-kind metrics for one group of faulted runs.

    clean_2d are clean-run scores for the same channel(s); they give the 'chance' rate:
    how often the detector alarms inside the same windows when nothing is wrong.
    """
    starts, ends, breach, manual = meta
    detected, first = _first_rising_in_window(score_2d > thr, starts, ends)
    chance, _ = _first_rising_in_window(clean_2d > thr, starts, ends)

    both = detected & (manual >= 0)                  # detector AND routine check flagged it
    det_breach = detected & (breach >= 0)
    return {
        "rate": detected.mean(),
        "chance": chance.mean(),
        "delay_min": _median((first[detected] - starts[detected]) * MIN_PER_SAMPLE),
        "lead_check_min": _median((manual[both] - first[both]) * MIN_PER_SAMPLE),
        "unseen_by_checks": (detected & (manual < 0)).mean(),   # caught, but checks never flagged
        "lead_breach_min": _median((breach[det_breach] - first[det_breach]) * MIN_PER_SAMPLE),
    }


def pool_meta(meta, kind, channels):
    return tuple(np.concatenate([meta[(kind, c)][i] for c in channels]) for i in range(4))


def evaluate_threshold(scores, meta, dname, thr, channels):
    """Pooled metrics over the given channels."""
    clean = np.vstack([scores[dname][("clean", c)] for c in channels])
    out = {"fa_per_hr": false_alarms_per_hour(clean, thr)}
    for kind in FAULT_KINDS:
        sc = np.vstack([scores[dname][(kind, c)] for c in channels])
        for k, v in detection_stats(sc, pool_meta(meta, kind, channels), thr, clean).items():
            out[f"{kind}_{k}"] = v
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


def manual_baseline(meta, runs_by_key):
    """What routine 4-hourly checks alone achieve (independent of any detector)."""
    rows = []
    for ch in CHANNELS:
        for kind in FAULT_KINDS:
            starts, _, breach, manual = meta[(kind, ch)]
            seen = manual >= 0
            rows.append({
                "channel": ch, "fault": kind,
                "left_range_pct": (breach >= 0).mean() * 100,
                "flagged_by_checks_pct": seen.mean() * 100,
                "check_delay_min": _median((manual[seen] - starts[seen]) * MIN_PER_SAMPLE),
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- reporting
def _fmt(v, digits=0, suffix=""):
    return "n/a" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{v:.{digits}f}{suffix}"


def markdown_summary(cal, base, target_fa, n_runs):
    lines = [f"# Results (false-alarm budget {target_fa}/hour per detector)", "",
             f"{n_runs} runs per channel per fault type. Each run is 12 h sampled every "
             f"{SAMPLE_SECONDS} s (1 h clean warm-up, fault starts between 1.25 h and 3.8 h, "
             "faults persist except spikes). All process ranges are illustrative.", "",
             "## Baseline: routine manual checks every 4 h (flag only if outside the range)", "",
             "| channel | fault | left the range (%) | flagged by a check (%) | check delay (min) |",
             "|---|---|---|---|---|"]
    for _, r in base.iterrows():
        lines.append(f"| {r['channel']} | {r['fault']} | {r['left_range_pct']:.0f} | "
                     f"{r['flagged_by_checks_pct']:.0f} | {_fmt(r['check_delay_min'])} |")

    def table(title, note, cell):
        out = ["", f"## {title}", "", note, "",
               "| channel | detector | thr | FA/hr | " + " | ".join(FAULT_KINDS) + " |",
               "|" + "---|" * (4 + len(FAULT_KINDS))]
        for _, r in cal.iterrows():
            flag = "" if r["meets_budget"] else " (!)"
            out.append(f"| {r['channel']} | {r['detector']} | {r['threshold']:.1f}{flag} | "
                       f"{r['fa_per_hr']:.2f} | " + " | ".join(cell(r, k) for k in FAULT_KINDS) + " |")
        return out

    lines += table("Detection: % of faults caught (% by chance) / median delay (min)",
                   "'By chance' = how often the same detector alarms inside the same time window "
                   "on CLEAN data, so a catch rate near the chance rate means no real skill. "
                   "Windows: spike ~1 min after it ends; drift until the reading leaves the range; "
                   "offset and stuck 2 h. Delay is in minutes; 'n/a' = never caught.",
                   lambda r, k: f"{r[f'{k}_rate']*100:.0f}% ({r[f'{k}_chance']*100:.0f}%) / "
                                f"{_fmt(r[f'{k}_delay_min'])}")
    lines += table("Drift: warning ahead of the range limit and of the routine check",
                   "Median minutes before the reading left the range / before the 4-hourly check "
                   "flagged it. Only drift runs, only where the detector fired before the breach.",
                   lambda r, k: (f"{_fmt(r['drift_lead_breach_min'])} / {_fmt(r['drift_lead_check_min'])}"
                                 if k == "drift" else ""))
    lines += ["", "(!) = no threshold in the grid met the false-alarm budget; the least sensitive "
              "one is shown."]
    return "\n".join(lines)


def plot_tradeoff(pooled, path):
    fig, axes = plt.subplots(2, len(FAULT_KINDS), figsize=(16, 7), sharex=True)
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
    """One example run per fault, with the standard range, manual checks and calibrated alarms."""
    combos = [("pH", "drift"), ("temperature", "offset"), ("weight", "stuck"), ("DO", "spike")]
    fig, axes = plt.subplots(len(combos), 1, figsize=(12, 11))
    colors = {"zscore": "tab:red", "ewma": "tab:green", "cusum": "tab:purple", "flatline": "tab:orange"}
    for ax, (ch_name, kind) in zip(axes, combos):
        run = runs_by_key[(kind, ch_name)][0]
        ch = run.channel
        mu, sg = estimate_baseline(run.x, run.warmup)
        h = run.minutes / 60.0
        span = ch.high - ch.low
        ax.plot(h, run.x, lw=0.8, color="0.25")
        ax.axhline(ch.low, color="crimson", ls="--", lw=1)
        ax.axhline(ch.high, color="crimson", ls="--", lw=1, label="standard range")
        for c in check_indices(len(run.x)):
            ax.axvline(c * MIN_PER_SAMPLE / 60.0, color="steelblue", ls=":", lw=1.2)
        ax.axvline(run.fault.start * MIN_PER_SAMPLE / 60.0, color="goldenrod", lw=2, label="fault starts")
        for k, det in enumerate(DETECTORS):
            row = cal[(cal.channel == ch_name) & (cal.detector == det)].iloc[0]
            alarm = np.where(DETECTORS[det][0](run.x, mu, sg) > row["threshold"])[0]
            alarm = alarm[alarm >= run.warmup]
            if len(alarm):
                ax.scatter(h[alarm], np.full(len(alarm), ch.low - (0.06 + 0.04 * k) * span),
                           s=5, color=colors[det], label=det)
        ax.set_ylim(ch.low - 0.25 * span, ch.high + 0.08 * span)
        ax.set_title(f"{ch_name}: {kind}   (blue dotted = 4-hourly manual checks)")
        ax.set_ylabel(ch.unit)
        ax.legend(loc="upper left", ncol=6, fontsize=8)
    axes[-1].set_xlabel("hours")
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
    pd.concat([per_channel, pooled], ignore_index=True).to_csv(
        os.path.join(args.outdir, "tradeoff.csv"), index=False)

    cal = calibrate(per_channel, args.target_fa)
    cal.to_csv(os.path.join(args.outdir, "calibrated.csv"), index=False)
    base = manual_baseline(meta, runs_by_key)
    base.to_csv(os.path.join(args.outdir, "manual_baseline.csv"), index=False)

    md = markdown_summary(cal, base, args.target_fa, args.runs)
    with open(os.path.join(args.outdir, "summary.md"), "w") as f:
        f.write(md + "\n")
    print(md)

    plot_tradeoff(pooled, os.path.join(args.outdir, "tradeoff.png"))
    plot_examples(cal, runs_by_key, os.path.join(args.outdir, "example_runs.png"))
    print(f"\nWrote tradeoff.csv, calibrated.csv, manual_baseline.csv, summary.md, "
          f"tradeoff.png, example_runs.png to {args.outdir}/")


if __name__ == "__main__":
    main()
