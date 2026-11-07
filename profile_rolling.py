"""Measure, don't guess: time the rolling z-score written four ways.

  loop_naive    plain Python, recomputes every window from scratch  O(n * window)
  loop_running  plain Python, running sums (smarter algorithm)      O(n)
  pandas        s.rolling(...).mean()/.std()
  numpy         sliding_window_view + mean/std along an axis

Also times the EWMA as a loop vs pandas' ewm.

Each result is checked against the naive loop, so a fast-but-wrong version cannot win.

Usage:  python profile_rolling.py [--window 60] [--outdir results]
"""
import argparse
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from detectors import (EPS, ewma_score, ewma_score_loop, rolling_zscore,
                       rolling_zscore_loop, rolling_zscore_loop_running,
                       rolling_zscore_numpy)

SIZES = [1_000, 10_000, 100_000, 1_000_000]
NAIVE_MAX = 100_000        # the naive loop at 1M samples would take minutes


def best_of(fn, repeat):
    best = float("inf")
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", type=int, default=60)
    ap.add_argument("--outdir", default="results")
    args = ap.parse_args()
    w = args.window
    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(0)

    impls = {
        "loop_naive": rolling_zscore_loop,
        "loop_running": rolling_zscore_loop_running,
        "pandas": rolling_zscore,
        "numpy": rolling_zscore_numpy,
    }

    rows = []
    for n in SIZES:
        x = rng.normal(7.0, 0.02, n)
        ref = rolling_zscore_loop(x, w) if n <= NAIVE_MAX else rolling_zscore(x, w)
        for name, fn in impls.items():
            if name == "loop_naive" and n > NAIVE_MAX:
                continue
            out = fn(x, w)
            ok = bool(np.allclose(out, ref, equal_nan=True, atol=1e-6))
            repeat = 3 if n <= 10_000 else 1 if name.startswith("loop") else 3
            sec = best_of(lambda: fn(x, w), repeat)
            rows.append({"task": "rolling_zscore", "n": n, "impl": name,
                         "seconds": sec, "matches_reference": ok})
            print(f"rolling_zscore n={n:>9,} {name:13s} {sec*1000:10.2f} ms   matches={ok}")

    # EWMA: recursive filter, loop vs pandas
    for n in SIZES:
        x = rng.normal(7.0, 0.02, n)
        ref = ewma_score_loop(x, 7.0, 0.02)
        for name, fn in (("loop", ewma_score_loop), ("pandas_ewm", ewma_score)):
            out = fn(x, 7.0, 0.02)
            ok = bool(np.allclose(out, ref, atol=1e-6))
            sec = best_of(lambda: fn(x, 7.0, 0.02), 3 if n <= 100_000 else 1)
            rows.append({"task": "ewma", "n": n, "impl": name, "seconds": sec,
                         "matches_reference": ok})
            print(f"ewma           n={n:>9,} {name:13s} {sec*1000:10.2f} ms   matches={ok}")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(args.outdir, "profile_timings.csv"), index=False)

    # Speed-up table for the rolling z-score
    z = df[df.task == "rolling_zscore"].pivot(index="n", columns="impl", values="seconds")
    z["pandas_speedup_vs_naive"] = z["loop_naive"] / z["pandas"]
    z["numpy_speedup_vs_naive"] = z["loop_naive"] / z["numpy"]
    print("\nSeconds per implementation (rolling z-score, window =", w, ")")
    print(z.round(4).to_string())

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for name in impls:
        g = df[(df.task == "rolling_zscore") & (df.impl == name)]
        ax.plot(g["n"], g["seconds"], marker="o", label=name)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("number of samples")
    ax.set_ylabel("seconds")
    ax.set_title(f"Rolling z-score timing (window={w})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "profile_timings.png"), dpi=130)
    print(f"\nWrote profile_timings.csv and profile_timings.png to {args.outdir}/")


if __name__ == "__main__":
    main()
