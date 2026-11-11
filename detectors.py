"""
Simple statistical fault detectors.

Each detector turns a series of readings into a score (bigger = more suspicious).
An alarm is raised when score > threshold. Keeping score and threshold separate lets
evaluate.py sweep the threshold and plot the fast-vs-false-alarm trade-off.

Detectors
---------
rolling_zscore : how far is this reading from the mean of the previous `window` readings?
ewma_score     : exponentially weighted moving average of the deviation from baseline
cusum_score    : cumulative sum of deviations (sensitive to small persistent shifts)
flatline_score : is the recent variation suspiciously close to zero? (stuck sensor)

The rolling z-score exists in four forms (plain loop, running-sum loop, pandas, NumPy)
and the tests check they all agree.
"""
import numpy as np
import pandas as pd
from numpy.lib.stride_tricks import sliding_window_view

EPS = 1e-9   # floor on a standard deviation, to avoid dividing by zero


# --------------------------------------------------------------------------
# Baseline from the fault-free warm-up period
# --------------------------------------------------------------------------
def estimate_baseline(x, warmup):
    """Robust mean and sigma from the first `warmup` samples (median / MAD)."""
    seg = np.asarray(x[:warmup], dtype=float)
    mu = float(np.median(seg))
    sigma = float(1.4826 * np.median(np.abs(seg - mu)))
    return mu, max(sigma, EPS)


# --------------------------------------------------------------------------
# Rolling z-score: four implementations of the same calculation
# z[t] = (x[t] - mean(x[t-w : t])) / std(x[t-w : t], ddof=1)   (window excludes x[t])
# --------------------------------------------------------------------------
def rolling_zscore_loop(x, window):
    """Plain Python loop. Recomputes every window from scratch: O(n * window)."""
    xs = [float(v) for v in x]
    n = len(xs)
    z = np.full(n, np.nan)
    for t in range(window, n):
        win = xs[t - window:t]
        mean = sum(win) / window
        var = sum((v - mean) ** 2 for v in win) / (window - 1)
        sd = max(var ** 0.5, EPS)
        z[t] = (xs[t] - mean) / sd
    return z


def rolling_zscore(x, window):
    """Vectorised pandas version (the one the detector uses)."""
    s = pd.Series(np.asarray(x, dtype=float))
    mean = s.rolling(window).mean().shift(1)
    sd = s.rolling(window).std().shift(1)
    return ((s - mean) / np.maximum(sd, EPS)).to_numpy()


def rolling_zscore_numpy(x, window):
    """Vectorised NumPy version using sliding windows (views, no copying)."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    z = np.full(n, np.nan)
    if n <= window:
        return z
    win = sliding_window_view(x, window)[:-1]    # window ending at t-1, for t = window..n-1
    mean = win.mean(axis=1)
    sd = np.maximum(win.std(axis=1, ddof=1), EPS)
    z[window:] = (x[window:] - mean) / sd
    return z


def rolling_zscore_loop_running(x, window):
    """Plain loop, but O(n): keep running sums instead of re-summing each window."""
    xs = [float(v) for v in x]
    shift = xs[0]
    xs = [v - shift for v in xs]          # z-score is shift-invariant. It reduces cancellation error
    n = len(xs)
    z = np.full(n, np.nan)
    s = sum(xs[:window])
    s2 = sum(v * v for v in xs[:window])
    for t in range(window, n):
        mean = s / window
        var = max((s2 - window * mean * mean) / (window - 1), 0.0)
        z[t] = (xs[t] - mean) / max(var ** 0.5, EPS)
        old, new = xs[t - window], xs[t]
        s += new - old
        s2 += new * new - old * old
    return z


def zscore_score(x, mu0, sigma0, window=60):
    """Score = |rolling z|. (mu0, sigma0 unused: this detector learns its own baseline.)"""
    return np.abs(rolling_zscore(x, window))


# --------------------------------------------------------------------------
# EWMA
# --------------------------------------------------------------------------
def ewma_score_loop(x, mu0, sigma0, lam=0.1):
    """EWMA control chart, plain loop. Score is in units of the EWMA's own std."""
    y = 0.0
    out = np.empty(len(x))
    for i, v in enumerate(x):
        y = lam * (float(v) - mu0) / sigma0 + (1 - lam) * y
        out[i] = abs(y)
    return out / np.sqrt(lam / (2 - lam))


def ewma_score(x, mu0, sigma0, lam=0.1):
    """Same thing using pandas' built-in exponentially weighted mean."""
    d = (np.asarray(x, dtype=float) - mu0) / sigma0
    y = pd.Series(np.concatenate([[0.0], d])).ewm(alpha=lam, adjust=False).mean().to_numpy()[1:]
    return np.abs(y) / np.sqrt(lam / (2 - lam))


# --------------------------------------------------------------------------
# CUSUM (two-sided). Inherently sequential, so it stays a loop.
# --------------------------------------------------------------------------
def cusum_score(x, mu0, sigma0, k=0.5):
    """Score = max(S+, S-) in sigma units. k is the slack (half the shift to detect)."""
    z = (np.asarray(x, dtype=float) - mu0) / sigma0
    sp = sn = 0.0
    out = np.empty(len(z))
    for i, v in enumerate(z):
        sp = max(0.0, sp + v - k)
        sn = max(0.0, sn - v - k)
        out[i] = max(sp, sn)
    return out


# --------------------------------------------------------------------------
# Flatline (stuck sensor)
# --------------------------------------------------------------------------
def flatline_score(x, mu0, sigma0, window=12):
    """Score = sigma0 / (rolling std + small floor). Huge when readings stop varying."""
    sd = pd.Series(np.asarray(x, dtype=float)).rolling(window).std().to_numpy()
    return sigma0 / (sd + 1e-3 * sigma0)


# --------------------------------------------------------------------------
# Registry used by evaluate.py: name -> (score function, threshold grid)
# --------------------------------------------------------------------------
DETECTORS = {
    "zscore":   (zscore_score,   np.arange(3.0, 10.01, 0.5)),
    "ewma":     (ewma_score,     np.arange(2.0, 30.01, 0.5)),
    "cusum":    (cusum_score,    np.arange(2.0, 120.01, 2.0)),
    "flatline": (flatline_score, np.array([2, 3, 5, 10, 20, 50, 100, 500], dtype=float)),
}
