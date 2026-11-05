"""Run with:  python -m pytest -q"""
import numpy as np

from detectors import (cusum_score, estimate_baseline, ewma_score, ewma_score_loop,
                       flatline_score, rolling_zscore, rolling_zscore_loop,
                       rolling_zscore_numpy, zscore_score)
from profile_rolling import rolling_zscore_loop_running
from simulate import CHANNELS, FAULT_DURATION, WARMUP, simulate_run


# ---------------- simulator ----------------
def test_same_seed_is_reproducible():
    a = simulate_run(CHANNELS["pH"], seed=5, fault_kind="step")
    b = simulate_run(CHANNELS["pH"], seed=5, fault_kind="step")
    assert np.array_equal(a.x, b.x)


def test_clean_and_faulted_runs_match_before_fault():
    clean = simulate_run(CHANNELS["weight"], seed=7)
    bad = simulate_run(CHANNELS["weight"], seed=7, fault_kind="step")
    s = bad.fault.start
    assert np.array_equal(clean.x[:s], bad.x[:s])
    assert not np.array_equal(clean.x[s:], bad.x[s:])


def test_fault_labels_are_sane():
    for kind in ("drift", "step", "stuck"):
        run = simulate_run(CHANNELS["temperature"], seed=3, fault_kind=kind)
        assert run.fault.kind == kind
        assert run.fault.start >= WARMUP
        assert run.fault.end - run.fault.start == FAULT_DURATION
        assert run.fault.end <= len(run.x)


def test_stuck_sensor_repeats_one_value():
    run = simulate_run(CHANNELS["pH"], seed=11, fault_kind="stuck")
    seg = run.x[run.fault.start:run.fault.end]
    assert np.all(seg == seg[0])


def test_warmup_is_fault_free():
    for seed in range(20):
        run = simulate_run(CHANNELS["pH"], seed=seed, fault_kind="drift")
        assert run.fault.start > run.warmup


# ---------------- rolling z-score: all versions agree ----------------
def test_zscore_versions_agree():
    x = np.random.default_rng(0).normal(7, 0.02, 800)
    ref = rolling_zscore_loop(x, 60)
    for fn in (rolling_zscore, rolling_zscore_numpy, rolling_zscore_loop_running):
        assert np.allclose(fn(x, 60), ref, equal_nan=True, atol=1e-6), fn.__name__


def test_zscore_excludes_current_sample():
    x = np.zeros(100) + np.random.default_rng(1).normal(0, 1, 100)
    x[80] = 50.0                                   # huge spike must show up as a huge z
    assert rolling_zscore(x, 30)[80] > 20


def test_zscore_first_window_is_nan():
    z = rolling_zscore(np.arange(100.0), 20)
    assert np.all(np.isnan(z[:20])) and not np.isnan(z[20])


# ---------------- EWMA / CUSUM ----------------
def test_ewma_vectorised_matches_loop():
    x = np.random.default_rng(2).normal(10, 1, 500)
    assert np.allclose(ewma_score(x, 10, 1, 0.1), ewma_score_loop(x, 10, 1, 0.1))


def test_ewma_and_cusum_catch_a_step_quickly():
    rng = np.random.default_rng(3)
    x = rng.normal(0, 1, 600)
    x[400:] += 3.0
    for fn, thr in ((ewma_score, 4.0), (cusum_score, 8.0)):
        sc = fn(x, 0.0, 1.0)
        assert not (sc[:400] > thr).any(), fn.__name__          # quiet on clean data
        first = np.argmax(sc[400:] > thr)
        assert (sc[400:] > thr).any() and first < 30, fn.__name__


def test_cusum_is_two_sided():
    x = np.random.default_rng(4).normal(0, 1, 600)
    x[300:] -= 2.0
    assert (cusum_score(x, 0.0, 1.0)[300:] > 10).any()


def test_zscore_misses_a_slow_drift_that_cusum_catches():
    # The rolling baseline follows the drift, so z never gets big. CUSUM remembers.
    rng = np.random.default_rng(5)
    x = rng.normal(0, 1, 1000)
    x[400:] += np.minimum(np.arange(600) / 240.0, 1.0) * 4.0
    assert not (zscore_score(x, 0, 1)[400:] > 5).any()
    assert (cusum_score(x, 0.0, 1.0)[400:] > 15).any()


# ---------------- flatline ----------------
def test_flatline_flags_constant_readings_only():
    x = np.random.default_rng(6).normal(5, 0.1, 300)
    x[200:260] = x[199]
    sc = flatline_score(x, 5.0, 0.1, window=12)
    assert (sc[215:260] > 100).all()
    assert not (sc[:190] > 20).any()


def test_baseline_estimate_is_close():
    x = np.random.default_rng(8).normal(37, 0.1, 400)
    mu, sg = estimate_baseline(x, 300)
    assert abs(mu - 37) < 0.05 and abs(sg - 0.1) < 0.03
