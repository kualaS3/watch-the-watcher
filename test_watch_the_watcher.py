"""Run with:  python -m pytest -q"""
import numpy as np

from detectors import (cusum_score, estimate_baseline, ewma_score, ewma_score_loop,
                       flatline_score, rolling_zscore, rolling_zscore_loop,
                       rolling_zscore_loop_running, rolling_zscore_numpy, zscore_score)
from simulate import (CHANNELS, CHECK_EVERY, FAULT_KINDS, FLAT_MIN_SAMPLES, SPIKE_GRACE,
                      WARMUP, breach_index, check_indices, manual_check_index, simulate_run)


# ---------------- channels ----------------
def test_ranges_are_the_ones_requested_and_centred():
    c = CHANNELS
    assert (c["pH"].low, c["pH"].high) == (6.8, 7.2)
    assert (c["temperature"].low, c["temperature"].high) == (36.5, 37.5)
    assert (c["weight"].low, c["weight"].high) == (2300.0, 2600.0)
    assert (c["DO"].low, c["DO"].high) == (30.0, 50.0)
    for ch in c.values():
        assert ch.mean == (ch.low + ch.high) / 2
        assert abs(ch.sigma * 3 - ch.wobble) < 1e-12


def test_clean_readings_stay_inside_the_range_and_match_resolution():
    for name, ch in CHANNELS.items():
        for seed in range(5):
            x = simulate_run(ch, seed=seed).x
            assert x.min() >= ch.low and x.max() <= ch.high, name
            assert np.allclose(x / ch.resolution, np.round(x / ch.resolution), atol=1e-6)


def test_ph_has_three_decimals_and_temperature_two():
    ph = simulate_run(CHANNELS["pH"], seed=1).x
    assert np.allclose(ph, np.round(ph, 3)) and not np.allclose(ph, np.round(ph, 2))
    t = simulate_run(CHANNELS["temperature"], seed=1).x
    assert np.allclose(t, np.round(t, 2))


# ---------------- simulator ----------------
def test_same_seed_is_reproducible():
    a = simulate_run(CHANNELS["pH"], seed=5, fault_kind="offset")
    b = simulate_run(CHANNELS["pH"], seed=5, fault_kind="offset")
    assert np.array_equal(a.x, b.x)


def test_clean_and_faulted_runs_match_before_fault():
    for kind in FAULT_KINDS:
        clean = simulate_run(CHANNELS["weight"], seed=7)
        bad = simulate_run(CHANNELS["weight"], seed=7, fault_kind=kind)
        s = bad.fault.start
        assert np.array_equal(clean.x[:s], bad.x[:s]), kind


def test_fault_labels_are_sane():
    for kind in FAULT_KINDS:
        run = simulate_run(CHANNELS["temperature"], seed=3, fault_kind=kind)
        assert run.fault.kind == kind
        assert run.fault.start > WARMUP
        assert run.fault.end <= len(run.x)


def test_spike_is_short_and_then_recovers():
    for seed in range(10):
        clean = simulate_run(CHANNELS["DO"], seed=seed)
        run = simulate_run(CHANNELS["DO"], seed=seed, fault_kind="spike")
        f = run.fault
        assert 1 <= f.length <= 3 and f.end == f.start + f.length + SPIKE_GRACE
        assert np.array_equal(clean.x[f.start + f.length:], run.x[f.start + f.length:])
        assert not np.array_equal(clean.x[f.start:f.start + f.length], run.x[f.start:f.start + f.length])


def test_offset_is_constant_after_the_start():
    clean = simulate_run(CHANNELS["pH"], seed=2)
    run = simulate_run(CHANNELS["pH"], seed=2, fault_kind="offset")
    diff = (run.x - clean.x)[run.fault.start:]
    assert np.allclose(diff, diff[0], atol=0.0011)       # equal up to rounding
    assert abs(abs(diff[0]) - 0.3 * CHANNELS["pH"].half_range) < 0.0011


def test_drift_ramps_towards_and_past_the_limit():
    for name, ch in CHANNELS.items():
        run = simulate_run(ch, seed=4, fault_kind="drift")
        assert breach_index(run) > run.fault.start, name


def test_stuck_sensor_repeats_one_value_until_the_end():
    run = simulate_run(CHANNELS["pH"], seed=11, fault_kind="stuck")
    seg = run.x[run.fault.start:]
    assert np.all(seg == seg[0])


def test_warmup_is_fault_free():
    for seed in range(20):
        run = simulate_run(CHANNELS["pH"], seed=seed, fault_kind="drift")
        assert run.fault.start > run.warmup


# ---------------- manual checks ----------------
def test_check_times_are_every_four_hours():
    assert list(check_indices(4321)) == [CHECK_EVERY, 2 * CHECK_EVERY, 3 * CHECK_EVERY]
    assert CHECK_EVERY * 10 / 3600 == 4.0


def test_checks_see_drift_and_stuck_but_not_offset_or_spike():
    # Window-review checks: a drift leaves the range inside the window, and a stuck
    # sensor shows >= 1 h of flat trace, so both get flagged at a check. Offsets stay
    # inside the range, and the simulated spikes (50% of half-range) never leave it
    # either, so neither is ever flagged.
    for seed in range(5):
        ch = CHANNELS["pH"]
        assert manual_check_index(simulate_run(ch, seed, "drift")) > 0
        assert manual_check_index(simulate_run(ch, seed, "stuck")) > 0
        for kind in ("offset", "spike"):
            assert manual_check_index(simulate_run(ch, seed, kind)) == -1, kind
    assert manual_check_index(simulate_run(CHANNELS["pH"], 0)) == -1      # clean run: no fault


def test_window_review_flags_a_breach_between_checks():
    # A single out-of-range reading in the middle of a window is seen at that
    # window's check, even though the reading at the check itself is in range.
    run = simulate_run(CHANNELS["pH"], seed=9, fault_kind="offset")
    run.x[CHECK_EVERY - 10] = CHANNELS["pH"].high + 0.05
    assert manual_check_index(run) == CHECK_EVERY


def test_stuck_is_flagged_only_once_an_hour_of_flat_is_visible():
    for seed in range(5):
        run = simulate_run(CHANNELS["pH"], seed, "stuck")
        if run.fault.start + FLAT_MIN_SAMPLES <= CHECK_EVERY:
            assert manual_check_index(run) == CHECK_EVERY
        else:
            assert manual_check_index(run) == 2 * CHECK_EVERY


def test_manual_check_flags_a_reading_outside_the_range():
    run = simulate_run(CHANNELS["pH"], seed=9, fault_kind="offset")
    run.x[CHECK_EVERY] = CHANNELS["pH"].high + 0.05
    run.x[run.fault.start:CHECK_EVERY] = 7.0
    assert manual_check_index(run) == CHECK_EVERY


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
