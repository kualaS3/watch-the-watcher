"""Simulate bioreactor-style process readings and inject known faults.

ALL NUMBERS HERE ARE ILLUSTRATIVE. They are made-up values chosen for this exercise,
not real process limits from any employer.

Every reading is:  centre of range + autocorrelated noise (+ optional slow cycle),
                   then a fault is added, then the value is rounded to sensor resolution.

Faults (all sized as a fraction of the channel's half-range, so "0.3" means 30% of the way
from the centre to a limit):
  spike   a short (1-3 sample) jump, e.g. from a test or a process step
  drift   a slow creep that ramps towards (and past) the range limit and then holds
  offset  a constant calibration error: reads too high or too low from then on
  stuck   the sensor repeats its last reported value from then on

Faults other than spikes persist until the end of the run, like a real fault nobody has fixed.

Routine manual checks are simulated too: every 4 hours someone reviews the whole window of
readings since the previous check (that is what the real routine check looks like) and flags
the fault if any reading in the window is outside the standard range, or if the trace has been
flat for at least an hour (a stuck sensor is only obvious once enough frozen data is on the
graph). The flag time is the check itself: nothing is "seen" until someone looks. That is the
baseline the detectors are compared with.
"""
from dataclasses import dataclass
from typing import Optional

import numpy as np

SAMPLE_SECONDS = 10                # one reading every 10 s
N_SAMPLES = 4321                   # 12 h, with the final sample landing exactly on the 12 h mark
WARMUP = 360                       # first hour is fault-free. It is used to learn the baseline
FAULT_START_MAX = 1380             # faults begin between 1.25 h and 3.8 h (before the 4 h check)
DRIFT_RAMP = 1440                  # drift takes 4 h to reach its full size
SPIKE_GRACE = 6                    # a spike counts as caught if alarmed within 1 min of ending
CHECK_EVERY = 1440                 # manual check every 4 h (samples)
FLAT_MIN_SAMPLES = 360             # 1 h of flat trace before a human calls a sensor stuck

FAULT_KINDS = ("spike", "drift", "offset", "stuck")
# Default fault sizes, as a fraction of the half-range (signed at random per run)
DEFAULT_MAGNITUDE = {"spike": 0.5, "drift": 1.2, "offset": 0.3, "stuck": 0.0}


@dataclass(frozen=True)
class Channel:
    name: str
    unit: str
    low: float             # standard range, lower limit
    high: float            # standard range, upper limit
    wobble: float          # largest swing of a healthy reading around the centre (~3 sigma)
    phi: float             # AR(1) coefficient: 0 = white noise, ~1 = slow wander
    resolution: float      # sensor quantisation step
    cycle_amp: float = 0.0     # slow cycle amplitude, in sigmas
    cycle_period: int = 120    # samples

    @property
    def mean(self) -> float:           # readings are centred in the range
        return (self.low + self.high) / 2.0

    @property
    def half_range(self) -> float:
        return (self.high - self.low) / 2.0

    @property
    def sigma(self) -> float:          # stationary std of the noise
        return self.wobble / 3.0


CHANNELS = {
    "pH": Channel("pH", "pH", 6.8, 7.2, wobble=0.055, phi=0.90, resolution=0.001),
    "temperature": Channel("temperature", "degC", 36.5, 37.5, wobble=0.12, phi=0.80,
                           resolution=0.01, cycle_amp=0.6, cycle_period=120),
    # weight shifts slowly (high phi), with a wobble of 20 kg
    "weight": Channel("weight", "kg", 2300.0, 2600.0, wobble=20.0, phi=0.97, resolution=1.0),
    # DO wobble was not specified: assumed +/- 2 % saturation
    "DO": Channel("DO", "% sat", 30.0, 50.0, wobble=2.0, phi=0.70, resolution=0.1),
}


@dataclass(frozen=True)
class Fault:
    kind: str
    start: int
    end: int           # end of the detection window (n for persistent faults)
    magnitude: float   # signed, fraction of half-range (0 for stuck)
    length: int        # samples the disturbance lasts


@dataclass
class Run:
    channel: Channel
    x: np.ndarray            # what the operator sees
    fault: Optional[Fault]
    warmup: int
    seed: int

    @property
    def minutes(self) -> np.ndarray:
        return np.arange(len(self.x)) * SAMPLE_SECONDS / 60.0


def check_indices(n: int):
    """Sample indices of the routine manual checks (4 h, 8 h, 12 h)."""
    return np.arange(CHECK_EVERY, n, CHECK_EVERY)


def out_of_range(run: Run) -> np.ndarray:
    return (run.x < run.channel.low) | (run.x > run.channel.high)


def breach_index(run: Run) -> int:
    """First sample at/after the fault start where the reading leaves the range, or -1."""
    if run.fault is None:
        return -1
    bad = out_of_range(run)
    bad[:run.fault.start] = False
    return int(np.argmax(bad)) if bad.any() else -1


def _run_starts(x: np.ndarray) -> np.ndarray:
    """Index where the run of identical values ending at each sample begins."""
    change = np.empty(len(x), dtype=bool)
    change[0] = True
    change[1:] = x[1:] != x[:-1]
    idx = np.arange(len(x))
    return np.maximum.accumulate(np.where(change, idx, 0))


def manual_check_index(run: Run) -> int:
    """First routine check that flags the fault, or -1.

    The check reviews the whole window of readings since the previous check, the way the
    real routine check does. It flags the fault at a check if either:
      * any reading in the window (at/after the fault start) is outside the range, or
      * the trace has been flat for at least FLAT_MIN_SAMPLES inside the window, counting
        only flat that belongs to the fault. A human needs about an hour of frozen data
        on the graph before calling a sensor stuck.
    The returned index is the check itself: the fault is only seen when someone looks,
    so a fault at hour 5 is "detected" by the human at the 8 h check.
    """
    if run.fault is None:
        return -1
    bad = out_of_range(run)
    bad[:run.fault.start] = False
    starts = _run_starts(run.x)
    prev = 0
    for c in check_indices(len(run.x)):
        c = int(c)
        if c >= run.fault.start:
            if bad[prev:c + 1].any():
                return c
            ends = np.arange(prev, c + 1)
            # A flat run only counts for the part of it that lies inside this window
            # and at/after the fault start. (A run can begin a sample or two before
            # the fault start: the frozen value repeats the last healthy reading, and
            # rounded healthy readings sometimes repeat too.)
            visible = ends - np.maximum.reduce(
                [starts[prev:c + 1], np.full(c + 1 - prev, prev),
                 np.full(c + 1 - prev, run.fault.start)]) + 1
            if (visible >= FLAT_MIN_SAMPLES).any():
                return c
        prev = c
    return -1


def simulate_run(channel: Channel,
                 seed: int,
                 fault_kind: Optional[str] = None,
                 magnitude: Optional[float] = None,
                 n: int = N_SAMPLES,
                 warmup: int = WARMUP,
                 spike_prob: float = 0.0) -> Run:
    """Simulate one run. fault_kind=None gives a clean run.

    spike_prob adds random background spikes to every run. The default is 0 because spikes
    are treated as faults to detect, not as noise.
    """
    if fault_kind is not None and fault_kind not in FAULT_KINDS:
        raise ValueError(f"unknown fault kind: {fault_kind}")

    rng = np.random.default_rng(seed)
    sigma, phi = channel.sigma, channel.phi

    # --- All random draws happen in a fixed order, whatever the fault kind ---
    eps = rng.normal(0.0, sigma * np.sqrt(1.0 - phi ** 2), n)
    first = rng.normal(0.0, sigma)
    phase = rng.uniform(0, 2 * np.pi)
    spike_mask = rng.random(n) < spike_prob
    spike_size = rng.uniform(4.0, 6.0, n) * rng.choice([-1.0, 1.0], n) * sigma
    start = int(rng.integers(warmup + 90, FAULT_START_MAX))
    sign = float(rng.choice([-1.0, 1.0]))
    spike_len = int(rng.integers(1, 4))

    # --- AR(1) noise: slow-ish wander like a real process ---
    noise = np.empty(n)
    noise[0] = first
    for i in range(1, n):
        noise[i] = phi * noise[i - 1] + eps[i]

    t = np.arange(n)
    cycle = channel.cycle_amp * sigma * np.sin(2 * np.pi * t / channel.cycle_period + phase)
    raw = channel.mean + noise + cycle + np.where(spike_mask, spike_size, 0.0)

    fault = None
    if fault_kind in ("spike", "drift", "offset"):
        mag = DEFAULT_MAGNITUDE[fault_kind] if magnitude is None else magnitude
        mag *= sign
        size = mag * channel.half_range
        offset = np.zeros(n)
        if fault_kind == "spike":
            offset[start:start + spike_len] = size
            fault = Fault("spike", start, start + spike_len + SPIKE_GRACE, mag, spike_len)
        elif fault_kind == "drift":
            ramp = np.clip((t[start:] - start) / DRIFT_RAMP, 0.0, 1.0)
            offset[start:] = size * ramp
            fault = Fault("drift", start, n, mag, n - start)
        else:
            offset[start:] = size
            fault = Fault("offset", start, n, mag, n - start)
        raw = raw + offset

    # --- Sensor quantisation ---
    x = np.round(raw / channel.resolution) * channel.resolution

    # --- Stuck sensor: keeps repeating its last reported value ---
    if fault_kind == "stuck":
        x[start:] = x[start - 1]
        fault = Fault("stuck", start, n, 0.0, n - start)

    return Run(channel=channel, x=x, fault=fault, warmup=warmup, seed=seed)


if __name__ == "__main__":
    for name, ch in CHANNELS.items():
        run = simulate_run(ch, seed=1, fault_kind="drift")
        print(f"{name:12s} range {ch.low}-{ch.high}  sigma {ch.sigma:.4f}  "
              f"first readings {np.round(run.x[:4], 3)}  fault {run.fault.kind} "
              f"breach@{breach_index(run)} manual@{manual_check_index(run)}")
