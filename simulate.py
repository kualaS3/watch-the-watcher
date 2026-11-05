"""Simulate process readings (pH, temperature, weight) and inject known faults.

Every reading is:  baseline + autocorrelated noise + (optional) slow cycle
                   + rare spikes, then quantised to the sensor resolution.

Faults are injected AFTER the clean signal is built, so the same seed gives a
clean run and a faulted run that are identical up to the fault start. That
pairing makes comparisons fair.

Fault sizes are given in units of the channel's noise sigma, so "a 3-sigma step"
means the same thing for pH, temperature and weight.
"""
from dataclasses import dataclass
from typing import Optional

import numpy as np

SAMPLE_SECONDS = 10          # one reading every 10 s
N_SAMPLES = 1440             # 4 hours
WARMUP = 360                 # first hour is fault-free; used to learn the baseline
FAULT_DURATION = 360         # a fault lasts 1 hour (the detection window)
DRIFT_RAMP = 240             # drift reaches full size after 40 min

FAULT_KINDS = ("drift", "step", "stuck")
DEFAULT_MAGNITUDE = {"drift": 4.0, "step": 3.0, "stuck": 0.0}  # in sigmas


@dataclass(frozen=True)
class Channel:
    name: str
    unit: str
    mean: float
    sigma: float           # stationary std of the noise
    phi: float             # AR(1) coefficient: 0 = white noise, ~1 = slow wander
    resolution: float      # sensor quantisation step
    cycle_amp: float = 0.0     # slow cycle amplitude, in sigmas
    cycle_period: int = 120    # samples


CHANNELS = {
    "pH": Channel("pH", "pH", 7.00, 0.02, 0.90, 0.01),
    "temperature": Channel("temperature", "degC", 37.0, 0.10, 0.80, 0.01,
                           cycle_amp=0.6, cycle_period=120),
    "weight": Channel("weight", "g", 500.0, 0.50, 0.20, 0.1),
}


@dataclass(frozen=True)
class Fault:
    kind: str
    start: int
    end: int
    magnitude: float   # in sigmas (signed); unused for "stuck"


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


def simulate_run(channel: Channel,
                 seed: int,
                 fault_kind: Optional[str] = None,
                 magnitude: Optional[float] = None,
                 n: int = N_SAMPLES,
                 warmup: int = WARMUP,
                 spike_prob: float = 0.0005) -> Run:
    """Simulate one run. fault_kind=None gives a clean run."""
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
    start = int(rng.integers(warmup + 90, n - FAULT_DURATION - 1))
    sign = float(rng.choice([-1.0, 1.0]))

    # --- AR(1) noise: slow-ish wander like a real process ---
    noise = np.empty(n)
    noise[0] = first
    for i in range(1, n):
        noise[i] = phi * noise[i - 1] + eps[i]

    t = np.arange(n)
    cycle = channel.cycle_amp * sigma * np.sin(2 * np.pi * t / channel.cycle_period + phase)
    raw = channel.mean + noise + cycle + np.where(spike_mask, spike_size, 0.0)

    fault = None
    end = start + FAULT_DURATION
    if fault_kind in ("drift", "step"):
        mag = DEFAULT_MAGNITUDE[fault_kind] if magnitude is None else magnitude
        mag *= sign
        offset = np.zeros(n)
        if fault_kind == "step":
            offset[start:end] = mag * sigma
        else:
            ramp = np.clip((t[start:end] - start) / DRIFT_RAMP, 0.0, 1.0)
            offset[start:end] = mag * sigma * ramp
        raw = raw + offset
        fault = Fault(fault_kind, start, end, mag)

    # --- Sensor quantisation ---
    x = np.round(raw / channel.resolution) * channel.resolution

    # --- Stuck sensor: keeps repeating its last reported value ---
    if fault_kind == "stuck":
        x[start:end] = x[start - 1]
        fault = Fault("stuck", start, end, 0.0)

    return Run(channel=channel, x=x, fault=fault, warmup=warmup, seed=seed)


if __name__ == "__main__":
    # Quick look: print a few readings from each channel
    for name, ch in CHANNELS.items():
        run = simulate_run(ch, seed=1, fault_kind="step")
        print(f"{name:12s} first readings: {np.round(run.x[:5], 3)}   fault: {run.fault}")
