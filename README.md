# Watch the watcher

Hello! This is my first project that I'm endeavouring on and the inspiration stems from my current job as a biotechnologist working in the biopharma industry. A part of my job is to monitor the bioreactors and tanks for deviations that may negatively impact our production pace and the product quality. So, biotechs like myself do routine checks about every 4 hours and log any deviations in the process readings (pH, temperature, weight, dissolved oxygen etc) that are outside our standard/expected range of values. 

And well, I thought to myself, logging the deviations is good but what if there's a way to either automate or do some form of early detection? 

From my experience and understanding, there are 2 main mechanisms of deviations that I've seen at work. The first would be sudden spikes that occur once, usually caused by a test or a process that suddenly increases one of the reading metrics (i mean, hopefully not due to a sudden catastrophic failure). And the second would be gradual increases or decreases that slowly creeps towards the min/max edge of the standard/expected range limits.

Then there's of course deviations caused my instrument failure or error in calibration. Not frequent, but they do occur. 

So in this project, my aim is to explore the use of statistical detectors and explore early detection ideas. Lets see how it goes? :)

---

## Disclaimer? Values are only meant for simulation

Every range, noise level and fault size below is an illustrative, made-up value list that I came up with for this project. Nothing here is derived or taken from any real life process or actual data. They are simply ranges that I think should be typical for these variables.

| Channel | Standard range | Healthy wobble (max swing) | Resolution | Behaviour |
|---|---|---|---|---|
| pH | 6.8 to 7.2 | about 0.055 | 3 dp | wanders moderately slowly |
| Temperature (C) | 36.5 to 37.5 | about 0.12 | 2 dp | thermal inertia plus a slow cycle |
| Weight (kg) | 2300 to 2600 | about ~20 | 1 kg | shifts very slowly |
| Dissolved oxygen (% sat) | 30 to 50 | about 2 (assumed) | 0.1 | fairly quick wobble |

Each reading is centred in its range. "Wobble" is treated as roughly 3 sigma of the noise.

## What the simulation does

Every run is 12 hours of readings, one every 10 seconds. The first hour is fault-free so the detectors get a clean baseline to learn from (real life is rarely this kind, but one step at a time). After that, I inject exactly one fault per run. Fault sizes are measured as a fraction of the half-range, so "30%" means 30% of the way from the centre of the range to one of its limits:

| Fault | What it models | Size | Behaviour |
|---|---|---|---|
| **spike** | a test or process step that briefly bumps a reading | 50% of half-range | lasts 1 to 3 samples, then recovers |
| **drift** | slow creep towards the limit | ramps to 120% of half-range over 4 h | passes the limit, then holds |
| **offset** | calibration error | 30% of half-range | constant once it starts, usually stays inside the range |
| **stuck** | instrument failure: reading freezes | n/a | repeats its last value until the end |

I simulate the routine checks too, because that's what I'm comparing against: at 4, 8 and 12 hours someone reviews the readings since the last check and flags the fault if any of them is outside the standard range, or if the trace has been flat for at least an hour (a frozen sensor is only obvious once enough flat data is on the graph). That is the baseline the detectors have to beat. There's not much point building a detector that does worse than a human reviewing the log three times a shift.

## Detectors

So for different fault modes and mechanisms, there won't be a one-size-fits-all ideal detector. This is a brief summary of the detectors that I've tried to implement:

- **rolling z-score**: how far is this reading from the mean of the previous 10 minutes?
- **Exponentially Weighted Moving Average (EWMA)**: weighted average that exponentially decreases weight to older data.
- **Cumulative Sum (CUSUM)**: running total (accumulation) of deviations. This should be good at small persistent shifts.
- **flatline**: is the recent variation suspiciously close to zero? (stuck sensors)

I want to point out that although the detectors are mostly complementary, where each one specialises in a certain fault modality, I wanted to try both CUSUM and EWMA as I've never implemented such "memory-based" control charting before and I thought it might be interesting to try, especially in the case where we need to detect persistent or gradual drifts in the data. Then the obvious (to me) follow up was "which is better". And so the next step: the enemy of joy, comparison.

## Comparing the detectors and measuring them

Here's how I score them:

- **False alarms per hour** on clean data, counted as alarm events, not individual samples. (One alarm that stays on for an hour is still one alarm.)
- **Caught**: I count a fault as caught if an alarm starts inside a window that depends on the fault: for a spike, about 1 minute after it ends; for a drift, any time before the reading actually leaves the range; for offsets and stuck sensors, 2 hours.
-  The window sizes are also based on the worst case timing of the routine check. If a fault starts right after a check, nobody is looking at the process for almost 4 hours, and a drift still needs to actually leave the range before the next check can flag it. So I think the windows are really asking whether the detector managed to say something while the human was guaranteed to be blind.
- **Chance rate**: how often the same detector alarms inside the same window on clean data. This matters a lot, and I found it out the hard way: at 0.5 false alarms per hour, a 2-hour window contains a random alarm roughly half the time, so a catch rate near the chance rate means the detector has no real skill. An earlier version of mine used long windows and made the z-score look good at drifts when it was basically guessing.
- **Delay**: median minutes from when the fault starts to the first alarm.
- **Warning time** (drift only): how much warning I get before the reading leaves the range, and before the routine check would have flagged it.

And to keep the comparison honest, I tune every detector to the same false-alarm budget (default 0.5/hour, tuned per channel). A detector that alarms constantly will "catch" everything and be completely useless, so the budget is kind of the whole game.

## Run it (Google Colab or any Python 3.9+)

```
pip install -r requirements.txt
python -m pytest -q          # 25 tests
python evaluate.py           # ~10 s, writes results/
```

The options I actually use: `python evaluate.py --runs 200 --target-fa 1.0 --seed 42`, if I want more runs or a looser alarm budget. For clicking through: the notebook `watch_the_watcher.ipynb` runs everything in Colab.

## Results snapshot (seed 0, 100 runs per cell, 0.5 false alarms/hour budget)

These are my results from one run of the evaluation. The full tables are in `results/summary.md`, and the plots are `results/example_runs.png` and `results/tradeoff.png`. If you run it yourself and get slightly different numbers, trust yours. This is a snapshot, not a gospel.

This table compares the human routine check versus the performance of the best detector. In practice, the tested best detector would be used so this comparison makes the most sense. I think it sums up the project's results pretty nicely.

| Fault | Routine check | Best detector |
|---|---|---|
| Drift | 100% caught, but ~5.5h late | EWMA: 100%, ~140–175 min *before* the reading leaves the range, ~300 min before the check |
| Spike | 0% — never leaves the range | Rolling z-score: 94–100%, instantly |
| Offset | 0% — sits inside the range forever | EWMA: 98–100%, in ~1–3 min |
| Stuck | 100%, but ~2h late | Flatline: 100%, in ~2 min |

So the detectors don't just catch more faults than the standard process monitoring routine checks, but they are also warn the user hours before the problem reaches criticality, instead of hours after it starts. Studying the tradeoffs shows that no single detector is a one-size-fits-all, and using different detectors for different fault modalities is the way to go.

1. **The routine checks only catch faults that leave the range, plus frozen sensors once enough flat data has piled up.** In my simulation they flagged 100% of drifts (about 5.5 h after the drift began, at the 8 h check) and 100% of stuck sensors, but only after about 2 hours, because a human needs at least an hour of flat trace before calling a sensor stuck. Spikes and offsets still scored 0%: spikes in my simulation never actually leave the range, and a calibration offset can sit inside the range indefinitely, so the routine check will never see it. This was honestly the result that sold me on the whole project.

The ~5.5 h figure is not an estimate I picked. It is the median delay the evaluation measured for the routine check on drifts (320-330 min across the four channels, in results/manual_baseline.csv): the drift first has to ramp far enough to leave the range, then it still has to wait for the next 4-hourly check to be seen.

2. **Drift: EWMA and CUSUM give hours of warning.** EWMA caught 100% of drifts on every channel, roughly 140 to 175 minutes before the reading left the range and about 300 minutes before the routine check would have flagged it. CUSUM did similarly on pH, temperature and DO. The rolling z-score was much weaker: on pH it caught 32% of drifts, no better than its 39% chance rate, because its baseline follows the drift along with it.
3. **Spikes: only the rolling z-score.** It caught 94 to 100% with a 0% chance rate. EWMA caught 15 to 45% and CUSUM 4 to 14%, because a single spike barely moves an averaged score. Which makes sense when I think about it: the memory that helps them with drifts is exactly what makes them blind here.
4. **Offsets: all of the level-based detectors react within about 3 minutes**, though EWMA (98 to 100%) was more reliable than the z-score (64 to 100%), whose baseline adapts to the new level.
5. **Stuck sensor: only the flatline detector works** (100% within about 2 minutes). EWMA, CUSUM and z-score are at or below their chance rates. The routine check does eventually flag a stuck sensor too, but only about 2 hours later, once a full hour of flat trace is sitting in the window it reviews. So no single detector covers every fault type, and I need one per fault type. A bit annoying, but that's the answer.
6. **Slow-moving signals break textbook thresholds.** Weight shifts very slowly, so the CUSUM threshold I needed to stay under the false-alarm budget was very high. It then caught only 63% of drifts (chance rate 65%) and 60% of offsets (chance 55%): no real skill. Control-chart limits assume independent noise, and weight does not give me independent noise.

I think a practical and real world benefit that can be derived from the results is the gain in time that the process can be in the optimal range. Because with earlier and more consistent detection, and 24/7 monitoring, faults are fixed earlier and 1-2 hours on each fault would translate into substantial gains to the bioreactor's production output. Especially when accounting for the high sensitivity of these process factors in relation to their productive efficiency. This is what gives me the motivation to explore this topic in these directions.

## Limitations

- This project uses simulated data from ranges, noise levels and fault shapes that I came up with. This is only a representation, and is not real data.
- The faults are not attest to real faults, theres only one fault per run and overlapping faults were not tested. Spikes are also modelled as just 1 to 3 samples.
- detector baselines are derived from the clean first hour, which sometimes real processes may not have the luxury of giving.
- the thresholds currently are tuned and tested on the same data, which is prime territory for textbook overfitting. I've done it this way for now for the simplicity of just using one set of data. It's definitely not hard to have a separate set of seeds to tune first and then test on a fresh set.

## Files

If you're reading the code: `simulate.py` and `detectors.py` are the heart of it, and everything else is measuring or testing them.

| File | What it does |
|---|---|
| `simulate.py` | Channels, noise model, fault injection, manual-check helpers |
| `detectors.py` | Rolling z-score (4 implementations), EWMA, CUSUM, flatline |
| `evaluate.py` | Detection, delay, false alarms, chance rate, calibration, plots |
| `test_watch_the_watcher.py` | 25 tests |
| `watch_the_watcher.ipynb` | Colab notebook that runs everything |
