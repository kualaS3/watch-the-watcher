# Results (false-alarm budget 0.5/hour per detector)

100 runs per channel per fault type. Each run is 12 h sampled every 10 s (1 h clean warm-up, fault starts between 1.25 h and 3.8 h, faults persist except spikes). All process ranges are illustrative.

## Baseline: routine manual checks every 4 h (each check reviews the readings since the last one; flags a fault if any reviewed reading is outside the range, or the trace has been flat for at least 1 h)

| channel | fault | left the range (%) | flagged by a check (%) | check delay (min) |
|---|---|---|---|---|
| pH | spike | 0 | 0 | n/a |
| pH | drift | 100 | 100 | 320 |
| pH | offset | 0 | 0 | n/a |
| pH | stuck | 0 | 100 | 138 |
| temperature | spike | 0 | 0 | n/a |
| temperature | drift | 100 | 100 | 325 |
| temperature | offset | 0 | 0 | n/a |
| temperature | stuck | 0 | 100 | 122 |
| weight | spike | 0 | 0 | n/a |
| weight | drift | 100 | 100 | 320 |
| weight | offset | 0 | 0 | n/a |
| weight | stuck | 0 | 100 | 144 |
| DO | spike | 0 | 0 | n/a |
| DO | drift | 100 | 100 | 327 |
| DO | offset | 0 | 0 | n/a |
| DO | stuck | 0 | 100 | 133 |

## Detection: % of faults caught (% by chance) / median delay (min)

'By chance' = how often the same detector alarms inside the same time window on CLEAN data, so a catch rate near the chance rate means no real skill. Windows: spike ~1 min after it ends; drift until the reading leaves the range; offset and stuck 2 h. Delay is in minutes; 'n/a' = never caught.

| channel | detector | thr | FA/hr | spike | drift | offset | stuck |
|---|---|---|---|---|---|---|---|
| DO | cusum | 28.0 | 0.48 | 14% (1%) / 0 | 99% (60%) / 15 | 99% (45%) / 1 | 59% (41%) / 6 |
| DO | ewma | 6.5 | 0.44 | 45% (1%) / 0 | 100% (58%) / 16 | 100% (46%) / 0 | 17% (47%) / 2 |
| DO | flatline | 5.0 | 0.01 | 0% (0%) / n/a | 3% (2%) / 110 | 1% (1%) / 49 | 100% (2%) / 2 |
| DO | zscore | 3.5 | 0.46 | 100% (0%) / 0 | 83% (76%) / 49 | 99% (59%) / 0 | 0% (65%) / n/a |
| pH | cusum | 60.0 | 0.49 | 4% (0%) / 0 | 95% (58%) / 27 | 92% (49%) / 3 | 47% (55%) / 11 |
| pH | ewma | 9.0 | 0.43 | 17% (0%) / 0 | 100% (55%) / 28 | 100% (54%) / 1 | 3% (49%) / 2 |
| pH | flatline | 10.0 | 0.01 | 0% (0%) / n/a | 0% (0%) / n/a | 0% (0%) / n/a | 100% (0%) / 2 |
| pH | zscore | 4.0 | 0.16 | 94% (0%) / 0 | 32% (39%) / 68 | 64% (30%) / 0 | 0% (29%) / n/a |
| temperature | cusum | 40.0 | 0.45 | 4% (0%) / 0 | 97% (52%) / 19 | 100% (43%) / 2 | 48% (42%) / 9 |
| temperature | ewma | 7.5 | 0.40 | 15% (0%) / 0 | 100% (53%) / 20 | 100% (52%) / 1 | 10% (43%) / 2 |
| temperature | flatline | 5.0 | 0.19 | 0% (0%) / n/a | 23% (38%) / 74 | 13% (27%) / 49 | 100% (30%) / 2 |
| temperature | zscore | 4.0 | 0.18 | 97% (0%) / 0 | 39% (47%) / 67 | 76% (29%) / 0 | 0% (33%) / n/a |
| weight | cusum | 70.0 | 0.50 | 6% (1%) / 0 | 63% (65%) / 19 | 60% (55%) / 1 | 38% (56%) / 14 |
| weight | ewma | 12.0 | 0.45 | 33% (1%) / 0 | 100% (53%) / 20 | 98% (43%) / 1 | 3% (46%) / 1 |
| weight | flatline | 10.0 | 0.33 | 0% (0%) / n/a | 24% (36%) / 41 | 32% (29%) / 42 | 100% (25%) / 1 |
| weight | zscore | 4.0 | 0.24 | 100% (0%) / 0 | 62% (55%) / 78 | 100% (33%) / 0 | 0% (36%) / n/a |

## Drift: warning ahead of the range limit and of the routine check

Median minutes before the reading left the range / before the 4-hourly check flagged it. Only drift runs, only where the detector fired before the breach.

| channel | detector | thr | FA/hr | spike | drift | offset | stuck |
|---|---|---|---|---|---|---|---|
| DO | cusum | 28.0 | 0.48 |  | 161 / 316 |  |  |
| DO | ewma | 6.5 | 0.44 |  | 160 / 314 |  |  |
| DO | flatline | 5.0 | 0.01 |  | 64 / 236 |  |  |
| DO | zscore | 3.5 | 0.46 |  | 120 / 266 |  |  |
| pH | cusum | 60.0 | 0.49 |  | 142 / 295 |  |  |
| pH | ewma | 9.0 | 0.43 |  | 139 / 295 |  |  |
| pH | flatline | 10.0 | 0.01 |  | n/a / n/a |  |  |
| pH | zscore | 4.0 | 0.16 |  | 91 / 254 |  |  |
| temperature | cusum | 40.0 | 0.45 |  | 153 / 301 |  |  |
| temperature | ewma | 7.5 | 0.40 |  | 151 / 302 |  |  |
| temperature | flatline | 5.0 | 0.19 |  | 99 / 238 |  |  |
| temperature | zscore | 4.0 | 0.18 |  | 87 / 236 |  |  |
| weight | cusum | 70.0 | 0.50 |  | 179 / 310 |  |  |
| weight | ewma | 12.0 | 0.45 |  | 173 / 299 |  |  |
| weight | flatline | 10.0 | 0.33 |  | 150 / 251 |  |  |
| weight | zscore | 4.0 | 0.24 |  | 115 / 241 |  |  |

(!) = no threshold in the grid met the false-alarm budget; the least sensitive one is shown.
