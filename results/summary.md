# Detector comparison at a false-alarm budget of 0.5 per hour

100 runs per channel per fault type; each run is 4 h sampled every 10 s (1 h warm-up, 1 h fault window). Delay = median minutes from fault start to first alarm, over detected faults only. 'n/a' = never detected.

| channel | detector | thr | FA/hr | drift detected / delay (min) | step detected / delay (min) | stuck detected / delay (min) |
|---|---|---|---|---|---|---|
| pH | cusum | 84.00 | 0.44 | 98% / 18.2 | 91% / 4.3 | 35% / 15.5 |
| pH | ewma | 11.00 | 0.39 | 99% / 18.2 | 99% / 2.0 | 8% / 2.7 |
| pH | flatline | 10.00 | 0.02 | 1% / 42.7 | 5% / 42.5 | 100% / 1.7 |
| pH | zscore | 4.00 | 0.30 | 38% / 20.8 | 56% / 0.0 | 0% / n/a |
| temperature | cusum | 44.00 | 0.48 | 100% / 15.6 | 97% / 3.0 | 50% / 7.2 |
| temperature | ewma | 8.00 | 0.35 | 100% / 15.8 | 100% / 1.5 | 10% / 2.0 |
| temperature | flatline | 5.00 | 0.13 | 14% / 21.2 | 12% / 47.3 | 100% / 1.5 |
| temperature | zscore | 4.00 | 0.26 | 27% / 30.0 | 50% / 0.0 | 0% / n/a |
| weight | cusum | 12.00 | 0.25 | 100% / 9.2 | 100% / 0.7 | 51% / 3.0 |
| weight | ewma | 4.50 | 0.21 | 100% / 9.2 | 100% / 0.7 | 27% / 1.8 |
| weight | flatline | 3.00 | 0.08 | 0% / n/a | 3% / 13.3 | 100% / 1.5 |
| weight | zscore | 3.50 | 0.48 | 44% / 30.3 | 80% / 0.3 | 0% / n/a |

(!) = no threshold in the grid met the false-alarm budget; the least sensitive one is shown.
