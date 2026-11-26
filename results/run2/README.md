# Run 2 – result of the pre-registered TRL 4 retrospective validation

| | |
|---|---|
| Run | 2 (run 1 invalid, see `../run1_invalid_et0/`) |
| Commit | `ccb84fd` (fieldcard 0.2.3; protocol v1, thresholds unchanged since tag `prereg-v1`) |
| Run finished | 2026-11-26 11:58 UTC |
| Field-seasons | 360 (30 maize + 30 winter wheat × 6 seasons); 300 used in A–C |
| Seasons | drought 2018, 2019, 2022 · normal 2017, 2021 · excluded 2020 |

## Pre-registered criteria

| Criterion | Threshold | Result | |
|---|---|---|---|
| A – share flagged, drought minus normal | ≥ 30 pp | **67.5 pp** (68.3 % vs 0.8 %) | pass |
| B – share flagged in normal seasons | ≤ 20 % | **0.8 %** | pass |
| C – median lead, first water warning → greenness loss | ≥ 7 days | **0 days** | **fail** |

C population: 135 drought field-seasons with both events; 51 negative leads,
0 left-censored, 2 with a warning only, 19 with greenness loss only.

By season and crop (share flagged): 2018 maize 93.3 %, wheat 90.0 % · 2019
maize 93.3 %, wheat 0.0 % · 2022 maize 93.3 %, wheat 40.0 % · 2017 maize and
wheat 0.0 % · 2021 maize 0.0 %, wheat 3.3 %.

## Weather

FAO-56 ET0 came from the Open-Meteo default reanalysis on all 2,191 days:
station 00445 has no wind record for 2017–2022. Station rain and temperature
were used; 12 rain days were filled from the same reanalysis. See
`weather_provenance.json`.

The notebook copy used for this run still had the 0.2.2 cells, so the run's own
summary lacks the weather-source block and its printed Section 2.4 sentence says
"DWD station weather". The model code was 0.2.3 (commit above); the sources were
read afterwards from the metadata of the run's weather file.

## Exploratory (not pre-registered)

Lead times of the 135 field-seasons in C: 51 negative, 39 zero, 45 positive
(all 45 ≥ 7 days); quartiles −7, 0, +7 days. In 106 of 165 field-seasons with a
water-deficit alert, the first such alert was a week already showing greenness
loss, i.e. the water label mostly arrives with the visible symptom rather than
before it at weekly resolution.

## Files

- `trl4_summary.json`, `trl4_field_seasons.csv`: outputs of run 2, unchanged.
- `trl4_weekly.csv`: the weekly cards of run 2 (one row per field and week), converted unchanged from the run's `trl4_weekly.parquet`.
- `weather_provenance.json`: ET0 and rain sources, read from the run's weather file.
- `run2_notebook_with_outputs.ipynb`: the executed notebook (Colab user name and
  ID removed from the cell metadata; timestamps kept).
