# Changelog

## 0.2.3 (2026-11-26), post-registration fix after run 1

Run 1 (tag `prereg-v1`, 0.2.2) is invalid and kept unchanged in
`results/run1_invalid_et0/`. Thresholds, rules, fields, seeds and satellite
processing are unchanged.

- **Weather:** ET0 was missing on every day of run 1. Station 00445 has no wind
  record for 2017–2022, so no station ET0 could be computed, and ERA5-Land
  returned empty ET0. Empty ERA5-Land values are now asked again from the
  Open-Meteo default reanalysis (0.2.2 did this only on an HTTP error). The
  weather table records how many days each source filled and stops with an
  error if ET0 or rain is still missing.
- **Water balance:** refuses missing rain, ET0 or Kcb instead of treating a gap
  as zero demand.
- **Notebook:** prints and saves the ET0 sources; the Section 2.4 sentence names
  them.

## 0.2.2 (2026-11-26), before the pre-registration tag

- **Criterion C counter:** a lead is flagged left-censored when the first water
  warning falls in the earliest week an alert can occur. Alerts need
  `persist_weeks` consecutive weeks, so when the window opens with an assessed
  week the earliest alert is the second week, not the first; 0.2.1 compared
  with the first assessed week and so almost never set the flag. A, B and the
  C median are unchanged; only the reported count changes. New column
  `earliest_possible_alert` in the season summary.

## 0.2.1 (2026-11-26), before the pre-registration commit

- **Water signal:** a drop of Ks below the usual value for the week now needs a
  canopy-moisture decline as well (NDMI z ≤ `ks_branch_ndmi_z`, −0.5). Before,
  the weather-driven branch alone could flag every field sharing one station in
  a dry spell of a normal season. New parameter in `Params` and in
  `configs/preregistered_trl4.yaml`.
- **Thünen maps:** the file name of each year is checked against the Zenodo API
  (falls back to the known name if the API cannot be reached); if the legend
  gives no class attribute, or the file lacks it, the attribute is found from
  the values, which must be legend codes. Codes compare as numbers or text.
- **Sentinel-2:** each day is read with one parallel `.compute()` instead of one
  read per band; `dask[array]` added to the `data` extra.
- **Notebook:** the Section 2.4 sentence counts drought + normal field-seasons
  and the classified seasons only, and lists the excluded seasons.
- **Criterion C:** field-seasons whose first water warning falls in their first
  assessed week are flagged as left-censored (lead is a lower bound); they stay
  in C and their number is reported.

## 0.2.0

- Leave-one-out reference, season classes, TRL 4 criteria, data sources,
  pre-registered configuration and Colab notebook.
