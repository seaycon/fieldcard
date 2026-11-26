# Changelog

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
