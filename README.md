# FIELDCARD

Explainable GeoAI for weekly crop-stress attribution and irrigation timing.

FIELDCARD issues one card per field and week. Each card states the crop status, the most likely cause of stress (water deficit, heat, waterlogging, disease-favourable weather or "unexplained") with a confidence value, the evidence behind it, a recommended action and the cheapest next observation when the evidence is weak. The method is deliberately transparent: a daily FAO-56 root-zone water balance driven by weather and an NDVI-based crop coefficient, robust anomalies against a same-crop reference, and rule-based cause attribution. A root-zone drop counts as a water signal only when canopy moisture (NDMI) also declines, so station weather alone never flags a field.

Status: TRL 4 core (`fieldcard.core`) with an automated test suite, and a pre-registered retrospective validation on open data around DLG IPZ Bernburg, Germany.

## Install

```bash
pip install "fieldcard[data,test] @ git+https://github.com/seaycon/fieldcard.git"
python -m pytest -q
```

`fieldcard.core` needs only NumPy and pandas. The `data` extra adds the readers for open data (`fieldcard.sources`).

## Use

```python
import fieldcard as fc
weekly = fc.run_pipeline(obs, wx, fc.Params(taw_mm=250), baseline_years=[2017, 2021],
                         reference_mode="leave_one_out")
print(fc.cards_json(weekly, week="2018-06-04"))
```

`obs`: one row per field and cloud-free observation (`parcel_id, date, ndvi, ndmi, ndre, crop`).
`wx`: one row per day (`date, precip, et0, tmax, tmean, rh`).

## TRL 4 retrospective validation (pre-registered)

The protocol is fixed in [`configs/preregistered_trl4.yaml`](configs/preregistered_trl4.yaml) before any satellite data are read, and run once with [`notebooks/trl4_retrospective_validation.ipynb`](notebooks/trl4_retrospective_validation.ipynb) (Google Colab, CPU):

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/seaycon/fieldcard/blob/main/notebooks/trl4_retrospective_validation.ipynb)

- **Area and seasons:** winter-wheat and maize fields ≥ 5 ha around Bernburg (Saxony-Anhalt), seasons 2017–2022.
- **Season classes:** from station rain only (April–August rain against the 1991–2020 mean): below 80 % drought, from 90 % normal, in between excluded.
- **Reference:** normal seasons only, never including the season (or, with one normal season, the field) being evaluated.
- **Criteria:** (A) share of fields flagged with water deficit in drought seasons minus normal seasons ≥ 30 percentage points; (B) share flagged in normal seasons ≤ 20 %; (C) median lead from the first water-deficit warning to visible greenness loss ≥ 7 days. Every count behind A–C is reported, including leads that are left-censored (warning in the earliest week an alert can occur; kept as lower bounds).

**Results:** run 2 (fieldcard 0.2.3): A 67.5 pp (pass), B 0.8 % (pass), C median lead 0 days (fail). Run 1 was invalid (no ET0); both are kept in [`results/`](results/).

The run writes `trl4_summary.json`, which records the commit it used. Results are reported as they come out.

## Data and attribution

- Weather: Deutscher Wetterdienst (DWD), Climate Data Center, station 00445 Bernburg/Saale (Nord). Gap filling: ERA5-Land via Open-Meteo (contains modified Copernicus Climate Change Service information).
- Fields and crop types: Thünen Institute, national crop-type maps of Germany v201 (CC BY 4.0), doi:10.5281/zenodo.10619783 and doi:10.5281/zenodo.10621629.
- Satellite: Copernicus Sentinel-2 L2A data, accessed through the Microsoft Planetary Computer.
- Soil: usable field capacity of the loess black earth at DLG IPZ (250 mm), NOSTRADAMUS D2.1.

## Licence

Apache License 2.0. Copyright Seaycon Technology Inc.
