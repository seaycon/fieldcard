# Run 1 – invalid (no reference evapotranspiration)

| | |
|---|---|
| Run | 1 of the pre-registered TRL 4 retrospective validation |
| Commit | `7f8ab95` (tag `prereg-v1`, fieldcard 0.2.2) |
| Tag published | 2026-11-26 10:52 UTC |
| Run finished | 2026-11-26 11:25 UTC |
| Status | **Invalid – implementation defect.** Kept here unchanged as evidence. |

## What went wrong

The protocol computes FAO-56 reference evapotranspiration (ET0) from the DWD
station and fills missing days from ERA5-Land. In this run ET0 was missing on
all 2,191 days, so the root-zone water balance had no evaporative demand:
Ks stayed at 1.0 in every week, and the water-deficit cause could never be
assigned. Criteria A, B and C were therefore not tested.

1. Station 00445 has no wind record (FM) for 2017–2022, and sunshine duration
   (SDK) is missing on 40 % of the days, so no station ET0 could be computed.
2. The ERA5-Land request returned rain but empty ET0. The code fell back to
   another source only on an HTTP error, not on empty values, and passed the
   gaps on.
3. The water balance treated the missing ET0 as no demand instead of stopping.

The notebook output shows it: "days filled from ERA5-Land -> ET0: 2191" and an
April–August ET0 sum of 0 mm in every season.

## Outcome as printed (not a test result)

A: 0.0 pp (fail) · B: 0.0 % flagged (pass) · C: no field-season with both a
water warning and greenness loss (fail). 154 drought field-seasons showed
greenness loss; it was labelled "unexplained", as the water balance was inactive.

## Fix

fieldcard 0.2.3 (see `CHANGELOG.md` and `post_registration_fixes` in
`configs/preregistered_trl4.yaml`). Thresholds, rules, fields, seeds and
satellite processing are unchanged. Run 2 uses the same cached fields and
Sentinel-2 observations as run 1.

## Files

- `trl4_summary.json`, `trl4_field_seasons.csv`: outputs of run 1, unchanged.
- `run1_notebook_with_outputs.ipynb`: the executed notebook (Colab user name
  and ID removed from the cell metadata; timestamps kept).
