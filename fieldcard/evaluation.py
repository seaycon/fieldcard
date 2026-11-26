"""Season classification and the pre-registered TRL 4 criteria.

Both functions are written so that every number in the evidence can be traced
back to one rule in configs/preregistered_trl4.yaml.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def classify_seasons(precip: pd.Series, seasons, months=(4, 8), normal_period=(1991, 2020),
                     drought_below=0.80, normal_at_or_above=0.90, max_missing_frac=0.05,
                     min_normal_years=25) -> pd.DataFrame:
    """Label each season drought / normal / excluded from station rainfall.

    precip : daily precipitation (mm) indexed by date.
    The ratio is the season's rainfall in `months` (inclusive) divided by the
    mean of the same months over `normal_period`. Seasons below `drought_below`
    are drought seasons, seasons at or above `normal_at_or_above` are normal,
    the rest are excluded. If no season qualifies as normal, the season with the
    highest ratio among the non-drought seasons is used as normal; if every
    season is a drought season, the wettest season is used (label
    'normal_fallback'). Seasons with too many missing days are excluded.
    """
    s = precip.copy()
    s.index = pd.to_datetime(s.index)
    s = s[(s.index.month >= months[0]) & (s.index.month <= months[1])]
    days = s.groupby(s.index.year).size()
    missing = s.isna().groupby(s.index.year).mean()
    total = s.groupby(s.index.year).sum(min_count=1)
    expected = pd.Series({y: len(pd.date_range(f"{y}-{months[0]:02d}-01",
                                               pd.Timestamp(f"{y}-{months[1]:02d}-01") + pd.offsets.MonthEnd(0)))
                          for y in total.index})
    complete = (missing <= max_missing_frac) & (days >= expected * (1 - max_missing_frac))
    ref_years = [y for y in range(normal_period[0], normal_period[1] + 1) if y in total.index and complete.get(y, False)]
    if len(ref_years) < min_normal_years:
        raise ValueError(f"only {len(ref_years)} complete years in the normal period; need {min_normal_years}")
    normal_mean = float(total.loc[ref_years].mean())

    rows = []
    for y in seasons:
        if y not in total.index or not complete.get(y, False):
            rows.append({"season": y, "rain_mm": np.nan, "ratio": np.nan, "class": "excluded_missing_data"})
            continue
        r = float(total.loc[y]) / normal_mean
        cls = "drought" if r < drought_below else ("normal" if r >= normal_at_or_above else "excluded")
        rows.append({"season": y, "rain_mm": round(float(total.loc[y]), 1), "ratio": round(r, 3), "class": cls})
    out = pd.DataFrame(rows)
    if not (out["class"] == "normal").any():
        cand = out[out["class"].isin(["excluded", "drought"])].copy()
        non_drought = cand[cand["class"] == "excluded"]
        pick_from = non_drought if len(non_drought) else cand
        if len(pick_from):
            i = pick_from["ratio"].idxmax()
            out.loc[i, "class"] = "normal_fallback"
    out.attrs["normal_mean_mm"] = round(normal_mean, 1)
    out.attrs["normal_years_used"] = len(ref_years)
    return out


def _pct(x: pd.Series) -> float:
    return float(100.0 * x.mean()) if len(x) else float("nan")


def trl4_criteria(summary: pd.DataFrame, classes: pd.DataFrame, min_weeks: int = 2,
                  a_min_pp: float = 30.0, b_max_pct: float = 20.0, c_min_days: float = 7.0) -> dict:
    """Criteria A, B and C on field-season summaries (fieldcard.core.season_summary).

    flagged  : at least `min_weeks` alert weeks attributed to water deficit.
    A        : % flagged in drought seasons minus % flagged in normal seasons (pp).
    B        : % flagged in normal seasons.
    C        : median days from the first water-deficit warning to the first
               greenness-loss week, drought seasons, field-seasons with both
               events; negative values are kept. Field-seasons whose warning
               falls in the earliest week an alert can occur are left-censored
               (their lead is a lower bound); they stay in C and are counted.
    """
    cls = classes.set_index("season")["class"]
    s = summary.copy()
    s["season_class"] = s["year"].map(cls)
    s["flagged"] = s["weeks_water_deficit"] >= min_weeks
    drought = s[s["season_class"] == "drought"]
    normal = s[s["season_class"].isin(["normal", "normal_fallback"])]

    pct_d, pct_n = _pct(drought["flagged"]), _pct(normal["flagged"])
    a = pct_d - pct_n
    with_both = drought.dropna(subset=["lead_water_days"])
    censored = with_both["lead_left_censored"].astype(bool) if "lead_left_censored" in with_both else \
        pd.Series(False, index=with_both.index)
    both = with_both
    c = float(both["lead_water_days"].median()) if len(both) else float("nan")

    def breakdown(df):
        return [{"season": int(y), "crop": cr, "fields": int(len(g)), "pct_flagged": round(_pct(g["flagged"]), 1),
                 "median_lead_days": (None if g["lead_water_days"].isna().all()
                                      else float(g["lead_water_days"].median()))}
                for (y, cr), g in df.groupby(["year", "crop"])]

    return {
        "fields": {"drought_field_seasons": int(len(drought)), "normal_field_seasons": int(len(normal)),
                   "excluded_field_seasons": int(len(s) - len(drought) - len(normal))},
        "A_drought_separation_pp": round(a, 1),
        "A_pass": bool(a >= a_min_pp),
        "pct_flagged_drought": round(pct_d, 1),
        "B_normal_flag_pct": round(pct_n, 1),
        "B_pass": bool(pct_n <= b_max_pct),
        "C_median_lead_days": None if np.isnan(c) else round(c, 1),
        "C_pass": bool(c >= c_min_days) if not np.isnan(c) else False,
        "C_counts": {
            "with_warning_and_loss": int(len(with_both)),
            "left_censored_lower_bound": int(censored.sum()),
            "warning_only": int((drought["first_water_warning"].notna() & drought["first_greenness_loss"].isna()).sum()),
            "loss_only": int((drought["first_water_warning"].isna() & drought["first_greenness_loss"].notna()).sum()),
            "negative_leads": int((both["lead_water_days"] < 0).sum()),
        },
        "thresholds": {"A_min_pp": a_min_pp, "B_max_pct": b_max_pct, "C_min_days": c_min_days,
                       "flagged_min_weeks": min_weeks},
        "by_season_and_crop": breakdown(s[s["season_class"].isin(["drought", "normal", "normal_fallback"])]),
    }
