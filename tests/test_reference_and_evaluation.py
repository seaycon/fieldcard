"""Tests for the leave-one-out reference and the pre-registered criteria."""
import numpy as np
import pandas as pd
import pytest

import fieldcard as fc
from fieldcard import core
from fieldcard.evaluation import classify_seasons, trl4_criteria


def _daily(rows):
    df = pd.DataFrame(rows, columns=["parcel_id", "year", "crop", "doy", "ndvi"])
    df["date"] = pd.to_datetime(df["year"].astype(str) + "-01-01") + pd.to_timedelta(df["doy"] - 1, unit="D")
    for c in ("ndmi", "ndre"):
        df[c] = df["ndvi"]
    for c in ("ndvi", "ndmi", "ndre"):
        df[f"{c}_obs"] = df[c]
    return df


def test_fold_labels():
    df = pd.DataFrame({"parcel_id": ["a", "b", "c"], "year": [2017, 2018, 2021]})
    lab = core._fold_labels(df, [2017, 2021], "leave_one_out")
    assert lab.tolist() == ["y2017", "all", "y2021"]
    lab1 = core._fold_labels(df, [2021], "leave_one_out")
    assert lab1.tolist() == ["all", "all", "pc"]
    assert core._fold_labels(df, [2017], "pooled").eq("all").all()
    with pytest.raises(ValueError):
        core._fold_labels(df, [2017], "other")


def test_leave_one_year_out_excludes_own_season():
    # two normal seasons with different values; each must be compared with the other one
    rows = [(f"{y}_{i}", y, "wheat", 150, v) for y, v in ((2017, 0.80), (2021, 0.60)) for i in range(3)]
    rows += [(f"2018_{i}", 2018, "wheat", 150, 0.50) for i in range(3)]
    d = core._add_reference(_daily(rows), [2017, 2021], fc.Params(), "leave_one_out")
    ref = d.groupby("year")["ndvi_ref"].first()
    assert ref[2017] == pytest.approx(0.60)       # reference = the other normal season
    assert ref[2021] == pytest.approx(0.80)
    assert ref[2018] == pytest.approx(0.70)       # drought season uses both normal seasons


def test_leave_one_field_out_with_single_normal_season():
    rows = [("A", 2021, "wheat", 150, 0.80), ("B", 2021, "wheat", 150, 0.60), ("C", 2021, "wheat", 150, 0.70)]
    d = core._add_reference(_daily(rows), [2021], fc.Params(), "leave_one_out").set_index("parcel_id")
    assert d.loc["A", "ndvi_ref"] == pytest.approx(0.65)   # median of B and C
    assert d.loc["B", "ndvi_ref"] == pytest.approx(0.75)


def test_pipeline_runs_with_leave_one_out():
    obs, wx, truth = fc.synthetic_dataset(years=(2017, 2018, 2021), drought_years=(2018,), n_parcels=12, seed=7)
    obs["parcel_id"] = obs["date"].dt.year.astype(str) + "_" + obs["parcel_id"]
    w = fc.run_pipeline(obs, wx, baseline_years=[2017, 2021], reference_mode="leave_one_out")
    s = fc.season_summary(w)
    assert {"lead_water_days", "first_water_warning", "first_greenness_loss"} <= set(s.columns)
    weeks = s.groupby("year")["weeks_water_deficit"].mean()
    assert weeks[2018] > weeks[2017] and weeks[2018] > weeks[2021]


def test_classify_seasons_rules_and_fallback():
    days = pd.date_range("1991-01-01", "2022-12-31", freq="D")
    p = pd.Series(1.0, index=days)
    summer = (days.month >= 4) & (days.month <= 8)
    p[summer] = 2.0
    p[(days.year == 2018) & summer] = 1.0       # ratio ~0.5 -> drought
    p[(days.year == 2019) & summer] = 1.7       # ratio ~0.85 -> excluded
    c = classify_seasons(p, [2018, 2019, 2021]).set_index("season")
    assert c.loc[2018, "class"] == "drought"
    assert c.loc[2019, "class"] == "excluded"
    assert c.loc[2021, "class"] == "normal"
    # no normal season: the wettest non-drought season becomes the fallback
    c2 = classify_seasons(p, [2018, 2019]).set_index("season")
    assert c2.loc[2019, "class"] == "normal_fallback"


def test_trl4_criteria():
    summ = pd.DataFrame({
        "parcel_id": list("abcdefgh"), "year": [2018] * 4 + [2021] * 4, "crop": ["wheat"] * 8,
        "weeks_water_deficit": [3, 2, 0, 4, 0, 2, 0, 0],
        "lead_water_days": [7, 14, np.nan, -7, np.nan, np.nan, np.nan, np.nan],
        "first_water_warning": [pd.Timestamp("2018-05-01")] * 4 + [pd.NaT] * 4,
        "first_greenness_loss": [pd.Timestamp("2018-05-08")] * 2 + [pd.NaT] + [pd.Timestamp("2018-04-24")] + [pd.NaT] * 4,
    })
    classes = pd.DataFrame({"season": [2018, 2021], "class": ["drought", "normal"]})
    r = trl4_criteria(summ, classes)
    assert r["pct_flagged_drought"] == 75.0
    assert r["B_normal_flag_pct"] == 25.0
    assert r["A_drought_separation_pp"] == 50.0
    assert r["C_median_lead_days"] == 7.0
    assert r["C_counts"]["negative_leads"] == 1
    assert r["A_pass"] and not r["B_pass"] and r["C_pass"]


def _weekly(ks_min, ndmi_z):
    n = len(ks_min)
    return pd.DataFrame({
        "parcel_id": "f1", "year": 2021, "crop": "winter_wheat",
        "week": pd.date_range("2021-05-03", periods=n, freq="7D"),
        "ndvi_z": 0.0, "ndre_z": 0.0, "ndmi_z": ndmi_z, "ndvi_ref": 0.8,
        "ks_min": ks_min, "ks_ref": 1.0, "hot_days": 0, "sat_days": 0, "excess": 0.0,
        "humid_warm_days": 0, "gap_weeks": 0, "n_obs": 1,
    })


def test_weather_alone_never_raises_a_water_alert():
    # root zone clearly below the usual level, canopy moisture normal -> no alert
    w = fc.attribute(_weekly([0.6] * 4, [0.0] * 4))
    assert (w["status"] == "normal").all()
    # same root zone with a mild canopy-moisture decline -> early warning, cause water deficit
    w2 = fc.attribute(_weekly([0.6] * 4, [-0.6] * 4))
    assert (w2["status"].iloc[1:] == "early_warning").all()
    assert (w2["cause"].iloc[1:] == "water_deficit").all()


def test_left_censored_leads_are_kept_and_counted():
    summ = pd.DataFrame({
        "parcel_id": ["a", "b", "c"], "year": [2018] * 3, "crop": ["maize"] * 3,
        "weeks_water_deficit": [3, 3, 3], "lead_water_days": [0.0, 14.0, 21.0],
        "lead_left_censored": [True, False, False],
        "first_water_warning": [pd.Timestamp("2018-06-04")] * 3,
        "first_greenness_loss": [pd.Timestamp("2018-06-04"), pd.Timestamp("2018-06-18"), pd.Timestamp("2018-06-25")],
    })
    classes = pd.DataFrame({"season": [2018], "class": ["drought"]})
    r = trl4_criteria(summ, classes)
    assert r["C_median_lead_days"] == 14.0
    assert r["C_counts"]["left_censored_lower_bound"] == 1 and r["C_counts"]["with_warning_and_loss"] == 3


def test_left_censored_flag_uses_earliest_possible_alert_week():
    # signal from the first week: the warning comes in week 2, the earliest possible -> censored
    w = fc.attribute(_weekly([0.6] * 4, [-0.6] * 4))
    s = fc.season_summary(w)
    assert bool(s["lead_left_censored"].iloc[0]) is True
    # signal starts in week 3: the warning in week 4 shows real lead -> not censored
    w = fc.attribute(_weekly([1.0, 1.0, 0.6, 0.6, 0.6], [0.0, 0.0, -0.6, -0.6, -0.6]))
    s = fc.season_summary(w)
    assert bool(s["lead_left_censored"].iloc[0]) is False
