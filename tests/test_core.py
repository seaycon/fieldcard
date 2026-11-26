"""Tests for fieldcard. Run:  python -m pytest -q

These check that the code does what the method says. They use synthetic data
generated with the same water-balance model, so they are NOT scientific
validation; that comes from the real-data run in the notebook.
"""
import json

import numpy as np
import pandas as pd
import pytest

import fieldcard as fc


def test_indices():
    assert fc.ndvi(0.5, 0.1) == pytest.approx(0.6667, abs=1e-3)
    assert np.isnan(fc.ndvi(0.0, 0.0))


def test_kcb_bounds_and_monotonic():
    p = fc.Params()
    k = fc.kcb_from_ndvi(np.linspace(-0.2, 1.0, 50), p)
    assert k.min() == 0 and k.max() == pytest.approx(p.kcb_max)
    assert np.all(np.diff(k) >= 0)


def test_dry_down_reaches_full_stress():
    p = fc.Params(taw_mm=100)
    n = 200
    wb = fc.water_balance(np.zeros(n), np.full(n, 5.0), np.full(n, 1.0), p)
    assert np.all(np.diff(wb["dr"]) >= -1e-9)
    assert wb["dr"].max() <= p.taw_mm + 1e-9
    assert wb["ks"].iloc[0] == 1.0 and wb["ks"].iloc[-1] < 0.05


def test_rain_refills_and_excess_is_reported():
    p = fc.Params(taw_mm=100, dr0_mm=80)
    wb = fc.water_balance(np.array([200.0, 0.0]), np.array([3.0, 3.0]), np.array([1.0, 1.0]), p)
    assert wb["dr"].iloc[0] == 0.0
    assert wb["excess"].iloc[0] > 100
    assert wb["ks"].iloc[1] == 1.0


def test_mass_balance_closes():
    rng = np.random.default_rng(1)
    n = 365
    precip = np.where(rng.random(n) < 0.3, rng.exponential(6, n), 0)
    p = fc.Params(dr0_mm=20)
    wb = fc.water_balance(precip, rng.uniform(0.5, 5, n), rng.uniform(0.1, 1.1, n), p)
    lhs = wb["dr"].iloc[-1]
    rhs = p.dr0_mm - precip.sum() + wb["eta"].sum() + wb["excess"].sum()
    assert lhs == pytest.approx(rhs, abs=1e-6)


@pytest.fixture(scope="module")
def result():
    obs, wx, truth = fc.synthetic_dataset(seed=3)
    weekly = fc.run_pipeline(obs, wx, baseline_years=[2021])
    return weekly, fc.season_summary(weekly), truth


def test_drought_years_separate_from_normal_year(result):
    _, summary, _ = result
    yt = fc.year_table(summary).set_index("year")
    assert yt.loc[2018, "pct_parcels_water_deficit_2w"] >= 80
    assert yt.loc[2022, "pct_parcels_water_deficit_2w"] >= 80
    assert yt.loc[2021, "pct_parcels_water_deficit_2w"] <= 20
    assert yt.loc[2018, "median_ks_min"] < yt.loc[2021, "median_ks_min"]


def test_early_warning_precedes_ndvi_anomaly(result):
    _, summary, _ = result
    lead = summary.loc[summary["year"].isin([2018, 2022]), "lead_days"].dropna()
    assert len(lead) > 10
    assert lead.median() >= 7


def test_local_loss_in_normal_year_is_not_blamed_on_weather(result):
    weekly, _, truth = result
    local = truth[(truth["year"] == 2021) & truth["local_stress"]]["parcel_id"]
    assert len(local) >= 2
    w = weekly[(weekly["year"] == 2021) & weekly["parcel_id"].isin(local) & (weekly["status"] == "stress")]
    assert len(w) > 0
    assert (w["cause"] == "unexplained").mean() >= 0.6
    assert (w["next_observation"].str.startswith("field visit")).mean() >= 0.6


def test_bare_or_harvested_weeks_are_not_assessed(result):
    weekly, _, _ = result
    low = weekly[weekly["ndvi_ref"] < fc.Params().ndvi_active]
    assert len(low) > 0 and set(low["status"]) == {"not_assessed"}
    assert set(low["cause"]) == {"none"}


def test_cards_are_valid_json(result):
    weekly, _, _ = result
    week = weekly["week"].iloc[len(weekly) // 2]
    cards = json.loads(fc.cards_json(weekly, week))
    assert cards and {"status", "cause", "scores", "next_observation"} <= set(cards[0])
    assert all(c["status"] in {"normal", "early_warning", "stress", "not_assessed"} for c in cards)


def test_missing_crop_column_and_short_series():
    obs, wx, _ = fc.synthetic_dataset(years=(2021,), drought_years=(), n_parcels=3, seed=5)
    obs = obs.drop(columns="crop")
    extra = pd.DataFrame({"parcel_id": ["X"], "date": [pd.Timestamp("2021-05-01")],
                          "ndvi": [0.5], "ndmi": [0.2], "ndre": [0.3]})
    weekly = fc.run_pipeline(pd.concat([obs, extra], ignore_index=True), wx)
    assert "X" not in set(weekly["parcel_id"])          # too few observations, skipped
    assert set(weekly["crop"]) == {"unknown"}
