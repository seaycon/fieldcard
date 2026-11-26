"""fieldcard.core - field health card with stress-cause attribution (core logic).

The module is deliberately independent of any data source. It takes plain
pandas tables, so the same code runs on public STAC data (prototype), on the
NOSTRADAMUS Data Cube / IoT Observatory (project), or on synthetic data (tests).

Inputs
------
obs : one row per parcel and cloud-free observation
      columns: parcel_id, date, ndvi, ndmi, ndre, [crop]
wx  : one row per day (site weather)
      columns: date, precip (mm), et0 (mm, FAO-56 reference ET),
               tmax (degC), tmean (degC), rh (%, daily mean)

Output
------
weekly : one row per parcel and week with indicators, cause scores, status,
         attributed cause, confidence and a suggested next observation.

Method (v0, rule based and transparent on purpose)
--------------------------------------------------
1. Basal crop coefficient from NDVI (linear cover scaling).
2. FAO-56 daily root-zone water balance -> stress coefficient Ks.
3. Each real observation is scored against a same-crop reference curve for its
   own day of year (robust z-score); alerts need two observed weeks in a row.
   A root-zone drop raises a water signal only together with a canopy-moisture
   decline, so station weather alone never flags a field.
4. Cause scores: water deficit, heat, waterlogging, disease-conducive weather.
   A vegetation anomaly that none of these explains is reported as
   "unexplained" and routed to a field check instead of being guessed.

Simplifications that must be calibrated on real sites are marked ASSUMPTION.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import json

import numpy as np
import pandas as pd

__version__ = "0.2.3"

CAUSES = ["water_deficit", "heat", "waterlogging", "disease_conducive"]


# --------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------
@dataclass
class Params:
    # soil / crop (ASSUMPTION: replace with site soil map and crop tables)
    taw_mm: float = 180.0          # total available water in the root zone
    p_depletion: float = 0.55      # FAO-56 depletion fraction, cereals
    kcb_max: float = 1.10          # basal crop coefficient at full cover
    ndvi_soil: float = 0.15        # NDVI of bare soil
    ndvi_full: float = 0.90        # NDVI at full cover
    ke_bare: float = 0.10          # soil evaporation coefficient, bare soil
    dr0_mm: float = 0.0            # depletion on 1 January (0 = field capacity)
    # season analysed
    season_start: str = "03-01"
    season_end: str = "08-15"
    # weather thresholds
    heat_tmax_c: float = 30.0
    heat_doy: tuple = (140, 190)   # sensitive window (flowering, grain fill)
    humid_rh: float = 80.0
    humid_t: tuple = (10.0, 25.0)
    sat_dr_mm: float = 2.0         # root zone counted as saturated below this
    waterlog_excess_mm: float = 40.0   # weekly excess water that counts as heavy
    # decision thresholds
    z_alert: float = -1.0
    ks_drop_alert: float = 0.2     # Ks this far below the usual value for the week
    ks_branch_ndmi_z: float = -0.5 # ...and only if canopy moisture (NDMI z) is at least this low
    persist_weeks: int = 2         # an alert must hold this many weeks in a row
    ndvi_active: float = 0.35      # below this reference NDVI the crop is not assessed
    cause_min_score: float = 0.4
    low_confidence: float = 0.2
    mad_floor: dict = field(default_factory=lambda: {"ndvi": 0.05, "ndmi": 0.05, "ndre": 0.04})


# --------------------------------------------------------------------------
# Spectral indices
# --------------------------------------------------------------------------
def _nd(a, b):
    a = np.asarray(a, dtype="float64")
    b = np.asarray(b, dtype="float64")
    with np.errstate(divide="ignore", invalid="ignore"):
        out = (a - b) / (a + b)
    return np.where(np.isfinite(out), out, np.nan)


def ndvi(nir, red):
    return _nd(nir, red)


def ndmi(nir, swir1):
    return _nd(nir, swir1)


def ndre(nir, rededge):
    return _nd(nir, rededge)


# --------------------------------------------------------------------------
# Crop coefficient and water balance
# --------------------------------------------------------------------------
def kcb_from_ndvi(ndvi_values, p: Params = Params()):
    """Basal crop coefficient from NDVI by linear cover scaling (ASSUMPTION)."""
    cover = (np.asarray(ndvi_values, dtype="float64") - p.ndvi_soil) / (p.ndvi_full - p.ndvi_soil)
    return p.kcb_max * np.clip(cover, 0.0, 1.0)


def water_balance(precip, et0, kcb, p: Params = Params()) -> pd.DataFrame:
    """Daily FAO-56 root-zone depletion.

    Ks is evaluated on the depletion at the start of the day. Water that the
    root zone cannot hold is returned as `excess` (deep percolation / runoff).
    """
    precip = np.asarray(precip, dtype="float64")
    et0 = np.asarray(et0, dtype="float64")
    kcb = np.asarray(kcb, dtype="float64")
    for name, arr in (("precip", precip), ("et0", et0), ("kcb", kcb)):
        if np.isnan(arr).any():  # a gap would silently freeze the depletion at its last value
            raise ValueError(f"water_balance: {int(np.isnan(arr).sum())} missing {name} values")
    n = len(precip)
    taw, raw = p.taw_mm, p.p_depletion * p.taw_mm
    dr = np.empty(n)
    ks = np.empty(n)
    eta = np.empty(n)
    etc = np.empty(n)
    excess = np.empty(n)
    prev = float(np.clip(p.dr0_mm, 0.0, taw))
    for i in range(n):
        k = 1.0 if prev <= raw else max(0.0, (taw - prev) / (taw - raw))
        ke = p.ke_bare * (1.0 - min(1.0, kcb[i] / p.kcb_max))
        etc[i] = (kcb[i] + ke) * et0[i]
        eta[i] = (k * kcb[i] + ke) * et0[i]
        cand = prev - precip[i] + eta[i]
        excess[i] = max(0.0, -cand)
        # at the wilting point the crop cannot extract more water
        if cand > taw:
            eta[i] -= cand - taw
            cand = taw
        prev = max(0.0, cand)
        dr[i], ks[i] = prev, k
    return pd.DataFrame({"dr": dr, "ks": ks, "eta": eta, "etc": etc, "excess": excess})


# --------------------------------------------------------------------------
# Pipeline: observations + weather -> weekly table
# --------------------------------------------------------------------------
INDEX_COLS = ("ndvi", "ndmi", "ndre")


def _week_start(dates: pd.Series) -> pd.Series:
    d = pd.to_datetime(dates).dt.normalize()
    return d - pd.to_timedelta(d.dt.weekday, unit="D")


def _mad(x) -> float:
    x = np.asarray(x, dtype="float64")
    x = x[np.isfinite(x)]
    return float(np.median(np.abs(x - np.median(x)))) if len(x) else np.nan


def _parcel_year_daily(o: pd.DataFrame, wx_year: pd.DataFrame, p: Params) -> pd.DataFrame:
    """Daily table for one parcel and calendar year: indices, water balance, weather flags."""
    days = pd.DatetimeIndex(wx_year["date"])
    o = o.sort_values("date").drop_duplicates("date").set_index("date")
    daily = pd.DataFrame(index=days)
    for col in INDEX_COLS:
        s = o[col].reindex(days.union(o.index)).interpolate(method="time").reindex(days)
        daily[col] = s.bfill().ffill()              # gap-filled, drives the water balance
        daily[f"{col}_obs"] = o[col].reindex(days)  # real observations only, drive the alerts
    wb = water_balance(wx_year["precip"].to_numpy(), wx_year["et0"].to_numpy(),
                       kcb_from_ndvi(daily["ndvi"].to_numpy(), p), p)
    wb.index = days
    daily = daily.join(wb)
    doy = days.dayofyear.to_numpy()
    tmax, tmean, rh = (wx_year[c].to_numpy() for c in ("tmax", "tmean", "rh"))
    daily["precip"] = wx_year["precip"].to_numpy()
    daily["tmax"] = tmax
    daily["hot"] = ((tmax >= p.heat_tmax_c) & (doy >= p.heat_doy[0]) & (doy <= p.heat_doy[1])).astype(int)
    daily["humid_warm"] = ((rh >= p.humid_rh) & (tmean >= p.humid_t[0]) & (tmean <= p.humid_t[1])).astype(int)
    daily["sat"] = (daily["dr"] <= p.sat_dr_mm).astype(int)
    daily["doy"] = doy
    daily["week"] = _week_start(pd.Series(days)).to_numpy()
    return daily.rename_axis("date").reset_index()


REFERENCE_MODES = ("pooled", "leave_one_out")


def _fold_labels(df: pd.DataFrame, baseline_years, mode: str) -> pd.Series:
    """Which reference each row is compared with.

    pooled        : every row uses all baseline rows.
    leave_one_out : rows of a baseline season use the baseline without their own
                    season; with a single baseline season, without their own
                    parcel. Rows outside the baseline use the whole baseline.
                    This keeps a season from being compared with itself.
    """
    if mode not in REFERENCE_MODES:
        raise ValueError(f"reference mode must be one of {REFERENCE_MODES}")
    base_years = sorted(set(baseline_years))
    lab = pd.Series("all", index=df.index, dtype="object")
    if mode == "pooled":
        return lab
    in_base = df["year"].isin(base_years)
    if len(base_years) >= 2:
        lab[in_base] = "y" + df.loc[in_base, "year"].astype(str)
    else:
        lab[in_base] = "p" + df.loc[in_base, "parcel_id"].astype(str)
    return lab


def _fold_base(base: pd.DataFrame, label: str) -> pd.DataFrame:
    if label == "all":
        return base
    if label.startswith("y"):
        return base[base["year"] != int(label[1:])]
    return base[base["parcel_id"].astype(str) != label[1:]]


def _ref_table(base: pd.DataFrame, col: str) -> pd.DataFrame:
    g = base.groupby(["crop", "doy"])[col]
    ref = pd.concat([g.median().rename("med"), g.agg(_mad).rename("mad")], axis=1).reset_index()
    ref = ref.sort_values(["crop", "doy"])
    for c in ("med", "mad"):  # light smoothing of the day-to-day reference
        ref[c] = ref.groupby("crop")[c].transform(lambda x: x.rolling(7, center=True, min_periods=1).median())
    return ref.rename(columns={"med": f"{col}_ref", "mad": f"{col}_mad"})


def _add_reference(daily: pd.DataFrame, baseline_years, p: Params, mode: str = "pooled") -> pd.DataFrame:
    """Reference curve per crop and day of year, and the z-score of each real observation.

    Each observation is compared with the reference on its own day, so sparse
    sampling during fast green-up does not create false anomalies. Calendar
    days are compared across years; a year with an early or late spring shifts
    the whole curve, and aligning on thermal time is the planned fix.
    The reference should describe an unstressed crop, so pass only normal
    (non-drought) seasons as `baseline_years`.
    """
    base = daily[daily["year"].isin(list(baseline_years))]
    labels = _fold_labels(daily, baseline_years, mode)
    parts = []
    for label, rows in daily.groupby(labels, sort=False):
        b = _fold_base(base, label)
        rows = rows.copy()
        for col in INDEX_COLS:
            rows = rows.merge(_ref_table(b, col), on=["crop", "doy"], how="left")
            scale = np.maximum(1.4826 * rows[f"{col}_mad"].to_numpy(), p.mad_floor[col])
            rows[f"{col}_zobs"] = (rows[f"{col}_obs"] - rows[f"{col}_ref"]) / scale
        parts.append(rows)
    out = pd.concat(parts, ignore_index=True)
    return out.sort_values(["parcel_id", "year", "date"]).reset_index(drop=True)


def _add_ks_reference(weekly: pd.DataFrame, baseline_years, mode: str = "pooled") -> pd.DataFrame:
    """Usual water status for this crop and week, with the same folds as the index reference."""
    base = weekly[weekly["year"].isin(list(baseline_years))]
    labels = _fold_labels(weekly, baseline_years, mode)
    parts = []
    for label, rows in weekly.groupby(labels, sort=False):
        b = _fold_base(base, label)
        ks_ref = b.groupby(["crop", "woy"])["ks_min"].median().rename("ks_ref").reset_index()
        parts.append(rows.merge(ks_ref, on=["crop", "woy"], how="left"))
    out = pd.concat(parts, ignore_index=True)
    return out.sort_values(["parcel_id", "year", "week"]).reset_index(drop=True)


def run_pipeline(obs: pd.DataFrame, wx: pd.DataFrame, p: Params = Params(),
                 baseline_years=None, reference_mode: str = "pooled") -> pd.DataFrame:
    """Full chain.

    `baseline_years` defaults to all years present. The reference should
    describe an unstressed crop: pass normal (non-drought) seasons only.
    `reference_mode="leave_one_out"` compares each baseline season with the
    other baseline seasons (or, with one baseline season, each parcel with the
    other parcels), so no season is compared with itself.

    The water balance runs over the whole calendar year from 1 January
    (`p.dr0_mm`, default 0 = field capacity after winter); only the weeks
    between `p.season_start` and `p.season_end` are assessed.
    """
    obs = obs.copy()
    wx = wx.copy()
    obs["date"] = pd.to_datetime(obs["date"]).dt.normalize()
    wx["date"] = pd.to_datetime(wx["date"]).dt.normalize()
    if "crop" not in obs.columns:
        obs["crop"] = "unknown"
    obs["year"] = obs["date"].dt.year
    wx["year"] = wx["date"].dt.year

    parts = []
    for (pid, year), o in obs.groupby(["parcel_id", "year"], sort=True):
        wx_year = wx[wx["year"] == year].sort_values("date")
        if len(wx_year) < 200 or o["ndvi"].notna().sum() < 3:
            continue  # not enough weather or observations for this parcel-year
        d = _parcel_year_daily(o[["date", *INDEX_COLS]], wx_year, p)
        d.insert(0, "parcel_id", pid)
        d["year"] = year
        d["crop"] = o["crop"].iloc[0]
        parts.append(d)
    if not parts:
        raise ValueError("no parcel-year has enough observations and weather data")
    daily = pd.concat(parts, ignore_index=True)
    years = sorted(daily["year"].unique()) if baseline_years is None else list(baseline_years)
    daily = _add_reference(daily, years, p, reference_mode)

    weekly = daily.groupby(["parcel_id", "year", "crop", "week"], sort=True).agg(
        ndvi=("ndvi", "mean"), ndmi=("ndmi", "mean"), ndre=("ndre", "mean"),
        ndvi_z=("ndvi_zobs", "mean"), ndmi_z=("ndmi_zobs", "mean"), ndre_z=("ndre_zobs", "mean"),
        ndvi_ref=("ndvi_ref", "mean"),
        ks_min=("ks", "min"), dr_max=("dr", "max"), eta=("eta", "sum"), etc=("etc", "sum"),
        precip=("precip", "sum"), excess=("excess", "sum"), tmax=("tmax", "max"),
        hot_days=("hot", "sum"), humid_warm_days=("humid_warm", "sum"),
        sat_days=("sat", "sum"), n_obs=("ndvi_obs", "count"),
    ).reset_index()
    weekly["rel_depletion"] = weekly["dr_max"] / p.taw_mm
    start = pd.to_datetime(weekly["year"].astype(str) + "-" + p.season_start)
    end = pd.to_datetime(weekly["year"].astype(str) + "-" + p.season_end)
    weekly = weekly[(weekly["week"] >= start) & (weekly["week"] <= end)].reset_index(drop=True)
    weekly["woy"] = weekly["week"].dt.isocalendar().week.astype(int)

    # consecutive weeks without any optical observation (data gap)
    has = weekly["n_obs"] > 0
    block = has.groupby([weekly["parcel_id"], weekly["year"]]).cumsum()
    weekly["gap_weeks"] = (~has).groupby([weekly["parcel_id"], weekly["year"], block]).cumsum().astype(int)

    # usual water status for this crop and week: late-season depletion is normal
    weekly = _add_ks_reference(weekly, years, reference_mode)
    return attribute(weekly, p)


# --------------------------------------------------------------------------
# Attribution
# --------------------------------------------------------------------------
def _persist(flag: pd.Series, keys, n: int) -> pd.Series:
    """True where the flag held for n consecutive rows of the same parcel-year."""
    out = flag.copy()
    for k in range(1, n):
        out &= flag.groupby(keys).shift(k, fill_value=False).astype(bool)
    return out


def _persist_observed(z: pd.Series, keys: pd.DataFrame, alert: float, n: int, carry: int = 2) -> pd.Series:
    """Anomaly confirmed by n consecutive *observed* weeks, carried over short cloud gaps."""
    seen = z.notna()
    flag = z[seen] <= alert
    pers = _persist(flag, [keys.loc[seen, c] for c in keys.columns], n)
    out = pd.Series(np.nan, index=z.index, dtype="float64")
    out.loc[seen] = pers.astype("float64")
    out = out.groupby([keys[c] for c in keys.columns]).ffill(limit=carry)
    return out.fillna(0.0).astype(bool)


def attribute(weekly: pd.DataFrame, p: Params = Params()) -> pd.DataFrame:
    w = weekly.sort_values(["parcel_id", "year", "week"]).reset_index(drop=True)
    keys = w[["parcel_id", "year"]]
    glist = [w["parcel_id"], w["year"]]
    clip = lambda x: np.clip(x, 0.0, 1.0)
    ndmi_recent = w["ndmi_z"].groupby(glist).ffill(limit=2).fillna(0.0)
    ks_drop = (w["ks_ref"].fillna(1.0) - w["ks_min"]).clip(lower=0.0)
    w["s_water_deficit"] = 0.7 * clip(ks_drop / 0.4) + 0.3 * clip(-ndmi_recent / 2.0)
    w["s_heat"] = clip(w["hot_days"] / 3.0)
    w["s_waterlogging"] = clip(w["sat_days"] / 5.0) * clip(w["excess"] / p.waterlog_excess_mm)
    w["s_disease_conducive"] = clip(w["humid_warm_days"] / 5.0)

    assessed = w["ndvi_ref"].fillna(0.0) >= p.ndvi_active
    veg_anom = (_persist_observed(w["ndvi_z"], keys, p.z_alert, p.persist_weeks)
                | _persist_observed(w["ndre_z"], keys, p.z_alert, p.persist_weeks)) & assessed
    # water signal: a clear canopy-moisture anomaly, or a root-zone drop confirmed by at
    # least a mild canopy-moisture decline (weather alone never raises an alert)
    ks_branch = (ks_drop >= p.ks_drop_alert) & (ndmi_recent <= p.ks_branch_ndmi_z)
    water_sig = (_persist_observed(w["ndmi_z"], keys, p.z_alert, p.persist_weeks)
                 | _persist(ks_branch, glist, p.persist_weeks)) & assessed
    early = water_sig & ~veg_anom
    w["status"] = np.select([~assessed, veg_anom, early],
                            ["not_assessed", "stress", "early_warning"], default="normal")

    s = w[[f"s_{c}" for c in CAUSES]].to_numpy(dtype="float64")
    order = np.argsort(-s, axis=1)
    top = np.take_along_axis(s, order[:, :1], axis=1)[:, 0]
    second = np.take_along_axis(s, order[:, 1:2], axis=1)[:, 0]
    top_name = np.array(CAUSES)[order[:, 0]]
    explained = top >= p.cause_min_score
    alert = w["status"].isin(["stress", "early_warning"]).to_numpy()
    w["cause"] = np.where(~alert, "none", np.where(explained, top_name, "unexplained"))
    w["confidence"] = np.where(alert & explained, top - second, 0.0)
    w["second_cause"] = np.array(CAUSES)[order[:, 1]]
    # weather risk flags are reported even when the crop still looks normal
    w["flag_disease_weather"] = w["s_disease_conducive"] >= 0.6
    w["flag_heat"] = w["s_heat"] >= 0.6
    w["next_observation"] = [next_observation(r, p) for r in w.to_dict("records")]
    return w


def next_observation(row: dict, p: Params = Params()) -> str:
    """Cheapest observation that would reduce the uncertainty of this card."""
    if row["status"] in ("normal", "not_assessed"):
        return "none"
    if row.get("gap_weeks", 0) >= 2:
        return "no optical data for 2+ weeks: check radar series or request a new acquisition"
    if row["cause"] == "unexplained":
        return "field visit to the anomalous zone with photos (pest, disease, nutrient or management cause not visible from weather)"
    if row["confidence"] < p.low_confidence:
        pair = {row["cause"], row["second_cause"]}
        if pair == {"water_deficit", "heat"}:
            return "soil moisture reading at root depth to separate drought from heat"
        if "disease_conducive" in pair:
            return "leaf photos from 3 points to confirm or rule out disease"
        if "waterlogging" in pair:
            return "check standing water and soil structure in low spots"
        return "field check at the weakest zone"
    return "none"


_ACTION = {
    "none": "No action needed.",
    "water_deficit": "Water deficit is the most likely cause. Prioritise irrigation where available and expect yield impact if it persists in the sensitive stage.",
    "heat": "Heat during the sensitive stage is the most likely cause. Little can be corrected now; record for yield expectation.",
    "waterlogging": "Saturated soil is the most likely cause. Check drainage and delay field traffic.",
    "disease_conducive": "Weather has favoured fungal disease. Scout the crop before deciding on treatment.",
    "unexplained": "Vegetation is below normal and weather does not explain it. A field check is needed.",
}


def make_card(row: dict) -> dict:
    """One JSON-serialisable card for one parcel and week."""
    f = lambda v, n=2: None if v is None or (isinstance(v, float) and np.isnan(v)) else round(float(v), n)
    return {
        "parcel_id": str(row["parcel_id"]),
        "week_start": pd.Timestamp(row["week"]).strftime("%Y-%m-%d"),
        "crop": row.get("crop", "unknown"),
        "status": row["status"],
        "cause": row["cause"],
        "confidence": f(row["confidence"]),
        "scores": {c: f(row[f"s_{c}"]) for c in CAUSES},
        "indicators": {
            "ndvi": f(row["ndvi"]), "ndvi_z": f(row["ndvi_z"]),
            "ndmi": f(row["ndmi"]), "ndmi_z": f(row["ndmi_z"]),
            "ks_min": f(row["ks_min"]), "ks_usual": f(row["ks_ref"]), "root_zone_depletion": f(row["rel_depletion"]),
            "precip_mm": f(row["precip"], 1), "eta_mm": f(row["eta"], 1),
        },
        "flags": {"disease_weather": bool(row["flag_disease_weather"]), "heat": bool(row["flag_heat"])},
        "data": {"optical_observations": int(row["n_obs"]), "weeks_without_observation": int(row["gap_weeks"])},
        "action": _ACTION[row["cause"]],
        "next_observation": row["next_observation"],
        "method_version": __version__,
    }


def cards_json(weekly: pd.DataFrame, week=None) -> str:
    w = weekly if week is None else weekly[weekly["week"] == pd.Timestamp(week)]
    return json.dumps([make_card(r) for r in w.to_dict("records")], indent=2)


# --------------------------------------------------------------------------
# Season and year summaries (evidence tables)
# --------------------------------------------------------------------------
def season_summary(weekly: pd.DataFrame, p: Params = Params()) -> pd.DataFrame:
    """One row per parcel-year: how much stress, which cause, how early."""
    rows = []
    for (pid, year), g in weekly.groupby(["parcel_id", "year"], sort=True):
        g = g.sort_values("week")
        first_early = g.loc[g["status"] == "early_warning", "week"].min()
        first_ndvi = g.loc[g["status"] == "stress", "week"].min()
        lead = (first_ndvi - first_early).days if pd.notna(first_early) and pd.notna(first_ndvi) else np.nan
        is_alert = g["status"].isin(["stress", "early_warning"])
        counts = g.loc[is_alert, "cause"].value_counts()
        # first week the card names water deficit (early warning or stress)
        first_water = g.loc[is_alert & (g["cause"] == "water_deficit"), "week"].min()
        lead_water = (first_ndvi - first_water).days if pd.notna(first_water) and pd.notna(first_ndvi) else np.nan
        # a warning at the earliest week an alert can occur cannot show its full lead
        # time: the signal was already present when assessment began (left-censored).
        # An alert needs persist_weeks consecutive weeks, so on a parcel-year whose
        # first row is assessed the earliest alert is row persist_weeks - 1.
        first_assessed = g.loc[g["status"] != "not_assessed", "week"].min()
        k = min(max(p.persist_weeks - 1, 0), len(g) - 1)
        earliest_alert = max(first_assessed, g["week"].iloc[k]) if pd.notna(first_assessed) else pd.NaT
        censored = bool(pd.notna(first_water) and pd.notna(earliest_alert) and first_water <= earliest_alert)
        rows.append({
            "parcel_id": pid, "year": year, "crop": g["crop"].iloc[0],
            "ks_min": g["ks_min"].min(),
            "weeks_alert": int(is_alert.sum()),
            "weeks_water_deficit": int(counts.get("water_deficit", 0)),
            "weeks_unexplained": int(counts.get("unexplained", 0)),
            "main_cause": counts.index[0] if len(counts) else "none",
            "lead_days": lead,
            "first_water_warning": first_water,
            "first_greenness_loss": first_ndvi,
            "lead_water_days": lead_water,
            "first_assessed": first_assessed,
            "earliest_possible_alert": earliest_alert,
            "lead_left_censored": censored,
        })
    return pd.DataFrame(rows)


def year_table(summary: pd.DataFrame) -> pd.DataFrame:
    """One row per year: the comparison used as proof-of-concept evidence."""
    g = summary.groupby("year")
    return pd.DataFrame({
        "parcels": g.size(),
        "median_ks_min": g["ks_min"].median().round(2),
        "pct_parcels_water_deficit_2w": (g["weeks_water_deficit"].apply(lambda x: (x >= 2).mean()) * 100).round(0),
        "median_weeks_water_deficit": g["weeks_water_deficit"].median(),
        "median_lead_days": g["lead_days"].median(),
    }).reset_index()


# --------------------------------------------------------------------------
# Synthetic data (tests and offline demo only - not scientific evidence)
# --------------------------------------------------------------------------
def synthetic_dataset(years=(2018, 2021, 2022), drought_years=(2018, 2022), n_parcels=30,
                      local_stress_frac=0.15, seed=0):
    """Weather + per-parcel index series with known causes.

    Drought years get little spring/summer rain. A share of parcels gets a
    local vegetation loss that weather does not explain (e.g. disease, pest).
    Returns (obs, wx, truth).
    """
    rng = np.random.default_rng(seed)
    obs_rows, wx_parts, truth = [], [], []
    for year in years:
        days = pd.date_range(f"{year}-01-01", f"{year}-12-31", freq="D")
        doy = days.dayofyear.to_numpy()
        season = np.clip(np.sin((doy - 80) / 365 * 2 * np.pi), 0, None)  # 0 winter, 1 midsummer
        dry = year in drought_years
        spring_summer = (doy >= 95) & (doy <= 215)
        wet_prob = np.where(spring_summer & dry, 0.10, 0.36)
        wet = rng.random(len(days)) < wet_prob
        precip = np.where(wet, rng.exponential(5.0, len(days)), 0.0)
        et0 = (0.5 + 4.2 * season) * (1.15 if dry else 1.0) * rng.normal(1, 0.08, len(days)).clip(0.6)
        tmean = 2 + 17 * season + (2.5 if dry else 0) * spring_summer + rng.normal(0, 2.5, len(days))
        tmax = tmean + 6 + rng.normal(0, 1.5, len(days))
        rh = np.clip(88 - 22 * season - (8 if dry else 0) * spring_summer + 10 * wet + rng.normal(0, 4, len(days)), 35, 100)
        wx_parts.append(pd.DataFrame({"date": days, "precip": precip, "et0": et0,
                                      "tmax": tmax, "tmean": tmean, "rh": rh}))
        rise = 1 / (1 + np.exp(-(doy - 105) / 9.0))
        fall = 1 / (1 + np.exp((doy - 200) / 7.0))
        pot = 0.22 + 0.63 * rise * fall                      # unstressed NDVI of a winter cereal
        for i in range(n_parcels):
            prng = np.random.default_rng(seed * 100003 + year * 1009 + i)
            pp = Params(taw_mm=float(prng.uniform(150, 210)))
            wb = water_balance(precip, et0, kcb_from_ndvi(pot, pp), pp)
            deficit = 1 - wb["ks"].to_numpy()
            fast = pd.Series(deficit).ewm(span=5).mean().to_numpy()               # canopy water reacts fast
            slow = pd.Series(deficit).ewm(span=14).mean().shift(12).fillna(0).to_numpy()  # greenness lags
            local = prng.random() < local_stress_frac
            loss = np.where(doy >= 130, 0.16, 0.0) * local
            v = pot - 0.35 * slow * rise - loss
            m = 0.05 + 0.40 * rise * fall - 0.30 * fast * rise - 0.3 * loss
            e = 0.10 + 0.38 * rise * fall - 0.22 * slow * rise - loss
            pick = (np.arange(len(days)) % 5 == (i % 5)) & (prng.random(len(days)) > 0.45)
            for j in np.flatnonzero(pick):
                obs_rows.append((f"P{i:03d}", days[j], v[j] + prng.normal(0, 0.015),
                                 m[j] + prng.normal(0, 0.015), e[j] + prng.normal(0, 0.012), "winter_wheat"))
            truth.append({"parcel_id": f"P{i:03d}", "year": year, "drought_year": dry, "local_stress": bool(local)})
    obs = pd.DataFrame(obs_rows, columns=["parcel_id", "date", "ndvi", "ndmi", "ndre", "crop"])
    return obs, pd.concat(wx_parts, ignore_index=True), pd.DataFrame(truth)
