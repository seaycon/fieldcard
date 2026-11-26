"""Readers for the open data used in the TRL 4 retrospective validation.

Nothing here is needed by fieldcard.core. Heavy libraries are imported inside
the functions, so `import fieldcard` stays light. Install with
`pip install "fieldcard[data]"`.

Sources
-------
DWD Climate Data Center   daily station data (KL), station 00445 Bernburg/Saale (Nord)
Open-Meteo archive         ERA5-Land, used only to fill missing station days
Thuenen Institute          national crop-type maps v201 (vector, FlatGeobuf), Zenodo
Microsoft Planetary Comp.  Sentinel-2 L2A (STAC + cloud-optimised GeoTIFF)
"""
from __future__ import annotations

import io
import re
import zipfile
from datetime import timedelta

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# DWD daily station data
# --------------------------------------------------------------------------
DWD_KL = "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/"


def _dwd_zip_urls(station_id: str, include_recent: bool):
    import requests
    urls = []
    html = requests.get(DWD_KL + "historical/", timeout=60).text
    urls += [DWD_KL + "historical/" + m for m in
             sorted(set(re.findall(rf"tageswerte_KL_{station_id}_\d{{8}}_\d{{8}}_hist\.zip", html)))]
    if include_recent:
        html = requests.get(DWD_KL + "recent/", timeout=60).text
        urls += [DWD_KL + "recent/" + m for m in sorted(set(re.findall(rf"tageswerte_KL_{station_id}_akt\.zip", html)))]
    if not urls:
        raise RuntimeError(f"no DWD KL files found for station {station_id}")
    return urls


def _read_dwd_zip(content: bytes) -> pd.DataFrame:
    with zipfile.ZipFile(io.BytesIO(content)) as z:
        name = next(n for n in z.namelist() if n.startswith("produkt_klima_tag"))
        df = pd.read_csv(z.open(name), sep=";", skipinitialspace=True, na_values=[-999, "-999"])
    df.columns = [c.strip() for c in df.columns]
    df["date"] = pd.to_datetime(df["MESS_DATUM"].astype(str), format="%Y%m%d")
    return df


def dwd_daily_kl(station_id: str = "00445", include_recent: bool = False) -> pd.DataFrame:
    """Daily KL table of one DWD station; quality-checked historical values win over recent ones."""
    import requests
    parts = []
    for url in _dwd_zip_urls(station_id, include_recent):
        r = requests.get(url, timeout=120)
        r.raise_for_status()
        d = _read_dwd_zip(r.content)
        d["source_file"] = url.rsplit("/", 1)[-1]
        parts.append(d)
    df = pd.concat(parts, ignore_index=True)
    df = df.drop_duplicates("date", keep="first").sort_values("date").reset_index(drop=True)
    return df


# --------------------------------------------------------------------------
# FAO-56 reference evapotranspiration (daily)
# --------------------------------------------------------------------------
def _e0(t):
    return 0.6108 * np.exp(17.27 * t / (t + 237.3))


def et0_fao56(tmax, tmin, rh_mean, u_z, sunshine_h, doy, lat_deg, elev_m,
              wind_height_m=10.0, ea=None, a_s=0.25, b_s=0.50):
    """Daily FAO-56 Penman-Monteith ET0 (mm/day), Allen et al. (1998).

    Solar radiation is estimated from sunshine duration with the Angstrom
    formula (a_s, b_s). Wind is reduced from `wind_height_m` to 2 m (eq. 47).
    `ea` (kPa) overrides the vapour pressure derived from mean humidity.
    """
    tmax, tmin = np.asarray(tmax, float), np.asarray(tmin, float)
    doy, n = np.asarray(doy, float), np.asarray(sunshine_h, float)
    u2 = np.asarray(u_z, float) * 4.87 / np.log(67.8 * wind_height_m - 5.42)
    tmean = (tmax + tmin) / 2.0
    pres = 101.3 * ((293.0 - 0.0065 * elev_m) / 293.0) ** 5.26
    gamma = 0.000665 * pres
    es = (_e0(tmax) + _e0(tmin)) / 2.0
    ea = np.asarray(rh_mean, float) / 100.0 * es if ea is None else np.asarray(ea, float)
    delta = 4098.0 * _e0(tmean) / (tmean + 237.3) ** 2
    phi = np.deg2rad(lat_deg)
    dr = 1 + 0.033 * np.cos(2 * np.pi * doy / 365.0)
    dec = 0.409 * np.sin(2 * np.pi * doy / 365.0 - 1.39)
    ws = np.arccos(np.clip(-np.tan(phi) * np.tan(dec), -1, 1))
    ra = 24 * 60 / np.pi * 0.0820 * dr * (ws * np.sin(phi) * np.sin(dec) + np.cos(phi) * np.cos(dec) * np.sin(ws))
    daylen = 24 / np.pi * ws
    rs = (a_s + b_s * np.clip(n / daylen, 0, 1)) * ra
    rso = (0.75 + 2e-5 * elev_m) * ra
    rns = 0.77 * rs
    sigma = 4.903e-9
    rnl = sigma * ((tmax + 273.16) ** 4 + (tmin + 273.16) ** 4) / 2 * (0.34 - 0.14 * np.sqrt(ea)) \
        * (1.35 * np.clip(rs / rso, 0, 1) - 0.35)
    rn = rns - rnl
    et0 = (0.408 * delta * rn + gamma * 900.0 / (tmean + 273.0) * u2 * (es - ea)) / (delta + gamma * (1 + 0.34 * u2))
    return np.maximum(et0, 0.0)


def open_meteo_daily(lat, lon, start, end, variables=("et0_fao_evapotranspiration", "precipitation_sum")):
    """ERA5-Land daily values from the Open-Meteo archive (used only to fill gaps)."""
    import requests
    base = "https://archive-api.open-meteo.com/v1/archive"
    params = {"latitude": lat, "longitude": lon, "start_date": str(start), "end_date": str(end),
              "daily": ",".join(variables), "timezone": "Europe/Berlin", "models": "era5_land"}
    r = requests.get(base, params=params, timeout=120)
    if r.status_code != 200:  # model not available for a variable: fall back to the default reanalysis
        params.pop("models")
        r = requests.get(base, params=params, timeout=120)
    r.raise_for_status()
    d = r.json()["daily"]
    out = pd.DataFrame(d)
    out["date"] = pd.to_datetime(out.pop("time"))
    return out


def weather_table(dwd: pd.DataFrame, lat, lon, elev_m, start, end, wind_height_m=10.0,
                  fill_from_era5=True) -> pd.DataFrame:
    """Daily weather in the form fieldcard.core.run_pipeline expects.

    columns: date, precip, et0, tmax, tmean, rh, filled_precip, filled_et0
    """
    days = pd.date_range(start, end, freq="D")
    d = dwd.set_index("date").reindex(days)
    et0 = et0_fao56(d["TXK"], d["TNK"], d["UPM"], d["FM"], d["SDK"], days.dayofyear, lat, elev_m, wind_height_m)
    out = pd.DataFrame({"date": days, "precip": d["RSK"].to_numpy(), "et0": et0,
                        "tmax": d["TXK"].to_numpy(), "tmean": d["TMK"].to_numpy(), "rh": d["UPM"].to_numpy()})
    out["filled_precip"] = out["precip"].isna()
    out["filled_et0"] = out["et0"].isna()
    if fill_from_era5 and (out["filled_precip"].any() or out["filled_et0"].any()):
        era = open_meteo_daily(lat, lon, days[0].date(), days[-1].date()).set_index("date").reindex(days)
        out.loc[out["filled_et0"].to_numpy(), "et0"] = era["et0_fao_evapotranspiration"].to_numpy()[out["filled_et0"].to_numpy()]
        out.loc[out["filled_precip"].to_numpy(), "precip"] = era["precipitation_sum"].to_numpy()[out["filled_precip"].to_numpy()]
    for c in ("et0", "tmax", "tmean", "rh"):
        out[c] = out[c].interpolate(limit=5, limit_direction="both")
    out["precip"] = out["precip"].fillna(0.0)
    return out


# --------------------------------------------------------------------------
# Thuenen national crop-type maps (vector, FlatGeobuf on Zenodo)
# --------------------------------------------------------------------------
THUENEN_FILES = {
    2017: ("10619783", "CTM_GER_2017_seg_v201.fgb"),
    2018: ("10619783", "CTM_GER_2018_seg_v201.fgb"),
    2019: ("10619783", "CTM_GER_2019_seg_v201.fgb"),
    2020: ("10619783", "CTM_GER_2020_seg_v201.fgb"),
    2021: ("10619783", "CTM_GER_2021_seg_v201.fgb"),
    2022: ("10621629", "CTM_GER_2022_seg_v201.fgb"),
}
THUENEN_LEGEND = "https://zenodo.org/records/10619783/files/CTM_GER_LegendEN_seg_v201.sld"


def thuenen_url(year: int) -> str:
    rec, name = THUENEN_FILES[year]
    return f"https://zenodo.org/records/{rec}/files/{name}"


def _pick_thuenen_key(keys, year: int):
    """File name of the vector map of one year among a record's file names."""
    hits = sorted(k for k in keys if re.fullmatch(rf"CTM_GER_{year}_seg.*\.fgb", k))
    return hits[0] if hits else None


def resolve_thuenen_url(year: int, log=print) -> str:
    """URL of the year's map, with the file name checked against the Zenodo API.

    Falls back to the known name if the API cannot be reached.
    """
    import requests
    rec, name = THUENEN_FILES[year]
    try:
        r = requests.get(f"https://zenodo.org/api/records/{rec}", timeout=60)
        r.raise_for_status()
        keys = [f["key"] for f in r.json().get("files", [])]
        key = _pick_thuenen_key(keys, year)
        if key is None:
            raise FileNotFoundError(f"no CTM_GER_{year}_seg*.fgb in Zenodo record {rec}: {keys}")
        if key != name and log:
            log(f"  {year}: Zenodo file name is {key} (expected {name})")
        return f"https://zenodo.org/records/{rec}/files/{key}"
    except (requests.RequestException, ValueError) as e:
        if log:
            log(f"  {year}: Zenodo API not reachable ({e}); using {name}")
        return thuenen_url(year)


def _norm_code(v) -> str:
    """'4', 4, 4.0 -> '4' so that legend codes and attribute values compare equal."""
    try:
        f = float(v)
        return str(int(f)) if f.is_integer() else str(f)
    except (TypeError, ValueError):
        return str(v).strip()


def detect_class_column(df, codes, min_share=0.8):
    """Attribute whose values are legend codes (the column name may change between years)."""
    codes = {_norm_code(c) for c in codes}
    best, best_share = None, 0.0
    for col in df.columns:
        if col == "geometry":
            continue
        vals = df[col].dropna()
        if not len(vals):
            continue
        share = float(vals.map(_norm_code).isin(codes).mean())
        if share > best_share:
            best, best_share = col, share
    if best is None or best_share < min_share:
        raise KeyError(f"no attribute holds the legend codes (best: {best!r}, {best_share:.0%}); "
                       f"columns: {list(df.columns)}")
    return best


def parse_sld(text: str):
    """(attribute name, table of class code and label) from an SLD legend."""
    rules = re.findall(r"<(?:\w+:)?Rule\b[^>]*>(.*?)</(?:\w+:)?Rule>", text, flags=re.S)
    rows, prop = [], None
    for r in rules:
        title = re.search(r"<(?:\w+:)?Title>(.*?)</(?:\w+:)?Title>", r, flags=re.S) \
            or re.search(r"<(?:\w+:)?Name>(.*?)</(?:\w+:)?Name>", r, flags=re.S)
        pn = re.search(r"<(?:\w+:)?PropertyName>(.*?)</(?:\w+:)?PropertyName>", r, flags=re.S)
        lit = re.search(r"<(?:\w+:)?Literal>(.*?)</(?:\w+:)?Literal>", r, flags=re.S)
        if title and lit:
            rows.append({"code": lit.group(1).strip(), "label": re.sub(r"\s+", " ", title.group(1)).strip()})
            prop = prop or (pn.group(1).strip() if pn else None)
    if not rows:  # raster-style legend
        for m in re.finditer(r"<(?:\w+:)?ColorMapEntry\b([^>]*)/?>", text):
            attrs = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
            if "quantity" in attrs and "label" in attrs:
                rows.append({"code": attrs["quantity"].strip(), "label": attrs["label"].strip()})
    return prop, pd.DataFrame(rows, columns=["code", "label"])


def crop_codes(legend: pd.DataFrame, labels: dict) -> dict:
    """Class codes per crop, matched on the legend label (case-insensitive substring)."""
    out = {}
    low = legend["label"].str.lower()
    for crop, patterns in labels.items():
        hit = np.zeros(len(legend), dtype=bool)
        for pat in patterns:
            hit |= low.str.contains(pat.lower(), regex=False).to_numpy()
        out[crop] = legend.loc[hit, "code"].astype(str).tolist()
    return out


def read_thuenen_fields(year, bbox_lonlat, prop, codes_by_crop, min_area_ha=5.0, inner_buffer_m=20.0,
                        path=None, legend_codes=None, log=print):
    """Fields of the wanted crops inside the box, in EPSG:3035 (equal area).

    The box is transformed to the file's own CRS before the spatial query.
    Fields touching the box edge are dropped, so every field is complete.
    If `prop` is missing from the file, the class attribute is found from the
    values (they must be legend codes); `legend_codes` makes that check strict.
    """
    import geopandas as gpd
    import pyogrio
    from pyproj import Transformer
    from shapely.geometry import box

    path = path or "/vsicurl/" + resolve_thuenen_url(year, log)
    info = pyogrio.read_info(path)
    crs = info["crs"]
    tr = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    native = tr.transform_bounds(*bbox_lonlat, densify_pts=21)
    gdf = pyogrio.read_dataframe(path, bbox=native, fid_as_index=True)
    if prop is None or prop not in gdf.columns:
        wanted = [c for cs in codes_by_crop.values() for c in cs]
        prop = detect_class_column(gdf, legend_codes if legend_codes is not None else wanted,
                                   min_share=0.8 if legend_codes is not None else 0.01)
        if log:
            log(f"  {year}: class attribute detected: {prop}")
    code_to_crop = {_norm_code(c): crop for crop, codes in codes_by_crop.items() for c in codes}
    cls = gdf[prop].map(_norm_code)
    gdf = gdf[cls.isin(code_to_crop)].copy()
    gdf["crop"] = cls[gdf.index].map(code_to_crop)
    gdf = gdf[gdf.within(box(*native))]
    gdf = gdf.set_crs(crs, allow_override=True).to_crs(3035)
    gdf["area_ha"] = gdf.area / 1e4
    gdf = gdf[gdf["area_ha"] >= min_area_ha].copy()
    gdf["geometry"] = gdf.buffer(-abs(inner_buffer_m))
    gdf = gdf[~gdf.is_empty & gdf.is_valid]
    gdf["field_id"] = [f"{year}_{i}" for i in gdf.index]
    gdf["year"] = int(year)
    return gpd.GeoDataFrame(gdf[["field_id", "year", "crop", "area_ha", "geometry"]], crs=3035).reset_index(drop=True)


def sample_fields(gdf, n_per_crop: int, seed: int):
    """Reproducible sample: for each season and crop, numpy default_rng(seed + season)
    draws n fields without replacement from the candidates sorted by field_id."""
    parts = []
    for (year, crop), g in gdf.groupby(["year", "crop"], sort=True):
        g = g.sort_values("field_id").reset_index(drop=True)
        rng = np.random.default_rng(int(seed) + int(year))
        k = min(n_per_crop, len(g))
        parts.append(g.iloc[np.sort(rng.choice(len(g), size=k, replace=False))])
    return pd.concat(parts, ignore_index=True) if parts else gdf.iloc[0:0]


# --------------------------------------------------------------------------
# Sentinel-2 L2A on Microsoft Planetary Computer
# --------------------------------------------------------------------------
PC_STAC = "https://planetarycomputer.microsoft.com/api/stac/v1"
S2_BANDS = ["B04", "B05", "B08", "B11", "SCL"]


def s2_search(bbox_lonlat, start, end, max_cloud=80):
    import planetary_computer
    import pystac_client
    cat = pystac_client.Client.open(PC_STAC, modifier=planetary_computer.sign_inplace)
    search = cat.search(collections=["sentinel-2-l2a"], bbox=bbox_lonlat, datetime=f"{start}/{end}",
                        query={"eo:cloud_cover": {"lt": max_cloud}})
    return list(search.items())


def s2_offset(item) -> int:
    """Digital-number offset of one item, read from its own metadata (never from the date).

    Processing baseline 04.00 and later store reflectance with +1000; the
    returned value (-1000) is added to the digital numbers before dividing by
    10000. If the baseline is missing, the raster:bands offset of B04 is used.
    """
    pb = (getattr(item, "properties", {}) or {}).get("s2:processing_baseline")
    try:
        if pb is not None:
            return -1000 if float(pb) >= 4.0 else 0
    except (TypeError, ValueError):
        pass
    try:
        rb = item.assets["B04"].extra_fields.get("raster:bands", [{}])[0]
        off, scale = rb.get("offset"), rb.get("scale", 1e-4)
        if off is not None:
            return int(round(off / scale))
    except Exception:
        pass
    return 0


def _solar_day(item):
    """Local solar date of an acquisition (as odc-stac groups by 'solar_day')."""
    from shapely.geometry import shape
    lon = shape(item.geometry).centroid.x if item.geometry else 0.0
    return (item.datetime + timedelta(hours=float(lon) / 15.0)).date()


def field_medians(arr: dict, lab, total, off: int, valid_scl=(4, 5), min_valid_frac=0.6) -> pd.DataFrame:
    """Median indices per field label for one day.

    arr   : band -> digital numbers of the pixels inside fields (same order as lab)
    lab   : field label (1..n) of each pixel; total: pixel count per label
    off   : digital-number offset of the day (see s2_offset)
    Fields with less than `min_valid_frac` clear pixels (SCL vegetation or bare
    soil) are dropped.
    """
    valid = np.isin(arr["SCL"], valid_scl) & (arr["B08"] > 0)
    if not valid.any():
        return pd.DataFrame(columns=["ndvi", "ndmi", "ndre", "valid_frac"])
    refl = {b: np.clip((np.asarray(arr[b])[valid].astype("float32") + off) / 10000.0, 1e-4, None)
            for b in ("B04", "B05", "B08", "B11")}
    nir = refl["B08"]
    px = pd.DataFrame({"label": np.asarray(lab)[valid],
                       "ndvi": (nir - refl["B04"]) / (nir + refl["B04"]),
                       "ndmi": (nir - refl["B11"]) / (nir + refl["B11"]),
                       "ndre": (nir - refl["B05"]) / (nir + refl["B05"])})
    med = px.groupby("label").median()
    frac = np.bincount(np.asarray(lab)[valid], minlength=len(total)) / np.maximum(total, 1)
    med["valid_frac"] = frac[med.index.to_numpy()]
    return med[med["valid_frac"] >= min_valid_frac]


def extract_field_indices(items, fields, resolution=20, crs="EPSG:32632", valid_scl=(4, 5),
                          min_valid_frac=0.6, log=print):
    """Per-field median NDVI, NDMI and NDRE for every acquisition day.

    The box around all fields is read once per day (not once per field); fields
    are rasterised onto that grid and summarised in memory.
    """
    import odc.stac
    import rasterio.features

    if not items:
        return pd.DataFrame(columns=["parcel_id", "date", "ndvi", "ndmi", "ndre", "crop", "valid_frac"])
    f = fields.to_crs(crs).reset_index(drop=True)
    bbox = tuple(fields.to_crs(4326).total_bounds)
    offsets = {}
    for it in items:
        offsets.setdefault(_solar_day(it), []).append(s2_offset(it))
    offsets = {k: max(set(v), key=v.count) for k, v in offsets.items()}

    ds = odc.stac.load(items, bands=S2_BANDS, bbox=bbox, crs=crs, resolution=resolution,
                       groupby="solar_day", chunks={"x": 2048, "y": 2048},
                       resampling={"*": "average", "SCL": "nearest"})
    gbox = ds.odc.geobox
    labels = rasterio.features.rasterize(((g, i + 1) for i, g in enumerate(f.geometry)),
                                         out_shape=gbox.shape, transform=gbox.affine, fill=0, dtype="int32")
    sel = labels > 0
    lab = labels[sel]
    total = np.bincount(lab, minlength=len(f) + 1)
    rows = []
    for t in range(ds.sizes["time"]):
        day = pd.Timestamp(ds.time.values[t]).normalize()
        off = offsets.get(day.date(), 0)
        sub = ds.isel(time=t).compute()          # all bands of the day in one parallel read
        arr = {b: sub[b].values[sel] for b in S2_BANDS}
        med = field_medians(arr, lab, total, off, valid_scl, min_valid_frac)
        for lb, r in med.iterrows():
            rows.append((f.at[lb - 1, "field_id"], day, r["ndvi"], r["ndmi"], r["ndre"],
                         f.at[lb - 1, "crop"], r["valid_frac"]))
        if log and t % 10 == 0:
            log(f"  {day.date()}  offset {off:+d}  fields kept {len(med)}/{len(f)}")
    return pd.DataFrame(rows, columns=["parcel_id", "date", "ndvi", "ndmi", "ndre", "crop", "valid_frac"])
