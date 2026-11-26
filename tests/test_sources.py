"""Offline tests for the data readers (no network)."""
from types import SimpleNamespace

import numpy as np
import pytest

from fieldcard.sources import crop_codes, et0_fao56, parse_sld, s2_offset


def test_et0_fao56_example_18():
    # FAO-56 Example 18: Brussels, 6 July, 50 deg 48' N, 100 m; ET0 = 3.9 mm/day
    et0 = et0_fao56(tmax=21.5, tmin=12.3, rh_mean=None, u_z=10 / 3.6, sunshine_h=9.25, doy=187,
                    lat_deg=50.8, elev_m=100, wind_height_m=10, ea=1.409)
    assert float(et0) == pytest.approx(3.9, abs=0.1)


def test_s2_offset_from_metadata_not_date():
    new = SimpleNamespace(properties={"s2:processing_baseline": "05.00"}, assets={})
    old = SimpleNamespace(properties={"s2:processing_baseline": "02.07"}, assets={})
    assert s2_offset(new) == -1000      # reprocessed archive scenes also carry the offset
    assert s2_offset(old) == 0
    band = SimpleNamespace(extra_fields={"raster:bands": [{"offset": -0.1, "scale": 0.0001}]})
    fallback = SimpleNamespace(properties={}, assets={"B04": band})
    assert s2_offset(fallback) == -1000


SLD_RULES = """<StyledLayerDescriptor><NamedLayer><UserStyle><FeatureTypeStyle>
<se:Rule><se:Name>Winter wheat</se:Name><ogc:Filter><ogc:PropertyIsEqualTo>
<ogc:PropertyName>class</ogc:PropertyName><ogc:Literal>101</ogc:Literal></ogc:PropertyIsEqualTo></ogc:Filter></se:Rule>
<se:Rule><se:Name>Winter barley</se:Name><ogc:Filter><ogc:PropertyIsEqualTo>
<ogc:PropertyName>class</ogc:PropertyName><ogc:Literal>103</ogc:Literal></ogc:PropertyIsEqualTo></ogc:Filter></se:Rule>
<se:Rule><se:Name>Maize</se:Name><ogc:Filter><ogc:PropertyIsEqualTo>
<ogc:PropertyName>class</ogc:PropertyName><ogc:Literal>4</ogc:Literal></ogc:PropertyIsEqualTo></ogc:Filter></se:Rule>
</FeatureTypeStyle></UserStyle></NamedLayer></StyledLayerDescriptor>"""


def test_parse_sld_and_crop_codes():
    prop, legend = parse_sld(SLD_RULES)
    assert prop == "class"
    codes = crop_codes(legend, {"winter_wheat": ["winter wheat"], "maize": ["maize"]})
    assert codes == {"winter_wheat": ["101"], "maize": ["4"]}


def test_parse_sld_colormap():
    txt = '<ColorMap><ColorMapEntry color="#ff0" quantity="4" label="Maize"/></ColorMap>'
    prop, legend = parse_sld(txt)
    assert prop is None and legend.iloc[0].to_dict() == {"code": "4", "label": "Maize"}


def test_field_medians_offset_and_cloud_mask():
    from fieldcard.sources import field_medians
    # field 1: 4 clear pixels; field 2: 3 of 4 pixels cloudy (SCL 9) -> dropped at 0.6
    lab = np.array([1, 1, 1, 1, 2, 2, 2, 2])
    total = np.bincount(lab, minlength=3)
    dn = lambda refl: np.full(8, refl * 10000 + 1000)        # baseline >= 04.00: +1000 stored
    arr = {"B04": dn(0.05), "B05": dn(0.15), "B08": dn(0.45), "B11": dn(0.20),
           "SCL": np.array([4, 4, 4, 4, 4, 9, 9, 9])}
    med = field_medians(arr, lab, total, off=-1000)
    assert list(med.index) == [1]
    assert med.loc[1, "ndvi"] == pytest.approx((0.45 - 0.05) / (0.45 + 0.05), abs=1e-4)
    assert med.loc[1, "valid_frac"] == 1.0
    wrong = field_medians(arr, lab, total, off=0)            # forgetting the offset biases NDVI
    assert wrong.loc[1, "ndvi"] < med.loc[1, "ndvi"] - 0.1


def test_read_and_sample_fields_from_local_flatgeobuf(tmp_path):
    gpd = pytest.importorskip("geopandas")
    from shapely.geometry import box
    from fieldcard.sources import read_thuenen_fields, sample_fields
    # fields in EPSG:3035 near Bernburg; the box query must be transformed to the file CRS
    x0, y0 = 4_400_000, 3_190_000
    geoms, cls = [], []
    for i in range(6):
        geoms.append(box(x0 + i * 400, y0, x0 + i * 400 + 300, y0 + 300))   # 9 ha each
        cls.append(101 if i % 2 == 0 else 4)
    geoms.append(box(x0, y0 + 1000, x0 + 100, y0 + 1100))                    # 1 ha: too small
    cls.append(101)
    gdf = gpd.GeoDataFrame({"class": cls}, geometry=geoms, crs=3035)
    path = str(tmp_path / "ctm.fgb")
    gdf.to_file(path, driver="FlatGeobuf")
    bbox = tuple(gdf.to_crs(4326).total_bounds + np.array([-0.05, -0.05, 0.05, 0.05]))
    f = read_thuenen_fields(2018, bbox, "class", {"winter_wheat": ["101"], "maize": ["4"]},
                            min_area_ha=5, inner_buffer_m=20, path=path)
    assert len(f) == 6 and set(f["crop"]) == {"winter_wheat", "maize"}
    assert f["area_ha"].round(1).eq(9.0).all()
    s1 = sample_fields(f, n_per_crop=2, seed=42)
    s2 = sample_fields(f, n_per_crop=2, seed=42)
    assert s1["field_id"].tolist() == s2["field_id"].tolist() and len(s1) == 4


def test_pick_thuenen_key():
    from fieldcard.sources import _pick_thuenen_key
    keys = ["CTM_GER_2017-21_vector_accuracy_v201.pdf", "CTM_GER_2018_seg_v201.fgb", "CTM_GER_LegendEN_seg_v201.sld"]
    assert _pick_thuenen_key(keys, 2018) == "CTM_GER_2018_seg_v201.fgb"
    assert _pick_thuenen_key(keys, 2019) is None


def test_class_column_detected_when_name_differs(tmp_path):
    gpd = pytest.importorskip("geopandas")
    from shapely.geometry import box
    from fieldcard.sources import read_thuenen_fields
    x0, y0 = 4_400_000, 3_190_000
    geoms = [box(x0 + i * 400, y0, x0 + i * 400 + 300, y0 + 300) for i in range(4)]
    gdf = gpd.GeoDataFrame({"OBJ_NAME": ["a", "b", "c", "d"], "CTM_2018": [101.0, 4.0, 101.0, 60.0]},
                           geometry=geoms, crs=3035)
    path = str(tmp_path / "ctm.fgb")
    gdf.to_file(path, driver="FlatGeobuf")
    bbox = tuple(gdf.to_crs(4326).total_bounds + np.array([-0.05, -0.05, 0.05, 0.05]))
    f = read_thuenen_fields(2018, bbox, None, {"winter_wheat": ["101"], "maize": ["4"]}, path=path,
                            legend_codes=["101", "4", "60"], log=None)
    assert sorted(f["crop"]) == ["maize", "winter_wheat", "winter_wheat"]
