"""Recompute Supplementary Fig. S6 using repaired 0.083-degree future projections.

The calculation preserves the equal-area allocation and joint climate-member plus
5-degree spatial-block bootstrap specified in M5.3. It does not alter raw inputs.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from scipy.spatial import cKDTree
from rasterio.transform import from_origin

PROJECT = Path(r"D:\ecoMarine")
REPAIR = PROJECT / "output" / "02_experiments" / "20260816_future_coverage_repair"
PACKAGE = REPAIR / "figureS6_protection_configuration_repaired_20260821"
SOURCE, TABLES, QA = PACKAGE / "source", PACKAGE / "tables", PACKAGE / "qa"
REEF_MPA = PROJECT / "data" / "processed" / "conservation_screening_20260731" / "grid_0083" / "reef_mpa_area_by_hsi_cell_0083_20260731.parquet"
OR10 = PROJECT / "output" / "01_mainline" / "inputs" / "sdm_global_0083_corrected_20260722" / "tables" / "or10_constrained_temporal_mean_1993_2022_20260724.csv"
STATIC_HWEI = PROJECT / "output" / "01_mainline" / "future_risk_ensemble_20260726" / "inputs" / "ptracer_static_hwei_1993_2019" / "hwei_surface_mean_1993_2019_ptracers.nc"
HSI_ROOT = REPAIR / "derived" / "hsi_member_projections"
MESS_ROOT = REPAIR / "derived" / "mess" / "rasters"

MEMBERS = ("gfdl_esm4_r1i1p1f1_gr", "cnrm_esm2_1_r1i1p1f2_gn", "ipsl_cm6a_lr_r1i1p1f1_gn")
TAXA = ("Acropora", "Lobophora", "Scarus", "Cephalopholis")
SCENARIOS = ("ssp126_2050", "ssp245_2050", "ssp585_2050")
REGIONS = {"Global": None, "Caribbean": (-98.0375, -55.0435, 8.0805, 33.0635), "GBR": (141.9985, 156.0255, -24.9535, -9.9305), "Southeast Asia": (89.9575, 159.9265, -10.0135, 25.0125)}
SCENARIO_LABELS = {"ssp126_2050": "SSP126", "ssp245_2050": "SSP245", "ssp585_2050": "SSP585"}
SCHEME_LABELS = {"existing_mpa": "Existing MPA", "climate_priority": "Co-suitability priority", "high_exposure_priority": "High-exposure priority"}
SCORE_EPSILON = 1e-4
SEED = 20260804
N_BOOTSTRAP = 1000
BLOCK_DEG = 5.0
DISPLAY_RES = 0.5
LON_MIN, LON_MAX, LAT_MIN, LAT_MAX = -180., 180., -40., 40.
NODATA = -9999.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def hsi_path(member: str, taxon: str, scenario: str) -> Path:
    return HSI_ROOT / member / "maxent" / "projection" / f"{taxon}_global" / f"{scenario}_hsi_mean.tif"


def mess_path(member: str, taxon: str, scenario: str) -> Path:
    return MESS_ROOT / member / taxon / f"{scenario}_mess_supported_binary.tif"


def exact_budget_selection(scores: np.ndarray, area: np.ndarray, budget: float) -> np.ndarray:
    """Allocate highest scores to the exact reef-area budget; score ties are fractional."""
    selection = np.zeros(scores.size, dtype=np.float64)
    eligible = np.isfinite(scores) & np.isfinite(area) & (area > 0)
    if budget <= 0 or not eligible.any():
        return selection
    available = float(area[eligible].sum(dtype=np.float64))
    if budget >= available * (1 - 1e-12):
        selection[eligible] = 1.
        return selection
    idx = np.flatnonzero(eligible)
    order = idx[np.argsort(-scores[idx], kind="mergesort")]
    cumulative = np.cumsum(area[order], dtype=np.float64)
    cutoff = int(np.searchsorted(cumulative, budget, side="left"))
    score = scores[order[cutoff]]
    above = eligible & (scores > score)
    tied = eligible & (scores == score)
    selected_above = float(area[above].sum(dtype=np.float64))
    tied_area = float(area[tied].sum(dtype=np.float64))
    fraction = (budget - selected_above) / tied_area
    if not (-1e-12 <= fraction <= 1 + 1e-12):
        raise RuntimeError("Invalid fractional allocation at score threshold.")
    selection[above] = 1.
    selection[tied] = np.clip(fraction, 0, 1)
    return selection


def regional_mask(data: pd.DataFrame, bounds: tuple[float, float, float, float] | None) -> np.ndarray:
    if bounds is None:
        return np.ones(len(data), dtype=bool)
    west, east, south, north = bounds
    return ((data.longitude.to_numpy() >= west) & (data.longitude.to_numpy() <= east) & (data.latitude.to_numpy() >= south) & (data.latitude.to_numpy() <= north))


def nearest(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    right = np.searchsorted(source, target)
    right = np.clip(right, 1, len(source) - 1)
    left = right - 1
    return np.where(np.abs(target - source[left]) <= np.abs(target - source[right]), left, right)


def lonlat_to_xyz(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    """Unit-sphere coordinates used only for documented coastal gap filling."""
    lon_rad, lat_rad = np.deg2rad(lon), np.deg2rad(lat)
    cos_lat = np.cos(lat_rad)
    return np.column_stack((cos_lat * np.cos(lon_rad), cos_lat * np.sin(lon_rad), np.sin(lat_rad)))


def static_hwei(data: pd.DataFrame) -> tuple[np.ndarray, dict[str, float | int | str], np.ndarray]:
    """Sample the 1993-2019 index and fill source-land coastal gaps from nearest positive source-ocean cell."""
    with xr.open_dataset(STATIC_HWEI, decode_times=False) as ds:
        field = ds["hwei_surface_mean_1993_2019"].values.astype(float)
        lat = ds.lat.values.astype(float)
        lon = ds.lon.values.astype(float)
    values = field[nearest(lat, data.latitude.to_numpy(float)), nearest(lon, np.mod(data.longitude.to_numpy(float), 360.0))]
    valid_source = np.isfinite(field) & (field > 0)
    source_rows, source_cols = np.nonzero(valid_source)
    source_lon, source_lat = lon[source_cols], lat[source_rows]
    missing = ~np.isfinite(values) | (values <= 0)
    distance_deg = np.zeros(len(data), dtype=float)
    if missing.any():
        tree = cKDTree(lonlat_to_xyz(source_lon, source_lat))
        chord, match = tree.query(lonlat_to_xyz(data.longitude.to_numpy(float)[missing], data.latitude.to_numpy(float)[missing]), k=1)
        values[missing] = field[source_rows[match], source_cols[match]]
        distance_deg[missing] = np.rad2deg(2 * np.arcsin(np.clip(chord / 2, 0, 1)))
    if not np.isfinite(values).all() or (values <= 0).any():
        raise ValueError("The static exposure field could not be completed at all mapped reef cells.")
    return values, {
        "method": "nearest positive 2.5-degree source-ocean cell on a sphere for source-land coastal gaps",
        "filled_cell_count": int(missing.sum()),
        "filled_cell_fraction": float(missing.mean()),
        "filled_reef_area_km2": float(data.loc[missing, "reef_area_km2"].sum()),
        "filled_reef_area_fraction": float(data.loc[missing, "reef_area_km2"].sum() / data.reef_area_km2.sum()),
        "median_fill_distance_degrees": float(np.median(distance_deg[missing])) if missing.any() else 0.0,
        "p95_fill_distance_degrees": float(np.quantile(distance_deg[missing], .95)) if missing.any() else 0.0,
        "max_fill_distance_degrees": float(distance_deg.max()),
    }, missing


def read_values(path: Path, rows: np.ndarray, cols: np.ndarray) -> tuple[np.ndarray, tuple]:
    with rasterio.open(path) as src:
        values = src.read(1).astype(float)
        nodata = src.nodata
        if nodata is not None and np.isfinite(nodata):
            values[np.isclose(values, nodata)] = np.nan
        signature = (src.shape, tuple(src.transform), str(src.crs))
    values[values <= -9990] = np.nan
    return values[rows, cols], signature


def calculate_metrics(area: np.ndarray, mpa: np.ndarray, p: np.ndarray, e_star: np.ndarray, hwei: np.ndarray) -> dict[str, float]:
    budget = float(np.sum(area * mpa, dtype=np.float64))
    total_p = float(np.sum(area * p, dtype=np.float64))
    climate = exact_budget_selection(p, area, budget)
    high = exact_budget_selection(p + SCORE_EPSILON * e_star, area, budget)
    output: dict[str, float] = {"candidate_reef_area_km2": float(area.sum(dtype=np.float64)), "mpa_budget_km2": budget, "total_area_weighted_support_km2": total_p}
    for name, selection in {"existing_mpa": mpa, "climate_priority": climate, "high_exposure_priority": high}.items():
        selected_area = float(np.sum(area * selection, dtype=np.float64))
        captured = float(np.sum(area * p * selection, dtype=np.float64))
        output[f"{name}_selected_area_km2"] = selected_area
        output[f"{name}_captured_co_suitability_fraction"] = captured / total_p if total_p > 0 else np.nan
        output[f"{name}_mean_hwei"] = float(np.sum(area * p * selection * hwei, dtype=np.float64)) / captured if captured > 0 else np.nan
    overlap = float(np.sum(area * np.minimum(climate, high), dtype=np.float64))
    output["high_pressure_vs_climate_area_changed_fraction"] = 1 - overlap / budget if budget > 0 else np.nan
    output["maximum_area_budget_closure_error_km2"] = max(abs(output[f"{key}_selected_area_km2"] - budget) for key in ("existing_mpa", "climate_priority", "high_exposure_priority"))
    return output


def build_current_fields(data: pd.DataFrame) -> tuple[dict[str, np.ndarray], np.ndarray, list[dict[str, str]]]:
    rows, cols = data.row.to_numpy(int), data.col.to_numpy(int)
    thresholds = pd.read_csv(OR10).set_index("species").loc[list(TAXA), "or10"].to_numpy(float)
    masks: dict[str, np.ndarray] = {}
    continuous_585: list[np.ndarray] = []
    inventory: list[dict[str, str]] = []
    for scenario in SCENARIOS:
        member_masks = []
        for member in MEMBERS:
            hsi_taxa, support_taxa, signatures = [], [], []
            for taxon, threshold in zip(TAXA, thresholds, strict=True):
                hp = hsi_path(member, taxon, scenario)
                mp = mess_path(member, taxon, scenario)
                hsi, sig_h = read_values(hp, rows, cols)
                support, sig_m = read_values(mp, rows, cols)
                if sig_h != sig_m:
                    raise ValueError(f"HSI/MESS grid mismatch: {member}/{taxon}/{scenario}")
                hsi_taxa.append(hsi); support_taxa.append(support); signatures += [sig_h]
                inventory += [{"kind": "future_constrained_hsi", "member": member, "taxon": taxon, "scenario": scenario, "path": str(hp)}, {"kind": "mess_support_binary", "member": member, "taxon": taxon, "scenario": scenario, "path": str(mp)}]
            if len(set(signatures)) != 1:
                raise ValueError("Future HSI grids differ within a member/scenario.")
            hsi_stack = np.vstack(hsi_taxa)
            support_stack = np.vstack(support_taxa)
            finite = np.all(np.isfinite(hsi_stack), axis=0) & np.all(np.isfinite(support_stack), axis=0)
            strict = finite & np.all(hsi_stack >= thresholds[:, None], axis=0) & np.all(support_stack == 1, axis=0)
            member_mask = np.full(len(data), np.nan)
            member_mask[finite] = strict[finite].astype(float)
            member_masks.append(member_mask)
            if scenario == "ssp585_2050":
                co = np.full(len(data), np.nan)
                hsi_complete = np.all(np.isfinite(hsi_stack), axis=0)
                co[hsi_complete] = np.prod(np.clip(hsi_stack[:, hsi_complete], 0, 1), axis=0) ** 0.25
                continuous_585.append(co)
        masks[scenario] = np.vstack(member_masks)
    matrix = np.vstack(continuous_585)
    has_member = np.any(np.isfinite(matrix), axis=0)
    median = np.full(len(data), np.nan)
    median[has_member] = np.nanmedian(matrix[:, has_member], axis=0)
    return masks, median, inventory


def bootstrap(data: pd.DataFrame, masks: dict[str, np.ndarray], hwei: np.ndarray, e_star: np.ndarray) -> pd.DataFrame:
    area_all = data.reef_area_km2.to_numpy(float)
    mpa_all = np.clip(data.mpa_fraction_with_point_buffers.to_numpy(float), 0, 1)
    rng = np.random.default_rng(SEED)
    cache: dict[str, dict[str, np.ndarray | int]] = {}
    for region, bounds in REGIONS.items():
        idx = np.flatnonzero(regional_mask(data, bounds))
        lon, lat = data.longitude.iloc[idx].to_numpy(float), data.latitude.iloc[idx].to_numpy(float)
        keys = np.floor((lat + 90) / BLOCK_DEG).astype(int) * 1000 + np.floor((lon + 180) / BLOCK_DEG).astype(int)
        _, inverse = np.unique(keys, return_inverse=True)
        cache[region] = {"idx": idx, "inverse": inverse, "n_blocks": int(inverse.max() + 1)}
    records: list[dict[str, float | int | str]] = []
    for rep in range(N_BOOTSTRAP):
        member_draw = rng.integers(0, len(MEMBERS), len(MEMBERS))
        for region, values in cache.items():
            idx = values["idx"]; inverse = values["inverse"]; n_blocks = values["n_blocks"]
            counts = np.bincount(rng.integers(0, n_blocks, n_blocks), minlength=n_blocks).astype(float)
            multiplicity = counts[inverse]
            for scenario in SCENARIOS:
                matrix = masks[scenario][:, idx]
                p = np.mean(matrix[member_draw], axis=0)
                candidate = np.all(np.isfinite(matrix), axis=0) & (p > 0)
                metrics = calculate_metrics(area_all[idx][candidate] * multiplicity[candidate], mpa_all[idx][candidate], p[candidate], e_star[idx][candidate], hwei[idx][candidate])
                records.append({"replicate": rep, "region": region, "scenario": scenario, **metrics})
        if (rep + 1) % 100 == 0:
            print(f"bootstrap {rep + 1}/{N_BOOTSTRAP}", flush=True)
    return pd.DataFrame(records)


def summarize(point: pd.DataFrame, boot: pd.DataFrame) -> pd.DataFrame:
    keys = ["region", "scenario"]
    metrics = [c for c in point.columns if c not in keys]
    records = []
    for row in point.itertuples(index=False):
        sample = boot[(boot.region == row.region) & (boot.scenario == row.scenario)]
        for metric in metrics:
            values = pd.to_numeric(sample[metric], errors="coerce").dropna().to_numpy(float)
            estimate = float(getattr(row, metric))
            lower, upper = (np.quantile(values, [0.025, 0.975]) if len(values) else (np.nan, np.nan))
            unit = "km2" if metric.endswith("km2") else ("index_units" if metric.endswith("mean_hwei") else "fraction_or_ratio")
            records.append({"region": row.region, "scenario": row.scenario, "metric": metric, "estimate": estimate, "lower_95": float(lower), "upper_95": float(upper), "n_bootstrap": int(len(values)), "unit": unit})
    return pd.DataFrame(records)


def aggregate_display(data: pd.DataFrame, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    display = ((data.longitude >= LON_MIN) & (data.longitude < LON_MAX) & (data.latitude >= LAT_MIN) & (data.latitude <= LAT_MAX)).to_numpy()
    lon, lat = data.longitude.to_numpy(float)[display], data.latitude.to_numpy(float)[display]
    area = data.reef_area_km2.to_numpy(float)[display]
    values = values[display]
    nlon, nlat = int((LON_MAX-LON_MIN)/DISPLAY_RES), int((LAT_MAX-LAT_MIN)/DISPLAY_RES)
    col = np.clip(np.floor((lon-LON_MIN)/DISPLAY_RES).astype(int), 0, nlon-1)
    row = np.clip(np.floor((LAT_MAX-lat)/DISPLAY_RES).astype(int), 0, nlat-1)
    bin_id = row*nlon + col
    valid = np.isfinite(values) & np.isfinite(area) & (area > 0)
    numerator = np.bincount(bin_id[valid], weights=area[valid]*values[valid], minlength=nlat*nlon)
    denominator = np.bincount(bin_id[valid], weights=area[valid], minlength=nlat*nlon)
    output = np.full(nlat*nlon, np.nan); output[denominator > 0] = numerator[denominator > 0]/denominator[denominator > 0]
    reef = np.bincount(bin_id, weights=area, minlength=nlat*nlon).reshape(nlat,nlon)
    return output.reshape(nlat,nlon), reef


def write_tif(path: Path, values: np.ndarray) -> None:
    profile = {"driver":"GTiff", "height":values.shape[0], "width":values.shape[1], "count":1, "dtype":"float32", "crs":"EPSG:4326", "transform":from_origin(LON_MIN,LAT_MAX,DISPLAY_RES,DISPLAY_RES), "nodata":NODATA, "compress":"LZW"}
    with rasterio.open(path,"w",**profile) as dst:
        dst.write(np.where(np.isfinite(values), values, NODATA).astype("float32"),1)


def main() -> None:
    for directory in (SOURCE, TABLES, QA): directory.mkdir(parents=True, exist_ok=True)
    data = pd.read_parquet(REEF_MPA)
    required = {"cell_id","row","col","longitude","latitude","reef_area_km2","mpa_fraction_with_point_buffers"}
    if missing := required.difference(data.columns): raise ValueError(f"Missing reef/MPA fields: {sorted(missing)}")
    data = data.loc[data.reef_area_km2.gt(0)].copy().reset_index(drop=True)
    mpa = np.clip(data.mpa_fraction_with_point_buffers.to_numpy(float), 0, 1)
    if not np.isfinite(mpa).all(): raise ValueError("MPA fractions contain missing values.")
    hwei, hwei_gap_fill, hwei_gap_filled = static_hwei(data)
    if not np.isfinite(hwei).all(): raise ValueError("Mapped reef cells have missing or zero static MPEI-based exposure values.")
    e_star = (hwei - hwei.min()) / (hwei.max() - hwei.min())
    masks, continuous_585, inventory = build_current_fields(data)
    inventory += [{"kind":"reef_mpa_overlay","member":"","taxon":"","scenario":"","path":str(REEF_MPA)}, {"kind":"or10_thresholds","member":"","taxon":"","scenario":"","path":str(OR10)}, {"kind":"static_mpei_based_exposure_1993_2019","member":"","taxon":"","scenario":"","path":str(STATIC_HWEI)}]
    point_records = []
    p_records = []
    for region, bounds in REGIONS.items():
        region_mask_values = regional_mask(data, bounds)
        for scenario in SCENARIOS:
            matrix = masks[scenario]
            valid = np.all(np.isfinite(matrix),axis=0)
            p = np.full(len(data),np.nan); p[valid]=matrix[:,valid].mean(axis=0)
            candidate = region_mask_values & valid & (p > 0)
            if not candidate.any(): raise ValueError(f"No positive strict-support candidates: {region}/{scenario}")
            metrics = calculate_metrics(data.reef_area_km2.to_numpy(float)[candidate],mpa[candidate],p[candidate],e_star[candidate],hwei[candidate])
            if metrics["maximum_area_budget_closure_error_km2"] > 1e-6: raise RuntimeError("Equal-area budget did not close.")
            point_records.append({"region":region,"scenario":scenario,**metrics})
            p_records.append({"region":region,"scenario":scenario,"candidate_cells":int(candidate.sum()),"candidate_reef_area_km2":float(data.reef_area_km2.to_numpy(float)[candidate].sum()),"coastal_gap_filled_candidate_cells":int(hwei_gap_filled[candidate].sum()),"coastal_gap_filled_candidate_reef_area_km2":float(data.reef_area_km2.to_numpy(float)[candidate & hwei_gap_filled].sum()),"coastal_gap_filled_candidate_reef_area_pct":float(100 * data.reef_area_km2.to_numpy(float)[candidate & hwei_gap_filled].sum()/data.reef_area_km2.to_numpy(float)[candidate].sum()),"support_area_weighted_mean":float(np.average(p[candidate],weights=data.reef_area_km2.to_numpy(float)[candidate]) )})
    point = pd.DataFrame(point_records)
    boot = bootstrap(data,masks,hwei,e_star)
    summary = summarize(point,boot)
    if (summary.n_bootstrap < 990).any(): raise RuntimeError("Too few finite bootstrap replicates for a protection metric.")
    point.to_csv(TABLES / "figS6_protection_point_estimates_20260821.csv",index=False,float_format="%.12g")
    boot.to_csv(TABLES / "figS6_protection_joint_bootstrap_20260821.csv",index=False,float_format="%.12g")
    summary.to_csv(TABLES / "figS6_protection_joint_bootstrap_summary_20260821.csv",index=False,float_format="%.12g")
    pd.DataFrame(p_records).to_csv(TABLES / "figS6_candidate_domain_summary_20260821.csv",index=False,float_format="%.12g")
    p585=masks["ssp585_2050"]; valid585=np.all(np.isfinite(p585),axis=0); support585=np.full(len(data),np.nan); support585[valid585]=p585[:,valid585].mean(axis=0)
    mpa_map, reef = aggregate_display(data,mpa)
    hsi_map,_=aggregate_display(data,continuous_585)
    p_map,_=aggregate_display(data,support585)
    write_tif(SOURCE / "figS6a_mpa_reef_fraction_05deg_20260821.tif",mpa_map)
    write_tif(SOURCE / "figS6b_continuous_four_taxon_hsi_ssp585_05deg_20260821.tif",hsi_map)
    write_tif(SOURCE / "figS6_strict_support_fraction_ssp585_05deg_20260821.tif",p_map)
    native=pd.DataFrame({"cell_id":data.cell_id,"longitude":data.longitude,"latitude":data.latitude,"reef_area_km2":data.reef_area_km2,"mpa_fraction":mpa,"static_mpei_based_exposure_index":hwei,"static_mpei_coastal_gap_filled":hwei_gap_filled,"static_exposure_normalised_e_star":e_star,"ssp585_support_fraction":support585,"ssp585_continuous_four_taxon_hsi_median":continuous_585})
    native.to_csv(SOURCE / "figS6_native_display_and_priority_inputs_20260821.csv",index=False,float_format="%.12g")
    pd.DataFrame(inventory).to_csv(SOURCE / "figS6_input_inventory_20260821.csv",index=False)
    qc={"status":"PASS","source_period_static_exposure":"1993-2019","static_exposure_units":"index units (source physical concentration units unavailable)","static_exposure_minmax_index_units":[float(hwei.min()),float(hwei.max())],"static_exposure_coastal_gap_fill":hwei_gap_fill,"mapped_reef_cells":int(len(data)),"bootstrap_replicates":N_BOOTSTRAP,"bootstrap_seed":SEED,"spatial_block_degrees":BLOCK_DEG,"all_metrics_min_finite_replicates":int(summary.n_bootstrap.min()),"maximum_area_budget_closure_error_km2":float(point.maximum_area_budget_closure_error_km2.max()),"selection_score":"P or P + 0.0001 * E_star; P first, exposure breaks ties","raw_inputs_modified":False,"input_hashes":{"reef_mpa":sha256(REEF_MPA),"or10":sha256(OR10),"static_exposure":sha256(STATIC_HWEI)}}
    (QA / "figS6_source_qc_20260821.json").write_text(json.dumps(qc,indent=2),encoding="utf-8")
    print(json.dumps(qc,indent=2))

if __name__ == "__main__": main()

