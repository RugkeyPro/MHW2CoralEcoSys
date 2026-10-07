#!/usr/bin/env python3
"""Recompute equal-area reef priorities with mean 2045-2055 36-tracer MPEI."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree


WORKTREE = Path(__file__).resolve().parents[4]
OUTPUT_ROOT = Path(__file__).resolve().parents[1]

MPEI_FILE = (
    WORKTREE
    / "output"
    / "01_mainline"
    / "inputs"
    / "mpei_36tracer_annual_1950_2060_20260901"
    / "derived"
    / "mpei_36tracer_0_100m_mean_2045_2055_20260901.nc"
)
SCENARIOS = ("ssp126_2050", "ssp245_2050", "ssp585_2050")
MEMBERS = (
    "gfdl_esm4_r1i1p1f1_gr",
    "cnrm_esm2_1_r1i1p1f2_gn",
    "ipsl_cm6a_lr_r1i1p1f1_gn",
)
REGIONS = {
    "Global": None,
    "Caribbean": (-98.0375, -55.0435, 8.0805, 33.0635),
    "GBR": (141.9985, 156.0255, -24.9535, -9.9305),
    "Southeast Asia": (89.9575, 159.9265, -10.0135, 25.0125),
}
SCORE_EPSILON = 1e-4
N_BOOTSTRAP = 1000
SPATIAL_BLOCK_DEG = 5.0
SEED = 20260804


def lonlat_to_xyz(longitude: np.ndarray, latitude: np.ndarray) -> np.ndarray:
    lon_rad = np.deg2rad(longitude)
    lat_rad = np.deg2rad(latitude)
    cos_lat = np.cos(lat_rad)
    return np.column_stack(
        [cos_lat * np.cos(lon_rad), cos_lat * np.sin(lon_rad), np.sin(lat_rad)]
    )


def read_native_mpei() -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray], dict]:
    with Dataset(MPEI_FILE, "r") as source:
        latitude = np.asarray(source.variables["lat"][:], dtype=np.float64)
        longitude = np.asarray(source.variables["lon"][:], dtype=np.float64)
        fields = {
            "mpei_total": np.asarray(
                np.ma.filled(source.variables["mpei_0_100m"][:], np.nan),
                dtype=np.float64,
            )
        }
        for group, coordinate, variable in (
            ("polymer", "polymer", "mpei_by_polymer_0_100m"),
            ("size", "size_class", "mpei_by_size_class_0_100m"),
            ("state", "state", "mpei_by_transport_state_0_100m"),
        ):
            labels = [str(value) for value in source.variables[coordinate][:].tolist()]
            values = np.asarray(
                np.ma.filled(source.variables[variable][:], np.nan), dtype=np.float64
            )
            for index, label in enumerate(labels):
                fields[f"mpei_{group}_{label}"] = values[index]
        metadata = {
            "completed_source_year_start": int(source.completed_source_year_start),
            "completed_source_year_end": int(source.completed_source_year_end),
            "annual_endpoint_count": int(source.annual_endpoint_count),
            "temporal_aggregation": str(source.temporal_aggregation),
            "physical_forcing_note": str(source.physical_forcing_note),
            "source_sha256": str(source.source_sha256),
        }
    return latitude, longitude, fields, metadata


def interpolate_to_reef_cells(
    reef: pd.DataFrame,
    latitude: np.ndarray,
    longitude: np.ndarray,
    fields: dict[str, np.ndarray],
) -> tuple[pd.DataFrame, dict]:
    target_latitude = reef["latitude"].to_numpy(dtype=np.float64)
    target_longitude = np.mod(reef["longitude"].to_numpy(dtype=np.float64), 360.0)
    points = np.column_stack([target_latitude, target_longitude])
    total = fields["mpei_total"]
    native_valid = np.isfinite(total)
    extended_longitude = np.concatenate(
        ([longitude[-1] - 360.0], longitude, [longitude[0] + 360.0])
    )
    extended_valid = np.concatenate(
        [native_valid[:, -1:], native_valid, native_valid[:, :1]], axis=1
    ).astype(np.float64)
    valid_interpolator = RegularGridInterpolator(
        (latitude, extended_longitude),
        extended_valid,
        method="linear",
        bounds_error=False,
        fill_value=0.0,
    )
    interpolated_valid_weight = valid_interpolator(points)
    output = pd.DataFrame(index=reef.index)
    for name, field in fields.items():
        extended_field = np.concatenate(
            [field[:, -1:], field, field[:, :1]], axis=1
        )
        numerator_interpolator = RegularGridInterpolator(
            (latitude, extended_longitude),
            np.where(np.isfinite(extended_field), extended_field, 0.0),
            method="linear",
            bounds_error=False,
            fill_value=0.0,
        )
        numerator = numerator_interpolator(points)
        values = np.full(len(reef), np.nan, dtype=np.float64)
        np.divide(
            numerator,
            interpolated_valid_weight,
            out=values,
            where=interpolated_valid_weight > 1e-12,
        )
        output[name] = values

    initially_missing = ~np.isfinite(output["mpei_total"].to_numpy())
    fill_flag = np.zeros(len(reef), dtype=bool)
    unsupported = np.zeros(len(reef), dtype=bool)
    fill_distance = np.zeros(len(reef), dtype=np.float64)
    native_grid_diagonal_deg = float(
        np.hypot(np.max(np.diff(latitude)), np.max(np.diff(longitude)))
    )
    if np.any(initially_missing):
        native_rows, native_cols = np.nonzero(native_valid)
        tree = cKDTree(
            lonlat_to_xyz(longitude[native_cols], latitude[native_rows])
        )
        chord, nearest_index = tree.query(
            lonlat_to_xyz(
                reef.loc[initially_missing, "longitude"].to_numpy(dtype=float),
                reef.loc[initially_missing, "latitude"].to_numpy(dtype=float),
            ),
            k=1,
        )
        angular_distance = np.rad2deg(
            2.0 * np.arcsin(np.clip(chord / 2.0, 0.0, 1.0))
        )
        fill_distance[initially_missing] = angular_distance
        acceptable = angular_distance <= native_grid_diagonal_deg
        missing_indices = np.flatnonzero(initially_missing)
        fill_indices = missing_indices[acceptable]
        unsupported_indices = missing_indices[~acceptable]
        fill_flag[fill_indices] = True
        unsupported[unsupported_indices] = True
        selected_rows = native_rows[nearest_index[acceptable]]
        selected_cols = native_cols[nearest_index[acceptable]]
        for name, field in fields.items():
            output.loc[fill_indices, name] = field[selected_rows, selected_cols]

    supported = np.isfinite(output["mpei_total"].to_numpy(dtype=float))
    if np.any(~np.isfinite(output.loc[supported].to_numpy(dtype=float))):
        raise RuntimeError("Supported future MPEI reef cells contain non-finite components")
    if np.nanmin(output.to_numpy(dtype=float)) < -1e-15:
        raise RuntimeError("Future MPEI reef interpolation contains negative values")

    component_groups = {
        "polymer": [name for name in output if name.startswith("mpei_polymer_")],
        "size": [name for name in output if name.startswith("mpei_size_")],
        "state": [name for name in output if name.startswith("mpei_state_")],
    }
    closure = {}
    for group, columns in component_groups.items():
        difference = (
            output.loc[supported, columns].sum(axis=1)
            - output.loc[supported, "mpei_total"]
        )
        closure[group] = float(np.max(np.abs(difference)))
    if max(closure.values()) > 5e-14:
        raise RuntimeError(f"Interpolated component closure failed: {closure}")

    interpolation_qc = {
        "method": (
            "coordinate-aware bilinear interpolation on the irregular IGSM latitude "
            "axis with cyclic longitude; missing coastal reef cells filled from the "
            "nearest valid native ocean-cell centre on a sphere"
        ),
        "reef_cell_count": int(len(reef)),
        "bilinear_valid_cell_count": int((~initially_missing).sum()),
        "filled_cell_count": int(fill_flag.sum()),
        "filled_cell_fraction": float(fill_flag.mean()),
        "filled_reef_area_fraction": float(
            reef.loc[fill_flag, "reef_area_km2"].sum() / reef["reef_area_km2"].sum()
        ),
        "maximum_permitted_fill_distance_degrees": native_grid_diagonal_deg,
        "maximum_accepted_fill_distance_degrees": float(
            fill_distance[fill_flag].max() if fill_flag.any() else 0.0
        ),
        "p95_fill_distance_degrees": float(
            np.quantile(fill_distance[fill_flag], 0.95) if fill_flag.any() else 0.0
        ),
        "unsupported_cell_count": int(unsupported.sum()),
        "unsupported_cell_fraction": float(unsupported.mean()),
        "unsupported_reef_area_fraction": float(
            reef.loc[unsupported, "reef_area_km2"].sum()
            / reef["reef_area_km2"].sum()
        ),
        "unsupported_rule": (
            "exclude cells whose nearest valid native MPEI cell is farther than one "
            "native-grid diagonal"
        ),
        "component_closure_max_abs_kg_m3": closure,
    }
    output["mpei_fill_flag"] = fill_flag
    output["mpei_fill_distance_deg"] = fill_distance
    output["mpei_supported"] = supported
    output["mpei_unsupported"] = unsupported
    return output, interpolation_qc


def exact_budget_selection(scores: np.ndarray, area: np.ndarray, budget: float) -> np.ndarray:
    selection = np.zeros(scores.size, dtype=np.float64)
    eligible = np.isfinite(scores) & np.isfinite(area) & (area > 0.0)
    available = float(np.sum(area[eligible], dtype=np.float64))
    if budget <= 0.0 or available <= 0.0:
        return selection
    if budget >= available * (1.0 - 1e-12):
        selection[eligible] = 1.0
        return selection
    indices = np.flatnonzero(eligible)
    ordered = indices[np.argsort(-scores[indices], kind="mergesort")]
    cumulative = np.cumsum(area[ordered], dtype=np.float64)
    cutoff = int(np.searchsorted(cumulative, budget, side="left"))
    cutoff_score = scores[ordered[cutoff]]
    above = eligible & (scores > cutoff_score)
    tied = eligible & (scores == cutoff_score)
    selected_above = float(np.sum(area[above], dtype=np.float64))
    tied_area = float(np.sum(area[tied], dtype=np.float64))
    selection[above] = 1.0
    selection[tied] = np.clip((budget - selected_above) / tied_area, 0.0, 1.0)
    return selection


def ratio(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator > 0.0 else np.nan


def calculate_metrics(
    area: np.ndarray,
    existing_mpa: np.ndarray,
    suitability_support: np.ndarray,
    exposure_score: np.ndarray,
    mpei: np.ndarray,
    mpei_supported: np.ndarray,
) -> tuple[dict[str, float], dict[str, np.ndarray]]:
    budget = float(np.sum(area * existing_mpa, dtype=np.float64))
    total_suitability = float(np.sum(area * suitability_support, dtype=np.float64))
    climate = exact_budget_selection(suitability_support, area, budget)
    exposure = exact_budget_selection(
        suitability_support + SCORE_EPSILON * exposure_score, area, budget
    )
    selections = {
        "existing_mpa": existing_mpa,
        "climate_priority": climate,
        "future_mpei_priority": exposure,
    }
    output = {
        "candidate_reef_area_km2": float(np.sum(area, dtype=np.float64)),
        "mpa_budget_km2": budget,
        "total_area_weighted_suitability_support": total_suitability,
    }
    for name, selection in selections.items():
        selected_area = float(np.sum(area * selection, dtype=np.float64))
        captured = float(
            np.sum(area * suitability_support * selection, dtype=np.float64)
        )
        supported_selected_area = float(
            np.sum(area * selection * mpei_supported, dtype=np.float64)
        )
        supported_captured = float(
            np.sum(
                area * suitability_support * selection * mpei_supported,
                dtype=np.float64,
            )
        )
        exposure_numerator = float(
            np.sum(
                area
                * suitability_support
                * selection
                * np.where(mpei_supported, mpei, 0.0),
                dtype=np.float64,
            )
        )
        output[f"{name}_selected_area_km2"] = selected_area
        output[f"{name}_captured_co_suitability_fraction"] = ratio(
            captured, total_suitability
        )
        output[f"{name}_mean_future_mpei_kg_m3"] = ratio(
            exposure_numerator, supported_captured
        )
        output[f"{name}_mpei_supported_selected_area_fraction"] = ratio(
            supported_selected_area, selected_area
        )
        output[f"{name}_mpei_supported_captured_co_suitability_fraction"] = ratio(
            supported_captured, captured
        )
    overlap = float(np.sum(area * np.minimum(climate, exposure), dtype=np.float64))
    output["future_mpei_vs_climate_co_suitability_retained_fraction"] = ratio(
        output["future_mpei_priority_captured_co_suitability_fraction"],
        output["climate_priority_captured_co_suitability_fraction"],
    )
    output["future_mpei_vs_climate_mean_mpei_change_pct"] = 100.0 * ratio(
        output["future_mpei_priority_mean_future_mpei_kg_m3"]
        - output["climate_priority_mean_future_mpei_kg_m3"],
        output["climate_priority_mean_future_mpei_kg_m3"],
    )
    output["future_mpei_vs_existing_mean_mpei_change_pct"] = 100.0 * ratio(
        output["future_mpei_priority_mean_future_mpei_kg_m3"]
        - output["existing_mpa_mean_future_mpei_kg_m3"],
        output["existing_mpa_mean_future_mpei_kg_m3"],
    )
    output["future_mpei_vs_climate_area_changed_fraction"] = 1.0 - ratio(
        overlap, budget
    )
    output["maximum_area_budget_closure_error_km2"] = float(
        max(abs(output[f"{name}_selected_area_km2"] - budget) for name in selections)
    )
    return output, selections


def region_mask(reef: pd.DataFrame, bounds) -> np.ndarray:
    if bounds is None:
        return np.ones(len(reef), dtype=bool)
    west, east, south, north = bounds
    return (
        (reef["longitude"].to_numpy() >= west)
        & (reef["longitude"].to_numpy() <= east)
        & (reef["latitude"].to_numpy() >= south)
        & (reef["latitude"].to_numpy() <= north)
    )


def point_estimates(
    reef: pd.DataFrame,
    component_columns: list[str],
    supported_domain_only: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    records = []
    components = []
    area_all = reef["reef_area_km2"].to_numpy(dtype=float)
    mpa_all = np.clip(
        reef["mpa_fraction_with_point_buffers"].to_numpy(dtype=float), 0.0, 1.0
    )
    mpei_all = reef["mpei_total"].to_numpy(dtype=float)
    exposure_all = reef["mpei_e_star_global_reef_minmax"].to_numpy(dtype=float)
    mpei_supported_all = reef["mpei_supported"].to_numpy(dtype=bool)
    for region_name, bounds in REGIONS.items():
        region = region_mask(reef, bounds)
        for scenario in SCENARIOS:
            support = reef[f"P_{scenario}"].to_numpy(dtype=float)
            valid = region & reef[f"P_valid_{scenario}"].to_numpy(dtype=bool)
            candidate = valid & (support > 0.0)
            if supported_domain_only:
                candidate &= mpei_supported_all
            metrics, selections = calculate_metrics(
                area_all[candidate],
                mpa_all[candidate],
                support[candidate],
                exposure_all[candidate],
                mpei_all[candidate],
                mpei_supported_all[candidate],
            )
            records.append({"region": region_name, "scenario": scenario, **metrics})
            for scheme, selection in selections.items():
                denominator = np.sum(
                    area_all[candidate]
                    * support[candidate]
                    * selection
                    * np.where(mpei_supported_all[candidate], mpei_all[candidate], 0.0),
                    dtype=np.float64,
                )
                for column in component_columns:
                    numerator = np.sum(
                        area_all[candidate]
                        * support[candidate]
                        * selection
                        * np.where(
                            mpei_supported_all[candidate],
                            reef[column].to_numpy(dtype=float)[candidate],
                            0.0,
                        ),
                        dtype=np.float64,
                    )
                    components.append(
                        {
                            "region": region_name,
                            "scenario": scenario,
                            "scheme": scheme,
                            "component_type": column.split("_")[1],
                            "component": "_".join(column.split("_")[2:]),
                            "share_of_selected_future_mpei_percent": 100.0
                            * ratio(float(numerator), float(denominator)),
                        }
                    )
    return pd.DataFrame(records), pd.DataFrame(components)


def bootstrap(reef: pd.DataFrame, member_masks: dict[str, np.ndarray]) -> pd.DataFrame:
    area_all = reef["reef_area_km2"].to_numpy(dtype=float)
    mpa_all = np.clip(
        reef["mpa_fraction_with_point_buffers"].to_numpy(dtype=float), 0.0, 1.0
    )
    mpei_all = reef["mpei_total"].to_numpy(dtype=float)
    exposure_all = reef["mpei_e_star_global_reef_minmax"].to_numpy(dtype=float)
    mpei_supported_all = reef["mpei_supported"].to_numpy(dtype=bool)
    rng = np.random.default_rng(SEED)
    region_cache = {}
    for region_name, bounds in REGIONS.items():
        indices = np.flatnonzero(region_mask(reef, bounds))
        longitude = reef.loc[indices, "longitude"].to_numpy(dtype=float)
        latitude = reef.loc[indices, "latitude"].to_numpy(dtype=float)
        block_key = (
            np.floor((latitude + 90.0) / SPATIAL_BLOCK_DEG).astype(int) * 1000
            + np.floor((longitude + 180.0) / SPATIAL_BLOCK_DEG).astype(int)
        )
        _, inverse = np.unique(block_key, return_inverse=True)
        region_cache[region_name] = {
            "indices": indices,
            "inverse": inverse,
            "n_blocks": int(inverse.max() + 1),
        }
    records = []
    for replicate in range(N_BOOTSTRAP):
        member_draw = rng.integers(0, len(MEMBERS), len(MEMBERS))
        for region_name, cache in region_cache.items():
            sampled_blocks = rng.integers(
                0, cache["n_blocks"], cache["n_blocks"]
            )
            counts = np.bincount(
                sampled_blocks, minlength=cache["n_blocks"]
            ).astype(float)
            multiplicity = counts[cache["inverse"]]
            indices = cache["indices"]
            for scenario in SCENARIOS:
                matrix = member_masks[scenario][:, indices]
                support = np.mean(matrix[member_draw], axis=0)
                valid = np.all(np.isfinite(matrix), axis=0)
                candidate = valid & (support > 0.0)
                metrics, _ = calculate_metrics(
                    area_all[indices][candidate] * multiplicity[candidate],
                    mpa_all[indices][candidate],
                    support[candidate],
                    exposure_all[indices][candidate],
                    mpei_all[indices][candidate],
                    mpei_supported_all[indices][candidate],
                )
                records.append(
                    {
                        "replicate": replicate,
                        "region": region_name,
                        "scenario": scenario,
                        **metrics,
                    }
                )
        if (replicate + 1) % 100 == 0:
            print(f"bootstrap {replicate + 1}/{N_BOOTSTRAP}", flush=True)
    return pd.DataFrame(records)


def summarize(point: pd.DataFrame, draws: pd.DataFrame) -> pd.DataFrame:
    metrics = [column for column in point if column not in {"region", "scenario"}]
    records = []
    for row in point.itertuples(index=False):
        subset = draws.loc[
            (draws["region"] == row.region) & (draws["scenario"] == row.scenario)
        ]
        for metric in metrics:
            values = pd.to_numeric(subset[metric], errors="coerce").dropna().to_numpy()
            lower, upper = np.quantile(values, [0.025, 0.975])
            records.append(
                {
                    "region": row.region,
                    "scenario": row.scenario,
                    "metric": metric,
                    "estimate": float(getattr(row, metric)),
                    "lower_95": float(lower),
                    "upper_95": float(upper),
                    "n_bootstrap": int(values.size),
                }
            )
    return pd.DataFrame(records)


def main() -> None:
    # The retained command now runs the manuscript-matched repaired habitat line.
    from run_repaired_analysis_20260906 import main as run_repaired
    run_repaired()


if __name__ == "__main__":
    main()
