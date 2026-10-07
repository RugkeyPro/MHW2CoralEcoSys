#!/usr/bin/env python3
"""Recalculate Figure 4 risk with period-matched 36-tracer MPEI fields.

Panels a-d retain the validated fixed-historical-domain habitat estimand. The
exposure analysis compares 1993-2022 historical habitat plus the mean
1993-2022 MPEI field with 2045-2055 future habitat plus the mean 2045-2055
MPEI field. Future exposure uses the four-taxon OR10 co-suitable domain and no
MESS screen, matching the current Figure 4 exposure definition.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator


PROJECT_ROOT = Path(__file__).resolve().parents[4]
AUTHORITY_PROJECT = Path(r"D:\ecoMarine")
PACKAGE = Path(__file__).resolve().parents[1]
SOURCE = PACKAGE / "source"
QA = PACKAGE / "qa"

HISTORICAL_ROOT = (
    AUTHORITY_PROJECT
    / "output"
    / "01_mainline"
    / "inputs"
    / "sdm_global_0083_corrected_20260722"
)
HISTORICAL_HSI_ROOT = HISTORICAL_ROOT / "maxent" / "temporal"
OR10_TABLE = (
    HISTORICAL_ROOT
    / "tables"
    / "or10_constrained_temporal_mean_1993_2022_20260724.csv"
)
REPAIR_ROOT = (
    AUTHORITY_PROJECT
    / "output"
    / "02_experiments"
    / "20260816_future_coverage_repair"
)
FUTURE_HSI_ROOT = REPAIR_ROOT / "derived" / "hsi_member_projections"
QUALITY_REFERENCE_TABLE = (
    REPAIR_ROOT
    / "results"
    / "continuous_cosuitability_capacity_20260820"
    / "tables"
    / "continuous_cosuitability_capacity_by_member_20260820.csv"
)
MPEI_PACKAGE = (
    PROJECT_ROOT
    / "output"
    / "01_mainline"
    / "inputs"
    / "mpei_36tracer_annual_1950_2060_20260901"
)
HISTORICAL_MPEI_FILE = (
    MPEI_PACKAGE
    / "derived"
    / "mpei_36tracer_0_100m_mean_1993_2022_20260901.nc"
)
FUTURE_MPEI_FILE = (
    MPEI_PACKAGE
    / "derived"
    / "mpei_36tracer_0_100m_mean_2045_2055_20260901.nc"
)

SPECIES = ("Acropora", "Lobophora", "Scarus", "Cephalopholis")
MEMBERS = (
    "gfdl_esm4_r1i1p1f1_gr",
    "cnrm_esm2_1_r1i1p1f2_gn",
    "ipsl_cm6a_lr_r1i1p1f1_gn",
)
SCENARIOS = ("ssp126_2050", "ssp245_2050", "ssp585_2050")
REGIONS = {
    "Global": None,
    "Caribbean": (-98.0375, -55.0435, 8.0805, 33.0635),
    "GBR": (141.9985, 156.0255, -24.9535, -9.9305),
    "Southeast_Asia": (89.9575, 159.9265, -10.0135, 25.0125),
}
EARTH_RADIUS_M = 6_371_008.8


def sha256(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def canonical_path(path: Path) -> str:
    """Record repository outputs at their canonical D:/ecoMarine location."""
    resolved = path.resolve()
    try:
        relative = resolved.relative_to(PROJECT_ROOT.resolve())
    except ValueError:
        return str(resolved)
    return str(AUTHORITY_PROJECT / relative)


def read_hsi(path: Path, expected_grid: dict | None = None) -> tuple[np.ndarray, dict, dict]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with rasterio.open(path) as dataset:
        values = dataset.read(1).astype(np.float32)
        if dataset.nodata is not None and np.isfinite(dataset.nodata):
            values[values == dataset.nodata] = np.nan
        values[~np.isfinite(values)] = np.nan
        grid = {
            "shape": dataset.shape,
            "transform": tuple(dataset.transform),
            "crs": dataset.crs.to_string() if dataset.crs else None,
        }
        profile = dataset.profile.copy()
    finite = values[np.isfinite(values)]
    if finite.size == 0 or finite.min() < -1e-6 or finite.max() > 1.0 + 1e-6:
        raise ValueError(f"HSI outside [0, 1] or entirely missing: {path}")
    if expected_grid is not None and grid != expected_grid:
        raise ValueError(f"HSI grid mismatch: {path}")
    return values, grid, profile


def historical_hsi_path(species: str) -> Path:
    return HISTORICAL_HSI_ROOT / f"{species}_global" / "temporal_mean.tif"


def future_hsi_path(member: str, species: str, scenario: str) -> Path:
    return (
        FUTURE_HSI_ROOT
        / member
        / "maxent"
        / "projection"
        / f"{species}_global"
        / f"{scenario}_hsi_mean.tif"
    )


def grid_coordinates(grid: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    height, width = grid["shape"]
    transform = rasterio.Affine(*grid["transform"])
    longitude = transform.c + (np.arange(width) + 0.5) * transform.a
    latitude = transform.f + (np.arange(height) + 0.5) * transform.e
    latitude_edges = transform.f + np.arange(height + 1) * transform.e
    row_area_km2 = EARTH_RADIUS_M**2 * abs(transform.a) * np.pi / 180.0
    row_area_km2 *= np.abs(
        np.sin(np.deg2rad(latitude_edges[1:]))
        - np.sin(np.deg2rad(latitude_edges[:-1]))
    )
    return longitude, latitude, row_area_km2 / 1e6


def regional_indices(
    longitude: np.ndarray, latitude: np.ndarray
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    output: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for name, bounds in REGIONS.items():
        if bounds is None:
            output[name] = (np.arange(latitude.size), np.arange(longitude.size))
            continue
        west, east, south, north = bounds
        rows = np.flatnonzero((latitude >= south) & (latitude <= north))
        columns = np.flatnonzero((longitude >= west) & (longitude <= east))
        if rows.size == 0 or columns.size == 0:
            raise ValueError(f"No HSI cells in {name}")
        output[name] = (rows, columns)
    return output


def nearest(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    position = np.searchsorted(source, target)
    position = np.clip(position, 1, len(source) - 1)
    left = position - 1
    return np.where(
        np.abs(target - source[left]) <= np.abs(target - source[position]),
        left,
        position,
    )


def read_mpei(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with Dataset(path, "r") as dataset:
        latitude = np.asarray(dataset.variables["lat"][:], dtype=float)
        longitude = np.asarray(dataset.variables["lon"][:], dtype=float)
        field = np.asarray(
            np.ma.filled(dataset.variables["mpei_0_100m"][:], np.nan),
            dtype=np.float64,
        )
        metadata = {
            "path": canonical_path(path),
            "sha256": sha256(path),
            "source_sha256": str(dataset.source_sha256),
            "completed_source_year_start": int(dataset.completed_source_year_start),
            "completed_source_year_end": int(dataset.completed_source_year_end),
            "annual_endpoint_count": int(dataset.annual_endpoint_count),
            "temporal_aggregation": str(dataset.temporal_aggregation),
            "physical_forcing_note": str(dataset.physical_forcing_note),
        }
    if not np.all(np.diff(latitude) > 0) or not np.all(np.diff(longitude) > 0):
        raise ValueError(f"MPEI coordinates must be strictly ascending: {path}")
    finite = field[np.isfinite(field)]
    if finite.size == 0 or finite.min() < 0:
        raise ValueError(f"MPEI is empty or contains negative values: {path}")
    return latitude, longitude, field, metadata


def map_nearest(
    source_latitude: np.ndarray,
    source_longitude: np.ndarray,
    source_field: np.ndarray,
    target_latitude: np.ndarray,
    target_longitude: np.ndarray,
) -> tuple[np.ndarray, dict]:
    cyclic_longitude = np.mod(target_longitude, 360.0)
    source_rows = nearest(source_latitude, target_latitude)
    source_columns = nearest(source_longitude, cyclic_longitude)
    mapped = source_field[source_rows[:, None], source_columns[None, :]].astype(
        np.float32
    )
    diagnostics = {
        "maximum_latitude_centre_distance_degrees": float(
            np.max(np.abs(target_latitude - source_latitude[source_rows]))
        ),
        "maximum_longitude_centre_distance_degrees": float(
            np.max(np.abs(cyclic_longitude - source_longitude[source_columns]))
        ),
        "finite_target_cells": int(np.isfinite(mapped).sum()),
        "target_cells": int(mapped.size),
    }
    return mapped, diagnostics


def map_bilinear(
    source_latitude: np.ndarray,
    source_longitude: np.ndarray,
    source_field: np.ndarray,
    target_latitude: np.ndarray,
    target_longitude: np.ndarray,
) -> np.ndarray:
    valid = np.isfinite(source_field)
    extended_longitude = np.concatenate(
        ([source_longitude[-1] - 360.0], source_longitude, [source_longitude[0] + 360.0])
    )
    extended_field = np.concatenate(
        [source_field[:, -1:], source_field, source_field[:, :1]], axis=1
    )
    extended_valid = np.concatenate(
        [valid[:, -1:], valid, valid[:, :1]], axis=1
    ).astype(float)
    numerator = RegularGridInterpolator(
        (source_latitude, extended_longitude),
        np.where(np.isfinite(extended_field), extended_field, 0.0),
        method="linear",
        bounds_error=False,
        fill_value=0.0,
    )
    denominator = RegularGridInterpolator(
        (source_latitude, extended_longitude),
        extended_valid,
        method="linear",
        bounds_error=False,
        fill_value=0.0,
    )
    output = np.full((target_latitude.size, target_longitude.size), np.nan, dtype=np.float32)
    cyclic_longitude = np.mod(target_longitude, 360.0)
    for start in range(0, target_latitude.size, 100):
        stop = min(start + 100, target_latitude.size)
        lon2, lat2 = np.meshgrid(cyclic_longitude, target_latitude[start:stop])
        points = np.column_stack([lat2.ravel(), lon2.ravel()])
        numerator_values = numerator(points)
        denominator_values = denominator(points)
        values = np.full(points.shape[0], np.nan, dtype=np.float64)
        np.divide(
            numerator_values,
            denominator_values,
            out=values,
            where=denominator_values > 1e-12,
        )
        output[start:stop] = values.reshape(stop - start, target_longitude.size)
    return output


def joint_hsi(stack: np.ndarray, domain: np.ndarray) -> np.ndarray:
    output = np.zeros(domain.shape, dtype=np.float32)
    output[domain] = np.prod(stack[:, domain], axis=0, dtype=np.float64) ** 0.25
    return output


def weighted_sum(
    values: np.ndarray,
    mask: np.ndarray,
    rows: np.ndarray,
    columns: np.ndarray,
    row_area_km2: np.ndarray,
) -> float:
    subset = np.where(mask[np.ix_(rows, columns)], values[np.ix_(rows, columns)], 0.0)
    return float(np.dot(subset.sum(axis=1, dtype=np.float64), row_area_km2[rows]))


def exposure_metrics(
    habitat_weight: np.ndarray,
    habitat_domain: np.ndarray,
    mpei: np.ndarray,
    rows: np.ndarray,
    columns: np.ndarray,
    row_area_km2: np.ndarray,
) -> dict[str, float | int]:
    valid = habitat_domain & np.isfinite(mpei) & (mpei >= 0.0)
    weighted_area = weighted_sum(habitat_weight, valid, rows, columns, row_area_km2)
    burden = weighted_sum(habitat_weight * mpei, valid, rows, columns, row_area_km2)
    area = weighted_sum(np.ones_like(habitat_weight), valid, rows, columns, row_area_km2)
    cell_count = int(valid[np.ix_(rows, columns)].sum())
    if weighted_area <= 0.0 or burden < 0.0 or cell_count == 0:
        raise ValueError("Invalid habitat-weighted exposure state")
    return {
        "habitat_area_km2": area,
        "hsi_weighted_area_km2": weighted_area,
        "mpei_burden_kg_m3_km2": burden,
        "mpei_pressure_kg_m3": burden / weighted_area,
        "n_habitat_cells_with_mpei": cell_count,
    }


def write_raster(path: Path, values: np.ndarray, profile: dict) -> None:
    output = np.where(np.isfinite(values), values, -9999.0).astype(np.float32)
    output_profile = profile.copy()
    output_profile.update(
        driver="GTiff",
        count=1,
        dtype="float32",
        nodata=-9999.0,
        compress="LZW",
    )
    with rasterio.open(path, "w", **output_profile) as dataset:
        dataset.write(output, 1)


def percentage_change(future: float, historical: float) -> float:
    if historical <= 0.0:
        raise ValueError("Historical reference must be positive")
    return 100.0 * (future - historical) / historical


def main() -> None:
    SOURCE.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)

    thresholds = (
        pd.read_csv(OR10_TABLE)
        .set_index("species")
        .loc[list(SPECIES), "or10"]
        .to_numpy(dtype=float)
    )
    historical_layers: list[np.ndarray] = []
    grid: dict | None = None
    raster_profile: dict | None = None
    manifest: list[dict[str, str]] = [
        {
            "kind": "fixed_global_or10",
            "path": canonical_path(OR10_TABLE),
            "sha256": sha256(OR10_TABLE),
        }
    ]
    for species in SPECIES:
        path = historical_hsi_path(species)
        layer, candidate_grid, candidate_profile = read_hsi(path, grid)
        if grid is None:
            grid, raster_profile = candidate_grid, candidate_profile
        historical_layers.append(layer)
        manifest.append(
            {
                "kind": "historical_hsi",
                "species": species,
                "path": canonical_path(path),
                "sha256": sha256(path),
            }
        )
    assert grid is not None and raster_profile is not None
    historical_stack = np.stack(historical_layers)
    historical_finite = np.all(np.isfinite(historical_stack), axis=0)
    historical_suitable = historical_finite & np.all(
        historical_stack >= thresholds[:, None, None], axis=0
    )
    historical_joint = joint_hsi(historical_stack, historical_suitable)
    longitude, latitude, row_area_km2 = grid_coordinates(grid)
    regions = regional_indices(longitude, latitude)

    hist_lat, hist_lon, historical_native_mpei, historical_mpei_metadata = read_mpei(
        HISTORICAL_MPEI_FILE
    )
    fut_lat, fut_lon, future_native_mpei, future_mpei_metadata = read_mpei(FUTURE_MPEI_FILE)
    if not np.array_equal(hist_lat, fut_lat) or not np.array_equal(hist_lon, fut_lon):
        raise ValueError("Historical and future MPEI coordinates differ")
    if not np.array_equal(np.isfinite(historical_native_mpei), np.isfinite(future_native_mpei)):
        raise ValueError("Historical and future MPEI finite masks differ")
    historical_mpei, nearest_diagnostics = map_nearest(
        hist_lat, hist_lon, historical_native_mpei, latitude, longitude
    )
    future_mpei, future_nearest_diagnostics = map_nearest(
        fut_lat, fut_lon, future_native_mpei, latitude, longitude
    )
    historical_mpei_bilinear = map_bilinear(
        hist_lat, hist_lon, historical_native_mpei, latitude, longitude
    )
    future_mpei_bilinear = map_bilinear(
        fut_lat, fut_lon, future_native_mpei, latitude, longitude
    )
    if not np.array_equal(np.isfinite(historical_mpei), np.isfinite(future_mpei)):
        raise ValueError("Nearest-mapped historical and future MPEI finite masks differ")

    historical_states: dict[str, dict[str, float | int]] = {}
    historical_states_bilinear: dict[str, dict[str, float | int]] = {}
    for region, (rows, columns) in regions.items():
        historical_states[region] = exposure_metrics(
            historical_joint,
            historical_suitable,
            historical_mpei,
            rows,
            columns,
            row_area_km2,
        )
        historical_states_bilinear[region] = exposure_metrics(
            historical_joint,
            historical_suitable,
            historical_mpei_bilinear,
            rows,
            columns,
            row_area_km2,
        )

    quality_rows: list[dict[str, float | str]] = []
    exposure_rows: list[dict[str, float | str | int]] = []
    state_rows: list[dict[str, float | str | int]] = []
    factorial_rows: list[dict[str, float | str]] = []
    mapping_sensitivity_rows: list[dict[str, float | str | bool]] = []
    map_statistics: list[dict[str, float | str | int]] = []

    for scenario in SCENARIOS:
        future_stacks: dict[str, np.ndarray] = {}
        common_finite_masks = [historical_finite]
        for member in MEMBERS:
            layers: list[np.ndarray] = []
            for species in SPECIES:
                path = future_hsi_path(member, species, scenario)
                layer, _, _ = read_hsi(path, grid)
                layers.append(layer)
                manifest.append(
                    {
                        "kind": "future_hsi",
                        "member": member,
                        "scenario": scenario,
                        "species": species,
                        "path": canonical_path(path),
                        "sha256": sha256(path),
                    }
                )
            stack = np.stack(layers)
            future_stacks[member] = stack
            common_finite_masks.append(np.all(np.isfinite(stack), axis=0))

        common_finite_support = np.all(np.stack(common_finite_masks), axis=0)
        quality_domain = historical_suitable & common_finite_support
        if not quality_domain.any():
            raise ValueError(f"Empty fixed historical quality domain: {scenario}")
        historical_quality_joint = joint_hsi(historical_stack, quality_domain)
        ensemble_future_quality = np.mean(
            [joint_hsi(stack, quality_domain) for stack in future_stacks.values()],
            axis=0,
        )
        quality_reduction = np.full(historical_quality_joint.shape, np.nan, dtype=np.float32)
        quality_reduction[quality_domain] = 100.0 * (
            1.0
            - ensemble_future_quality[quality_domain]
            / historical_quality_joint[quality_domain]
        )
        map_path = SOURCE / f"{scenario}_ensemble_mean_quality_reduction_pct_20260902.tif"
        write_raster(map_path, quality_reduction, raster_profile)
        manifest.append(
            {
                "kind": "derived_quality_reduction_map",
                "scenario": scenario,
                "path": canonical_path(map_path),
                "sha256": sha256(map_path),
            }
        )
        finite_reduction = quality_reduction[np.isfinite(quality_reduction)]
        map_statistics.append(
            {
                "scenario": scenario,
                "finite_cells": int(finite_reduction.size),
                "quality_reduction_min_pct": float(finite_reduction.min()),
                "quality_reduction_p01_pct": float(np.quantile(finite_reduction, 0.01)),
                "quality_reduction_p99_pct": float(np.quantile(finite_reduction, 0.99)),
                "quality_reduction_max_pct": float(finite_reduction.max()),
            }
        )

        for member, future_stack in future_stacks.items():
            future_quality_joint = joint_hsi(future_stack, quality_domain)
            future_finite = np.all(np.isfinite(future_stack), axis=0)
            future_suitable = future_finite & np.all(
                future_stack >= thresholds[:, None, None], axis=0
            )
            future_joint = joint_hsi(future_stack, future_suitable)
            for region, (rows, columns) in regions.items():
                historical_quality = weighted_sum(
                    historical_quality_joint,
                    quality_domain,
                    rows,
                    columns,
                    row_area_km2,
                )
                future_quality = weighted_sum(
                    future_quality_joint,
                    quality_domain,
                    rows,
                    columns,
                    row_area_km2,
                )
                quality_rows.append(
                    {
                        "member": member,
                        "scenario": scenario,
                        "region": region,
                        "historical_continuous_quality_km2_hsi": historical_quality,
                        "future_continuous_quality_km2_hsi": future_quality,
                        "continuous_quality_reduction_pct": 100.0
                        * (historical_quality - future_quality)
                        / historical_quality,
                        "continuous_quality_retention_pct": 100.0
                        * future_quality
                        / historical_quality,
                    }
                )

                configurations = {
                    "historical_habitat_historical_mpei": (
                        historical_joint,
                        historical_suitable,
                        historical_mpei,
                    ),
                    "historical_habitat_future_mpei": (
                        historical_joint,
                        historical_suitable,
                        future_mpei,
                    ),
                    "future_habitat_historical_mpei": (
                        future_joint,
                        future_suitable,
                        historical_mpei,
                    ),
                    "future_habitat_future_mpei": (
                        future_joint,
                        future_suitable,
                        future_mpei,
                    ),
                }
                states: dict[str, dict[str, float | int]] = {}
                for configuration, (weight, domain, field) in configurations.items():
                    metrics = exposure_metrics(
                        weight, domain, field, rows, columns, row_area_km2
                    )
                    states[configuration] = metrics
                    state_rows.append(
                        {
                            "member": member,
                            "scenario": scenario,
                            "region": region,
                            "configuration": configuration,
                            "habitat_period": (
                                "1993-2022"
                                if configuration.startswith("historical_habitat")
                                else "2045-2055"
                            ),
                            "mpei_period": (
                                "1993-2022"
                                if configuration.endswith("historical_mpei")
                                else "2045-2055"
                            ),
                            **metrics,
                        }
                    )

                baseline = states["historical_habitat_historical_mpei"]
                mpei_only = states["historical_habitat_future_mpei"]
                habitat_only = states["future_habitat_historical_mpei"]
                joint = states["future_habitat_future_mpei"]
                exposure_record: dict[str, float | str | int] = {
                    "member": member,
                    "scenario": scenario,
                    "region": region,
                    "historical_mpei_burden_kg_m3_km2": baseline[
                        "mpei_burden_kg_m3_km2"
                    ],
                    "future_mpei_burden_kg_m3_km2": joint[
                        "mpei_burden_kg_m3_km2"
                    ],
                    "historical_mpei_pressure_kg_m3": baseline["mpei_pressure_kg_m3"],
                    "future_mpei_pressure_kg_m3": joint["mpei_pressure_kg_m3"],
                    "dynamic_mpei_burden_change_pct": percentage_change(
                        float(joint["mpei_burden_kg_m3_km2"]),
                        float(baseline["mpei_burden_kg_m3_km2"]),
                    ),
                    "dynamic_mpei_pressure_change_pct": percentage_change(
                        float(joint["mpei_pressure_kg_m3"]),
                        float(baseline["mpei_pressure_kg_m3"]),
                    ),
                    "future_habitat_area_km2": joint["habitat_area_km2"],
                    "future_hsi_weighted_area_km2": joint["hsi_weighted_area_km2"],
                }
                exposure_rows.append(exposure_record)

                factorial_record: dict[str, float | str] = {
                    "member": member,
                    "scenario": scenario,
                    "region": region,
                }
                for metric in ("mpei_burden_kg_m3_km2", "mpei_pressure_kg_m3"):
                    base_value = float(baseline[metric])
                    habitat_value = float(habitat_only[metric])
                    mpei_value = float(mpei_only[metric])
                    joint_value = float(joint[metric])
                    habitat_effect = habitat_value - base_value
                    mpei_effect = mpei_value - base_value
                    interaction = joint_value - habitat_value - mpei_value + base_value
                    factorial_record.update(
                        {
                            f"baseline_{metric}": base_value,
                            f"habitat_only_{metric}": habitat_value,
                            f"mpei_only_{metric}": mpei_value,
                            f"joint_{metric}": joint_value,
                            f"joint_change_pct_{metric}": percentage_change(
                                joint_value, base_value
                            ),
                            f"habitat_contribution_pct_of_baseline_{metric}": 100.0
                            * habitat_effect
                            / base_value,
                            f"mpei_contribution_pct_of_baseline_{metric}": 100.0
                            * mpei_effect
                            / base_value,
                            f"interaction_contribution_pct_of_baseline_{metric}": 100.0
                            * interaction
                            / base_value,
                            f"factorial_closure_error_{metric}": joint_value
                            - base_value
                            - habitat_effect
                            - mpei_effect
                            - interaction,
                        }
                    )
                factorial_rows.append(factorial_record)

                future_bilinear = exposure_metrics(
                    future_joint,
                    future_suitable,
                    future_mpei_bilinear,
                    rows,
                    columns,
                    row_area_km2,
                )
                baseline_bilinear = historical_states_bilinear[region]
                for metric, output_name in (
                    ("mpei_burden_kg_m3_km2", "dynamic_mpei_burden_change_pct"),
                    ("mpei_pressure_kg_m3", "dynamic_mpei_pressure_change_pct"),
                ):
                    nearest_change = float(exposure_record[output_name])
                    bilinear_change = percentage_change(
                        float(future_bilinear[metric]), float(baseline_bilinear[metric])
                    )
                    mapping_sensitivity_rows.append(
                        {
                            "member": member,
                            "scenario": scenario,
                            "region": region,
                            "metric": output_name,
                            "nearest_change_pct": nearest_change,
                            "bilinear_change_pct": bilinear_change,
                            "bilinear_minus_nearest_percentage_points": bilinear_change
                            - nearest_change,
                            "change_sign_matches": bool(
                                np.sign(bilinear_change) == np.sign(nearest_change)
                            ),
                        }
                    )

    quality_by_member = pd.DataFrame(quality_rows)
    exposure_by_member = pd.DataFrame(exposure_rows)
    states = pd.DataFrame(state_rows)
    factorial = pd.DataFrame(factorial_rows)
    mapping_sensitivity = pd.DataFrame(mapping_sensitivity_rows)

    quality_measures = [
        column
        for column in quality_by_member.columns
        if column not in {"member", "scenario", "region"}
    ]
    quality_summary = quality_by_member.groupby(
        ["scenario", "region"], sort=False
    )[quality_measures].agg(["mean", "min", "max"])
    quality_summary.columns = [
        f"{metric}_{stat}" for metric, stat in quality_summary.columns
    ]
    quality_summary = quality_summary.reset_index()

    exposure_measures = [
        "dynamic_mpei_burden_change_pct",
        "dynamic_mpei_pressure_change_pct",
        "historical_mpei_burden_kg_m3_km2",
        "future_mpei_burden_kg_m3_km2",
        "historical_mpei_pressure_kg_m3",
        "future_mpei_pressure_kg_m3",
        "future_habitat_area_km2",
        "future_hsi_weighted_area_km2",
    ]
    exposure_summary = exposure_by_member.groupby(
        ["scenario", "region"], sort=False
    )[exposure_measures].agg(["mean", "min", "max"])
    exposure_summary.columns = [
        f"{metric}_{stat}" for metric, stat in exposure_summary.columns
    ]
    exposure_summary = exposure_summary.reset_index()

    factorial_measures = [
        column
        for column in factorial.columns
        if column not in {"member", "scenario", "region"}
    ]
    factorial_summary = factorial.groupby(["scenario", "region"], sort=False)[
        factorial_measures
    ].agg(["mean", "min", "max"])
    factorial_summary.columns = [
        f"{metric}_{stat}" for metric, stat in factorial_summary.columns
    ]
    factorial_summary = factorial_summary.reset_index()

    output_tables = {
        "quality_by_member": SOURCE
        / "figure4_quality_reduction_by_member_20260902.csv",
        "quality_summary": SOURCE
        / "figure4_quality_reduction_ensemble_20260902.csv",
        "exposure_states": SOURCE
        / "figure4_period_matched_mpei_states_by_member_20260902.csv",
        "exposure_by_member": SOURCE
        / "figure4_period_matched_mpei_change_by_member_20260902.csv",
        "exposure_summary": SOURCE
        / "figure4_period_matched_mpei_change_ensemble_20260902.csv",
        "factorial_by_member": SOURCE
        / "figure4_period_matched_mpei_factorial_by_member_20260902.csv",
        "factorial_summary": SOURCE
        / "figure4_period_matched_mpei_factorial_ensemble_20260902.csv",
        "mapping_sensitivity": SOURCE
        / "figure4_period_matched_mpei_mapping_sensitivity_20260902.csv",
        "map_statistics": SOURCE
        / "figure4_quality_reduction_map_statistics_20260902.csv",
        "input_manifest": SOURCE / "figure4_input_manifest_20260902.csv",
    }
    quality_by_member.to_csv(output_tables["quality_by_member"], index=False)
    quality_summary.to_csv(output_tables["quality_summary"], index=False)
    states.to_csv(output_tables["exposure_states"], index=False)
    exposure_by_member.to_csv(output_tables["exposure_by_member"], index=False)
    exposure_summary.to_csv(output_tables["exposure_summary"], index=False)
    factorial.to_csv(output_tables["factorial_by_member"], index=False)
    factorial_summary.to_csv(output_tables["factorial_summary"], index=False)
    mapping_sensitivity.to_csv(output_tables["mapping_sensitivity"], index=False)
    pd.DataFrame(map_statistics).to_csv(output_tables["map_statistics"], index=False)
    manifest.extend(
        [
            {
                "kind": "historical_mpei_1993_2022_mean",
                "path": canonical_path(HISTORICAL_MPEI_FILE),
                "sha256": historical_mpei_metadata["sha256"],
            },
            {
                "kind": "future_mpei_2045_2055_mean",
                "path": canonical_path(FUTURE_MPEI_FILE),
                "sha256": future_mpei_metadata["sha256"],
            },
        ]
    )
    pd.DataFrame(manifest).to_csv(output_tables["input_manifest"], index=False)

    reference = pd.read_csv(QUALITY_REFERENCE_TABLE)
    quality_audit = quality_by_member.merge(
        reference[
            [
                "member",
                "scenario",
                "region",
                "historical_or10_domain_capacity_retention_pct",
            ]
        ],
        on=["member", "scenario", "region"],
        validate="one_to_one",
    )
    maximum_quality_reference_difference = float(
        np.max(
            np.abs(
                quality_audit["continuous_quality_reduction_pct"].to_numpy(float)
                - (
                    100.0
                    - quality_audit[
                        "historical_or10_domain_capacity_retention_pct"
                    ].to_numpy(float)
                )
            )
        )
    )
    closure_columns = [
        column for column in factorial.columns if column.startswith("factorial_closure_error_")
    ]
    maximum_factorial_closure = float(
        max(factorial[column].abs().max() for column in closure_columns)
    )
    maximum_factorial_scale = float(
        max(
            factorial[column].abs().max()
            for column in factorial.columns
            if column.startswith("joint_") and not column.startswith("joint_change")
        )
    )
    factorial_tolerance = float(
        64.0 * np.finfo(np.float64).eps * maximum_factorial_scale
    )
    expected_member_rows = len(MEMBERS) * len(SCENARIOS) * len(REGIONS)
    expected_summary_rows = len(SCENARIOS) * len(REGIONS)
    numerical_frames = [quality_by_member, exposure_by_member, states, factorial]
    nonfinite_counts = {
        str(index): int(
            np.size(frame.select_dtypes(include=[np.number]).to_numpy())
            - np.isfinite(frame.select_dtypes(include=[np.number]).to_numpy()).sum()
        )
        for index, frame in enumerate(numerical_frames)
    }
    status = "PASS"
    if (
        len(quality_by_member) != expected_member_rows
        or len(exposure_by_member) != expected_member_rows
        or len(quality_summary) != expected_summary_rows
        or len(exposure_summary) != expected_summary_rows
        or maximum_quality_reference_difference > 1e-6
        or maximum_factorial_closure > factorial_tolerance
        or any(nonfinite_counts.values())
    ):
        status = "FAIL"

    qc = {
        "status": status,
        "analysis_definition": (
            "Panels a-d reproduce the fixed 1993-2022 historical OR10-domain "
            "continuous-quality analysis with repaired 2045-2055 HSI. Panels e-f "
            "compare historical habitat with mean 1993-2022 MPEI against future "
            "four-taxon OR10 habitat with mean 2045-2055 MPEI; no MESS screen."
        ),
        "future_mpei_scenario_rule": (
            "One common 2045-2055 MPEI field is paired with all three SSP habitat "
            "projections; the MPEI source does not provide SSP-specific trajectories."
        ),
        "mpei_field_semantics": (
            "Arithmetic means of annual instantaneous endpoints after each completed "
            "source year; not within-year concentration means."
        ),
        "historical_mpei_metadata": historical_mpei_metadata,
        "future_mpei_metadata": future_mpei_metadata,
        "primary_grid_mapping": "nearest native IGSM cell centre",
        "nearest_mapping_diagnostics": nearest_diagnostics,
        "future_nearest_mapping_diagnostics": future_nearest_diagnostics,
        "bilinear_sensitivity_rows": int(len(mapping_sensitivity)),
        "bilinear_sensitivity_sign_match_fraction": float(
            mapping_sensitivity["change_sign_matches"].mean()
        ),
        "maximum_absolute_bilinear_minus_nearest_percentage_points": float(
            mapping_sensitivity[
                "bilinear_minus_nearest_percentage_points"
            ].abs().max()
        ),
        "quality_member_rows": int(len(quality_by_member)),
        "exposure_member_rows": int(len(exposure_by_member)),
        "exposure_state_rows": int(len(states)),
        "factorial_rows": int(len(factorial)),
        "summary_rows": int(len(exposure_summary)),
        "maximum_quality_reference_difference_percentage_points": maximum_quality_reference_difference,
        "maximum_factorial_closure_native_units": maximum_factorial_closure,
        "factorial_closure_tolerance_native_units": factorial_tolerance,
        "nonfinite_counts": nonfinite_counts,
        "output_tables": {
            key: canonical_path(value) for key, value in output_tables.items()
        },
    }
    (QA / "figure4_period_matched_mpei_source_qc_20260902.json").write_text(
        json.dumps(qc, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    if status != "PASS":
        raise RuntimeError(json.dumps(qc, ensure_ascii=False, indent=2))
    print(
        json.dumps(
            {
                "status": status,
                "quality_rows": len(quality_by_member),
                "exposure_rows": len(exposure_by_member),
                "mapping_sign_match_fraction": qc[
                    "bilinear_sensitivity_sign_match_fraction"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
