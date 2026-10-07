"""Scientific calculation functions copied unchanged from the accepted September 6 workflow.

Only source IO and plotting are outside this module; see provenance.json for origins.
"""

from __future__ import annotations

import numpy as np

import pandas as pd

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
