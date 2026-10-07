"""Prepare current-test-branch inputs for Supplementary Fig. S1.

This script reads only the current 0.083-degree mainline inputs. It does not
read archive directories or legacy supplementary branches. Raw records and
rasters are never modified; all products are written below OUTPUT_DIR.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import Resampling, reproject
from sklearn.metrics import roc_auc_score


PROJECT_ROOT = Path(r"D:\ecoMarine")
OUTPUT_DIR = PROJECT_ROOT / "output" / "01_mainline" / "figureS1_sdm_support_20260806"
TABLE_DIR = OUTPUT_DIR / "tables"
RASTER_DIR = OUTPUT_DIR / "rasters"
DOC_DIR = OUTPUT_DIR / "docs"

TAXA = ("Acropora", "Lobophora", "Scarus", "Cephalopholis")
OR10 = {
    "Acropora": 0.3185778,
    "Lobophora": 0.3281061,
    "Scarus": 0.2463830,
    "Cephalopholis": 0.3564393,
}
IUCN_TAXON = {
    "Acropora": "Acropora",
    "Lobophora": "Turbinaria",
    "Scarus": "Scarus",
    "Cephalopholis": "Cephalopholis",
}
IUCN_ROLE = {
    "Acropora": "direct",
    "Lobophora": "proxy_supplementary",
    "Scarus": "direct",
    "Cephalopholis": "direct",
}

MAINLINE = PROJECT_ROOT / "output" / "01_mainline"
SDM_ROOT = MAINLINE / "inputs" / "sdm_global_0083_corrected_20260722"
OCCURRENCE_DIR = MAINLINE / "inputs" / "occurrence_083global_20260727"
IUCN_DIR = PROJECT_ROOT / "validation" / "IUCN" / "rasterized"
IUCN_AUTHORITY = (
    MAINLINE
    / "supplementary_support"
    / "iucn_range_validation_083global_20260804"
    / "tables"
    / "iucn_range_validation_083global_20260804.csv"
)
SURVEY_DIR = MAINLINE / "supplementary_support" / "external_survey_validation_083global_20260804" / "tables"
SURVEY_EVENT_HSI = SURVEY_DIR / "caribbean_survey_event_hsi_0083_complete.csv"
METHODS_FILE = MAINLINE / "methods_083global_20260801" / "supplementary_methods_083global_20260801.md"
PROJECTION_SCRIPT = SDM_ROOT / "scripts" / "branch_05_maxent_projection.R"
ALTERNATIVE_AGREEMENT = (
    MAINLINE
    / "supplementary_support"
    / "alternative_sdm_083global_20260806"
    / "tables"
    / "alternative_sdm_spatial_agreement_083global_20260806.csv"
)


def relative_path(path: Path) -> str:
    """Return a checked project-relative source path."""
    resolved = path.resolve()
    root = PROJECT_ROOT.resolve()
    if root not in resolved.parents and resolved != root:
        raise ValueError(f"Source is outside project root: {resolved}")
    relative = resolved.relative_to(root).as_posix()
    forbidden = ("archive/", "output/02_supplementary/")
    if any(token in relative.lower() for token in forbidden):
        raise ValueError(f"Forbidden non-mainline source: {relative}")
    return relative


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite_values(array: np.ndarray, nodata: float | None) -> tuple[np.ndarray, np.ndarray]:
    values = array.astype("float64", copy=True)
    if nodata is not None and np.isfinite(nodata):
        values[np.isclose(values, nodata)] = np.nan
    return values, np.isfinite(values)


def align_iucn(iucn_path: Path, destination_source: rasterio.DatasetReader) -> np.ndarray:
    """Match IUCN binary ranges to the current 0.083-degree HSI grid."""
    with rasterio.open(iucn_path) as source:
        aligned = np.full((destination_source.height, destination_source.width), np.nan, dtype="float32")
        reproject(
            source=rasterio.band(source, 1),
            destination=aligned,
            src_transform=source.transform,
            src_crs=source.crs,
            dst_transform=destination_source.transform,
            dst_crs=destination_source.crs,
            resampling=Resampling.nearest,
            src_nodata=source.nodata,
            dst_nodata=np.nan,
        )
    return aligned


def binary_metrics(prediction: np.ndarray, observed: np.ndarray) -> dict[str, float | int]:
    predicted = prediction.astype(bool)
    observed_bool = observed.astype(bool)
    tp = int(np.sum(predicted & observed_bool))
    fp = int(np.sum(predicted & ~observed_bool))
    fn = int(np.sum(~predicted & observed_bool))
    tn = int(np.sum(~predicted & ~observed_bool))
    return {
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "TN": tn,
        "specificity": tn / (tn + fp),
        "precision": tp / (tp + fp),
        "sensitivity": tp / (tp + fn),
    }


def prepare_occurrences() -> tuple[pd.DataFrame, pd.DataFrame, list[Path]]:
    records: list[pd.DataFrame] = []
    sources: list[Path] = []
    for taxon in TAXA:
        path = OCCURRENCE_DIR / f"occurrence_thinned_0p083deg_year_month_{taxon}.csv"
        if not path.exists():
            raise FileNotFoundError(path)
        data = pd.read_csv(path)
        required = {"focal_taxon", "longitude", "latitude", "year", "month", "source", "grid_row", "grid_col"}
        if not required.issubset(data.columns):
            raise ValueError(f"Unexpected occurrence schema in {path.name}")
        data = data.loc[:, ["focal_taxon", "longitude", "latitude", "year", "month", "source", "grid_row", "grid_col"]].copy()
        if data["focal_taxon"].nunique() != 1 or data["focal_taxon"].iat[0] != taxon:
            raise ValueError(f"Taxon mismatch in {path.name}")
        if not data["longitude"].between(-180, 180).all() or not data["latitude"].between(-90, 90).all():
            raise ValueError(f"Invalid coordinates in {path.name}")
        records.append(data)
        sources.append(path)
    combined = pd.concat(records, ignore_index=True)
    summary = (
        combined.groupby("focal_taxon", as_index=False)
        .agg(
            displayed_records=("focal_taxon", "size"),
            year_min=("year", "min"),
            year_max=("year", "max"),
            source_count=("source", "nunique"),
        )
        .rename(columns={"focal_taxon": "taxon"})
    )
    occupied_cells = (
        combined.loc[:, ["focal_taxon", "grid_row", "grid_col"]]
        .drop_duplicates()
        .groupby("focal_taxon", as_index=False, sort=False)
        .size()
        .rename(columns={"focal_taxon": "taxon", "size": "occupied_grid_cells"})
    )
    summary = summary.merge(occupied_cells, on="taxon", validate="one_to_one")
    return combined, summary, sources


def prepare_physiology() -> pd.DataFrame:
    """Export the stated physiological multipliers from the active methods and projection equations."""
    temperature = np.round(np.arange(10.0, 40.01, 0.1), 2)
    outer_gamma = 0.5
    rows: list[pd.DataFrame] = []

    omega = np.round(np.arange(1.0, 5.01, 0.05), 2)
    temperature_grid, omega_grid = np.meshgrid(temperature, omega, indexing="ij")
    sigma = np.where(temperature_grid <= 27.0, 5.0, 2.0)
    coral_temperature = np.exp(-((temperature_grid - 27.0) ** 2) / (2.0 * sigma**2))
    coral_omega = 1.0 / (1.0 + np.exp(-4.0 * (omega_grid - 3.0)))
    rows.append(
        pd.DataFrame(
            {
                "taxon": "Acropora",
                "response_geometry": "temperature_by_aragonite_saturation",
                "temperature_c": temperature_grid.ravel(),
                "omega_arag": omega_grid.ravel(),
                "oxygen_mmol_m3": np.nan,
                "physiological_multiplier": (coral_temperature * coral_omega).ravel(),
                "effective_hsi_multiplier": (coral_temperature * coral_omega).ravel() ** outer_gamma,
            }
        )
    )

    kb = 8.617333262e-5
    temperature_kelvin = temperature + 273.15
    activation = np.exp((0.55 / kb) * (1.0 / 293.15 - 1.0 / temperature_kelvin))
    deactivation = 1.0 + np.exp((2.5 / kb) * (1.0 / 308.15 - 1.0 / temperature_kelvin))
    rate = activation / deactivation
    optimum_kelvin = 308.15 * (2.5 - 0.55) / (2.5 - 0.55 + kb * 308.15 * np.log(2.5 / 0.55))
    rate_optimum = (
        np.exp((0.55 / kb) * (1.0 / 293.15 - 1.0 / optimum_kelvin))
        / (1.0 + np.exp((2.5 / kb) * (1.0 / 308.15 - 1.0 / optimum_kelvin)))
    )
    algae_multiplier = np.clip(rate / rate_optimum, 0.0, 1.0)
    rows.append(
        pd.DataFrame(
            {
                "taxon": "Lobophora",
                "response_geometry": "temperature",
                "temperature_c": temperature,
                "omega_arag": np.nan,
                "oxygen_mmol_m3": np.nan,
                "physiological_multiplier": algae_multiplier,
                "effective_hsi_multiplier": algae_multiplier**outer_gamma,
            }
        )
    )

    oxygen = np.arange(0.0, 300.1, 2.0)
    temperature_grid, oxygen_grid = np.meshgrid(temperature, oxygen, indexing="ij")
    temperature_kelvin = temperature_grid + 273.15
    smr = np.exp((0.63 / kb) * (1.0 / 293.15 - 1.0 / temperature_kelvin))
    for taxon, optimum_temperature, aerobic_scope_requirement, response_shape in (
        ("Scarus", 30.0, 4.0, 1.5),
        ("Cephalopholis", 29.0, 3.0, 0.8),
    ):
        smr_at_optimum = np.exp((0.63 / kb) * (1.0 / 293.15 - 1.0 / (optimum_temperature + 273.15)))
        maximum_metabolic_rate = (
            aerobic_scope_requirement
            * smr_at_optimum
            * np.exp(-((temperature_grid - optimum_temperature) ** 2) / (2.0 * 5.0**2))
            * (oxygen_grid / (50.0 + oxygen_grid))
        )
        aerobic_scope = np.maximum(maximum_metabolic_rate - smr, 0.0)
        aerobic_scope_maximum = (aerobic_scope_requirement - 1.0) * smr_at_optimum
        physiological_multiplier = np.clip(aerobic_scope / aerobic_scope_maximum, 0.0, 1.0) ** response_shape
        rows.append(
            pd.DataFrame(
                {
                    "taxon": taxon,
                    "response_geometry": "temperature_by_oxygen",
                    "temperature_c": temperature_grid.ravel(),
                    "omega_arag": np.nan,
                    "oxygen_mmol_m3": oxygen_grid.ravel(),
                    "physiological_multiplier": physiological_multiplier.ravel(),
                    "effective_hsi_multiplier": physiological_multiplier.ravel() ** outer_gamma,
                }
            )
        )
    return pd.concat(rows, ignore_index=True)


def prepare_iucn_constraint_effect() -> tuple[pd.DataFrame, np.ndarray, rasterio.profiles.Profile, list[Path]]:
    """Compare constrained and unconstrained HSI on their common finite domain."""
    authority = pd.read_csv(IUCN_AUTHORITY).set_index("species")
    rows: list[dict[str, float | int | str]] = []
    shared_domain: np.ndarray | None = None
    profile: rasterio.profiles.Profile | None = None
    sources: list[Path] = [IUCN_AUTHORITY]
    for taxon in TAXA:
        temporal_dir = SDM_ROOT / "maxent" / "temporal" / f"{taxon}_global"
        constrained_path = temporal_dir / "temporal_mean.tif"
        unconstrained_path = temporal_dir / "temporal_mean_unconstrained.tif"
        iucn_path = IUCN_DIR / f"{IUCN_TAXON[taxon]}_global_iucn_range.tif"
        for path in (constrained_path, unconstrained_path, iucn_path):
            if not path.exists():
                raise FileNotFoundError(path)
        sources.extend([constrained_path, unconstrained_path, iucn_path])
        with rasterio.open(constrained_path) as constrained_source, rasterio.open(unconstrained_path) as unconstrained_source:
            if (
                constrained_source.width != unconstrained_source.width
                or constrained_source.height != unconstrained_source.height
                or constrained_source.transform != unconstrained_source.transform
                or constrained_source.crs != unconstrained_source.crs
            ):
                raise ValueError(f"Constrained and unconstrained grids differ for {taxon}")
            constrained, constrained_valid = finite_values(constrained_source.read(1), constrained_source.nodata)
            unconstrained, unconstrained_valid = finite_values(unconstrained_source.read(1), unconstrained_source.nodata)
            observed_range = align_iucn(iucn_path, constrained_source)
            if profile is None:
                profile = constrained_source.profile.copy()
        valid = constrained_valid & unconstrained_valid & np.isfinite(observed_range)
        observed = observed_range[valid] > 0.5
        if not valid.any():
            raise ValueError(f"No common finite validation pixels for {taxon}")
        for state, values in (("Unconstrained MaxEnt", unconstrained), ("Physiologically constrained HSI", constrained)):
            metrics = binary_metrics(values[valid] >= OR10[taxon], observed)
            rows.append(
                {
                    "taxon": taxon,
                    "iucn_comparator": IUCN_TAXON[taxon],
                    "validation_role": IUCN_ROLE[taxon],
                    "state": state,
                    "threshold_or10": OR10[taxon],
                    "common_valid_pixels": int(valid.sum()),
                    **metrics,
                }
            )
        constrained_row = rows[-1]
        reference = authority.loc[taxon]
        delta = abs(int(constrained_row["TP"]) - int(reference["TP"])) + abs(int(constrained_row["FP"]) - int(reference["FP"]))
        if delta > 2:
            raise ValueError(f"Current constrained validation mismatch exceeds two pixels for {taxon}")
        shared_domain = valid if shared_domain is None else (shared_domain | valid)
    if profile is None or shared_domain is None:
        raise RuntimeError("No IUCN constraint data were prepared")
    return pd.DataFrame(rows), shared_domain.astype("uint8"), profile, sources


def bootstrap_auc_interval(presence: np.ndarray, hsi: np.ndarray, seed: int = 42, n_bootstrap: int = 1000) -> tuple[float, float, float]:
    """Estimate an event-level ROC AUC and a stratified percentile interval."""
    presence = np.asarray(presence, dtype=np.int8)
    hsi = np.asarray(hsi, dtype=float)
    positive = hsi[presence == 1]
    negative = hsi[presence == 0]
    if not len(positive) or not len(negative):
        raise ValueError("Caribbean event-level validation requires both presence states")
    auc = float(roc_auc_score(presence, hsi))
    generator = np.random.default_rng(seed)
    bootstrap = np.empty(n_bootstrap, dtype=float)
    for index in range(n_bootstrap):
        resampled_positive = positive[generator.integers(0, len(positive), len(positive))]
        resampled_negative = negative[generator.integers(0, len(negative), len(negative))]
        bootstrap[index] = roc_auc_score(
            np.r_[np.ones(len(resampled_positive), dtype=np.int8), np.zeros(len(resampled_negative), dtype=np.int8)],
            np.r_[resampled_positive, resampled_negative],
        )
    lower, upper = np.quantile(bootstrap, [0.025, 0.975])
    return auc, float(lower), float(upper)


def prepare_survey() -> tuple[pd.DataFrame, pd.DataFrame, list[Path]]:
    """Prepare direct, event-level Caribbean external-validation results.

    A survey event is the independent unit. Repeated taxon rows within a
    survey are reduced to a single presence indicator and mean sampled HSI.
    This avoids comparing the survey series with manuscript trend products.
    """
    if not SURVEY_EVENT_HSI.exists():
        raise FileNotFoundError(SURVEY_EVENT_HSI)
    records = pd.read_csv(SURVEY_EVENT_HSI)
    required = {"survey_id", "year", "taxon", "present", "hsi_0083"}
    if not required.issubset(records.columns):
        raise ValueError("Unexpected Caribbean survey-event validation schema")
    records = records.loc[records["taxon"].isin(TAXA)].copy()
    records["present"] = pd.to_numeric(records["present"], errors="raise").astype("int8")
    if not records["present"].isin((0, 1)).all():
        raise ValueError("Caribbean survey presence must be binary")
    records["hsi_0083"] = pd.to_numeric(records["hsi_0083"], errors="coerce")
    finite = records.dropna(subset=["hsi_0083"]).copy()
    events = (
        finite.groupby(["taxon", "survey_id", "year"], as_index=False)
        .agg(present=("present", "max"), hsi_0083=("hsi_0083", "mean"), n_source_rows=("present", "size"))
        .sort_values(["taxon", "year", "survey_id"])
        .reset_index(drop=True)
    )
    rows = []
    for taxon in TAXA:
        subset = events.loc[events["taxon"].eq(taxon)].copy()
        auc, auc_low, auc_high = bootstrap_auc_interval(subset["present"].to_numpy(), subset["hsi_0083"].to_numpy())
        rows.append(
            {
                "taxon": taxon,
                "n_unique_survey_events": int(len(subset)),
                "n_present": int(subset["present"].sum()),
                "n_absent": int((subset["present"] == 0).sum()),
                "n_source_rows": int(subset["n_source_rows"].sum()),
                "auc": auc,
                "auc_ci95_low": auc_low,
                "auc_ci95_high": auc_high,
                "median_hsi_present": float(subset.loc[subset["present"].eq(1), "hsi_0083"].median()),
                "median_hsi_absent": float(subset.loc[subset["present"].eq(0), "hsi_0083"].median()),
                "year_min": int(subset["year"].min()),
                "year_max": int(subset["year"].max()),
            }
        )
    return events, pd.DataFrame(rows), [SURVEY_EVENT_HSI]


def prepare_alternative_sdm() -> tuple[pd.DataFrame, list[Path]]:
    """Read the completed current-branch BRT/RF/GLM agreement table."""
    if not ALTERNATIVE_AGREEMENT.exists():
        raise FileNotFoundError(ALTERNATIVE_AGREEMENT)
    agreement = pd.read_csv(ALTERNATIVE_AGREEMENT)
    required = {
        "taxon",
        "algorithm",
        "comparison_domain",
        "n_common_finite_cells",
        "spearman_rho",
    }
    if not required.issubset(agreement.columns):
        raise ValueError("Unexpected alternative-SDM agreement schema")
    expected = pd.MultiIndex.from_product([TAXA, ("BRT", "RF", "GLM")], names=["taxon", "algorithm"])
    found = pd.MultiIndex.from_frame(agreement.loc[:, ["taxon", "algorithm"]])
    if not found.is_unique or set(found) != set(expected):
        raise ValueError("Alternative-SDM agreement table does not contain one result per taxon and algorithm")
    if not agreement["spearman_rho"].between(-1.0, 1.0).all():
        raise ValueError("Alternative-SDM Spearman correlations are outside [-1, 1]")
    if not agreement["comparison_domain"].eq("Common finite support across MaxEnt, BRT, RF, and GLM").all():
        raise ValueError("Alternative-SDM results do not use the strict common comparison domain")
    return agreement.sort_values(["taxon", "algorithm"]).reset_index(drop=True), [ALTERNATIVE_AGREEMENT]


def write_study_domain(mask: np.ndarray, profile: rasterio.profiles.Profile) -> Path:
    output = RASTER_DIR / "figS1_iucn_common_validation_domain_0083.tif"
    profile.update(count=1, dtype="uint8", nodata=0, compress="lzw")
    with rasterio.open(output, "w", **profile) as destination:
        destination.write(mask, 1)
    return output


def write_availability_audit() -> None:
    audit = """# Supplementary Fig. S1 current-test-branch data audit

## Available and prepared

- Thinned 0.083-degree occurrence coordinates for all four focal taxa.
- Current physiological equations and stated parameter values from the 0.083-degree projection script and Supplementary Methods.
- Constrained and unconstrained 1993-2022 mean HSI TIFFs, paired with IUCN ranges on their common finite domain.
- Caribbean external-survey event data for 2013-2022, with sampled annual 0.083-degree HSI and binary observed occurrence.
- Completed current-branch BRT, random-forest and binomial-GLM models, their 0.083-degree predictions and a strict common-domain MaxEnt agreement table.

## Alternative-SDM concordance

The three alternative algorithms were refitted using the restored corrected-20260722 presence and target-group-background SWD tables, the taxon-specific current selected predictors, and a cellwise 1993-2022 mean predictor baseline. Continuous baseline predictors were bilinearly aligned from their native 0.5-degree grid to the current 0.083-degree MaxEnt template before prediction. The panel uses Spearman correlations over the finite intersection of MaxEnt, BRT, RF and GLM predictions for each taxon; it does not mix support domains among algorithms.

## Interpretation boundary

The Lobophora comparison uses the Turbinaria IUCN range as a proxy and remains supplementary only. The IUCN constraint-effect comparison applies each taxon's fixed constrained OR10 threshold to both layers on their common finite domain. Caribbean validation uses survey events as the independent unit and reports presence-versus-absence discrimination; it is observational predictive validation, not a causal test.
"""
    (DOC_DIR / "figS1_data_availability_audit_20260806.md").write_text(audit, encoding="utf-8")


def write_figure_spec() -> None:
    spec = """# Supplementary Fig. S1 figure specification

## Scientific question

Does the current 0.083-degree SDM workflow provide traceable record coverage, explicit physiological constraints and geographically or field-based predictive checks for the four focal taxa?

## Layout and inputs

- Canvas: 210 mm wide by 195 mm high. Every data rectangle is 40 x 40 mm; labels, legends and colourbars are outside these rectangles.
- First row: four global occurrence maps of the thinned 0.083-degree records, ordered Acropora, Lobophora, Scarus and Cephalopholis from left to right.
- Second row: the matching four physiological-response panels. Acropora is temperature by aragonite saturation; Lobophora is the temperature response; Scarus and Cephalopholis are temperature by dissolved oxygen. Heatmaps use the shared 0 to 1 effective multiplier scale.
- Third row left: paired direct IUCN range-agreement values for specificity and precision before and after the physiological constraint. The Lobophora-Turbinaria comparison remains a proxy.
- Third row centre: taxon-specific point ranges of Spearman spatial agreement between MaxEnt and BRT, RF and binomial GLM, calculated on the strict common finite support.
- Third row right: the direct Caribbean external-validation result: distributions of 0.083-degree HSI sampled at unique survey events, grouped by observed presence and absence. This panel shows validation-event data rather than a comparison with manuscript time-series products.

## Colour and scale rules

- Taxa: Acropora `#115FA4`, Lobophora `#F79015`, Scarus `#DAA81C`, Cephalopholis `#AB2428`.
- Physiological multiplier: a perceptually ordered scicomap palette with fixed limits 0 to 1 for all heatmaps.
- IUCN state: unconstrained `#7F9DC3`, constrained `#115FA4`; proxy comparison uses a hollow marker and dashed connector.
- Survey-event absence and presence: `#C7E5F2` and `#115FA4`, respectively. BRT, RF and binomial GLM use `#115FA4`, `#F79015` and `#AB2428`.
- No gridlines, panel letters, titles, fitted trend lines, numerical annotations or significance stars. The caption carries the fixed column order, sample sizes, proxy status and statistic definitions.
"""
    (DOC_DIR / "figS1_figure_spec_20260806.md").write_text(spec, encoding="utf-8")


def main() -> None:
    for directory in (TABLE_DIR, RASTER_DIR, DOC_DIR):
        directory.mkdir(parents=True, exist_ok=True)

    occurrences, occurrence_summary, occurrence_sources = prepare_occurrences()
    physiology = prepare_physiology()
    constraint_effect, common_domain, profile, iucn_sources = prepare_iucn_constraint_effect()
    survey_events, survey_summary, survey_sources = prepare_survey()
    alternative_sdm, alternative_sources = prepare_alternative_sdm()
    domain_path = write_study_domain(common_domain, profile)

    occurrences.to_csv(TABLE_DIR / "figS1_occurrence_records_thinned_0083.csv", index=False)
    occurrence_summary.to_csv(TABLE_DIR / "figS1_occurrence_summary.csv", index=False)
    physiology.to_csv(TABLE_DIR / "figS1_physiological_response_curves.csv", index=False)
    constraint_effect.to_csv(TABLE_DIR / "figS1_iucn_constraint_effect_common_domain.csv", index=False)
    survey_events.to_csv(TABLE_DIR / "figS1_caribbean_event_hsi_validation.csv", index=False)
    survey_summary.to_csv(TABLE_DIR / "figS1_caribbean_event_hsi_validation_summary.csv", index=False)
    (TABLE_DIR / "figS1_caribbean_validation_summary.csv").unlink(missing_ok=True)
    (TABLE_DIR / "figS1_caribbean_validation_annual.csv").unlink(missing_ok=True)
    alternative_sdm.to_csv(TABLE_DIR / "figS1_alternative_sdm_agreement.csv", index=False)
    (TABLE_DIR / "figS1_alternative_sdm_agreement_status.csv").unlink(missing_ok=True)

    source_paths = [*occurrence_sources, METHODS_FILE, PROJECTION_SCRIPT, *iucn_sources, *survey_sources, *alternative_sources]
    manifest_rows = []
    for source in source_paths:
        if not source.exists():
            raise FileNotFoundError(source)
        manifest_rows.append(
            {
                "source_path": relative_path(source),
                "sha256": sha256(source),
                "bytes": source.stat().st_size,
            }
        )
    manifest_rows.append(
        {
            "source_path": relative_path(domain_path),
            "sha256": sha256(domain_path),
            "bytes": domain_path.stat().st_size,
        }
    )
    pd.DataFrame(manifest_rows).drop_duplicates("source_path").to_csv(TABLE_DIR / "figS1_input_manifest_20260806.csv", index=False)
    write_availability_audit()
    write_figure_spec()
    print(f"Prepared data under {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
