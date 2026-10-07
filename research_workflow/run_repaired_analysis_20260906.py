"""Use repaired August habitat inputs with the mean 2045-2055 MPEI field.

Reuse the validated habitat builder and period-matched exposure implementation.
No habitat thresholds, selection rules, or bootstrap design are changed.
"""
from pathlib import Path
import hashlib
import importlib.util
import json

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
WORKTREE = ROOT.parents[2]
REPAIR = Path(r"D:\ecoMarine\output\02_experiments\20260816_future_coverage_repair\figureS6_protection_configuration_repaired_20260821")
EXPOSURE = WORKTREE / "output/01_mainline/future_conservation_period_matched_mpei_20260903"
MPEI = Path(r"D:\ecoMarine_worktrees\20260903_mainline_publish\output\01_mainline\inputs\mpei_36tracer_annual_1950_2060_20260901\derived\mpei_36tracer_0_100m_mean_2045_2055_20260901.nc")


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def compare_habitat(new, old, keys):
    """Verify that only MP-dependent metrics changed, including every bootstrap draw."""
    metrics = ["candidate_reef_area_km2", "mpa_budget_km2",
               "existing_mpa_captured_co_suitability_fraction",
               "climate_priority_captured_co_suitability_fraction"]
    paired = new.merge(old, on=keys, validate="one_to_one", suffixes=("_new", "_old"))
    assert len(paired) == len(new) == len(old)
    errors = {}
    for metric in metrics:
        a, b = paired[metric + "_new"], paired[metric + "_old"]
        np.testing.assert_allclose(a, b, rtol=1e-9, atol=1e-10, equal_nan=True)
        errors[metric] = float((a - b).abs().max())
    return errors


def main():
    for name in ("tables", "derived", "qa", "docs", "source"):
        (ROOT / name).mkdir(parents=True, exist_ok=True)
    habitat_script = REPAIR / "scripts/prepare_figureS6_repaired_20260821.py"
    exposure_script = EXPOSURE / "scripts/process_future_mpei_conservation_20260903.py"
    habitat = load_module("repaired_habitat", habitat_script)
    exposure = load_module("period_matched_exposure", exposure_script)
    # Explicit run configuration: retain the August bootstrap seed, while using
    # the validated future-MPEI source and interpolation support rules.
    exposure.MPEI_FILE = MPEI
    exposure.SEED = habitat.SEED
    assert exposure.SCENARIOS == habitat.SCENARIOS
    assert exposure.MEMBERS == habitat.MEMBERS
    assert exposure.REGIONS == habitat.REGIONS
    assert exposure.SCORE_EPSILON == habitat.SCORE_EPSILON
    assert exposure.N_BOOTSTRAP == habitat.N_BOOTSTRAP == 1000
    assert exposure.SPATIAL_BLOCK_DEG == habitat.BLOCK_DEG

    reef = pd.read_parquet(habitat.REEF_MPA)
    reef = reef.loc[reef.reef_area_km2.gt(0)].copy().reset_index(drop=True)
    assert reef.cell_id.is_unique
    assert np.isfinite(reef.mpa_fraction_with_point_buffers).all()
    masks, continuous_585, inventory = habitat.build_current_fields(reef)
    for name, field in (("figS8a_mpa_reef_fraction_05deg_20260903.tif",
                         np.clip(reef.mpa_fraction_with_point_buffers.to_numpy(float), 0, 1)),
                        ("figS8b_continuous_four_taxon_hsi_ssp585_05deg_20260903.tif", continuous_585)):
        grid, _ = habitat.aggregate_display(reef, field)
        habitat.write_tif(ROOT / "source" / name, grid)
    score_table = reef[["cell_id", "row", "col", "longitude", "latitude"]].copy()
    for scenario, matrix in masks.items():
        valid = np.all(np.isfinite(matrix), axis=0)
        reef[f"P_valid_{scenario}"] = valid
        reef[f"P_{scenario}"] = np.mean(matrix, axis=0)
        for index, member in enumerate(habitat.MEMBERS):
            score_table[f"{scenario}__{member}"] = matrix[index]
        score_table[f"P_{scenario}"] = reef[f"P_{scenario}"]
    lat, lon, fields, metadata = exposure.read_native_mpei()
    assert metadata["completed_source_year_start"] == 2045
    assert metadata["completed_source_year_end"] == 2055
    assert metadata["annual_endpoint_count"] == 11
    values, interpolation_qc = exposure.interpolate_to_reef_cells(reef, lat, lon, fields)
    for column in values:
        reef[column] = values[column].to_numpy()
    supported = reef.mpei_supported.to_numpy(bool)
    minimum, maximum = reef.loc[supported, "mpei_total"].agg(["min", "max"])
    assert np.isfinite([minimum, maximum]).all() and 0 <= minimum < maximum
    reef["mpei_e_star_global_reef_minmax"] = 0.0
    reef.loc[supported, "mpei_e_star_global_reef_minmax"] = (
        reef.loc[supported, "mpei_total"] - minimum
    ) / (maximum - minimum)
    components = [c for c in reef if c.startswith(("mpei_polymer_", "mpei_size_", "mpei_state_"))]
    point, shares = exposure.point_estimates(reef, components)
    sensitivity, _ = exposure.point_estimates(reef, [], supported_domain_only=True)
    old_point = pd.read_csv(REPAIR / "tables/figS6_protection_point_estimates_20260821.csv")
    point_errors = compare_habitat(point, old_point, ["region", "scenario"])
    print("PASS: repaired habitat point estimates match August manuscript source", flush=True)
    draws = exposure.bootstrap(reef, masks)
    old_draws = pd.read_csv(REPAIR / "tables/figS6_protection_joint_bootstrap_20260821.csv")
    draw_errors = compare_habitat(draws, old_draws, ["replicate", "region", "scenario"])
    summary = exposure.summarize(point, draws)
    summary["outside_percentile_interval"] = (
        (summary.estimate < summary.lower_95 - 1e-10)
        | (summary.estimate > summary.upper_95 + 1e-10)
    )
    table_metrics = {
        "existing_mpa_captured_co_suitability_fraction": 100,
        "climate_priority_captured_co_suitability_fraction": 100,
        "future_mpei_priority_captured_co_suitability_fraction": 100,
        "future_mpei_vs_climate_area_changed_fraction": 100,
        "future_mpei_vs_climate_mean_mpei_change_pct": 1,
        "future_mpei_vs_existing_mean_mpei_change_pct": 1,
    }
    s5 = summary.loc[summary.metric.isin(table_metrics)].copy()
    for col in ["estimate", "lower_95", "upper_95"]:
        s5[col] *= s5.metric.map(table_metrics)
    s5["unit"] = "%"
    # Independently check exported percentage intervals against stored draws.
    for row in s5.itertuples():
        values = draws.loc[(draws.region == row.region)
                           & (draws.scenario == row.scenario), row.metric].dropna()
        assert np.isfinite(values).all() and len(values) == row.n_bootstrap
        np.testing.assert_allclose(
            np.quantile(values, [0.025, 0.975]) * table_metrics[row.metric],
            [row.lower_95, row.upper_95], atol=1e-9,
        )
    for part, metrics in (("a", list(table_metrics)[:4]),
                          ("b", list(table_metrics)[4:])):
        wide = s5.loc[s5.metric.isin(metrics)].pivot(
            index=["region", "scenario"], columns="metric",
            values=["estimate", "lower_95", "upper_95"],
        )
        wide.columns = [f"{metric}__{stat}" for stat, metric in wide.columns]
        wide.to_csv(ROOT / f"tables/Table_S5{part}_wide_percent.csv", encoding="utf-8-sig")
    s5.to_csv(ROOT / "tables/Table_S5_percent.csv", index=False)
    # Replace established consumer paths, so old table/figure commands cannot read stale results.
    for label, frame in {"point_estimates": point, "bootstrap": draws, "summary": summary,
                         "component_shares": shares, "supported_domain_sensitivity": sensitivity}.items():
        frame.to_csv(ROOT / f"tables/future_mpei_conservation_{label}_20260903.csv", index=False)
    summary.to_csv(ROOT / "tables/Supplementary_Table_S5_all_metrics_20260903.csv", index=False)
    reef.to_parquet(ROOT / "derived/future_mpei_mean_2045_2055_reef_cells_20260903.parquet", index=False)
    score_table.to_parquet(ROOT / "derived/repaired_member_support.parquet", index=False)
    assert len(point) == 12 and len(draws) == 12000 and len(s5) == 72
    np.testing.assert_allclose(
        draws.climate_priority_captured_co_suitability_fraction,
        draws.future_mpei_priority_captured_co_suitability_fraction,
        atol=1e-12, rtol=1e-12, equal_nan=True,
    )
    assert draws.maximum_area_budget_closure_error_km2.max() < 1e-6
    assert summary.n_bootstrap.min() >= 990
    assert np.allclose(shares.groupby(["region", "scenario", "scheme", "component_type"])
                       .share_of_selected_future_mpei_percent.sum(), 100, atol=1e-5)
    support_columns = [f"{s}_mpei_supported_captured_co_suitability_fraction"
                       for s in ("existing_mpa", "climate_priority", "future_mpei_priority")]
    min_support = float(point[support_columns].min().min())
    assert min_support >= 0.8
    qc = {
        "numerical_checks": "PASS",
        "intervals_require_review": bool(s5.outside_percentile_interval.any()),
        "table_s5_outside_intervals": s5.loc[s5.outside_percentile_interval].to_dict("records"),
        "point_habitat_difference_vs_august": point_errors,
        "bootstrap_habitat_difference_vs_august": draw_errors,
        "bootstrap_replicates": 1000, "bootstrap_seed": exposure.SEED,
        "minimum_effective_replicates": int(summary.n_bootstrap.min()),
        "minimum_selected_P_weight_with_MP_support_fraction": min_support,
        "mpei_metadata": metadata, "interpolation": interpolation_qc,
        "score": "P + 1e-4 * E_star; P is repaired three-member strict OR10/MESS support",
        "candidate_and_budget": "valid P>0 reef cells; B=sum(reef_area*MPA_fraction) in that region/scenario candidate domain",
        "all_mapped_reef_MPA_fraction": float(np.average(
            np.clip(reef.mpa_fraction_with_point_buffers, 0, 1), weights=reef.reef_area_km2)),
    }
    for path in (habitat_script, exposure_script, habitat.REEF_MPA, habitat.OR10, MPEI):
        inventory.append({"kind": "input_or_implementation", "path": str(path),
                          "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    pd.DataFrame(inventory).to_csv(ROOT / "qa/input_manifest.csv", index=False)
    (ROOT / "qa/validation.json").write_text(json.dumps(qc, indent=2, ensure_ascii=False), encoding="utf-8")
    (ROOT / "qa/future_mpei_conservation_qc_20260903.json").write_text(
        json.dumps(qc, indent=2, ensure_ascii=False), encoding="utf-8")
    from prepare_supplementary_table_s5_20260903 import main as prepare_table
    prepare_table()
    print(json.dumps(qc, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
