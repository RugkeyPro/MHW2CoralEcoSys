"""Reproduce reef priority selection and joint model/spatial-block uncertainty."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from . import conservation_core as core
except ImportError:
    import conservation_core as core

MODULE = Path(__file__).resolve().parents[1]
REPO = MODULE.parents[2]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_inputs(path: Path) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    """Validate required physical fields and original per-member support contracts."""
    frame = pd.read_csv(path, float_precision='round_trip')
    required = ['cell_id', 'longitude', 'latitude', 'reef_area_km2',
                'mpa_fraction_with_point_buffers', 'mpei_total', 'mpei_supported',
                'mpei_e_star_global_reef_minmax']
    for scenario in core.SCENARIOS:
        required += [f'P_valid_{scenario}', f'P_{scenario}']
        required += [f'{scenario}__{member}' for member in core.MEMBERS]
    missing = set(required) - set(frame)
    if missing:
        raise ValueError(f'Missing input columns: {sorted(missing)}')
    if frame.empty or not frame.cell_id.is_unique:
        raise ValueError('Input requires nonempty, unique reef-cell identifiers')
    for field, low, high in [('longitude', -180, 180), ('latitude', -90, 90),
                             ('mpa_fraction_with_point_buffers', 0, 1),
                             ('mpei_e_star_global_reef_minmax', 0, 1)]:
        values = frame[field].to_numpy(float)
        if not np.isfinite(values).all() or np.any(values < low - 1e-10) or np.any(values > high + 1e-10):
            raise ValueError(f'Invalid input range or missing value: {field}')
    area = frame.reef_area_km2.to_numpy(float)
    if not np.isfinite(area).all() or np.any(area <= 0):
        raise ValueError('Reef area must be finite and strictly positive, in km2')
    boolean_fields = ['mpei_supported', *[f'P_valid_{scenario}' for scenario in core.SCENARIOS]]
    for field in boolean_fields:
        if not pd.api.types.is_bool_dtype(frame[field]):
            raise ValueError(f'{field} must contain explicit True/False values; missing flags are not allowed')
    supported = frame.mpei_supported.to_numpy(bool)
    mpei = frame.mpei_total.to_numpy(float)
    if not np.isfinite(mpei[supported]).all() or np.any(mpei[supported] < 0):
        raise ValueError('Supported MPEI must be finite and nonnegative, in kg/m3')
    masks = {}
    for scenario in core.SCENARIOS:
        matrix = frame[[f'{scenario}__{member}' for member in core.MEMBERS]].to_numpy(float).T
        finite = np.isfinite(matrix)
        if not np.isin(matrix[finite], [0.0, 1.0]).all():
            raise ValueError('Member support must be 0, 1 or NaN for missing environmental support')
        valid = np.all(finite, axis=0)
        if not np.array_equal(valid, frame[f'P_valid_{scenario}'].to_numpy(bool)):
            raise ValueError('Validity flag disagrees with complete member support')
        np.testing.assert_allclose(np.mean(matrix, axis=0), frame[f'P_{scenario}'].to_numpy(float),
                                   rtol=0, atol=1e-15, equal_nan=True)
        masks[scenario] = matrix
    return frame, masks


def compare_frames(actual: pd.DataFrame, reference: pd.DataFrame, keys: list[str]) -> dict:
    """Compare aligned source records, preserving missingness and metric units."""
    if actual.duplicated(keys).any() or reference.duplicated(keys).any():
        raise ValueError('Duplicate comparison keys')
    left, right = actual.sort_values(keys).reset_index(drop=True), reference.sort_values(keys).reset_index(drop=True)
    if len(left) != len(right) or not left[keys].equals(right[keys]):
        raise ValueError('Recomputed record keys differ from frozen reference')
    fields = [column for column in right if column not in keys]
    max_errors = {}
    finite_comparisons = 0
    for field in fields:
        a, b = left[field].to_numpy(float), right[field].to_numpy(float)
        if not np.array_equal(np.isnan(a), np.isnan(b)):
            raise ValueError(f'Changed missingness in {field}')
        atol = 1e-18 if field.endswith('kg_m3') else 1e-8
        np.testing.assert_allclose(a, b, rtol=1e-9, atol=atol, equal_nan=True,
                                   err_msg=f'Frozen reference mismatch: {field}')
        finite = np.isfinite(a) & np.isfinite(b)
        finite_comparisons += int(finite.sum())
        max_errors[field] = float(np.max(np.abs(a[finite] - b[finite]))) if finite.any() else None
    return {'status': 'passed', 'records': len(left), 'numeric_fields': len(fields),
            'finite_values_compared': finite_comparisons, 'max_absolute_error_by_field': max_errors}


def compact_table(summary: pd.DataFrame) -> pd.DataFrame:
    scenario_names = dict(zip(core.SCENARIOS, ['SSP1-2.6', 'SSP2-4.5', 'SSP5-8.5']))
    records = []
    for region in core.REGIONS:
        for scenario in core.SCENARIOS:
            rows = summary.loc[(summary.region == region) & (summary.scenario == scenario)].set_index('metric')
            reallocation = rows.loc['future_mpei_vs_climate_area_changed_fraction']
            exposure = rows.loc['future_mpei_vs_climate_mean_mpei_change_pct']
            records.append({'region': region, 'scenario': scenario_names[scenario],
                'future_mpei_vs_climate_area_changed_pct': float(reallocation.estimate) * 100,
                'area_changed_lower_95': float(reallocation.lower_95) * 100,
                'area_changed_upper_95': float(reallocation.upper_95) * 100,
                'future_mpei_vs_climate_mean_mpei_change_pct': float(exposure.estimate),
                'mpei_change_lower_95': float(exposure.lower_95),
                'mpei_change_upper_95': float(exposure.upper_95)})
    return pd.DataFrame(records)


def export_priorities(reef: pd.DataFrame, output: Path) -> pd.DataFrame:
    records = []
    for scenario in core.SCENARIOS:
        candidate = reef[f'P_valid_{scenario}'] & reef[f'P_{scenario}'].gt(0)
        cells = reef.loc[candidate].copy()
        _, selected = core.calculate_metrics(cells.reef_area_km2.to_numpy(float),
            np.clip(cells.mpa_fraction_with_point_buffers.to_numpy(float), 0, 1),
            cells[f'P_{scenario}'].to_numpy(float), cells.mpei_e_star_global_reef_minmax.to_numpy(float),
            cells.mpei_total.to_numpy(float), cells.mpei_supported.to_numpy(bool))
        subset = cells[['cell_id', 'longitude', 'latitude', 'reef_area_km2', 'mpei_total', 'mpei_supported']].copy()
        subset['scenario'] = scenario
        subset['model_support_fraction'] = cells[f'P_{scenario}'].to_numpy(float)
        for scheme, values in selected.items():
            subset[f'{scheme}_selected_fraction'] = values
        subset['selection_fraction_change'] = selected['future_mpei_priority'] - selected['climate_priority']
        records.append(subset)
    table = pd.concat(records, ignore_index=True)
    table.to_csv(output / 'selected_reef_cells.csv.gz', index=False, float_format='%.17g', compression='gzip')
    return table


def render_figures(table: pd.DataFrame, cells: pd.DataFrame, output: Path) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'svg.fonttype': 'none',
                         'axes.spines.top': False, 'axes.spines.right': False})
    names = ['Global', 'Caribbean', 'GBR', 'Southeast Asia']
    colors = ['#377b68', '#bb9750', '#c56750']
    scenarios = ['SSP1-2.6', 'SSP2-4.5', 'SSP5-8.5']
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    fig.subplots_adjust(left=.11, right=.98, bottom=.14, top=.76, wspace=.34)
    metrics = [('future_mpei_vs_climate_area_changed_pct', 'area_changed_lower_95', 'area_changed_upper_95', 'Priority area reallocated (%)'),
               ('future_mpei_vs_climate_mean_mpei_change_pct', 'mpei_change_lower_95', 'mpei_change_upper_95', 'Selected-area mean MPEI change (%)')]
    for ax, (metric, lower, upper, title) in zip(axes, metrics):
        for sindex, scenario in enumerate(scenarios):
            subset = table.loc[table.scenario == scenario].set_index('region').loc[names]
            y = np.arange(len(names)) + (sindex - 1) * .20
            ax.hlines(y, subset[lower], subset[upper], color=colors[sindex], linewidth=1.5)
            ax.scatter(subset[metric], y, color=colors[sindex], s=24, zorder=3)
        ax.axvline(0, color='#a4aaa0', linewidth=.8, linestyle='--')
        ax.set_yticks(np.arange(4), names)
        ax.invert_yaxis()
        ax.set_xlabel(title)
        ax.grid(axis='x', alpha=.15)
    handles = [Line2D([0], [0], color=color, marker='o', label=scenario) for color, scenario in zip(colors, scenarios)]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .91), ncol=3, frameon=False)
    fig.suptitle('Recomputed equal-area reef priorities: estimates and 95% bootstrap intervals', y=.98, fontsize=12)
    for suffix in ('png', 'svg', 'pdf'):
        fig.savefig(output / f'conservation_comparison.{suffix}', dpi=180, bbox_inches='tight')
    plt.close(fig)
    subset = cells.loc[cells.scenario == 'ssp585_2050']
    fig, ax = plt.subplots(figsize=(12, 4), layout='constrained')
    colors_map = ax.scatter(subset.longitude, subset.latitude, c=subset.selection_fraction_change,
                           cmap='BrBG', vmin=-1, vmax=1, s=3, linewidths=0)
    ax.set(xlim=(-180, 180), ylim=(-45, 45), xlabel='Longitude (degrees)', ylabel='Latitude (degrees)',
           title='Actual candidate reef cells: exposure-aware minus habitat-only selection (SSP5-8.5)')
    fig.colorbar(colors_map, ax=ax, label='Change in selected cell fraction')
    ax.grid(alpha=.15)
    fig.savefig(output / 'reef_priority_reallocation.png', dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict:
    start = time.perf_counter()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((MODULE / 'provenance.json').read_text(encoding='utf-8'))
    # References are used only after computation, never as analysis inputs.
    for item in manifest['files']:
        if sha(REPO / item['file']) != item['sha256']:
            raise ValueError(f'Packaged file hash mismatch: {item["file"]}')
    input_path = args.input.resolve()
    original_input = sha(input_path) == next(item['sha256'] for item in manifest['files'] if item['role'] == 'runtime_input')
    print('1/5 Loading actual reef-cell inputs', flush=True)
    reef, masks = load_inputs(input_path)
    print(f'    {len(reef):,} cells; {len(core.MEMBERS)} climate members; 3 scenarios', flush=True)
    print('2/5 Computing equal-area priorities and exposure metrics', flush=True)
    point, _ = core.point_estimates(reef, [])
    point.to_csv(output / 'point_estimates.csv', index=False, float_format='%.17g')
    if point.maximum_area_budget_closure_error_km2.max() > 1e-7:
        raise ValueError('Selected-area budget does not close')
    np.testing.assert_allclose(point.future_mpei_priority_captured_co_suitability_fraction,
                               point.climate_priority_captured_co_suitability_fraction, rtol=0, atol=1e-12)
    priorities = export_priorities(reef, output)
    print(f'3/5 Running {args.replicates} joint climate-member / 5-degree spatial-block bootstrap draws', flush=True)
    core.N_BOOTSTRAP, core.SEED = args.replicates, args.seed
    draws = core.bootstrap(reef, masks)
    draws.to_csv(output / 'bootstrap_draws.csv.gz', index=False, float_format='%.17g', compression='gzip')
    summary = core.summarize(point, draws)
    # The accepted promotion wrapper records this diagnostic after summarizing.
    summary['outside_percentile_interval'] = (
        (summary.estimate < summary.lower_95 - 1e-10)
        | (summary.estimate > summary.upper_95 + 1e-10)
    )
    summary.to_csv(output / 'summary.csv', index=False, float_format='%.17g')
    compact = compact_table(summary)
    compact.to_csv(output / 'Table_S5_recomputed.csv', index=False, float_format='%.17g')
    print('4/5 Checking numerical outputs against the frozen accepted references', flush=True)
    checks = {'schema': 'passed', 'area_budget_closure': 'passed', 'support_tier_retention': 'passed',
              'input_rows': len(reef), 'point_rows': len(point), 'draw_rows': len(draws),
              'requested_replicates': args.replicates, 'seed': args.seed, 'original_packaged_input': original_input,
              'full_frozen_design': original_input and args.replicates == 1000 and args.seed == 20260804,
              'source_hashes_checked': len(manifest['files']), 'comparisons': {}}
    if original_input:
        checks['comparisons']['point_estimates'] = compare_frames(point, pd.read_csv(MODULE / 'reference/point_estimates.csv', float_precision='round_trip'), ['region','scenario'])
    if original_input and args.seed == 20260804 and args.replicates <= 1000:
        frozen = pd.read_csv(MODULE / 'reference/bootstrap_draws.csv.gz', float_precision='round_trip')
        checks['comparisons']['bootstrap_draws'] = compare_frames(draws, frozen.loc[frozen.replicate < args.replicates], ['replicate','region','scenario'])
    if checks['full_frozen_design']:
        checks['comparisons']['summary'] = compare_frames(summary, pd.read_csv(MODULE / 'reference/summary.csv', float_precision='round_trip'), ['region','scenario','metric'])
        checks['comparisons']['Table_S5'] = compare_frames(compact, pd.read_csv(MODULE / 'reference/Table_S5.csv', float_precision='round_trip'), ['region','scenario'])
    checks['effective_replicates_by_region'] = {region: int(summary.loc[summary.region == region].n_bootstrap.min()) for region in core.REGIONS}
    checks['reference_scope'] = 'Full accepted calculation' if checks['full_frozen_design'] else 'Customized/shortened run; frozen full uncertainty intervals are not reproduced'
    print('5/5 Writing calculated tables, selected reef cells and scientific figures', flush=True)
    if not args.no_plots:
        render_figures(compact, priorities, output)
    checks.update(status='passed', elapsed_seconds=round(time.perf_counter() - start, 3),
                  environment={'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__},
                  raw_inputs_modified=False)
    (output / 'run_checks.json').write_text(json.dumps(checks, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps({key: checks[key] for key in ['status','input_rows','draw_rows','full_frozen_design','effective_replicates_by_region','elapsed_seconds']}, indent=2), flush=True)
    return checks


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=MODULE / 'data/reef_cells.csv.gz')
    parser.add_argument('--output', type=Path, default=REPO / 'outputs/conservation_demo')
    parser.add_argument('--replicates', type=int, default=1000, help='Default reproduces the full accepted 1,000-replicate design')
    parser.add_argument('--seed', type=int, default=20260804)
    parser.add_argument('--no-plots', action='store_true')
    args = parser.parse_args()
    if args.replicates < 1 or args.seed < 0:
        parser.error('Replicates must be positive and seed must be nonnegative')
    try:
        run(args)
    except Exception as error:
        args.output.mkdir(parents=True, exist_ok=True)
        (args.output / 'run_checks.json').write_text(json.dumps({'status':'failed','error':str(error)}, indent=2) + '\n', encoding='utf-8')
        raise


if __name__ == '__main__':
    main()
