"""Package the accepted reef-cell inputs and unchanged scientific calculation core."""
import argparse
import ast
import gzip
import hashlib
import io
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / 'analysis/modules/conservation_demo'
FUNCTIONS = ('exact_budget_selection', 'ratio', 'calculate_metrics', 'region_mask', 'point_estimates', 'bootstrap', 'summarize')
SCENARIOS = ('ssp126_2050', 'ssp245_2050', 'ssp585_2050')
MEMBERS = ('gfdl_esm4_r1i1p1f1_gr', 'cnrm_esm2_1_r1i1p1f2_gn', 'ipsl_cm6a_lr_r1i1p1f1_gn')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def write_gzip_csv(frame, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open('wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as compressed:
        with io.TextIOWrapper(compressed, encoding='utf-8', newline='') as text:
            frame.to_csv(text, index=False, float_format='%.17g', lineterminator='\n')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', type=Path, required=True)
    args = parser.parse_args()
    project = args.project_root
    mainline = project / '解压整理/01_ecoMarine/03_mainline/ecoMarine_0907_20260913_03_mainline_009/filesystem/D/ecoMarine_worktrees/20260903_mainline_publish/output/01_mainline'
    parent = mainline / 'future_conservation_period_matched_mpei_20260903'
    cells_file = parent / 'derived/future_mpei_mean_2045_2055_reef_cells_20260903.parquet'
    masks_file = parent / 'derived/repaired_member_support.parquet'
    source_script = parent / 'scripts/process_future_mpei_conservation_20260903.py'
    source_hashes_before = {f.name: sha(f) for f in (cells_file, masks_file, source_script)}
    reef = pd.read_parquet(cells_file)
    masks = pd.read_parquet(masks_file)
    if not reef.cell_id.is_unique or not masks.cell_id.is_unique:
        raise ValueError('Duplicate source cell identifiers')
    if not np.array_equal(reef.cell_id.to_numpy(), masks.cell_id.to_numpy()):
        raise ValueError('Source row order differs; bootstrap order must remain unchanged')
    fields = ['cell_id', 'row', 'col', 'longitude', 'latitude', 'reef_area_km2', 'mpa_fraction_with_point_buffers',
              *[f'P_valid_{s}' for s in SCENARIOS], *[f'P_{s}' for s in SCENARIOS],
              'mpei_total', 'mpei_supported', 'mpei_e_star_global_reef_minmax']
    selected = reef[fields].copy()
    mask_fields = [f'{s}__{m}' for s in SCENARIOS for m in MEMBERS]
    for column in mask_fields:
        selected[column] = masks[column].to_numpy()
    input_file = MODULE / 'data/reef_cells.csv.gz'
    write_gzip_csv(selected, input_file)
    restored = pd.read_csv(input_file, float_precision='round_trip')
    for column in selected:
        if pd.api.types.is_bool_dtype(selected[column]):
            if not np.array_equal(selected[column].to_numpy(), restored[column].to_numpy()):
                raise ValueError(f'Boolean round-trip mismatch: {column}')
        else:
            np.testing.assert_array_equal(selected[column].to_numpy(), restored[column].to_numpy())
    parent_text = source_script.read_text(encoding='utf-8-sig')
    tree = ast.parse(parent_text)
    definitions = {node.name: node for node in tree.body if isinstance(node, ast.FunctionDef)}
    constants = [node for node in tree.body if isinstance(node, ast.Assign)
                 and any(isinstance(target, ast.Name) and target.id in {'SCENARIOS','MEMBERS','REGIONS','SCORE_EPSILON','N_BOOTSTRAP','SPATIAL_BLOCK_DEG','SEED'} for target in node.targets)]
    core = ['"""Scientific calculation functions copied unchanged from the accepted September 6 workflow.\n\nOnly source IO and plotting are outside this module; see provenance.json for origins.\n"""',
            'from __future__ import annotations', 'import numpy as np', 'import pandas as pd']
    core += [ast.get_source_segment(parent_text, node) for node in constants]
    core += [ast.get_source_segment(parent_text, definitions[name]) for name in FUNCTIONS]
    core_file = MODULE / 'code/conservation_core.py'
    core_file.write_text('\n\n'.join(core) + '\n', encoding='utf-8', newline='\n')
    files = []
    def register(path, role, source):
        files.append({'file': path.relative_to(ROOT).as_posix(), 'role': role, 'sha256': sha(path),
                      'bytes': path.stat().st_size, 'source': source})
    register(input_file, 'runtime_input', 'Accepted derived reef cells and member masks; numeric columns only, original row order and exact float round-trip')
    register(core_file, 'unchanged_function_bodies', 'process_future_mpei_conservation_20260903.py: seven independent calculation functions and original constants')
    for source_name, target_name in [('future_mpei_conservation_point_estimates_20260903.csv', 'point_estimates.csv'),
                                     ('future_mpei_conservation_bootstrap_20260903.csv', 'bootstrap_draws.csv.gz'),
                                     ('future_mpei_conservation_summary_20260903.csv', 'summary.csv'),
                                     ('Supplementary_Table_S5_20260906.csv', 'Table_S5.csv')]:
        source = parent / 'tables' / source_name
        target = MODULE / 'reference' / target_name
        if target.suffix == '.gz':
            with target.open('wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as stream:
                stream.write(source.read_bytes())
            if gzip.decompress(target.read_bytes()) != source.read_bytes():
                raise ValueError('Frozen draw reference compression changed bytes')
        else:
            shutil.copyfile(source, target)
        register(target, 'frozen_reference_not_runtime_input', source.relative_to(mainline).as_posix())
    original_specs = [
        (source_script, 'process_future_mpei_conservation_20260903.py'),
        (parent / 'scripts/run_repaired_analysis_20260906.py', 'run_repaired_analysis_20260906.py'),
        (mainline / 'figure4_future_period_matched_mpei_20260902/scripts/analyze_future_risk_period_matched_mpei_20260902.py', 'analyze_future_risk_period_matched_mpei_20260902.py'),
        (project / '解压整理/01_ecoMarine/03_mainline/ecoMarine_0907_20260913_03_mainline_008/filesystem/D/ecoMarine/output/02_experiments/20260816_future_coverage_repair/figureS6_protection_configuration_repaired_20260821/scripts/prepare_figureS6_repaired_20260821.py', 'prepare_figureS6_repaired_20260821.py'),
    ]
    for source, name in original_specs:
        target = ROOT / 'research_workflow' / name
        shutil.copyfile(source, target)
        register(target, 'original_full_workflow_requires_external_inputs', name)
    metadata = {
        'curated_date': '2026-10-07', 'accepted_lineage': '2026-09-06 repaired conservation',
        'input_rows': len(selected), 'input_columns': list(selected), 'coordinate_system': 'EPSG:4326',
        'source_records': [{'name': f.name, 'sha256': sha(f), 'bytes': f.stat().st_size} for f in (cells_file, masks_file, source_script)],
        'row_selection': 'All 29,746 archived positive-area reef cells; no result-dependent filtering',
        'column_selection': 'Only numeric/boolean runtime fields; WKB geometry and unused polymer/size component fields omitted',
        'float_conversion': '17 significant digits, pandas round_trip; every source column tested for exact equality after read-back',
        'normalization': 'Retain the original global-reef min-max MPEI score; do not renormalize per reporting region',
        'source_hashes_before_and_after_unchanged': True, 'function_bodies_copied_unchanged': list(FUNCTIONS), 'files': files,
    }
    for source in (cells_file, masks_file, source_script):
        if sha(source) != source_hashes_before[source.name]:
            raise ValueError(f'Original source changed during curation: {source.name}')
    (MODULE / 'provenance.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(f'Curated {len(selected):,} real cells, {len(selected.columns)} fields; exact round-trip checked.')

if __name__ == '__main__':
    main()
