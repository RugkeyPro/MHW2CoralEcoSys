"""Independently reproject real raster-cell samples and audit full output masks/formulas."""
import json
from pathlib import Path
import sys
import argparse
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
import rasterio
from analysis.modules.physiology_maxent_demo import core, full_grid

ROOT = Path(__file__).resolve().parents[1]
module = ROOT / 'analysis/modules/physiology_maxent_demo'
parser = argparse.ArgumentParser()
parser.add_argument('--project', type=Path)
args = parser.parse_args()
original_hashes = 0
if args.project:
    source_manifest = json.loads((module / 'full_grid_sources.json').read_text(encoding='utf-8'))
    for record in source_manifest['files']:
        if full_grid.digest(args.project / record['original_source']) != record['sha256']:
            raise ValueError('An original full-grid source changed')
        original_hashes += 1
outputs = ROOT / 'outputs/physiology_maxent_demo'
audit_output = outputs / 'independent_full_grid_audit'
audit_output.mkdir(exist_ok=True)
records = []
for scenario in full_grid.SCENARIOS:
    profile, environment = full_grid.load_environment(module / 'data/full_grid', scenario)
    valid = np.logical_and.reduce([np.isfinite(environment[v]) for v in core.VARIABLES])
    folder = outputs / 'full_grid' / scenario
    with rasterio.open(folder / 'hsi_unconstrained.tif') as src:
        unconstrained = src.read(1, masked=True).filled(np.nan)
        assert src.transform == profile['transform'] and src.crs == profile['crs']
    with rasterio.open(folder / 'hsi_constrained.tif') as src:
        constrained = src.read(1, masked=True).filled(np.nan)
    np.testing.assert_array_equal(np.isfinite(unconstrained), valid)
    multiplier = core.coral_multiplier(environment['thetao'], environment['omega_arag']).astype('float32')
    np.testing.assert_array_equal(constrained, unconstrained * multiplier)
    indices = np.random.default_rng(20261008).choice(np.flatnonzero(valid), size=256, replace=False)
    rows, cols = np.unravel_index(indices, valid.shape)
    longitude, latitude = rasterio.transform.xy(profile['transform'], rows, cols)
    frame = pd.DataFrame({'longitude': longitude, 'latitude': latitude})
    for variable in core.VARIABLES:
        frame[variable] = environment[variable].ravel()[indices]
    _, native = core.project(module / 'vendor/maxent.jar', outputs / 'model/Acropora.lambdas',
                             frame, audit_output, scenario)
    columns = [c for c in native if c not in ['longitude', 'latitude', 'species']]
    actual = native[columns[0]].to_numpy().astype('float32')
    np.testing.assert_array_equal(actual, unconstrained.ravel()[indices])
    records.append({'scenario': scenario, 'entire_grid_cells_mask_checked': int(valid.size),
                    'valid_native_cells': int(valid.sum()), 'independent_native_reprojections': len(indices),
                    'maximum_sample_error': float(np.max(abs(actual - unconstrained.ravel()[indices]))),
                    'full_constrained_formula_and_missingness': 'exact float32 match'})
    print(f'{scenario}: independent native reprojection and full masks passed', flush=True)
report = {'status': 'passed', 'scope': 'Full grid masks/formula plus independent native reprojections',
          'spatial_subsampling_of_delivered_projection': False, 'audit_sample_seed': 20261008,
          'original_full_grid_source_hashes_rechecked': original_hashes,
          'scenarios': records}
(ROOT / 'docs/full_grid_projection_validation.json').write_text(json.dumps(report, indent=2) + '\n')
