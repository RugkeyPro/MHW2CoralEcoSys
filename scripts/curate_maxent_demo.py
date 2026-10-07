"""Curate a seeded actual-input subset without modifying parent research files."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument('--project', type=Path, required=True)
args = parser.parse_args()
project = args.project
module = ROOT / 'analysis/modules/physiology_maxent_demo'
data = module / 'data'
data.mkdir(parents=True, exist_ok=True)
sdm = project / '解压整理/01_ecoMarine/03_mainline/ecoMarine_0907_20260913_03_mainline_006/filesystem/D/ecoMarine/output/01_mainline/inputs/sdm_global_0083_corrected_20260722'
member = project / '解压整理/01_ecoMarine/02_data/ecoMarine_0907_20260913_02_data_004/filesystem/D/ecoMarine/data/processed/scenarios_shared_cmip6_ensemble/members/gfdl_esm4_r1i1p1f1_gr'
variables = (sdm / 'vif/Acropora_global/final_variables.txt').read_text().split()
columns = ['longitude', 'latitude'] + variables
sources = []

def register(path, role):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    sources.append({'file': path.relative_to(project).as_posix(), 'role': role,
                    'bytes': path.stat().st_size, 'sha256': digest.hexdigest()})

tables = {}
eligibility = {}
for label, count in [('presence', 1500), ('background', 10000)]:
    path = sdm / f'swd/{label}_raw_Acropora_global.csv'
    frame = pd.read_csv(path, usecols=columns, float_precision='round_trip')
    valid = frame.dropna(subset=columns).copy()
    eligibility[label] = {'original_rows': len(frame), 'complete_rows': len(valid),
                          'excluded_missing_environment': len(frame) - len(valid)}
    subset = valid.sample(n=min(count, len(valid)), random_state=20261007).sort_index()
    subset.insert(0, 'parent_row', subset.index)
    subset.to_csv(data / f'{label}.csv.gz', index=False, float_format='%.17g')
    tables[label] = subset
    register(path, label)

locations = tables['background'].sample(n=2000, random_state=20261008).copy()
locations.insert(0, 'location_id', np.arange(len(locations)))
projection = []
for scenario in ['baseline', 'ssp126_2050', 'ssp245_2050', 'ssp585_2050']:
    frame = locations[['location_id', 'parent_row', 'longitude', 'latitude', 'bathymetry', 'rugosity']].copy()
    frame.insert(0, 'scenario', scenario)
    for variable in ['bottomT', 'thetao', 'si', 'so', 'zos', 'omega_arag']:
        path = member / f'{scenario}_regional/{variable}.tif'
        with rasterio.open(path) as src:
            if src.crs.to_epsg() != 4326:
                raise ValueError(f'Unexpected CRS: {path}')
            rows, cols = rasterio.transform.rowcol(src.transform, frame.longitude, frame.latitude)
            band = src.read(1, masked=True)
            values = band[np.asarray(rows), np.asarray(cols)].filled(np.nan).astype(float)
            values[~np.isfinite(values)] = np.nan
            frame[variable] = values
        register(path, f'{scenario}:{variable}')
        print(f'Sampled and hashed {scenario}:{variable}', flush=True)
    frame['complete_environment'] = frame[variables + ['omega_arag']].notna().all(axis=1)
    projection.append(frame)
pd.concat(projection, ignore_index=True).to_csv(data / 'projection.csv.gz', index=False, float_format='%.17g')
parent_scripts = module / 'parent_sources'
parent_scripts.mkdir(exist_ok=True)
for filename in ['branch_04_maxent_formal.R', 'branch_05_maxent_projection.R']:
    path = sdm / 'scripts' / filename
    (parent_scripts / filename).write_bytes(path.read_bytes())
    register(path, 'parent_code')
authority = project / '解压整理/01_ecoMarine/03_mainline/ecoMarine_0907_20260913_03_mainline_009/filesystem/D/ecoMarine_worktrees/20260903_mainline_publish'
for path in [authority / 'scripts/prepare_figS1_sdm_support_data_20260806.py',
             authority / 'output/01_mainline/inputs/sdm_global_0083_corrected_20260722/maxent/optimization/Acropora_global/Acropora_best_model.csv']:
    (parent_scripts / path.name).write_bytes(path.read_bytes())
    register(path, 'parameter_authority')
manifest = {'curation_seed': 20261007, 'taxon': 'Acropora', 'variables': variables,
            'member': 'gfdl_esm4_r1i1p1f1_gr', 'presence_rows': len(tables['presence']),
            'background_rows': len(tables['background']), 'projection_locations': len(locations),
            'projection_sampling': 'Seeded global background locations; no outcome-dependent selection',
            'static_fields': 'bathymetry and rugosity copied from actual background SWD; other fields sampled from member rasters',
            'eligible_training_rows': eligibility, 'sources': sources, 'packaged_files': [],
            'native_engine': {'version': '3.4.4', 'upstream_commit': '963d0114e05a39f92f832dfbe800c22442d9d8a6',
                              'source': 'https://github.com/mrmaxent/Maxent/tree/963d0114e05a39f92f832dfbe800c22442d9d8a6/ArchivedReleases/3.4.4'}}
for path in sorted(list(data.glob('*')) + list(parent_scripts.glob('*')) + list((module / 'vendor').glob('*'))):
    manifest['packaged_files'].append({'file': path.relative_to(ROOT).as_posix(),
                                      'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
(module / 'provenance.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print(json.dumps({k: v for k, v in manifest.items() if k not in ['sources', 'packaged_files']}, indent=2))
