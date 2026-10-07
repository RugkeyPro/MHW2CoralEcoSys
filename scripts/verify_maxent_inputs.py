"""Independently verify subset values, raster samples and unchanged source hashes."""
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
module = ROOT / 'analysis/modules/physiology_maxent_demo'
manifest = json.loads((module / 'provenance.json').read_text(encoding='utf-8'))
source_values = 0
for record in manifest['sources']:
    path = args.project / record['file']
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != record['sha256']:
        raise ValueError(f"Original source changed: {record['role']}")
    if record['role'] in ['presence', 'background']:
        subset = pd.read_csv(module / f"data/{record['role']}.csv.gz", float_precision='round_trip')
        columns = list(subset.columns.drop('parent_row'))
        parent = pd.read_csv(path, usecols=columns, float_precision='round_trip')
        np.testing.assert_array_equal(subset[columns].to_numpy(), parent.loc[subset.parent_row, columns].to_numpy())
        source_values += subset[columns].size
projection = pd.read_csv(module / 'data/projection.csv.gz', float_precision='round_trip')
raster_values = 0
for record in manifest['sources']:
    if ':' not in record['role']:
        continue
    scenario, variable = record['role'].split(':')
    frame = projection.loc[projection.scenario == scenario].iloc[::167]
    with rasterio.open(args.project / record['file']) as src:
        samples = np.ma.vstack(list(src.sample(zip(frame.longitude, frame.latitude), masked=True)))[:, 0].filled(np.nan)
    np.testing.assert_allclose(samples, frame[variable], rtol=0, atol=0, equal_nan=True)
    raster_values += len(frame)
report = {'status': 'passed', 'original_source_hashes_rechecked': len(manifest['sources']),
          'subset_numeric_values_exactly_matched': source_values,
          'independently_sampled_raster_values_matched': raster_values,
          'sampling_comparison': 'Independent rasterio.sample versus full-band indexing',
          'raw_inputs_modified': False}
(ROOT / 'docs/maxent_input_validation.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
print(json.dumps(report, indent=2))
