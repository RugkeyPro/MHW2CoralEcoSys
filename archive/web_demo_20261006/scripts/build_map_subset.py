"""Aggregate valid original Figure 4 grid cells into one-degree display bins."""
import json
from pathlib import Path
import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
target = ROOT / 'public/data'
target.mkdir(parents=True, exist_ok=True)
report = {'display_resolution_degrees': 1, 'aggregation': 'Unweighted arithmetic mean of valid source cell values within each 1-degree display bin. Not used for regional headline statistics.', 'scenarios': {}}
for scenario in ('ssp126_2050', 'ssp245_2050', 'ssp585_2050'):
    name = f'{scenario}_ensemble_mean_quality_reduction_pct_20260902.tif'
    with rasterio.open(ROOT / 'data/source' / name) as dataset:
        if str(dataset.crs) != 'EPSG:4326':
            raise ValueError('Expected geographic source raster')
        grid = dataset.read(1, masked=True)
        row, col = np.where(~np.ma.getmaskarray(grid) & np.isfinite(grid.data))
        values = grid.data[row, col].astype(float)
        lon, lat = rasterio.transform.xy(dataset.transform, row, col)
        lon, lat = np.asarray(lon), np.asarray(lat)
        valid = (lon >= -180) & (lon < 180) & (lat >= -90) & (lat < 90)
        x, y = np.floor(lon[valid]).astype(int), np.floor(lat[valid]).astype(int)
        keys = (y + 90) * 360 + x + 180
        unique, inverse = np.unique(keys, return_inverse=True)
        sums = np.bincount(inverse, weights=values[valid])
        counts = np.bincount(inverse)
        cells = [[int(k % 360) - 179.5, int(k // 360) - 89.5, round(float(total / count), 5), int(count)] for k, total, count in zip(unique, sums, counts)]
        (target / f'map_{scenario}.json').write_text(json.dumps(cells, separators=(',', ':')) + '\n', encoding='utf-8')
        report['scenarios'][scenario] = {'source_file': name, 'original_shape': list(dataset.shape), 'source_valid_cells': int(len(values)), 'display_bins': len(cells), 'source_range_pct': [float(values.min()), float(values.max())], 'displayed_source_cells': int(counts.sum())}
(target / 'map_manifest.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
print(json.dumps(report, indent=2))
