"""Full native-grid projection; missing source environments remain NoData."""
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window

from . import core

SCENARIOS = ['baseline', 'ssp126_2050', 'ssp245_2050', 'ssp585_2050']
NODATA = -9999.0


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def load_environment(folder, scenario):
    """Crop aligned static rasters to the complete native member grid, without interpolation."""
    with rasterio.open(folder / scenario / 'thetao.tif') as src:
        profile = src.profile.copy()
    arrays = {}
    for variable in core.VARIABLES + ['omega_arag']:
        path = folder / ('static' if variable in ['bathymetry', 'rugosity'] else scenario) / f'{variable}.tif'
        with rasterio.open(path) as src:
            if src.crs != profile['crs'] or src.transform != profile['transform']:
                raise ValueError(f'Grid alignment differs: {path.name}')
            if src.width != profile['width'] or src.height < profile['height']:
                raise ValueError(f'Insufficient source coverage: {path.name}')
            array = src.read(1, window=Window(0, 0, profile['width'], profile['height']), masked=True)
            arrays[variable] = array.filled(np.nan)
    return profile, arrays


def write_grid(path, values, profile):
    output_profile = profile.copy()
    output_profile.update(driver='GTiff', dtype='float32', count=1, nodata=NODATA,
                          compress='deflate', tiled=True, blockxsize=256, blockysize=256)
    with rasterio.open(path, 'w', **output_profile) as dst:
        dst.write(np.where(np.isfinite(values), values, NODATA).astype('float32'), 1)


def materialize_native_clamping(output):
    """Retain the engine's clamping diagnostic as a validated GeoTIFF, not a large temporary CSV."""
    source = output / 'native_predictions_clamping.csv'
    target = output / 'native_clamping.tif'
    if not source.is_file():
        if target.is_file():
            return
        raise ValueError('Native clamping diagnostic is missing')
    with rasterio.open(output / 'hsi_unconstrained.tif') as src:
        profile = src.profile.copy()
        valid = ~np.ma.getmaskarray(src.read(1, masked=True))
    indices = np.flatnonzero(valid)
    values = np.full(valid.shape, np.nan, dtype='float32')
    count = 0
    for frame in pd.read_csv(source, chunksize=100000):
        stop = count + len(frame)
        rows, cols = np.unravel_index(indices[count:stop], valid.shape)
        longitude, latitude = rasterio.transform.xy(profile['transform'], rows, cols)
        np.testing.assert_allclose(frame.longitude, longitude, rtol=0, atol=1e-5)
        np.testing.assert_allclose(frame.latitude, latitude, rtol=0, atol=1e-5)
        if not np.isfinite(frame.Clamping).all():
            raise ValueError('Non-finite native clamping values')
        values.ravel()[indices[count:stop]] = frame.Clamping.to_numpy()
        count = stop
    if count != len(indices):
        raise ValueError('Clamping coverage differs from HSI')
    write_grid(target, values, profile)
    with rasterio.open(target) as src:
        np.testing.assert_array_equal(src.read(1, masked=True).filled(np.nan), values)
    source.unlink()


def project_full_grid(jar, model, folder, output, scenario, training_environment):
    """Every valid grid cell is passed to native density.Project, not interpolated from samples."""
    started = time.perf_counter()
    output = output / scenario
    output.mkdir(parents=True, exist_ok=True)
    profile, arrays = load_environment(folder, scenario)
    shape = (profile['height'], profile['width'])
    valid = np.ones(shape, dtype=bool)
    for variable in core.VARIABLES:
        valid &= np.isfinite(arrays[variable])
    indices = np.flatnonzero(valid)
    rows, cols = np.unravel_index(indices, shape)
    longitude, latitude = rasterio.transform.xy(profile['transform'], rows, cols)
    environment_file = output / 'native_environment.csv'
    prediction_file = output / 'native_predictions.csv'
    # The intermediate SWD contains every complete native cell; it is not a sampled grid.
    with environment_file.open('w', encoding='utf-8', newline='') as stream:
        for start in range(0, len(indices), 100000):
            stop = min(start + 100000, len(indices))
            data = {'species': 'projection', 'longitude': np.asarray(longitude[start:stop]),
                    'latitude': np.asarray(latitude[start:stop])}
            data.update({variable: arrays[variable].ravel()[indices[start:stop]] for variable in core.VARIABLES})
            pd.DataFrame(data).to_csv(stream, index=False, header=start == 0, float_format='%.17g')
    print(f'{scenario}: submitting all {len(indices):,} complete grid cells to native MaxEnt', flush=True)
    core.run_java(jar, ['density.Project', str(model), str(environment_file), str(prediction_file),
                       'doclamp=true', 'extrapolate=false', 'fadebyclamping=true', 'outputformat=logistic'],
                  output / 'native_projection.log', timeout=1800)
    unconstrained = np.full(shape, np.nan, dtype='float32')
    count = 0
    for frame in pd.read_csv(prediction_file, chunksize=100000):
        prediction_columns = [column for column in frame.columns if column not in ['longitude', 'latitude', 'species']]
        if len(prediction_columns) != 1:
            raise ValueError('Unexpected native grid prediction schema')
        stop = count + len(frame)
        np.testing.assert_allclose(frame.longitude, longitude[count:stop], rtol=0, atol=1e-5)
        np.testing.assert_allclose(frame.latitude, latitude[count:stop], rtol=0, atol=1e-5)
        values = frame[prediction_columns[0]].to_numpy()
        if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
            raise ValueError('Invalid native HSI values')
        unconstrained.ravel()[indices[count:stop]] = values
        count = stop
    if count != len(indices):
        raise ValueError('Native prediction coverage is incomplete')
    multiplier = core.coral_multiplier(arrays['thetao'], arrays['omega_arag']).astype('float32')
    constrained = unconstrained * multiplier
    if np.any(constrained[np.isfinite(constrained)] > unconstrained[np.isfinite(constrained)] + 1e-7):
        raise ValueError('Physiological constraint increased coral HSI')
    novelty = np.zeros(shape, dtype='float32')
    for variable in core.VARIABLES:
        minimum, maximum = training_environment[variable].min(), training_environment[variable].max()
        novelty += (arrays[variable] < minimum) | (arrays[variable] > maximum)
    novelty[~valid] = np.nan
    for name, values in [('hsi_unconstrained', unconstrained), ('physiology_effective_multiplier', multiplier),
                         ('hsi_constrained', constrained), ('outside_training_range_variables', novelty)]:
        write_grid(output / f'{name}.tif', values, profile)
    complete = np.isfinite(constrained)
    # Spherical cell area at the source-grid latitude; avoids equating simple pixel means to area means.
    transform = profile['transform']
    latitude_edges = transform.f + np.arange(shape[0] + 1) * transform.e
    row_area = 6371.0088 ** 2 * np.deg2rad(abs(transform.a)) * np.abs(np.diff(np.sin(np.deg2rad(latitude_edges))))
    area = np.broadcast_to(row_area[:, None], shape)
    weights = area[complete]
    result = {'scenario': scenario, 'width': shape[1], 'height': shape[0],
              'total_grid_cells': int(np.prod(shape)), 'native_projected_cells': count,
              'physiology_complete_cells': int(complete.sum()),
              'missing_training_environment_cells': int((~valid).sum()),
              'missing_physiology_after_maxent_cells': int((valid & ~complete).sum()),
              'area_weighted_mean_unconstrained_on_common_domain': float(np.average(unconstrained[complete], weights=weights)),
              'area_weighted_mean_constrained': float(np.average(constrained[complete], weights=weights)),
              'outside_training_range_cells': int(((novelty > 0) & valid).sum()),
              'elapsed_seconds': round(time.perf_counter() - started, 3),
              'native_prediction_complete': True, 'spatial_subsampling': False,
              'source_predictor_precision_preserved': True,
              'maxent_projection_flags': ['doclamp=true', 'extrapolate=false', 'fadebyclamping=true', 'outputformat=logistic']}
    # Read each GeoTIFF fully and check values/NoData, not merely file existence or length.
    for name, expected in [('hsi_unconstrained', unconstrained), ('hsi_constrained', constrained)]:
        with rasterio.open(output / f'{name}.tif') as src:
            actual = src.read(1, masked=True).filled(np.nan)
            np.testing.assert_array_equal(actual, expected)
    (output / 'projection_checks.json').write_text(json.dumps(result, indent=2) + '\n')
    if not (output / 'native_predictions_clamping.csv').is_file():
        raise ValueError('The current native projection did not emit its clamping diagnostic')
    materialize_native_clamping(output)
    environment_file.unlink()
    prediction_file.unlink()
    print(f'{scenario}: verified full-grid outputs', flush=True)
    return result
