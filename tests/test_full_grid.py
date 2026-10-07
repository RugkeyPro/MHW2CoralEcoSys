"""Explicitly synthetic raster fixtures test alignment and preservation of missing inputs."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
import rasterio
from rasterio.transform import from_origin
from analysis.modules.physiology_maxent_demo import core, full_grid


class FullGridTests(unittest.TestCase):
    def make_sources(self, directory):
        root = Path(directory)
        for variable in core.VARIABLES + ['omega_arag']:
            static = variable in ['bathymetry', 'rugosity']
            folder = root / ('static' if static else 'baseline')
            folder.mkdir(exist_ok=True)
            values = np.ones((3 if static else 2, 2), dtype='float32')
            if variable == 'si':
                values[0, 0] = -9999
            with rasterio.open(folder / f'{variable}.tif', 'w', driver='GTiff', width=2,
                               height=values.shape[0], count=1, dtype='float32', crs='EPSG:4326',
                               transform=from_origin(-180, 80, 0.083, 0.083), nodata=-9999) as dst:
                dst.write(values, 1)
        return root

    def test_aligned_static_crop_does_not_resample(self):
        with tempfile.TemporaryDirectory() as folder:
            root = self.make_sources(folder)
            profile, arrays = full_grid.load_environment(root, 'baseline')
            self.assertEqual((profile['height'], profile['width']), (2, 2))
            np.testing.assert_array_equal(arrays['bathymetry'], np.ones((2, 2)))
            self.assertTrue(np.isnan(arrays['si'][0, 0]))

    def test_geo_alignment_difference_is_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root = self.make_sources(folder)
            with rasterio.open(root / 'static/rugosity.tif', 'r+') as dst:
                dst.transform = from_origin(-179.9, 80, 0.083, 0.083)
            with self.assertRaises(ValueError):
                full_grid.load_environment(root, 'baseline')

    def test_source_float64_precision_is_preserved(self):
        with tempfile.TemporaryDirectory() as folder:
            root = self.make_sources(folder)
            value = 1.0 + 2.0 ** -30
            with rasterio.open(root / 'static/bathymetry.tif', 'w', driver='GTiff', width=2,
                               height=3, count=1, dtype='float64', crs='EPSG:4326',
                               transform=from_origin(-180, 80, 0.083, 0.083), nodata=-9999) as dst:
                dst.write(np.full((3, 2), value, dtype='float64'), 1)
            _, arrays = full_grid.load_environment(root, 'baseline')
            self.assertEqual(arrays['bathymetry'].dtype, np.dtype('float64'))
            self.assertEqual(arrays['bathymetry'][0, 0], value)
