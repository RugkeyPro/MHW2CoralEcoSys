"""Selection, physical input and scientific regression checks; fixtures are test-only."""
import argparse
import ast
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.modules.conservation_demo.code import conservation_core as core
from analysis.modules.conservation_demo.code.run_all import MODULE, REPO, load_inputs, compare_frames


class ConservationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reef, cls.masks = load_inputs(MODULE / 'data/reef_cells.csv.gz')

    def test_unmodified_scientific_function_bodies(self):
        original = ast.parse((REPO / 'research_workflow/process_future_mpei_conservation_20260903.py').read_text(encoding='utf-8-sig'))
        packaged = ast.parse((MODULE / 'code/conservation_core.py').read_text(encoding='utf-8'))
        def functions(tree):
            return {node.name: ast.dump(node, include_attributes=False) for node in tree.body if isinstance(node, ast.FunctionDef)}
        old, new = functions(original), functions(packaged)
        for name in new:
            self.assertEqual(new[name], old[name])

    def test_source_manifest_byte_hashes(self):
        manifest = json.loads((MODULE / 'provenance.json').read_text())
        for item in manifest['files']:
            self.assertEqual(hashlib.sha256((REPO / item['file']).read_bytes()).hexdigest(), item['sha256'])

    def test_real_input_coverage_and_masks(self):
        self.assertEqual(len(self.reef), 29746)
        self.assertTrue(self.reef.cell_id.is_unique)
        self.assertEqual(set(self.masks), set(core.SCENARIOS))
        for matrix in self.masks.values():
            self.assertEqual(matrix.shape, (3, 29746))

    def test_equal_score_fractional_tie_budget(self):
        # Explicit synthetic unit-test fixture, not scientific observations.
        area = np.array([2., 3.])
        selection = core.exact_budget_selection(np.array([1., 1.]), area, 2.5)
        np.testing.assert_array_equal(selection, [0.5, 0.5])
        self.assertEqual(float(np.sum(selection * area)), 2.5)

    def test_exposure_breaks_ties_without_changing_support_capture(self):
        # Synthetic unit-test fixture isolates the primary/secondary ranking rule.
        metrics, selections = core.calculate_metrics(np.ones(3), np.array([1., 0., 0.]),
            np.array([2/3, 2/3, 1/3]), np.array([0., 1., 1.]), np.array([1e-8, 2e-8, 3e-8]), np.ones(3, bool))
        np.testing.assert_array_equal(selections['climate_priority'], [0.5, 0.5, 0.])
        np.testing.assert_array_equal(selections['future_mpei_priority'], [0., 1., 0.])
        self.assertAlmostEqual(metrics['future_mpei_vs_climate_co_suitability_retained_fraction'], 1.)
        self.assertAlmostEqual(metrics['future_mpei_vs_climate_area_changed_fraction'], .5)

    def test_real_points_match_accepted_table(self):
        point, _ = core.point_estimates(self.reef, [])
        reference = pd.read_csv(MODULE / 'reference/point_estimates.csv', float_precision='round_trip')
        result = compare_frames(point, reference, ['region','scenario'])
        self.assertEqual(result['status'], 'passed')
        self.assertEqual(result['records'], 12)

    def test_invalid_physical_input_is_not_imputed(self):
        with tempfile.TemporaryDirectory() as directory:
            file = Path(directory) / 'invalid.csv'
            subset = self.reef.iloc[:10].copy()
            subset.loc[subset.index[0], 'reef_area_km2'] = -1
            subset.to_csv(file, index=False, float_format='%.17g')
            with self.assertRaisesRegex(ValueError, 'area'):
                load_inputs(file)


if __name__ == '__main__':
    unittest.main()
