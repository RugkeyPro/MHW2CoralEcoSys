"""Meaningful source and numerical-reproduction checks for the selected dataset."""
import json
import unittest
from pathlib import Path
from analysis.run_demo import assemble, close

ROOT = Path(__file__).resolve().parents[1]

class DemoTests(unittest.TestCase):
    def test_reproduction(self):
        data = assemble()
        self.assertEqual(data['checks']['numericComparisons'], 252)
        row = next(r for r in data['future'] if r['scenario'] == 'ssp585_2050' and r['region'] == 'Global')
        self.assertAlmostEqual(row['pressure']['mean'], 47.870631310467964, places=10)

    def test_frozen_browser_data_matches_recalculation(self):
        published = json.loads((ROOT / 'public/data/demo.json').read_text(encoding='utf-8'))
        self.assertEqual(published, assemble())

    def test_no_duplicate_historical_keys(self):
        rows = assemble()['historical']
        self.assertEqual(len({(r['region'],r['taxon'],r['group']) for r in rows}), 60)

    def test_invalid_numbers_fail_loudly(self):
        with self.assertRaises(ValueError):
            close(float('nan'), 0, 'Missing input')
        with self.assertRaises(ValueError):
            close(12, 3, 'Changed source')

    def test_map_samples_preserve_source_coverage(self):
        manifest = json.loads((ROOT / 'public/data/map_manifest.json').read_text())
        for scenario, metadata in manifest['scenarios'].items():
            points = json.loads((ROOT / f'public/data/map_{scenario}.json').read_text())
            self.assertEqual(sum(p[3] for p in points), metadata['source_valid_cells'])
            self.assertEqual(len(points), 2479)
            self.assertTrue(any(p[2] < 0 for p in points))
            self.assertTrue(all(-180 <= p[0] <= 180 and -90 <= p[1] <= 90 for p in points))

if __name__ == '__main__':
    unittest.main()
