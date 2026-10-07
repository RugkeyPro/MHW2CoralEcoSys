"""Synthetic formula tests are separate from the actual-input MaxEnt demo."""
import unittest
import numpy as np
import pandas as pd
from analysis.modules.physiology_maxent_demo import core


class PhysiologyTests(unittest.TestCase):
    def test_coral_equation_at_optimum(self):
        self.assertAlmostEqual(float(core.coral_multiplier(27.0, 3.0)), np.sqrt(0.5))

    def test_hot_side_is_steeper(self):
        cold, hot = core.coral_multiplier([24.0, 30.0], 4.0)
        self.assertGreater(cold, hot)

    def test_missing_environment_is_not_filled(self):
        self.assertTrue(np.isnan(core.coral_multiplier(np.nan, 3.0)))

    def test_5degree_block_is_not_split(self):
        frame = pd.DataFrame({'longitude': [10.1, 10.9], 'latitude': [5.1, 5.9]})
        self.assertEqual(*core.split_spatial(frame))

    def test_rank_auc_with_ties(self):
        self.assertEqual(core.auc([1, 1], [0, 0]), 1.0)
        self.assertEqual(core.auc([0.5], [0.5]), 0.5)
