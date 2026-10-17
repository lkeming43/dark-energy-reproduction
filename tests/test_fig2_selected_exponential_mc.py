"""Numerical and interpretation checks for Notebook 11d."""

import unittest

import numpy as np
import pandas as pd

from src.fig2_selected_exponential_mc import _outside_segments, _raw_curve


class SelectedExponentialMCTests(unittest.TestCase):
    def test_raw_exponential_distance_and_w_are_finite(self):
        grid = np.linspace(0, 2.33, 31)
        frame = _raw_curve(None, 29.321836856709833, -0.2502321164835667,
                           0.3, 101.0, grid)
        np.testing.assert_allclose(frame.distance,
                                   29.321836856709833 * grid *
                                   np.exp(-0.2502321164835667 * grid))
        self.assertTrue(np.isfinite(frame.w).all())
        self.assertEqual(frame.distance.iloc[0], 0.0)

    def test_outside_segments_have_grid_endpoints(self):
        rows = pd.DataFrame({"center": ["raw"] * 4,
                             "selection": ["pole_free"] * 4,
                             "quantity": ["w"] * 4,
                             "z": [0., .1, .2, .3],
                             "q0p135": [-.9, -1.1, -1.2, -.9],
                             "q99p865": [-.8, -.8, -.7, -.8]})
        self.assertEqual(_outside_segments(rows, "raw", "pole_free"),
                         [[0.0, 0.0], [0.3, 0.3]])


if __name__ == "__main__":
    unittest.main()
