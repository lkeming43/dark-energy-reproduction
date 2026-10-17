"""Checks for error-version isolation and full-covariance whitening."""
import unittest
from pathlib import Path

import numpy as np

from src.data_loading import load_standard_dataset, with_desi_error_version


ROOT = Path(__file__).resolve().parents[1]


class PySRInputTests(unittest.TestCase):
    def test_error_switch_changes_only_six_desi_variances(self):
        for scenario, n in (("fig2_pantheon", 1596), ("fig3_desy5", 1760)):
            original = load_standard_dataset(
                ROOT / "results/08_fig23_data_preparation" / f"{scenario}.npz")
            public = with_desi_error_version(original, "public_release")
            paper = with_desi_error_version(original, "paper_table")
            self.assertEqual(len(paper["z"]), n)
            np.testing.assert_array_equal(public["z"], paper["z"])
            np.testing.assert_array_equal(public["y"], paper["y"])
            np.testing.assert_array_equal(public["covariance"][6:, 6:],
                                          paper["covariance"][6:, 6:])
            np.testing.assert_array_equal(original["covariance"], public["covariance"])
            np.testing.assert_allclose(np.sqrt(np.diag(paper["covariance"][:6, :6])),
                                       [.167, .177, .152, .318, .760, .531])
            self.assertEqual(original["metadata"].get("desi_error_version"), None)

    def test_column_major_whitening_reproduces_chi_square(self):
        cov = np.array([[4.0, 1.0], [1.0, 9.0]])
        residual = np.array([2.0, -3.0])
        w = np.linalg.inv(np.linalg.cholesky(cov))
        julia_reshape_equivalent = w.T.ravel().reshape(2, 2, order="F")
        np.testing.assert_array_equal(julia_reshape_equivalent, w)
        chi2 = float((w @ residual) @ (w @ residual))
        np.testing.assert_allclose(chi2, residual @ np.linalg.solve(cov, residual))


if __name__ == "__main__":
    unittest.main()
