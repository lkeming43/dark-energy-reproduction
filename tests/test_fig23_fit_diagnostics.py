"""Small independent numerical checks for Notebook 12 diagnostics."""

import unittest

import numpy as np

from src.fig23_fit_diagnostics import covariance_diagnostics, fixed_prediction


class FitDiagnosticsTests(unittest.TestCase):
    def test_full_covariance_chi2_and_whitening(self):
        observed = np.array([1.0, 3.0])
        predicted = np.array([0.5, 2.5])
        covariance = np.array([[2.0, 0.8], [0.8, 1.0]])
        result = covariance_diagnostics(observed, predicted, covariance)
        residual = observed - predicted
        expected = residual @ np.linalg.solve(covariance, residual)
        self.assertAlmostEqual(result["chi2"], expected)
        self.assertAlmostEqual(result["chi2"], result["whitened"] @ result["whitened"])
        np.testing.assert_allclose(result["marginal_sigma"],
                                   residual / np.sqrt(np.diag(covariance)))

    def test_subset_chi2_adds_for_block_diagonal_matrix(self):
        covariance = np.array([[2.0, 0.0, 0.0],
                               [0.0, 1.0, 0.4],
                               [0.0, 0.4, 2.0]])
        observed = np.array([1.0, 2.0, 3.0])
        predicted = np.zeros(3)
        joint = covariance_diagnostics(observed, predicted, covariance)["chi2"]
        first = covariance_diagnostics(observed[:1], predicted[:1], covariance[:1, :1])["chi2"]
        second = covariance_diagnostics(observed[1:], predicted[1:], covariance[1:, 1:])["chi2"]
        self.assertAlmostEqual(joint, first + second)

    def test_frozen_prediction_has_no_refit(self):
        record = {"printed_paper_expression": "33.3*z*exp(-0.274*z)",
                  "branches": {"independent": {"refitted_template": "z*(p0*z + p1)",
                                               "refitted_parameters": [-5.7, 29.0]}}}
        z = np.array([0.0, 1.0, 2.0])
        np.testing.assert_allclose(fixed_prediction(z, "independent", record),
                                   z * (29.0 - 5.7 * z))
        np.testing.assert_allclose(fixed_prediction(z, "printed_paper", record),
                                   33.3 * z * np.exp(-0.274 * z))

    def test_invalid_covariance_fails(self):
        with self.assertRaises(np.linalg.LinAlgError):
            covariance_diagnostics([1.0], [0.0], [[-1.0]])


if __name__ == "__main__":
    unittest.main()
