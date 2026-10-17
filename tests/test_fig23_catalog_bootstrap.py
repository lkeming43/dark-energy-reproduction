"""Checks for the separate Notebook 11b catalogue-level MC experiment."""

from unittest import TestCase

import numpy as np
from scipy.optimize import minimize_scalar

from src.fig23_catalog_bootstrap import WhitenedGLS


def synthetic_data():
    rng = np.random.default_rng(928)
    z = np.linspace(.03, 2.1, 18)
    sn_factor = rng.normal(size=(12, 12))
    sn_cov = .08 * np.eye(12) + .002 * sn_factor @ sn_factor.T
    covariance = np.zeros((18, 18))
    covariance[:6, :6] = np.diag(np.linspace(.03, .08, 6) ** 2)
    covariance[6:, 6:] = sn_cov
    y = 29.3 * z * np.exp(-.25 * z) + rng.normal(0, .04, len(z))
    return {"z": z, "y": y, "covariance": covariance,
            "source": np.array(["DESI_DR2"] * 6 + ["SN"] * 12)}


class CatalogBootstrapTests(TestCase):
    def test_quadratic_whitened_draw_matches_explicit_covariance_gls(self):
        data = synthetic_data()
        problem = WhitenedGLS(data, "quadratic")
        rd0, rd = 101.0, 104.0
        seed = 527
        fitted = problem.draw_and_fit(data["y"], np.array([rd]), rd0,
                                      np.random.default_rng(seed))[0]
        replay = np.random.default_rng(seed)
        mock = data["y"].copy()
        mock[:6] += problem.lower_bao @ replay.standard_normal(6)
        mock[6:] += problem.lower_sn @ replay.standard_normal(12)
        scale = rd0 / rd
        mock[6:] *= scale
        covariance = data["covariance"].copy()
        covariance[6:, 6:] *= scale ** 2
        design = np.column_stack((data["z"], data["z"] ** 2))
        inverse_design = np.linalg.solve(covariance, design)
        expected = np.linalg.solve(design.T @ inverse_design,
                                   inverse_design.T @ mock)
        residual = mock - design @ expected
        chi2 = residual @ np.linalg.solve(covariance, residual)
        np.testing.assert_allclose(fitted[:2], expected, atol=1e-9, rtol=0)
        np.testing.assert_allclose(fitted[2], chi2, atol=1e-7, rtol=0)

    def test_exponential_taylor_profile_matches_direct_gls(self):
        data = synthetic_data()
        problem = WhitenedGLS(data, "exponential")
        result = problem.fit_observed()
        inverse_y = np.linalg.solve(data["covariance"], data["y"])

        def direct(slope):
            x = data["z"] * np.exp(slope * data["z"])
            inverse_x = np.linalg.solve(data["covariance"], x)
            a = (x @ inverse_y) / (x @ inverse_x)
            residual = data["y"] - a * x
            return residual @ np.linalg.solve(data["covariance"], residual)

        optimum = minimize_scalar(direct, bounds=(-2, 0), method="bounded",
                                  options={"xatol": 1e-11})
        self.assertAlmostEqual(result[1], optimum.x, delta=1e-6)
        self.assertAlmostEqual(result[2], optimum.fun, delta=1e-6)
