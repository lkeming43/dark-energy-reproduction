"""Mathematical and selection checks for Notebook 11 (no PySR needed)."""

from unittest import TestCase

import numpy as np
import pandas as pd

from src.fig23_reconstruction import (
    C_KM_S, audit_and_select, candidate_diagnostics, gamma_value,
    lcdm_distance, w_from_derivatives,
)


class Figure23ReconstructionTests(TestCase):
    def test_lcdm_analytic_derivatives_give_minus_one(self):
        z = np.linspace(0, 2.33, 100)
        omega_m, rd = 0.3, 101.0
        e2 = omega_m * (1 + z) ** 3 + 1 - omega_m
        amplitude = C_KM_S / (100 * rd)
        d1 = amplitude / np.sqrt(e2)
        d2 = -1.5 * amplitude * omega_m * (1 + z) ** 2 / e2 ** 1.5
        w, factor = w_from_derivatives(z, d1, d2, gamma_value(omega_m, rd))
        np.testing.assert_allclose(w, -1, atol=2e-14)
        self.assertTrue((factor > 0).all())

    def test_lcdm_distance_integral_and_origin(self):
        z = np.array([0.0, .1, .5, 1.0])
        result = lcdm_distance(z, .3, 101.0)
        self.assertEqual(result[0], 0.0)
        self.assertTrue((np.diff(result) > 0).all())

    def test_exclude_nonphysical_and_select_without_w_preference(self):
        grid = np.linspace(0, 2.33, 301)
        self.assertFalse(candidate_diagnostics("30 + z", grid)["eligible"])
        rows = pd.DataFrame([
            {"scenario": "fig2_pantheon", "desi_error_version": "public_release",
             "seed": 1, "complexity": 5, "python_loss": 1.0,
             "sympy_expression": "30*z*exp(-0.25*z)"},
            {"scenario": "fig2_pantheon", "desi_error_version": "public_release",
             "seed": 2, "complexity": 4, "python_loss": 1.2,
             "sympy_expression": "z*(29-5.7*z)"},
        ])
        audit, chosen, decision = audit_and_select(rows, "fig2_pantheon", grid)
        self.assertEqual(int(chosen.seed), 2)
        self.assertEqual(decision["eligible_count"], 2)
        self.assertEqual(int(audit.independent_primary.sum()), 1)
