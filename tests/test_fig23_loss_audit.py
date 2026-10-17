import unittest

import numpy as np

from src.loss_audit import (
    amplitude_loss,
    amplitude_quadratic,
    block_chi2,
    evaluate_formula_bank,
    magnitude_space_diagnostic,
)


class LossAuditTests(unittest.TestCase):
    def test_amplitude_quadratic_matches_direct_loss(self):
        z = np.array([0.1, 0.5, 1.0])
        y = np.array([2.0, 8.0, 14.0])
        cov = np.array([[0.4, 0.03, 0.0], [0.03, 0.6, 0.02], [0.0, 0.02, 0.8]])
        dataset = {"z": z, "y": y}
        shape = z * 0.76 ** z
        quadratic = amplitude_quadratic(dataset, cov, shape)
        for amplitude in (25.0, 30.0, 33.0):
            residual = amplitude * shape - y
            direct = residual @ np.linalg.solve(cov, residual) / len(z)
            self.assertAlmostEqual(float(amplitude_loss(quadratic, amplitude)), direct)
        self.assertLessEqual(
            float(amplitude_loss(quadratic, quadratic["best_amplitude"])),
            float(amplitude_loss(quadratic, 30.0)),
        )

    def test_formula_bank_and_block_decomposition(self):
        z = np.arange(1.0, 9.0)
        y = 2 * z + 0.2
        cov = np.diag(np.linspace(0.1, 0.8, 8))
        cov[6, 7] = cov[7, 6] = 0.03
        dataset = {"z": z, "y": y}
        formulas = [{"id": "linear", "expression": "2*z", "paper_location": "test",
                     "paper_reported_loss": 0.0}]
        result = evaluate_formula_bank(dataset, formulas, cov)[0]
        prediction = 2 * z
        parts = block_chi2(dataset, prediction, cov)
        self.assertAlmostEqual(result["our_chi2"], parts["total_chi2"])
        self.assertAlmostEqual(parts["total_loss"], result["our_loss"])
        self.assertGreater(parts["sn_chi2"], 0)

    def test_magnitude_space_zero_residual(self):
        z = np.linspace(0.1, 0.7, 7)
        prediction = 20 * z
        calibration = {"M": -19.253, "rd_mpc": 138.28}
        observable = 25 + 5 * np.log10(prediction[6:] * (1 + z[6:])
                                         * calibration["rd_mpc"]) + calibration["M"]
        dataset = {"z": z, "sn_observable": observable,
                   "sn_observable_covariance": np.array([[0.1]])}
        result = magnitude_space_diagnostic(dataset, prediction, calibration)
        self.assertAlmostEqual(result["sn_observable_chi2"], 0.0)
        self.assertAlmostEqual(result["profiled_common_offset_mag"], 0.0)


if __name__ == "__main__":
    unittest.main()
