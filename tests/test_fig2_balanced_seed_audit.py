"""Selection-rule checks for Notebook 11c."""

import unittest

import pandas as pd

from src.fig2_balanced_seed_audit import select_balanced


class BalancedSeedAuditTests(unittest.TestCase):
    def test_worst_subset_beats_lower_complexity(self):
        scores = pd.DataFrame([
            {"seed": 1, "eligible": True, "within_notebook11_loss_gate": True,
             "worst_subset_chi2_per_point": 7.0,
             "combined_chi2_per_point": .92, "complexity": 7},
            {"seed": 2, "eligible": True, "within_notebook11_loss_gate": True,
             "worst_subset_chi2_per_point": 1.7,
             "combined_chi2_per_point": .88, "complexity": 9},
            {"seed": 3, "eligible": False, "within_notebook11_loss_gate": True,
             "worst_subset_chi2_per_point": .5,
             "combined_chi2_per_point": .7, "complexity": 5},
        ])
        self.assertEqual(int(select_balanced(scores).seed), 2)

    def test_loss_gate_excludes_candidate(self):
        scores = pd.DataFrame([
            {"seed": 1, "eligible": True, "within_notebook11_loss_gate": False,
             "worst_subset_chi2_per_point": .1,
             "combined_chi2_per_point": 2.0, "complexity": 5},
            {"seed": 2, "eligible": True, "within_notebook11_loss_gate": True,
             "worst_subset_chi2_per_point": 1.7,
             "combined_chi2_per_point": .88, "complexity": 9},
        ])
        self.assertEqual(int(select_balanced(scores).seed), 2)


if __name__ == "__main__":
    unittest.main()
