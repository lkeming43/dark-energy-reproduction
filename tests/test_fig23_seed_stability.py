"""Fast tests for the resumable Figure 2/3 seed ensemble (no Julia required)."""

from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase
from unittest.mock import patch

import numpy as np
import pandas as pd

from src.seed_stability import (
    classify_expression, collect_runs, run_seed_ensemble, save_summary_plots,
    setup_experiment, summarize_runs,
)


class SeedStabilityTests(TestCase):
    def test_expression_families(self):
        exponential = classify_expression("30.1*z/1.31**z")
        self.assertEqual(exponential["family"], "a*z*exp(b*z)")
        self.assertAlmostEqual(exponential["exp_a"], 30.1)
        self.assertAlmostEqual(exponential["exp_b"], -np.log(1.31))
        self.assertEqual(classify_expression("z*(29-5*z)")["family"],
                         "polynomial_degree_2")
        self.assertEqual(classify_expression("sqrt(z)-z**3+24*z")["family"],
                         "fractional_power_or_radical")

    def test_resume_extension_and_changed_inputs_rejected(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            prepared = root / "data.npz"
            prepared.write_bytes(b"fixed data")
            kwargs = dict(output_root=root / "results", config={"search": {"niterations": 30}},
                          prepared_paths={"fig2": prepared}, scenarios=["fig2"],
                          error_versions=["public_release"], paper_expressions={"fig2": "z"})
            directory = setup_experiment(**kwargs, seeds=[1, 2])
            self.assertEqual(len(collect_runs(directory)[0]), 2)
            self.assertEqual(setup_experiment(**kwargs, seeds=[1, 2, 3],
                                              resume_directory=directory), directory)
            self.assertEqual(len(collect_runs(directory)[0]), 3)
            with self.assertRaises(ValueError):
                setup_experiment(**{**kwargs, "config": {"search": {"niterations": 31}}},
                                 seeds=[1, 2, 3], resume_directory=directory)
            prepared.write_bytes(b"changed data")
            with self.assertRaises(ValueError):
                setup_experiment(**kwargs, seeds=[1, 2, 3], resume_directory=directory)

    def test_partial_attempt_is_preserved_and_completed_seed_is_skipped(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            prepared = root / "data.npz"
            prepared.write_bytes(b"fixed data")
            config = {"search": {"niterations": 1}}
            paper = {"fig2": "z"}
            directory = setup_experiment(root / "results", config, {"fig2": prepared},
                                         ["fig2"], ["public_release"], [7], paper)
            data = {"z": np.array([.5, 1.0]), "y": np.array([.5, 1.0]),
                    "covariance": np.eye(2), "metadata": {"scenario": "fig2",
                                                     "desi_error_version": "public_release"}}
            variants = {("fig2", "public_release"): data}

            def fake_fit(dataset, fit_config, seed, path):
                path.mkdir(parents=True, exist_ok=False)
                return {"summary": {"scenario": "fig2",
                                    "desi_error_version": "public_release", "seed": seed,
                                    "sympy_expression": "2*z*exp(-0.2*z)",
                                    "python_loss": .9, "elapsed_seconds": .01,
                                    "finite_on_grid": True,
                                    "positive_for_z_gt_zero": True,
                                    "satisfies_dm_zero_constraint": True},
                        "z_grid": np.array([0.0, 1.0]),
                        "selected_grid": np.array([0.0, 2*np.exp(-.2)])}

            def fake_export(result, path, *_args):
                pd.DataFrame({"z": result["z_grid"],
                              "selected": result["selected_grid"]}).to_csv(
                                  path / "distance_grid.csv", index=False)

            with patch("src.seed_stability.fit_exponential_reference",
                       return_value={"a": 1.0, "b": 0.0, "loss": 0.0}), \
                 patch("src.seed_stability.fit_symbolic_distance",
                       side_effect=RuntimeError("simulated failure")):
                with self.assertRaises(RuntimeError):
                    run_seed_ensemble(directory, variants, config, paper)
            manifest, _ = collect_runs(directory)
            self.assertEqual(manifest.status.iloc[0], "partial")
            first = directory / "runs/fig2/public_release/seed_7/attempt_001"
            self.assertTrue((first / "failure.json").is_file())

            with patch("src.seed_stability.fit_exponential_reference",
                       return_value={"a": 1.0, "b": 0.0, "loss": 0.0}), \
                 patch("src.seed_stability.fit_symbolic_distance", side_effect=fake_fit) as fit, \
                 patch("src.seed_stability.export_fit_result", side_effect=fake_export):
                run_seed_ensemble(directory, variants, config, paper)
                fit.assert_called_once()
                run_seed_ensemble(directory, variants, config, paper)
                fit.assert_called_once()
            second = directory / "runs/fig2/public_release/seed_7/attempt_002"
            self.assertTrue((second / "run_complete.json").is_file())
            self.assertTrue((first / "failure.json").is_file())
            summary, families, curves, selected = summarize_runs(directory)
            self.assertEqual(summary.n_completed.iloc[0], 1)
            self.assertEqual(families.family.iloc[0], "a*z*exp(b*z)")
            self.assertEqual(len(curves), 2)
            self.assertEqual(len(selected), 1)
            loss, distance = save_summary_plots(directory, paper)
            self.assertTrue(loss.is_file())
            self.assertTrue(distance.is_file())
