"""Exact data-level Gaussian bootstrap for the fixed Figure 2/3 model families.

Gaussian catalogue draws are represented in whitened coordinates.  Projecting
those draws onto a model basis gives exactly the same GLS fit as materializing
each dense, correlated mock catalogue, while avoiding enormous mock arrays.
The exponential basis uses a verified high-order Taylor representation.
"""

from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4
import hashlib
import json
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular
from scipy.optimize import minimize_scalar

from .data_loading import load_standard_dataset
from .fig23_reconstruction import (
    SCENARIOS, _sha256, compile_curve, gamma_value, lcdm_distance,
    w_from_derivatives,
)


CENTERS = ("best_fit_expression", "observed_data")
QUANTILES = (.001349898, .022750132, .025, .158655254, .5,
             .841344746, .975, .977249868, .998650102)
QUANTILE_NAMES = ("q0p135", "q2p275", "q025", "q15p866", "median",
                  "q84p134", "q975", "q97p725", "q99p865")
SLOPE_BOUNDS = (-2.0, 0.0)  # Same optimization bounds as Notebook 07.
EXP_ORDER = 26


def _write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")
    temporary.replace(path)


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True,
                                  separators=(",", ":")).encode()).hexdigest()


def _anti_diagonal_sum(gram):
    result = np.zeros(gram.shape[0] * 2 - 1)
    for row in range(gram.shape[0]):
        for col in range(gram.shape[1]):
            result[row + col] += gram[row, col]
    return result


def _exp_taylor_columns(z, order=EXP_ORDER):
    base = z * np.exp(-z)
    columns = [base]
    for degree in range(1, order + 1):
        columns.append(columns[-1] * z / degree)
    return np.column_stack(columns)


class WhitenedGLS:
    """Block-GLS sufficient statistics for exact quadratic or exponential refits."""

    def __init__(self, data, family):
        if family not in ("quadratic", "exponential"):
            raise ValueError(f"Unsupported family: {family}")
        self.family = family
        self.z = np.asarray(data["z"], dtype=float)
        self.y = np.asarray(data["y"], dtype=float)
        source = np.asarray(data["source"])
        self.bao = np.flatnonzero(source == "DESI_DR2")
        self.sn = np.flatnonzero(source != "DESI_DR2")
        covariance = np.asarray(data["covariance"], dtype=float)
        if len(self.bao) != 6 or not np.array_equal(self.bao, np.arange(6)):
            raise ValueError("Expected the six DESI measurements first")
        if np.any(covariance[np.ix_(self.bao, self.sn)] != 0):
            raise ValueError("The factorized bootstrap requires zero BAO/SN cross-covariance")
        self.lower_bao = np.linalg.cholesky(covariance[np.ix_(self.bao, self.bao)])
        self.lower_sn = np.linalg.cholesky(covariance[np.ix_(self.sn, self.sn)])
        if family == "quadratic":
            design = np.column_stack((self.z, self.z ** 2))
        else:
            design = _exp_taylor_columns(self.z)
        self.v_bao = solve_triangular(self.lower_bao, design[self.bao],
                                      lower=True, check_finite=False)
        self.v_sn = solve_triangular(self.lower_sn, design[self.sn],
                                     lower=True, check_finite=False)
        self.gram_bao = self.v_bao.T @ self.v_bao
        self.gram_sn = self.v_sn.T @ self.v_sn
        self.observed_u_bao = self.whiten(self.y[self.bao], "bao")
        self.observed_u_sn = self.whiten(self.y[self.sn], "sn")
        if family == "exponential":
            self.den_bao = _anti_diagonal_sum(self.gram_bao)
            self.den_sn = _anti_diagonal_sum(self.gram_sn)
            self._check_taylor_accuracy()

    def whiten(self, values, source):
        lower = self.lower_bao if source == "bao" else self.lower_sn
        return solve_triangular(lower, values, lower=True, check_finite=False)

    def predict(self, theta, z=None):
        z = self.z if z is None else np.asarray(z, dtype=float)
        a, b = np.asarray(theta, dtype=float)
        if self.family == "quadratic":
            return a * z + b * z ** 2
        return a * z * np.exp(b * z)

    def _check_taylor_accuracy(self):
        # Direct whitening is used only for validation, never for every MC fit.
        for slope in (-2.0, -1.25, -.5, -.26, 0.0):
            powers = np.power(slope + 1.0, np.arange(EXP_ORDER + 1))
            for indices, lower, basis in ((self.bao, self.lower_bao, self.v_bao),
                                          (self.sn, self.lower_sn, self.v_sn)):
                direct = solve_triangular(lower,
                                          self.z[indices] * np.exp(slope * self.z[indices]),
                                          lower=True, check_finite=False)
                approximate = basis @ powers
                np.testing.assert_allclose(approximate, direct, rtol=2e-10, atol=2e-10)

    def _exponential_fit_one(self, p_bao, p_sn, t, norm_squared):
        numerator_coeff = p_bao + t * p_sn
        denominator_coeff = self.den_bao + t * t * self.den_sn

        def minus_gain(slope):
            delta = slope + 1.0
            numerator = np.polynomial.polynomial.polyval(delta, numerator_coeff)
            denominator = np.polynomial.polynomial.polyval(delta, denominator_coeff)
            return -numerator * numerator / denominator if denominator > 0 else np.inf

        optimum = minimize_scalar(minus_gain, bounds=SLOPE_BOUNDS,
                                  method="bounded", options={"xatol": 1e-11,
                                                             "maxiter": 200})
        if not optimum.success:
            raise RuntimeError("Exponential profile fit failed")
        slope = float(optimum.x)
        delta = slope + 1.0
        numerator = np.polynomial.polynomial.polyval(delta, numerator_coeff)
        denominator = np.polynomial.polynomial.polyval(delta, denominator_coeff)
        amplitude = numerator / denominator
        return amplitude, slope, norm_squared - numerator * numerator / denominator

    def fit_observed(self):
        p_bao = self.observed_u_bao @ self.v_bao
        p_sn = self.observed_u_sn @ self.v_sn
        norm = self.observed_u_bao @ self.observed_u_bao + self.observed_u_sn @ self.observed_u_sn
        if self.family == "exponential":
            return np.asarray(self._exponential_fit_one(p_bao, p_sn, 1.0, norm))
        gram = self.gram_bao + self.gram_sn
        theta = np.linalg.solve(gram, p_bao + p_sn)
        chi2 = norm - theta @ (p_bao + p_sn)
        return np.array([theta[0], theta[1], chi2])

    def draw_and_fit(self, center, rd, rd0, rng):
        """One batch of Gaussian mock catalogues and exact fixed-family GLS fits."""
        center = np.asarray(center, dtype=float)
        batch = len(rd)
        u_bao = self.whiten(center[self.bao], "bao")[None, :] + rng.standard_normal((batch, len(self.bao)))
        u_sn = self.whiten(center[self.sn], "sn")[None, :] + rng.standard_normal((batch, len(self.sn)))
        p_bao = u_bao @ self.v_bao
        p_sn = u_sn @ self.v_sn
        norm = np.sum(u_bao * u_bao, axis=1) + np.sum(u_sn * u_sn, axis=1)
        t = rd / rd0
        if self.family == "quadratic":
            gram = self.gram_bao[None, :, :] + t[:, None, None] ** 2 * self.gram_sn[None, :, :]
            p = p_bao + t[:, None] * p_sn
            theta = np.linalg.solve(gram, p[:, :, None])[:, :, 0]
            chi2 = norm - np.sum(theta * p, axis=1)
            return np.column_stack((theta, chi2))
        fitted = np.empty((batch, 3))
        for index in range(batch):
            fitted[index] = self._exponential_fit_one(p_bao[index], p_sn[index],
                                                       t[index], norm[index])
        return fitted


def reconstruct_samples(family, amplitude, slope, rd, z, omega_m=.3):
    z2 = z[None, :]
    a = amplitude[:, None]
    b = slope[:, None]
    if family == "quadratic":
        dm = a * z2 + b * z2 ** 2
        d1 = a + 2 * b * z2
        d2 = np.broadcast_to(2 * b, dm.shape)
    else:
        dm = a * z2 * np.exp(b * z2)
        d1 = a * np.exp(b * z2) * (1 + b * z2)
        d2 = a * np.exp(b * z2) * (2 * b + b * b * z2)
    w, factor = w_from_derivatives(z2, d1, d2,
                                   gamma_value(omega_m, rd)[:, None])
    finite = (np.isfinite(dm).all(axis=1) & np.isfinite(d1).all(axis=1)
              & np.isfinite(d2).all(axis=1) & np.isfinite(w).all(axis=1)
              & np.isfinite(rd) & (rd > 0) & np.isfinite(amplitude) & (amplitude > 0))
    pole = np.any(factor[:, :-1] * factor[:, 1:] <= 0, axis=1)
    pole |= np.any(factor == 0, axis=1)
    monotonic = (d1 > 0).all(axis=1)
    return dm, w, finite, pole, monotonic


def _bands(z, values, quantity, center, selection):
    frame = pd.DataFrame({"z": z, "quantity": quantity, "center": center,
                          "selection": selection,
                          "n_finite": np.isfinite(values).sum(axis=0)})
    for quantile, name in zip(QUANTILES, QUANTILE_NAMES):
        frame[name] = np.nanquantile(values, quantile, axis=0)
    return frame


def _figure(path, data, grid, family, center, bands, paper_expression, omega_m, rd0):
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), sharex="col")
    records = bands.loc[bands.center == center]
    paper_dm, paper_d1, paper_d2 = compile_curve(paper_expression)(grid)
    paper_w, _ = w_from_derivatives(grid, paper_d1, paper_d2,
                                    gamma_value(omega_m, rd0))
    for row, selection in enumerate(("all", "pole_free")):
        ax_d, ax_w = axes[row]
        distance = records.loc[(records.selection == selection) & (records.quantity == "distance")]
        w = records.loc[(records.selection == selection) & (records.quantity == "w")]
        z = distance.z.to_numpy()
        ax_d.fill_between(z, distance.q0p135, distance.q99p865, color="tab:blue", alpha=.12, label="3 sigma")
        ax_d.fill_between(z, distance.q2p275, distance.q97p725, color="tab:blue", alpha=.2, label="2 sigma")
        ax_d.fill_between(z, distance.q15p866, distance.q84p134, color="tab:blue", alpha=.3, label="1 sigma")
        ax_d.plot(z, distance["median"], color="tab:blue", label="Bootstrap median")
        desi = data["source"] == "DESI_DR2"
        ax_d.scatter(data["z"][~desi], data["y"][~desi], s=3, alpha=.06,
                     color="gray", label="Supernovae")
        ax_d.errorbar(data["z"][desi], data["y"][desi], yerr=data["sigma"][desi],
                      fmt="o", color="black", markersize=3, label="DESI DR2")
        ax_d.plot(z, lcdm_distance(z, omega_m, rd0), "--", color="darkgreen", label="LambdaCDM")
        ax_d.plot(z, paper_dm, ":", color="darkorange", label="Printed paper expression")
        ax_d.set(ylabel=f"{selection}: $D_M/r_d$", title=f"{family} · {center}")
        ax_w.fill_between(z, w.q0p135, w.q99p865, color="tab:blue", alpha=.12)
        ax_w.fill_between(z, w.q2p275, w.q97p725, color="tab:blue", alpha=.2)
        ax_w.fill_between(z, w.q15p866, w.q84p134, color="tab:blue", alpha=.3)
        ax_w.plot(z, w["median"], color="tab:blue", label="Bootstrap median")
        ax_w.axhline(-1, ls="--", color="darkgreen", label="w=-1")
        ax_w.plot(z, paper_w, ":", color="darkorange", label="Printed paper expression")
        ax_w.set(ylabel=f"{selection}: w(z)", ylim=(-4, 2), title="Pointwise intervals")
        ax_d.legend(fontsize=7, loc="best")
        ax_w.legend(fontsize=7, loc="best")
    axes[-1, 0].set_xlabel("Redshift z")
    axes[-1, 1].set_xlabel("Redshift z")
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def _compare_linearized(bootstrap_bands, linearized_file, scenario, branch):
    prior = pd.read_csv(linearized_file)
    out = []
    for selection in ("all", "pole_free"):
        boot = bootstrap_bands.loc[(bootstrap_bands.center == "best_fit_expression") &
                                   (bootstrap_bands.selection == selection) &
                                   (bootstrap_bands.quantity == "w")].reset_index(drop=True)
        approx = prior.loc[(prior.selection == selection) &
                           (prior.quantity == "w")].reset_index(drop=True)
        np.testing.assert_allclose(boot.z, approx.z, atol=1e-12, rtol=0)
        for index in (0, 64, 129, 193, 257, 300):
            row_b, row_a = boot.iloc[index], approx.iloc[index]
            out.append({"scenario": scenario, "branch": branch,
                        "selection": selection, "z": float(row_b.z),
                        "bootstrap_median": row_b["median"],
                        "linearized_median": row_a["median"],
                        "bootstrap_95_low": row_b.q025,
                        "bootstrap_95_high": row_b.q975,
                        "linearized_95_low": row_a.q025,
                        "linearized_95_high": row_a.q975})
    return pd.DataFrame(out)


def run_catalog_bootstrap(project_root, prior_experiment, n_draws=20000,
                          batch_size=500, base_seed=20261014,
                          resume_directory=None):
    """Run/resume each branch in a fresh experiment directory, saving batches."""
    root = Path(project_root).resolve()
    prior = Path(prior_experiment).resolve()
    if n_draws < 100 or batch_size < 1:
        raise ValueError("Need at least 100 draws and a positive batch size")
    prior_manifest = json.loads((prior / "experiment.json").read_text())
    if prior_manifest.get("status") != "complete":
        raise ValueError("Notebook 11 selection experiment is not complete")
    core = {"schema_version": 1, "prior_path": str(prior),
            "prior_sha256": _sha256(prior / "experiment.json"),
            "n_draws_per_center": int(n_draws), "batch_size": int(batch_size),
            "base_seed": int(base_seed), "centers": CENTERS,
            "scenario_branches": {"fig2_pantheon": ("independent", "refitted_exponential"),
                                  "fig3_desy5": ("independent",)},
            "rd_policy": "per-draw shared rd; exact 1/rd SN mean and covariance rescaling in distance-space GLS",
            "fit_policy": "full-covariance Gaussian mock catalogue; fixed family refitted per draw",
            "slope_bounds": SLOPE_BOUNDS, "exp_taylor_order": EXP_ORDER,
            "selection_policy": "same fixed structures as completed Notebook 11; no reselection per mock"}
    core_hash = _digest(core)
    root_out = root / "results/11b_fig23_catalog_bootstrap"
    root_out.mkdir(parents=True, exist_ok=True)
    if resume_directory is None:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        destination = root_out / f"fig23_catalog_bootstrap_{stamp}_{uuid4().hex[:8]}"
        destination.mkdir(exist_ok=False)
        manifest = {"core": core, "core_sha256": core_hash,
                    "created_utc": datetime.now(timezone.utc).isoformat(),
                    "status": "in_progress", "scenarios": {}}
        _write_json(destination / "experiment.json", manifest)
    else:
        destination = Path(resume_directory).resolve()
        manifest = json.loads((destination / "experiment.json").read_text())
        if manifest.get("core_sha256") != core_hash or _digest(manifest.get("core")) != core_hash:
            raise ValueError("Cannot resume with changed inputs or settings")
    grid = np.linspace(0, 2.33, 301)
    overview = []
    comparisons = []
    for scenario_index, scenario in enumerate(SCENARIOS):
        data_file = root / f"results/08_fig23_data_preparation/{scenario}.npz"
        data = load_standard_dataset(data_file)
        prior_scenario = json.loads((prior / scenario / "scenario.json").read_text())
        if _sha256(data_file) != prior_scenario["data_sha256"]:
            raise ValueError(f"Prepared data changed since Notebook 11: {scenario}")
        cosmos = data["metadata"]["config"]["cosmology"]
        rd0 = float(cosmos["rd_hinv_mpc"])
        rd_sigma = float(cosmos["rd_sigma_hinv_mpc"])
        omega_m = float(cosmos["Omega_m"])
        for branch_index, branch in enumerate(core["scenario_branches"][scenario]):
            family = ("quadratic" if scenario == "fig2_pantheon" and branch == "independent"
                      else "exponential")
            problem = WhitenedGLS(data, family)
            theta_hat = problem.fit_observed()
            previous = prior_scenario["branches"][branch]
            if abs(theta_hat[2] / len(data["y"]) - previous["loss"]) > 1e-7:
                raise ValueError("Deterministic fit does not match Notebook 11")
            centers = {"best_fit_expression": problem.predict(theta_hat[:2]),
                       "observed_data": problem.y}
            branch_dir = destination / scenario / branch
            branch_dir.mkdir(parents=True, exist_ok=True)
            branch_record = {"scenario": scenario, "branch": branch, "family": family,
                             "n_observations": len(data["y"]), "data_sha256": _sha256(data_file),
                             "coefficient_order": ["a", "b"],
                             "deterministic_fit": {"a": float(theta_hat[0]),
                                                   "b": float(theta_hat[1]),
                                                   "chi2": float(theta_hat[2]),
                                                   "loss": float(theta_hat[2] / len(data["y"]))},
                             "source_expression": previous["source_expression"],
                             "rd_mean": rd0, "rd_sigma": rd_sigma, "omega_m": omega_m,
                             "n_draws_per_center": n_draws, "centers": {}}
            all_bands = []
            for center_index, center in enumerate(CENTERS):
                center_dir = branch_dir / center
                batches_dir = center_dir / "batches"
                batches_dir.mkdir(parents=True, exist_ok=True)
                count_batches = math.ceil(n_draws / batch_size)
                for batch_index in range(count_batches):
                    batch_path = batches_dir / f"batch_{batch_index:04d}.npz"
                    if batch_path.exists():
                        continue
                    count = min(batch_size, n_draws - batch_index * batch_size)
                    sequence = np.random.SeedSequence(base_seed, spawn_key=(
                        scenario_index, branch_index, center_index, batch_index))
                    rng = np.random.default_rng(sequence)
                    rd = rng.normal(rd0, rd_sigma, size=count)
                    fits = problem.draw_and_fit(centers[center], rd, rd0, rng)
                    np.savez_compressed(batch_path, a=fits[:, 0], b=fits[:, 1],
                                        chi2=fits[:, 2], rd_hinv_mpc=rd)
                    if (batch_index + 1) % 10 == 0 or batch_index == count_batches - 1:
                        print(f"{scenario}/{branch}/{center}: {batch_index + 1}/{count_batches} batches", flush=True)
                chunks = [np.load(batches_dir / f"batch_{index:04d}.npz")
                          for index in range(count_batches)]
                columns = {key: np.concatenate([chunk[key] for chunk in chunks])
                           for key in ("a", "b", "chi2", "rd_hinv_mpc")}
                for chunk in chunks:
                    chunk.close()
                dm, w, finite, pole, monotonic = reconstruct_samples(
                    family, columns["a"], columns["b"], columns["rd_hinv_mpc"],
                    grid, omega_m)
                samples = pd.DataFrame({"draw_index": np.arange(n_draws), **columns,
                                        "valid": finite, "pole_crossing": pole,
                                        "distance_increasing": monotonic})
                samples.to_csv(center_dir / "fit_samples.csv.gz", index=False,
                               compression="gzip")
                tables = []
                for selection, mask in (("all", finite), ("pole_free", finite & ~pole)):
                    if mask.sum() < min(100, math.ceil(n_draws / 2)):
                        raise ValueError("Too few retained draws for pointwise bands")
                    tables.append(_bands(grid, dm[mask], "distance", center, selection))
                    tables.append(_bands(grid, w[mask], "w", center, selection))
                bands = pd.concat(tables, ignore_index=True)
                bands.to_csv(center_dir / "bands.csv", index=False)
                all_bands.append(bands)
                diagnostics = {"n_draws": n_draws, "n_valid": int(finite.sum()),
                               "n_pole_crossing": int((finite & pole).sum()),
                               "n_pole_free": int((finite & ~pole).sum()),
                               "n_nonincreasing": int((finite & ~monotonic).sum()),
                               "pointwise_intervals": True}
                _write_json(center_dir / "diagnostics.json", diagnostics)
                branch_record["centers"][center] = diagnostics
                overview.append({"scenario": scenario, "branch": branch,
                                 "family": family, "center": center,
                                 "deterministic_loss": branch_record["deterministic_fit"]["loss"],
                                 **diagnostics})
                _figure(center_dir / "reconstruction.png", data, grid, family,
                        center, bands, prior_scenario["printed_paper_expression"],
                        omega_m, rd0)
            full_bands = pd.concat(all_bands, ignore_index=True)
            full_bands.to_csv(branch_dir / "all_bands.csv", index=False)
            prior_bands = prior / scenario / f"{branch}_bands.csv"
            comparison = _compare_linearized(full_bands, prior_bands, scenario, branch)
            comparison.to_csv(branch_dir / "comparison_to_linearized.csv", index=False)
            comparisons.append(comparison)
            _write_json(branch_dir / "branch.json", branch_record)
            manifest["scenarios"].setdefault(scenario, {})[branch] = branch_record
            _write_json(destination / "experiment.json", manifest)
    pd.DataFrame(overview).to_csv(destination / "overview.csv", index=False)
    pd.concat(comparisons, ignore_index=True).to_csv(
        destination / "comparison_to_linearized.csv", index=False)
    manifest["status"] = "complete"
    _write_json(destination / "experiment.json", manifest)
    return destination, pd.DataFrame(overview)
