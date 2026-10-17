"""Notebook 11: auditable Figure 2/3 selection and w(z) reconstruction.

The Monte Carlo uses a first-order GLS coefficient approximation.  It does
not re-run PySR or refit every synthetic 1600/1760-point catalogue.
"""

from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sympy as sp
from scipy.integrate import cumulative_trapezoid
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares

from .data_loading import load_standard_dataset


C_KM_S = 299792.458
SCENARIOS = ("fig2_pantheon", "fig3_desy5")
QUANTILES = (0.00135, 0.025, 0.16, 0.5, 0.84, 0.975, 0.99865)
QUANTILE_NAMES = ("q00135", "q025", "q16", "median", "q84", "q975", "q99865")


def _sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False,
                                     default=str) + "\n")


def parse_expression(expression):
    z = sp.Symbol("z", real=True)
    expr = sp.sympify(str(expression), locals={"z": z})
    if expr.free_symbols - {z}:
        raise ValueError(f"Unexpected symbols in {expression}")
    return z, expr


def _array(func, z, *params):
    with np.errstate(all="ignore"):
        value = np.asarray(func(z, *params), dtype=float)
    shape = np.broadcast_shapes(np.shape(z), value.shape)
    return np.broadcast_to(value, shape).copy()


def compile_curve(expression):
    z, expr = parse_expression(expression)
    functions = [sp.lambdify(z, part, modules="numpy")
                 for part in (expr, sp.diff(expr, z), sp.diff(expr, z, 2))]
    return lambda grid: tuple(_array(function, grid) for function in functions)


def gamma_value(omega_m, rd_hinv_mpc):
    return float(omega_m) * (100.0 * np.asarray(rd_hinv_mpc) / C_KM_S) ** 2


def w_from_derivatives(z, d1, d2, gamma):
    """Paper Eq. (5), with d1=d(D_M/r_d)/dz and d2=d2(D_M/r_d)/dz2."""
    matter = gamma * (1.0 + z) ** 3 * d1 ** 3
    factor = 1.0 - gamma * (1.0 + z) ** 3 * d1 ** 2
    with np.errstate(all="ignore"):
        w = -1.0 - (2.0 * (1.0 + z) * d2 + 3.0 * matter) / (3.0 * (d1 - matter))
    return w, factor


def candidate_diagnostics(expression, grid, omega_m=0.3, rd_hinv_mpc=101.0,
                          origin_tolerance=1e-8):
    try:
        dm, d1, d2 = compile_curve(expression)(grid)
        w, factor = w_from_derivatives(grid, d1, d2,
                                       gamma_value(omega_m, rd_hinv_mpc))
        finite = bool(np.isfinite(dm).all() and np.isfinite(d1).all()
                      and np.isfinite(d2).all() and np.isfinite(w).all())
        origin = bool(np.isfinite(dm[0]) and abs(dm[0]) <= origin_tolerance)
        positive = bool(np.isfinite(dm[1:]).all() and (dm[1:] > 0).all())
        increasing = bool(np.isfinite(d1).all() and (d1 > 0).all())
        no_pole = bool(np.isfinite(factor).all() and (factor > 0).all())
        return {"finite_derivatives_and_w": finite, "zero_at_origin": origin,
                "positive_distance": positive, "increasing_distance": increasing,
                "positive_pole_factor": no_pole,
                "minimum_pole_factor": float(np.nanmin(factor)),
                "eligible": finite and origin and positive and increasing and no_pole}
    except (ValueError, TypeError, OverflowError, ZeroDivisionError,
            NameError, AttributeError) as error:
        return {"finite_derivatives_and_w": False, "zero_at_origin": False,
                "positive_distance": False, "increasing_distance": False,
                "positive_pole_factor": False, "minimum_pole_factor": np.nan,
                "eligible": False, "diagnostic_error": str(error)}


def audit_and_select(selected, scenario, grid, loss_factor=1.5):
    """Predeclared selection: eligible; within 1.5x best; min complexity/loss/seed."""
    rows = selected.loc[(selected.scenario == scenario) &
                        (selected.desi_error_version == "public_release")].copy()
    if rows.empty:
        raise ValueError(f"No public-release Notebook 10 fits for {scenario}")
    diagnostics = pd.DataFrame([candidate_diagnostics(expr, grid)
                                for expr in rows.sympy_expression])
    rows = pd.concat([rows.reset_index(drop=True), diagnostics], axis=1)
    eligible = rows.loc[rows.eligible].copy()
    if eligible.empty:
        raise ValueError(f"No physically eligible candidate for {scenario}")
    threshold = float(loss_factor * eligible.python_loss.min())
    rows["within_loss_gate"] = rows.eligible & (rows.python_loss <= threshold)
    finalists = rows.loc[rows.within_loss_gate].sort_values(
        ["complexity", "python_loss", "seed"], kind="mergesort")
    chosen = finalists.iloc[0]
    rows["independent_primary"] = rows.seed.eq(chosen.seed)
    return rows, chosen, {"minimum_eligible_loss": float(eligible.python_loss.min()),
                          "loss_gate_factor": float(loss_factor),
                          "loss_gate_threshold": threshold,
                          "eligible_count": int(len(eligible)),
                          "finalist_count": int(len(finalists)),
                          "selected_seed": int(chosen.seed)}


def parameterize_expression(expression):
    z, expr = parse_expression(expression)
    floats = sorted(expr.atoms(sp.Float), key=lambda value: (float(value), str(value)))
    if not floats:
        raise ValueError("Chosen expression has no floating coefficients to propagate")
    parameters = sp.symbols(f"p0:{len(floats)}", real=True)
    template = expr.xreplace(dict(zip(floats, parameters)))
    initial = np.array([float(value) for value in floats], dtype=float)
    return z, template, parameters, initial


def _curve_functions(z, template, parameters):
    parts = [template, sp.diff(template, z), sp.diff(template, z, 2)]
    return [sp.lambdify((z, *parameters), part, modules="numpy") for part in parts]


def _fitted_model(data, expression):
    z, template, parameters, initial = parameterize_expression(expression)
    funcs = _curve_functions(z, template, parameters)
    observations = data["y"]
    redshifts = data["z"]
    covariance = data["covariance"]
    lower = np.linalg.cholesky(covariance)

    def predict(theta, grid=redshifts):
        return _array(funcs[0], grid, *theta)

    def whitened_residual(theta):
        values = predict(theta)
        if not np.isfinite(values).all():
            return np.full(len(values), 1e12)
        return solve_triangular(lower, values - observations, lower=True,
                                check_finite=False)

    fit = least_squares(whitened_residual, initial, max_nfev=150,
                        x_scale="jac", ftol=1e-10, xtol=1e-10, gtol=1e-10)
    if not fit.success or not np.isfinite(fit.x).all():
        raise RuntimeError(f"GLS coefficient refit failed: {fit.message}")
    fitted = fit.x
    derivative_functions = [sp.lambdify((z, *parameters), sp.diff(template, p),
                                        modules="numpy") for p in parameters]
    jacobian = np.column_stack([_array(function, redshifts, *fitted)
                                for function in derivative_functions])
    whitened_jacobian = solve_triangular(lower, jacobian, lower=True,
                                         check_finite=False)
    fisher = whitened_jacobian.T @ whitened_jacobian
    coefficient_covariance = np.linalg.inv(fisher)
    if not np.isfinite(coefficient_covariance).all():
        raise ValueError("Non-finite linearized coefficient covariance")
    np.linalg.cholesky(coefficient_covariance)
    # A single r_d realization rescales all SN distance measurements together.
    rd0 = float(data["metadata"]["config"]["cosmology"]["rd_hinv_mpc"])
    sn_response = np.where(data["source"] == "DESI_DR2", 0.0,
                           -observations / rd0)
    whitened_response = solve_triangular(lower, sn_response, lower=True,
                                         check_finite=False)
    rd_coefficient_response = coefficient_covariance @ (
        whitened_jacobian.T @ whitened_response)
    residual = whitened_residual(fitted)
    return {"template": str(template), "parameter_count": len(parameters),
            "parameters": fitted, "parameter_covariance": coefficient_covariance,
            "rd_coefficient_response": rd_coefficient_response,
            "chi2": float(residual @ residual),
            "loss": float(residual @ residual / len(redshifts)),
            "functions": funcs}


def _quantile_frame(z, values, label):
    frame = pd.DataFrame({"z": z, "quantity": label,
                          "n_finite": np.isfinite(values).sum(axis=0)})
    for quantile, name in zip(QUANTILES, QUANTILE_NAMES):
        frame[name] = np.nanquantile(values, quantile, axis=0)
    return frame


def _evaluate_samples(model, z, parameters, rd_values, omega_m):
    shape = (len(parameters), len(z))
    z2 = z[None, :]
    args = [parameters[:, index, None] for index in range(parameters.shape[1])]
    dm, d1, d2 = tuple(_array(func, z2, *args) for func in model["functions"])
    dm = np.broadcast_to(dm, shape).copy()
    d1 = np.broadcast_to(d1, shape).copy()
    d2 = np.broadcast_to(d2, shape).copy()
    gamma = gamma_value(omega_m, rd_values)[:, None]
    w, factor = w_from_derivatives(z2, d1, d2, gamma)
    valid_finite = np.isfinite(dm).all(axis=1) & np.isfinite(d1).all(axis=1) \
        & np.isfinite(d2).all(axis=1) & np.isfinite(w).all(axis=1)
    crossing = np.any(factor[:, :-1] * factor[:, 1:] <= 0, axis=1)
    crossing |= np.any(factor == 0, axis=1)
    return dm, w, factor, valid_finite, crossing


def monte_carlo_model(model, z, omega_m, rd_mean, rd_sigma, n_draws, seed):
    """First-order GLS coefficient draws plus one shared r_d draw per realization."""
    generator = np.random.default_rng(seed)
    rd_draws = generator.normal(rd_mean, rd_sigma, size=n_draws)
    noise = generator.multivariate_normal(np.zeros(model["parameter_count"]),
                                           model["parameter_covariance"], size=n_draws)
    coefficient_draws = (model["parameters"][None, :] + noise +
                         (rd_draws - rd_mean)[:, None] *
                         model["rd_coefficient_response"][None, :])
    dm, w, factor, finite, crossing = _evaluate_samples(
        model, z, coefficient_draws, rd_draws, omega_m)
    all_w = np.where(finite[:, None], w, np.nan)
    all_dm = np.where(finite[:, None], dm, np.nan)
    retained = finite & ~crossing
    bands = []
    for selection, mask in (("all", finite), ("pole_free", retained)):
        if not mask.any():
            raise ValueError(f"No finite Monte Carlo samples for {selection}")
        bands.append(_quantile_frame(z, dm[mask], "distance").assign(selection=selection))
        bands.append(_quantile_frame(z, w[mask], "w").assign(selection=selection))
    bands = pd.concat(bands, ignore_index=True)
    diagnostics = {"n_draws": int(n_draws), "n_finite": int(finite.sum()),
                   "n_pole_crossing": int((finite & crossing).sum()),
                   "n_pole_free": int(retained.sum()),
                   "pole_free_fraction_of_finite": float(retained.sum() / finite.sum()),
                   "minimum_absolute_pole_factor_q01": float(np.nanquantile(
                       np.min(np.abs(factor[finite]), axis=1), .01)),
                   "three_sigma_is_pointwise": True}
    return bands, diagnostics, {"coefficients": coefficient_draws,
                                 "rd_hinv_mpc": rd_draws,
                                 "finite": finite, "pole_crossing": crossing,
                                 "central_w_draws": all_w[:, 0],
                                 "central_dm_draws": all_dm[:, 0]}


def lcdm_distance(z, omega_m, rd_hinv_mpc):
    dense = np.linspace(0, max(z), max(4000, 10 * len(z)))
    h_over_h0 = np.sqrt(omega_m * (1 + dense) ** 3 + 1 - omega_m)
    integral = cumulative_trapezoid(1 / h_over_h0, dense, initial=0)
    return C_KM_S / (100 * rd_hinv_mpc) * np.interp(z, dense, integral)


def _seed_model_spread(candidate_rows, z, omega_m, rd_mean):
    curves = []
    for row in candidate_rows.loc[candidate_rows.eligible].itertuples():
        try:
            _, d1, d2 = compile_curve(row.sympy_expression)(z)
            w, factor = w_from_derivatives(z, d1, d2,
                                           gamma_value(omega_m, rd_mean))
            if np.isfinite(w).all() and (factor > 0).all():
                curves.append(w)
        except (ValueError, TypeError, OverflowError, NameError, AttributeError):
            continue
    if not curves:
        raise ValueError("No eligible seed-selected w curves")
    matrix = np.stack(curves)
    return _quantile_frame(z, matrix, "w_seed_spread"), len(curves)


def _plot_scenario(path, scenario, data, z, branch_bands, selected_curve,
                   printed_curve, seed_spread, lcdm_curve):
    fig, (left, right) = plt.subplots(1, 2, figsize=(14, 5.2))
    desi = data["source"] == "DESI_DR2"
    sn = ~desi
    left.scatter(data["z"][sn], data["y"][sn], s=2, alpha=.09,
                 color="gray", label="Supernovae")
    left.errorbar(data["z"][desi], data["y"][desi],
                  yerr=data["sigma"][desi], fmt="o", color="black",
                  markersize=4, label="DESI DR2")
    for branch, color, label in (("independent", "tab:blue", "Independent selection"),
                                 ("refitted_exponential", "tab:orange", "Refitted exponential")):
        bands = branch_bands[branch]
        distance = bands.loc[(bands.selection == "all") &
                             (bands.quantity == "distance")]
        left.plot(z, distance["median"], color=color, label=label)
        left.fill_between(z, distance.q025, distance.q975, color=color, alpha=.12)
    left.plot(z, printed_curve[0], ":", color="firebrick", label="Printed paper expression")
    left.plot(z, lcdm_curve, "--", color="darkgreen", label=r"$\Lambda$CDM")
    left.set(xlabel="Redshift z", ylabel=r"$D_M/r_d$", title=f"{scenario}: distance")
    left.legend(fontsize=7, loc="upper left")

    right.fill_between(z, seed_spread.q16, seed_spread.q84, color="gray", alpha=.16,
                       label="16–84% across eligible seeds (not CI)")
    for branch, color, label in (("independent", "tab:blue", "Independent selection"),
                                 ("refitted_exponential", "tab:orange", "Refitted exponential")):
        wbands = branch_bands[branch]
        wband = wbands.loc[(wbands.selection == "all") & (wbands.quantity == "w")]
        right.plot(z, wband["median"], color=color, label=label)
        right.fill_between(z, wband.q025, wband.q975, color=color, alpha=.14)
    right.plot(z, selected_curve[1], color="tab:blue", lw=.8, alpha=.6)
    right.plot(z, printed_curve[1], ":", color="firebrick",
               label="Printed paper expression")
    right.axhline(-1, color="darkgreen", ls="--", label=r"$w=-1$")
    right.set(xlabel="Redshift z", ylabel="w(z)", title="Pointwise 95% bands; model dependent",
              ylim=(-4, 2))
    right.legend(fontsize=7, loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def run_notebook11(project_root, upstream, n_mc=10000, mc_seed=20261011,
                   output_root=None):
    """Run both scenarios into one new, non-overwriting experiment directory."""
    root = Path(project_root).resolve()
    upstream = Path(upstream).resolve()
    record = json.loads((upstream / "experiment.json").read_text())
    if record.get("status") != "complete":
        raise ValueError("Notebook 10 upstream experiment is not marked complete")
    source_table = upstream / "selected_runs.csv"
    selected = pd.read_csv(source_table)
    grid = np.linspace(0, 2.33, 301)
    output_root = Path(output_root or root / "results/11_fig23_mc_and_reconstruction")
    output_root.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    destination = output_root / f"fig23_mc_{stamp}_{uuid4().hex[:8]}"
    destination.mkdir(exist_ok=False)
    manifest = {"created_utc": datetime.now(timezone.utc).isoformat(),
                "status": "in_progress",
                "upstream": str(upstream), "upstream_experiment_sha256": _sha256(upstream / "experiment.json"),
                "selected_runs_sha256": _sha256(source_table), "n_mc": int(n_mc),
                "mc_seed": int(mc_seed), "grid": {"start": 0.0, "end": 2.33, "points": len(grid)},
                "selection_rule": "public_release; finite D,D',D'',w; D(0)=0; D>0 for z>0; D'>0; positive pole factor; raw loss <=1.5x minimum eligible; lowest complexity then loss then seed",
                "method": "first-order GLS coefficient-covariance Monte Carlo; one shared rd draw per realization; SN calibration response linearized; no PySR refit per draw",
                "branches": {"independent": "Notebook 10 selected structure, GLS-refitted coefficients",
                             "refitted_exponential": "separate a*z*exp(b*z) comparison",
                             "printed_paper": "fixed printed expression, no fitted uncertainty band"},
                "limitations": ["Conditional coefficient bands exclude symbolic model-selection uncertainty.",
                                "Seed-spread bands are not observational confidence intervals.",
                                "3-sigma bands are pointwise, not a global hypothesis test.",
                                "The coefficient covariance and rd response use a first-order approximation."],
                "scenarios": {}}
    _json(destination / "experiment.json", manifest)
    overview = []
    for index, scenario in enumerate(SCENARIOS):
        data_file = root / f"results/08_fig23_data_preparation/{scenario}.npz"
        data = load_standard_dataset(data_file)
        expected = record["core"]["prepared_data"][scenario]["sha256"]
        if _sha256(data_file) != expected:
            raise ValueError(f"Prepared input changed since Notebook 10: {scenario}")
        scenario_dir = destination / scenario
        scenario_dir.mkdir()
        audit, primary, decision = audit_and_select(selected, scenario, grid)
        audit.to_csv(scenario_dir / "candidate_diagnostics.csv", index=False)
        cosmos = data["metadata"]["config"]["cosmology"]
        omega_m = float(cosmos["Omega_m"])
        rd_mean = float(cosmos["rd_hinv_mpc"])
        rd_sigma = float(cosmos["rd_sigma_hinv_mpc"])
        paper = data["metadata"]["config"]["references"]["paper_expression"]
        if paper != record["core"]["paper_expressions"][scenario]:
            raise ValueError("Paper reference changed between notebooks")
        branches = {"independent": str(primary.sympy_expression),
                    "refitted_exponential": "30.0*z*exp(-0.25*z)"}
        branch_bands, branch_metadata = {}, {}
        for branch, expression in branches.items():
            model = _fitted_model(data, expression)
            bands, diagnostics, draws = monte_carlo_model(
                model, grid, omega_m, rd_mean, rd_sigma, n_mc, mc_seed + index * 100 +
                (0 if branch == "independent" else 1))
            branch_bands[branch] = bands
            bands.to_csv(scenario_dir / f"{branch}_bands.csv", index=False)
            np.savez_compressed(scenario_dir / f"{branch}_mc_samples.npz", **draws)
            branch_metadata[branch] = {"source_expression": expression,
                                       "refitted_template": model["template"],
                                       "refitted_parameters": model["parameters"].tolist(),
                                       "coefficient_covariance": model["parameter_covariance"].tolist(),
                                       "rd_coefficient_response": model["rd_coefficient_response"].tolist(),
                                       "chi2": model["chi2"], "loss": model["loss"],
                                       "diagnostics": diagnostics}
        selected_d, selected_d1, selected_d2 = compile_curve(primary.sympy_expression)(grid)
        selected_w, selected_factor = w_from_derivatives(
            grid, selected_d1, selected_d2, gamma_value(omega_m, rd_mean))
        paper_d, paper_d1, paper_d2 = compile_curve(paper)(grid)
        paper_w, paper_factor = w_from_derivatives(
            grid, paper_d1, paper_d2, gamma_value(omega_m, rd_mean))
        seed_spread, seed_count = _seed_model_spread(audit, grid, omega_m, rd_mean)
        seed_spread.to_csv(scenario_dir / "eligible_seed_w_spread.csv", index=False)
        lcdm = lcdm_distance(grid, omega_m, rd_mean)
        pd.DataFrame({"z": grid, "raw_selected_distance": selected_d,
                      "raw_selected_w": selected_w, "raw_selected_pole_factor": selected_factor,
                      "printed_paper_distance": paper_d, "printed_paper_w": paper_w,
                      "printed_paper_pole_factor": paper_factor,
                      "lcdm_distance": lcdm, "lcdm_w": -np.ones_like(grid)}).to_csv(
                          scenario_dir / "reference_curves.csv", index=False)
        _plot_scenario(scenario_dir / "figure_reconstruction.png", scenario,
                       data, grid, branch_bands, (selected_d, selected_w),
                       (paper_d, paper_w), seed_spread, lcdm)
        metadata = {"scenario": scenario, "data_file": str(data_file),
                    "data_sha256": _sha256(data_file), "n_observations": len(data["z"]),
                    "selection": decision, "raw_selected_expression": str(primary.sympy_expression),
                    "raw_selected_loss": float(primary.python_loss),
                    "raw_selected_complexity": int(primary.complexity),
                    "eligible_seed_spread_count": seed_count,
                    "printed_paper_expression": paper, "branches": branch_metadata}
        _json(scenario_dir / "scenario.json", metadata)
        manifest["scenarios"][scenario] = metadata
        overview.append({"scenario": scenario, "selected_seed": int(primary.seed),
                         "selected_expression": str(primary.sympy_expression),
                         "selected_raw_loss": float(primary.python_loss),
                         "selected_refit_loss": branch_metadata["independent"]["loss"],
                         "exponential_refit_loss": branch_metadata["refitted_exponential"]["loss"],
                         "eligible_seed_count": seed_count,
                         "primary_pole_free_fraction": branch_metadata["independent"]["diagnostics"]["pole_free_fraction_of_finite"]})
        _json(destination / "experiment.json", manifest)
    pd.DataFrame(overview).to_csv(destination / "overview.csv", index=False)
    manifest["status"] = "complete"
    _json(destination / "experiment.json", manifest)
    return destination, pd.DataFrame(overview)
