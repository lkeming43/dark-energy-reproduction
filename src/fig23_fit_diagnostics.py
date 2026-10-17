"""Notebook 12: frozen-model, separated-data fit diagnostics for Figures 2/3.

No structure search or coefficient refit occurs here.  This audit uses the
completed Notebook 11 coefficients and the immutable Notebook 08 data files.
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
from scipy.linalg import solve_triangular

from .data_loading import load_standard_dataset
from .fig23_reconstruction import compile_curve


SCENARIOS = ("fig2_pantheon", "fig3_desy5")
BRANCHES = ("independent", "refitted_exponential", "printed_paper")
LABELS = {"independent": "Independent selection",
          "refitted_exponential": "Refitted exponential",
          "printed_paper": "Printed paper expression"}
COLORS = {"independent": "tab:blue", "refitted_exponential": "tab:orange",
          "printed_paper": "firebrick"}


def _sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False,
                                     default=str) + "\n")


def _markdown_table(frame):
    columns = list(frame.columns)
    lines = ["| " + " | ".join(columns) + " |",
             "| " + " | ".join("---" for _ in columns) + " |"]
    for values in frame.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(str(value) for value in values) + " |")
    return "\n".join(lines)


def fixed_prediction(z, branch, scenario_record):
    """Evaluate the coefficients saved by Notebook 11, without refitting."""
    z = np.asarray(z, dtype=float)
    if branch == "printed_paper":
        return compile_curve(scenario_record["printed_paper_expression"])(z)[0]
    record = scenario_record["branches"][branch]
    parameters = np.asarray(record["refitted_parameters"], dtype=float)
    symbols = sp.symbols(f"p0:{len(parameters)}", real=True)
    z_symbol = sp.Symbol("z", real=True)
    template = sp.sympify(record["refitted_template"],
                          locals={"z": z_symbol, **{str(p): p for p in symbols}})
    if template.free_symbols - {z_symbol, *symbols}:
        raise ValueError("Unexpected symbol in frozen Notebook 11 template")
    function = sp.lambdify((z_symbol, *symbols), template, modules="numpy")
    with np.errstate(all="ignore"):
        prediction = np.asarray(function(z, *parameters), dtype=float)
    prediction = np.broadcast_to(prediction, z.shape).copy()
    if not np.isfinite(prediction).all():
        raise ValueError(f"Non-finite {branch} prediction")
    return prediction


def covariance_diagnostics(observed, predicted, covariance):
    """Exact marginal-subset chi-square and two *different* residual views.

    Marginal standardized residuals are interpretable point-by-point but are
    correlated for supernovae.  Whitened residuals sum in quadrature to chi2,
    but individual whitened coordinates depend on the supplied row ordering.
    """
    observed, predicted = np.asarray(observed), np.asarray(predicted)
    covariance = np.asarray(covariance)
    if observed.ndim != 1 or observed.shape != predicted.shape:
        raise ValueError("Observed and predicted vectors must have matching 1D shape")
    if covariance.shape != (len(observed), len(observed)):
        raise ValueError("Covariance dimension mismatch")
    lower = np.linalg.cholesky(covariance)
    residual = observed - predicted
    whitened = solve_triangular(lower, residual, lower=True, check_finite=False)
    marginal = residual / np.sqrt(np.diag(covariance))
    return {"residual": residual, "marginal_sigma": marginal,
            "whitened": whitened, "chi2": float(whitened @ whitened),
            "n": len(observed)}


def _model_table(data, scenario_record):
    desi = np.flatnonzero(data["source"] == "DESI_DR2")
    sn = np.flatnonzero(data["source"] != "DESI_DR2")
    if len(desi) != 6 or desi.tolist() != list(range(6)):
        raise ValueError("Notebook 08 contract changed: expected DESI in first six rows")
    if np.max(np.abs(data["covariance"][np.ix_(desi, sn)])) != 0:
        raise ValueError("DESI/SN cross-covariance is no longer zero")
    table, residuals = [], {}
    for branch in BRANCHES:
        predicted = fixed_prediction(data["z"], branch, scenario_record)
        parts = {}
        for name, indices in (("DESI", desi), ("SN", sn)):
            part = covariance_diagnostics(data["y"][indices], predicted[indices],
                                          data["covariance"][np.ix_(indices, indices)])
            parts[name] = part
            table.append({"scenario": scenario_record["scenario"], "branch": branch,
                          "subset": name, "n": part["n"], "chi2": part["chi2"],
                          "chi2_per_observation": part["chi2"] / part["n"],
                          "marginal_residual_rms": float(np.sqrt(np.mean(part["marginal_sigma"] ** 2))),
                          "largest_absolute_marginal_sigma": float(np.max(np.abs(part["marginal_sigma"])))})
        total = parts["DESI"]["chi2"] + parts["SN"]["chi2"]
        table.append({"scenario": scenario_record["scenario"], "branch": branch,
                      "subset": "combined", "n": len(data["z"]), "chi2": total,
                      "chi2_per_observation": total / len(data["z"]),
                      "marginal_residual_rms": np.nan,
                      "largest_absolute_marginal_sigma": np.nan})
        if branch != "printed_paper":
            expected = scenario_record["branches"][branch]["chi2"]
            if not np.isclose(total, expected, rtol=1e-9, atol=1e-6):
                raise ValueError(f"Frozen {branch} chi2 no longer matches Notebook 11")
        residuals[branch] = {"prediction": predicted, "DESI": parts["DESI"],
                             "SN": parts["SN"]}
    return pd.DataFrame(table), residuals


def _point_table(data, scenario_record, residuals):
    rows = []
    for branch in BRANCHES:
        for index in range(len(data["z"])):
            subset = "DESI" if index < 6 else "SN"
            local_index = index if index < 6 else index - 6
            part = residuals[branch][subset]
            rows.append({"scenario": scenario_record["scenario"], "branch": branch,
                         "subset": subset, "id": data["id"][index], "row_index": int(data["row_index"][index]),
                         "z": float(data["z"][index]), "observed": float(data["y"][index]),
                         "predicted": float(residuals[branch]["prediction"][index]),
                         "marginal_sigma": float(part["marginal_sigma"][local_index]),
                         "whitened_residual": float(part["whitened"][local_index])})
    return pd.DataFrame(rows)


def _highest_desi_check(data, scenario_record, residuals):
    index = int(np.argmax(data["z"][:6]))
    keep = np.array([i for i in range(6) if i != index])
    rows = []
    for branch in BRANCHES:
        prediction = residuals[branch]["prediction"]
        reduced = covariance_diagnostics(data["y"][keep], prediction[keep],
                                         data["covariance"][np.ix_(keep, keep)])
        desi = residuals[branch]["DESI"]
        rows.append({"scenario": scenario_record["scenario"], "branch": branch,
                     "excluded_id": str(data["id"][index]), "excluded_z": float(data["z"][index]),
                     "excluded_observed": float(data["y"][index]),
                     "excluded_prediction": float(prediction[index]),
                     "excluded_marginal_sigma": float(desi["marginal_sigma"][index]),
                     "desi_chi2_all_six": desi["chi2"],
                     "desi_chi2_other_five_fixed_coefficients": reduced["chi2"],
                     "no_refit": True})
    return pd.DataFrame(rows)


def _binned_median(z, values, bins=24):
    edges = np.linspace(float(np.min(z)), float(np.max(z)), bins + 1)
    centers, medians = [], []
    for left, right in zip(edges[:-1], edges[1:]):
        inside = (z >= left) & (z < right if right < edges[-1] else z <= right)
        if inside.sum() >= 4:
            centers.append(float(np.median(z[inside])))
            medians.append(float(np.median(values[inside])))
    return np.asarray(centers), np.asarray(medians)


def _plot(path, data, scenario, residuals, summary):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    desi_z, sn_z = data["z"][:6], data["z"][6:]
    dense = np.linspace(0, 2.33, 400)
    axes[0, 0].errorbar(desi_z, data["y"][:6], yerr=data["sigma"][:6],
                        fmt="ko", ms=4, label="DESI DR2")
    for branch in BRANCHES:
        label, color = LABELS[branch], COLORS[branch]
        axes[0, 0].plot(dense, residuals[branch]["dense_prediction"],
                        color=color, lw=1.5, ls=":" if branch == "printed_paper" else "-", label=label)
        axes[0, 1].plot(desi_z, residuals[branch]["DESI"]["marginal_sigma"],
                        marker="o", ms=4, color=color, lw=1.2, label=label)
    axes[0, 0].set(xlabel="Redshift z", ylabel=r"$D_M/r_d$", title="Six DESI measurements")
    axes[0, 1].set(xlabel="Redshift z", ylabel="Residual / marginal SD",
                   title="DESI residuals (six-point covariance diagonal here)")
    axes[0, 1].axhline(0, color="black", lw=.7)
    axes[0, 1].axhline(3, color="gray", lw=.7, ls="--")
    axes[0, 1].axhline(-3, color="gray", lw=.7, ls="--")
    base = residuals["independent"]["SN"]["marginal_sigma"]
    axes[1, 0].scatter(sn_z, base, s=3, alpha=.10, color="gray",
                       label="Individual SN: independent model")
    for branch in BRANCHES:
        zbin, median = _binned_median(sn_z, residuals[branch]["SN"]["marginal_sigma"])
        axes[1, 0].plot(zbin, median, marker=".", ms=3, lw=1.2,
                        color=COLORS[branch], label=f"Binned median: {LABELS[branch]}")
        axes[1, 1].hist(residuals[branch]["SN"]["whitened"], bins=45, range=(-4, 4),
                        density=True, histtype="step", lw=1.2,
                        color=COLORS[branch], label=LABELS[branch])
    axes[1, 0].axhline(0, color="black", lw=.7)
    axes[1, 0].set(xlabel="Redshift z", ylabel="Residual / marginal SD",
                   title="SN marginal residuals (correlated; visual diagnostic)")
    x = np.linspace(-4, 4, 300)
    axes[1, 1].plot(x, np.exp(-x*x/2)/np.sqrt(2*np.pi), "k--", lw=1, label="N(0,1) reference")
    axes[1, 1].set(xlabel="Whitened residual coordinate", ylabel="Density",
                   title="SN covariance-whitened residual distribution")
    axes[0, 0].legend(fontsize=7)
    axes[1, 0].legend(fontsize=6)
    axes[1, 1].legend(fontsize=7)
    fig.suptitle(f"{scenario}: frozen Notebook 11 fits; nominal calibration", fontsize=12)
    fig.savefig(path, dpi=170)
    plt.close(fig)


def _plot_group_scores(path, summary):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for ax, scenario in zip(axes, SCENARIOS):
        table = summary[(summary.scenario == scenario) & (summary.subset != "combined")]
        xs = np.arange(len(BRANCHES))
        for offset, subset in ((-.18, "DESI"), (.18, "SN")):
            values = [float(table[(table.branch == b) & (table.subset == subset)]
                            .chi2_per_observation.iloc[0]) for b in BRANCHES]
            ax.bar(xs + offset, values, width=.34, label=subset)
        ax.set_xticks(xs, ["Independent", "Refit exp", "Paper print"], rotation=15)
        ax.set(title=scenario, ylabel=r"Subset $\chi^2/N$ (descriptive)")
        ax.legend(fontsize=8)
    fig.savefig(path, dpi=170)
    plt.close(fig)


def _write_analysis(path, summary, high_z, prior, bootstrap_summary):
    def score(scenario, branch, subset):
        chosen = summary[(summary.scenario == scenario) &
                         (summary.branch == branch) & (summary.subset == subset)]
        return float(chosen.chi2_per_observation.iloc[0])

    lines = ["# Notebook 12 — Figure 2/3 frozen-fit diagnostic", "",
             f"Completed {datetime.now(timezone.utc).date().isoformat()} (UTC).",
             "The Notebook 11 model structures and coefficients are frozen. No PySR, coefficient refit, or new Monte Carlo run was performed.",
             "The nominal Notebook 08 DESI and SN data/covariance are used; the DESI and SN blocks have zero cross-covariance.",
             "", "## Subset scores", "",
             "The table reports exact covariance-weighted chi-square for each marginal subset. Chi-square per observation is descriptive, not reduced chi-square or a p-value. Both fitted branches were optimized jointly, so no unique per-subset fitted degrees of freedom are assigned.",
             "", _markdown_table(summary[["scenario", "branch", "subset", "n", "chi2", "chi2_per_observation"]].round(5)),
             "", "## Highest-redshift DESI point", "",
             "This is a single-point marginal-standard-deviation residual, not a global significance test. The five-point score below merely drops that point while holding each published/fitted curve fixed; it does not refit.",
             "", _markdown_table(high_z.round(5)), "",
             "## Notebook 11b stability context", "",
             "These are existing full-catalogue Monte Carlo results, not newly generated by Notebook 12. Pole-free counts and nonincreasing-curve counts are diagnostics, not posterior model probabilities.",
             "", _markdown_table(bootstrap_summary), "",
             "## Evidence-based assessment", "",
             "**Technical reproducibility:** the fixed-branch combined chi-squared values exactly reproduce the completed Notebook 11 values within numerical tolerance. Dataset and scenario hashes are recorded in the manifest. This verifies arithmetic/provenance, not the scientific interpretation.",
             f"**Figure 2 independent quadratic:** combined loss {score('fig2_pantheon', 'independent', 'combined'):.3f} conceals DESI chi-squared per point {score('fig2_pantheon', 'independent', 'DESI'):.3f} versus SN {score('fig2_pantheon', 'independent', 'SN'):.3f}. The highest-redshift DESI point has a +4.61 marginal-SD residual. This branch cannot be described as a uniformly good joint fit without qualification.",
             f"**Figure 2 refitted exponential comparator:** DESI chi-squared per point is {score('fig2_pantheon', 'refitted_exponential', 'DESI'):.3f} and combined loss is {score('fig2_pantheon', 'refitted_exponential', 'combined'):.3f}. This is a better balanced fit on these diagnostics, but it is a separately labelled refit, not a replacement for the independently selected PySR result.",
             f"**Figure 3 independent exponential:** combined loss {score('fig3_desy5', 'independent', 'combined'):.3f}, DESI chi-squared per point {score('fig3_desy5', 'independent', 'DESI'):.3f}, SN {score('fig3_desy5', 'independent', 'SN'):.3f}. It fits the two subsets more evenly than the Figure 2 quadratic, while some individual DESI residuals remain near 2.7 marginal SD.",
             f"**Paper comparison:** the printed Figure 2 expression has combined loss {score('fig2_pantheon', 'printed_paper', 'combined'):.3f} on our prepared public data, so the published-loss discrepancy remains unresolved. The printed Figure 3 expression gives {score('fig3_desy5', 'printed_paper', 'combined'):.3f}.",
             "**Cosmological claim:** the existing Notebook 11b pointwise pole-free w(z) interval at z=0 includes -1 for the Figure 2 quadratic but excludes it for the exponential branches. Consequently the endpoint statement is model-family/conditioning dependent and does not establish a global exclusion of LambdaCDM. This validation notebook is complete, but scientific reproduction of all published losses and inference is not closed.", "",
             "## Interpretation and decision guardrails", "",
             "The combined loss weights all observations through the stated covariance but can obscure a poor fit in the six-point DESI subset beside thousands of supernovae. Compare both subgroup scores before describing an expression as a good joint fit.",
             "The SN per-point residual plot uses marginal standard deviations; those points are correlated. The SN histogram uses Cholesky-whitened coordinates, whose individual values depend on the fixed catalogue ordering. Only the sum of their squares is used as the full-covariance chi-square.",
             "All scores are in prepared distance space at nominal r_d. They do not resolve the paper's published-loss mismatch or constitute a test that globally rejects LambdaCDM. Notebook 11b's pole-free/all-draw w(z) sensitivity remains a separate, conditional inference.",
             "Do not change the Notebook 11 primary structure after seeing this audit. If a different selection criterion is proposed, declare it first and run it as a new experiment, with new conditional MC bands.",
             "", "## Provenance", "",
             f"Frozen Notebook 11 run: `{prior}`", "",
             "The experiment manifest records hashes of every prepared dataset and frozen scenario manifest. Figures and CSVs are in this new run folder.", ""]
    path.write_text("\n".join(lines))


def run_fit_diagnostics(project_root, prior, output_root=None):
    """Create one non-overwriting Notebook 12 run containing both scenarios."""
    root, prior = Path(project_root).resolve(), Path(prior).resolve()
    prior_manifest = json.loads((prior / "experiment.json").read_text())
    if prior_manifest.get("status") != "complete":
        raise ValueError("Notebook 11 prior is not complete")
    output_root = Path(output_root or root / "results/12_fig23_fit_diagnostics")
    output_root.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    destination = output_root / f"fig23_fit_check_{stamp}_{uuid4().hex[:8]}"
    destination.mkdir(exist_ok=False)
    manifest = {"status": "in_progress", "created_utc": datetime.now(timezone.utc).isoformat(),
                "method": "frozen Notebook 11 coefficients; separated full-covariance subset chi2; no new fit",
                "prior": str(prior), "prior_experiment_sha256": _sha256(prior / "experiment.json"),
                "scenarios": {}}
    _write_json(destination / "experiment.json", manifest)
    summaries, high_z_frames = [], []
    for scenario in SCENARIOS:
        data_path = root / "results/08_fig23_data_preparation" / f"{scenario}.npz"
        record_path = prior / scenario / "scenario.json"
        record = json.loads(record_path.read_text())
        if record["scenario"] != scenario or _sha256(data_path) != record["data_sha256"]:
            raise ValueError(f"Frozen data provenance mismatch for {scenario}")
        data = load_standard_dataset(data_path)
        summary, residuals = _model_table(data, record)
        high_z = _highest_desi_check(data, record, residuals)
        points = _point_table(data, record, residuals)
        dense = np.linspace(0, 2.33, 400)
        for branch in BRANCHES:
            residuals[branch]["dense_prediction"] = fixed_prediction(dense, branch, record)
        scenario_dir = destination / scenario
        scenario_dir.mkdir()
        summary.to_csv(scenario_dir / "subset_scores.csv", index=False)
        points.to_csv(scenario_dir / "point_residuals.csv.gz", index=False,
                      compression="gzip", float_format="%.10g")
        high_z.to_csv(scenario_dir / "highest_z_desi_check.csv", index=False)
        _plot(scenario_dir / "residual_diagnostics.png", data, scenario, residuals, summary)
        manifest["scenarios"][scenario] = {"data_sha256": _sha256(data_path),
                                            "scenario_json_sha256": _sha256(record_path),
                                            "n_total": len(data["z"]), "n_desi": 6,
                                            "n_sn": len(data["z"]) - 6,
                                            "branches": list(BRANCHES)}
        _write_json(destination / "experiment.json", manifest)
        summaries.append(summary)
        high_z_frames.append(high_z)
    all_summary = pd.concat(summaries, ignore_index=True)
    all_high_z = pd.concat(high_z_frames, ignore_index=True)
    all_summary.to_csv(destination / "subset_scores.csv", index=False)
    all_high_z.to_csv(destination / "highest_z_desi_check.csv", index=False)
    _plot_group_scores(destination / "subset_score_comparison.png", all_summary)
    bootstrap_dir = root / "results/11b_fig23_catalog_bootstrap/fig23_catalog_bootstrap_20261014T033440741851Z_c88d56d5"
    bootstrap = pd.read_csv(bootstrap_dir / "overview.csv")
    bootstrap_summary = bootstrap.loc[bootstrap["center"] == "best_fit_expression",
                                      ["scenario", "branch", "n_draws", "n_pole_free", "n_nonincreasing"]].copy()
    bootstrap_summary.to_csv(destination / "notebook11b_stability_context.csv", index=False)
    manifest["notebook11b_overview"] = str(bootstrap_dir / "overview.csv")
    manifest["notebook11b_overview_sha256"] = _sha256(bootstrap_dir / "overview.csv")
    _write_analysis(destination / "analysis.md", all_summary, all_high_z, prior, bootstrap_summary)
    manifest["status"] = "complete"
    _write_json(destination / "experiment.json", manifest)
    return destination, all_summary
