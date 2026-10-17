"""Notebook 11c: post-Notebook-12, subset-balanced audit of fixed PySR seeds.

This does not rerun PySR, refit coefficients, or replace Notebook 11's result.
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
from scipy.linalg import solve_triangular

from .data_loading import load_standard_dataset
from .fig23_reconstruction import compile_curve


def _hash(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def evaluate_fixed_candidates(candidates, data):
    """Use the two marginal covariance blocks; coefficients remain untouched."""
    if len(candidates) != 100 or candidates.seed.nunique() != 100:
        raise ValueError("Expected exactly 100 distinct Notebook 10 seed selections")
    if not (candidates.scenario.eq("fig2_pantheon") &
            candidates.desi_error_version.eq("public_release")).all():
        raise ValueError("The candidate file is not the expected Figure 2 public sample")
    if candidates.independent_primary.sum() != 1:
        raise ValueError("Expected one frozen Notebook 11 primary candidate")
    desi = np.flatnonzero(data["source"] == "DESI_DR2")
    sn = np.flatnonzero(data["source"] != "DESI_DR2")
    if len(desi) != 6 or not np.array_equal(desi, np.arange(6)):
        raise ValueError("Expected six DESI points at the start of the data vector")
    cov = data["covariance"]
    if np.any(cov[np.ix_(desi, sn)]):
        raise ValueError("Unexpected DESI/SN cross covariance")
    factors = {"DESI": np.linalg.cholesky(cov[np.ix_(desi, desi)]),
               "SN": np.linalg.cholesky(cov[np.ix_(sn, sn)])}
    indices = {"DESI": desi, "SN": sn}
    rows = []
    for row in candidates.itertuples(index=False):
        item = {"seed": int(row.seed), "complexity": int(row.complexity),
                "family": row.family, "expression": row.sympy_expression,
                "original_combined_loss": float(row.python_loss),
                "eligible": bool(row.eligible),
                "within_notebook11_loss_gate": bool(row.within_loss_gate),
                "notebook11_primary": bool(row.independent_primary)}
        if item["eligible"]:
            predicted = compile_curve(item["expression"])(data["z"])[0]
            if not np.isfinite(predicted).all():
                raise ValueError(f"Nonfinite data prediction for seed {row.seed}")
            chi2 = {}
            for name, ix in indices.items():
                residual = data["y"][ix] - predicted[ix]
                white = solve_triangular(factors[name], residual,
                                         lower=True, check_finite=False)
                chi2[name] = float(white @ white)
                item[f"{name.lower()}_chi2"] = chi2[name]
                item[f"{name.lower()}_chi2_per_point"] = chi2[name] / len(ix)
            item["combined_chi2_per_point"] = sum(chi2.values()) / len(data["z"])
            if not np.isclose(item["combined_chi2_per_point"],
                              item["original_combined_loss"], rtol=2e-6, atol=2e-6):
                raise ValueError(f"Covariance loss mismatch for seed {row.seed}")
            item["worst_subset_chi2_per_point"] = max(
                item["desi_chi2_per_point"], item["sn_chi2_per_point"])
        rows.append(item)
    return pd.DataFrame(rows)


def select_balanced(scores):
    """Exploratory rule disclosed before this run, after Notebook 12's audit."""
    finalists = scores.loc[scores.eligible &
                           scores.within_notebook11_loss_gate].copy()
    if finalists.empty:
        raise ValueError("No physically eligible candidates inside the loss gate")
    finalists = finalists.sort_values(
        ["worst_subset_chi2_per_point", "combined_chi2_per_point",
         "complexity", "seed"], kind="mergesort")
    return finalists.iloc[0]


def _plot_scores(path, scores, chosen):
    plot = scores.loc[scores.eligible].copy()
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharey=True)
    old = plot.loc[plot.notebook11_primary].iloc[0]
    for ax in axes:
        for family, group in plot.groupby("family"):
            ax.scatter(group.desi_chi2_per_point, group.sn_chi2_per_point,
                       label=f"{family} (n={len(group)})", alpha=.7, s=46)
        ax.scatter([old.desi_chi2_per_point], [old.sn_chi2_per_point],
                   marker="s", s=160, facecolors="none", edgecolors="black",
                   linewidths=1.7, label="Notebook 11 primary")
        ax.scatter([chosen.desi_chi2_per_point], [chosen.sn_chi2_per_point],
                   marker="*", s=240, facecolors="gold", edgecolors="black",
                   label="Notebook 11c balanced")
        ax.set(xlabel=r"DESI $\chi^2/6$")
        ax.set_xlim(left=0)
        ax.grid(alpha=.2)
    axes[0].set(xlim=(0, 12), ylabel=r"SN $\chi^2/1590$",
                title="Eligible candidates: useful score range")
    axes[1].set(title="All eligible candidates")
    axes[1].legend(fontsize=7, loc="upper center")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_data(path, data, scores, chosen):
    old = scores.loc[scores.notebook11_primary].iloc[0]
    grid = np.linspace(0, 2.33, 401)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].errorbar(data["z"][:6], data["y"][:6],
                     yerr=np.sqrt(np.diag(data["covariance"])[:6]),
                     fmt="ko", capsize=3, label="DESI DR2")
    for row, label, color in ((old, "Notebook 11 quadratic", "tab:blue"),
                              (chosen, "Notebook 11c selected", "tab:orange")):
        func = compile_curve(row.expression)
        axes[0].plot(grid, func(grid)[0], color=color, label=label)
        residual = ((data["y"][:6] - func(data["z"][:6])[0]) /
                    np.sqrt(np.diag(data["covariance"])[:6]))
        axes[1].plot(data["z"][:6], residual, "o-", color=color, label=label)
    axes[0].set(xlabel="Redshift z", ylabel=r"$D_M/r_d$", title="Six DESI measurements")
    axes[1].axhline(0, color="black", linewidth=.8)
    axes[1].set(xlabel="Redshift z", ylabel="Residual / marginal SD",
                title="DESI pointwise residuals")
    for ax in axes:
        ax.legend(fontsize=8)
        ax.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run_notebook11c(project_root, prior=None, output_root=None):
    root = Path(project_root).resolve()
    prior = Path(prior or root / "results/11_fig23_mc_and_reconstruction" /
                 "fig23_mc_20261010T044029855521Z_5cf83467").resolve()
    candidate_path = prior / "fig2_pantheon/candidate_diagnostics.csv"
    data_path = root / "results/08_fig23_data_preparation/fig2_pantheon.npz"
    upstream_path = root / "results/10_fig23_seed_stability" / \
        "fig23_seeds_100_20261004T071816359765Z_101dd054/selected_runs.csv"
    if not all(path.is_file() for path in (candidate_path, data_path, upstream_path)):
        raise FileNotFoundError("Notebook 08, 10 or 11 input is missing")
    candidates = pd.read_csv(candidate_path)
    upstream = pd.read_csv(upstream_path)
    source = upstream.loc[(upstream.scenario == "fig2_pantheon") &
                          (upstream.desi_error_version == "public_release")]
    if set(source.seed) != set(candidates.seed):
        raise ValueError("Notebook 11 candidate seeds differ from Notebook 10")
    data = load_standard_dataset(data_path)
    destination_root = Path(output_root or root / "results/11c_fig2_balanced_seed_audit")
    destination_root.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    destination = destination_root / f"fig2_balanced_{stamp}_{uuid4().hex[:8]}"
    destination.mkdir(exist_ok=False)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "in_progress", "scope": "Figure 2 only",
        "study_type": "exploratory post-Notebook-12 re-selection; not preregistered",
        "inputs": {str(path): _hash(path) for path in
                   (candidate_path, data_path, upstream_path)},
        "selection_rule": [
            "Use only the 100 fixed public-covariance Notebook 10 seed-selected expressions.",
            "Retain Notebook 11 physical eligibility and <=1.5x original combined-loss gate.",
            "No coefficient refit. Use separate marginal full-covariance DESI and SN chi2/N.",
            "Minimize max(DESI chi2/6, SN chi2/1590); break ties by combined loss, complexity, seed."],
        "limitations": [
            "Subset chi2/N is descriptive; six DESI points give it substantial sampling variance.",
            "One selected structure from 100 seeds does not demonstrate search stability.",
            "This audit does not compute new Monte Carlo confidence bands or reproduce Figure 2 completely."]}
    _json(destination / "experiment.json", manifest)
    try:
        scores = evaluate_fixed_candidates(candidates, data)
        chosen = select_balanced(scores)
        old = scores.loc[scores.notebook11_primary].iloc[0]
        scores["notebook11c_selected"] = scores.seed.eq(chosen.seed)
        scores.to_csv(destination / "candidate_subset_scores.csv", index=False)
        _plot_scores(destination / "candidate_subset_scores.png", scores, chosen)
        _plot_data(destination / "desi_fit_comparison.png", data, scores, chosen)
        summary = {"candidate_count": len(scores),
                   "eligible_count": int(scores.eligible.sum()),
                   "loss_gate_count": int((scores.eligible & scores.within_notebook11_loss_gate).sum()),
                   "notebook11_primary": {key: old[key] for key in
                      ("seed", "expression", "complexity", "desi_chi2_per_point",
                       "sn_chi2_per_point", "combined_chi2_per_point")},
                   "notebook11c_selected": {key: chosen[key] for key in
                      ("seed", "expression", "complexity", "desi_chi2_per_point",
                       "sn_chi2_per_point", "combined_chi2_per_point")}}
        for branch in ("notebook11_primary", "notebook11c_selected"):
            summary[branch]["seed"] = int(summary[branch]["seed"])
            summary[branch]["complexity"] = int(summary[branch]["complexity"])
        _json(destination / "summary.json", summary)
        lines = ["# Notebook 11c — Figure 2 balanced seed audit", "",
                 "Exploratory re-selection **after** Notebook 12 exposed the DESI/SN imbalance. "
                 "The original Notebook 11 result remains unchanged. No PySR rerun, coefficient "
                 "refit, or new Monte Carlo is performed.", "",
                 "The fixed rule is recorded in `experiment.json`. Separate subgroup scores "
                 "use the corresponding marginal covariance blocks; the joint score is their "
                 "sum divided by 1596.", "",
                 "| Branch | Seed | Complexity | DESI chi2/6 | SN chi2/1590 | Combined chi2/1596 |",
                 "| --- | ---: | ---: | ---: | ---: | ---: |"]
        for label, item in (("Notebook 11 primary", old), ("Notebook 11c selected", chosen)):
            lines.append(f"| {label} | {int(item.seed)} | {int(item.complexity)} | "
                         f"{item.desi_chi2_per_point:.5f} | {item.sn_chi2_per_point:.5f} | "
                         f"{item.combined_chi2_per_point:.5f} |")
        lines += ["", f"New fixed PySR expression: `{chosen.expression}`", "",
                  "The new choice is a seed-selected PySR expression, not the separately "
                  "postulated exponential comparison. Its appearance in only one of 100 "
                  "seed-selected rows does not show stable selection.", "",
                  "The second panel uses marginal DESI standard deviations for visual "
                  "residuals; the subgroup chi-squared uses its full covariance.", "",
                  "Next: if this exploratory choice is retained, compute fresh conditional "
                  "Monte Carlo bands under the same uncertainty protocol as Notebook 11/11b. "
                  "Do not relabel the existing comparison-model bands as independent bands.", ""]
        (destination / "analysis.md").write_text("\n".join(lines))
        manifest["status"] = "complete"
        manifest["selected_seed"] = int(chosen.seed)
        _json(destination / "experiment.json", manifest)
        return destination, summary
    except Exception:
        manifest["status"] = "failed"
        _json(destination / "experiment.json", manifest)
        raise
