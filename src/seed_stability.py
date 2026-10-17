"""Resumable, auditable seed ensembles for the Figure 2/3 PySR fits."""

from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4
import hashlib
import json
import re

import numpy as np
import pandas as pd
import sympy as sp
import matplotlib.pyplot as plt

from .symbolic_regression import (
    chi_square, evaluate_sympy, export_fit_result, fit_exponential_reference,
    fit_symbolic_distance,
)


def _utc_now():
    return datetime.now(timezone.utc).isoformat()


def _json_digest(value):
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _file_digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")
    temporary.replace(path)


def classify_expression(expression):
    """Broad structural label; exact a*z*exp(b*z) gets comparable a and b."""
    z = sp.Symbol("z")
    try:
        expr = sp.sympify(str(expression), locals={"z": z})
        if not expr.has(z):
            return {"family": "constant", "exp_a": None, "exp_b": None}
        quotient = sp.simplify(expr / z)
        slope = sp.simplify(sp.diff(quotient, z) / quotient)
        if not slope.has(z):
            intercept = sp.limit(quotient, z, 0)
            a, b = float(intercept), float(slope)
            if np.isfinite([a, b]).all() and a != 0:
                probes = np.array([0.1, 0.4, 1.0, 2.0])
                observed = evaluate_sympy(expr, probes)
                target = a * probes * np.exp(b * probes)
                if np.isfinite(observed).all() and np.allclose(
                        observed, target, rtol=1e-8, atol=1e-9):
                    return {"family": "a*z*exp(b*z)", "exp_a": a, "exp_b": b}
        try:
            degree = sp.Poly(expr, z).degree()
            return {"family": f"polynomial_degree_{degree}",
                    "exp_a": None, "exp_b": None}
        except sp.PolynomialError:
            pass
        if expr.has(sp.log):
            family = "logarithmic_or_mixed"
        elif expr.has(sp.exp) or any(
                power.exp.has(z) for power in expr.atoms(sp.Pow)):
            family = "exponential_or_variable_power_mixed"
        elif any(power.exp.is_Rational and power.exp.q != 1
                 for power in expr.atoms(sp.Pow)):
            family = "fractional_power_or_radical"
        elif expr.is_rational_function(z):
            family = "rational"
        else:
            family = "other_mixed"
        return {"family": family, "exp_a": None, "exp_b": None}
    except (TypeError, ValueError, sp.SympifyError, OverflowError):
        return {"family": "unclassified", "exp_a": None, "exp_b": None}


def _validated_seeds(seeds):
    numbers = [int(seed) for seed in seeds]
    if not numbers or len(set(numbers)) != len(numbers) or any(seed < 0 for seed in numbers):
        raise ValueError("Provide a nonempty list of unique, nonnegative seeds")
    return numbers


def setup_experiment(output_root, config, prepared_paths, scenarios,
                     error_versions, seeds, paper_expressions,
                     resume_directory=None):
    """Create an experiment or resume it; an existing seed list may be extended."""
    scenarios, error_versions = list(scenarios), list(error_versions)
    seeds = _validated_seeds(seeds)
    if len(set(scenarios)) != len(scenarios) or len(set(error_versions)) != len(error_versions):
        raise ValueError("Scenarios and error versions must be unique")
    if set(scenarios) != set(prepared_paths) or set(scenarios) != set(paper_expressions):
        raise ValueError("Each scenario needs one prepared dataset and paper expression")
    core = {
        "schema_version": 1,
        "config": config,
        "scenarios": scenarios,
        "error_versions": error_versions,
        "prepared_data": {name: {"path": str(Path(prepared_paths[name]).resolve()),
                                  "sha256": _file_digest(prepared_paths[name])}
                          for name in scenarios},
        "paper_expressions": dict(paper_expressions),
    }
    core_hash = _json_digest(core)
    if resume_directory is None:
        root = Path(output_root)
        root.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        directory = (root / f"fig23_seeds_{len(seeds)}_{stamp}_{uuid4().hex[:8]}").resolve()
        directory.mkdir(parents=False, exist_ok=False)
        record = {"core": core, "core_sha256": core_hash, "seeds": seeds,
                  "created_utc": _utc_now(), "updated_utc": _utc_now(),
                  "status": "in_progress", "seed_extensions": []}
    else:
        directory = Path(resume_directory).resolve()
        record_path = directory / "experiment.json"
        if not record_path.is_file():
            raise FileNotFoundError(f"No experiment manifest: {record_path}")
        record = json.loads(record_path.read_text())
        if record.get("core_sha256") != core_hash or record.get("core") != core:
            raise ValueError("Cannot resume: configuration, inputs, or paper references changed")
        prior = record["seeds"]
        if seeds[:len(prior)] != prior:
            raise ValueError("Resume seeds must preserve the existing seed-list prefix")
        if len(seeds) > len(prior):
            record["seed_extensions"].append(
                {"utc": _utc_now(), "previous_count": len(prior), "new_count": len(seeds)})
            record["seeds"] = seeds
            record["status"] = "in_progress"
            record["updated_utc"] = _utc_now()
    _write_json(directory / "experiment.json", record)
    return directory


def _seed_path(directory, scenario, error_version, seed):
    return Path(directory) / "runs" / scenario / error_version / f"seed_{seed}"


def _completed_attempt(seed_path):
    completed = sorted(path for path in seed_path.glob("attempt_*")
                       if (path / "run_complete.json").is_file())
    if len(completed) > 1:
        raise ValueError(f"More than one completed attempt for {seed_path}")
    return completed[0] if completed else None


def collect_runs(directory):
    """Return a planned-run manifest and one row per completed seed."""
    directory = Path(directory)
    record = json.loads((directory / "experiment.json").read_text())
    manifest, selected = [], []
    for scenario in record["core"]["scenarios"]:
        for version in record["core"]["error_versions"]:
            for seed in record["seeds"]:
                seed_path = _seed_path(directory, scenario, version, seed)
                done = _completed_attempt(seed_path)
                attempts = sorted(seed_path.glob("attempt_*"))
                status = "complete" if done else "partial" if attempts else "pending"
                manifest.append({"scenario": scenario, "desi_error_version": version,
                                 "seed": seed, "status": status,
                                 "attempts": len(attempts),
                                 "completed_path": str(done.relative_to(directory)) if done else ""})
                if done:
                    marker = json.loads((done / "run_complete.json").read_text())
                    if marker["scenario"] != scenario or marker["desi_error_version"] != version or marker["seed"] != seed:
                        raise ValueError(f"Corrupt completion marker: {done}")
                    classification = marker.get("classification")
                    if classification is None:
                        classification = classify_expression(
                            marker["summary"]["sympy_expression"])
                    selected.append({**marker["summary"], **classification,
                                     "run_path": str(done.relative_to(directory))})
    return pd.DataFrame(manifest), pd.DataFrame(selected)


def save_progress(directory):
    manifest, selected = collect_runs(directory)
    manifest.to_csv(Path(directory) / "seed_manifest.csv", index=False)
    selected.to_csv(Path(directory) / "selected_runs.csv", index=False)
    return manifest, selected


def run_seed_ensemble(directory, variants, config, paper_expressions):
    """Save every fit immediately; retry partial fits in a new attempt folder."""
    directory = Path(directory)
    record = json.loads((directory / "experiment.json").read_text())
    scenarios = record["core"]["scenarios"]
    versions = record["core"]["error_versions"]
    if set(variants) != {(scenario, version) for scenario in scenarios for version in versions}:
        raise ValueError("Variants do not match the experiment manifest")
    if record["core"]["config"] != config or record["core"]["paper_expressions"] != paper_expressions:
        raise ValueError("Current fit configuration or paper expressions differ from the manifest")
    save_progress(directory)
    for scenario in scenarios:
        for version in versions:
            data = variants[scenario, version]
            benchmark = None
            paper_grid = None
            for seed in record["seeds"]:
                seed_path = _seed_path(directory, scenario, version, seed)
                if _completed_attempt(seed_path) is not None:
                    continue
                if benchmark is None:
                    benchmark = fit_exponential_reference(data)
                    paper_data = evaluate_sympy(paper_expressions[scenario], data["z"])
                    benchmark["paper_expression_loss"] = (
                        chi_square(data["y"], paper_data, data["covariance"]) / len(data["y"]))
                seed_path.mkdir(parents=True, exist_ok=True)
                attempt_numbers = [int(match.group(1)) for path in seed_path.glob("attempt_*")
                                   if (match := re.fullmatch(r"attempt_(\d+)", path.name))]
                attempt = seed_path / f"attempt_{max(attempt_numbers, default=0) + 1:03d}"
                print(f"RUN {scenario} / {version} / seed {seed} -> {attempt}", flush=True)
                try:
                    result = fit_symbolic_distance(data, config, seed, attempt)
                    if paper_grid is None:
                        paper_grid = evaluate_sympy(paper_expressions[scenario], result["z_grid"])
                    result["summary"].update({
                        "paper_expression": paper_expressions[scenario],
                        "paper_expression_loss": benchmark["paper_expression_loss"],
                        "benchmark_a": benchmark["a"], "benchmark_b": benchmark["b"],
                        "benchmark_loss": benchmark["loss"],
                    })
                    export_fit_result(result, attempt, config, paper_grid,
                                      benchmark, data["metadata"])
                    _write_json(attempt / "run_complete.json", {
                        "scenario": scenario, "desi_error_version": version,
                        "seed": seed, "completed_utc": _utc_now(),
                        "summary": result["summary"],
                        "classification": classify_expression(
                            result["summary"]["sympy_expression"]),
                    })
                    save_progress(directory)
                except Exception as error:
                    attempt.mkdir(parents=True, exist_ok=True)
                    _write_json(attempt / "failure.json", {
                        "utc": _utc_now(), "error_type": type(error).__name__,
                        "message": str(error),
                    })
                    save_progress(directory)
                    raise
    manifest, selected = save_progress(directory)
    record["status"] = "complete" if manifest.status.eq("complete").all() else "in_progress"
    record["updated_utc"] = _utc_now()
    _write_json(directory / "experiment.json", record)
    return manifest, selected


def summarize_runs(directory):
    """Write family frequencies, loss/physical diagnostics, and curve quantiles."""
    directory = Path(directory)
    manifest, selected = save_progress(directory)
    if selected.empty:
        raise ValueError("No completed seed runs to summarize")
    group = ["scenario", "desi_error_version"]
    family = selected.groupby(group + ["family"], dropna=False).size().rename("count").reset_index()
    family["fraction_of_completed"] = family["count"] / family.groupby(group)["count"].transform("sum")
    family.to_csv(directory / "structure_counts.csv", index=False)
    planned = manifest.groupby(group).size().rename("n_planned")
    counts = selected.groupby(group).agg(
        n_completed=("seed", "size"),
        median_loss=("python_loss", "median"),
        loss_q10=("python_loss", lambda x: x.quantile(.10)),
        loss_q90=("python_loss", lambda x: x.quantile(.90)),
        exponential_fraction=("family", lambda x: x.eq("a*z*exp(b*z)").mean()),
        finite_fraction=("finite_on_grid", "mean"),
        positive_fraction=("positive_for_z_gt_zero", "mean"),
        zero_origin_fraction=("satisfies_dm_zero_constraint", "mean"),
        median_elapsed_seconds=("elapsed_seconds", "median"),
    ).join(planned).reset_index()
    counts["completion_fraction"] = counts.n_completed / counts.n_planned
    counts.to_csv(directory / "stability_summary.csv", index=False)
    curve_rows = []
    for (scenario, version), rows in selected.groupby(group):
        grids = []
        common_z = None
        for _, row in rows.iterrows():
            grid = pd.read_csv(directory / row.run_path / "distance_grid.csv")
            z = grid.z.to_numpy(dtype=float)
            if common_z is None:
                common_z = z
            else:
                np.testing.assert_array_equal(common_z, z)
            values = grid.selected.to_numpy(dtype=float)
            grids.append(np.where(np.isfinite(values), values, np.nan))
        stack = np.stack(grids)
        for j, z in enumerate(common_z):
            finite = stack[:, j][np.isfinite(stack[:, j])]
            quantiles = np.quantile(finite, [.025, .16, .5, .84, .975]) if len(finite) else [np.nan] * 5
            curve_rows.append({"scenario": scenario, "desi_error_version": version,
                               "z": z, "n_finite": len(finite),
                               "q025": quantiles[0], "q16": quantiles[1],
                               "median": quantiles[2], "q84": quantiles[3],
                               "q975": quantiles[4]})
    curves = pd.DataFrame(curve_rows)
    curves.to_csv(directory / "distance_curve_quantiles.csv", index=False)
    return counts, family, curves, selected


def save_summary_plots(directory, paper_expressions):
    """Plot loss distributions and pointwise between-seed distance envelopes."""
    directory = Path(directory)
    _, _, curves, selected = summarize_runs(directory)
    groups = list(selected.groupby(["scenario", "desi_error_version"], sort=False))
    figure, ax = plt.subplots(figsize=(max(8, len(groups) * 2), 5))
    for position, ((scenario, version), rows) in enumerate(groups, start=1):
        ax.boxplot(rows.python_loss.to_numpy(), positions=[position], widths=.5)
        ax.scatter(position, rows.paper_expression_loss.iloc[0], marker="x",
                   color="firebrick", zorder=3,
                   label="Printed paper formula" if position == 1 else None)
    ax.set_xticks(range(1, len(groups) + 1),
                  [f"{scenario}\n{version}" for (scenario, version), _ in groups])
    ax.set_ylabel(r"Full-covariance $\chi^2/N$")
    ax.set_title("Selected PySR loss across seeds")
    ax.legend(frameon=False)
    figure.tight_layout()
    loss_path = directory / "loss_across_seeds.png"
    figure.savefig(loss_path, dpi=180)
    plt.close(figure)

    scenarios = list(dict.fromkeys(key[0] for key, _ in groups))
    versions = list(dict.fromkeys(key[1] for key, _ in groups))
    figure, axes = plt.subplots(len(scenarios), len(versions),
                                figsize=(6 * len(versions), 4 * len(scenarios)),
                                squeeze=False, sharex=True, sharey=True)
    for i, scenario in enumerate(scenarios):
        for j, version in enumerate(versions):
            ax = axes[i, j]
            subset = curves.loc[(curves.scenario == scenario) &
                                (curves.desi_error_version == version)]
            if subset.empty:
                ax.set_visible(False)
                continue
            z = subset.z.to_numpy()
            ax.fill_between(z, subset.q16.to_numpy(), subset.q84.to_numpy(),
                            color="steelblue", alpha=.25, label="16–84% of seeds")
            ax.plot(z, subset["median"].to_numpy(), color="steelblue",
                    label="Median selected fit")
            ax.plot(z, evaluate_sympy(paper_expressions[scenario], z), "--",
                    color="firebrick", label="Printed paper formula")
            ax.set(title=f"{scenario} · {version}", xlabel="Redshift z",
                   ylabel=r"$D_M/r_d$")
            ax.legend(frameon=False, fontsize=8)
    figure.tight_layout()
    distance_path = directory / "distance_seed_envelopes.png"
    figure.savefig(distance_path, dpi=180)
    plt.close(figure)
    return loss_path, distance_path
