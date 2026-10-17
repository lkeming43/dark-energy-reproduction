"""Cross-figure audit of the losses printed in Tables III--V.

This experiment uses only public/local inputs.  It compares several objective
definitions for every legible printed candidate and performs a deterministic
bounded search in the rounding interval implied by every printed decimal.
It does not rerun PySR and does not tune the data or covariance.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ModuleNotFoundError:  # The numerical audit does not require plotting.
    plt = None


ROOT = Path(__file__).resolve().parent
RESULTS_ROOT = ROOT / "results/09b_cross_figure_loss_objective_audit"
PAPER = ROOT / "arXiv-2502.10506v2.pdf"


def parameter(value, decimals, label):
    half_width = 0.5 * 10.0 ** (-decimals)
    return {"value": float(value), "lower": float(value - half_width),
            "upper": float(value + half_width), "label": label,
            "decimals": int(decimals)}


def candidate(figure, complexity, paper_loss, equation, fn, params=(), status="literal"):
    return {"figure": figure, "complexity": int(complexity),
            "formula_id": f"{figure}_c{complexity}",
            "paper_loss": float(paper_loss), "equation": equation,
            "fn": fn, "params": list(params), "parse_status": status}


def candidate_bank():
    p = parameter
    return [
        # Table III: source-verified. C13/C17 contain a double minus. C15 is
        # malformed in the official LaTeX source and uses a documented repair.
        candidate("fig1", 1, 2519.06, "2.8069465",
                  lambda z, a: np.full_like(z, a), [p(2.8069465, 7, "a")]),
        candidate("fig1", 2, 120.73, "sqrt(276.74304)",
                  lambda z, a: np.full_like(z, np.sqrt(a)), [p(276.74304, 5, "a")]),
        candidate("fig1", 3, 91.5, "17.705526-z",
                  lambda z, a: a-z, [p(17.705526, 6, "a")]),
        candidate("fig1", 4, 1.22, "(5.127707-z)^2",
                  lambda z, a: (a-z)**2, [p(5.127707, 6, "a")]),
        candidate("fig1", 5, .30, "28.42/1.68^z",
                  lambda z, a, b: a/b**z, [p(28.42, 2, "a"), p(1.68, 2, "b")]),
        candidate("fig1", 7, .25, "(3.77-z)^2.32+6.24",
                  lambda z, a, b, c: (a-z)**b+c,
                  [p(3.77, 2, "a"), p(2.32, 2, "b"), p(6.24, 2, "c")]),
        candidate("fig1", 8, .24, "((z-3.18)^2+3.47)*2.06",
                  lambda z, a, b, c: ((z-a)**2+b)*c,
                  [p(3.18, 2, "a"), p(3.47, 2, "b"), p(2.06, 2, "c")]),
        candidate("fig1", 9, .21, "2.70^z-11.8*z+26.0",
                  lambda z, a, b, c: a**z-b*z+c,
                  [p(2.70, 2, "a"), p(11.8, 1, "b"), p(26.0, 1, "c")]),
        candidate("fig1", 11, .17, "4.19*(1-0.787*z)^6+(4.94-z)^2.1",
                  lambda z, a, b, c, d: a*(1-b*z)**6+(c-z)**d,
                  [p(4.19, 2, "a"), p(.787, 3, "b"), p(4.94, 2, "c"), p(2.1, 1, "d")]),
        candidate("fig1", 13, .020, "(3.68-z)^2.4+6.56- -0.02/(z-0.76)",
                  lambda z, a, b, c, d, e: (a-z)**b+c+d/(z-e),
                  [p(3.68, 2, "a"), p(2.4, 1, "b"), p(6.56, 2, "c"),
                   p(.02, 2, "d"), p(.76, 2, "e")]),
        candidate("fig1", 15, .016, "minimal repair: (3.68-z)^2.37+6.56- -0.0169/[0.862*(z-0.763)]",
                  lambda z, a, b, c, d, e, f: (a-z)**b+c+d/(e*(z-f)),
                  [p(3.68, 2, "a"), p(2.37, 2, "b"), p(6.56, 2, "c"),
                   p(.0169, 4, "d"), p(.862, 3, "e"), p(.763, 3, "f")],
                  status="minimal_repair_of_malformed_source"),
        candidate("fig1", 16, .015, "21.9*(1-0.272*z)^2.37+6.56+0.0169/[sqrt(z)*(z-0.763)]",
                  lambda z, a, b, c, d, e, f: a*(1-b*z)**c+d+e/(np.sqrt(z)*(z-f)),
                  [p(21.9, 1, "a"), p(.272, 3, "b"), p(2.37, 2, "c"),
                   p(6.56, 2, "d"), p(.0169, 4, "e"), p(.763, 3, "f")]),
        candidate("fig1", 17, .013, "(3.68-z)^2.37+6.56- -0.0169/[z*(z-0.763-0.0216)]",
                  lambda z, a, b, c, d, e, f: (a-z)**b+c+d/(z*(z-e-f)),
                  [p(3.68, 2, "a"), p(2.37, 2, "b"), p(6.56, 2, "c"),
                   p(.0169, 4, "d"), p(.763, 3, "e"), p(.0216, 4, "f")]),

        # Table IV: DESI + Pantheon+.
        candidate("fig2", 1, 124.5, "1.24", lambda z, a: np.full_like(z, a), [p(1.24, 2, "a")]),
        candidate("fig2", 2, 116.74, "exp(z)", lambda z: np.exp(z)),
        candidate("fig2", 3, 4.65, "27.1*z", lambda z, a: a*z, [p(27.1, 1, "a")]),
        candidate("fig2", 5, 2.69, "(z/0.034)^0.915", lambda z, a, b: (z/a)**b,
                  [p(.034, 3, "a"), p(.915, 3, "b")]),
        candidate("fig2", 6, 2.39, "z*30.20*(0.18*z-1)^2",
                  lambda z, a, b: z*a*(b*z-1)**2, [p(30.20, 2, "a"), p(.18, 2, "b")]),
        candidate("fig2", 7, 1.70, "Table IV literal: 0.76^2*z/0.03",
                  lambda z, a, b: a**2*z/b, [p(.76, 2, "a"), p(.03, 2, "b")],
                  status="source_verified_table_literal_conflicts_with_appendix"),
        candidate("fig2", 9, 1.68, "0.77^z*(z*29.96)^0.99",
                  lambda z, a, b, c: a**z*(z*b)**c,
                  [p(.77, 2, "a"), p(29.96, 2, "b"), p(.99, 2, "c")]),
        candidate("fig2", 10, 1.68, "z^3+z*(29.8-z*8.24)",
                  lambda z, a, b: z**3+z*(a-z*b), [p(29.8, 1, "a"), p(8.24, 2, "b")]),
        candidate("fig2", 13, 1.68, "(z^3-z^2*7.60+z*27.1)/0.907",
                  lambda z, a, b, c: (z**3-z**2*a+z*b)/c,
                  [p(7.60, 2, "a"), p(27.1, 1, "b"), p(.907, 3, "c")]),
        candidate("fig2", 14, 1.68, "0.739^z*z*(z*(sqrt(z)-0.528)+30.0)",
                  lambda z, a, b, c: a**z*z*(z*(np.sqrt(z)-b)+c),
                  [p(.739, 3, "a"), p(.528, 3, "b"), p(30.0, 1, "c")]),
        candidate("fig2", 17, 1.68, "z*30.0/(0.196*0.465+exp(0.266))^z+z^2.13",
                  lambda z, a, b, c, d, e: z*a/(b*c+np.exp(d))**z+z**e,
                  [p(30.0, 1, "a"), p(.196, 3, "b"), p(.465, 3, "c"),
                   p(.266, 3, "d"), p(2.13, 2, "e")]),
        candidate("fig2", 19, 1.68, "30.0*z/sqrt(1.95^z)+z^2.13",
                  lambda z, a, b, c: a*z/np.sqrt(b**z)+z**c,
                  [p(30.0, 1, "a"), p(1.95, 2, "b"), p(2.13, 2, "c")],
                  status="literal_scope_from_rendered_pdf"),

        # Table V: DESI + DESY5.
        candidate("fig3", 1, 796., "z", lambda z: z),
        candidate("fig3", 2, 499., "8.65", lambda z, a: np.full_like(z, a), [p(8.65, 2, "a")]),
        candidate("fig3", 3, 22.5, "21.5*z", lambda z, a: a*z, [p(21.5, 1, "a")]),
        candidate("fig3", 5, 4.82, "(z*54.8)^0.77", lambda z, a, b: (z*a)**b,
                  [p(54.8, 1, "a"), p(.77, 2, "b")]),
        candidate("fig3", 6, 4.82, "22.1*z^0.774", lambda z, a, b: a*z**b,
                  [p(22.1, 1, "a"), p(.774, 3, "b")]),
        candidate("fig3", 7, 2.23, "30.1*z/1.31^z", lambda z, a, b: a*z/b**z,
                  [p(30.1, 1, "a"), p(1.31, 2, "b")]),
        candidate("fig3", 8, 2.10, "29.44^z*z", lambda z, a: a**z*z, [p(29.44, 2, "a")]),
        candidate("fig3", 9, 2.05, "z*(z+30.1/1.37^z)",
                  lambda z, a, b: z*(z+a/b**z), [p(30.1, 1, "a"), p(1.37, 2, "b")]),
        candidate("fig3", 11, 2.05, "z*(z+30.0/1.37^z)+0.0361",
                  lambda z, a, b, c: z*(z+a/b**z)+c,
                  [p(30.0, 1, "a"), p(1.37, 2, "b"), p(.0361, 4, "c")]),
        candidate("fig3", 13, 2.05, "z*(z^2+31.7)/sqrt(1.88^z)-z",
                  lambda z, a, b: z*(z**2+a)/np.sqrt(b**z)-z,
                  [p(31.7, 1, "a"), p(1.88, 2, "b")]),
        candidate("fig3", 14, 2.04, "0.0129+z*((z^2)^0.923+30.5)/1.38^z",
                  lambda z, a, b, c, d: a+z*((z**2)**b+c)/d**z,
                  [p(.0129, 4, "a"), p(.923, 3, "b"), p(30.5, 1, "c"), p(1.38, 2, "d")]),
        candidate("fig3", 18, 2.04, "31.1*z*fourth_root(z^(0.2*z))/sqrt(1.85^z)",
                  lambda z, a, b, c: a*z*(z**(b*z))**.25/np.sqrt(c**z),
                  [p(31.1, 1, "a"), p(.2, 1, "b"), p(1.85, 2, "c")]),
        candidate("fig3", 20, 2.04, "z*(z+30^(-z))-0.0218/(z-1.61)",
                  lambda z, a, b, c: z*(z+a**(-z))-b/(z-c),
                  [p(30., 0, "a"), p(.0218, 4, "b"), p(1.61, 2, "c")]),
    ]


def selected_expression_variants():
    """Printed representations of the selected equations in each paper location."""
    p = parameter
    specs = [
        ("fig1_table_iii", "fig1", 5, .30, "Table III", "28.42/1.68^z",
         lambda z, a, b: a/b**z, [p(28.42, 2, "a"), p(1.68, 2, "b")]),
        ("fig1_eq10", "fig1", 5, .30, "Eq. (10)", "28.4*exp(-0.519*z)",
         lambda z, a, b: a*np.exp(-b*z), [p(28.4, 1, "a"), p(.519, 3, "b")]),
        ("fig2_table_iv_literal", "fig2", 7, 1.70, "Table IV official source",
         "0.76^2*z/0.03", lambda z, a, b: a**2*z/b,
         [p(.76, 2, "a"), p(.03, 2, "b")]),
        ("fig2_appendix_a2", "fig2", 7, 1.70, "Appendix Eq. (A2)",
         "0.76^z*z/0.03", lambda z, a, b: a**z*z/b,
         [p(.76, 2, "a"), p(.03, 2, "b")]),
        ("fig2_eq12", "fig2", 7, 1.70, "Eq. (12)", "33.3*z*exp(-0.274*z)",
         lambda z, a, b: a*z*np.exp(-b*z), [p(33.3, 1, "a"), p(.274, 3, "b")]),
        ("fig3_table_v", "fig3", 7, 2.23, "Table V", "30.1*z/1.31^z",
         lambda z, a, b: a*z/b**z, [p(30.1, 1, "a"), p(1.31, 2, "b")]),
        ("fig3_appendix_a3", "fig3", 7, 2.23, "Appendix Eq. (A3)", "30.1*z/1.30^z",
         lambda z, a, b: a*z/b**z, [p(30.1, 1, "a"), p(1.30, 2, "b")]),
        ("fig3_eq13", "fig3", 7, 2.23, "Eq. (13)", "30.1*z*exp(-0.262*z)",
         lambda z, a, b: a*z*np.exp(-b*z), [p(30.1, 1, "a"), p(.262, 3, "b")]),
    ]
    result = []
    for formula_id, figure, complexity, loss, location, equation, fn, params in specs:
        item = candidate(figure, complexity, loss, equation, fn, params)
        item["formula_id"] = formula_id
        item["location"] = location
        result.append(item)
    return result


def load_inputs():
    files = {"fig1": ROOT / "results/01_desi_data_check/desi_dr2_dh_over_rd.npz",
             "fig2": ROOT / "results/08_fig23_data_preparation/fig2_pantheon.npz",
             "fig3": ROOT / "results/08_fig23_data_preparation/fig3_desy5.npz"}
    result = {}
    for figure, path in files.items():
        d = np.load(path, allow_pickle=False)
        covariance = d["covariance"].astype(float)
        result[figure] = {"path": path, "z": d["z"].astype(float),
                          "y": d["y"].astype(float), "covariance": covariance,
                          "precision": np.linalg.inv(covariance),
                          "variance": np.diag(covariance).copy()}
    return result


def predictions(item, z, theta_rows):
    values = []
    for theta in theta_rows:
        with np.errstate(all="ignore"):
            value = np.asarray(item["fn"](z, *theta), dtype=float)
        if value.ndim == 0:
            value = np.full_like(z, float(value))
        values.append(value)
    return np.vstack(values)


def full_losses(item, data, theta_rows):
    pred = predictions(item, data["z"], theta_rows)
    residual = pred - data["y"][None, :]
    finite = np.isfinite(residual).all(axis=1)
    result = np.full(len(theta_rows), np.inf)
    if finite.any():
        r = residual[finite]
        result[finite] = np.einsum("bi,ij,bj->b", r, data["precision"], r,
                                   optimize=True) / len(data["y"])
    return result


def objective_row(item, data):
    theta = [p["value"] for p in item["params"]]
    pred = predictions(item, data["z"], [theta])[0]
    residual = pred - data["y"]
    n = len(residual)
    full_chi2 = float(residual @ data["precision"] @ residual)
    diagonal_chi2 = float(np.sum(residual**2 / data["variance"]))
    sse = float(residual @ residual)
    metrics = {"full_cov_chi2_over_n": full_chi2/n,
               "diagonal_cov_chi2_over_n": diagonal_chi2/n,
               "unweighted_sse": sse, "unweighted_mse": sse/n,
               "full_cov_chi2": full_chi2}
    paper = item["paper_loss"]
    closest = min(metrics, key=lambda key: abs(np.log10(max(metrics[key], 1e-300)/paper)))
    return {"figure": item["figure"], "formula_id": item["formula_id"],
            "complexity": item["complexity"], "equation": item["equation"],
            "parse_status": item["parse_status"], "n": n,
            "paper_loss": paper, **metrics,
            "full_cov_to_paper_ratio": metrics["full_cov_chi2_over_n"]/paper,
            "closest_objective": closest,
            "closest_ratio": metrics[closest]/paper}


def rounding_scan(item, data):
    params = item["params"]
    center = np.array([p["value"] for p in params], dtype=float)
    if item["parse_status"].startswith("minimal_repair"):
        return {"center_loss": np.nan, "sampled_min_loss": np.nan,
                "sampled_max_loss": np.nan, "nearest_loss": np.nan,
                "nearest_abs_difference": np.nan,
                "paper_loss_bracketed": False, "evaluations": 0,
                "parameter_box": json.dumps(params),
                "scan_status": "skipped_malformed_source_formula"}
    if not params:
        value = float(full_losses(item, data, [center])[0])
        return {"center_loss": value, "sampled_min_loss": value,
                "sampled_max_loss": value, "nearest_loss": value,
                "nearest_abs_difference": abs(value-item["paper_loss"]),
                "paper_loss_bracketed": bool(np.isclose(value, item["paper_loss"])),
                "evaluations": 1, "parameter_box": "[]", "scan_status": "complete"}
    lower = np.array([p["lower"] for p in params]); upper = np.array([p["upper"] for p in params])
    samples = [center]
    for bits in itertools.product([0, 1], repeat=len(params)):
        samples.append(np.where(np.asarray(bits, dtype=bool), upper, lower))
    # Deterministic interior coverage; every sample remains strictly in the
    # rounding box.  This is a bounded numerical search, not a global proof
    # when the reported loss is not reached.
    rng = np.random.default_rng(20260922 + 100*{"fig1":1,"fig2":2,"fig3":3}[item["figure"]]
                                + item["complexity"])
    for row in rng.uniform(lower, upper, size=(32, len(params))):
        samples.append(row)
    samples = np.unique(np.asarray(samples), axis=0)
    losses = full_losses(item, data, samples)
    finite = np.isfinite(losses)
    samples, losses = samples[finite], losses[finite]
    if not len(losses):
        return {"center_loss": np.inf, "sampled_min_loss": np.inf,
                "sampled_max_loss": np.inf, "nearest_loss": np.inf,
                "nearest_abs_difference": np.inf, "paper_loss_bracketed": False,
                "evaluations": 0, "parameter_box": json.dumps(params),
                "scan_status": "no_finite_evaluation"}
    # Two coordinate sweeps refine both the minimum and maximum found above.
    all_losses = list(losses); eval_count = len(losses)
    for mode in ("min", "max"):
        theta = samples[int(np.argmin(losses) if mode == "min" else np.argmax(losses))].copy()
        for _ in range(2):
            for j in range(len(params)):
                trial = np.tile(theta, (9, 1)); trial[:, j] = np.linspace(lower[j], upper[j], 9)
                values = full_losses(item, data, trial); eval_count += len(values)
                valid = np.isfinite(values)
                if valid.any():
                    ids = np.flatnonzero(valid)
                    chosen = ids[np.argmin(values[valid]) if mode == "min" else np.argmax(values[valid])]
                    theta, value = trial[chosen].copy(), float(values[chosen])
                    all_losses.append(value)
    all_losses = np.asarray(all_losses)
    paper = item["paper_loss"]
    nearest = float(all_losses[np.argmin(np.abs(all_losses-paper))])
    lo, hi = float(np.min(all_losses)), float(np.max(all_losses))
    center_loss = float(full_losses(item, data, [center])[0])
    return {"center_loss": center_loss, "sampled_min_loss": lo,
            "sampled_max_loss": hi, "nearest_loss": nearest,
            "nearest_abs_difference": abs(nearest-paper),
            "paper_loss_bracketed": bool(lo <= paper <= hi),
            "evaluations": int(eval_count), "parameter_box": json.dumps(params),
            "scan_status": "complete"}


def markdown_table(frame, columns):
    view = frame[columns].copy()
    for column in view.select_dtypes(include=[np.number]).columns:
        view[column] = view[column].map(lambda value: f"{value:.6g}")
    lines = ["| " + " | ".join(columns) + " |",
             "| " + " | ".join("---" for _ in columns) + " |"]
    lines.extend("| " + " | ".join(map(str, row)) + " |"
                 for row in view.itertuples(index=False, name=None))
    return "\n".join(lines)


def write_plot(objectives, destination):
    if plt is None:
        return False
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.3))
    for ax, figure in zip(axes, ("fig1", "fig2", "fig3")):
        group = objectives[objectives.figure.eq(figure)]
        ax.axhline(1, color="black", ls="--", lw=1)
        ax.plot(group.complexity, group.full_cov_to_paper_ratio, "o-")
        ax.set_yscale("log"); ax.set_title(figure.upper())
        ax.set_xlabel("Complexity"); ax.set_ylabel("Recomputed / paper loss")
        ax.grid(True, which="both", alpha=.25)
    fig.tight_layout(); fig.savefig(destination, dpi=200); plt.close(fig)
    return True


def main():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    experiment_id = f"cross_figure_loss_audit_{stamp}_{uuid4().hex[:8]}"
    output = RESULTS_ROOT / experiment_id
    output.mkdir(parents=True, exist_ok=False)
    data = load_inputs(); bank = candidate_bank()
    objective_rows, rounding_rows, parameter_rows = [], [], []
    for item in bank:
        objective_rows.append(objective_row(item, data[item["figure"]]))
        scan = rounding_scan(item, data[item["figure"]])
        rounding_rows.append({"figure": item["figure"], "formula_id": item["formula_id"],
                              "complexity": item["complexity"], "paper_loss": item["paper_loss"],
                              **scan})
        for par in item["params"]:
            parameter_rows.append({"figure": item["figure"], "formula_id": item["formula_id"],
                                   "complexity": item["complexity"], **par})
    objectives = pd.DataFrame(objective_rows)
    rounding = pd.DataFrame(rounding_rows)
    parameters = pd.DataFrame(parameter_rows)
    selected_rows = []
    for item in selected_expression_variants():
        row = objective_row(item, data[item["figure"]])
        scan = rounding_scan(item, data[item["figure"]])
        selected_rows.append({"location": item["location"], **row, **scan})
    selected = pd.DataFrame(selected_rows)
    conclusive = objectives[~objectives.parse_status.str.startswith("minimal_repair")]
    summary = conclusive.groupby("figure").agg(
        candidates=("formula_id", "size"), ratio_min=("full_cov_to_paper_ratio", "min"),
        ratio_max=("full_cov_to_paper_ratio", "max"),
        exact_at_5pct=("full_cov_to_paper_ratio", lambda x: int(np.sum(np.abs(x-1)<=.05))))
    summary["rounding_boxes_bracketing_paper_loss"] = rounding.groupby("figure").paper_loss_bracketed.sum()
    objectives.to_csv(output / "objective_comparison.csv", index=False)
    rounding.to_csv(output / "rounding_box_scan.csv", index=False)
    parameters.to_csv(output / "rounding_parameter_intervals.csv", index=False)
    selected.to_csv(output / "selected_expression_crosschecks.csv", index=False)
    summary.to_csv(output / "summary_by_figure.csv")
    plot_written = write_plot(objectives, output / "loss_ratio_by_complexity.png")
    inputs = {key: {"path": str(value["path"].relative_to(ROOT)),
                    "sha256": hashlib.sha256(value["path"].read_bytes()).hexdigest(),
                    "n": len(value["z"])} for key, value in data.items()}
    record = {"experiment_id": experiment_id, "created_utc": datetime.now(timezone.utc).isoformat(),
              "paper": str(PAPER.relative_to(ROOT)), "paper_sha256": hashlib.sha256(PAPER.read_bytes()).hexdigest(),
              "inputs": inputs, "candidate_count": len(bank), "status": "complete",
              "rounding_search": "center + every corner + 32 deterministic interior samples + two 9-point coordinate sweeps for min and max; all values constrained to printed rounding intervals",
              "plot_written": plot_written}
    (output / "experiment.json").write_text(json.dumps(record, indent=2)+"\n", encoding="utf-8")
    closest_counts = conclusive.groupby(["figure", "closest_objective"]).size().rename("count").reset_index()
    closest_counts.to_csv(output / "closest_objective_counts.csv", index=False)
    metric_names = ["full_cov_chi2_over_n", "diagonal_cov_chi2_over_n",
                    "unweighted_sse", "unweighted_mse", "full_cov_chi2"]
    agreement_rows = []
    for figure, group in conclusive.groupby("figure"):
        for metric in metric_names:
            ratio = group[metric].to_numpy()/group.paper_loss.to_numpy()
            agreement_rows.append({"figure": figure, "objective": metric,
                                   "candidates": len(group),
                                   "within_5_percent": int(np.sum(np.abs(ratio-1) <= .05)),
                                   "within_10_percent": int(np.sum(np.abs(ratio-1) <= .10)),
                                   "within_factor_2": int(np.sum((ratio >= .5) & (ratio <= 2))),
                                   "median_log10_abs_ratio": float(np.median(np.abs(np.log10(ratio))))})
    agreement = pd.DataFrame(agreement_rows)
    agreement.to_csv(output / "objective_agreement_summary.csv", index=False)
    lines = ["# Cross-figure loss-objective and rounding audit", "",
             f"Experiment: `{experiment_id}`.", "", "## Figure-level summary", "",
             markdown_table(summary.reset_index(), list(summary.reset_index().columns)), "",
             "## Closest tested objective by candidate count", "",
             markdown_table(closest_counts, list(closest_counts.columns)), "",
             "## Findings", "",
             "- No single multiplicative normalization explains any of the three printed loss tables: the full-covariance loss ratio changes direction within Figures 2 and 3 and spans orders of magnitude in Figure 1.",
             "- Figure 1 complexities 2--4 remain close to unweighted SSE, whereas complexities 8--9 are close to the stated full-covariance loss. No one tested objective reproduces Table III uniformly.",
             "- The official arXiv LaTeX source has double minus signs in Figure 1 complexities 13 and 17. Complexity 15 is syntactically incomplete and is excluded from the rounding conclusion.",
             "- The official source for Table IV complexity 7 prints `0.76^2*z/0.03`, while Appendix A2 prints `0.76^z*z/0.03`; these are different formulas and must not be silently merged. Their separate checks are in `selected_expression_crosschecks.csv`.",
             "- A rounding box marked `paper_loss_bracketed=True` means the bounded numerical search found losses on both sides of the paper value while respecting every displayed decimal. A False result is failure to reach the value in this search, not a mathematical impossibility proof.",
             "- Remaining explanations that require author material are full-precision hall-of-fame expressions, a different fit-ready vector/covariance, a different loss implementation, or formula/loss rows taken from different runs.", "",
             "## Artifacts", "",
             "`objective_comparison.csv`, `objective_agreement_summary.csv`, `rounding_box_scan.csv`, `rounding_parameter_intervals.csv`, `selected_expression_crosschecks.csv`, `closest_objective_counts.csv`, `summary_by_figure.csv`, optional `loss_ratio_by_complexity.png`, and `experiment.json`."]
    (output / "findings.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
    print(output)
    print(summary.to_string())
    print("\nRounding boxes bracketing paper loss:")
    print(rounding.groupby("figure").paper_loss_bracketed.agg(["sum", "count"]).to_string())
    return output


if __name__ == "__main__":
    main()
