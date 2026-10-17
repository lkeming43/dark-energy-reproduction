"""Notebook 11d: catalogue-level MC for the raw Notebook 11c PySR exponent.

The fixed expression is used as one mock-catalogue center. Each synthetic
catalogue is refitted within its exponential family, as in Notebook 11b.
"""

from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4
import hashlib
import json

import numpy as np
import pandas as pd

from .data_loading import load_standard_dataset
from .fig23_catalog_bootstrap import (
    WhitenedGLS, _bands, _figure, reconstruct_samples,
)
from .fig23_reconstruction import (
    compile_curve, gamma_value, w_from_derivatives,
)


def _hash(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def _selected_expression(root, selection_dir):
    summary_path = selection_dir / "summary.json"
    score_path = selection_dir / "candidate_subset_scores.csv"
    if not summary_path.is_file() or not score_path.is_file():
        raise FileNotFoundError("Completed Notebook 11c summary/scores are required")
    summary = json.loads(summary_path.read_text())
    selected = summary["notebook11c_selected"]
    scores = pd.read_csv(score_path)
    rows = scores.loc[scores.notebook11c_selected]
    if len(rows) != 1 or int(rows.iloc[0].seed) != selected["seed"]:
        raise ValueError("Notebook 11c selection metadata is inconsistent")
    row = rows.iloc[0]
    if row.family != "a*z*exp(b*z)" or row.expression != selected["expression"]:
        raise ValueError("Notebook 11c did not select the expected exponential family")
    prior_candidates = (root / "results/11_fig23_mc_and_reconstruction" /
        "fig23_mc_20261010T044029855521Z_5cf83467" /
        "fig2_pantheon/candidate_diagnostics.csv")
    source = pd.read_csv(prior_candidates)
    source_rows = source.loc[source.seed.eq(selected["seed"])]
    if len(source_rows) != 1 or source_rows.iloc[0].sympy_expression != row.expression:
        raise ValueError("Original PySR expression differs from Notebook 11c")
    source_row = source_rows.iloc[0]
    return row.expression, int(row.seed), float(source_row.exp_a), float(source_row.exp_b), \
        (summary_path, score_path, prior_candidates)


def _raw_curve(data, amplitude, slope, omega_m, rd0, grid):
    z = np.asarray(grid)
    exp_part = np.exp(slope * z)
    dm = amplitude * z * exp_part
    d1 = amplitude * exp_part * (1 + slope * z)
    d2 = amplitude * exp_part * (2 * slope + slope * slope * z)
    w, factor = w_from_derivatives(z, d1, d2, gamma_value(omega_m, rd0))
    if not np.isfinite(w).all() or not (factor > 0).all():
        raise ValueError("Raw PySR curve has a nonfinite value or nominal pole")
    return pd.DataFrame({"z": z, "distance": dm, "distance_first_derivative": d1,
                         "distance_second_derivative": d2, "w": w,
                         "pole_factor": factor})


def _checkpoints(bands, raw_curve):
    redshifts = (0.0, 0.1, 0.2, 0.5, 1.0, 1.5, 2.33)
    records = []
    for center in ("raw_selected_expression", "observed_data"):
        for selection in ("all", "pole_free"):
            values = bands.loc[(bands.center == center) &
                               (bands.selection == selection) &
                               (bands.quantity == "w")].reset_index(drop=True)
            for z in redshifts:
                index = int(np.argmin(np.abs(values.z.to_numpy() - z)))
                row = values.iloc[index]
                raw = raw_curve.iloc[index]
                records.append({"center": center, "selection": selection,
                                "z": float(row.z), "raw_central_w": float(raw.w),
                                "median_w": float(row["median"]),
                                "pointwise_3sigma_low": float(row.q0p135),
                                "pointwise_3sigma_high": float(row.q99p865),
                                "minus_one_inside_pointwise_3sigma": bool(
                                    row.q0p135 <= -1 <= row.q99p865)})
    return pd.DataFrame(records)


def _outside_segments(bands, center, selection):
    """Grid intervals where -1 misses the pointwise 3-sigma quantile band."""
    rows = bands.loc[(bands.center == center) &
                     (bands.selection == selection) &
                     (bands.quantity == "w")].sort_values("z")
    redshifts = rows.z.to_numpy()
    outside = ((rows.q0p135 > -1) | (rows.q99p865 < -1)).to_numpy()
    segments = []
    start = None
    for index, flag in enumerate(outside):
        if flag and start is None:
            start = index
        if not flag and start is not None:
            segments.append([float(redshifts[start]), float(redshifts[index - 1])])
            start = None
    if start is not None:
        segments.append([float(redshifts[start]), float(redshifts[-1])])
    return segments


def run_notebook11d(project_root, selection_directory=None, n_draws=20000,
                    batch_size=500, base_seed=20261016, output_root=None):
    """Run two catalogue-MC centers in one new, non-overwriting experiment folder."""
    root = Path(project_root).resolve()
    if n_draws < 100 or batch_size < 1:
        raise ValueError("Require at least 100 draws and a positive batch size")
    selection_dir = Path(selection_directory or root /
        "results/11c_fig2_balanced_seed_audit" /
        "fig2_balanced_20261016T044854959271Z_1a017e27").resolve()
    selection_manifest = json.loads((selection_dir / "experiment.json").read_text())
    if selection_manifest.get("status") != "complete":
        raise ValueError("Notebook 11c selection is not marked complete")
    expression, seed, amplitude, slope, selection_files = _selected_expression(root, selection_dir)
    data_path = root / "results/08_fig23_data_preparation/fig2_pantheon.npz"
    data = load_standard_dataset(data_path)
    problem = WhitenedGLS(data, "exponential")
    grid = np.linspace(0, 2.33, 301)
    np.testing.assert_allclose(problem.predict((amplitude, slope), grid),
                               compile_curve(expression)(grid)[0],
                               rtol=1e-13, atol=1e-12)
    fit = problem.fit_observed()
    if not np.allclose((amplitude, slope), fit[:2], rtol=2e-6, atol=2e-6):
        raise ValueError("Raw PySR coefficients are not sufficiently close to the GLS fit")
    raw_prediction = problem.predict((amplitude, slope))
    raw_bao = problem.whiten(data["y"][:6] - raw_prediction[:6], "bao")
    raw_sn = problem.whiten(data["y"][6:] - raw_prediction[6:], "sn")
    raw_loss = float((raw_bao @ raw_bao + raw_sn @ raw_sn) / len(data["y"]))
    cosmos = data["metadata"]["config"]["cosmology"]
    rd0, rd_sigma, omega_m = (float(cosmos[key]) for key in
                              ("rd_hinv_mpc", "rd_sigma_hinv_mpc", "Omega_m"))
    raw_curve = _raw_curve(data, amplitude, slope, omega_m, rd0, grid)
    output_root = Path(output_root or root / "results/11d_fig2_selected_exponential_mc")
    output_root.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    destination = output_root / f"fig2_exp_mc_{stamp}_{uuid4().hex[:8]}"
    destination.mkdir(exist_ok=False)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "in_progress", "scope": "Figure 2 only",
        "selected_pysr_seed": seed, "raw_pysr_expression": expression,
        "equivalent_exponential_parameters": {"amplitude": amplitude, "slope": slope},
        "coefficient_policy": "raw PySR coefficients fixed as the model-centered mock mean; each mock catalogue refitted in the same exponential family",
        "method": "full-covariance Gaussian catalogue bootstrap, 20k draws per center by default; exact shared r_d SN rescaling and per-draw GLS refit",
        "centers": ["raw_selected_expression", "observed_data"],
        "selection_policy": "structure fixed from Notebook 11c; no per-draw PySR reselection",
        "pointwise_intervals": True,
        "input_sha256": {str(path): _hash(path) for path in
                         (data_path, selection_dir / "experiment.json", *selection_files)},
        "n_draws_per_center": int(n_draws), "batch_size": int(batch_size),
        "base_seed": int(base_seed),
        "limitations": [
            "The result is exploratory because 11c re-selection followed Notebook 12 diagnostics.",
            "Only one of 100 seed-selected expressions has this exponential family.",
            "Intervals are conditional on the fixed exponential structure and are pointwise, not simultaneous.",
            "Pole-free bands condition on discarding pole-crossing realizations."]}
    _write_json(destination / "experiment.json", manifest)
    try:
        raw_curve.to_csv(destination / "raw_central_curve.csv", index=False)
        centers = {"raw_selected_expression": raw_prediction,
                   "observed_data": problem.y}
        all_bands = []
        overview = []
        for center_index, (name, center_vector) in enumerate(centers.items()):
            rng_rows = []
            count_batches = (n_draws + batch_size - 1) // batch_size
            for batch_index in range(count_batches):
                count = min(batch_size, n_draws - batch_index * batch_size)
                sequence = np.random.SeedSequence(base_seed, spawn_key=(center_index, batch_index))
                rng = np.random.default_rng(sequence)
                rd = rng.normal(rd0, rd_sigma, size=count)
                fits = problem.draw_and_fit(center_vector, rd, rd0, rng)
                rng_rows.append(np.column_stack((fits, rd)))
                if (batch_index + 1) % 10 == 0 or batch_index + 1 == count_batches:
                    print(f"{name}: {batch_index + 1}/{count_batches} batches", flush=True)
            result = np.concatenate(rng_rows, axis=0)
            a, b, chi2, rd = result.T
            dm, w, finite, pole, monotonic = reconstruct_samples(
                "exponential", a, b, rd, grid, omega_m)
            records = pd.DataFrame({"draw_index": np.arange(n_draws),
                                    "amplitude": a, "slope": b,
                                    "chi2": chi2, "rd_hinv_mpc": rd,
                                    "finite": finite, "pole_crossing": pole,
                                    "distance_increasing": monotonic})
            records.to_csv(destination / f"{name}_samples.csv.gz", index=False,
                           compression="gzip")
            subsets = []
            for selection, mask in (("all", finite), ("pole_free", finite & ~pole)):
                if int(mask.sum()) < 100:
                    raise ValueError(f"Too few {name}/{selection} samples")
                subsets.append(_bands(grid, dm[mask], "distance", name, selection))
                subsets.append(_bands(grid, w[mask], "w", name, selection))
            bands = pd.concat(subsets, ignore_index=True)
            bands.to_csv(destination / f"{name}_bands.csv", index=False)
            all_bands.append(bands)
            diagnostics = {"center": name, "n_draws": int(n_draws),
                           "n_finite": int(finite.sum()),
                           "n_pole_crossing": int((finite & pole).sum()),
                           "n_pole_free": int((finite & ~pole).sum()),
                           "n_nonincreasing": int((finite & ~monotonic).sum())}
            overview.append(diagnostics)
            _write_json(destination / f"{name}_diagnostics.json", diagnostics)
            _figure(destination / f"{name}_reconstruction.png", data, grid,
                    "11c selected exponential", name, bands,
                    data["metadata"]["config"]["references"]["paper_expression"],
                    omega_m, rd0)
            del dm, w, result
        all_bands = pd.concat(all_bands, ignore_index=True)
        all_bands.to_csv(destination / "all_bands.csv", index=False)
        checkpoints = _checkpoints(all_bands, raw_curve)
        checkpoints.to_csv(destination / "w_checkpoints.csv", index=False)
        summary = {"selected_seed": seed, "expression": expression,
                   "raw_amplitude": amplitude, "raw_slope": slope,
                   "raw_combined_loss": raw_loss,
                   "gls_refitted_amplitude": float(fit[0]),
                   "gls_refitted_slope": float(fit[1]),
                   "gls_refitted_loss": float(fit[2] / len(data["y"])),
                   "rd_mean": rd0, "rd_sigma": rd_sigma, "omega_m": omega_m,
                   "centers": overview,
                   "pointwise_3sigma_outside_minus_one_segments": {
                       name: {selection: _outside_segments(all_bands, name, selection)
                              for selection in ("all", "pole_free")}
                       for name in centers}}
        _write_json(destination / "summary.json", summary)
        lines = ["# Notebook 11d — MC for the Figure 2 seed-20261093 exponential", "",
                 "The raw PySR expression is the source of the fixed model structure and "
                 "the model-centered mock-catalogue mean. Each Gaussian mock catalogue "
                 "is refitted within the exponential family using the full prepared covariance.", "",
                 f"Raw expression: `{expression}`", "",
                 f"Equivalent form: `D_M/r_d = {amplitude:.12f} z exp({slope:.12f} z)`", "",
                 f"Raw-expression combined loss: {raw_loss:.6f}; GLS refitted-family "
                 f"loss: {fit[2]/len(data['y']):.6f}.", "",
                 "| Mock centre | Draws | Finite | Pole crossing | Pole-free |",
                 "| --- | ---: | ---: | ---: | ---: |"]
        for item in overview:
            lines.append(f"| {item['center']} | {item['n_draws']} | "
                         f"{item['n_finite']} | {item['n_pole_crossing']} | "
                         f"{item['n_pole_free']} |")
        lines += ["", "Approximate grid ranges where `w=-1` falls outside the "
                  "pointwise 3-sigma interval:"]
        for name in centers:
            for selection in ("all", "pole_free"):
                segments = summary["pointwise_3sigma_outside_minus_one_segments"][name][selection]
                lines.append(f"- {name}, {selection}: {segments or 'none'}")
        lines += ["", "The two centres answer different questions: simulated data "
                  "around the raw selected expression versus around the observed catalogue. "
                  "Neither re-selects a PySR structure for each draw.", "",
                  "All bands are pointwise. The pole-free bands are conditional on the "
                  "filter and cannot replace the unfiltered result. Check `w_checkpoints.csv` "
                  "for selected redshifts; `all_bands.csv` contains the full grid.", "",
                  "At z below the lowest DESI point (0.51), DESI-only support is an "
                  "extrapolation; supernova data do extend to lower redshift.", "",
                  "This is an exploratory follow-up to Notebook 11c, not evidence that "
                  "the printed paper expression or reported loss has been reproduced.", ""]
        (destination / "analysis.md").write_text("\n".join(lines))
        manifest["status"] = "complete"
        _write_json(destination / "experiment.json", manifest)
        return destination, summary
    except Exception:
        manifest["status"] = "failed"
        _write_json(destination / "experiment.json", manifest)
        raise
