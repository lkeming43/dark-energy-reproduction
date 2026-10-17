"""Apply the arXiv v1 DESI pointwise objective to the v2 Table III formulas.

The v1 text defines the DESI loss per observation as

    abs(prediction - target) * weight,

with weight described as the inverse measurement uncertainty.  This audit
uses w_i = 1/sigma_i and reports both the mean aggregation normally used for
elementwise PySR losses and the unnormalised sum.  It is a diagnostic for a
possible legacy objective in the v2 tables, not a reproduction of the v1 DR1
analysis (which used different data).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

if __package__:
    from .run_fig1_candidate_loss_audit import (
        DATA_FILE, TABLE_I_SIGMA, TABLE_I_Y, candidate_bank,
    )
else:
    from run_fig1_candidate_loss_audit import (
        DATA_FILE, TABLE_I_SIGMA, TABLE_I_Y, candidate_bank,
    )


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = ROOT / "results/02_desi_pysr/table_iii_v1_objective_audit_20260926"


def multiplicative_distance(value: float, reference: float) -> float:
    return float(abs(np.log(value / reference)))


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    loaded = np.load(DATA_FILE)
    z = loaded["z"].astype(float)
    public_y = loaded["y"].astype(float)
    public_sigma = loaded["sigma"].astype(float)
    covariance = loaded["covariance"].astype(float)

    formulas = list(candidate_bank())
    formulas.append(
        {
            "formula_id": "figure_1_main_text_eq10",
            "complexity": 5,
            "paper_loss": 0.30,
            "equation": "28.4*exp(-0.519*z)",
            "evaluate": lambda x: 28.4 * np.exp(-0.519 * x),
            "parse_status": "main_text_equation; 0.30 is selected-row reference",
        }
    )

    rows = []
    for item in formulas:
        predicted = np.asarray(item["evaluate"](z), dtype=float)
        residual_public = predicted - public_y
        residual_rounded = predicted - TABLE_I_Y
        public_point_losses = np.abs(residual_public) / public_sigma
        rounded_point_losses = np.abs(residual_rounded) / TABLE_I_SIGMA
        covariance_chi2 = float(
            residual_public @ np.linalg.solve(covariance, residual_public)
        )

        paper_loss = float(item["paper_loss"])
        v1_mean = float(public_point_losses.mean())
        v1_sum = float(public_point_losses.sum())
        covariance_mean = covariance_chi2 / len(z)

        choices = {
            "v2_covariance_mean": covariance_mean,
            "v1_abs_inv_sigma_mean": v1_mean,
            "v1_abs_inv_sigma_sum": v1_sum,
        }
        closest_name = min(
            choices,
            key=lambda name: multiplicative_distance(choices[name], paper_loss),
        )

        rows.append(
            {
                "formula_id": item["formula_id"],
                "complexity": int(item["complexity"]),
                "equation": item["equation"],
                "parse_status": item.get("parse_status", "literal_transcription"),
                "paper_reported_loss": paper_loss,
                "v2_covariance_loss_mean": covariance_mean,
                "v1_public_abs_inv_sigma_mean": v1_mean,
                "v1_public_abs_inv_sigma_sum": v1_sum,
                "v1_rounded_table_abs_inv_sigma_mean": float(
                    rounded_point_losses.mean()
                ),
                "v1_mean_to_paper_ratio": v1_mean / paper_loss,
                "v1_sum_to_paper_ratio": v1_sum / paper_loss,
                "closest_of_tested_definitions": closest_name,
            }
        )

    result = pd.DataFrame(rows)
    result.to_csv(OUTPUT_DIR / "candidate_v1_objective_comparison.csv", index=False)

    table_rows = result[result["formula_id"].str.startswith("table_iii")]
    exact_rows = table_rows[
        table_rows["parse_status"] != "minimal_repair_of_malformed_source"
    ]
    v1_mean_closest_count = int(
        (exact_rows["closest_of_tested_definitions"] == "v1_abs_inv_sigma_mean").sum()
    )
    v1_sum_closest_count = int(
        (exact_rows["closest_of_tested_definitions"] == "v1_abs_inv_sigma_sum").sum()
    )
    v2_closest_count = int(
        (exact_rows["closest_of_tested_definitions"] == "v2_covariance_mean").sum()
    )

    metadata = {
        "date": "2026-09-26",
        "paper_versions": {
            "v1": "arXiv:2502.10506v1 (DESI DR1; submitted 2025-02-14)",
            "v2": "arXiv:2502.10506v2 (DESI DR2; submitted 2025-06-13)",
        },
        "diagnostic_question": (
            "Do the v2 Table III losses look as if the v1 DESI elementwise "
            "absolute-error objective was retained?"
        ),
        "v1_literal_point_loss": "abs(prediction-target)/sigma",
        "primary_aggregation": "mean over the six v2 D_H/r_d observations",
        "important_limitation": (
            "The v1 paper used DESI DR1 and a different Figure 1 analysis. This "
            "audit changes only the objective while retaining the v2 DR2 data and "
            "v2 Table III expressions."
        ),
        "exact_parseable_table_rows": int(len(exact_rows)),
        "closest_counts": {
            "v2_covariance_mean": v2_closest_count,
            "v1_abs_inv_sigma_mean": v1_mean_closest_count,
            "v1_abs_inv_sigma_sum": v1_sum_closest_count,
        },
    }
    (OUTPUT_DIR / "audit_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )

    columns = [
        "complexity",
        "paper_reported_loss",
        "v2_covariance_loss_mean",
        "v1_public_abs_inv_sigma_mean",
        "v1_public_abs_inv_sigma_sum",
        "closest_of_tested_definitions",
    ]
    lines = [
        "# Figure 1 audit using the arXiv v1 DESI objective",
        "",
        "The v1 pointwise definition is `abs(prediction-target)/sigma`. The mean "
        "and sum below apply that objective to the v2 DR2 data and Table III "
        "expressions. This is not a reproduction of the v1 DR1 analysis.",
        "",
        "```text",
        table_rows[columns].to_string(index=False),
        "```",
        "",
        "## Conclusion",
        "",
        f"Across {len(exact_rows)} exactly parseable Table III rows, the closest "
        f"tested definition is the v2 covariance mean for {v2_closest_count}, "
        f"the v1 absolute-error mean for {v1_mean_closest_count}, and the v1 "
        f"absolute-error sum for {v1_sum_closest_count}.",
        "",
        "A legacy v1 pointwise objective does not provide a single consistent "
        "explanation unless it aligns the reported values across the table, rather "
        "than improving isolated rows only.",
    ]
    (OUTPUT_DIR / "findings.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(result.to_string(index=False))
    print("\nClosest counts:", metadata["closest_counts"])
    print("Saved to", OUTPUT_DIR)


if __name__ == "__main__":
    main()
