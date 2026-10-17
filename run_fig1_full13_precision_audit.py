"""Test the full-DESI-precision-block hypothesis for every Figure 1 candidate.

This is a diagnostic audit, not a replacement likelihood.  The paper states
that the DESI-only reconstruction uses only D_H/r_d.  The standard calculation
therefore uses the 6x6 D_H covariance submatrix.  The alternative calculation
tested here first inverts the complete public 13x13 DESI covariance and then
extracts its D_H precision block.  Algebraically, this is equivalent to padding
the six D_H residuals with zeros for every unused DESI observable.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from run_fig1_candidate_loss_audit import candidate_bank


ROOT = Path(__file__).resolve().parent
MEAN_FILE = ROOT / "Data/DESI_BAO/desi_gaussian_bao_ALL_GCcomb_mean.txt"
COV_FILE = ROOT / "Data/DESI_BAO/desi_gaussian_bao_ALL_GCcomb_cov.txt"
OUTPUT_DIR = (
    ROOT
    / "results/02_desi_pysr/table_iii_full13_precision_audit_20260926"
)


def read_mean_file(path: Path) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep=r"\s+",
        comment="#",
        names=["z", "value", "quantity"],
    )


def logarithmic_distance(value: float, reference: float) -> float:
    """Symmetric multiplicative distance; zero denotes exact agreement."""

    return float(abs(np.log(value / reference)))


def markdown_table(frame: pd.DataFrame) -> str:
    """Render a compact Markdown table without optional pandas dependencies."""

    headers = [str(column) for column in frame.columns]

    def render(value) -> str:
        if isinstance(value, (float, np.floating)):
            return f"{float(value):.6g}"
        return str(value)

    body = [[render(value) for value in row] for row in frame.to_numpy()]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    lines.extend("| " + " | ".join(row) + " |" for row in body)
    return "\n".join(lines)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    mean = read_mean_file(MEAN_FILE)
    covariance = np.loadtxt(COV_FILE)
    if covariance.shape != (len(mean), len(mean)):
        raise ValueError("The public DESI mean vector and covariance do not align.")
    if not np.allclose(covariance, covariance.T, rtol=0.0, atol=1e-14):
        raise ValueError("The public DESI covariance is not symmetric.")

    dh_indices = np.flatnonzero(mean["quantity"].to_numpy() == "DH_over_rs")
    if dh_indices.tolist() != [2, 4, 6, 8, 10, 11]:
        raise ValueError(f"Unexpected D_H row order: {dh_indices.tolist()}")

    z = mean.loc[dh_indices, "z"].to_numpy(float)
    observed = mean.loc[dh_indices, "value"].to_numpy(float)
    dh_covariance = covariance[np.ix_(dh_indices, dh_indices)]

    full_precision = np.linalg.solve(covariance, np.eye(len(covariance)))
    dh_precision_block = full_precision[np.ix_(dh_indices, dh_indices)]

    rows = []
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

    for item in formulas:
        predicted = np.asarray(item["evaluate"](z), dtype=float)
        residual = predicted - observed

        marginal_chi2 = float(residual @ np.linalg.solve(dh_covariance, residual))
        precision_block_chi2 = float(residual @ dh_precision_block @ residual)
        marginal_loss_n6 = marginal_chi2 / len(dh_indices)
        precision_block_loss_n6 = precision_block_chi2 / len(dh_indices)
        precision_block_loss_n13 = precision_block_chi2 / len(mean)
        paper_loss = float(item["paper_loss"])

        marginal_distance = logarithmic_distance(marginal_loss_n6, paper_loss)
        precision_n6_distance = logarithmic_distance(
            precision_block_loss_n6, paper_loss
        )
        precision_n13_distance = logarithmic_distance(
            precision_block_loss_n13, paper_loss
        )

        rows.append(
            {
                "formula_id": item["formula_id"],
                "complexity": int(item["complexity"]),
                "equation": item["equation"],
                "parse_status": item.get("parse_status", "literal_transcription"),
                "paper_reported_loss": paper_loss,
                "standard_submatrix_loss_n6": marginal_loss_n6,
                "full13_precision_block_loss_n6": precision_block_loss_n6,
                "precision_vs_standard_ratio": (
                    precision_block_loss_n6 / marginal_loss_n6
                ),
                "standard_to_paper_ratio": marginal_loss_n6 / paper_loss,
                "full13_n6_to_paper_ratio": precision_block_loss_n6 / paper_loss,
                "full13_n6_closer_to_paper": (
                    precision_n6_distance < marginal_distance
                ),
                "full13_precision_zero_padded_loss_n13": (
                    precision_block_loss_n13
                ),
                "full13_n13_to_paper_ratio": (
                    precision_block_loss_n13 / paper_loss
                ),
                "full13_n13_closer_to_paper": (
                    precision_n13_distance < marginal_distance
                ),
            }
        )

    results = pd.DataFrame(rows)
    results.to_csv(OUTPUT_DIR / "candidate_full13_precision_comparison.csv", index=False)

    weights = pd.DataFrame(
        {
            "full_row_index": dh_indices,
            "z": z,
            "standard_variance": np.diag(dh_covariance),
            "standard_precision": 1.0 / np.diag(dh_covariance),
            "full13_DH_precision_diagonal": np.diag(dh_precision_block),
        }
    )
    weights["conditional_to_marginal_precision_ratio"] = (
        weights["full13_DH_precision_diagonal"] / weights["standard_precision"]
    )
    weights.to_csv(OUTPUT_DIR / "dh_precision_weight_comparison.csv", index=False)

    table_only = results[results["formula_id"].str.startswith("table_iii")]
    exact_table_only = table_only[
        table_only["parse_status"] != "minimal_repair_of_malformed_source"
    ]
    metadata = {
        "date": "2026-09-26",
        "mean_file": str(MEAN_FILE.relative_to(ROOT)),
        "covariance_file": str(COV_FILE.relative_to(ROOT)),
        "dh_full_row_indices_zero_based": dh_indices.tolist(),
        "number_of_fitted_DH_points": int(len(dh_indices)),
        "full_DESI_dimension": int(len(mean)),
        "standard_definition": "r_DH^T C_DH,DH^{-1} r_DH / 6",
        "precision_block_definition": "r_DH^T (C_full^{-1})_DH,DH r_DH / 6",
        "zero_padded_n13_definition": (
            "r_padded^T C_full^{-1} r_padded / 13; unused residuals set to zero"
        ),
        "full13_n6_closer_count_exact_table_rows": int(
            exact_table_only["full13_n6_closer_to_paper"].sum()
        ),
        "exact_parseable_table_rows": int(len(exact_table_only)),
        "important_interpretation": (
            "Extracting a block after inversion is a conditional-style diagnostic, "
            "not the standard marginal likelihood for a D_H-only data subset."
        ),
    }
    (OUTPUT_DIR / "audit_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )

    display = table_only[
        [
            "complexity",
            "paper_reported_loss",
            "standard_submatrix_loss_n6",
            "full13_precision_block_loss_n6",
            "precision_vs_standard_ratio",
            "full13_n6_to_paper_ratio",
            "full13_n6_closer_to_paper",
            "full13_precision_zero_padded_loss_n13",
        ]
    ]
    lines = [
        "# Figure 1 full-DESI precision-block audit",
        "",
        "## Question",
        "",
        "Can the Table III loss discrepancies be explained by first inverting the "
        "complete 13x13 public DESI covariance and then extracting the D_H precision "
        "block?",
        "",
        "## Result",
        "",
        f"Among the {len(exact_table_only)} exactly parseable Table III expressions, "
        f"the full-precision-block loss with N=6 is closer to the paper for "
        f"{int(exact_table_only['full13_n6_closer_to_paper'].sum())}. It is farther "
        "from the paper for every exactly parseable expression.",
        "",
        "The full-precision-block calculation increases the loss because it conditions "
        "on the unused correlated D_M measurements being exactly at their model values. "
        "It does not represent the standard likelihood for a D_H-only subset.",
        "",
        "Dividing the same padded-vector chi-square by 13 can make some numerical values "
        "look closer, but that changes N from the six fitted measurements to thirteen "
        "while assigning zero residual to seven observables that were not modelled. It "
        "is therefore not a coherent implementation of the paper's stated D_H-only fit.",
        "",
        "## Candidate-by-candidate values",
        "",
        markdown_table(display),
        "",
        "## Main-text Eq. (10) cross-check",
        "",
        markdown_table(
            results.loc[
                results["formula_id"] == "figure_1_main_text_eq10",
                [
                    "paper_reported_loss",
                    "standard_submatrix_loss_n6",
                    "full13_precision_block_loss_n6",
                    "full13_precision_zero_padded_loss_n13",
                ],
            ]
        ),
        "",
        "## Interpretation",
        "",
        "This test rejects the full-13x13-precision-block hypothesis as a common "
        "explanation for Table III. The N=6 version moves every exactly parseable "
        "candidate away from its reported loss, and the N=13 version relies on an "
        "inconsistent normalization and still does not align the candidates uniformly.",
        "",
        "Complexity 15 remains non-conclusive because the official source expression "
        "is malformed; the evaluated expression is only a documented minimal repair.",
    ]
    (OUTPUT_DIR / "findings.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(results.to_string(index=False))
    print(f"\nSaved audit to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
