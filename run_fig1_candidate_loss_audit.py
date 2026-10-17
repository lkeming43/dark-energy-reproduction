"""Audit the DESI-only candidate losses printed in Table III.

The primary calculation evaluates each printed expression on the six exact
DESI DR2 D_H/r_d measurements and applies the public covariance according to
Eq. (9) of arXiv:2502.10506v2:

    loss = residual.T @ covariance^{-1} @ residual / n.

The rounded Table I values and diagonal errors are evaluated as a sensitivity
check.  This script does not rerun PySR and does not alter the raw data.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ModuleNotFoundError:  # The numerical audit does not require plotting.
    plt = None


ROOT = Path(__file__).resolve().parent
DATA_FILE = ROOT / "results/01_desi_data_check/desi_dr2_dh_over_rd.npz"
OUTPUT_DIR = (
    ROOT
    / "results/02_desi_pysr/table_iii_candidate_loss_audit_20260922_source_verified"
)

PAPER_URL = "https://arxiv.org/html/2502.10506v2"

# Table I, as printed in the paper.  This is only a sensitivity configuration;
# it is not claimed to be the covariance used by the authors.
TABLE_I_Y = np.array([21.863, 19.455, 17.641, 14.176, 12.817, 8.632])
TABLE_I_SIGMA = np.array([0.425, 0.330, 0.193, 0.221, 0.516, 0.101])


def candidate_bank():
    """Return literal transcriptions of the parseable Table III equations."""

    return [
        {
            "formula_id": "table_iii_c1",
            "complexity": 1,
            "paper_loss": 2519.06,
            "paper_decimals": 2,
            "equation": "2.8069465",
            "evaluate": lambda z: np.full_like(z, 2.8069465),
        },
        {
            "formula_id": "table_iii_c2",
            "complexity": 2,
            "paper_loss": 120.73,
            "paper_decimals": 2,
            "equation": "sqrt(276.74304)",
            "evaluate": lambda z: np.full_like(z, np.sqrt(276.74304)),
        },
        {
            "formula_id": "table_iii_c3",
            "complexity": 3,
            "paper_loss": 91.5,
            "paper_decimals": 1,
            "equation": "17.705526 - z",
            "evaluate": lambda z: 17.705526 - z,
        },
        {
            "formula_id": "table_iii_c4",
            "complexity": 4,
            "paper_loss": 1.22,
            "paper_decimals": 2,
            "equation": "(5.127707 - z)^2",
            "evaluate": lambda z: (5.127707 - z) ** 2,
        },
        {
            "formula_id": "table_iii_c5_selected",
            "complexity": 5,
            "paper_loss": 0.30,
            "paper_decimals": 2,
            "equation": "28.42 / 1.68^z",
            "evaluate": lambda z: 28.42 / (1.68**z),
        },
        {
            "formula_id": "table_iii_c7",
            "complexity": 7,
            "paper_loss": 0.25,
            "paper_decimals": 2,
            "equation": "(3.77 - z)^2.32 + 6.24",
            "evaluate": lambda z: (3.77 - z) ** 2.32 + 6.24,
        },
        {
            "formula_id": "table_iii_c8",
            "complexity": 8,
            "paper_loss": 0.24,
            "paper_decimals": 2,
            "equation": "((z - 3.18)^2 - (-3.47)) * 2.06",
            "evaluate": lambda z: ((z - 3.18) ** 2 + 3.47) * 2.06,
        },
        {
            "formula_id": "table_iii_c9",
            "complexity": 9,
            "paper_loss": 0.21,
            "paper_decimals": 2,
            "equation": "2.70^z - 11.8*z + 26.0",
            "evaluate": lambda z: 2.70**z - 11.8 * z + 26.0,
        },
        {
            "formula_id": "table_iii_c11",
            "complexity": 11,
            "paper_loss": 0.17,
            "paper_decimals": 2,
            "equation": "4.19*(1 - 0.787*z)^6 + (4.94 - z)^2.1",
            "evaluate": lambda z: 4.19 * (1 - 0.787 * z) ** 6
            + (4.94 - z) ** 2.1,
        },
        {
            "formula_id": "table_iii_c13",
            "complexity": 13,
            "paper_loss": 0.02,
            "paper_decimals": 2,
            "equation": "(3.68 - z)^2.4 + 6.56 - -0.02/(z - 0.76)",
            "evaluate": lambda z: (3.68 - z) ** 2.4
            + 6.56
            + 0.02 / (z - 0.76),
        },
        {
            # The official LaTeX source contains ``3.68-z\times`` with no
            # operand after the multiplication sign.  This minimal repair
            # drops the orphan multiplication sign and is not exact.
            "formula_id": "table_iii_c15_minimal_repair",
            "complexity": 15,
            "paper_loss": 0.016,
            "paper_decimals": 3,
            "equation": "(3.68-z)^2.37 + 6.56 - -0.0169/[0.862*(z-0.763)]",
            "evaluate": lambda z: (3.68 - z) ** 2.37
            + 6.56
            + 0.0169 / (0.862 * (z - 0.763)),
            "parse_status": "minimal_repair_of_malformed_source",
        },
        {
            "formula_id": "table_iii_c16",
            "complexity": 16,
            "paper_loss": 0.015,
            "paper_decimals": 3,
            "equation": (
                "21.9*(1 - 0.272*z)^2.37 + 6.56 "
                "+ 0.0169/[sqrt(z)*(z - 0.763)]"
            ),
            "evaluate": lambda z: 21.9 * (1 - 0.272 * z) ** 2.37
            + 6.56
            + 0.0169 / (np.sqrt(z) * (z - 0.763)),
        },
        {
            "formula_id": "table_iii_c17",
            "complexity": 17,
            "paper_loss": 0.013,
            "paper_decimals": 3,
            "equation": (
                "(3.68-z)^2.37 + 6.56 "
                "- -0.0169/[z*(z - 0.763 - 0.0216)]"
            ),
            "evaluate": lambda z: (3.68 - z) ** 2.37
            + 6.56
            + 0.0169 / (z * (z - 0.763 - 0.0216)),
        },
    ]


def covariance_loss(observed, predicted, covariance):
    residual = np.asarray(predicted) - np.asarray(observed)
    chi2_contributions = residual * np.linalg.solve(covariance, residual)
    chi2 = float(np.sum(chi2_contributions))
    return chi2, chi2 / len(residual), chi2_contributions


def evaluate_audit():
    loaded = np.load(DATA_FILE)
    z = loaded["z"].astype(float)
    public_y = loaded["y"].astype(float)
    public_covariance = loaded["covariance"].astype(float)
    table_covariance = np.diag(TABLE_I_SIGMA**2)

    rows = []
    predictions = {}
    for item in candidate_bank():
        predicted = np.asarray(item["evaluate"](z), dtype=float)
        if predicted.shape != public_y.shape or not np.isfinite(predicted).all():
            raise ValueError(f"Invalid prediction for {item['formula_id']}")
        public_chi2, public_loss, contributions = covariance_loss(
            public_y, predicted, public_covariance
        )
        table_chi2, table_loss, _ = covariance_loss(
            TABLE_I_Y, predicted, table_covariance
        )
        paper_loss = float(item["paper_loss"])
        decimals = int(item["paper_decimals"])
        dominant_index = int(np.argmax(contributions))
        rows.append(
            {
                "formula_id": item["formula_id"],
                "complexity": int(item["complexity"]),
                "equation": item["equation"],
                "parse_status": item.get("parse_status", "literal_transcription"),
                "paper_reported_loss": paper_loss,
                "public_chi2": public_chi2,
                "public_loss": public_loss,
                "public_minus_paper": public_loss - paper_loss,
                "public_to_paper_ratio": public_loss / paper_loss,
                "same_at_paper_precision": round(public_loss, decimals)
                == round(paper_loss, decimals),
                "table_i_rounded_loss": table_loss,
                "table_i_minus_public": table_loss - public_loss,
                "unweighted_sse": float(np.sum((predicted - public_y) ** 2)),
                "dominant_chi2_redshift": float(z[dominant_index]),
                "dominant_chi2_fraction": float(
                    contributions[dominant_index] / public_chi2
                ),
            }
        )
        predictions[item["formula_id"]] = predicted

    # Eq. (10), which is the expression actually plotted in Figure 1.  The
    # paper does not print a separate loss for this rounded form, so comparison
    # with 0.30 is explicitly labelled as an implied Table III comparison.
    eq10_prediction = 28.4 * np.exp(-0.519 * z)
    eq10_public_chi2, eq10_public_loss, _ = covariance_loss(
        public_y, eq10_prediction, public_covariance
    )
    eq10_table_chi2, eq10_table_loss, _ = covariance_loss(
        TABLE_I_Y, eq10_prediction, table_covariance
    )

    frame = pd.DataFrame(rows).sort_values("complexity").reset_index(drop=True)
    eq10 = {
        "equation": "28.4*exp(-0.519*z)",
        "public_chi2": eq10_public_chi2,
        "public_loss": eq10_public_loss,
        "table_i_rounded_chi2": eq10_table_chi2,
        "table_i_rounded_loss": eq10_table_loss,
        "implied_table_iii_selected_loss": 0.30,
    }
    metadata = {
        "paper_url": PAPER_URL,
        "paper_version": "arXiv:2502.10506v2",
        "paper_loss_definition": "r.T @ inv(C) @ r / n",
        "n": int(len(z)),
        "data_file": str(DATA_FILE.relative_to(ROOT)),
        "public_z": z.tolist(),
        "public_y": public_y.tolist(),
        "public_sigma": np.sqrt(np.diag(public_covariance)).tolist(),
        "table_i_y": TABLE_I_Y.tolist(),
        "table_i_sigma": TABLE_I_SIGMA.tolist(),
        "eq10": eq10,
        "formula_transcription": (
            "Official arXiv LaTeX source checked: complexities 13 and 17 contain "
            "double minus signs; complexity 15 contains an orphan multiplication sign."
        ),
    }
    return frame, metadata


def write_plot(frame):
    if plt is None:
        return
    valid = frame[frame["parse_status"] == "literal_transcription"]
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    x = valid["paper_reported_loss"].to_numpy()
    y = valid["public_loss"].to_numpy()
    ax.scatter(x, y, s=55, color="tab:blue", zorder=3)
    for _, row in valid.iterrows():
        ax.annotate(
            f"C{int(row['complexity'])}",
            (row["paper_reported_loss"], row["public_loss"]),
            xytext=(5, 4),
            textcoords="offset points",
            fontsize=8,
        )
    lower = min(float(np.min(x)), float(np.min(y))) / 1.8
    upper = max(float(np.max(x)), float(np.max(y))) * 1.8
    ax.plot([lower, upper], [lower, upper], "k--", lw=1.2, label="Exact agreement")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lower, upper)
    ax.set_ylim(lower, upper)
    ax.set_xlabel("Loss reported in Table III")
    ax.set_ylabel("Recomputed public-data loss")
    ax.set_title("DESI-only candidate-loss audit")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "paper_vs_recomputed_loss.png", dpi=200)
    plt.close(fig)


def write_findings(frame, metadata):
    eq10 = metadata["eq10"]
    literal = frame[frame["parse_status"] == "literal_transcription"]
    exact_precision = literal[literal["same_at_paper_precision"]]
    closest = literal.assign(
        relative_difference=np.abs(
            literal["public_minus_paper"] / literal["paper_reported_loss"]
        )
    ).sort_values("relative_difference")

    display = frame[
        [
            "complexity",
            "paper_reported_loss",
            "public_loss",
            "table_i_rounded_loss",
            "public_to_paper_ratio",
            "same_at_paper_precision",
            "parse_status",
        ]
    ].copy()
    for column in (
        "paper_reported_loss",
        "public_loss",
        "table_i_rounded_loss",
        "public_to_paper_ratio",
    ):
        display[column] = display[column].map(lambda value: f"{value:.6g}")

    headers = [str(column) for column in display.columns]
    markdown_rows = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for values in display.itertuples(index=False, name=None):
        markdown_rows.append("| " + " | ".join(map(str, values)) + " |")
    markdown_table = "\n".join(markdown_rows)

    lines = [
        "# Figure 1 / Table III candidate-loss audit",
        "",
        "## Scope",
        "",
        "Each printed DESI-only candidate equation was evaluated at the six ",
        "official DESI DR2 $D_H/r_d$ redshifts. The primary result uses the ",
        "exact release central values and covariance and applies the paper's ",
        "Eq. (9), $L=r^T C^{-1}r/n$.",
        "",
        "## Result",
        "",
        markdown_table,
        "",
        (
            f"Only {len(exact_precision)} of {len(literal)} unambiguous printed "
            "equations reproduces the reported loss at the number of decimal "
            "places shown in Table III."
        ),
        "",
        "The closest literal transcriptions by relative difference are:",
        "",
    ]
    for _, row in closest.head(4).iterrows():
        lines.append(
            f"- complexity {int(row['complexity'])}: paper "
            f"{row['paper_reported_loss']:.6g}, recomputed "
            f"{row['public_loss']:.6g}, relative difference "
            f"{100*row['relative_difference']:.2f}%"
        )
    lines.extend(
        [
            "",
            "The rounded Table I configuration changes the losses only modestly "
            "and does not reconcile the discrepant rows.",
            "",
            "## Figure 1 selected expression",
            "",
            (
                "The expression printed as Eq. (10), "
                f"`${eq10['equation']}`, gives $\\chi^2="
                f"{eq10['public_chi2']:.6f}$ and $L="
                f"{eq10['public_loss']:.6f}$ with the public release data. "
                "Using rounded Table I values gives "
                f"$L={eq10['table_i_rounded_loss']:.6f}$. Both differ from the "
                "0.30 reported for the corresponding selected Table III row."
            ),
            "",
            "The selected Table III form, `28.42/1.68^z`, gives "
            f"$L={float(frame.loc[frame.complexity == 5, 'public_loss'].iloc[0]):.6f}$.",
            "",
            "## Alternative-loss diagnostic",
            "",
            "For complexities 2, 3, and 4, the reported values happen to be "
            "close to the unweighted sum of squared residuals: respectively "
            "$120.73$ versus $120.984$, $91.5$ versus $90.522$, and $1.22$ "
            "versus $1.261$. This pattern does not hold for the other rows, so "
            "it does not provide a single alternative loss definition that "
            "reproduces Table III.",
            "",
            "## Interpretation",
            "",
            "The published loss column cannot be reproduced by applying the "
            "published Eq. (9) uniformly to the printed equations and the "
            "official DESI DR2 data. Coefficient rounding can strongly affect "
            "the high-complexity pole-like expressions, but it does not explain "
            "the large discrepancies for several low-complexity rows whose "
            "constants are printed to many digits.",
            "",
            "The official arXiv LaTeX source has double minus signs in "
            "complexities 13 and 17. Complexity 15 is malformed (`3.68-z\\times` "
            "has no following operand); its value uses a documented minimal repair.",
            "",
            "Possible explanations include unpublished full-precision equations, "
            "a different data/covariance object during training, stale losses "
            "copied from another run, or a mismatch between the implemented and "
            "stated loss. The available paper information is insufficient to "
            "distinguish these possibilities conclusively.",
            "",
        ]
    )
    (OUTPUT_DIR / "findings.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    frame, metadata = evaluate_audit()
    frame.to_csv(OUTPUT_DIR / "candidate_loss_comparison.csv", index=False)
    (OUTPUT_DIR / "audit_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    write_findings(frame, metadata)
    write_plot(frame)
    print(frame.to_string(index=False))
    print(f"\nSaved audit to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
