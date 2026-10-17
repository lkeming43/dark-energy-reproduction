# Reproduction Log

## Project Goal

Independently reproduce the symbolic-regression reconstruction of cosmological distance functions and the dark-energy equation of state using DESI DR2, Pantheon+, and DESY5 data.

## Progress Checklist

- [ ] Verify Equations (1)-(6)
- [ ] Validate DESI DR2 Table I data
- [ ] Reproduce the DESI-only PySR fit
- [ ] Reproduce Figure 1
- [ ] Reproduce DESI + Pantheon+ and Figure 2
- [ ] Reproduce DESI + DESY5 and Figure 3
- [ ] Reconstruct confidence intervals
- [ ] Test multiple random seeds
- [ ] Document discrepancies and ambiguities
- [ ] Reproduce Appendix B parameter tests

---

## 2026-08-23 — Project Setup

### Objective

Create a reproducible project structure and define the first reproduction target.

### Current Data

- DESI DR2 BAO mean and covariance files
- Pantheon+ data and covariance file
- DESY5 Hubble-diagram data and covariance file

### First Target

Validate the DESI DR2 data and reproduce the left panel of Figure 1 before running symbolic regression.

### Next Step

Open `notebooks/01_desi_data_check.ipynb`, inspect the DESI data format, and compare the values with Table I of the paper.

---

## Experiment Template

### Date and Experiment ID

Date:

Experiment ID:

### Objective

What specific question does this run test?

### Input Data

- Dataset:
- Source file:
- Number of data points:
- Redshift range:
- Cleaning or conversion steps:
- Covariance treatment:

### Configuration

- Python version:
- PySR version:
- Random seed:
- Unary operators:
- Binary operators:
- Loss function:
- Model-selection criterion:
- Iterations and population settings:
- Physical constraints:
- Cosmological parameters:

### Changes From the Previous Run

What changed and why?

### Expected Result

What does the paper report?

### Actual Result

- Selected expression:
- Loss:
- Complexity:
- Runtime:
- Output files:

### Comparison and Interpretation

- What agrees with the paper?
- What differs?
- Is the difference scientifically meaningful?
- What are the possible causes?

### Validation

- [ ] Units checked
- [ ] Source data checked
- [ ] Covariance matrix checked
- [ ] Derivatives checked
- [ ] No unphysical poles in the fitted range
- [ ] Run can be reproduced from saved configuration

### AI Assistance

Record what AI helped with and how the result was manually verified.

### Conclusion and Next Step

Conclusion:

Next step:

---

## 2026-09-21 — Figure 1 / Table III Candidate-Loss Audit

### Objective

Check whether the losses printed for the DESI-only candidate equations in Table III can be recovered by evaluating those equations on the official DESI DR2 $D_H/r_d$ data with the loss stated in Eq. (9),

$$
\mathcal{L}=\frac{1}{n}\mathbf{r}^{T}\mathbf{C}^{-1}\mathbf{r}.
$$

### Method

- Used the six exact public DESI DR2 $D_H/r_d$ central values and the corresponding public covariance submatrix.
- Evaluated every unambiguous Table III expression literally at the six published redshifts.
- Repeated the calculation with the rounded Table I central values and diagonal errors as a sensitivity test.
- Did not rerun PySR.
- Treated the complexity-15 expression separately because the printed equation is syntactically incomplete.

### Main Results

- Only the complexity-9 expression reproduces the paper loss at the displayed two-decimal precision: paper $0.21$, recomputed $0.207807$.
- Complexity 8 is close but not identical: paper $0.24$, recomputed $0.248055$.
- The selected Table III expression, $28.42/1.68^z$, gives $\mathcal{L}=0.664026$, rather than the reported $0.30$.
- The rounded Figure 1 / Eq. (10) expression, $28.4\exp(-0.519z)$, gives $\chi^2=4.222491$ and $\mathcal{L}=0.703749$.
- Several rows differ by factors of approximately 4–65. The low-complexity 2–4 discrepancies are dominated by the $z=2.33$ Ly$\alpha$ measurement.
- Replacing the public covariance errors with the rounded Table I errors changes the results only modestly and does not reconcile the table.
- The reported losses for complexities 2–4 are close to the unweighted sums of squared residuals, but this pattern does not hold for the remaining rows and therefore is not a consistent alternative definition for the table.

### Interpretation

The Table III loss column is not reproducible by applying the published Eq. (9) uniformly to the printed candidate equations and the official DESI DR2 data. Rounded coefficients may explain part of the disagreement for pole-like high-complexity expressions, but they cannot explain all low-complexity discrepancies. Possible causes include unpublished full-precision expressions, a different training data/covariance object, stale values from another run, or a difference between the implemented and documented objective.

### Outputs

- `results/02_desi_pysr/table_iii_candidate_loss_audit_20260921/candidate_loss_comparison.csv`
- `results/02_desi_pysr/table_iii_candidate_loss_audit_20260921/findings.md`
- `results/02_desi_pysr/table_iii_candidate_loss_audit_20260921/paper_vs_recomputed_loss.png`
- `scripts/audits/run_fig1_candidate_loss_audit.py`

---

## 2026-09-21 — DESI DR2 Input-Data Provenance Check

### Result

- The paper's original v1 used DESI DR1. The DR2 analysis first appears in v2, submitted on 13 June 2025, after the DESI DR2 public release in March 2025. No evidence was found that the authors needed or used pre-release DR2 files.
- The six Table I central values exactly equal the current public DESI DR2 mean values rounded to three decimal places.
- The public mean and covariance files were added to `CobayaSampler/bao_data` in a single commit on 20 March 2025 and have no later changes in that repository. The project copies are byte-for-byte identical to the current public files.
- Table I exactly reproduces DESI DR2 Results II v2 Table 4. Those values are presentation-level marginalized-posterior means and standard deviations, not a complete covariance specification.
- The authors did not publish the exact DESI input file, hash, or code. Use of full-precision data is plausible but cannot be independently proved.
- Extra input precision or ordinary Table I rounding does not account for the previously measured Table III loss discrepancies.

### Detailed Record

- `results/02_desi_pysr/table_iii_candidate_loss_audit_20260921/data_provenance_check.md`

## 2026-09-21 — Table III Formula-Transcription Check

- Inspected the official arXiv v2 LaTeX source directly.
- Confirmed that repeated minus signs and the malformed complexity-15 `z\times` term are present in the author source and are not OCR artifacts.
- Verified that interpreting the repeated minus as subtraction of a negative is the only plausible reading; the opposite sign produces much worse losses.
- Confirmed that formula parsing cannot explain the selected complexity-5 loss discrepancy because that expression is unambiguous, and even the best refitted exponential has loss about 0.407 rather than the reported 0.30.
- Detailed record: `results/02_desi_pysr/table_iii_candidate_loss_audit_20260921/formula_transcription_check.md`.

---

## 2026-10-10 — Notebook 11 Figure 2/3 Reconstruction

- Implemented a documented independent selection rule for the completed Notebook 10 public-covariance 100-seed results. The paper's printed and refitted exponential expressions remain separate comparisons.
- Reconstructed $D_M/r_d$ and $w(z)$ for Figures 2 and 3, with 10,000 first-order GLS coefficient Monte Carlo draws and one shared $r_d$ draw per realization. This is a conditional, linearized uncertainty analysis, **not** a full per-draw PySR rerun or a full data-level bootstrap.
- The independent Figure 2 rule selected a quadratic expression (seed 20261064, original loss 0.920080); Figure 3 selected an exponential expression (seed 20261023, original loss 0.937176). There were 29 physically eligible seed-selected expressions in each scenario under the stated grid checks.
- The primary branches retained 9,700/10,000 and 9,772/10,000 pole-free realizations respectively. Wide high-redshift $w(z)$ tails and structure sensitivity mean these outputs cannot be used as a global rejection of $w=-1$.
- Full outputs and the analysis log: `results/11_fig23_mc_and_reconstruction/fig23_mc_20261010T044029855521Z_5cf83467/`.
- Verified the notebook schema, executed every notebook code cell through the runner, inspected both saved figures, and passed all 14 project tests. The source-loss mismatch with Tables III–V remains unresolved.

---

## 2026-10-14 — Separate Notebook 11b Full-Catalogue Monte Carlo Experiment

- Kept Notebook 11 untouched and froze its three selected/comparator model-family branches. No PySR expression was reselected based on its $w(z)$ outcome.
- Implemented a data-level Gaussian Monte Carlo: draw correlated mock DESI + supernova catalogues from the prepared full covariance; draw one shared $r_d$ per catalogue; refit the fixed quadratic or exponential family using GLS; reconstruct $w(z)$; report both all-valid and pole-free pointwise bands.
- Ran 20,000 realizations for each branch under two centring choices (best-fit expression and observed data), totalling 120,000 fixed-family refits. The run has its own timestamped experiment folder, with per-branch batch files, fit samples, diagnostics, bands, and plots.
- The catalogue-refit medians and 95% bands closely match the previous linearized result at representative redshifts. The earlier high-redshift pole and model-family sensitivity remain. At $z=0$, the pole-free 3-sigma interval includes $w=-1$ for the Figure 2 quadratic but excludes it for the Figure 2 and Figure 3 exponentials; this is not a global LambdaCDM test.
- The new run coherently rescales the SN covariance with $r_d$, whereas Notebook 11 used a first-order coefficient response at nominal covariance. Thus the comparison does not isolate linearization alone. The prepared SN covariance still relies on the Notebook 08 first-order Jacobian.
- Complete result and interpretation: `results/11b_fig23_catalog_bootstrap/fig23_catalog_bootstrap_20261014T033440741851Z_c88d56d5/analysis.md`. The separate 100-draw folder was a validation pilot only.
- Verified notebook inspection mode, direct-GLS equivalence tests, and all 16 Figure 2/3 tests. Published-loss discrepancies remain unresolved.
