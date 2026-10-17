# Cross-figure loss-objective and rounding audit

Experiment: `cross_figure_loss_audit_20260922T030440370758Z_e42b0972`.

## Figure-level summary

| figure | candidates | ratio_min | ratio_max | exact_at_5pct | rounding_boxes_bracketing_paper_loss |
| --- | --- | --- | --- | --- | --- |
| fig1 | 13 | 0.989556 | 134.815 | 2 | 3 |
| fig2 | 12 | 0.541733 | 54907.4 | 0 | 1 |
| fig3 | 13 | 0.215739 | 36012.3 | 0 | 0 |

## Closest tested objective by candidate count

| figure | closest_objective | count |
| --- | --- | --- |
| fig1 | diagonal_cov_chi2_over_n | 1 |
| fig1 | full_cov_chi2_over_n | 3 |
| fig1 | unweighted_mse | 5 |
| fig1 | unweighted_sse | 4 |
| fig2 | diagonal_cov_chi2_over_n | 3 |
| fig2 | full_cov_chi2_over_n | 7 |
| fig2 | unweighted_mse | 2 |
| fig3 | diagonal_cov_chi2_over_n | 2 |
| fig3 | full_cov_chi2_over_n | 3 |
| fig3 | unweighted_mse | 8 |

## Findings

- No single multiplicative normalization explains any of the three printed loss tables: the full-covariance loss ratio changes direction within Figures 2 and 3 and spans orders of magnitude in Figure 1.
- Figure 1 complexities 2--4 remain close to unweighted SSE, whereas complexities 8--9 are close to the stated full-covariance loss. No one tested objective reproduces Table III uniformly.
- The corrected rendered-PDF transcriptions of Figure 1 complexities 13, 15 and 17 use a single minus sign before the rational term.
- A rounding box marked `paper_loss_bracketed=True` means the bounded numerical search found losses on both sides of the paper value while respecting every displayed decimal. A False result is failure to reach the value in this search, not a mathematical impossibility proof.
- Remaining explanations that require author material are full-precision hall-of-fame expressions, a different fit-ready vector/covariance, a different loss implementation, or formula/loss rows taken from different runs.

## Artifacts

`objective_comparison.csv`, `rounding_box_scan.csv`, `rounding_parameter_intervals.csv`, `closest_objective_counts.csv`, `summary_by_figure.csv`, optional `loss_ratio_by_complexity.png`, and `experiment.json`.
