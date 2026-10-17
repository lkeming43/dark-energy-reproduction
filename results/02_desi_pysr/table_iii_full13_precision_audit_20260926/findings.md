# Figure 1 full-DESI precision-block audit

## Question

Can the Table III loss discrepancies be explained by first inverting the complete 13x13 public DESI covariance and then extracting the D_H precision block?

## Result

Among the 12 exactly parseable Table III expressions, the full-precision-block loss with N=6 is closer to the paper for 0. It is farther from the paper for every exactly parseable expression.

The full-precision-block calculation increases the loss because it conditions on the unused correlated D_M measurements being exactly at their model values. It does not represent the standard likelihood for a D_H-only subset.

Dividing the same padded-vector chi-square by 13 can make some numerical values look closer, but that changes N from the six fitted measurements to thirteen while assigning zero residual to seven observables that were not modelled. It is therefore not a coherent implementation of the paper's stated D_H-only fit.

## Candidate-by-candidate values

| complexity | paper_reported_loss | standard_submatrix_loss_n6 | full13_precision_block_loss_n6 | precision_vs_standard_ratio | full13_n6_to_paper_ratio | full13_n6_closer_to_paper | full13_precision_zero_padded_loss_n13 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 2519.06 | 2693.99 | 3206.08 | 1.19008 | 1.27273 | False | 1479.73 |
| 2 | 120.73 | 1115.29 | 1369.1 | 1.22758 | 11.3402 | False | 631.893 |
| 3 | 91.5 | 797.367 | 978.849 | 1.2276 | 10.6978 | False | 451.777 |
| 4 | 1.22 | 11.3074 | 13.8871 | 1.22815 | 11.3829 | False | 6.40945 |
| 5 | 0.3 | 0.664026 | 0.808646 | 1.21779 | 2.69549 | False | 0.373221 |
| 7 | 0.25 | 0.311022 | 0.370594 | 1.19154 | 1.48238 | False | 0.171043 |
| 8 | 0.24 | 0.248055 | 0.297704 | 1.20016 | 1.24043 | False | 0.137402 |
| 9 | 0.21 | 0.207807 | 0.246074 | 1.18415 | 1.17178 | False | 0.113573 |
| 11 | 0.17 | 9.43666 | 11.2364 | 1.19071 | 66.0963 | False | 5.18601 |
| 13 | 0.02 | 1.29987 | 1.55213 | 1.19407 | 77.6067 | False | 0.71637 |
| 15 | 0.016 | 0.0529975 | 0.0673983 | 1.27173 | 4.2124 | False | 0.0311069 |
| 16 | 0.015 | 0.0631307 | 0.0787332 | 1.24715 | 5.24888 | False | 0.0363384 |
| 17 | 0.013 | 0.0544846 | 0.0690594 | 1.2675 | 5.31226 | False | 0.0318736 |

## Main-text Eq. (10) cross-check

| paper_reported_loss | standard_submatrix_loss_n6 | full13_precision_block_loss_n6 | full13_precision_zero_padded_loss_n13 |
| --- | --- | --- | --- |
| 0.3 | 0.703749 | 0.856272 | 0.395202 |

## Interpretation

This test rejects the full-13x13-precision-block hypothesis as a common explanation for Table III. The N=6 version moves every exactly parseable candidate away from its reported loss, and the N=13 version relies on an inconsistent normalization and still does not align the candidates uniformly.

Complexity 15 remains non-conclusive because the official source expression is malformed; the evaluated expression is only a documented minimal repair.
