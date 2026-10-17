# Notebook 11c — Figure 2 balanced seed audit

Exploratory re-selection **after** Notebook 12 exposed the DESI/SN imbalance. The original Notebook 11 result remains unchanged. No PySR rerun, coefficient refit, or new Monte Carlo is performed.

The fixed rule is recorded in `experiment.json`. Separate subgroup scores use the corresponding marginal covariance blocks; the joint score is their sum divided by 1596.

| Branch | Seed | Complexity | DESI chi2/6 | SN chi2/1590 | Combined chi2/1596 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Notebook 11 primary | 20261064 | 7 | 6.98795 | 0.89718 | 0.92008 |
| Notebook 11c selected | 20261093 | 9 | 1.71371 | 0.88052 | 0.88365 |

New fixed PySR expression: `29.321836856709833*z/1.2843234947453503**z`

The new choice is a seed-selected PySR expression, not the separately postulated exponential comparison. Its appearance in only one of 100 seed-selected rows does not show stable selection.

The second panel uses marginal DESI standard deviations for visual residuals; the subgroup chi-squared uses its full covariance.

Next: if this exploratory choice is retained, compute fresh conditional Monte Carlo bands under the same uncertainty protocol as Notebook 11/11b. Do not relabel the existing comparison-model bands as independent bands.
