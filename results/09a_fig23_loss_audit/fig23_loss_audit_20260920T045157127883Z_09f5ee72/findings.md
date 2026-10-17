# Figure 2/3 published-loss audit

Experiment ID: `fig23_loss_audit_20260920T045157127883Z_09f5ee72`.

## Main numerical checks

- Printed Fig. 2 Eq. 12: loss 4.729393; Table IV reports about 1.70.
- Literal Table IV selected row with denominator 0.030: loss 4.788842.
- Illustrative denominator 0.032, which also displays as 0.03 to two decimals: loss 1.686104. This is not a recovered author constant.
- Independent Table IV complexity-10 expression: loss 0.910111; paper reports about 1.68.
- Printed Fig. 3 Eq. 13: loss 1.016472; Table V reports about 2.23.
- Fig. 2 printed Eq. 12: BAO chi-square 526.201; SN chi-square 7021.910.

## Interpretation

The selected Fig. 2 formula is highly sensitive to the denominator rounded to 0.03 in the paper. A value near 0.032 can reproduce its reported loss under our inputs, but this numerical coincidence does not prove the hidden PySR coefficient.
Other printed candidates also have loss discrepancies, so rounding alone has not established full reproduction of the authors' data vector and covariance.
Changing only six DESI variances or dropping SN off-diagonal correlations does not explain the printed Eq. 12 discrepancy.
Next: request or locate the authors' unrounded PySR hall-of-fame output, transformed SN vector, covariance, and exact loss implementation. Until then, keep the data/covariance and formula-precision hypotheses separate.

## Artifacts

`experiment.json`, `input_manifest.csv`, `formula_loss_comparison.csv`, `denominator_scan.csv`, `denominator_anchors.csv`, `fig2_denominator_scan.png`, and `covariance_diagnostics.csv`.

All calculations are read-only with respect to the prepared Notebook 08 data and the Notebook 09 fit results.
