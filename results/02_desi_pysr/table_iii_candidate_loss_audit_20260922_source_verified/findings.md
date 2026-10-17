# Figure 1 / Table III candidate-loss audit

## Scope

Each printed DESI-only candidate equation was evaluated at the six 
official DESI DR2 $D_H/r_d$ redshifts. The primary result uses the 
exact release central values and covariance and applies the paper's 
Eq. (9), $L=r^T C^{-1}r/n$.

## Result

| complexity | paper_reported_loss | public_loss | table_i_rounded_loss | public_to_paper_ratio | same_at_paper_precision | parse_status |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 2519.06 | 2693.99 | 2802.01 | 1.06944 | False | literal_transcription |
| 2 | 120.73 | 1115.29 | 1118.26 | 9.23788 | False | literal_transcription |
| 3 | 91.5 | 797.367 | 799.597 | 8.71439 | False | literal_transcription |
| 4 | 1.22 | 11.3074 | 11.3502 | 9.26836 | False | literal_transcription |
| 5 | 0.3 | 0.664026 | 0.677843 | 2.21342 | False | literal_transcription |
| 7 | 0.25 | 0.311022 | 0.324566 | 1.24409 | False | literal_transcription |
| 8 | 0.24 | 0.248055 | 0.255927 | 1.03356 | False | literal_transcription |
| 9 | 0.21 | 0.207807 | 0.216853 | 0.989556 | True | literal_transcription |
| 11 | 0.17 | 9.43666 | 9.79827 | 55.5097 | False | literal_transcription |
| 13 | 0.02 | 1.29987 | 1.35894 | 64.9936 | False | literal_transcription |
| 15 | 0.016 | 0.0529975 | 0.0540233 | 3.31234 | False | minimal_repair_of_malformed_source |
| 16 | 0.015 | 0.0631307 | 0.0646034 | 4.20871 | False | literal_transcription |
| 17 | 0.013 | 0.0544846 | 0.0555719 | 4.19113 | False | literal_transcription |

Only 1 of 12 unambiguous printed equations reproduces the reported loss at the number of decimal places shown in Table III.

The closest literal transcriptions by relative difference are:

- complexity 9: paper 0.21, recomputed 0.207807, relative difference 1.04%
- complexity 8: paper 0.24, recomputed 0.248055, relative difference 3.36%
- complexity 1: paper 2519.06, recomputed 2693.99, relative difference 6.94%
- complexity 7: paper 0.25, recomputed 0.311022, relative difference 24.41%

The rounded Table I configuration changes the losses only modestly and does not reconcile the discrepant rows.

## Figure 1 selected expression

The expression printed as Eq. (10), `$28.4*exp(-0.519*z)`, gives $\chi^2=4.222491$ and $L=0.703749$ with the public release data. Using rounded Table I values gives $L=0.718559$. Both differ from the 0.30 reported for the corresponding selected Table III row.

The selected Table III form, `28.42/1.68^z`, gives $L=0.664026$.

## Alternative-loss diagnostic

For complexities 2, 3, and 4, the reported values happen to be close to the unweighted sum of squared residuals: respectively $120.73$ versus $120.984$, $91.5$ versus $90.522$, and $1.22$ versus $1.261$. This pattern does not hold for the other rows, so it does not provide a single alternative loss definition that reproduces Table III.

## Interpretation

The published loss column cannot be reproduced by applying the published Eq. (9) uniformly to the printed equations and the official DESI DR2 data. Coefficient rounding can strongly affect the high-complexity pole-like expressions, but it does not explain the large discrepancies for several low-complexity rows whose constants are printed to many digits.

The official arXiv LaTeX source has double minus signs in complexities 13 and 17. Complexity 15 is malformed (`3.68-z\times` has no following operand); its value uses a documented minimal repair.

Possible explanations include unpublished full-precision equations, a different data/covariance object during training, stale losses copied from another run, or a mismatch between the implemented and stated loss. The available paper information is insufficient to distinguish these possibilities conclusively.
