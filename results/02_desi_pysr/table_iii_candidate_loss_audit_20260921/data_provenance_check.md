# DESI DR2 input-data provenance check

## Question-by-question verdict

| Question | Verdict | Evidence |
| --- | --- | --- |
| Did the authors use pre-release DESI DR2 files? | **No supporting evidence; unlikely for the published DR2 revision.** | The symbolic-regression paper v1 (14 February 2025) used DESI DR1. Its DR2 revision, v2, was submitted on 13 June 2025. DESI DR2 results were released on 19 March 2025, and the public Gaussian BAO files were added to `CobayaSampler/bao_data` on 20 March 2025. |
| Is the current public mean exactly the same file used by the authors? | **The displayed central values agree, but bit-for-bit identity cannot be proved.** | Every DESI-only Table I central value is the current public mean rounded to three decimal places. The authors did not publish their exact input file, hash, or analysis code. |
| Was the public covariance updated? | **Not in the public `CobayaSampler/bao_data` history.** | Both the mean and covariance files have only one commit (`b7b8a36e`, 20 March 2025), and neither differs from that initial commit at the current repository head. This does not rule out a private pre-release file or a different file source. |
| Are the Table I numbers display values? | **Yes.** | The central values are three-decimal renderings of the more precise public mean. Table I reproduces the DESI DR2 Results II v2 Table 4, whose caption defines the quoted values as marginalized-posterior means and standard deviations. The printed table is not a complete likelihood or covariance specification. |
| Did the authors use more precision than Table I? | **Plausible, but not documented and therefore not confirmable.** | A full-covariance loss normally requires machine-readable means and covariance, but the paper supplies neither the exact file nor a hash. Coefficient/data precision cannot explain the large Table III loss discrepancies found in the candidate-loss audit. |

## Timeline

- 14 February 2025: symbolic-regression paper v1 submitted; this version used DESI DR1.
- 19 March 2025: DESI DR2 papers released.
- 20 March 2025: the public Gaussian DESI DR2 mean and covariance files were added to `CobayaSampler/bao_data`.
- 26 March 2025: DESI DR2 Results II v2 submitted. Its Table 4 contains exactly the values reproduced in the symbolic-regression paper's Table I.
- 13 June 2025: symbolic-regression paper v2 submitted with the DESI DR2 analysis.

## File identity and history

Public-repository commit:

`b7b8a36e9bccb063081f811f323cada21ab5fbdd` (`add desi_bao_dr2`, 2025-03-20 22:10:01 UTC)

Local/public SHA-256 values:

- mean: `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585`
- covariance: `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509`

The local files are byte-for-byte identical to the current public repository files. The repository contains no later commit changing either file.

## Mean-value comparison

| tracer | current public mean | paper Table I | relation |
| --- | ---: | ---: | --- |
| LRG1 | 21.86294686 | 21.863 | rounded to 3 decimals |
| LRG2 | 19.45534918 | 19.455 | rounded to 3 decimals |
| LRG3+ELG1 | 17.64149464 | 17.641 | rounded to 3 decimals |
| ELG2 | 14.17602155 | 14.176 | rounded to 3 decimals |
| QSO | 12.81699964 | 12.817 | rounded to 3 decimals |
| Ly-alpha | 8.631545674846294 | 8.632 | rounded to 3 decimals |

## Why the Table I errors differ from the public covariance diagonal

The Table I uncertainties are copied from DESI DR2 Results II v2 Table 4 and are reported as marginalized-posterior standard deviations. The released covariance belongs to a Gaussian likelihood data product. These are not necessarily the same summary statistic; parameter-basis conversion and Gaussian approximation are possible reasons for a small difference, but the public material inspected here does not identify the exact cause.

For the six radial measurements, the public covariance gives approximately

`[0.428868, 0.333870, 0.201043, 0.224551, 0.518012, 0.101062]`,

whereas Table I prints

`[0.425, 0.330, 0.193, 0.221, 0.516, 0.101]`.

Therefore the difference does **not**, by itself, show that DESI later changed the public covariance.

## Consequence for the Table III loss discrepancy

The evidence does not support a later public covariance update as the main explanation. It also does not support ordinary Table I rounding as the explanation: using exact public inputs versus rounded Table I inputs changes the audited losses only modestly, while several Table III discrepancies are factors of about 4--65. The remaining plausible causes are an unpublished input or transformation, unpublished full-precision equations, stale table values from another run, or a mismatch between the implemented and documented loss.

## Primary sources

- Symbolic-regression paper v1 (DESI DR1): <https://arxiv.org/html/2502.10506v1>
- Symbolic-regression paper v2 (DESI DR2): <https://arxiv.org/html/2502.10506v2>
- DESI DR2 Results II v2, including Table 4: <https://arxiv.org/html/2503.14738v2>
- DESI DR2 publication page: <https://data.desi.lbl.gov/doc/papers/dr2/>
- Public-file addition commit: <https://github.com/CobayaSampler/bao_data/commit/b7b8a36e9bccb063081f811f323cada21ab5fbdd>
