# Research Learning Notebook

This notebook records concepts learned through research discussions, starting on 2026-09-16. Each entry distinguishes conceptual understanding, locally verified results, and unresolved questions.

## 2026-09-16 — DESI errors, correlations, and covariance

### 1. Are individual standard deviations sufficient?

The standard deviation \(\sigma_i\) describes the error scale of a single measurement. If a covariance matrix is diagonal, the individual standard deviations are sufficient to reconstruct it. If the matrix has nonzero off-diagonal elements, the covariance between measurements is also needed.

Zero covariance does not generally imply statistical independence. Under the multivariate Gaussian observational-error model used here, however, zero cross-covariance corresponds to independent components.

### 2. Is correlation the same as covariance?

They are related but not the same quantity. The correlation coefficient \(\rho_{ij}\) is dimensionless and normally ranges from −1 to 1. Covariance \(C_{ij}\) incorporates both error scales:

\[
C_{ij}=\rho_{ij}\sigma_i\sigma_j,\qquad C_{ii}=\sigma_i^2.
\]

If two observations are both above a model prediction and their errors are strongly positively correlated, they should not be counted as two independent upward fluctuations. The full covariance determines the appropriate weighting of residuals in the fit.

### 3. What does the full DESI dataset contain, and what do we use?

The local DESI mean file contains 13 BAO observables, with a 13 × 13 covariance matrix:

| z | Released observables | Figure 1 | Figures 2/3 |
|---:|---|---|---|
| 0.295 | \(D_V/r_d\) | not used | not used |
| 0.510 | \(D_M/r_d\), \(D_H/r_d\) | \(D_H/r_d\) | \(D_M/r_d\) |
| 0.706 | \(D_M/r_d\), \(D_H/r_d\) | \(D_H/r_d\) | \(D_M/r_d\) |
| 0.934 | \(D_M/r_d\), \(D_H/r_d\) | \(D_H/r_d\) | \(D_M/r_d\) |
| 1.321 | \(D_M/r_d\), \(D_H/r_d\) | \(D_H/r_d\) | \(D_M/r_d\) |
| 1.484 | \(D_M/r_d\), \(D_H/r_d\) | \(D_H/r_d\) | \(D_M/r_d\) |
| 2.330 | \(D_M/r_d\), \(D_H/r_d\) | \(D_H/r_d\) | \(D_M/r_d\) |

- \(D_M/r_d\): transverse-distance observable.
- \(D_H/r_d=c/[H(z)r_d]\): radial-distance observable.
- \(D_V/r_d\): isotropic distance observable combining transverse and radial information.
- The correlations come from the DESI-released covariance. Our symbolic-regression analysis does not re-estimate them.
- Figure 3 keeps the same DESI measurements as Figure 2, replacing Pantheon+ with DESY5.

### 4. What was verified today?

The zero-based source-row indices are:

- DM: [1, 3, 5, 7, 9, 12]
- DH: [2, 4, 6, 8, 10, 11]

The final pair has DH before DM in the source file, so odd/even indexing is not a valid selection rule.

Both 6 × 6 submatrices have exactly zero off-diagonal entries. The six nonzero cross-terms in the full matrix all couple DM and DH at the same redshift, with correlations of roughly −0.35 to −0.49. At \(z=0.934\), the correlation is approximately −0.347.

Therefore the six DESI DM points used in Figures 2/3 can be represented by their six standard deviations under this released covariance convention. The earlier concern that paper standard deviations might be insufficient due to unrecorded DM–DM correlations does not apply to the checked DM subset.

Selecting the DM submatrix does not condition on the DH measurements. DH values are not added as extra information to the DM fit. Zero cross-bin covariance in this release does not prove that every conceivable survey systematic is independent.

### 5. How should public errors and paper-table errors be compared?

Public DM errors are \(\sqrt{\mathrm{diag}(C_{\rm DM})}\). The paper lists six standard deviations, from which a diagonal matrix \(\mathrm{diag}(\sigma_{\rm paper}^2)\) can be constructed.

For LRG3+ELG1 at \(z=0.934\), the public standard deviation is 0.1617815861, while the paper gives 0.152: a difference of about 6.44%. The cause is not yet known, so this alone cannot establish an author error.

Replacing only the BAO variances while retaining the central values and full SN covariance isolates the effect of the error choice. This should be described as a “paper-table-error version,” not the authors’ original covariance.

At this update, the combined NPZ files use public DESI errors. The paper errors are stored in the comparison CSV. A switchable combined dataset using paper-table errors has not yet been implemented.

### 6. Related concepts reviewed today

#### zHD

zHD is the Hubble-diagram redshift after corrections to the CMB frame and for peculiar velocity. zHEL is heliocentric redshift; zCMB is corrected to the CMB frame. At low redshift, peculiar velocity has a larger relative effect, which motivates the paper’s \(z<0.01\) exclusion.

The paper explicitly specifies the threshold but not the column name. Selecting with zHD yields 1590 rows and a redshift range of 0.01016–2.26137, matching the paper and supporting this implementation choice.

#### Fixed \(r_d\) and a shared calibration uncertainty

\(r_d\) is a ruler shared by every SN. The current conversion uses its central value, and the covariance contains only observational error conditional on that fixed calibration. In a future Monte Carlo run, sample one common \(r_d\) per realization and use it for every SN conversion and the corresponding \(w(z)\) reconstruction. Do not independently draw \(r_d\) for each SN. Published BAO dimensionless ratios should not be divided by \(r_d\) again as though they were SN distances.

#### Pantheon+ covariance symmetrization

Theory requires \(C_{ij}=C_{ji}\). The local file has a maximum mismatch of \(3\times10^{-8}\) mag². For example, 0.00048592 and 0.00048595 are both replaced by their average, 0.000485935:

\[
C_{\rm sym}=(C+C^T)/2.
\]

This does not alter central values or diagonal variances; only the derived matrix is changed and the raw file remains intact. An explicit tolerance is enforced, with larger discrepancies raising an error. Its cause is not established and should not automatically be attributed to rounding. This issue is separate from the 6.44% DESI standard-deviation difference.

### 7. Evidence and follow-up questions

- [Executable notebook 08](../notebooks/08_fig23_data_preparation.ipynb)
- [DESI submatrix audit](../results/08_fig23_data_preparation/desi_submatrix_audit.csv)
- [Nonzero DESI correlated pairs](../results/08_fig23_data_preparation/desi_correlated_pairs.csv)
- [Paper versus public error comparison](../results/08_fig23_data_preparation/desi_table1_comparison.csv)
- [Pantheon+ official column and covariance documentation](https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/README)
- [Legacy DESY5 covariance documentation](https://raw.githubusercontent.com/des-science/DES-SN5YR/e3493cb/4_DISTANCES_COVMAT/README.md)

Follow-up work: add a public/paper-table error switch, identify the source of the release-versus-paper difference, and specify a consistent shared-calibration Monte Carlo scheme.



| 数据                        | 文件                                                  | Figure 2 | Figure 3 | 用途                                                         |
| --------------------------- | ----------------------------------------------------- | -------- | -------- | ------------------------------------------------------------ |
| DESI DR2 BAO mean vector    | `Data/DESI_BAO/desi_gaussian_bao_ALL_GCcomb_mean.txt` | 用       | 用       | 取六个横向 BAO 测量 \(D_M/r_d\)，作为距离–红移关系的高红移锚点 |
| DESI DR2 BAO covariance     | `Data/DESI_BAO/desi_gaussian_bao_ALL_GCcomb_cov.txt`  | 用       | 用       | 给出 DESI 六点的误差；在所选 \(D_M/r_d\) 子集内是对角矩阵    |
| Pantheon+ catalog           | `Data/Pantheon/Pantheon+SH0ES.dat`                    | 用       | 不用     | 提供 `zHD` 和校正视星等 `m_b_corr`；筛掉 \(z<0.01\)，保留 1590 条 |
| Pantheon+ covariance        | `Data/Pantheon/Pantheon+SH0ES_STAT+SYS.cov`           | 用       | 不用     | 给 Pantheon+ 星等误差及不同超新星之间的系统相关性            |
| DESY5 catalog               | `Data/DESY5/DES-SN5YR_HD.csv`                         | 不用     | 用       | 提供 `zHD`、距离模数 `MU`、误差 `MUERR_FINAL`；筛掉 `MUERR_FINAL > 1`，保留 1754 条 |
| DESY5 systematic covariance | `Data/DESY5/STAT+SYS.txt.gz`                          | 不用     | 用       | 给 DESY5 的系统误差相关性；需要加上 `diag(MUERR_FINAL²)` 才是总协方差 |

DESI 已经直接给出无量纲的 \(D_M/r_d\)，可直接进入拟合。

Pantheon+ 给的是星等 `m_b_corr`，需要用 \(M=-19.253\) 和 \(r_d=138.28\) Mpc 转成 \(D_M/r_d\)。

DESY5 已给距离模数 `MU`，以 \(H_0=70\) 为标定；再配合 \(r_d=101/0.7\) Mpc 转成 \(D_M/r_d\)。

Figure 2 用 `DESI + Pantheon+`。

Figure 3 用 `DESI + DESY5`；DESI 部分没有增加或改变，只替换超新星样本。
