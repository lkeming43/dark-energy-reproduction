# Loss Recalculation Workflow for Figures 1-3

**Paper:** *Evidence for dynamical dark energy from DESI-DR2 and SN data? A symbolic regression analysis*  
**Prepared by:** Liu Keming  
**Date:** 22 September 2026

## 1. Purpose

Before rerunning PySR for Figures 2 and 3, I checked whether the losses reported in Tables III-V can be recovered by evaluating the published candidate expressions directly on the public data.

## 2. Loss definition

I used the covariance-weighted loss stated in the paper:

$$
\mathcal{L}=\frac{1}{N}(\mathbf{y}-\mathbf{f})^{T}
\mathbf{C}^{-1}(\mathbf{y}-\mathbf{f}).
$$

where:

- $N$ is the total number of fitted data points;
- $\mathbf{y}$ is the prepared data vector;
- $\mathbf{f}$ is the published expression evaluated at the same redshifts;
- $\mathbf{C}$ is the covariance matrix in the same observable space as $\mathbf{y}$.

For numerical stability, I calculated the quadratic form by solving the linear system

$$
\mathbf{C}\mathbf{x}=\mathbf{y}-\mathbf{f}
$$

and then evaluating

$$
\chi^2=(\mathbf{y}-\mathbf{f})^T\mathbf{x},
\qquad
\mathcal{L}=\frac{\chi^2}{N}.
$$

I did not explicitly form $\mathbf{C}^{-1}$.

The sign convention for the residual does not affect the result because the loss is quadratic.

## 3. Data vectors

### Figure 1: DESI only

The Figure 1 data vector contains the six public DESI DR2 radial BAO measurements of

$$
\frac{D_H}{r_d}
$$

and the corresponding $6\times6$ submatrix of the public DESI covariance matrix.

Therefore:

$$
N=6.
$$

No supernova conversion or Jacobian propagation is involved in Figure 1.

### Figure 2: DESI plus Pantheon+

The Figure 2 vector contains:

- six public DESI DR2 measurements of $D_M/r_d$;
- 1,590 Pantheon+ supernovae remaining after applying $z_{\mathrm{HD}}>0.01$.

Therefore:

$$
N=6+1590=1596.
$$

The selected Pantheon+ sample spans

$$
0.01016\le z_{\mathrm{HD}}\le2.26137.
$$

The DESI points are placed first, followed by the selected supernovae in their original row order.

### Figure 3: DESI plus DES-SN5YR

The Figure 3 vector contains:

- six public DESI DR2 measurements of $D_M/r_d$;
- 1,754 DES-SN5YR supernovae remaining after applying $\mathrm{MUERR\_FINAL}<1$.

Therefore:

$$
N=6+1754=1760.
$$

The selected DES-SN5YR sample spans

$$
0.02509\le z_{\mathrm{HD}}\le1.12132.
$$

## 4. Conversion of supernova data into $D_M/r_d$

The published symbolic expressions for Figures 2 and 3 describe $D_M/r_d$, whereas the supernova releases provide magnitudes or distance moduli. I therefore converted the supernova observables into the same distance space before combining them with DESI.

### Pantheon+

Pantheon+ provides the corrected apparent magnitude $m_{b,\mathrm{corr}}$. I used

$$
\mu=m_{b,\mathrm{corr}}-M
$$

with

$$
M=-19.253,
\qquad
r_d=138.28\ \mathrm{Mpc}.
$$

### DES-SN5YR

DES-SN5YR directly provides the distance modulus $\mu$. I used

$$
H_0=70\ \mathrm{km\,s^{-1}\,Mpc^{-1}},
\qquad
r_d=144.285714\ \mathrm{Mpc}.
$$

### Common distance conversion

For both supernova samples, I calculated

$$
D_L=10^{(\mu-25)/5}\ \mathrm{Mpc}
$$

and then

$$
\frac{D_M}{r_d}=\frac{D_L}{(1+z_{\mathrm{HD}})r_d}.
$$

The redshift used in the fit and in the $(1+z)$ factor is $z_{\mathrm{HD}}$.

For this deterministic loss audit, $r_d$ is fixed at its central value. Its uncertainty is a common calibration uncertainty and is reserved for later Monte Carlo propagation rather than being added independently to every supernova point.

## 5. Covariance propagation

The supernova covariance is originally defined in magnitude or distance-modulus space. I propagated it to $D_M/r_d$ space using the first-order Jacobian.

For each supernova,

$$
y_i=\frac{D_M(z_i)}{r_d}
$$

and

$$
J_i=\frac{\partial y_i}{\partial\mu_i}
=\frac{\ln 10}{5}y_i.
$$

The propagated covariance is therefore

$$
\mathbf{C}_{\mathrm{distance}}
=\mathbf{J}\mathbf{C}_{\mathrm{magnitude}}\mathbf{J}^{T},
$$

where $\mathbf{J}$ is diagonal.

For Pantheon+, I used the public STAT+SYS covariance. Its very small numerical asymmetry, with a maximum difference of approximately $3\times10^{-8}$ between transposed elements, was removed by averaging the two matrix triangles:

$$
\mathbf{C}\leftarrow\frac{\mathbf{C}+\mathbf{C}^{T}}{2}.
$$

No eigenvalue regularisation was applied.

For DES-SN5YR, the released systematic covariance was combined with the statistical variances from $\mathrm{MUERR\_FINAL}$ before the same Jacobian propagation.

The final Figure 2 and Figure 3 covariance matrices are block diagonal:

$$
\mathbf{C}_{\mathrm{total}}=
\begin{pmatrix}
\mathbf{C}_{\mathrm{DESI}} & \mathbf{0}\\
\mathbf{0} & \mathbf{C}_{\mathrm{SN}}
\end{pmatrix}.
$$

This assumes no cross-covariance between DESI and the supernova survey.

Both prepared matrices are finite, symmetric, and positive definite.

## 6. Candidate-expression evaluation

For every sufficiently well-defined expression in Tables III-V, I followed the same procedure:

1. Parse the printed expression literally.
2. Evaluate it at every redshift in the corresponding prepared vector.
3. Check that every predicted value is finite and real.
4. Form the residual vector $\mathbf{r}=\mathbf{y}-\mathbf{f}(z)$.
5. Solve $\mathbf{C}\mathbf{x}=\mathbf{r}$.
6. Calculate $\chi^2=\mathbf{r}^{T}\mathbf{x}$.
7. Divide by the total number of points to obtain $\mathcal{L}=\chi^2/N$.
8. Compare this value with the loss printed in the paper.

Expressions that are syntactically incomplete in the official source were not used to draw strict conclusions. Any minimally repaired interpretation was labelled explicitly.

## 7. Alternative loss definitions tested

To check whether the discrepancy could arise from a different implementation of the objective, I also calculated:

1. Full-covariance $\chi^2/N$;
2. diagonal-covariance $\chi^2/N$;
3. unweighted sum of squared errors, $\mathrm{SSE}=\sum_i r_i^2$;
4. unweighted mean squared error, $\mathrm{MSE}=N^{-1}\sum_i r_i^2$;
5. full-covariance $\chi^2$ without division by $N$.

No single tested definition reproduces an entire published loss table.

For Figure 1, some low-complexity candidates are close to unweighted SSE, while other candidates are close to the stated covariance-weighted loss. This pattern is not consistent across Table III.

## 8. Printed-constant rounding test

I also tested whether the loss discrepancies could be explained by rounding of the constants displayed in the paper.

Each printed constant was allowed to vary only within the interval that would round back to its displayed value. For example:

$$
33.3\longrightarrow[33.25,33.35),
\qquad
0.274\longrightarrow[0.2735,0.2745).
$$

The search included the printed centre, interval boundaries, deterministic interior samples, and coordinate sweeps. This was a constrained rounding test, not an unrestricted refit.

For the selected main-text Figure 2 expression,

$$
f(z)=33.3z\exp(-0.274z),
$$

the allowed rounding box produces losses of approximately

$$
4.62\lesssim\mathcal{L}\lesssim4.84,
$$

which does not reach the paper value of $1.70$.

## 9. Representative results

| Figure | Paper location | Published expression | Reported selected loss | Recomputed full-covariance loss |
|---|---|---|---:|---:|
| Figure 1 | Table III selected row | $28.42/1.68^z$ | 0.30 | 0.664 |
| Figure 1 | Main-text Eq. (10) | $28.4\exp(-0.519z)$ | 0.30* | 0.704 |
| Figure 2 | Table IV selected row, literal source | $0.76^2z/0.03$ | 1.70 | 21.333 |
| Figure 2 | Appendix Eq. (A2) | $0.76^z\,z/0.03$ | 1.70* | 4.789 |
| Figure 2 | Main-text Eq. (12) | $33.3z\exp(-0.274z)$ | 1.70* | 4.729 |
| Figure 3 | Table V selected row | $30.1z/1.31^z$ | 2.23 | 0.990 |
| Figure 3 | Appendix Eq. (A3) | $30.1z/1.30^z$ | 2.23* | 1.015 |
| Figure 3 | Main-text Eq. (13) | $30.1z\exp(-0.262z)$ | 2.23* | 1.016 |

\* The loss is printed on the selected table row. The same value is shown here as the reference for the appendix and main-text representations of that selected expression; it is not reported as a separate new fit.

The discrepancy is not uniform across the candidate expressions. Some expressions agree reasonably well with the reported losses, while others differ by factors of approximately 9-65. A few source expressions produce even larger discrepancies because their literal printed forms appear inconsistent with the corresponding appendix or main-text expressions.

### Figure-level summary

| Figure | Conclusive printed candidates | Within 5% under the stated loss | Reached by the constrained rounding scan |
|---|---:|---:|---:|
| Figure 1 | 12 | 2 | 3 |
| Figure 2 | 12 | 0 | 0 |
| Figure 3 | 13 | 0 | 0 |

For Figure 2, the official Table IV source and Appendix A2 contain different versions of the selected expression:

$$
\begin{aligned}
\text{Table IV source:}\quad & \frac{0.76^2z}{0.03},\\
\text{Appendix A2:}\quad & \frac{0.76^z\,z}{0.03},\\
\text{Main-text Eq. (12):}\quad & 33.3z\exp(-0.274z).
\end{aligned}
$$

The appendix expression, with exponent $z$, is approximately equivalent to the main-text exponential expression. The exponent-$2$ form is not equivalent and was evaluated separately rather than silently corrected.

## 10. Sensitivity checks

I repeated the relevant calculations using:

- the public DESI covariance and the rounded errors printed in the paper;
- the full covariance and diagonal-only covariance;
- Pantheon+ STAT+SYS and STATONLY covariance;
- alternative Pantheon+ redshift columns;
- alternative reasonable values of $r_d$;
- magnitude-space and transformed-distance-space residual diagnostics;
- the alternative objective definitions listed above.

These choices change individual values, but none consistently reconciles Tables III-V.

For the selected Figure 2 expression, approximately $93\%$ of the total $\chi^2$ in the primary calculation comes from the supernova block. However, changing the tested supernova covariance and redshift conventions does not reduce the total loss from approximately $4.73$ to $1.70$ while also explaining the other candidate rows.

Because similar discrepancies occur in Figure 1, which contains no supernova conversion, the Jacobian transformation cannot be the common explanation across all three figures.

## 11. Current interpretation

I found that rounding of the coefficients printed in the paper may explain some of the loss discrepancies. This is particularly relevant for expressions whose losses are highly sensitive to small changes in their coefficients. However, varying every printed constant within its allowed rounding interval does not explain all of the differences across Tables III-V. Therefore, coefficient rounding may be a contributing factor, but it is not a complete explanation for the observed discrepancies.

## 12. Available supporting material

The following can be provided if useful:

- candidate-by-candidate loss comparisons for Tables III-V;
- full-covariance, diagonal-covariance, SSE, MSE, and unnormalised chi-squared values;
- constrained rounding intervals and scan results;
- data-preparation notebooks and validation summaries;
- the complete loss-audit notebook and executable scripts.
