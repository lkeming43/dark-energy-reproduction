# Notebook 11 — Figure 2/3 reconstruction (10 October 2026)

## Scope

The same public-covariance selection and reconstruction pipeline was applied to
DESI DR2 + Pantheon+ (Figure 2) and DESI DR2 + DESY5 (Figure 3). The completed
Notebook 10 100-seed experiment and immutable Notebook 08 prepared data were
used; their hashes are recorded in `experiment.json` and each `scenario.json`.
The Table I error variant is not used for the primary result.

## Predefined expression choice

On the fixed 0–2.33 grid, candidate expressions had to satisfy finite
distance and derivatives, zero distance at zero redshift, positive/increasing
distance, finite reconstructed w, and a positive Eq. (5) denominator factor.
The candidate with the lowest complexity among those within 1.5 times the
smallest eligible original loss was selected; loss and then seed break ties.
This is **our independent rule**, not an assertion about the paper's rule.

| Scenario | Eligible seed selections | Chosen seed | Independent structure | Original loss | Exponential GLS loss |
| --- | ---: | ---: | --- | ---: | ---: |
| Figure 2 | 29/100 | 20261064 | `z*(28.99729801609955 - 5.71445618323668*z)` | 0.920080 | 0.883650 |
| Figure 3 | 29/100 | 20261023 | `z*29.43862992267616/1.2951851937274546**z` | 0.937176 | 0.937176 |

The separately refitted exponential coefficients are approximately
`29.321837*z*exp(-0.250232*z)` for Figure 2 and
`29.438630*z*exp(-0.258654*z)` for Figure 3. The paper's printed expressions
are retained as fixed visual comparisons, not substituted for the independent
selected fit.

## Monte Carlo method and diagnostics

Each branch used 10,000 draws from a first-order GLS coefficient-covariance
approximation, plus one common sound-horizon draw per realization. A coherent
first-order response of all supernova distances to this shared calibration
draw was included. PySR was **not** rerun for each draw; the coefficient bands
are conditional on the stated structure and fixed Omega_m=0.3.

| Scenario | Branch | Finite draws | Pole-free draws |
| --- | --- | ---: | ---: |
| Figure 2 | Independent quadratic | 10,000 | 9,700 |
| Figure 2 | Refitted exponential | 10,000 | 9,257 |
| Figure 3 | Independent exponential | 10,000 | 9,772 |
| Figure 3 | Refitted exponential | 10,000 | 9,778 |

The unfiltered w(z) quantiles become very wide at higher redshift because some
draws approach the reconstruction denominator zero. The pole-free bands are
a conditional diagnostic, not a replacement for the unfiltered result. The
gray seed-to-seed spread in the figures is a separate structural-sensitivity
display and is **not** an observational confidence interval. All 1/2/3-sigma
bands are pointwise; none is a global test that rejects LambdaCDM.

## Files

- `overview.csv`: selected structures and summary losses.
- `<scenario>/candidate_diagnostics.csv`: all 100 public-covariance candidates and selection gates.
- `<scenario>/scenario.json`: source hashes, selected expression, fitted coefficients, covariance, and diagnostics.
- `<scenario>/*_bands.csv`: all-draw and pole-free pointwise distance/w quantiles.
- `<scenario>/*_mc_samples.npz`: coefficient, shared-rd, and pole flags for the draws.
- `<scenario>/eligible_seed_w_spread.csv`: between-seed central-curve quantiles.
- `<scenario>/reference_curves.csv`: raw selected, printed paper, and LambdaCDM curves.
- `<scenario>/figure_reconstruction.png`: two-panel figure.

The published candidate-loss discrepancy investigated in Notebook 09a remains
unresolved. Notebook 11 is an explicitly independent reconstruction, not a
claim of exact reproduction of the authors' confidence intervals.
