# Notebook 11b: full-catalogue Monte Carlo refits

Completed 2026-10-14. This experiment is separate from the unchanged Notebook 11 run. It fixes the same three model-family branches chosen there, draws correlated mock catalogues from the complete prepared covariance, and refits the family to every mock. No PySR structure search is repeated. The experiment used 20,000 draws for each of two centring choices in each branch: 120,000 refits total.

The branches are Figure 2 independent quadratic, Figure 2 refitted exponential comparator, and Figure 3 independent exponential. The deterministic losses exactly reproduce the completed Notebook 11 fits to the reported precision: 0.920080, 0.883650, and 0.937176, respectively. The data and frozen Notebook 11 manifest hashes are recorded in `experiment.json` and each `branch.json`.

Each draw uses one shared value of r_d. The supernova mean vector and covariance are coherently rescaled as 1/r_d and 1/r_d² in prepared distance space, respectively. DESI measurements remain fixed in D_M/r_d space. The initial supernova covariance transformation from native observables to D_M/r_d was already first-order (Notebook 08); this run does not regenerate native magnitudes. A whitened sufficient-statistic implementation is algebraically equivalent to materializing each correlated mock catalogue for GLS; the exponential basis approximation was checked against direct evaluation to better than 2e-10 in whitened coordinates. Unit tests independently compare both families against direct GLS.

## Main result

The new catalogue-refit medians and 95% pointwise bands are very close to the earlier first-order coefficient-Monte-Carlo bands at representative redshifts. For example, Figure 2 quadratic at z=0.497 has a new w median -1.348 versus -1.349 previously; Figure 3 exponential at z=0.497 has -1.067 versus -1.065. At high redshift, the quadratic w bands can still be extremely wide because of denominator poles. This experiment therefore does **not** remove the model-family and pole-selection sensitivity already seen in Notebook 11.

With fitted-expression centring, the pole-free counts out of 20,000 are 19,486 (Figure 2 quadratic), 18,533 (Figure 2 exponential), and 19,640 (Figure 3 exponential). The corresponding observed-data-centred counts are 19,534, 18,579, and 19,620. These are diagnostic counts, not posterior model probabilities.

The `pole_free` selection does not additionally discard distance curves that become nonincreasing on the grid; that condition is recorded separately in `fit_samples.csv.gz`. It occurs in 680 and 653 of the Figure 2 quadratic realizations under the two centres, and in none of the exponential branches. Any stronger physical admissibility cut should therefore be specified and reported as a *different* sensitivity analysis, not silently substituted into these bands.

At z=0, the pole-free 3-sigma pointwise w intervals are approximately [-1.066, -0.982] for Figure 2 quadratic, [-0.983, -0.848] for Figure 2 exponential, and [-0.981, -0.845] for Figure 3 exponential (fitted-expression centring). Thus w=-1 is inside the quadratic band but outside the two exponential bands at this endpoint. This is a model-selection and pole-conditioning sensitivity, not a global statistical rejection of LambdaCDM. Supernovae include lower-redshift data than DESI BAO, but w itself is inferred from derivatives of the fitted distance curve.

The comparison to Notebook 11 is not a perfectly isolated test of linearization: the new run coherently rescales the SN covariance with r_d whereas Notebook 11 uses a first-order parameter response at nominal covariance. The identical 2.5-97.5% quantiles are used in `comparison_to_linearized.csv` so the intervals themselves are comparable.

## Files

- `experiment.json`: run settings, input provenance, completion status, and branch summary.
- `overview.csv`: six branch/centering diagnostic rows.
- `comparison_to_linearized.csv`: selected-redshift comparison with Notebook 11.
- `<scenario>/<branch>/branch.json`: fixed expression and deterministic fit.
- `<scenario>/<branch>/<center>/batches/*.npz`: resumable MC batches.
- `<scenario>/<branch>/<center>/fit_samples.csv.gz`: fitted coefficients, chi-squared, shared r_d, and quality flags for every draw.
- `<scenario>/<branch>/<center>/bands.csv`: pointwise distance and w quantiles for all and pole-free draws.
- `<scenario>/<branch>/<center>/reconstruction.png`: distance and w plots, including paper and LambdaCDM references.

The 100-draw folder next to this one is a validation pilot, not a scientific result. Open `notebooks/11b_fig23_catalog_bootstrap.ipynb` in inspection mode to browse this completed run, or use `scripts/notebooks/run_notebook11b.py` to create a new timestamped experiment. Set `--resume` to an existing folder to continue interrupted batches without redrawing them.
