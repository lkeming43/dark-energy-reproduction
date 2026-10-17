# Notebook 11d — MC for the Figure 2 seed-20261093 exponential

The raw PySR expression is the source of the fixed model structure and the model-centered mock-catalogue mean. Each Gaussian mock catalogue is refitted within the exponential family using the full prepared covariance.

Raw expression: `29.321836856709833*z/1.2843234947453503**z`

Equivalent form: `D_M/r_d = 29.321836856710 z exp(-0.250232116484 z)`

Raw-expression combined loss: 0.883650; GLS refitted-family loss: 0.883650.

| Mock centre | Draws | Finite | Pole crossing | Pole-free |
| --- | ---: | ---: | ---: | ---: |
| raw_selected_expression | 20000 | 20000 | 1478 | 18522 |
| observed_data | 20000 | 20000 | 1466 | 18534 |

The raw central equation gives w(0) = -0.942 and w(2.33) = 2.142. For the raw-expression-centred MC, the **unfiltered** pointwise 3-sigma band includes w = -1 at all listed checkpoints (z = 0, 0.101, 0.202, 0.497, 1.002, 1.499, 2.330); it misses at only one isolated grid point near z = 1.577. That one-point feature should not be interpreted as a robust exclusion.

After **pole filtering**, w = -1 falls outside the pointwise 3-sigma band near z = 0–0.039 and z = 1.577–2.33. At z = 2.33, the conditional 3-sigma interval is approximately [2.07, 3.25]. The observed-data-centred MC has a similar pole-free pattern. These are descriptive, grid-dependent, conditional intervals—not a global test of LambdaCDM.

The two centres answer different questions: simulated data around the raw selected expression versus around the observed catalogue. Neither re-selects a PySR structure for each draw.

All bands are pointwise. The pole-free bands are conditional on the filter and cannot replace the unfiltered result. Check `w_checkpoints.csv` for selected redshifts; `all_bands.csv` contains the full grid.

At z below the lowest DESI point (0.51), DESI-only support is an extrapolation; supernova data do extend to lower redshift.

This is an exploratory follow-up to Notebook 11c, not evidence that the printed paper expression or reported loss has been reproduced.
