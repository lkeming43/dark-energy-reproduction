# Figure 1 audit using the arXiv v1 DESI objective

The v1 pointwise definition is `abs(prediction-target)/sigma`. The mean and sum below apply that objective to the v2 DR2 data and Table III expressions. This is not a reproduction of the v1 DR1 analysis.

```text
 complexity  paper_reported_loss  v2_covariance_loss_mean  v1_public_abs_inv_sigma_mean  v1_public_abs_inv_sigma_sum closest_of_tested_definitions
          1             2519.060              2693.989642                     49.278969                   295.673813            v2_covariance_mean
          2              120.730              1115.288699                     20.526956                   123.161735          v1_abs_inv_sigma_sum
          3               91.500               797.366807                     17.617385                   105.704308          v1_abs_inv_sigma_sum
          4                1.220                11.307394                      2.010984                    12.065906         v1_abs_inv_sigma_mean
          5                0.300                 0.664026                      0.716029                     4.296172            v2_covariance_mean
          7                0.250                 0.311022                      0.510818                     3.064907            v2_covariance_mean
          8                0.240                 0.248055                      0.404630                     2.427780            v2_covariance_mean
          9                0.210                 0.207807                      0.374495                     2.246971            v2_covariance_mean
         11                0.170                 9.436657                      2.943319                    17.659912         v1_abs_inv_sigma_mean
         13                0.020                 1.299872                      1.018406                     6.110433         v1_abs_inv_sigma_mean
         15                0.016                 0.052997                      0.171642                     1.029852            v2_covariance_mean
         16                0.015                 0.063131                      0.217570                     1.305417            v2_covariance_mean
         17                0.013                 0.054485                      0.187391                     1.124348            v2_covariance_mean
```

## Conclusion

Across 12 exactly parseable Table III rows, the closest tested definition is the v2 covariance mean for 7, the v1 absolute-error mean for 3, and the v1 absolute-error sum for 2.

A legacy v1 pointwise objective does not provide a single consistent explanation unless it aligns the reported values across the table, rather than improving isolated rows only.
