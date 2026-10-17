# Input data

The analysis reads the public survey files below from this directory. The raw files are kept locally and excluded from Git; place downloaded files at the exact paths shown. Do not edit them in place.

| Survey | Local path | Release |
| --- | --- | --- |
| DESI DR2 BAO mean vector | `DESI_BAO/desi_gaussian_bao_ALL_GCcomb_mean.txt` | [Cobaya BAO DR2 files](https://github.com/CobayaSampler/bao_data/tree/master/desi_bao_dr2) |
| DESI DR2 BAO covariance | `DESI_BAO/desi_gaussian_bao_ALL_GCcomb_cov.txt` | Same release |
| Pantheon+ catalog | `Pantheon/Pantheon+SH0ES.dat` | [Pantheon+ distances and covariance](https://github.com/PantheonPlusSH0ES/DataRelease/tree/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR) |
| Pantheon+ STAT+SYS covariance | `Pantheon/Pantheon+SH0ES_STAT+SYS.cov` | Same release |
| DES-SN5YR Hubble diagram | `DESY5/DES-SN5YR_HD.csv` | [DES-SN5YR distance release at e3493cb](https://github.com/des-science/DES-SN5YR/tree/e3493cb/4_DISTANCES_COVMAT) |
| DES-SN5YR STAT+SYS covariance | `DESY5/STAT+SYS.txt.gz` | Same release |

The exact filenames used by the code are recorded in [`configs/fig2_pantheon.yaml`](../configs/fig2_pantheon.yaml) and [`configs/fig3_desy5.yaml`](../configs/fig3_desy5.yaml). The data-preparation notebook checks the expected row counts and covariance dimensions. The DES-SN5YR link is pinned because its main branch now carries a different Hubble diagram and covariance format.

SHA-256 hashes of the local inputs used for this reproduction:

| File | SHA-256 |
| --- | --- |
| `DESI_BAO/desi_gaussian_bao_ALL_GCcomb_mean.txt` | `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585` |
| `DESI_BAO/desi_gaussian_bao_ALL_GCcomb_cov.txt` | `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509` |
| `Pantheon/Pantheon+SH0ES.dat` | `1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8` |
| `Pantheon/Pantheon+SH0ES_STAT+SYS.cov` | `abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc` |
| `DESY5/DES-SN5YR_HD.csv` | `95ab1e9f5460a056de37d2b4e28df73ce0c7fd471baba9fb37f1667f41b7901c` |
| `DESY5/STAT+SYS.txt.gz` | `20bc94a881893235a94fc5a1b1c995f4b35bfae333996c2a0ce2696cb05381f2` |

If you use a different upstream version, record its hash before comparing results.
