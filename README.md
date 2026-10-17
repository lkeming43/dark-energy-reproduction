# Dark Energy: independent symbolic-regression reproduction

An independent reproduction and audit of [*Evidence for dynamical dark energy from DESI-DR2 and SN data? A symbolic regression analysis* (arXiv:2502.10506v2)](https://arxiv.org/abs/2502.10506). The project reconstructs the distance–redshift relations and dark-energy equation of state \(w(z)\) in Figures 1–3 using public DESI DR2 BAO, Pantheon+, and DES-SN5YR data.

**Status (17 October 2026): partial reproduction.** The DESI-only distance trend and an exponential model family are recovered. The published uncertainty bands and several losses calculated from the paper's printed expressions are not reproduced exactly. The Figure 2/3 analyses show that \(w(z)\) intervals depend on the selected symbolic family and how pole-producing Monte Carlo draws are treated. These pointwise, conditional intervals are not a global rejection test of \(\Lambda\)CDM. See the [reproduction report](report/main.tex), [experiment log](notes/reproduction_log.md), and [latest Figure 2/3 diagnostics](results/12_fig23_fit_diagnostics/fig23_fit_check_20261015T025206917931Z_c0798219/analysis.md).

## Repository map

| Path | Contents |
| --- | --- |
| [`notebooks/`](notebooks/) | Numbered exploratory and reproduction notebooks, 01–12 |
| [`src/`](src/) | Reusable data loading, fitting, reconstruction, and diagnostics code |
| [`configs/`](configs/) | Scenario, PySR, and Monte Carlo settings |
| [`tests/`](tests/) | Numerical and provenance checks |
| [`results/`](results/) | Small experiment summaries, tables, logs, and figures retained in Git |
| [`figures/`](figures/) | Selected project figures |
| [`notes/`](notes/) | Methods, literature notes, experiment log, and interim findings |
| [`report/`](report/) | LaTeX report source, bibliography, and report figures |
| [`Data/README.md`](Data/README.md) | Required public input files and their release locations |
| `run_*.py` | Scripts that execute or audit specific notebooks and experiments |

The original paper PDF, downloaded survey files, large prepared arrays, PySR checkpoints, 100-seed run directories, temporary files, and report ZIP exports remain local and are excluded from Git. The largest local run directory is about 17 GB. Its compact summaries and plots are included. No local files are deleted by this Git layout.

## Set up

From the repository root, use Python 3 and a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Install the optional symbolic-search dependency before running PySR notebooks 02, 03, 09, or 10:

```bash
python -m pip install -r requirements-pysr.txt
```

PySR also needs its Julia runtime; follow the [PySR installation instructions](https://ai.damtp.cam.ac.uk/pysr/) for your platform. The repository does not pin a single environment or promise bitwise identical PySR searches across platforms. Runs record seeds and settings in their experiment folders.

Download the six public input files listed in [Data/README.md](Data/README.md) to their specified paths. The Git checkout contains the data guide, not the raw files. The audit script `run_cross_figure_loss_objective_audit.py` also expects a local copy of the [arXiv v2 PDF](https://arxiv.org/pdf/2502.10506v2) named `arXiv-2502.10506v2.pdf` in the repository root.

Check the environment and begin with the data validation notebooks:

```bash
python -m unittest discover -s tests
jupyter lab
```

Open `notebooks/01_desi_data_check.ipynb` for Figure 1 or `notebooks/08_fig23_data_preparation.ipynb` for Figures 2/3. Run notebooks from the repository checkout so their project-root discovery and relative paths resolve correctly. The later notebooks need outputs from earlier stages; rerunning the 100-seed PySR search and full-catalogue Monte Carlo can take substantial time and disk space.

## Analysis sequence

| Stage | Notebooks | Purpose |
| --- | --- | --- |
| DESI only | 01–07 | Validate radial BAO data, search symbolic expressions, examine seed stability, reconstruct \(w(z)\), and propagate uncertainty |
| Combined data preparation and loss audits | 08, 09, 09a, 09b | Convert Pantheon+/DES-SN5YR distances, fit Figure 2/3 expressions, and check printed losses |
| Figure 2/3 robustness | 10, 11, 11b | Search 100 seeds per scenario and error variant; reconstruct \(w(z)\) with fixed-family coefficient draws and separate mock-catalogue refits |
| Follow-up checks | 12, 11c, 11d | Diagnose DESI/SN fit balance, apply an exploratory alternative seed rule, then run new conditional Monte Carlo for that selected exponential |

The labels `11c` and `11d` are follow-ups to the frozen Notebook 11 result; they do not replace it. Notebook 12 is a diagnostic of frozen fits. See [Notebook 11c's analysis](results/11c_fig2_balanced_seed_audit/fig2_balanced_20261016T044854959271Z_1a017e27/analysis.md) and [Notebook 11d's analysis](results/11d_fig2_selected_exponential_mc/fig2_exp_mc_20261016T045826099551Z_1de02bcb/analysis.md) for their scope and limitations.

The Python runners at the repository root execute selected notebook workflows. For example, `python run_notebook12.py` performs a fresh diagnostic run **after** its prerequisite prepared data and frozen-fit outputs exist. Some runners create a new timestamped results directory; consult the notebook and its configuration before rerunning a large experiment.

## Interpretation and reproducibility

- Figures 1–3 use covariance-weighted losses on the public data. The project records the adopted data selections, calibrations, and covariance treatment in [`configs/`](configs/) and [the methods note](notes/Loss_Recalculation_Workflow_Figures_1_to_3.md).
- The selected expression can change across PySR seeds. A good combined loss may hide a worse fit to the six DESI points because the supernova sample contains thousands of observations; [Notebook 12](results/12_fig23_fit_diagnostics/fig23_fit_check_20261015T025206917931Z_c0798219/analysis.md) reports the subset scores separately.
- Monte Carlo bands are pointwise and conditional on fixed model families. Pole-free bands additionally condition on a numerical exclusion rule. Neither includes a full PySR search for every mock catalogue.
- Raw survey files are never edited in place. Regenerated arrays and checkpoints are local; commit-ready summaries, tables, logs, and figures are in `results/`. To reproduce late-stage notebooks from a fresh clone, regenerate the ignored intermediate arrays and search runs first.
- Published-loss discrepancies remain unresolved. The numerical audits and qualifications are recorded in [`notes/reproduction_log.md`](notes/reproduction_log.md) and the [report](report/main.tex).

The original paper and survey releases retain their own attribution and terms. This repository contains an independent analysis, not the authors' implementation.
