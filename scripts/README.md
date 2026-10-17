# Run scripts

Run these commands from the repository root. They use the checkout's local `Data/` and `results/` directories; see the main [README](../README.md) for setup and input data.

## Notebook workflows

| Script | Work performed |
| --- | --- |
| [`run_notebook09.py`](notebooks/run_notebook09.py) | Execute the Figure 2/3 PySR fit notebook; computationally expensive |
| [`run_notebook09a.py`](notebooks/run_notebook09a.py) | Execute the Figure 2/3 printed-loss audit notebook |
| [`run_notebook10.py`](notebooks/run_notebook10.py) | Run or resume the Figure 2/3 multi-seed PySR search; computationally expensive |
| [`run_notebook11.py`](notebooks/run_notebook11.py) | Reconstruct Figure 2/3 curves with fixed-family coefficient Monte Carlo |
| [`run_notebook11b.py`](notebooks/run_notebook11b.py) | Run or resume full-catalogue fixed-family Monte Carlo |
| [`run_notebook11c.py`](notebooks/run_notebook11c.py) | Audit an exploratory Figure 2 seed-selection rule |
| [`run_notebook11d.py`](notebooks/run_notebook11d.py) | Run conditional Monte Carlo for the selected Figure 2 exponential |
| [`run_notebook12.py`](notebooks/run_notebook12.py) | Diagnose the frozen Figure 2/3 fits |

For example, after generating the prerequisite data and results:

```bash
python scripts/notebooks/run_notebook12.py
```

## Independent audits

| Script | Work performed |
| --- | --- |
| [`run_fig1_candidate_loss_audit.py`](audits/run_fig1_candidate_loss_audit.py) | Recalculate the Figure 1 Table III candidate losses |
| [`run_fig1_full13_precision_audit.py`](audits/run_fig1_full13_precision_audit.py) | Check a full-DESI precision-matrix alternative |
| [`run_fig1_v1_objective_audit.py`](audits/run_fig1_v1_objective_audit.py) | Check whether the paper's earlier loss definition explains Table III |
| [`run_cross_figure_loss_objective_audit.py`](audits/run_cross_figure_loss_objective_audit.py) | Compare loss definitions and coefficient rounding across Tables III–V |

The audit scripts evaluate published expressions without rerunning PySR. The cross-figure audit needs the paper PDF at the root path described in the main README.
