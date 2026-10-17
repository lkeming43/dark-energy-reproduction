"""Run Notebook 11 into a fresh experiment folder.

Usage:
    /path/to/python scripts/notebooks/run_notebook11.py --draws 10000

The notebook defaults to inspection mode in VS Code; this runner changes only
RUN_ANALYSIS and optionally the Monte Carlo draw count and random seed.
"""

import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--draws", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261011)
    args = parser.parse_args()
    if args.draws < 100:
        parser.error("--draws must be at least 100")
    notebook_path = Path(__file__).resolve().parents[2] / "notebooks/11_fig23_mc_and_reconstruction.ipynb"
    notebook = json.loads(notebook_path.read_text())
    namespace = {"__name__": "__main__"}
    for index, cell in enumerate(notebook["cells"]):
        if cell["cell_type"] != "code":
            continue
        source = "".join(cell["source"])
        if "RUN_ANALYSIS = False" in source:
            source = source.replace("RUN_ANALYSIS = False", "RUN_ANALYSIS = True", 1)
            source = source.replace("N_MC = 10000", f"N_MC = {args.draws}", 1)
            source = source.replace("MC_SEED = 20261011", f"MC_SEED = {args.seed}", 1)
        exec(compile(source, f"{notebook_path.name}:cell_{index}", "exec"), namespace)
    print("Notebook 11 experiment complete:", namespace["EXPERIMENT_DIR"], flush=True)


if __name__ == "__main__":
    main()
