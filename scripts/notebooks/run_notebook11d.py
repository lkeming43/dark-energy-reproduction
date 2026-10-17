"""Run the Figure 2 selected-exponential catalogue MC in a fresh folder."""

import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--draws", type=int, default=20000,
                        help="draws per catalogue-centering method")
    parser.add_argument("--batch-size", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20261016)
    args = parser.parse_args()
    if args.draws < 100 or args.batch_size < 1:
        parser.error("--draws must be >=100 and --batch-size must be positive")
    notebook_path = Path(__file__).resolve().parents[2] / "notebooks/11d_fig2_selected_exponential_mc.ipynb"
    notebook = json.loads(notebook_path.read_text())
    namespace = {"__name__": "__main__"}
    for index, cell in enumerate(notebook["cells"]):
        if cell["cell_type"] != "code":
            continue
        source = "".join(cell["source"])
        if "RUN_ANALYSIS = False" in source:
            source = source.replace("RUN_ANALYSIS = False", "RUN_ANALYSIS = True", 1)
            source = source.replace("N_DRAWS = 20000", f"N_DRAWS = {args.draws}", 1)
            source = source.replace("BATCH_SIZE = 500", f"BATCH_SIZE = {args.batch_size}", 1)
            source = source.replace("MC_SEED = 20261016", f"MC_SEED = {args.seed}", 1)
        exec(compile(source, f"{notebook_path.name}:cell_{index}", "exec"), namespace)
    print("Notebook 11d experiment complete:", namespace["EXPERIMENT_DIR"], flush=True)


if __name__ == "__main__":
    main()
