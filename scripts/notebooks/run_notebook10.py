"""Run Notebook 10 as a foreground, resumable seed ensemble.

Examples:
    python scripts/notebooks/run_notebook10.py
    python scripts/notebooks/run_notebook10.py --resume /absolute/experiment/folder
    python scripts/notebooks/run_notebook10.py --resume /absolute/experiment/folder --seed-count 500

This executes the notebook's code cells in order, overriding only its safety
switch (RUN_FITS) and optional seed count/resume path. Per-seed outputs are
saved by src.seed_stability after each completed fit.
"""

import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resume", type=Path, help="Existing experiment folder")
    parser.add_argument("--seed-count", type=int, default=100,
                        help="Seeds per dataset/error variant (default: 100)")
    args = parser.parse_args()
    if args.seed_count < 1:
        parser.error("--seed-count must be positive")
    if args.resume is not None and not (args.resume / "experiment.json").is_file():
        parser.error("--resume must contain experiment.json")

    path = Path(__file__).resolve().parents[2] / "notebooks/10_fig23_seed_stability.ipynb"
    notebook = json.loads(path.read_text())
    namespace = {"__name__": "__main__"}
    for index, cell in enumerate(notebook["cells"]):
        if cell["cell_type"] != "code":
            continue
        source = "".join(cell["source"])
        if "RUN_FITS = False" in source:
            source = source.replace("RUN_FITS = False", "RUN_FITS = True", 1)
            source = source.replace("N_SEEDS = 100", f"N_SEEDS = {args.seed_count}", 1)
            if args.resume is not None:
                marker = "RESUME_DIRECTORY = globals().get('EXPERIMENT_DIR')"
                if marker not in source:
                    raise RuntimeError("Notebook resume control has changed")
                source = source.replace(marker,
                                        f"RESUME_DIRECTORY = Path({str(args.resume.resolve())!r})", 1)
        exec(compile(source, f"{path.name}:cell_{index}", "exec"), namespace)
    print("Notebook 10 experiment complete:", namespace["EXPERIMENT_DIR"], flush=True)


if __name__ == "__main__":
    main()
