"""Run Notebook 11b into a separate, resumable experiment directory.

Example:
    python run_notebook11b.py --draws 20000 --batch-size 500
    python run_notebook11b.py --draws 20000 --batch-size 500 --resume PATH
"""

import argparse
from pathlib import Path

from src.fig23_catalog_bootstrap import run_catalog_bootstrap


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--draws", type=int, default=20000,
                        help="draws per scenario/branch/centering method")
    parser.add_argument("--batch-size", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20261014)
    parser.add_argument("--resume", type=Path)
    args = parser.parse_args()
    if args.draws < 100:
        parser.error("--draws must be at least 100")
    if args.batch_size < 1:
        parser.error("--batch-size must be positive")
    root = Path(__file__).resolve().parent
    prior = (root / "results/11_fig23_mc_and_reconstruction" /
             "fig23_mc_20261010T044029855521Z_5cf83467")
    destination, overview = run_catalog_bootstrap(
        root, prior, n_draws=args.draws, batch_size=args.batch_size,
        base_seed=args.seed, resume_directory=args.resume)
    print(overview.to_string(index=False), flush=True)
    print("Experiment:", destination, flush=True)


if __name__ == "__main__":
    main()
