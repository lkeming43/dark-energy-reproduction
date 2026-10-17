"""Run Notebook 11c into a fresh experiment folder without changing Notebook 11."""

import json
from pathlib import Path


def main():
    notebook_path = Path(__file__).resolve().parents[2] / "notebooks/11c_fig2_balanced_seed_audit.ipynb"
    notebook = json.loads(notebook_path.read_text())
    namespace = {"__name__": "__main__"}
    for index, cell in enumerate(notebook["cells"]):
        if cell["cell_type"] != "code":
            continue
        source = "".join(cell["source"])
        if "RUN_ANALYSIS = False" in source:
            source = source.replace("RUN_ANALYSIS = False", "RUN_ANALYSIS = True", 1)
        exec(compile(source, f"{notebook_path.name}:cell_{index}", "exec"), namespace)
    print("Notebook 11c experiment complete:", namespace["EXPERIMENT_DIR"], flush=True)


if __name__ == "__main__":
    main()
