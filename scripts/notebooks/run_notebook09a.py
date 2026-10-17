"""Execute Notebook 09a and preserve outputs if an audit cell fails."""

from pathlib import Path

import nbformat
from nbclient import NotebookClient


path = Path(__file__).resolve().parents[2] / "notebooks/09a_fig23_loss_audit.ipynb"
notebook = nbformat.read(path, as_version=4)
nbformat.validate(notebook)
client = NotebookClient(
    notebook,
    timeout=1200,
    kernel_name="jupyter-ml",
    resources={"metadata": {"path": str(path.parent)}},
)
try:
    client.execute()
finally:
    nbformat.write(notebook, path)
    for cell in notebook.cells:
        for output in cell.get("outputs", []):
            if output.output_type == "stream":
                print(output.text, flush=True)
print("Notebook 09a executed and saved:", path, flush=True)
