"""Execute notebook 09 and retain cell outputs even if a run stops early."""
from pathlib import Path
import nbformat
from nbclient import NotebookClient


path = Path(__file__).resolve().parent / "notebooks/09_fig23_pysr_reproduction.ipynb"
notebook = nbformat.read(path, as_version=4)
for index, cell in enumerate(notebook.cells):
    cell.setdefault("id", f"fig23-n09-{index:02d}")
nbformat.validate(notebook)
client = NotebookClient(notebook, timeout=7200, kernel_name="jupyter-ml",
                        resources={"metadata":{"path":str(path.parent)}})
try:
    client.execute()
finally:
    nbformat.write(notebook, path)
    for cell in notebook.cells:
        for output in cell.get("outputs", []):
            if output.output_type == "stream":
                print(output.text, flush=True)
print("Notebook 09 executed and saved:", path, flush=True)
