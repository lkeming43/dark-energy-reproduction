"""Shared covariance-aware PySR fitting for Figures 2 and 3."""
from pathlib import Path
import json
import platform
import time

import numpy as np
import pandas as pd
import sympy as sp
import yaml
from scipy.optimize import least_squares


def load_pysr_config(config_path):
    config = yaml.safe_load(Path(config_path).read_text())
    if config["loss"]["name"] != "full_covariance_chi2_over_n":
        raise ValueError("Only the full-covariance objective is supported")
    return config


def chi_square(observed, predicted, covariance):
    residual = np.asarray(predicted) - np.asarray(observed)
    return float(residual @ np.linalg.solve(covariance, residual))


def fit_exponential_reference(dataset, initial=(32.0, -0.27)):
    """GLS fit of a*z*exp(b*z), for comparison only."""
    z, y, cov = dataset["z"], dataset["y"], dataset["covariance"]
    whitening = np.linalg.inv(np.linalg.cholesky(cov))
    def residual(params):
        a, b = params
        return whitening @ (a*z*np.exp(b*z) - y)
    fit = least_squares(residual, initial, max_nfev=200)
    a, b = fit.x
    chi2 = float(fit.fun @ fit.fun)
    return {"a": float(a), "b": float(b), "chi2": chi2,
            "loss": chi2/len(z), "success": bool(fit.success),
            "message": str(fit.message)}


def _julia_covariance_loss(binary_path, n):
    """Load an inverse Cholesky matrix once into a Julia closure."""
    path = json.dumps(str(Path(binary_path).resolve()))
    return f'''let
    W = Matrix{{Float64}}(undef, {n}, {n})
    open({path}, "r") do io
        read!(io, W)
    end
    function covariance_loss(tree, dataset::Dataset{{T,L}}, options)::L where {{T,L}}
        prediction, flag = eval_tree_array(tree, dataset.X, options)
        if !flag || any(value -> !isfinite(value), prediction)
            return L(Inf)
        end
        residual = prediction .- dataset.y
        whitened = W * residual
        return L(sum(abs2, whitened) / dataset.n)
    end
    covariance_loss
end'''


def build_pysr_model(config, seed, covariance, output_directory):
    """Create one model and a saved column-major whitening matrix."""
    from pysr import PySRRegressor
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=False)
    whitening = np.linalg.inv(np.linalg.cholesky(covariance))
    factor_path = output_directory / "whitening_matrix_f64_column_major.bin"
    # Transpose then write in C order to match Julia's column-major matrix.
    whitening.T.tofile(factor_path)
    kwargs = dict(config["search"])
    kwargs.update(binary_operators=config["operators"]["binary"],
                  unary_operators=config["operators"]["unary"],
                  constraints={k:tuple(v) for k,v in config["operators"]["constraints"].items()},
                  random_state=seed,
                  loss_function=_julia_covariance_loss(factor_path, len(covariance)),
                  output_directory=str(output_directory / "pysr_raw"),
                  run_id="fit")
    return PySRRegressor(**kwargs), factor_path


def evaluate_sympy(expression, z):
    symbol = sp.Symbol("z")
    expr = sp.sympify(str(expression), locals={"z":symbol})
    with np.errstate(all="ignore"):
        value = np.asarray(sp.lambdify(symbol, expr, modules="numpy")(z), dtype=float)
    return np.broadcast_to(value, np.shape(z)).copy()


def fit_symbolic_distance(dataset, config, seed, output_directory):
    """Run a selected scenario and error convention in a fresh folder."""
    path = Path(output_directory)
    path.parent.mkdir(parents=True, exist_ok=True)
    model, factor_path = build_pysr_model(config, seed, dataset["covariance"], path)
    z, y, cov = dataset["z"], dataset["y"], dataset["covariance"]
    started = time.monotonic()
    model.fit(z.reshape(-1,1), y, variable_names=["z"])
    elapsed = time.monotonic() - started
    selected = model.get_best()
    index = int(selected.name)
    prediction = np.asarray(model.predict(z.reshape(-1,1), index=index), dtype=float).reshape(-1)
    chi2 = chi_square(y, prediction, cov)
    expression = str(selected["sympy_format"])
    settings = config["diagnostics"]
    z_grid = np.linspace(settings["z_grid_min"], settings["z_grid_max"], settings["z_grid_size"])
    grid = evaluate_sympy(expression, z_grid)
    summary = {"scenario":dataset["metadata"]["scenario"],
               "desi_error_version":dataset["metadata"]["desi_error_version"],
               "seed":int(seed), "n":len(y), "selected_index":index,
               "complexity":int(selected["complexity"]), "pysr_loss":float(selected["loss"]),
               "python_chi2":chi2, "python_loss":chi2/len(y),
               "equation":str(selected["equation"]), "sympy_expression":expression,
               "elapsed_seconds":elapsed, "finite_on_grid":bool(np.isfinite(grid).all()),
               "positive_for_z_gt_zero":bool(np.all(grid[1:] > 0)),
               "dm_at_zero":float(grid[0]),
               "satisfies_dm_zero_constraint":bool(abs(grid[0]) <= settings["physical_tolerance"])}
    candidates = model.equations_.copy()
    for column in ("lambda_format", "jax_format", "torch_format"):
        candidates = candidates.drop(columns=column, errors="ignore")
    if "sympy_format" in candidates:
        candidates["sympy_format"] = candidates["sympy_format"].astype(str)
    candidates["selected_by_best"] = candidates.index == index
    return {"summary":summary, "candidates":candidates, "z_grid":z_grid,
            "selected_grid":grid, "observed_prediction":prediction,
            "whitening_file":str(factor_path)}


def export_fit_result(result, output_directory, config, paper_grid, benchmark, data_metadata):
    path = Path(output_directory)
    pd.DataFrame([result["summary"]]).to_csv(path / "summary.csv", index=False)
    result["candidates"].to_csv(path / "candidates.csv", index=True)
    pd.DataFrame({"z":result["z_grid"], "selected":result["selected_grid"],
                  "paper":paper_grid}).to_csv(path / "distance_grid.csv", index=False)
    np.savez_compressed(path / "fit_arrays.npz", z_grid=result["z_grid"],
                        selected_grid=result["selected_grid"],
                        observed_prediction=result["observed_prediction"])
    record = {"summary":result["summary"], "benchmark_azexpbz":benchmark,
              "config":config, "dataset_metadata":data_metadata,
              "software":{"python":platform.python_version(),
                          "numpy":np.__version__, "pandas":pd.__version__}}
    (path / "configuration.json").write_text(json.dumps(record,indent=2,ensure_ascii=False)+"\n")
    return path
