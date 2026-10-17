"""Read-only numerical diagnostics for the published Figure 2/3 loss values."""

import numpy as np
import sympy as sp


def evaluate_expression(expression, z):
    """Evaluate a printed paper expression without invoking PySR."""
    symbol = sp.Symbol("z")
    parsed = sp.sympify(expression, locals={"z": symbol})
    with np.errstate(all="ignore"):
        values = np.asarray(sp.lambdify(symbol, parsed, modules="numpy")(z), dtype=float)
    values = np.broadcast_to(values, np.shape(z)).copy()
    if not np.isfinite(values).all():
        raise ValueError(f"Non-finite published expression: {expression}")
    return values


def evaluate_formula_bank(dataset, formulas, covariance):
    """Compute all full-covariance losses with one linear solve."""
    z, y = dataset["z"], dataset["y"]
    residuals = np.column_stack(
        [evaluate_expression(item["expression"], z) - y for item in formulas]
    )
    solved = np.linalg.solve(covariance, residuals)
    chi2 = np.einsum("ij,ij->j", residuals, solved)
    return [
        {
            "formula_id": item["id"],
            "expression": item["expression"],
            "paper_location": item["paper_location"],
            "paper_reported_loss": float(item["paper_reported_loss"]),
            "our_chi2": float(value),
            "our_loss": float(value / len(z)),
            "loss_difference": float(value / len(z) - item["paper_reported_loss"]),
        }
        for item, value in zip(formulas, chi2)
    ]


def block_chi2(dataset, predicted, covariance):
    """Split chi-square only where the BAO/SN cross-block is exactly zero."""
    y = dataset["y"]
    if not np.allclose(covariance[:6, 6:], 0, atol=1e-14, rtol=0):
        raise ValueError("BAO/SN covariance is not block diagonal")
    residual = np.asarray(predicted) - y
    bao = float(residual[:6] @ np.linalg.solve(covariance[:6, :6], residual[:6]))
    sn = float(residual[6:] @ np.linalg.solve(covariance[6:, 6:], residual[6:]))
    return {"bao_chi2": bao, "sn_chi2": sn, "total_chi2": bao + sn,
            "total_loss": (bao + sn) / len(y)}


def diagonal_sn_loss(dataset, predicted, covariance):
    """Counterfactual: keep BAO covariance, discard SN off-diagonals."""
    residual = np.asarray(predicted) - dataset["y"]
    bao = float(residual[:6] @ np.linalg.solve(covariance[:6, :6], residual[:6]))
    sn = float(np.sum(residual[6:] ** 2 / np.diag(covariance)[6:]))
    return (bao + sn) / len(residual)


def amplitude_quadratic(dataset, covariance, shape):
    """For model a*shape(z), return its exact GLS loss quadratic."""
    x = np.asarray(shape, dtype=float)
    y = dataset["y"]
    if x.shape != y.shape:
        raise ValueError("Shape function and data vector differ")
    solved = np.linalg.solve(covariance, np.column_stack((x, y)))
    xx = float(x @ solved[:, 0])
    xy = float(x @ solved[:, 1])
    yy = float(y @ solved[:, 1])
    if xx <= 0:
        raise ValueError("Non-positive quadratic coefficient")
    return {"xx": xx, "xy": xy, "yy": yy, "n": len(y),
            "best_amplitude": xy / xx,
            "best_loss": (yy - xy * xy / xx) / len(y)}


def amplitude_loss(quadratic, amplitude):
    a = np.asarray(amplitude, dtype=float)
    return (quadratic["xx"] * a*a - 2*quadratic["xy"]*a
            + quadratic["yy"]) / quadratic["n"]


def magnitude_space_diagnostic(dataset, predicted, calibration):
    """Alternative residual domain, for diagnosis only, not the Notebook 09 loss."""
    z = dataset["z"][6:]
    distance = np.asarray(predicted)[6:] * (1 + z) * calibration["rd_mpc"]
    if np.any(distance <= 0):
        raise ValueError("Predicted luminosity distance must be positive")
    observable = 25 + 5 * np.log10(distance)
    if "M" in calibration:
        observable = observable + calibration["M"]
    residual = observable - dataset["sn_observable"]
    covariance = dataset["sn_observable_covariance"]
    solved = np.linalg.solve(covariance, np.column_stack((residual, np.ones(len(z)))))
    sn_chi2 = float(residual @ solved[:, 0])
    weighted_shift = -float(np.sum(solved[:, 0]) / np.sum(solved[:, 1]))
    profiled_residual = residual + weighted_shift
    profiled_sn_chi2 = float(profiled_residual @ np.linalg.solve(covariance, profiled_residual))
    return {"sn_observable_chi2": sn_chi2,
            "mean_observable_residual": float(np.mean(residual)),
            "profiled_common_offset_mag": weighted_shift,
            "profiled_sn_chi2": profiled_sn_chi2}
