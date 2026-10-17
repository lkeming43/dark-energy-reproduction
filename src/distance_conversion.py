"""Eq. (6), conditional on fixed redshift and calibration."""
import numpy as np


def convert_supernova_to_dm_over_rd(supernova_data, calibration):
    """Central distance, not the mean of the transformed lognormal distribution."""
    obs = np.asarray(supernova_data["observable"], dtype=float)
    kind = supernova_data["observable_kind"]
    if kind == "corrected_magnitude":
        mu = obs - calibration["M"]
    elif kind == "distance_modulus":
        mu = obs
    else:
        raise ValueError(f"Unknown observable: {kind}")
    rd = float(calibration["rd_mpc"])
    z = np.asarray(supernova_data["z"], dtype=float)
    if rd <= 0 or np.any(z <= -1):
        raise ValueError("Invalid distance scale or redshift")
    return 10.0 ** ((mu - 25.0) / 5.0) / ((1 + z) * rd)


def propagate_covariance(supernova_data, covariance, calibration):
    """First-order J C J^T; common r_d uncertainty is kept separate."""
    y = convert_supernova_to_dm_over_rd(supernova_data, calibration)
    covariance = np.asarray(covariance)
    if covariance.shape != (len(y), len(y)):
        raise ValueError("Covariance/row mismatch")
    jac = np.log(10.0) / 5.0 * y
    return jac[:, None] * covariance * jac[None, :]
