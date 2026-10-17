"""Common Figure 2/3 data contract; raw inputs are read only."""
from pathlib import Path
import gzip
import hashlib
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import yaml

from .distance_conversion import convert_supernova_to_dm_over_rd, propagate_covariance


def load_scenario_config(config_path):
    path = Path(config_path).resolve()
    config = yaml.safe_load(path.read_text())
    config["_project_root"] = str(path.parent.parent)
    return config


def _path(config, relative):
    return Path(config["_project_root"]) / relative


def _read_covariance(path, expected_n):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        n = int(handle.readline().strip())
        values = np.loadtxt(handle)
    if n != expected_n or values.size != n * n:
        raise ValueError(f"Covariance dimensions disagree: {path}")
    return values.reshape(n, n)


def _validate_covariance(cov):
    if not np.isfinite(cov).all() or not np.allclose(cov, cov.T, rtol=1e-10, atol=1e-12):
        raise ValueError("Non-finite or asymmetric covariance")
    # Fail rather than silently regularize scientific inputs.
    np.linalg.cholesky(cov)


def load_desi_dm(config):
    cfg = config["desi"]
    frame = pd.read_csv(_path(config, cfg["mean_file"]), sep=r"\s+", comment="#",
                        names=["z", "value", "quantity"])
    cov = np.loadtxt(_path(config, cfg["covariance_file"]))
    if cov.shape != (len(frame), len(frame)):
        raise ValueError("DESI mean/covariance dimension mismatch")
    _validate_covariance(cov)
    ix = np.flatnonzero(frame.quantity.eq(cfg["observable"]).to_numpy())
    if len(ix) != 6:
        raise ValueError("Expected six transverse DESI DR2 measurements")
    selected = frame.iloc[ix]
    np.testing.assert_allclose(selected.z, [.510, .706, .934, 1.321, 1.484, 2.330], atol=5e-7, rtol=0)
    # Table I includes rounding/truncation; retain full release precision.
    np.testing.assert_allclose(selected.value, [13.588,17.351,21.576,27.601,30.512,38.988], atol=.001, rtol=0)
    return {"z": selected.z.to_numpy(), "y": selected.value.to_numpy(),
            "covariance": cov[np.ix_(ix, ix)], "row_index": ix,
            "id": np.array(["LRG1","LRG2","LRG3+ELG1","ELG2","QSO","Ly-alpha"])}


def load_supernova(config):
    cfg = config["supernova"]
    frame = pd.read_csv(_path(config, cfg["data_file"]), sep=cfg["separator"])
    if len(frame) != cfg["expected_raw_rows"]:
        raise ValueError("Unexpected catalog size; check release version")
    col = cfg["columns"]
    if not np.isfinite(frame[[col["redshift"],col["observable"],col["error"]]].to_numpy(dtype=float)).all():
        raise ValueError("Non-finite raw values")
    cut = cfg["selection"]
    mask = np.ones(len(frame), dtype=bool)
    if "redshift_min" in cut:
        mask &= frame[col["redshift"]].to_numpy() >= cut["redshift_min"]
    if "error_max" in cut:
        mask &= frame[col["error"]].to_numpy() <= cut["error_max"]
    ix = np.flatnonzero(mask)
    if len(ix) != cfg["expected_selected_rows"]:
        raise ValueError("Selected sample size disagrees with paper")
    cov = _read_covariance(_path(config, cfg["covariance_file"]), len(frame))
    if cfg["covariance_kind"] == "systematics_only":
        cov[np.diag_indices_from(cov)] += frame[col["error"]].to_numpy() ** 2
    elif cfg["covariance_kind"] != "stat_plus_sys":
        raise ValueError("Unknown covariance convention")
    asymmetry = float(np.max(np.abs(cov - cov.T)))
    # Local Pantheon+ text matrix has <=3e-8 mag^2 asymmetry; origin unverified.
    # Average the two triangles only within an explicit release-level tolerance.
    if asymmetry > cfg.get("symmetry_atol", 0.0):
        raise ValueError(f"Source covariance asymmetry exceeds configured tolerance: {asymmetry}")
    cov = (cov + cov.T) / 2.0
    _validate_covariance(cov)
    selected = frame.iloc[ix]
    return {"z": selected[col["redshift"]].to_numpy(),
            "observable": selected[col["observable"]].to_numpy(),
            "observable_kind": cfg["observable_kind"],
            "covariance": cov[np.ix_(ix, ix)], "row_index": ix,
            "id": selected[col["id"]].astype(str).to_numpy(dtype=str),
            "raw_count": len(frame), "name": cfg["name"], "raw_covariance_max_asymmetry": asymmetry}


def _sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def prepare_combined_dataset(config):
    bao, sn = load_desi_dm(config), load_supernova(config)
    cal = config["supernova"]["calibration"]
    sn_y = convert_supernova_to_dm_over_rd(sn, cal)
    sn_cov = propagate_covariance(sn, sn["covariance"], cal)
    covariance = np.zeros((6 + len(sn_y), 6 + len(sn_y)))
    covariance[:6, :6] = bao["covariance"]
    covariance[6:, 6:] = sn_cov
    sources = [config[g][f] for g, fields in
               [("desi", ("mean_file","covariance_file")), ("supernova", ("data_file","covariance_file"))]
               for f in fields]
    result = {
        "z": np.concatenate([bao["z"],sn["z"]]),
        "y": np.concatenate([bao["y"],sn_y]),
        "covariance": covariance, "sigma": np.sqrt(np.diag(covariance)),
        "source": np.array(["DESI_DR2"]*len(bao["z"]) + [sn["name"]]*len(sn["z"])),
        "row_index": np.concatenate([bao["row_index"],sn["row_index"]]),
        "id": np.concatenate([bao["id"],sn["id"]]).astype(str),
        "sn_observable": sn["observable"], "sn_observable_covariance": sn["covariance"],
        "sn_z": sn["z"],
        "metadata": {
            "schema_version": 1, "scenario": config["scenario"]["id"],
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "raw_sn_count": sn["raw_count"], "selected_sn_count": len(sn["z"]), "desi_count": 6,
            "raw_covariance_max_asymmetry": sn["raw_covariance_max_asymmetry"],
            "symmetry_handling": "average triangles within explicit source tolerance; no eigenvalue regularization",
            "ordering": "DESI first, then selected SN in original row order; no deduplication",
            "observable": "D_M/r_d", "units": "dimensionless",
            "covariance_convention": "Independent BAO/SN blocks; conditional on fixed calibration; first-order Jacobian",
            "redshift_convention": "zHD for fit and 1+z denominator following paper Eq.6",
            "anchor": "none; enforce D_M(0)=0 in model space",
            "rd_uncertainty": "Stored separately in config for shared nuisance MC, not per-SN noise",
            "redshift_uncertainty": "fixed zHD; no duplicate velocity variance",
            "source_sha256": {p: _sha256(_path(config,p)) for p in sources},
            "config": {k:v for k,v in config.items() if not k.startswith("_")},
        },
    }
    result["metadata"]["validation"] = validate_standard_dataset(result)
    return result


def validate_standard_dataset(dataset):
    n = len(dataset["z"])
    for key in ("y","sigma","source","row_index","id"):
        if len(dataset[key]) != n:
            raise ValueError(f"Length mismatch: {key}")
    if dataset["covariance"].shape != (n,n):
        raise ValueError("Combined covariance dimensions disagree")
    for key in ("z","y","sigma"):
        if not np.isfinite(dataset[key]).all() or np.any(dataset[key] <= 0):
            raise ValueError(f"Invalid {key}")
    _validate_covariance(dataset["covariance"])
    np.testing.assert_allclose(dataset["sigma"]**2, np.diag(dataset["covariance"]))
    return {"n": n, "z_min":float(min(dataset["z"])), "z_max":float(max(dataset["z"])),
            "positive_definite":True, "finite":True, "symmetric":True}


def save_standard_dataset(dataset, output_directory):
    """NPZ contains plain arrays and JSON metadata; no pickle required."""
    destination = Path(output_directory)
    destination.mkdir(parents=True, exist_ok=True)
    name = dataset["metadata"]["scenario"]
    metadata = json.dumps(dataset["metadata"], ensure_ascii=False)
    path = destination / f"{name}.npz"
    np.savez_compressed(path, **{k:v for k,v in dataset.items() if k != "metadata"}, metadata_json=metadata)
    (destination / f"{name}.json").write_text(json.dumps(dataset["metadata"],indent=2,ensure_ascii=False)+"\n")
    return path


def load_standard_dataset(path):
    with np.load(path, allow_pickle=False) as archive:
        result = {k:archive[k] for k in archive.files if k != "metadata_json"}
        result["metadata"] = json.loads(str(archive["metadata_json"]))
    validate_standard_dataset(result)
    return result


def with_desi_error_version(dataset, version):
    """Select public errors or a constructed diagonal Table I DM error matrix."""
    if version not in ("public_release", "paper_table"):
        raise ValueError("Unknown DESI error version")
    result = {key: (value.copy() if isinstance(value, np.ndarray) else value)
              for key, value in dataset.items()}
    result["metadata"] = dict(dataset["metadata"])
    result["metadata"]["desi_error_version"] = version
    if version == "paper_table":
        sub = result["covariance"][:6, :6]
        if np.count_nonzero(sub - np.diag(np.diag(sub))):
            raise ValueError("Table I does not provide off-diagonal correlations")
        sigma = np.array([.167, .177, .152, .318, .760, .531])
        result["covariance"][:6, :6] = np.diag(sigma ** 2)
        result["sigma"][:6] = sigma
        result["metadata"]["desi_error_description"] = (
            "Constructed diagonal covariance from rounded Table I DM errors; "
            "not the authors' original covariance")
    else:
        result["metadata"]["desi_error_description"] = "DESI release covariance submatrix"
    validate_standard_dataset(result)
    return result
