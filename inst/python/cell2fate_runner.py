"""Run the official Cell2fate workflow for the SCOP R package."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import random
import re
import shutil
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd


BACKEND_COMMIT = "c03d1ca0bb963f550001c6070d4986a61ec8456a"
PRODUCER = "RunCell2fate"
RUNNER_SCHEMA_VERSION = 1


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _write_json(value: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _versions() -> dict[str, str]:
    packages = (
        "cell2fate",
        "scvi-tools",
        "anndata",
        "scanpy",
        "scvelo",
        "torch",
        "jax",
        "jaxlib",
        "numpy",
        "pandas",
        "scipy",
    )
    versions: dict[str, str] = {}
    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "missing"
    return versions


def _installed_backend_commit() -> str:
    distribution = importlib.metadata.distribution("cell2fate")
    direct_url = distribution.read_text("direct_url.json")
    if direct_url is None:
        raise RuntimeError(
            "Cell2fate installation has no Git provenance; recreate the "
            "cell2fate_env with scop::PrepareEnv(modules = 'cell2fate', force = TRUE)"
        )
    metadata = json.loads(direct_url)
    observed = metadata.get("vcs_info", {}).get("commit_id")
    if observed != BACKEND_COMMIT:
        raise RuntimeError(
            "Cell2fate backend commit mismatch: expected "
            f"{BACKEND_COMMIT}, observed {observed or 'unknown'}. Recreate the "
            "cell2fate_env with scop::PrepareEnv(modules = 'cell2fate', force = TRUE)"
        )
    return observed


def _stable_parameters(config: dict[str, Any]) -> dict[str, Any]:
    keys = (
        "cluster_column",
        "remove_clusters",
        "cells_per_cluster",
        "min_shared_counts",
        "n_var_genes",
        "n_modules",
        "model_params",
        "train_params",
        "posterior_params",
        "seed",
        "prefix",
        "store_velocity",
    )
    return {key: config.get(key) for key in keys}


def _output_paths(result_dir: Path) -> dict[str, Path]:
    paths = {
        "inputs": result_dir / "inputs",
        "model": result_dir / "model",
        "posterior": result_dir / "posterior",
        "tables": result_dir / "tables",
        "manifest": result_dir / "manifest.json",
        "complete": result_dir / ".complete",
        "owner": result_dir / ".cell2fate.json",
    }
    return paths


def _ensure_output_directories(paths: dict[str, Path]) -> None:
    for key in ("inputs", "posterior", "tables"):
        paths[key].mkdir(parents=True, exist_ok=True)


def _expected_outputs(
    paths: dict[str, Path],
    require_velocity: bool,
) -> tuple[Path, ...]:
    outputs = [
        paths["model"],
        paths["posterior"] / "cell2fate_posterior.h5ad",
        paths["tables"] / "cell_metadata.csv",
        paths["manifest"],
        paths["complete"],
        paths["owner"],
    ]
    if require_velocity:
        outputs.append(paths["tables"] / "velocity.csv")
    return tuple(outputs)


def _is_owned_result(paths: dict[str, Path]) -> bool:
    if not paths["owner"].exists():
        return False
    try:
        owner = _read_json(paths["owner"])
    except (OSError, ValueError, TypeError):
        return False
    return (
        owner.get("producer") == PRODUCER
        and owner.get("runner_schema_version") == RUNNER_SCHEMA_VERSION
        and owner.get("backend_commit") == BACKEND_COMMIT
    )


def _clear_outputs(paths: dict[str, Path]) -> None:
    for key in ("inputs", "model", "posterior", "tables"):
        target = paths[key]
        if target.exists():
            shutil.rmtree(target)
    for key in ("manifest", "complete"):
        paths[key].unlink(missing_ok=True)
    _ensure_output_directories(paths)


def _can_resume(
    paths: dict[str, Path],
    fingerprint: str,
    parameters: dict[str, Any],
    require_velocity: bool,
) -> bool:
    if not _is_owned_result(paths):
        return False
    if not all(
        path.exists()
        for path in _expected_outputs(paths, require_velocity=require_velocity)
    ):
        return False
    manifest = _read_json(paths["manifest"])
    return (
        manifest.get("producer") == PRODUCER
        and manifest.get("runner_schema_version") == RUNNER_SCHEMA_VERSION
        and manifest.get("input_sha256") == fingerprint
        and manifest.get("parameters") == parameters
        and manifest.get("backend_commit") == BACKEND_COMMIT
        and manifest.get("status") == "complete"
    )


def _set_seed(seed: int) -> None:
    import scvi
    import torch

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    scvi.settings.seed = seed


def _safe_metadata_name(name: str, prefix: str) -> str:
    if name == "Time (hours)":
        suffix = "time"
    elif name.startswith("Time Uncertainty"):
        suffix = "time_uncertainty"
    else:
        match = re.fullmatch(r"Module\s+(\d+)\s+(Activation|State)", name)
        if match is None:
            suffix = re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_").lower()
        else:
            suffix = f"module_{match.group(1)}_{match.group(2).lower()}"
    return f"{prefix}_{suffix}"


def _posterior_metadata(
    adata: ad.AnnData,
    original_columns: set[str],
    prefix: str,
) -> pd.DataFrame:
    selected = [
        column
        for column in adata.obs.columns
        if column not in original_columns
        and (
            column.startswith("Time ")
            or re.fullmatch(r"Module\s+\d+\s+(Activation|State)", column)
        )
    ]
    if "Time (hours)" not in selected:
        raise KeyError("Cell2fate posterior did not contain 'Time (hours)'")
    output = adata.obs.loc[:, selected].copy()
    output.columns = [_safe_metadata_name(column, prefix) for column in selected]
    if output.columns.duplicated().any():
        raise ValueError("Cell2fate posterior metadata names are not unique")
    for column in output.select_dtypes(include=[np.number]).columns:
        if not np.isfinite(output[column].to_numpy(dtype=float)).all():
            raise ValueError(
                f"Cell2fate posterior metadata {column!r} contains non-finite values"
            )
    output.index = output.index.astype(str)
    return output


def _velocity_table(adata: ad.AnnData) -> pd.DataFrame:
    if "velocity" not in adata.layers:
        raise KeyError("Cell2fate posterior did not contain a velocity layer")
    values = adata.layers["velocity"]
    if hasattr(values, "toarray"):
        values = values.toarray()
    values = np.asarray(values, dtype=float)
    if values.shape != (adata.n_obs, adata.n_vars):
        raise ValueError("Cell2fate velocity dimensions do not match the posterior")
    if not np.isfinite(values).all():
        raise ValueError("Cell2fate velocity contains non-finite values")
    return pd.DataFrame(
        values,
        index=adata.obs_names.astype(str),
        columns=adata.var_names.astype(str),
    )


def _add_velocity_from_posterior(adata: ad.AnnData, model: Any) -> None:
    means = model.samples.get("post_sample_means", {})
    required = ("mu_expression", "beta_g", "gamma_g")
    missing = [key for key in required if key not in means]
    if missing:
        raise KeyError(
            "Cell2fate posterior is missing values required for velocity: "
            + ", ".join(missing)
        )
    expression = np.asarray(means["mu_expression"], dtype=float)
    beta = np.asarray(means["beta_g"], dtype=float)
    gamma = np.asarray(means["gamma_g"], dtype=float)
    if expression.shape != (adata.n_obs, adata.n_vars, 2):
        raise ValueError(
            "Cell2fate posterior expression dimensions do not match the data"
        )
    velocity = beta * expression[..., 0] - gamma * expression[..., 1]
    if velocity.shape != (adata.n_obs, adata.n_vars):
        raise ValueError(
            "Cell2fate posterior kinetic parameters do not match the data"
        )
    adata.layers["velocity"] = velocity


def run(config: dict[str, Any]) -> None:
    import cell2fate as c2f

    observed_commit = _installed_backend_commit()
    result_dir = Path(config["result_dir"])
    input_path = Path(config["input_path"])
    paths = _output_paths(result_dir)
    fingerprint = _sha256(input_path)
    parameters = _stable_parameters(config)
    store_velocity = bool(config.get("store_velocity", False))

    if bool(config.get("resume", True)) and _can_resume(
        paths,
        fingerprint,
        parameters,
        require_velocity=store_velocity,
    ):
        print("Reuse matching Cell2fate result")
        return
    managed_paths = tuple(
        paths[key]
        for key in ("inputs", "model", "posterior", "tables", "manifest", "complete")
    )
    existing = any(path.exists() for path in managed_paths)
    if existing and not _is_owned_result(paths):
        raise RuntimeError(
            "result_dir contains paths managed by Cell2fate but has no valid "
            "RunCell2fate ownership marker; refusing to replace them"
        )
    if existing and not bool(config.get("overwrite", False)):
        raise RuntimeError(
            "Existing Cell2fate artifacts do not match the current input or "
            "parameters; set overwrite=TRUE to replace them"
        )
    if existing:
        _clear_outputs(paths)
    result_dir.mkdir(parents=True, exist_ok=True)
    _write_json(
        {
            "producer": PRODUCER,
            "runner_schema_version": RUNNER_SCHEMA_VERSION,
            "backend_commit": observed_commit,
        },
        paths["owner"],
    )
    _ensure_output_directories(paths)

    seed = int(config.get("seed", 1))
    _set_seed(seed)
    adata = ad.read_h5ad(input_path)
    if "unspliced" not in adata.layers:
        raise KeyError("Input AnnData is missing the 'unspliced' layer")
    adata.layers["spliced"] = adata.X.copy()
    original_columns = set(str(column) for column in adata.obs.columns)

    cluster_column = str(config["cluster_column"])
    if cluster_column not in adata.obs:
        raise KeyError(f"Input AnnData is missing cluster column {cluster_column!r}")
    cells_per_cluster = config.get("cells_per_cluster")
    if cells_per_cluster is None:
        cells_per_cluster = adata.n_obs
    n_var_genes = min(int(config["n_var_genes"]), adata.n_vars)
    adata = c2f.utils.get_training_data(
        adata,
        cells_per_cluster=int(cells_per_cluster),
        cluster_column=cluster_column,
        remove_clusters=[str(value) for value in config.get("remove_clusters", [])],
        min_shared_counts=int(config["min_shared_counts"]),
        n_var_genes=n_var_genes,
    )
    if adata.n_obs < 2 or adata.n_vars < 2:
        raise ValueError(
            "Cell2fate preprocessing retained fewer than two cells or genes"
        )

    n_modules = config.get("n_modules")
    if n_modules is None:
        n_modules = int(c2f.utils.get_max_modules(adata))
    if int(n_modules) < 1:
        raise ValueError("Cell2fate requires at least one module")

    c2f.Cell2fate_DynamicalModel.setup_anndata(
        adata,
        spliced_label="spliced",
        unspliced_label="unspliced",
    )
    model = c2f.Cell2fate_DynamicalModel(
        adata,
        n_modules=int(n_modules),
        **dict(config.get("model_params") or {}),
    )
    model.train(**dict(config.get("train_params") or {}))
    posterior_params = dict(config.get("posterior_params") or {})
    posterior_params.setdefault("num_samples", 30)
    posterior_params.setdefault("batch_size", None)
    posterior_params.setdefault("use_gpu", False)
    posterior_params.setdefault("return_samples", False)
    adata = model.export_posterior(
        adata,
        sample_kwargs=posterior_params,
    )
    adata = model.compute_module_summary_statistics(adata)
    _add_velocity_from_posterior(adata, model)

    paths["model"].parent.mkdir(parents=True, exist_ok=True)
    model.save(str(paths["model"]), overwrite=True, save_anndata=True)
    posterior_path = paths["posterior"] / "cell2fate_posterior.h5ad"
    adata.write_h5ad(posterior_path)
    shutil.copy2(input_path, paths["inputs"] / "input.h5ad")

    metadata = _posterior_metadata(
        adata,
        original_columns=original_columns,
        prefix=str(config.get("prefix", "Cell2fate")),
    )
    metadata.to_csv(paths["tables"] / "cell_metadata.csv")
    if store_velocity:
        velocity = _velocity_table(adata)
        velocity.to_csv(paths["tables"] / "velocity.csv")

    manifest = {
        "producer": PRODUCER,
        "runner_schema_version": RUNNER_SCHEMA_VERSION,
        "status": "complete",
        "backend_commit": observed_commit,
        "input_sha256": fingerprint,
        "parameters": parameters,
        "versions": _versions(),
        "n_modules": int(n_modules),
        "cells": adata.obs_names.astype(str).tolist(),
        "features": adata.var_names.astype(str).tolist(),
    }
    _write_json(manifest, paths["manifest"])
    paths["complete"].write_text("complete\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    run(_read_json(Path(args.config)))


if __name__ == "__main__":
    main()
