"""Run the official Cell2fate workflow for the SCOP R package."""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import importlib.metadata
import json
import os
import random
import re
import shutil
import tempfile
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path
import importlib.util

_LOG_MESSAGE_ENV = os.environ.get("SCOP_LOG_MESSAGE_PATH")
_LOG_MESSAGE_PATH = (
    Path(_LOG_MESSAGE_ENV)
    if _LOG_MESSAGE_ENV
    else Path(__file__).resolve().parent / "log_message.py"
)
if not _LOG_MESSAGE_PATH.exists():
    raise ImportError(f"Cannot load log_message module from {_LOG_MESSAGE_PATH}")
_LOG_MESSAGE_SPEC = importlib.util.spec_from_file_location(
    "scop_log_message", _LOG_MESSAGE_PATH
)
_LOG_MESSAGE_MODULE = importlib.util.module_from_spec(_LOG_MESSAGE_SPEC)
_LOG_MESSAGE_SPEC.loader.exec_module(_LOG_MESSAGE_MODULE)
log_message = _LOG_MESSAGE_MODULE.log_message


BACKEND_COMMIT = "c03d1ca0bb963f550001c6070d4986a61ec8456a"
PRODUCER = "RunCell2fate"
RUNNER_SCHEMA_VERSION = 1


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _write_json(value: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=path.name + ".",
        suffix=".tmp",
        dir=path.parent,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


@contextmanager
def _result_lock(result_dir: Path):
    result_dir.mkdir(parents=True, exist_ok=True)
    lock_path = result_dir / ".cell2fate.lock"
    try:
        descriptor = os.open(
            lock_path,
            os.O_CREAT | os.O_EXCL | os.O_WRONLY,
            0o600,
        )
    except FileExistsError as error:
        log_message(
            "Another Cell2fate run is using result_dir. If no run is active, "
            "remove the stale .cell2fate.lock file.",
          message_type="error")
    token = f"{os.getpid()}-{os.urandom(16).hex()}"
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(token + "\n")
        yield token
    finally:
        try:
            observed = lock_path.read_text(encoding="utf-8").strip()
        except OSError:
            observed = None
        if observed == token:
            lock_path.unlink(missing_ok=True)


def _assert_external_lock(result_dir: Path, token: Any) -> None:
    if not isinstance(token, str) or not token:
        log_message("Cell2fate external result lock token is missing", message_type="error")
    lock_path = result_dir / ".cell2fate.lock"
    try:
        observed = lock_path.read_text(encoding="utf-8").strip()
    except OSError as error:
        log_message("Cell2fate external result lock is missing", message_type="error")
    if observed != token:
        log_message("Cell2fate external result lock owner changed", message_type="error")


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
        log_message(
            "Cell2fate installation has no Git provenance; recreate the "
            "cell2fate_env with scop::PrepareEnv(modules = 'cell2fate', force = TRUE)",
          message_type="error")
    metadata = json.loads(direct_url)
    observed = metadata.get("vcs_info", {}).get("commit_id")
    if observed != BACKEND_COMMIT:
        log_message(
            "Cell2fate backend commit mismatch: expected "
            f"{BACKEND_COMMIT}, observed {observed or 'unknown'}. Recreate the "
            "cell2fate_env with scop::PrepareEnv(modules = 'cell2fate', force = TRUE)",
          message_type="error")
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


def _parameter_fingerprint(parameters: dict[str, Any]) -> str:
    encoded = json.dumps(
        parameters,
        allow_nan=False,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _output_paths(result_dir: Path) -> dict[str, Path]:
    paths = {
        "inputs": result_dir / "inputs",
        "model": result_dir / "model",
        "posterior": result_dir / "posterior",
        "tables": result_dir / "tables",
        "manifest": result_dir / "manifest.json",
        "complete": result_dir / ".complete",
        "owner": result_dir / ".cell2fate.json",
        "lock": result_dir / ".cell2fate.lock",
        "logs": result_dir / "logs",
    }
    return paths


def _assert_safe_result_paths(result_dir: Path) -> None:
    if result_dir.is_symlink():
        log_message("Cell2fate result_dir must not be a symbolic link", message_type="error")
    root = result_dir.resolve(strict=False)
    paths = _output_paths(result_dir)
    managed = tuple(paths.values()) + (
        paths["inputs"] / "input.h5ad",
        paths["model"] / "model.pt",
        paths["model"] / "adata.h5ad",
        paths["posterior"] / "cell2fate_posterior.h5ad",
        paths["tables"] / "cell_metadata.csv",
        paths["tables"] / "velocity.csv",
    )
    for path in managed:
        if path.is_symlink():
            log_message(
                f"Cell2fate managed path must not be a symbolic link: {path}",
              message_type="error")
        try:
            path.resolve(strict=False).relative_to(root)
        except ValueError as error:
            log_message(
                f"Cell2fate managed path escapes result_dir: {path}",
              message_type="error")


def _ensure_output_directories(paths: dict[str, Path]) -> None:
    for key in ("inputs", "posterior", "tables"):
        paths[key].mkdir(parents=True, exist_ok=True)


def _expected_outputs(
    paths: dict[str, Path],
    require_velocity: bool,
) -> tuple[Path, ...]:
    outputs = list(_artifact_paths(paths, require_velocity).values())
    outputs.extend([
        paths["manifest"],
        paths["complete"],
        paths["owner"],
    ])
    return tuple(outputs)


def _artifact_paths(
    paths: dict[str, Path],
    require_velocity: bool,
) -> dict[str, Path]:
    artifacts = {
        "inputs/input.h5ad": paths["inputs"] / "input.h5ad",
        "model/model.pt": paths["model"] / "model.pt",
        "model/adata.h5ad": paths["model"] / "adata.h5ad",
        "posterior/cell2fate_posterior.h5ad": (
            paths["posterior"] / "cell2fate_posterior.h5ad"
        ),
        "tables/cell_metadata.csv": paths["tables"] / "cell_metadata.csv",
    }
    if require_velocity:
        artifacts["tables/velocity.csv"] = paths["tables"] / "velocity.csv"
    return artifacts


def _artifact_records(
    paths: dict[str, Path],
    require_velocity: bool,
) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    for name, path in _artifact_paths(paths, require_velocity).items():
        size = path.stat().st_size
        if size <= 0:
            log_message(f"Cell2fate artifact is empty: {name}", message_type="error")
        records[name] = {"size": size, "sha256": _sha256(path)}
    return records


def _artifacts_match(
    manifest: dict[str, Any],
    paths: dict[str, Path],
    require_velocity: bool,
) -> bool:
    expected = _artifact_paths(paths, require_velocity)
    observed = manifest.get("artifacts")
    if not isinstance(observed, dict) or set(observed) != set(expected):
        return False
    for name, path in expected.items():
        record = observed.get(name)
        if not isinstance(record, dict) or not path.is_file():
            return False
        try:
            size = path.stat().st_size
        except OSError:
            return False
        try:
            digest = _sha256(path)
        except OSError:
            return False
        if record.get("size") != size or size <= 0 or record.get("sha256") != digest:
            return False
    return True


def _is_owned_result(paths: dict[str, Path]) -> bool:
    if not paths["owner"].exists():
        return False
    try:
        owner = _read_json(paths["owner"])
    except (OSError, ValueError, TypeError):
        return False
    if not isinstance(owner, dict):
        return False
    return (
        owner.get("producer") == PRODUCER
        and owner.get("runner_schema_version") == RUNNER_SCHEMA_VERSION
    )


def _remove_outputs(paths: dict[str, Path]) -> None:
    for key in ("inputs", "model", "posterior", "tables"):
        target = paths[key]
        if target.is_dir():
            shutil.rmtree(target)
        elif target.exists():
            target.unlink()
    for key in ("manifest", "complete"):
        paths[key].unlink(missing_ok=True)


@contextmanager
def _staged_output_paths(result_dir: Path):
    stage_dir = Path(
        tempfile.mkdtemp(prefix=".cell2fate-stage-", dir=result_dir)
    )
    try:
        yield _output_paths(stage_dir)
    finally:
        shutil.rmtree(stage_dir, ignore_errors=True)


def _stale_stage_paths(result_dir: Path) -> tuple[Path, ...]:
    return tuple(result_dir.glob(".cell2fate-stage-*"))


def _cleanup_stale_stages(result_dir: Path, lock_token: str) -> None:
    for stage_dir in _stale_stage_paths(result_dir):
        _assert_external_lock(result_dir, lock_token)
        if stage_dir.is_symlink() or stage_dir.is_file():
            stage_dir.unlink(missing_ok=True)
        elif stage_dir.is_dir():
            shutil.rmtree(stage_dir)
        _assert_external_lock(result_dir, lock_token)


def _publish_staged_outputs(
    staged: dict[str, Path],
    final: dict[str, Path],
    result_dir: Path,
    lock_token: str,
) -> None:
    _assert_external_lock(result_dir, lock_token)
    _assert_safe_result_paths(result_dir)
    _remove_outputs(final)
    _assert_external_lock(result_dir, lock_token)
    for key in ("inputs", "model", "posterior", "tables", "manifest"):
        os.replace(staged[key], final[key])
        _assert_external_lock(result_dir, lock_token)
    os.replace(staged["complete"], final["complete"])
    _assert_external_lock(result_dir, lock_token)


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
    try:
        manifest = _read_json(paths["manifest"])
    except (OSError, TypeError, ValueError):
        return False
    if not isinstance(manifest, dict):
        return False
    manifest_parameters = manifest.get("parameters")
    if not isinstance(manifest_parameters, dict):
        return False
    observed_parameter_fingerprint = manifest.get("parameters_sha256")
    if not isinstance(observed_parameter_fingerprint, str):
        try:
            observed_parameter_fingerprint = _parameter_fingerprint(
                manifest_parameters
            )
        except (TypeError, ValueError):
            return False
    return (
        manifest.get("producer") == PRODUCER
        and manifest.get("runner_schema_version") == RUNNER_SCHEMA_VERSION
        and manifest.get("input_sha256") == fingerprint
        and observed_parameter_fingerprint == _parameter_fingerprint(parameters)
        and manifest.get("backend_commit") == BACKEND_COMMIT
        and manifest.get("status") == "complete"
        and _artifacts_match(
            manifest,
            paths,
            require_velocity=require_velocity,
        )
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
    n_modules: int,
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
        log_message("Cell2fate posterior did not contain 'Time (hours)'", message_type="error")
    uncertainty = [
        column for column in selected if column.startswith("Time Uncertainty")
    ]
    if len(uncertainty) != 1:
        log_message(
            "Cell2fate posterior did not contain one time uncertainty column",
          message_type="error")
    expected_modules = {
        f"Module {module} {metric}"
        for module in range(n_modules)
        for metric in ("Activation", "State")
    }
    observed_modules = {
        column
        for column in selected
        if re.fullmatch(r"Module\s+\d+\s+(Activation|State)", column)
    }
    if observed_modules != expected_modules:
        log_message(
            "Cell2fate posterior module summaries do not match the fitted modules",
          message_type="error")
    expected_columns = {"Time (hours)", uncertainty[0], *expected_modules}
    if set(selected) != expected_columns:
        log_message("Cell2fate posterior contained unexpected summary columns", message_type="error")
    output = adata.obs.loc[:, selected].copy()
    continuous = [
        column
        for column in selected
        if column == "Time (hours)"
        or column.startswith("Time Uncertainty")
        or re.fullmatch(r"Module\s+\d+\s+Activation", column)
    ]
    for column in continuous:
        try:
            output[column] = pd.to_numeric(output[column], errors="raise").astype(float)
        except (TypeError, ValueError) as error:
            log_message(
                f"Cell2fate posterior metadata {column!r} must be numeric",
              message_type="error")
    output.columns = [_safe_metadata_name(column, prefix) for column in selected]
    if output.columns.duplicated().any():
        log_message("Cell2fate posterior metadata names are not unique", message_type="error")
    for column in output.select_dtypes(include=[np.number]).columns:
        if not np.isfinite(output[column].to_numpy(dtype=float)).all():
            log_message(
                f"Cell2fate posterior metadata {column!r} contains non-finite values",
              message_type="error")
    output.index = output.index.astype(str)
    return output


def _velocity_table(adata: ad.AnnData) -> pd.DataFrame:
    if "velocity" not in adata.layers:
        log_message("Cell2fate posterior did not contain a velocity layer", message_type="error")
    values = adata.layers["velocity"]
    if hasattr(values, "toarray"):
        values = values.toarray()
    values = np.asarray(values, dtype=float)
    if values.shape != (adata.n_obs, adata.n_vars):
        log_message("Cell2fate velocity dimensions do not match the posterior", message_type="error")
    if not np.isfinite(values).all():
        log_message("Cell2fate velocity contains non-finite values", message_type="error")
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
        log_message(
            "Cell2fate posterior is missing values required for velocity: "
            + ", ".join(missing),
          message_type="error")
    expression = np.asarray(means["mu_expression"], dtype=float)
    beta = np.asarray(means["beta_g"], dtype=float)
    gamma = np.asarray(means["gamma_g"], dtype=float)
    if expression.shape != (adata.n_obs, adata.n_vars, 2):
        log_message(
            "Cell2fate posterior expression dimensions do not match the data",
          message_type="error")
    velocity = beta * expression[..., 0] - gamma * expression[..., 1]
    if velocity.shape != (adata.n_obs, adata.n_vars):
        log_message(
            "Cell2fate posterior kinetic parameters do not match the data",
          message_type="error")
    if not np.isfinite(velocity).all():
        log_message("Cell2fate posterior velocity contains non-finite values", message_type="error")
    adata.layers["velocity"] = velocity


def run(config: dict[str, Any]) -> None:
    result_dir = Path(config["result_dir"])
    _assert_safe_result_paths(result_dir)
    external_token = config.get("lock_token")
    if external_token is not None:
        _assert_external_lock(result_dir, external_token)
        _run_locked(config, lock_token=external_token)
        _assert_external_lock(result_dir, external_token)
        return
    with _result_lock(result_dir) as lock_token:
        _assert_safe_result_paths(result_dir)
        _run_locked(config, lock_token=lock_token)


def _run_locked(config: dict[str, Any], lock_token: str) -> None:
    import cell2fate as c2f

    observed_commit = _installed_backend_commit()
    result_dir = Path(config["result_dir"])
    input_path = Path(config["input_path"])
    paths = _output_paths(result_dir)
    fingerprint = _sha256(input_path)
    parameters = _stable_parameters(config)
    parameter_fingerprint = _parameter_fingerprint(parameters)
    store_velocity = bool(config.get("store_velocity", False))
    owned_result = _is_owned_result(paths)
    stale_stages = _stale_stage_paths(result_dir)
    if stale_stages and not owned_result:
        log_message(
            "result_dir contains unowned Cell2fate staging paths; refusing "
            "to remove them",
          message_type="error")
    if owned_result:
        _cleanup_stale_stages(result_dir, lock_token)

    if bool(config.get("resume", True)) and _can_resume(
        paths,
        fingerprint,
        parameters,
        require_velocity=store_velocity,
    ):
        _assert_external_lock(result_dir, lock_token)
        manifest = _read_json(paths["manifest"])
        manifest["request_id"] = config.get("request_id")
        manifest["reused"] = True
        manifest["parameters_sha256"] = parameter_fingerprint
        _write_json(manifest, paths["manifest"])
        log_message("Reuse matching Cell2fate result")
        return
    managed_paths = tuple(
        paths[key]
        for key in ("inputs", "model", "posterior", "tables", "manifest", "complete")
    )
    existing = any(path.exists() for path in managed_paths)
    if existing and not owned_result:
        log_message(
            "result_dir contains paths managed by Cell2fate but has no valid "
            "RunCell2fate ownership marker; refusing to replace them",
          message_type="error")
    if existing and not bool(config.get("overwrite", False)):
        log_message(
            "Existing Cell2fate artifacts do not match the current input or "
            "parameters; set overwrite=TRUE to replace them",
          message_type="error")
    result_dir.mkdir(parents=True, exist_ok=True)
    _write_json(
        {
            "producer": PRODUCER,
            "runner_schema_version": RUNNER_SCHEMA_VERSION,
            "backend_commit": observed_commit,
        },
        paths["owner"],
    )
    with _staged_output_paths(result_dir) as staged:
        _ensure_output_directories(staged)

        seed = int(config.get("seed", 1))
        _set_seed(seed)
        adata = ad.read_h5ad(input_path)
        if "unspliced" not in adata.layers:
            log_message("Input AnnData is missing the 'unspliced' layer", message_type="error")
        adata.layers["spliced"] = adata.X.copy()
        original_columns = set(str(column) for column in adata.obs.columns)

        cluster_column = str(config["cluster_column"])
        if cluster_column not in adata.obs:
            log_message(
                f"Input AnnData is missing cluster column {cluster_column!r}",
              message_type="error")
        cells_per_cluster = config.get("cells_per_cluster")
        if cells_per_cluster is None:
            cells_per_cluster = adata.n_obs
        n_var_genes = min(int(config["n_var_genes"]), adata.n_vars)
        adata = c2f.utils.get_training_data(
            adata,
            cells_per_cluster=int(cells_per_cluster),
            cluster_column=cluster_column,
            remove_clusters=[
                str(value) for value in config.get("remove_clusters", [])
            ],
            min_shared_counts=int(config["min_shared_counts"]),
            n_var_genes=n_var_genes,
        )
        if adata.n_obs < 2 or adata.n_vars < 2:
            log_message(
                "Cell2fate preprocessing retained fewer than two cells or genes",
              message_type="error")

        n_modules = config.get("n_modules")
        if n_modules is None:
            n_modules = int(c2f.utils.get_max_modules(adata))
        if int(n_modules) < 1:
            log_message("Cell2fate requires at least one module", message_type="error")

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

        staged["model"].parent.mkdir(parents=True, exist_ok=True)
        model.save(str(staged["model"]), overwrite=True, save_anndata=True)
        posterior_path = staged["posterior"] / "cell2fate_posterior.h5ad"
        adata.write_h5ad(posterior_path)
        shutil.copy2(input_path, staged["inputs"] / "input.h5ad")

        metadata = _posterior_metadata(
            adata,
            original_columns=original_columns,
            prefix=str(config.get("prefix", "Cell2fate")),
            n_modules=int(n_modules),
        )
        metadata.to_csv(staged["tables"] / "cell_metadata.csv")
        if store_velocity:
            velocity = _velocity_table(adata)
            velocity.to_csv(staged["tables"] / "velocity.csv")

        manifest = {
            "producer": PRODUCER,
            "runner_schema_version": RUNNER_SCHEMA_VERSION,
            "status": "complete",
            "request_id": config.get("request_id"),
            "reused": False,
            "backend_commit": observed_commit,
            "input_sha256": fingerprint,
            "parameters": parameters,
            "parameters_sha256": parameter_fingerprint,
            "versions": _versions(),
            "n_modules": int(n_modules),
            "cells": adata.obs_names.astype(str).tolist(),
            "features": adata.var_names.astype(str).tolist(),
            "artifacts": _artifact_records(
                staged,
                require_velocity=store_velocity,
            ),
        }
        _write_json(manifest, staged["manifest"])
        staged["complete"].write_text("complete\n", encoding="utf-8")
        _publish_staged_outputs(
            staged,
            final=paths,
            result_dir=result_dir,
            lock_token=lock_token,
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    run(_read_json(Path(args.config)))


if __name__ == "__main__":
    main()
