"""Run the official cell2location workflow for the SCOP R package.

The runner deliberately communicates through h5ad/CSV/JSON files so R never
depends on live Python model objects. It supports stage-level resume guarded by
content hashes and parameter manifests.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import shutil
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd


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


def _normalise_optional_list(value: Any) -> list[str] | None:
    if value is None:
        return None
    if isinstance(value, str):
        return [value]
    return [str(item) for item in value]


def _normalise_train_params(value: Any) -> dict[str, Any]:
    params = dict(value or {})
    if "devices" in params:
        if "device" in params:
            raise ValueError("Use only one of 'device' or legacy 'devices' in cell2location train parameters")
        params["device"] = params.pop("devices")
    return params


def _versions() -> dict[str, str]:
    packages = (
        "cell2location",
        "scvi-tools",
        "scanpy",
        "anndata",
        "numpy",
        "pandas",
        "scipy",
        "torch",
    )
    versions: dict[str, str] = {}
    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "missing"
    return versions


def _stable_parameters(config: dict[str, Any]) -> dict[str, Any]:
    keys = (
        "reference_label",
        "spatial_batch",
        "reference_batch",
        "reference_covariates",
        "N_cells_per_location",
        "detection_alpha",
        "gene_filter_params",
        "reference_train_params",
        "spatial_train_params",
        "reference_posterior_params",
        "spatial_posterior_params",
    )
    return {key: config.get(key) for key in keys}


def _fingerprints(config: dict[str, Any]) -> dict[str, Any]:
    spatial_path = Path(config["spatial_path"])
    reference_path = config.get("reference_path")
    signatures_path = config.get("signatures_path")
    return {
        "spatial_sha256": _sha256(spatial_path),
        "reference_sha256": _sha256(Path(reference_path)) if reference_path else None,
        "signatures_sha256": _sha256(Path(signatures_path)) if signatures_path else None,
    }


def _manifest_matches(
    manifest: dict[str, Any] | None,
    fingerprints: dict[str, Any],
    parameters: dict[str, Any],
    *,
    reference_only: bool,
) -> bool:
    if not manifest:
        return False
    old_fingerprints = manifest.get("fingerprints", {})
    if reference_only:
        fingerprint_keys = ("reference_sha256", "signatures_sha256")
        parameter_keys = (
            "reference_label",
            "reference_batch",
            "reference_covariates",
            "gene_filter_params",
            "reference_train_params",
            "reference_posterior_params",
            "runtime_versions",
        )
    else:
        fingerprint_keys = tuple(fingerprints)
        parameter_keys = tuple(parameters)
    return all(old_fingerprints.get(key) == fingerprints.get(key) for key in fingerprint_keys) and all(
        manifest.get("parameters", {}).get(key) == parameters.get(key) for key in parameter_keys
    )


def _read_signatures(path: Path) -> pd.DataFrame:
    signatures = pd.read_csv(path, index_col=0)
    signatures.index = signatures.index.astype(str)
    signatures.columns = signatures.columns.astype(str)
    if signatures.empty:
        raise ValueError("Reference signatures are empty")
    if not signatures.index.is_unique or not signatures.columns.is_unique:
        raise ValueError("Reference signature gene and cell-type names must be unique")
    values = signatures.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Reference signatures must be finite and non-negative")
    return signatures


def _extract_reference_signatures(adata_ref: ad.AnnData) -> pd.DataFrame:
    factor_names = [str(value) for value in adata_ref.uns["mod"]["factor_names"]]
    keys = [f"means_per_cluster_mu_fg_{factor}" for factor in factor_names]
    if "means_per_cluster_mu_fg" in adata_ref.varm:
        stored = adata_ref.varm["means_per_cluster_mu_fg"]
        if isinstance(stored, pd.DataFrame):
            signatures = stored.loc[:, keys].copy()
        else:
            signatures = pd.DataFrame(stored, index=adata_ref.var_names, columns=keys)
    elif all(key in adata_ref.var.columns for key in keys):
        signatures = adata_ref.var.loc[:, keys].copy()
    else:
        raise KeyError("Reference posterior did not contain means_per_cluster_mu_fg signatures")
    signatures.columns = factor_names
    signatures.index = adata_ref.var_names.astype(str)
    return signatures


def _prepare_output_directories(result_dir: Path) -> dict[str, Path]:
    directories = {
        "inputs": result_dir / "inputs",
        "reference": result_dir / "reference",
        "spatial": result_dir / "spatial",
        "tables": result_dir / "tables",
        "logs": result_dir / "logs",
    }
    for directory in directories.values():
        directory.mkdir(parents=True, exist_ok=True)
    return directories


def _clear_stage(directories: dict[str, Path], stage: str) -> None:
    targets = [directories[stage]]
    if stage == "spatial":
        targets.append(directories["tables"])
    for target in targets:
        if target.exists():
            shutil.rmtree(target)
        target.mkdir(parents=True, exist_ok=True)


def _reference_stage(
    config: dict[str, Any],
    directories: dict[str, Path],
    fingerprints: dict[str, Any],
    parameters: dict[str, Any],
) -> tuple[pd.DataFrame, bool]:
    from cell2location.models import RegressionModel
    from cell2location.utils.filtering import filter_genes

    signatures_out = directories["reference"] / "signatures.csv"
    stage_manifest_path = directories["reference"] / "manifest.json"
    done_path = directories["reference"] / ".complete"
    resume = bool(config.get("resume", True))
    overwrite = bool(config.get("overwrite", False))
    old_manifest = _read_json(stage_manifest_path) if stage_manifest_path.exists() else None
    matches = _manifest_matches(old_manifest, fingerprints, parameters, reference_only=True)

    if resume and done_path.exists() and signatures_out.exists() and matches:
        print("Reuse matching cell2location reference signatures")
        return _read_signatures(signatures_out), True
    if done_path.exists() and not matches:
        if not overwrite:
            raise RuntimeError(
                "Existing cell2location reference artifacts do not match the current input or parameters; "
                "set overwrite=TRUE to replace them"
            )
        _clear_stage(directories, "reference")

    signatures_path = config.get("signatures_path")
    if signatures_path:
        signatures = _read_signatures(Path(signatures_path))
        signatures.to_csv(signatures_out)
        shutil.copy2(signatures_path, directories["inputs"] / "reference_signatures.csv")
    else:
        reference_path = Path(config["reference_path"])
        adata_ref = ad.read_h5ad(reference_path)
        selected = filter_genes(adata_ref, **config.get("gene_filter_params", {}))
        adata_ref = adata_ref[:, np.asarray(selected)].copy()
        if adata_ref.n_vars == 0:
            raise ValueError("Reference gene filtering removed every gene")

        setup_kwargs: dict[str, Any] = {
            "adata": adata_ref,
            "labels_key": config["reference_label"],
        }
        if config.get("reference_batch") is not None:
            setup_kwargs["batch_key"] = config["reference_batch"]
        covariates = _normalise_optional_list(config.get("reference_covariates"))
        if covariates:
            setup_kwargs["categorical_covariate_keys"] = covariates
        RegressionModel.setup_anndata(**setup_kwargs)
        model_ref = RegressionModel(adata_ref)
        model_ref.train(**_normalise_train_params(config.get("reference_train_params")))
        adata_ref = model_ref.export_posterior(
            adata_ref,
            sample_kwargs=config.get("reference_posterior_params", {}),
        )
        signatures = _extract_reference_signatures(adata_ref)
        signatures.to_csv(signatures_out)
        model_ref.save(str(directories["reference"] / "model"), overwrite=True)
        adata_ref.write_h5ad(directories["reference"] / "reference_posterior.h5ad")
        shutil.copy2(reference_path, directories["inputs"] / "reference.h5ad")

    stage_manifest = {
        "fingerprints": fingerprints,
        "parameters": parameters,
        "versions": _versions(),
        "stage": "reference",
    }
    _write_json(stage_manifest, stage_manifest_path)
    done_path.write_text("complete\n", encoding="utf-8")
    return signatures, False


def _spatial_stage(
    config: dict[str, Any],
    directories: dict[str, Path],
    signatures: pd.DataFrame,
    fingerprints: dict[str, Any],
    parameters: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, bool]:
    from cell2location.models import Cell2location

    abundance_path = directories["tables"] / "abundance_q05.csv"
    proportions_path = directories["tables"] / "proportions.csv"
    done_path = directories["spatial"] / ".complete"
    manifest_path = Path(config["result_dir"]) / "manifest.json"
    resume = bool(config.get("resume", True))
    overwrite = bool(config.get("overwrite", False))
    old_manifest = _read_json(manifest_path) if manifest_path.exists() else None
    matches = _manifest_matches(old_manifest, fingerprints, parameters, reference_only=False)

    if resume and done_path.exists() and abundance_path.exists() and proportions_path.exists() and matches:
        print("Reuse matching cell2location spatial posterior")
        return _read_signatures(abundance_path), _read_signatures(proportions_path), True
    if done_path.exists() and not matches:
        if not overwrite:
            raise RuntimeError(
                "Existing cell2location spatial artifacts do not match the current input or parameters; "
                "set overwrite=TRUE to replace them"
            )
        _clear_stage(directories, "spatial")

    spatial_path = Path(config["spatial_path"])
    adata_sp = ad.read_h5ad(spatial_path)
    shared = adata_sp.var_names.intersection(signatures.index, sort=False)
    if len(shared) == 0:
        raise ValueError("Spatial data and reference signatures have no shared genes")
    adata_sp = adata_sp[:, shared].copy()
    signatures = signatures.loc[shared, :].copy()

    setup_kwargs: dict[str, Any] = {"adata": adata_sp}
    if config.get("spatial_batch") is not None:
        setup_kwargs["batch_key"] = config["spatial_batch"]
    Cell2location.setup_anndata(**setup_kwargs)
    model_sp = Cell2location(
        adata_sp,
        cell_state_df=signatures,
        N_cells_per_location=float(config["N_cells_per_location"]),
        detection_alpha=float(config["detection_alpha"]),
    )
    model_sp.train(**_normalise_train_params(config.get("spatial_train_params")))
    adata_sp = model_sp.export_posterior(
        adata_sp,
        use_quantiles=True,
        add_to_obsm=["q05"],
        sample_kwargs=config.get("spatial_posterior_params", {}),
    )

    factor_names = [str(value) for value in adata_sp.uns["mod"]["factor_names"]]
    raw_abundance = adata_sp.obsm["q05_cell_abundance_w_sf"]
    if isinstance(raw_abundance, pd.DataFrame):
        abundance = raw_abundance.copy()
        abundance.columns = factor_names
    else:
        abundance = pd.DataFrame(raw_abundance, index=adata_sp.obs_names, columns=factor_names)
    abundance.index = abundance.index.astype(str)
    denominator = abundance.sum(axis=1).clip(lower=1e-12)
    proportions = abundance.div(denominator, axis=0)

    abundance.to_csv(abundance_path)
    proportions.to_csv(proportions_path)
    model_sp.save(str(directories["spatial"] / "model"), overwrite=True)
    adata_sp.obsm["q05_cell_abundance_w_sf"] = abundance
    adata_sp.obsm["prop_celltypes"] = proportions
    adata_sp.write_h5ad(directories["spatial"] / "spatial_posterior.h5ad")
    shutil.copy2(spatial_path, directories["inputs"] / "spatial.h5ad")
    done_path.write_text("complete\n", encoding="utf-8")
    return abundance, proportions, False


def run(config: dict[str, Any]) -> None:
    result_dir = Path(config["result_dir"]).expanduser().resolve()
    result_dir.mkdir(parents=True, exist_ok=True)
    config["result_dir"] = str(result_dir)
    directories = _prepare_output_directories(result_dir)
    fingerprints = _fingerprints(config)
    parameters = _stable_parameters(config)
    parameters["runtime_versions"] = _versions()

    signatures, reused_reference = _reference_stage(
        config,
        directories,
        fingerprints,
        parameters,
    )
    abundance, proportions, reused_spatial = _spatial_stage(
        config,
        directories,
        signatures,
        fingerprints,
        parameters,
    )
    if list(abundance.index) != list(proportions.index) or list(abundance.columns) != list(proportions.columns):
        raise RuntimeError("Abundance and proportion output dimensions do not match")

    manifest = {
        "fingerprints": fingerprints,
        "parameters": parameters,
        "versions": _versions(),
        "reused_reference": reused_reference,
        "reused_spatial": reused_spatial,
        "n_spots": int(abundance.shape[0]),
        "n_cell_types": int(abundance.shape[1]),
        "n_signature_genes": int(signatures.shape[0]),
        "status": "complete",
    }
    _write_json(manifest, result_dir / "manifest.json")
    print(
        "cell2location completed: "
        f"{abundance.shape[0]} locations, {abundance.shape[1]} cell types, "
        f"{signatures.shape[0]} signature genes"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    run(_read_json(args.config))


if __name__ == "__main__":
    main()
