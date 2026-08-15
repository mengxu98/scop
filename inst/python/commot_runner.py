"""Subprocess runner for the official COMMOT Python package.

The file deliberately does not use the package name ``commot.py`` because that
would shadow the installed official package when Python resolves imports.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmread


COMMUNICATION_COLUMNS = [
    "sender", "receiver", "ligand", "receptor", "interaction_name",
    "pathway_name", "score", "pvalue", "method", "score_type",
]
CLUSTER_COLUMNS = COMMUNICATION_COLUMNS + ["key"]
DIRECTION_COLUMNS = [
    "cell_id", "group", "key", "perspective", "x", "y", "dx", "dy",
    "magnitude",
]


def _read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError("COMMOT configuration must be a JSON object")
    return value


def _write_json(value: dict, path: Path) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False, allow_nan=False)
    temporary.replace(path)


def _read_single_column(path: Path, label: str) -> list[str]:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    if frame.shape[1] != 1:
        raise ValueError(f"{label} must contain exactly one column")
    values = frame.iloc[:, 0].astype(str).tolist()
    if not values or len(values) != len(set(values)) or any(not value for value in values):
        raise ValueError(f"{label} must contain unique non-empty IDs")
    return values


def _sha256_ids(values: list[str]) -> str:
    return hashlib.sha256("\n".join(values).encode("utf-8")).hexdigest()


def _json_args(value: object, label: str) -> dict:
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object")
    return dict(value)


def _normalize_api_args(value: dict) -> dict:
    value = dict(value)
    for key in ("lr_pair", "cot_weights", "pos_idx"):
        if key in value and isinstance(value[key], list):
            if key == "pos_idx":
                value[key] = np.asarray(value[key], dtype=int)
            else:
                value[key] = tuple(value[key])
    return value


def _build_adata(config: dict) -> anndata.AnnData:
    features = _read_single_column(Path(config["features_path"]), "feature table")
    metadata = pd.read_csv(config["metadata_path"], sep="\t", dtype=str)
    coordinates = pd.read_csv(config["coordinates_path"], sep="\t")
    if list(metadata.columns) != ["cell_id", "group"]:
        raise ValueError("metadata columns must be cell_id and group")
    if list(coordinates.columns) != ["cell_id", "x", "y"]:
        raise ValueError("coordinate columns must be cell_id, x, and y")
    cell_ids = metadata["cell_id"].astype(str).tolist()
    coordinate_ids = coordinates["cell_id"].astype(str).tolist()
    if not cell_ids or len(cell_ids) != len(set(cell_ids)):
        raise ValueError("metadata cell IDs must be unique")
    if cell_ids != coordinate_ids:
        raise ValueError("metadata and coordinate IDs are not identically ordered")
    groups = metadata["group"]
    if groups.isna().any() or (groups.astype(str).str.len() == 0).any():
        raise ValueError("group labels must be non-missing and non-empty")
    groups = groups.astype(str)
    xy = coordinates[["x", "y"]].to_numpy(dtype=float)
    if not np.isfinite(xy).all():
        raise ValueError("raw coordinates must be finite")

    matrix = mmread(config["matrix_path"])
    matrix = sparse.csr_matrix(matrix).transpose().tocsr()
    if matrix.shape != (len(cell_ids), len(features)):
        raise ValueError("MatrixMarket dimensions do not match metadata and features")
    if matrix.data.size and (not np.isfinite(matrix.data).all() or (matrix.data < 0).any()):
        raise ValueError("expression matrix must contain finite non-negative values")
    obs = pd.DataFrame(index=pd.Index(cell_ids, name="cell_id"))
    obs["_scop_group"] = groups.to_numpy()
    var = pd.DataFrame(index=pd.Index(features, name="feature"))
    adata = anndata.AnnData(X=matrix, obs=obs, var=var)
    adata.obsm["spatial"] = xy
    return adata


def _aggregate_transport(adata: anndata.AnnData, database_name: str) -> pd.DataFrame:
    info_key = f"commot-{database_name}-info"
    info = adata.uns.get(info_key, {})
    ligrec = info.get("df_ligrec")
    if ligrec is None:
        raise ValueError("COMMOT did not preserve its ligand-receptor database")
    ligrec = pd.DataFrame(ligrec).iloc[:, :3].copy()
    ligrec.columns = ["ligand", "receptor", "pathway_name"]
    ligrec = ligrec.drop_duplicates(subset=["ligand", "receptor"])
    groups = adata.obs["_scop_group"].astype(str).to_numpy()
    rows: list[pd.DataFrame] = []
    for item in ligrec.itertuples(index=False):
        key = f"commot-{database_name}-{item.ligand}-{item.receptor}"
        if key not in adata.obsp:
            continue
        matrix = sparse.coo_matrix(adata.obsp[key])
        if matrix.nnz == 0:
            continue
        edge = pd.DataFrame({
            "sender": groups[matrix.row],
            "receiver": groups[matrix.col],
            "score": matrix.data.astype(float),
        })
        edge = edge.groupby(["sender", "receiver"], as_index=False, sort=False)["score"].sum()
        edge["ligand"] = str(item.ligand)
        edge["receptor"] = str(item.receptor)
        edge["interaction_name"] = f"{item.ligand}_{item.receptor}"
        edge["pathway_name"] = str(item.pathway_name)
        edge["pvalue"] = np.nan
        edge["method"] = "COMMOT"
        edge["score_type"] = "transport_mass"
        rows.append(edge[COMMUNICATION_COLUMNS])
    if not rows:
        return pd.DataFrame(columns=COMMUNICATION_COLUMNS)
    result = pd.concat(rows, ignore_index=True)
    return result[np.isfinite(result["score"]) & (result["score"] > 0)].reset_index(drop=True)


def _selection_key(database_name: str, args: dict) -> tuple:
    pathway = args.get("pathway_name")
    lr_pair = args.get("lr_pair")
    if pathway is not None and lr_pair is not None:
        raise ValueError("pathway_name and lr_pair cannot both be selected")
    if pathway is not None:
        return f"{database_name}-{pathway}", None, None, str(pathway)
    if lr_pair is not None:
        if not isinstance(lr_pair, (list, tuple)) or len(lr_pair) != 2:
            raise ValueError("lr_pair must contain ligand and receptor")
        return f"{database_name}-{lr_pair[0]}-{lr_pair[1]}", str(lr_pair[0]), str(lr_pair[1]), None
    return f"{database_name}-total-total", None, None, None


def _cluster_table(adata: anndata.AnnData, database_name: str, args: dict) -> pd.DataFrame:
    selection, ligand, receptor, pathway = _selection_key(database_name, args)
    uns_key = f"commot_cluster-_scop_group-{selection}"
    if uns_key not in adata.uns:
        raise ValueError(f"COMMOT cluster output is missing {uns_key}")
    value = adata.uns[uns_key]
    score = pd.DataFrame(value["communication_matrix"])
    pvalue = pd.DataFrame(value["communication_pvalue"])
    if (
        score.shape != pvalue.shape
        or list(score.index.astype(str)) != list(pvalue.index.astype(str))
        or list(score.columns.astype(str)) != list(pvalue.columns.astype(str))
    ):
        raise ValueError("COMMOT cluster score and p-value matrices are misaligned")
    rows = []
    for sender in score.index:
        for receiver in score.columns:
            rows.append({
                "sender": str(sender), "receiver": str(receiver),
                "ligand": ligand, "receptor": receptor,
                "interaction_name": None if ligand is None else f"{ligand}_{receptor}",
                "pathway_name": pathway,
                "score": float(score.loc[sender, receiver]),
                "pvalue": float(pvalue.loc[sender, receiver]),
                "method": "COMMOT", "score_type": "cluster_communication",
                "key": selection,
            })
    return pd.DataFrame(rows, columns=CLUSTER_COLUMNS)


def _direction_table(adata: anndata.AnnData, database_name: str, args: dict) -> pd.DataFrame:
    selection, _, _, _ = _selection_key(database_name, args)
    xy = np.asarray(adata.obsm["spatial"], dtype=float)
    frames = []
    for perspective in ("sender", "receiver"):
        obsm_key = f"commot_{perspective}_vf-{selection}"
        if obsm_key not in adata.obsm:
            raise ValueError(f"COMMOT direction output is missing {obsm_key}")
        vector = np.asarray(adata.obsm[obsm_key], dtype=float)
        if vector.shape[0] != adata.n_obs or vector.shape[1] < 2:
            raise ValueError("COMMOT direction vectors have invalid dimensions")
        frames.append(pd.DataFrame({
            "cell_id": adata.obs_names.astype(str),
            "group": adata.obs["_scop_group"].astype(str).to_numpy(),
            "key": selection,
            "perspective": perspective,
            "x": xy[:, 0], "y": xy[:, 1],
            "dx": vector[:, 0], "dy": vector[:, 1],
            "magnitude": np.sqrt(np.square(vector[:, 0]) + np.square(vector[:, 1])),
        }))
    return pd.concat(frames, ignore_index=True)[DIRECTION_COLUMNS]


def run(config: dict) -> None:
    import commot as ct

    output_dir = Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    adata = _build_adata(config)
    database_name = str(config["database"])
    species = str(config["species"])
    communication_args = _normalize_api_args(_json_args(config.get("communication_args"), "communication_args"))
    cluster_args = _normalize_api_args(_json_args(config.get("cluster_args"), "cluster_args"))
    direction_args = _normalize_api_args(_json_args(config.get("direction_args"), "direction_args"))

    normalize = bool(communication_args.pop("normalize", True))
    target_sum = float(communication_args.pop("target_sum", 1e4))
    signaling_type = communication_args.pop("signaling_type", "Secreted Signaling")
    filter_args = _json_args(communication_args.pop("filter_args", {}), "communication_args.filter_args")
    if normalize:
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)

    df_ligrec = ct.pp.ligand_receptor_database(
        database=database_name,
        species=species,
        signaling_type=signaling_type,
    )
    filter_defaults = {"heteromeric": True, "filter_criteria": "min_cell_pct", "min_cell_pct": 0.05}
    filter_defaults.update(filter_args)
    df_ligrec = ct.pp.filter_lr_database(df_ligrec, adata, **filter_defaults)
    if df_ligrec.shape[0] == 0:
        raise ValueError("No ligand-receptor pairs remain after official COMMOT filtering")

    reserved = {"adata", "database_name", "df_ligrec", "dis_thr", "copy"}
    overlap = reserved.intersection(communication_args)
    if overlap:
        raise ValueError(f"communication_args cannot override reserved fields: {sorted(overlap)}")
    spatial_args = {"pathway_sum": True, "heteromeric": True}
    spatial_args.update(communication_args)
    ct.tl.spatial_communication(
        adata,
        database_name=database_name,
        df_ligrec=df_ligrec,
        dis_thr=config.get("distance_threshold"),
        **spatial_args,
    )
    communication = _aggregate_transport(adata, database_name)
    if communication.empty:
        raise ValueError("Official COMMOT returned no positive transport mass")

    cluster = pd.DataFrame(columns=CLUSTER_COLUMNS)
    if bool(config.get("cluster", False)):
        if set(cluster_args).intersection({"adata", "database_name", "clustering", "copy"}):
            raise ValueError("cluster_args contains reserved fields")
        ct.tl.cluster_communication(
            adata, database_name=database_name, clustering="_scop_group", **cluster_args
        )
        cluster = _cluster_table(adata, database_name, cluster_args)

    direction = pd.DataFrame(columns=DIRECTION_COLUMNS)
    if bool(config.get("direction", False)):
        if set(direction_args).intersection({"adata", "database_name", "copy"}):
            raise ValueError("direction_args contains reserved fields")
        ct.tl.communication_direction(adata, database_name=database_name, **direction_args)
        direction = _direction_table(adata, database_name, direction_args)

    communication_path = output_dir / "communication.csv"
    cluster_path = output_dir / "cluster.csv"
    direction_path = output_dir / "direction.csv"
    communication.to_csv(communication_path, index=False)
    cluster.to_csv(cluster_path, index=False)
    direction.to_csv(direction_path, index=False)

    h5ad_path = None
    if bool(config.get("store_h5ad", False)):
        h5ad_path = output_dir / "commot_result.h5ad"
        adata.write_h5ad(h5ad_path)

    cell_ids = adata.obs_names.astype(str).tolist()
    feature_ids = adata.var_names.astype(str).tolist()
    versions = {
        package: importlib.metadata.version(package)
        for package in ("commot", "anndata", "numpy", "pandas", "scanpy", "scipy")
    }
    manifest = {
        "schema_version": 1,
        "backend": "commot",
        "database": database_name,
        "species": species,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "cell_ids": cell_ids,
        "feature_ids": feature_ids,
        "feature_ids_sha256": _sha256_ids(feature_ids),
        "cell_ids_sha256": _sha256_ids(cell_ids),
        "n_communication_rows": int(communication.shape[0]),
        "n_cluster_rows": int(cluster.shape[0]),
        "n_direction_rows": int(direction.shape[0]),
        "coordinate_space": "raw",
        "distance_threshold": config.get("distance_threshold"),
        "normalized": normalize,
        "target_sum": target_sum if normalize else None,
        "versions": versions,
        "files": {
            "communication": communication_path.name,
            "cluster": cluster_path.name,
            "direction": direction_path.name,
            "h5ad": None if h5ad_path is None else h5ad_path.name,
        },
    }
    _write_json(manifest, output_dir / "manifest.json")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    run(_read_json(Path(args.config)))


if __name__ == "__main__":
    main()
