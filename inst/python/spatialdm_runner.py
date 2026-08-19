"""Subprocess runner for the pinned official SpatialDM backend.

The runner deliberately has a method-specific filename and exchanges only
plain files with R.  AnnData objects and the optional backend never enter the
Seurat object directly.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import spatialdm as sdm
from scipy import sparse
from scipy.io import mmread, mmwrite


def _read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError("SpatialDM configuration must be a JSON object")
    return value


def _write_json(value: dict, path: Path) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False, allow_nan=False)
    temporary.replace(path)


def _read_ids(path: Path, label: str) -> list[str]:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    if frame.shape[1] != 1:
        raise ValueError(f"{label} must contain exactly one column")
    values = frame.iloc[:, 0].astype(str).tolist()
    if not values or len(values) != len(set(values)) or any(not value for value in values):
        raise ValueError(f"{label} must contain unique non-empty IDs")
    return values


def _build_adata(config: dict) -> anndata.AnnData:
    cells = _read_ids(Path(config["cells_path"]), "cell table")
    features = _read_ids(Path(config["features_path"]), "feature table")
    coordinates = pd.read_csv(config["coordinates_path"], sep="\t", dtype={"cell_id": str})
    if list(coordinates.columns) != ["cell_id", "x", "y"]:
        raise ValueError("coordinate columns must be cell_id, x and y")
    if coordinates["cell_id"].astype(str).tolist() != cells:
        raise ValueError("coordinate IDs are not identically ordered with cell IDs")
    xy = coordinates[["x", "y"]].to_numpy(dtype=float)
    if not np.isfinite(xy).all():
        raise ValueError("raw coordinates must be finite")

    expression = sparse.csr_matrix(mmread(config["expression_path"])).transpose().tocsr()
    counts = sparse.csr_matrix(mmread(config["counts_path"])).transpose().tocsr()
    expected = (len(cells), len(features))
    if expression.shape != expected or counts.shape != expected:
        raise ValueError("expression dimensions do not match cell and feature IDs")
    for matrix, label in ((expression, "normalized expression"), (counts, "count expression")):
        if matrix.data.size and (not np.isfinite(matrix.data).all() or (matrix.data < 0).any()):
            raise ValueError(f"{label} must contain finite non-negative values")

    obs = pd.DataFrame(index=pd.Index(cells, name="cell_id"))
    var = pd.DataFrame(index=pd.Index(features, name="feature"))
    adata = anndata.AnnData(X=expression, obs=obs.copy(), var=var.copy())
    raw = anndata.AnnData(X=counts, obs=obs.copy(), var=var.copy())
    adata.raw = raw
    adata.obsm["spatial"] = xy
    return adata


def _custom_lr(adata: anndata.AnnData, path: Path, min_cell: int) -> None:
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"ligand", "receptor", "annotation"}
    if not required.issubset(table.columns) or table.empty:
        raise ValueError("custom LR table requires ligand, receptor and annotation columns")
    genes = set(adata.var_names.astype(str))

    def components(value: str) -> list[str]:
        return [item.strip() for item in str(value).replace(",", ";").split(";") if item.strip()]

    rows = []
    ligands = []
    receptors = []
    for idx, row in table.iterrows():
        ligand = [gene for gene in components(row["ligand"]) if gene in genes]
        receptor = [gene for gene in components(row["receptor"]) if gene in genes]
        if not ligand or not receptor:
            continue
        if any(np.asarray((adata[:, gene].X > 0).sum()).item() < min_cell for gene in ligand + receptor):
            continue
        interaction = str(row.get("interaction", "")) or "_".join(ligand + receptor)
        rows.append({
            "interaction_name": interaction,
            "annotation": str(row["annotation"]),
            "pathway_name": str(row.get("pathway_name", "")),
        })
        ligands.append(ligand)
        receptors.append(receptor)
    if not rows:
        raise ValueError("No custom LR pairs remain after feature and expression filtering")
    index = pd.Index([row["interaction_name"] for row in rows], name="interaction")
    if index.has_duplicates:
        raise ValueError("custom LR interaction names must be unique")
    max_l = max(map(len, ligands))
    max_r = max(map(len, receptors))
    ldf = pd.DataFrame([[genes[i] if i < len(genes) else None for i in range(max_l)] for genes in ligands], index=index)
    rdf = pd.DataFrame([[genes[i] if i < len(genes) else None for i in range(max_r)] for genes in receptors], index=index)
    ldf.columns = [f"Ligand{i}" for i in range(max_l)]
    rdf.columns = [f"Receptor{i}" for i in range(max_r)]
    adata.uns["mean"] = str("algebra")
    adata.uns["ligand"] = ldf
    adata.uns["receptor"] = rdf
    adata.uns["geneInter"] = pd.DataFrame(rows, index=index)
    adata.uns["num_pairs"] = len(rows)


def _write_matrix(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=True, na_rep="NaN")


def _write_empty_matrix(path: Path, rows: list[str], columns: list[str]) -> None:
    _write_matrix(pd.DataFrame(index=pd.Index(rows, name="interaction"), columns=columns, dtype=float), path)


def _weight_summary(matrix: sparse.spmatrix, xy: np.ndarray) -> dict:
    coo = sparse.coo_matrix(matrix)
    if coo.nnz:
        distance = np.sqrt(np.square(xy[coo.row] - xy[coo.col]).sum(axis=1))
        finite = distance[np.isfinite(distance)]
        max_distance = float(finite.max()) if finite.size else None
        median_degree = float(np.median(np.bincount(coo.row, minlength=xy.shape[0])))
    else:
        max_distance = None
        median_degree = 0.0
    return {"nnz": int(coo.nnz), "median_degree": median_degree, "effective_range": max_distance}


def run(config: dict) -> None:
    np.random.seed(int(config.get("seed", 1)))
    output_dir = Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    adata = _build_adata(config)

    l = config.get("l")
    eff_dist = config.get("eff_dist")
    if l is None and eff_dist is None:
        raise ValueError("At least one of l and eff_dist must be supplied")
    n_neighbors = config.get("n_neighbors")
    if n_neighbors is None:
        n_neighbors = min(6 * 31, max(2, adata.n_obs))
    n_neighbors = int(n_neighbors)
    if n_neighbors < 2 or n_neighbors > adata.n_obs:
        raise ValueError("n_neighbors must be between 2 and the number of spots")
    n_neighbor_layers = min(int(config.get("n_neighbor_layers", 6)), adata.n_obs)
    if n_neighbor_layers < 1:
        raise ValueError("n_neighbor_layers must be positive")

    sdm.weight_matrix(
        adata,
        l=None if l is None else float(l),
        cutoff=float(config.get("cutoff", 0.1)),
        n_neighbors=n_neighbors,
        n_nearest_neighbors=n_neighbor_layers,
        single_cell=bool(config.get("single_cell", False)),
        eff_dist=None if eff_dist is None else float(eff_dist),
    )
    if "weight" not in adata.obsp or "nearest_neighbors" not in adata.obsp:
        raise ValueError("SpatialDM did not produce both secreted and contact weight matrices")
    spatial_w = adata.obsp["weight"]
    nearest = adata.obsp["nearest_neighbors"]
    adata.obsp["weight"] = spatial_w
    adata.obsp["nearest_neighbors"] = nearest

    custom_path = config.get("custom_lr_path")
    if custom_path:
        _custom_lr(adata, Path(custom_path), int(config.get("min_cell", 3)))
    else:
        sdm.extract_lr(
            adata,
            str(config["species"]),
            mean=str(config.get("complex_mean", "algebra")),
            min_cell=int(config.get("min_cell", 3)),
            datahost="package",
        )

    method = str(config.get("method", "z-score"))
    n_perm = int(config.get("n_perm", 1000))
    nproc = int(config.get("nproc", 1))
    sdm.spatialdm_global(adata, n_perm=n_perm, method=method, nproc=nproc)
    sdm.sig_pairs(
        adata,
        method=method,
        fdr=bool(config.get("global_fdr", True)),
        threshold=float(config.get("global_threshold", 0.1)),
    )

    global_result = adata.uns["global_res"].copy()
    gene_inter = pd.DataFrame(adata.uns["geneInter"])
    for column in ("annotation", "pathway_name"):
        if column in gene_inter.columns and column not in global_result.columns:
            global_result[column] = gene_inter.loc[global_result.index, column].to_numpy()
    global_result.insert(0, "interaction", global_result.index.astype(str))
    global_result["moran_r"] = np.asarray(adata.uns["global_I"], dtype=float)
    p_col = "z_pval" if method == "z-score" else "perm_pval"
    if p_col in global_result.columns:
        global_result["pvalue"] = pd.to_numeric(global_result[p_col], errors="coerce")
    else:
        raise ValueError(f"SpatialDM global output is missing {p_col}")
    global_result["z_score"] = pd.to_numeric(global_result.get("z", np.nan), errors="coerce")
    global_result["score_type"] = "spatial_association"
    global_result["n_spots"] = 0

    selected = global_result.loc[global_result["selected"].astype(bool), "interaction"].tolist()
    if bool(config.get("run_local", True)) and selected:
        sdm.spatialdm_local(adata, n_perm=n_perm, method=method, nproc=nproc)
        sdm.sig_spots(
            adata,
            method=method,
            fdr=bool(config.get("local_fdr", False)),
            threshold=float(config.get("local_threshold", 0.1)),
        )
        local_p = adata.uns["local_z_p"] if method == "z-score" else adata.uns["local_perm_p"]
        local_z = adata.uns.get("local_z")
        local_i = np.asarray(adata.uns["local_stat"]["local_I"], dtype=float).T
        local_i_r = np.asarray(adata.uns["local_stat"]["local_I_R"], dtype=float).T
        _write_matrix(pd.DataFrame(local_i, index=selected, columns=adata.obs_names), output_dir / "local_i.tsv")
        _write_matrix(pd.DataFrame(local_i_r, index=selected, columns=adata.obs_names), output_dir / "local_i_r.tsv")
        _write_matrix(pd.DataFrame(np.asarray(local_z, dtype=float), index=selected, columns=adata.obs_names), output_dir / "local_z.tsv")
        _write_matrix(pd.DataFrame(np.asarray(local_p, dtype=float), index=selected, columns=adata.obs_names), output_dir / "local_p.tsv")
        selected_spots = pd.DataFrame(np.asarray(adata.uns["selected_spots"], dtype=int), index=selected, columns=adata.obs_names)
        _write_matrix(selected_spots, output_dir / "selected_spots.tsv")
        n_spots = selected_spots.sum(axis=1).rename("n_spots")
        n_spots.to_csv(output_dir / "n_spots.tsv", sep="\t", header=True)
        global_result.loc[global_result["interaction"].isin(n_spots.index), "n_spots"] = global_result.loc[
            global_result["interaction"].isin(n_spots.index), "interaction"
        ].map(n_spots).fillna(0).astype(int)
    else:
        _write_empty_matrix(output_dir / "local_i.tsv", selected, list(adata.obs_names))
        _write_empty_matrix(output_dir / "local_i_r.tsv", selected, list(adata.obs_names))
        _write_empty_matrix(output_dir / "local_z.tsv", selected, list(adata.obs_names))
        _write_empty_matrix(output_dir / "local_p.tsv", selected, list(adata.obs_names))
        _write_empty_matrix(output_dir / "selected_spots.tsv", selected, list(adata.obs_names))
        pd.DataFrame({"n_spots": []}, index=pd.Index([], name="interaction")).to_csv(output_dir / "n_spots.tsv", sep="\t")

    global_result.to_csv(output_dir / "global.csv", index=False)

    mmwrite(output_dir / "weight.mtx", sparse.coo_matrix(spatial_w))
    mmwrite(output_dir / "nearest_neighbors.mtx", sparse.coo_matrix(nearest))
    pd.DataFrame({"cell_id": adata.obs_names.astype(str), "x": adata.obsm["spatial"][:, 0], "y": adata.obsm["spatial"][:, 1]}).to_csv(output_dir / "coordinates.tsv", sep="\t", index=False)
    versions = {}
    for package in ("SpatialDM", "anndata", "numpy", "pandas", "scipy", "scanpy", "SparseAEH"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    manifest = {
        "schema_version": 1,
        "backend": "SpatialDM",
        "species": str(config.get("species")),
        "method": method,
        "database_name": "custom" if custom_path else "CellChatDB",
        "database_source": "input table" if custom_path else "SpatialDM package",
        "cell_ids": adata.obs_names.astype(str).tolist(),
        "feature_ids": adata.var_names.astype(str).tolist(),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "n_candidate_pairs": int(len(adata.uns["geneInter"])),
        "n_selected_pairs": int(len(selected)),
        "n_local_spots": int(sum(pd.read_csv(output_dir / "n_spots.tsv", sep="\t").iloc[:, -1])) if selected else 0,
        "coordinate_space": "raw",
        "parameters": {key: config.get(key) for key in ("l", "eff_dist", "cutoff", "n_neighbors", "n_neighbor_layers", "single_cell", "min_cell", "global_fdr", "global_threshold", "local_fdr", "local_threshold", "seed")},
        "weight_summary": {"secreted": _weight_summary(spatial_w, adata.obsm["spatial"]), "contact": _weight_summary(nearest, adata.obsm["spatial"])},
        "versions": versions,
    }
    _write_json(manifest, output_dir / "manifest.json")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    run(_read_json(Path(args.config)))


if __name__ == "__main__":
    main()
