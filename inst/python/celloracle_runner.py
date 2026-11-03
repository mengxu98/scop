"""Object adaptation for the unmodified CellOracle 0.20.0 workflow."""

from pathlib import Path
import inspect
import os
import random

import numpy as np
import pandas as pd
from scipy import sparse


def run_celloracle(adata, group_by, embedding_name, base_grn, perturbation,
                   stage_args, cores=1, seed=1234, model_dir=None,
                   output_dir=None, store_expression=False):
    import celloracle as co

    if co.__version__ != "0.20.0":
        raise ValueError("This adapter requires official celloracle==0.20.0.")
    stage_args = dict(stage_args)
    # Restore caller RNG state; do not replace any upstream numerical routine.
    numpy_state, python_state = np.random.get_state(), random.getstate()
    try:
        np.random.seed(seed)
        random.seed(seed)
        oracle = co.Oracle()
        calls = {}

        def call(stage, method, defaults=None):
            args = dict(defaults or {})
            args.update(stage_args.get(stage, {}) or {})
            inspect.signature(method).bind(**args)
            calls[stage] = args.copy()
            return method(**args)

        if model_dir is None:
            oracle.import_anndata_as_raw_count(
                adata=adata, cluster_column_name=group_by,
                embedding_name=embedding_name)
            if isinstance(base_grn, str):
                oracle.import_TF_data(TF_info_matrix_path=base_grn)
            elif isinstance(base_grn, pd.DataFrame):
                oracle.import_TF_data(TF_info_matrix=base_grn)
            else:
                oracle.import_TF_data(TFdict={
                    target: list(np.atleast_1d(tfs))
                    for target, tfs in base_grn.items()})
            call("pca", oracle.perform_PCA)
            call("imputation", oracle.knn_imputation, {"n_jobs": cores})
            links = call("links", oracle.get_links, {"n_jobs": cores})
            call("filter", links.filter_links)
            oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
            call("fit", oracle.fit_GRN_for_simulation,
                 {"use_cluster_specific_TFdict": True})
        else:
            oracle = co.load_hdf5(str(Path(model_dir) / "model.celloracle.oracle"))
            links = co.load_hdf5(str(Path(model_dir) / "network.celloracle.links"))
            if oracle.celloracle_version_used != co.__version__:
                raise ValueError("Saved model was fitted with another CellOracle version.")
            if (list(oracle.adata.obs_names) != list(adata.obs_names)
                    or list(oracle.adata.var_names) != list(adata.var_names)
                    or oracle.cluster_column_name != group_by
                    or oracle.embedding_name != embedding_name
                    or not np.array_equal(oracle.adata.obs[group_by].astype(str),
                                          adata.obs[group_by].astype(str))
                    or not np.array_equal(oracle.embedding, adata.obsm[embedding_name])):
                raise ValueError("Saved model cells, genes, groups or embedding differ from input.")
            difference = sparse.csr_matrix(oracle.adata.layers["raw_count"]) - sparse.csr_matrix(adata.X)
            if difference.count_nonzero():
                raise ValueError("Saved model raw counts differ from input.")

        call("simulation", oracle.simulate_shift,
             {"perturb_condition": perturbation})
        call("transition", oracle.estimate_transition_prob,
             {"n_jobs": cores, "threads": cores, "random_seed": seed})
        call("embedding", oracle.calculate_embedding_shift)
        call("grid", oracle.calculate_p_mass)
        call("mass_filter", oracle.calculate_mass_filter, {"plot": False})

        coordinates = np.asarray(oracle.embedding)
        shift = np.asarray(oracle.delta_embedding)
        shift_random = np.asarray(oracle.delta_embedding_random)
        cell_shift = pd.DataFrame({
            "cell": list(oracle.adata.obs_names),
            "group": oracle.adata.obs[group_by].astype(str).to_numpy(),
            "x": coordinates[:, 0], "y": coordinates[:, 1],
            "dx": shift[:, 0], "dy": shift[:, 1],
            "dx_random": shift_random[:, 0], "dy_random": shift_random[:, 1]})
        grid = pd.DataFrame({
            "x": oracle.flow_grid[:, 0], "y": oracle.flow_grid[:, 1],
            "dx": oracle.flow[:, 0], "dy": oracle.flow[:, 1],
            "dx_random": oracle.flow_rndm[:, 0], "dy_random": oracle.flow_rndm[:, 1],
            "mass": oracle.total_p_mass, "filtered": oracle.mass_filter})
        delta = np.asarray(oracle.adata.layers["delta_X"])
        effect = pd.DataFrame({"gene": list(oracle.adata.var_names),
                               "mean_delta": delta.mean(axis=0),
                               "mean_abs_delta": np.abs(delta).mean(axis=0)})
        result = {"cell_shift": cell_shift, "grid_flow": grid,
                  "gene_effect": effect,
                  "links": {str(k): v.copy() for k, v in links.links_dict.items()},
                  "filtered_links": {str(k): v.copy() for k, v in links.filtered_links.items()},
                  "parameters": {"celloracle_version": co.__version__,
                                 "python_hash_seed": os.environ.get("PYTHONHASHSEED"),
                                 "seed": seed, "cores": cores, "calls": calls,
                                 "model_dir": model_dir},
                  "cells": list(oracle.adata.obs_names),
                  "genes": list(oracle.adata.var_names)}
        if store_expression:
            result["expression"] = {key: np.asarray(oracle.adata.layers[key]).copy()
                                    for key in ("imputed_count", "simulated_count", "delta_X")}
        if output_dir is not None:
            destination = Path(output_dir)
            destination.mkdir(parents=True, exist_ok=False)
            oracle.to_hdf5(str(destination / "model.celloracle.oracle"))
            links.to_hdf5(str(destination / "network.celloracle.links"))
            result["model_dir"] = str(destination.resolve())
        return result
    finally:
        np.random.set_state(numpy_state)
        random.setstate(python_state)


def main():
    """File transport for R sessions using an already prepared Linux container."""
    import json
    import sys
    import anndata as ad
    from scipy.io import mmread

    directory = Path(sys.argv[1])
    request = json.loads((directory / "request.json").read_text())
    adata = ad.AnnData(mmread(directory / "counts.mtx").T.tocsr())
    adata.obs_names = request["cells"]
    adata.var_names = request["genes"]
    adata.obs[request["group_by"]] = pd.Categorical(
        request["groups"], categories=np.atleast_1d(request["levels"]).tolist(),
        ordered=request["ordered"])
    adata.obsm[request["embedding_name"]] = np.loadtxt(
        directory / "embedding.csv", delimiter=",", ndmin=2)
    prior_type = request["prior_type"]
    if prior_type == "matrix":
        prior = pd.read_csv(directory / "prior.csv", dtype={"peak_id": str, "gene_short_name": str})
    elif prior_type == "parquet":
        prior = str(directory / "prior.parquet")
    else:
        prior = request["prior"]
    output = directory / "output"
    output.mkdir()
    result = run_celloracle(
        adata, request["group_by"], request["embedding_name"], prior,
        request["perturbation"], request["stage_args"],
        cores=request["cores"], seed=request["seed"],
        model_dir=str(directory / "model") if request["reuse_model"] else None,
        output_dir=str(output / "model") if request["save_model"] else None,
        store_expression=request["store_expression"])
    manifest = {key: result[key] for key in ("cells", "genes", "parameters")}
    for name in ("cell_shift", "grid_flow", "gene_effect"):
        result[name].to_csv(output / (name + ".csv"), index=False)
    for name in ("links", "filtered_links"):
        manifest[name] = {}
        for index, (cluster, frame) in enumerate(result[name].items()):
            filename = f"{name}_{index}.csv"
            frame.to_csv(output / filename, index=False)
            manifest[name][cluster] = filename
    manifest["expression"] = {}
    for name, matrix in result.get("expression", {}).items():
        filename = name + ".csv"
        np.savetxt(output / filename, matrix, delimiter=",", fmt="%.17g")
        manifest["expression"][name] = filename
    (output / "manifest.json").write_text(json.dumps(manifest))


if __name__ == "__main__":
    main()
