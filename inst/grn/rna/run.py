import importlib.metadata
import inspect
import json
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def infer(expression, request):
    method = request["method"]
    parameters = request["method_args"] or {}
    seed = int(request["seed"])
    np.random.seed(seed)
    random.seed(seed)
    regulators = request["regulators"]
    if method == "grnboost2":
        from arboreto.algo import grnboost2
        from distributed import Client, LocalCluster

        with LocalCluster(n_workers=int(request["cores"]), threads_per_worker=1,
                          processes=False, dashboard_address=None) as cluster:
            with Client(cluster) as client:
                edges = grnboost2(
                    expression_data=expression.T, tf_names=regulators,
                    client_or_address=client, seed=seed,
                    early_stop_window_length=int(parameters.get("early_stop_window_length", 25)),
                    verbose=bool(request["verbose"]),
                )
        parameters = dict(early_stop_window_length=int(parameters.get("early_stop_window_length", 25)))
        return edges.rename(columns={"TF": "regulator", "importance": "weight"}), "arboreto::grnboost2", ["arboreto", "dask", "distributed"], parameters
    if method == "gniplr":
        from gniplr import infer_gniplr

        parameters = dict(lasso_degree=int(parameters.get("lasso_degree", 30)), lasso_alpha=float(parameters.get("lasso_alpha", 0.04)))
        if expression.shape[1] < 8 or parameters["lasso_degree"] < 1 or parameters["lasso_alpha"] < 0:
            raise ValueError("GNIPLR requires at least eight cells, positive lasso_degree and nonnegative lasso_alpha")
        output = Path(request["work_dir"]) / "gniplr-native.csv"
        infer_gniplr(expression, output,
                     lasso_degree=int(parameters.get("lasso_degree", 30)),
                     lasso_alpha=float(parameters.get("lasso_alpha", 0.04)))
        return pd.read_csv(output, dtype={"regulator": str, "target": str}, keep_default_na=False), "multiCSN-manuscript::GNIPLR", ["statsmodels", "scipy", "scikit-learn"], parameters
    if method == "scsgl":
        sys.path.insert(0, str(Path(__file__).parent / "vendor" / "SCSGL"))
        import run_scSGL

        parameters = dict(pos_density=float(parameters.get("pos_density", 0.01)), neg_density=float(parameters.get("neg_density", 0.01)))
        if any(not np.isfinite(value) or value < 0 or value >= 1 for value in parameters.values()):
            raise ValueError("SCSGL densities must be finite, nonnegative and smaller than one")
        if not any(value > 0 for value in parameters.values()):
            raise ValueError("SCSGL requires at least one positive edge density")
        output = Path(request["work_dir"]) / "scsgl-native.tsv"
        run_scSGL.run_inference(
            expression, output,
            pos_density=float(parameters.get("pos_density", 0.01)),
            neg_density=float(parameters.get("neg_density", 0.01)),
        )
        edges = pd.read_csv(output, sep="\t", dtype={"Gene1": str, "Gene2": str}, keep_default_na=False)
        edges = edges.rename(columns={"Gene1": "regulator", "Gene2": "target", "EdgeWeight": "weight"})
        return edges.sort_values("weight", ascending=False), "scSGL::signed_graph_learning", ["scipy", "numba"], parameters
    if method == "regdiffusion":
        import torch
        import regdiffusion as rd

        torch.set_num_threads(int(request["cores"]))
        torch.manual_seed(seed)
        arguments = dict(parameters)
        percentile = int(arguments.pop("top_gene_percentile", 50))
        workers = int(arguments.pop("workers", 1))
        arguments.setdefault("device", "cpu")
        if "hidden_dims" in arguments and not isinstance(arguments["hidden_dims"], list):
            arguments["hidden_dims"] = [arguments["hidden_dims"]]
        timesteps = int(arguments.get("T", 5000))
        start_noise = float(arguments.get("start_noise", 0.0001))
        end_noise = float(arguments.get("end_noise", 0.02))
        if timesteps < 1 or not 0 < start_noise <= end_noise < timesteps / 1000:
            raise ValueError("RegDiffusion requires 0 < start_noise <= end_noise < T/1000 for a finite diffusion schedule")
        trainer = rd.RegDiffusionTrainer(np.asarray(expression.T, dtype=np.float32), **arguments)
        trainer.train()
        graph = trainer.get_grn(list(expression.index), top_gene_percentile=percentile)
        edges = graph.extract_edgelist(k=-1, workers=workers)
        edges.columns = ["regulator", "target", "weight"]
        parameters = {name: value.default for name, value in inspect.signature(rd.RegDiffusionTrainer).parameters.items()
                      if name != "exp_array" and value.default is not inspect.Parameter.empty}
        parameters.update(arguments)
        parameters.update(top_gene_percentile=percentile, workers=workers)
        return edges, "regdiffusion::RegDiffusionTrainer", ["regdiffusion", "torch"], parameters
    raise ValueError(f"Unsupported RNA Python method: {method}")


def main(request_path):
    request = json.loads(Path(request_path).read_text())
    expression = pd.read_csv(request["expression"], sep="\t", dtype={"gene": str}, keep_default_na=False).set_index("gene")
    if expression.index.has_duplicates or expression.columns.has_duplicates:
        raise ValueError("Expression gene and cell identifiers must be unique")
    if not np.isfinite(expression.to_numpy(dtype=float)).all():
        raise ValueError("Expression values must be finite")
    edges, implementation, packages, parameters = infer(expression, request)
    edges = edges.loc[(edges["regulator"] != edges["target"]) &
                      edges["regulator"].isin(request["regulators"]) &
                      edges["target"].isin(request["targets"]), :]
    if not np.isfinite(edges["weight"].to_numpy(dtype=float)).all():
        raise ValueError("Native inference returned non-finite edge scores")
    edges.to_csv(request["output"], sep="\t", index=False)
    versions = {name: importlib.metadata.version(name) for name in ["numpy", "pandas", *packages]}
    receipt = dict(implementation=implementation, package_versions=versions,
                   parameters=parameters, python=sys.executable,
                   seed=request["seed"], cores=request["cores"], rows=len(edges))
    Path(request["receipt"]).write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    main(sys.argv[1])
