import importlib.metadata
import importlib.util
import json
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def load_upstream(request):
    path = Path(request["upstream_source"])
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(infer, request_path):
    request = json.loads(Path(request_path).read_text())
    expression = pd.read_csv(request["expression"], sep="\t", dtype={"gene": str}, keep_default_na=False).set_index("gene")
    seed = int(request["seed"])
    np.random.seed(seed)
    random.seed(seed)
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
