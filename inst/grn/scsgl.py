import sys
from pathlib import Path

import numpy as np
import pandas as pd

from _rna import load_upstream, main


def infer(expression, request):
    parameters = request["method_args"] or {}
    parameters = dict(pos_density=float(parameters.get("pos_density", 0.01)), neg_density=float(parameters.get("neg_density", 0.01)))
    if any(not np.isfinite(value) or value < 0 or value >= 1 for value in parameters.values()):
        raise ValueError("SCSGL densities must be finite, nonnegative and smaller than one")
    if not any(value > 0 for value in parameters.values()):
        raise ValueError("SCSGL requires at least one positive edge density")
    output = Path(request["work_dir"]) / "scsgl-native.tsv"
    load_upstream(request).run_inference(expression, output, **parameters)
    edges = pd.read_csv(output, sep="\t", dtype={"Gene1": str, "Gene2": str}, keep_default_na=False)
    edges = edges.rename(columns={"Gene1": "regulator", "Gene2": "target", "EdgeWeight": "weight"})
    return edges.sort_values("weight", ascending=False), "scSGL::signed_graph_learning", ["scipy", "numba"], parameters


if __name__ == "__main__":
    main(infer, sys.argv[1])
