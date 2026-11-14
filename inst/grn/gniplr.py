import ast
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from _rna import main


def infer(expression, request):
    args = request["method_args"] or {}
    degree = int(args.get("lasso_degree", 30))
    alpha = float(args.get("lasso_alpha", 0.04))
    if expression.shape[1] < 8 or degree < 1 or alpha < 0:
        raise ValueError("GNIPLR requires at least eight cells, positive lasso_degree and nonnegative lasso_alpha")
    path = Path(request["upstream_source"])
    tree = ast.parse(path.read_text())
    tree.body = [node for node in tree.body if isinstance(node, ast.FunctionDef) or
                 isinstance(node, ast.ImportFrom) or
                 (isinstance(node, ast.Import) and all(alias.name != "matplotlib.pyplot" for alias in node.names))]
    upstream = {"samples1": expression.columns.to_list()}
    exec(compile(tree, str(path), "exec"), upstream)
    vectors = {gene: row.to_list() for gene, row in expression.iterrows()}
    rows = []
    for regulator in request["regulators"]:
        for target in request["targets"]:
            if regulator == target:
                continue
            left, right = upstream["lasso_regress"](vectors[regulator], vectors[target], degree, alpha)
            p = upstream["granger"](left, right)
            rows.append((regulator, target, -np.log10(max(p, 1e-300))))
    edges = pd.DataFrame(rows, columns=["regulator", "target", "weight"])
    return edges.sort_values("weight", ascending=False, kind="stable"), "GNIPLR::gniplr.py", ["statsmodels", "scipy", "scikit-learn"], dict(lasso_degree=degree, lasso_alpha=alpha)


if __name__ == "__main__":
    main(infer, sys.argv[1])
