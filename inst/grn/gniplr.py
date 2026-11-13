import inspect
import sys
from pathlib import Path

import pandas as pd

from _rna import load_upstream, main


def infer(expression, request):
    parameters = request["method_args"] or {}
    parameters = dict(lasso_degree=int(parameters.get("lasso_degree", 30)), lasso_alpha=float(parameters.get("lasso_alpha", 0.04)))
    if expression.shape[1] < 8 or parameters["lasso_degree"] < 1 or parameters["lasso_alpha"] < 0:
        raise ValueError("GNIPLR requires at least eight cells, positive lasso_degree and nonnegative lasso_alpha")
    infer_gniplr = load_upstream(request).infer_gniplr
    available = inspect.signature(infer_gniplr).parameters
    for name, default in (("lasso_degree", 30), ("lasso_alpha", 0.04)):
        if name not in available and parameters[name] != default:
            raise ValueError(f"The selected GNIPLR implementation fixes {name} at {default}; use an upstream implementation that exposes this argument.")
    output = Path(request["work_dir"]) / "gniplr-native.csv"
    infer_gniplr(expression, output, **{name: value for name, value in parameters.items() if name in available})
    return pd.read_csv(output, dtype={"regulator": str, "target": str}, keep_default_na=False), "multiCSN-manuscript::GNIPLR", ["statsmodels", "scipy", "scikit-learn"], parameters


if __name__ == "__main__":
    main(infer, sys.argv[1])
