import sys
from pathlib import Path

import numpy as np

from _rna import main


def infer(expression, request):
    sys.path.insert(0, str(Path(request["upstream_source"]).parents[2]))
    from pysrc.graphlearning import learn_signed_graph, convert_df, unsigned, _find_bs_upper_bound

    args = request["method_args"] or {}
    pos = float(args.get("pos_density", 0.01))
    neg = float(args.get("neg_density", 0.01))
    if any(not np.isfinite(value) or value < 0 or value >= 1 for value in (pos, neg)) or not (pos or neg):
        raise ValueError("SCSGL needs at least one positive edge density; both must be in [0, 1)")
    values = expression.to_numpy(dtype=float)
    genes = expression.index.to_numpy()
    keep = np.count_nonzero(values, axis=1) != 0
    if np.count_nonzero(keep) < 2:
        raise ValueError("SCSGL needs at least two nonzero genes")
    if pos and neg:
        edges = learn_signed_graph(values, pos, neg, assoc="correlation", gene_names=genes)
    else:
        from pysrc.associations import correlation

        kernel = correlation.calc(values[keep])
        k = kernel[np.triu_indices_from(kernel, k=1)]
        k /= np.max(np.abs(k))
        diag = np.diag(kernel)
        k = k if pos else -k
        density = pos or neg
        lower, upper = 0, _find_bs_upper_bound(k, diag, density)
        for _ in range(50):
            threshold = (lower + upper) / 2
            weights = unsigned.learn(k, diag, threshold)
            estimated = np.count_nonzero(weights) / len(weights)
            if abs(estimated - density) < 0.01:
                break
            if estimated < density:
                lower = threshold
            else:
                upper = threshold
        edges = convert_df(genes[keep], weights if pos else -weights, np.zeros_like(weights))
    edges = edges.rename(columns={"Gene1": "regulator", "Gene2": "target", "EdgeWeight": "weight"})
    return edges.sort_values("weight", ascending=False), "scSGL::learn_signed_graph", ["scipy", "numba"], dict(pos_density=pos, neg_density=neg)


if __name__ == "__main__":
    main(infer, sys.argv[1])
