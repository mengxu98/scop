import inspect
import sys
from pathlib import Path

import numpy as np

from _rna import main


def infer(expression, request):
    import torch
    import regdiffusion as rd

    torch.set_num_threads(int(request["cores"]))
    torch.manual_seed(int(request["seed"]))
    arguments = dict(request["method_args"] or {})
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


if __name__ == "__main__":
    script_dir = Path(__file__).resolve().parent
    sys.path = [path for path in sys.path if Path(path).resolve() != script_dir]
    main(infer, sys.argv[1])
