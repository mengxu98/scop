import sys

from _rna import main


def infer(expression, request):
    from arboreto.algo import grnboost2
    from distributed import Client, LocalCluster

    parameters = request["method_args"] or {}

    with LocalCluster(n_workers=int(request["cores"]), threads_per_worker=1,
                      processes=False, dashboard_address=None) as cluster:
        with Client(cluster) as client:
            edges = grnboost2(
                expression_data=expression.T, tf_names=request["regulators"],
                client_or_address=client, seed=int(request["seed"]),
                early_stop_window_length=int(parameters.get("early_stop_window_length", 25)),
                verbose=bool(request["verbose"]),
            )
    parameters = dict(early_stop_window_length=int(parameters.get("early_stop_window_length", 25)))
    return edges.rename(columns={"TF": "regulator", "importance": "weight"}), "arboreto::grnboost2", ["arboreto", "dask", "distributed"], parameters


if __name__ == "__main__":
    main(infer, sys.argv[1])
