import sys

from _rna import main


def infer(expression, request):
    from arboreto.algo import SGBM_KWARGS, grnboost2
    from distributed import Client, LocalCluster

    arguments = request["method_args"] or {}
    parameters = dict(early_stop_window_length=int(arguments.get("early_stop_window_length", 25)))
    for name, setting, cast in (
        ("n_rounds", "n_estimators", int),
        ("learning_rate", "learning_rate", float),
        ("max_depth", "max_depth", int),
        ("max_features", "max_features", float),
        ("subsample", "subsample", float),
    ):
        if name in arguments:
            SGBM_KWARGS[setting] = cast(arguments[name])
            parameters[name] = SGBM_KWARGS[setting]

    with LocalCluster(n_workers=int(request["cores"]), threads_per_worker=1,
                      processes=False, dashboard_address=None) as cluster:
        with Client(cluster) as client:
            edges = grnboost2(
                expression_data=expression.T, tf_names=request["regulators"],
                client_or_address=client, seed=int(request["seed"]),
                early_stop_window_length=parameters["early_stop_window_length"],
                verbose=bool(request["verbose"]),
            )
    return edges.rename(columns={"TF": "regulator", "importance": "weight"}), "arboreto::grnboost2", ["arboreto", "dask", "distributed"], parameters


if __name__ == "__main__":
    main(infer, sys.argv[1])
