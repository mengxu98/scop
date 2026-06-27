import os
from pathlib import Path


def _expand_path(path):
    if path is None:
        return None
    return str(Path(str(path)).expanduser())


def _expand_pathlike(value):
    if isinstance(value, (str, os.PathLike)):
        return _expand_path(value)
    return value


def _copy_obsm_value(value):
    try:
        return value.copy()
    except AttributeError:
        return value


def _as_list(value):
    if value is None:
        return None
    if isinstance(value, str):
        return [value]
    return [str(x) for x in list(value)]


def _read_anndata(test_input):
    if isinstance(test_input, (str, os.PathLike)):
        import scanpy as sc

        return sc.read_h5ad(_expand_path(test_input))
    return test_input


def _run_aucell(adata, gmt_file, norm_type=False):
    from scMalignantFinder import utils

    gmt_file = _expand_path(gmt_file)
    try:
        return utils.aucell_cal(adata, gmt_file, norm_type=bool(norm_type))
    except TypeError:
        return utils.aucell_cal(adata, gmt_file)


def _obs_frame(adata, columns=None, before=None):
    if columns is None:
        columns = []
    columns = [c for c in columns if c in adata.obs.columns]
    if not columns and before is not None:
        columns = [c for c in adata.obs.columns if c not in before]
    if not columns:
        raise ValueError("No expected scMalignantFinder result columns were returned")
    return adata.obs.loc[:, columns].copy()


def run_scmalignantfinder(
    test_input,
    pretrain_dir=None,
    train_h5ad_path=None,
    feature_path=None,
    model_method="LogisticRegression",
    norm_type=False,
    use_raw=False,
    n_thread=1,
    return_obs=True,
    verbose=True,
):
    from scMalignantFinder import classifier

    kwargs = {
        "test_input": _expand_pathlike(test_input),
        "pretrain_dir": _expand_path(pretrain_dir),
        "train_h5ad_path": _expand_path(train_h5ad_path),
        "feature_path": _expand_path(feature_path),
        "model_method": model_method,
        "norm_type": bool(norm_type),
        "use_raw": bool(use_raw),
        "n_thread": int(n_thread),
    }
    kwargs = {key: value for key, value in kwargs.items() if value is not None}

    try:
        model = classifier.scMalignantFinder(**kwargs)
    except TypeError as exc:
        if "use_raw" not in kwargs:
            raise
        kwargs.pop("use_raw", None)
        try:
            model = classifier.scMalignantFinder(**kwargs)
        except TypeError:
            raise exc

    model.load()
    result_adata = model.predict()

    if return_obs:
        return _obs_frame(
            result_adata,
            columns=["scMalignantFinder_prediction", "malignancy_probability"],
        )
    return result_adata


def run_scmalignant_region(
    test_input,
    signature_gmt,
    features=None,
    nclus=3,
    define_feature="Malignant_up",
    spatial_nn=True,
    spatial_coordinates=None,
    spatial_key="spatial",
    image=False,
    norm_type=False,
    return_obs=True,
    verbose=True,
):
    import numpy as np
    from scMalignantFinder import spatial as smf_spatial

    adata = _read_anndata(test_input)
    spatial_key = "spatial" if spatial_key is None else str(spatial_key)
    if spatial_coordinates is not None:
        spatial_coordinates = np.asarray(spatial_coordinates, dtype=float)
        if spatial_coordinates.shape[0] != adata.n_obs:
            raise ValueError(
                "spatial_coordinates must have one row for each observation"
            )
        adata.obsm[spatial_key] = spatial_coordinates

    adata = _run_aucell(adata, signature_gmt, norm_type=norm_type)
    if image:
        adata = smf_spatial.image_cal(adata)

    features = _as_list(features)
    if features is None:
        features = ["malignancy_probability", "Malignant_up"]
        if image:
            features.append("image_score")

    alias_spatial = bool(spatial_nn) and spatial_key != "spatial"
    had_spatial = False
    old_spatial = None
    if alias_spatial:
        if spatial_key not in adata.obsm:
            raise ValueError(f"spatial_key '{spatial_key}' was not found in adata.obsm")
        had_spatial = "spatial" in adata.obsm
        if had_spatial:
            old_spatial = _copy_obsm_value(adata.obsm["spatial"])
        adata.obsm["spatial"] = _copy_obsm_value(adata.obsm[spatial_key])

    try:
        adata = smf_spatial.region_identification(
            adata,
            features=features,
            nclus=int(nclus),
            define_feature=define_feature,
            spatial_nn=bool(spatial_nn),
        )
    finally:
        if alias_spatial:
            if had_spatial:
                adata.obsm["spatial"] = old_spatial
            elif "spatial" in adata.obsm:
                del adata.obsm["spatial"]

    if return_obs:
        preferred = ["cluster", "region_prediction", "Malignant_up", "image_score"]
        preferred.extend(features)
        return _obs_frame(adata, columns=list(dict.fromkeys(preferred)))
    return adata


def run_scmalignant_states(
    test_input,
    gene_sets,
    norm_type=False,
    return_obs=True,
    verbose=True,
):
    adata = _read_anndata(test_input)
    before = set(adata.obs.columns)
    adata = _run_aucell(adata, gene_sets, norm_type=norm_type)

    if return_obs:
        preferred = [c for c in adata.obs.columns if c not in before]
        if not preferred:
            preferred = [c for c in adata.obs.columns if str(c).startswith("MP")]
        return _obs_frame(adata, columns=preferred)
    return adata
