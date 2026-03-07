import importlib.util
import os
from pathlib import Path

_LOG_MESSAGE_PATH = Path(__file__).resolve().parent / "log_message.py"
_LOG_MESSAGE_SPEC = importlib.util.spec_from_file_location(
    "scop_log_message", _LOG_MESSAGE_PATH
)
if _LOG_MESSAGE_SPEC is None or _LOG_MESSAGE_SPEC.loader is None:
    raise ImportError(f"Cannot load log_message module from {_LOG_MESSAGE_PATH}")

_LOG_MESSAGE_MODULE = importlib.util.module_from_spec(_LOG_MESSAGE_SPEC)
_LOG_MESSAGE_SPEC.loader.exec_module(_LOG_MESSAGE_MODULE)

log_message = _LOG_MESSAGE_MODULE.log_message


def configure_apple_silicon_env(
    scanpy_settings=False,
    scanpy_verbosity=False,
    numba_threading=False,
    numba_disable_jit=False,
    configure_numba_runtime=False,
    numba_runtime_optional=False,
    verbose=True,
):
    log_message(
        "Apple silicon detected: Applying specific configurations",
        message_type="info",
        verbose=verbose,
    )
    os.environ["PYTHONHASHSEED"] = "0"
    os.environ["PYTHONUNBUFFERED"] = "1"
    os.environ["MPLBACKEND"] = "Agg"
    os.environ["DISPLAY"] = ""

    if scanpy_settings:
        os.environ["SCANPY_SETTINGS"] = "scanpy_settings"
    if scanpy_verbosity:
        os.environ["SCANPY_SETTINGS_VERBOSITY"] = "1"
    if numba_disable_jit:
        os.environ["NUMBA_DISABLE_JIT"] = "1"
    if numba_threading:
        os.environ["NUMBA_NUM_THREADS"] = "1"
        os.environ["NUMBA_THREADING_LAYER"] = "tbb"
        os.environ["NUMBA_DEFAULT_NUM_THREADS"] = "1"

    if configure_numba_runtime:
        try:
            import numba

            numba.config.DISABLE_JIT = True
            numba.set_num_threads(1)
            log_message(
                "{.pkg NUMBA} configured for Apple silicon",
                message_type="success",
                verbose=verbose,
            )
        except ImportError:
            if not numba_runtime_optional:
                raise

    try:
        import ctypes
        import sys

        prefixes = []
        conda_prefix = os.environ.get("CONDA_PREFIX")
        if conda_prefix:
            prefixes.append(conda_prefix)
        if sys.prefix:
            prefixes.append(sys.prefix)

        prefixes = list(dict.fromkeys(prefixes))
        lib_dirs = [str(Path(p) / "lib") for p in prefixes if p]

        existing_fallback = os.environ.get("DYLD_FALLBACK_LIBRARY_PATH", "")
        existing_library = os.environ.get("DYLD_LIBRARY_PATH", "")
        merged_fallback = ":".join(
            lib_dirs + ([existing_fallback] if existing_fallback else [])
        )
        merged_library = ":".join(
            lib_dirs + ([existing_library] if existing_library else [])
        )
        if merged_fallback:
            os.environ["DYLD_FALLBACK_LIBRARY_PATH"] = merged_fallback
        if merged_library:
            os.environ["DYLD_LIBRARY_PATH"] = merged_library

        loaded_omp = False
        for prefix in prefixes:
            lib_dir = Path(prefix) / "lib"
            for lib_name in ("libomp.dylib", "libiomp5.dylib"):
                lib_path = lib_dir / lib_name
                if lib_path.exists():
                    try:
                        ctypes.CDLL(str(lib_path), mode=ctypes.RTLD_GLOBAL)
                        log_message(
                            "Loaded OpenMP runtime from {.val {%s}}" % lib_path,
                            message_type="info",
                            verbose=verbose,
                        )
                        loaded_omp = True
                        break
                    except Exception:
                        pass
            if loaded_omp:
                break
    except Exception:
        pass


def SCVELO(
    adata=None,
    h5ad=None,
    group_by=None,
    palette=None,
    linear_reduction=None,
    nonlinear_reduction=None,
    basis=None,
    mode=["deterministic", "stochastic", "dynamical"],
    fitting_by="stochastic",
    magic_impute=False,
    knn=5,
    t=2,
    min_shared_counts=30,
    n_pcs=30,
    n_neighbors=30,
    filter_genes=True,
    min_counts=3,
    min_counts_u=3,
    normalize_per_cell=True,
    log_transform=True,
    use_raw=False,
    diff_kinetics=False,
    stream_smooth=None,
    stream_density=2,
    arrow_length=5,
    arrow_size=5,
    arrow_density=0.5,
    denoise=False,
    denoise_topn=3,
    kinetics=False,
    kinetics_topn=100,
    calculate_velocity_genes=False,
    compute_velocity_confidence=True,
    compute_terminal_states=True,
    compute_pseudotime=True,
    compute_paga=True,
    top_n=6,
    n_jobs=1,
    show_plot=True,
    save_plot=False,
    plot_format="png",
    plot_dpi=600,
    plot_prefix="scvelo",
    dirpath="./scvelo",
    save=False,
    dpi=300,
    fileprefix="",
    verbose=True,
    legend_loc="on data",
):
    import os
    import platform

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["KMP_WARNINGS"] = "0"
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

    is_apple_silicon = platform.system() == "Darwin" and platform.machine() == "arm64"

    if is_apple_silicon:
        configure_apple_silicon_env(
            scanpy_settings=True,
            configure_numba_runtime=True,
            verbose=verbose,
        )

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scvelo as scv
    import scanpy as sc
    import pandas as pd
    import numpy as np
    from scipy import sparse
    import warnings

    warnings.simplefilter("ignore", category=UserWarning)
    warnings.simplefilter("ignore", category=FutureWarning)
    warnings.simplefilter("ignore", category=DeprecationWarning)

    prevdir = os.getcwd()

    if save_plot:
        expanded_path = os.path.expanduser(dirpath)
        if not os.path.exists(expanded_path):
            os.makedirs(expanded_path, exist_ok=True)
            log_message(
                "Created directory: {.val {expanded_path}}",
                message_type="info",
                verbose=verbose,
            )
        os.chdir(expanded_path)
        sc.settings.figdir = "."
    else:
        sc.settings.figdir = "."
    sc.settings.file_format_figs = plot_format
    sc.settings.autosave = save_plot

    log_message(
        "Starting {.pkg scVelo} analysis...", message_type="running", verbose=verbose
    )
    try:
        if adata is None and h5ad is None:
            raise ValueError("Either {.arg adata} or {.arg h5ad} must be provided")

        if adata is None:
            adata = scv.read(h5ad)

        if group_by is None:
            raise ValueError("{.arg group_by} must be provided")

        if linear_reduction is None and nonlinear_reduction is None:
            raise ValueError(
                "At least one of {.arg linear_reduction} or {.arg nonlinear_reduction} must be provided"
            )

        if basis is None:
            if nonlinear_reduction is not None:
                if nonlinear_reduction in adata.obsm:
                    basis = nonlinear_reduction
                else:
                    log_message(
                        "{.val {nonlinear_reduction}} not found in adata.obsm. Available keys: {.val {list(adata.obsm.keys())}}",
                        message_type="error",
                        verbose=verbose,
                    )
                    raise ValueError(
                        f"nonlinear_reduction '{nonlinear_reduction}' not found in adata.obsm"
                    )
            else:
                if linear_reduction in adata.obsm:
                    basis = linear_reduction
                else:
                    log_message(
                        "{.val {linear_reduction}} not found in adata.obsm. Available keys: {.val {list(adata.obsm.keys())}}",
                        message_type="error",
                        verbose=verbose,
                    )
                    raise ValueError(
                        f"linear_reduction '{linear_reduction}' not found in adata.obsm"
                    )

        if basis not in adata.obsm:
            log_message(
                "basis {.val {basis}} not found in adata.obsm. Available keys: {.val {list(adata.obsm.keys())}}",
                message_type="error",
                verbose=verbose,
            )
            if linear_reduction in adata.obsm:
                adata.obsm["basis"] = adata.obsm[linear_reduction][:, 0:2]
                basis = "basis"
            else:
                raise ValueError(
                    f"Cannot find suitable basis. Available obsm keys: {list(adata.obsm.keys())}"
                )

        log_message("Using basis: {.val {basis}}", message_type="info", verbose=verbose)
        log_message(
            "Available embeddings in adata.obsm: {.val {list(adata.obsm.keys())}}",
            message_type="info",
            verbose=verbose,
        )

        adata.obs[group_by] = adata.obs[group_by].astype("category")

        log_message("Starting preprocessing", message_type="running", verbose=verbose)

        if filter_genes:
            log_message("Filtering genes...", message_type="info", verbose=verbose)
            scv.pp.filter_genes(adata, min_counts=min_counts)
            scv.pp.filter_genes(adata, min_counts_u=min_counts_u)

        if normalize_per_cell:
            log_message("Normalizing per cell...", message_type="info", verbose=verbose)
            scv.pp.normalize_per_cell(adata)

        if log_transform:
            log_message("Log transforming...", message_type="info", verbose=verbose)
            sc.pp.log1p(adata)

        if magic_impute:
            log_message(
                "Performing {.pkg magic-impute} imputation...",
                message_type="info",
                verbose=verbose,
            )
            try:
                import magic

                magic_operator = magic.MAGIC(knn=knn, t=t, random_state=42)
                adata.layers["spliced_raw"] = adata.layers["spliced"].copy()
                adata.layers["unspliced_raw"] = adata.layers["unspliced"].copy()
                adata.layers["spliced"] = magic_operator.fit_transform(
                    adata.layers["spliced"]
                )
                adata.layers["unspliced"] = magic_operator.transform(
                    adata.layers["unspliced"]
                )
            except ImportError:
                log_message(
                    "{.pkg magic-impute} not installed. Skipping imputation",
                    message_type="warning",
                    verbose=verbose,
                )

        log_message(
            "Computing {.pkg scvelo} neighbors and moments...",
            message_type="info",
            verbose=verbose,
        )

        use_rep = None
        if linear_reduction in adata.obsm:
            rep_dims = adata.obsm[linear_reduction].shape[1]
            if rep_dims >= n_pcs:
                use_rep = linear_reduction
                log_message(
                    "Using {.val {linear_reduction}} with {.val {rep_dims}} dimensions",
                    message_type="info",
                    verbose=verbose,
                )
            else:
                log_message(
                    "{.val {linear_reduction}} has only {.val {rep_dims}} dimensions, need {.val {n_pcs}}",
                    message_type="warning",
                    verbose=verbose,
                )

        if use_rep is None:
            pca_candidates = [
                key
                for key in adata.obsm.keys()
                if "pca" in key.lower() and "umap" not in key.lower()
            ]
            for candidate in pca_candidates:
                rep_dims = adata.obsm[candidate].shape[1]
                if rep_dims >= n_pcs:
                    use_rep = candidate
                    log_message(
                        "Using PCA representation {.val {candidate}} with {.val {rep_dims}} dimensions",
                        message_type="info",
                        verbose=verbose,
                    )
                    break

        if use_rep is None:
            max_dims = 0
            for key in adata.obsm.keys():
                if "pca" in key.lower() and "umap" not in key.lower():
                    dims = adata.obsm[key].shape[1]
                    if dims > max_dims:
                        max_dims = dims
                        use_rep = key

            if use_rep and max_dims > 0:
                n_pcs = min(n_pcs, max_dims)
                log_message(
                    "Reducing {.arg n_pcs} to {.val {n_pcs}} to match available dimensions in {.val {use_rep}}",
                    message_type="info",
                    verbose=verbose,
                )
            else:
                use_rep = None
                n_pcs = min(n_pcs, adata.X.shape[1])
                log_message(
                    "Using raw data with {.arg n_pcs}={.val {n_pcs}}",
                    message_type="info",
                    verbose=verbose,
                )

        try:
            scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors=n_neighbors, use_rep=use_rep)
        except Exception as e:
            log_message(
                "{.pkg scvelo} moments failed ({.val {e}}), using manual computation...",
                message_type="warning",
                verbose=verbose,
            )
            sc.pp.neighbors(
                adata, n_pcs=n_pcs, n_neighbors=n_neighbors, use_rep=use_rep
            )

            connectivities = adata.obsp["connectivities"]
            if sparse.issparse(adata.layers["spliced"]):
                Ms = connectivities @ adata.layers["spliced"]
                Mu = connectivities @ adata.layers["unspliced"]
            else:
                Ms = connectivities @ sparse.csr_matrix(adata.layers["spliced"])
                Mu = connectivities @ sparse.csr_matrix(adata.layers["unspliced"])

            adata.layers["Ms"] = Ms
            adata.layers["Mu"] = Mu

        log_message(
            "Starting {.pkg velocity} estimation",
            message_type="running",
            verbose=verbose,
        )

        for m in mode:
            log_message(
                "Processing mode: {.val {m}}", message_type="info", verbose=verbose
            )

            if m == "dynamical":
                log_message(
                    "Performing {.pkg velocity} dynamical modeling...",
                    message_type="info",
                    verbose=verbose,
                )
                gene_subset = (
                    adata.var[fitting_by + "_genes"]
                    if (fitting_by + "_genes") in adata.var.columns
                    else None
                )

                if gene_subset is not None and gene_subset.sum() > 0:
                    scv.tl.recover_dynamics(
                        adata,
                        var_names=gene_subset,
                        use_raw=use_raw,
                        n_jobs=n_jobs,
                    )
                else:
                    log_message(
                        "No genes found for dynamical modeling. Using all genes",
                        message_type="warning",
                        verbose=verbose,
                    )
                    scv.tl.recover_dynamics(adata, use_raw=use_raw, n_jobs=n_jobs)

            log_message(
                "Computing {.pkg velocity} {.val {m}} velocity...",
                message_type="info",
                verbose=verbose,
            )
            scv.tl.velocity(adata, mode=m, diff_kinetics=diff_kinetics)

            log_message(
                "Computing {.pkg velocity} {.val {m}} graph...",
                message_type="info",
                verbose=verbose,
            )
            scv.tl.velocity_graph(adata, vkey=m, n_neighbors=n_neighbors, n_jobs=n_jobs)

            log_message(
                "Downstream analysis for {.pkg velocity} {.val {m}}",
                message_type="running",
                verbose=verbose,
            )

            log_message(
                "Computing {.pkg velocity} embedding...",
                message_type="info",
                verbose=verbose,
            )
            if basis not in adata.obsm:
                log_message(
                    "Basis {.val {basis}} not found in adata.obsm. Available keys: {.val {list(adata.obsm.keys())}}",
                    message_type="error",
                    verbose=verbose,
                )
                raise ValueError(
                    f"Basis '{basis}' not found in adata.obsm. Available keys: {list(adata.obsm.keys())}"
                )

            basis_embedding = adata.obsm[basis]
            if basis_embedding.shape[1] < 2:
                log_message(
                    "Basis {.val {basis}} has only {.val {basis_embedding.shape[1]}} dimensions, need at least 2 for velocity embedding",
                    message_type="error",
                    verbose=verbose,
                )
                raise ValueError(
                    f"Basis '{basis}' must have at least 2 dimensions for velocity embedding"
                )

            x_basis_key = f"X_{basis}"
            if x_basis_key not in adata.obsm:
                adata.obsm[x_basis_key] = basis_embedding
                log_message(
                    "Created {.val {x_basis_key}} in adata.obsm for scvelo compatibility",
                    message_type="info",
                    verbose=verbose,
                )

            scv.tl.velocity_embedding(adata, basis=basis, vkey=m)

            if compute_velocity_confidence:
                log_message(
                    "Computing {.pkg velocity} confidence...",
                    message_type="info",
                    verbose=verbose,
                )
                try:
                    scv.tl.velocity_confidence(adata, vkey=m)
                except Exception as e:
                    log_message(
                        "velocity confidence failed ({.val {e}}), using default values",
                        message_type="warning",
                        verbose=verbose,
                    )
                    n_obs = adata.n_obs
                    adata.obs[m + "_length"] = np.ones(n_obs) * 0.5
                    adata.obs[m + "_confidence"] = np.ones(n_obs) * 0.5

            if compute_terminal_states:
                log_message(
                    "Computing terminal states...", message_type="info", verbose=verbose
                )
                try:
                    scv.tl.terminal_states(adata, vkey=m)
                    for term in ["root_cells", "end_points"]:
                        if term in adata.obs.columns:
                            adata.obs[m + "_" + term] = adata.obs[term]
                            adata.obs.drop(term, axis=1, inplace=True)
                except Exception as e:
                    log_message(
                        "Terminal states computation failed: {.val {e}}",
                        message_type="warning",
                        verbose=verbose,
                    )

            if compute_pseudotime:
                log_message(
                    "Computing {.pkg velocity} pseudotime...",
                    message_type="info",
                    verbose=verbose,
                )
                try:
                    root_key = (
                        m + "_root_cells"
                        if (m + "_root_cells") in adata.obs.columns
                        else None
                    )
                    end_key = (
                        m + "_end_points"
                        if (m + "_end_points") in adata.obs.columns
                        else None
                    )
                    scv.tl.velocity_pseudotime(
                        adata, vkey=m, root_key=root_key, end_key=end_key
                    )
                except Exception as e:
                    log_message(
                        "Pseudotime computation failed: {.val {e}}",
                        message_type="warning",
                        verbose=verbose,
                    )

            if compute_paga:
                log_message(
                    "Computing {.pkg PAGA}...", message_type="info", verbose=verbose
                )
                try:
                    if "neighbors" not in adata.uns:
                        adata.uns["neighbors"] = {}
                    adata.uns["neighbors"]["distances"] = adata.obsp["distances"]
                    adata.uns["neighbors"]["connectivities"] = adata.obsp[
                        "connectivities"
                    ]

                    if m + "_graph" in adata.uns:
                        adata.uns["velocity_graph"] = adata.uns[m + "_graph"]

                    adata.obs[group_by] = adata.obs[group_by].astype(dtype="category")

                    sc.tl.paga(adata, groups=group_by, use_rna_velocity=True)
                except Exception as e:
                    log_message(
                        "{.pkg PAGA} computation failed ({.val {e}})",
                        message_type="warning",
                        verbose=verbose,
                    )

            if calculate_velocity_genes:
                log_message(
                    "Ranking {.pkg velocity} genes...",
                    message_type="info",
                    verbose=verbose,
                )
                try:
                    if m != "dynamical":
                        scv.tl.rank_velocity_genes(adata, vkey=m, groupby=group_by)
                        if "spearmans_score" in adata.var.columns:
                            adata.var[m + "_score"] = adata.var["spearmans_score"]
                    else:
                        scv.tl.rank_dynamical_genes(adata, groupby=group_by)
                except Exception as e:
                    log_message(
                        "{.pkg velocity} genes ranking failed ({.val {e}})",
                        message_type="warning",
                        verbose=verbose,
                    )

            if show_plot:
                log_message(
                    "Generating plots for {.pkg {m}}...",
                    message_type="info",
                    verbose=verbose,
                )

                groups = (
                    adata.obs[group_by].cat.categories
                    if hasattr(adata.obs[group_by], "cat")
                    else adata.obs[group_by].unique()
                )
                if palette is None:
                    palette = dict(
                        zip(groups, plt.cm.tab10(np.linspace(0, 1, len(groups))))
                    )

                def get_scvelo_save_path(name):
                    if save_plot:
                        ext = (
                            plot_format
                            if plot_format in ["png", "pdf", "svg"]
                            else "png"
                        )
                        save_path = os.path.abspath(f"{plot_prefix}_{m}_{name}.{ext}")
                        return save_path
                    return False

                try:
                    scv.pl.velocity_embedding_stream(
                        adata,
                        vkey=m,
                        basis=basis,
                        title=f"{m} velocity",
                        color=group_by,
                        palette=palette,
                        smooth=stream_smooth,
                        density=stream_density,
                        legend_loc=legend_loc,
                        save=False,
                        show=show_plot,
                        dpi=plot_dpi,
                    )
                    if save_plot:
                        save_path = get_scvelo_save_path("stream")
                        try:
                            plt.savefig(save_path, dpi=plot_dpi, bbox_inches="tight")
                        except Exception as pdf_error:
                            if os.path.exists(save_path):
                                try:
                                    os.remove(save_path)
                                except Exception:
                                    pass
                            save_path_png = save_path.replace(".pdf", ".png")
                            plt.savefig(
                                save_path_png, dpi=plot_dpi, bbox_inches="tight"
                            )
                            log_message(
                                "PDF save failed, saved as PNG instead: {.val {save_path_png}}",
                                message_type="warning",
                                verbose=verbose,
                            )
                    plt.close()
                except Exception as e:
                    log_message(
                        "Stream plot failed ({.val {e}})",
                        message_type="warning",
                        verbose=verbose,
                    )
                    if plt is not None:
                        plt.close()

                try:
                    scv.pl.velocity_embedding(
                        adata,
                        vkey=m,
                        basis=basis,
                        title=f"{m} velocity",
                        color=group_by,
                        palette=palette,
                        arrow_length=arrow_length,
                        arrow_size=arrow_size,
                        density=arrow_density,
                        linewidth=0.3,
                        legend_loc=legend_loc,
                        save=False,
                        show=show_plot,
                        dpi=plot_dpi,
                    )
                    if save_plot:
                        save_path = get_scvelo_save_path("arrow")
                        try:
                            plt.savefig(save_path, dpi=plot_dpi, bbox_inches="tight")
                        except Exception as pdf_error:
                            if os.path.exists(save_path):
                                try:
                                    os.remove(save_path)
                                except Exception:
                                    pass
                            save_path_png = save_path.replace(".pdf", ".png")
                            plt.savefig(
                                save_path_png, dpi=plot_dpi, bbox_inches="tight"
                            )
                            log_message(
                                "PDF save failed, saved as PNG instead: {.val {save_path_png}}",
                                message_type="warning",
                                verbose=verbose,
                            )
                    plt.close()
                except Exception as e:
                    log_message(
                        "Arrow plot failed ({.val {e}})",
                        message_type="warning",
                        verbose=verbose,
                    )
                    if plt is not None:
                        plt.close()

                if compute_velocity_confidence:
                    for metric in ["length", "confidence"]:
                        color_key = f"{m}_{metric}"
                        if color_key in adata.obs.columns:
                            try:
                                scv.pl.scatter(
                                    adata,
                                    basis=basis,
                                    color=color_key,
                                    title=f"{m} {metric}",
                                    cmap="viridis",
                                    legend_loc="right margin",
                                    save=get_scvelo_save_path(metric),
                                    show=show_plot,
                                )
                            except Exception as e:
                                log_message(
                                    "{.pkg {metric}} plot failed ({.val {e}})",
                                    message_type="warning",
                                    verbose=verbose,
                                )

        log_message(
            "{.pkg scVelo} analysis completed", message_type="success", verbose=verbose
        )

    except Exception as e:
        log_message(
            "Error in {.pkg scVelo} analysis: {.val {e}}",
            message_type="error",
            verbose=verbose,
        )
        raise
    finally:
        try:
            figures_dir = os.path.join(os.getcwd(), "figures")
            if os.path.exists(figures_dir) and os.path.isdir(figures_dir):
                if not os.listdir(figures_dir):
                    os.rmdir(figures_dir)
                    log_message(
                        "Removed empty figures directory: {.val {figures_dir}}",
                        message_type="info",
                        verbose=verbose,
                    )
        except Exception:
            pass

        if save_plot:
            os.chdir(prevdir)

    try:
        if hasattr(adata, "_raw") and adata._raw is not None:
            if hasattr(adata._raw, "_var"):
                adata._raw._var = adata._raw._var.rename(columns={"_index": "features"})
    except Exception:
        pass

    return adata


def compute_transition_matrix(kernel, verbose=True):
    """
    Compute transition matrix

    Parameters
    ----------
    kernel : cellrank.kernels.Kernel
        The kernel object
    verbose : bool
        Whether to log messages

    Returns
    -------
    bool
        True if matrix was computed, False otherwise
    """
    try:
        kernel.compute_transition_matrix()
        return fix_transition_matrix(kernel, verbose=verbose)
    except ValueError as e:
        if "not row stochastic" in str(e):
            if verbose:
                log_message(
                    "Transition matrix validation failed: {.val {e}}. Attempting to fix...",
                    message_type="warning",
                    verbose=verbose,
                )
            try:
                import numpy as np
                from scipy.sparse import issparse, csr_matrix

                if hasattr(kernel, "connectivity"):
                    conn = kernel.connectivity
                elif hasattr(kernel, "_conn"):
                    conn = kernel._conn
                else:
                    if (
                        hasattr(kernel, "adata")
                        and "connectivities" in kernel.adata.obsp
                    ):
                        conn = kernel.adata.obsp["connectivities"]
                    else:
                        raise ValueError("Cannot access connectivity matrix")

                if issparse(conn):
                    row_sums = np.array(conn.sum(axis=1)).flatten()
                    row_sums_safe = np.maximum(row_sums, 1e-10)
                    from scipy.sparse import diags

                    tmat = diags(1.0 / row_sums_safe) @ conn
                else:
                    row_sums = conn.sum(axis=1, keepdims=True)
                    row_sums_safe = np.maximum(row_sums, 1e-10)
                    tmat = conn / row_sums_safe

                kernel._transition_matrix = tmat

                fix_transition_matrix(kernel, verbose=verbose)

                if verbose:
                    log_message(
                        "Transition matrix fixed successfully",
                        message_type="success",
                        verbose=verbose,
                    )
                return True
            except Exception as fix_error:
                if verbose:
                    log_message(
                        "Failed to fix transition matrix: {.val {fix_error}}",
                        message_type="error",
                        verbose=verbose,
                    )
                raise
        else:
            raise


def fix_transition_matrix(kernel, verbose=True):
    """
    Fix transition matrix to be row stochastic (rows sum to 1).

    This function ensures the transition matrix is valid for CellRank estimators
    by cleaning NaN/Inf values, clipping negatives, adding self-loops where needed,
    and normalizing rows to sum to 1.

    Parameters
    ----------
    kernel : cellrank.kernels.Kernel
        The kernel object with a transition_matrix attribute
    verbose : bool
        Whether to log messages

    Returns
    -------
    bool
        True if matrix was modified, False otherwise
    """
    try:
        import numpy as np
        from scipy.sparse import diags, issparse, csr_matrix

        tmat = kernel.transition_matrix
        matrix_modified = False

        if verbose:
            log_message(
                "Validating and fixing transition matrix...",
                message_type="info",
                verbose=verbose,
            )

        if issparse(tmat):
            if np.any(np.isnan(tmat.data)) or np.any(np.isinf(tmat.data)):
                if verbose:
                    log_message(
                        "Cleaning NaN/Inf values...",
                        message_type="warning",
                        verbose=verbose,
                    )
                tmat.data[np.isnan(tmat.data)] = 0
                tmat.data[np.isinf(tmat.data)] = 0
                tmat.eliminate_zeros()
                matrix_modified = True
        else:
            if np.any(np.isnan(tmat)) or np.any(np.isinf(tmat)):
                if verbose:
                    log_message(
                        "Cleaning NaN/Inf values...",
                        message_type="warning",
                        verbose=verbose,
                    )
                tmat[np.isnan(tmat)] = 0
                tmat[np.isinf(tmat)] = 0
                matrix_modified = True

        if issparse(tmat):
            if np.any(tmat.data < 0):
                if verbose:
                    log_message(
                        "Clipping {.val {(tmat.data < 0).sum()}} negative values to zero...",
                        message_type="warning",
                        verbose=verbose,
                    )
                tmat.data = np.maximum(tmat.data, 0)
                matrix_modified = True
        else:
            if np.any(tmat < 0):
                if verbose:
                    log_message(
                        "Clipping negative values to zero...",
                        message_type="warning",
                        verbose=verbose,
                    )
                tmat = np.maximum(tmat, 0)
                matrix_modified = True

        if issparse(tmat):
            diag_values = tmat.diagonal()
            min_self_loop = 0.01
            needs_self_loop = diag_values < min_self_loop

            if np.any(needs_self_loop):
                if verbose:
                    log_message(
                        "Adding self-loops to {.val {needs_self_loop.sum()}} cells for matrix primitivity...",
                        message_type="info",
                        verbose=verbose,
                    )

                tmat = tmat.tolil()
                for i in np.where(needs_self_loop)[0]:
                    tmat[i, i] = min_self_loop
                tmat = tmat.tocsr()
                matrix_modified = True

        row_sums = np.array(tmat.sum(axis=1)).flatten()
        zero_rows = row_sums < 1e-10

        if np.any(zero_rows):
            if verbose:
                log_message(
                    "Found {.val {zero_rows.sum()}} zero rows. Setting to self-loops...",
                    message_type="warning",
                    verbose=verbose,
                )

            if issparse(tmat):
                tmat = tmat.tolil()
                for i in np.where(zero_rows)[0]:
                    tmat[i, i] = 1.0
                tmat = tmat.tocsr()
            else:
                for i in np.where(zero_rows)[0]:
                    tmat[i, :] = 0
                    tmat[i, i] = 1.0

            row_sums = np.array(tmat.sum(axis=1)).flatten()
            matrix_modified = True

        if not np.allclose(row_sums, 1.0, rtol=1e-6, atol=1e-8):
            log_message(
                "Normalizing rows (current range: {.val {format(row_sums.min(), '.8f')}} - {.val {format(row_sums.max(), '.8f')}})...",
                message_type="info",
                verbose=verbose,
            )

            row_sums_safe = np.maximum(row_sums, 1e-10)

            if issparse(tmat):
                tmat = diags(1.0 / row_sums_safe) @ tmat
            else:
                tmat = tmat / row_sums_safe[:, np.newaxis]

            matrix_modified = True

        final_row_sums = np.array(tmat.sum(axis=1)).flatten()

        if issparse(tmat):
            has_negative = np.any(tmat.data < -1e-10)
        else:
            has_negative = np.any(tmat < -1e-10)

        rows_sum_to_one = np.allclose(final_row_sums, 1.0, rtol=1e-6, atol=1e-8)

        if has_negative and verbose:
            log_message(
                "Warning: Matrix still has negative values after fixing",
                message_type="warning",
                verbose=verbose,
            )

        if not rows_sum_to_one and verbose:
            log_message(
                "Warning: Matrix rows still don't sum to 1 (range: {.val {format(final_row_sums.min(), '.8f')}} - {.val {format(final_row_sums.max(), '.8f')}})",
                message_type="warning",
                verbose=verbose,
            )

        if matrix_modified:
            kernel._transition_matrix = tmat
            if verbose:
                log_message(
                    "Matrix fixed and validated (row sums: {.val {format(final_row_sums.min(), '.8f')}} - {.val {format(final_row_sums.max(), '.8f')}})",
                    message_type="success",
                    verbose=verbose,
                )
        elif verbose:
            log_message(
                "Matrix validation passed (row sums: {.val {format(final_row_sums.min(), '.8f')}} - {.val {format(final_row_sums.max(), '.8f')}})",
                message_type="success",
                verbose=verbose,
            )

        return matrix_modified

    except Exception as e:
        if verbose:
            log_message(
                "Matrix validation encountered error: {.val {e}}. Proceeding with original matrix...",
                message_type="warning",
                verbose=verbose,
            )
        return False


def CellRank(
    adata=None,
    h5ad=None,
    group_by=None,
    palette=None,
    linear_reduction=None,
    nonlinear_reduction=None,
    basis=None,
    mode=["deterministic", "stochastic", "dynamical"],
    fitting_by="stochastic",
    magic_impute=False,
    knn=5,
    t=2,
    min_shared_counts=30,
    n_pcs=30,
    n_neighbors=30,
    stream_smooth=None,
    stream_density=2,
    arrow_size=5,
    arrow_length=5,
    arrow_density=0.5,
    calculate_velocity_genes=False,
    denoise=False,
    kinetics=False,
    n_jobs=1,
    show_plot=True,
    dpi=300,
    save=False,
    dirpath="./cellrank",
    fileprefix="",
    verbose=True,
    kernel_type="velocity",
    time_key="dpt_pseudotime",
    estimator_type="GPCCA",
    use_connectivity_kernel=True,
    velocity_weight=0.8,
    connectivity_weight=0.2,
    softmax_scale=4,
    n_macrostates=None,
    schur_method="krylov",
    n_cells_terminal=10,
    save_plot=False,
    plot_format="png",
    plot_dpi=600,
    plot_prefix="cellrank",
    legend_loc="on data",
):
    import os
    import platform

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["KMP_WARNINGS"] = "0"
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

    is_apple_silicon = platform.system() == "Darwin" and platform.machine() == "arm64"

    if is_apple_silicon:
        configure_apple_silicon_env(
            scanpy_settings=True,
            configure_numba_runtime=True,
            verbose=verbose,
        )

    import matplotlib.pyplot as plt
    import scvelo as scv
    import cellrank as cr
    import pandas as pd
    import numpy as np
    import scanpy as sc

    import warnings

    warnings.simplefilter("ignore", category=UserWarning)
    warnings.simplefilter("ignore", category=FutureWarning)
    warnings.simplefilter("ignore", category=DeprecationWarning)

    prevdir = os.getcwd()

    if save_plot:
        expanded_path = os.path.expanduser(dirpath)
        if not os.path.exists(expanded_path):
            os.makedirs(expanded_path, exist_ok=True)
            log_message(
                "Created directory: {.val {expanded_path}}",
                message_type="info",
                verbose=verbose,
            )
        os.chdir(expanded_path)
        sc.settings.figdir = "."
        cr.settings.figdir = "."
    else:
        sc.settings.figdir = "."
        cr.settings.figdir = "."
    sc.settings.file_format_figs = plot_format
    sc.settings.autosave = save_plot
    log_message(
        "{.pkg CellRank} figdir set to: {.val {cr.settings.figdir}}",
        message_type="info",
        verbose=verbose,
    )

    if platform.system() == "Windows":
        import sys, multiprocessing, re

        if re.match(pattern=".*pythonw.exe$", string=sys.executable):
            pythonw = sys.executable
        else:
            pythonw = sys.executable.replace("python.exe", "pythonw.exe")
        sys.executable = pythonw
        sys._base_executable = pythonw
        multiprocessing.set_executable(pythonw)

    try:
        if adata is None and h5ad is None:
            log_message(
                "{.arg adata} or {.arg h5ad} must be provided", message_type="error"
            )
            exit()

        if adata is None:
            adata = scv.read(h5ad)

        if group_by is None:
            log_message("{.arg group_by} must be provided", message_type="error")
            exit()

        if linear_reduction is None and nonlinear_reduction is None:
            log_message(
                "{.arg linear_reduction} or {.arg nonlinear_reduction} must be provided at least one",
                message_type="error",
            )
            exit()

        if basis is None:
            if nonlinear_reduction is not None:
                basis = nonlinear_reduction
            else:
                basis = "basis"
                adata.obsm["basis"] = adata.obsm[linear_reduction][:, 0:2]

        mode.append(fitting_by)
        if kinetics is True or denoise is True:
            mode.append("dynamical")

        mode = list(set(mode))
        if "dynamical" in mode:
            mode.sort(key="dynamical".__eq__)

        if not fitting_by in ["deterministic", "stochastic"]:
            log_message(
                "{.arg fitting_by} must be one of {.val deterministic} and {.val stochastic}.",
                message_type="error",
            )
            exit()

        if not all([m in ["deterministic", "stochastic", "dynamical"] for m in mode]):
            log_message(
                "Invalid mode name! Must be one of {.val deterministic}, {.val stochastic} or {.val dynamical}.",
                message_type="error",
            )
            exit()

        adata.obs[group_by] = adata.obs[group_by].astype(dtype="category")

        log_message(
            "Using {.arg kernel_type}: {.val {kernel_type}}",
            message_type="info",
            verbose=verbose,
        )

        use_velocity = False
        use_pseudotime = False
        use_cytotrace = False

        if kernel_type == "velocity":
            has_velocity_data = (
                "spliced" in adata.layers and "unspliced" in adata.layers
            )

            if not has_velocity_data:
                log_message(
                    "No spliced/unspliced data found. Consider using {.arg kernel_type}={.val cytotrace} or {.arg kernel_type}={.val pseudotime} for RNA-only data",
                    message_type="warning",
                    verbose=verbose,
                )
                log_message(
                    "Falling back to {.pkg ConnectivityKernel} only",
                    message_type="info",
                    verbose=verbose,
                )
                use_velocity = False
            else:
                use_velocity = True

                if mode[-1] + "_graph" not in adata.obs.keys():
                    log_message(
                        "Running {.pkg scVelo} to compute RNA velocity...",
                        message_type="info",
                        verbose=verbose,
                    )
                    adata = SCVELO(
                        adata=adata,
                        group_by=group_by,
                        n_jobs=n_jobs,
                        linear_reduction=linear_reduction,
                        nonlinear_reduction=nonlinear_reduction,
                        basis=basis,
                        mode=mode,
                        fitting_by=fitting_by,
                        magic_impute=magic_impute,
                        knn=knn,
                        t=t,
                        min_shared_counts=min_shared_counts,
                        n_pcs=n_pcs,
                        n_neighbors=n_neighbors,
                        stream_smooth=stream_smooth,
                        stream_density=stream_density,
                        arrow_size=arrow_size,
                        arrow_length=arrow_length,
                        arrow_density=arrow_density,
                        denoise=denoise,
                        kinetics=kinetics,
                        calculate_velocity_genes=calculate_velocity_genes,
                        show_plot=show_plot,
                        save_plot=save_plot,
                        plot_format=plot_format,
                        plot_dpi=plot_dpi,
                        plot_prefix="scvelo",
                        dirpath=".",
                        legend_loc=legend_loc,
                        dpi=dpi,
                        save=save,
                        fileprefix=fileprefix,
                    )
                adata.layers["velocity"] = adata.layers[mode[-1]]

        elif kernel_type == "pseudotime":
            use_pseudotime = True
            if time_key not in adata.obs:
                log_message(
                    "Pseudotime {.val {time_key}} not found. Computing DPT pseudotime...",
                    message_type="info",
                    verbose=verbose,
                )
                if linear_reduction and linear_reduction in adata.obsm:
                    rep_key = linear_reduction
                else:
                    log_message(
                        "{.arg linear_reduction} {.val {linear_reduction}} not found in adata.obsm",
                        message_type="error",
                        verbose=verbose,
                    )
                    exit()

                if "connectivities" not in adata.obsp:
                    sc.pp.neighbors(
                        adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=rep_key
                    )

                sc.tl.diffmap(adata)
                adata.uns["iroot"] = 0
                sc.tl.dpt(adata)
                time_key = "dpt_pseudotime"
                adata.obs["cellrank_pseudotime"] = adata.obs["dpt_pseudotime"]
                log_message(
                    "{.pkg DPT} pseudotime computed and stored in {.val adata.obs['cellrank_pseudotime']}",
                    message_type="success",
                    verbose=verbose,
                )
            else:
                log_message(
                    "Using existing pseudotime from {.val adata.obs[{time_key}]}",
                    message_type="info",
                    verbose=verbose,
                )

        elif kernel_type == "cytotrace":
            use_cytotrace = True
            log_message(
                "Using {.pkg CytoTRACEKernel} for RNA-only data...",
                message_type="info",
                verbose=verbose,
            )

        log_message(
            "Using {.pkg CellRank} kernel-estimator architecture...",
            message_type="info",
            verbose=verbose,
        )

        main_kernel = None

        def get_rep_key():
            if linear_reduction and linear_reduction in adata.obsm:
                return linear_reduction
            else:
                log_message(
                    "{.arg linear_reduction} {.val {linear_reduction}} not found in adata.obsm",
                    message_type="error",
                    verbose=verbose,
                )
                exit()

        if use_pseudotime:
            log_message(
                "Creating {.pkg PseudotimeKernel} with {.arg time_key}={.val {time_key}}...",
                message_type="info",
                verbose=verbose,
            )
            try:
                pk = cr.kernels.PseudotimeKernel(adata, time_key=time_key)
                compute_transition_matrix(pk, verbose=verbose)
                main_kernel = pk
                log_message(
                    "{.pkg PseudotimeKernel} created successfully",
                    message_type="success",
                    verbose=verbose,
                )
            except Exception as e:
                log_message(
                    "{.pkg PseudotimeKernel} failed: {.val {e}}. Falling back to ConnectivityKernel.",
                    message_type="warning",
                    verbose=verbose,
                )
                use_pseudotime = False

        elif use_cytotrace:
            log_message(
                "Creating {.pkg CytoTRACEKernel} and computing {.pkg CytoTRACE} score...",
                message_type="info",
                verbose=verbose,
            )
            try:
                if "Ms" not in adata.layers:
                    log_message(
                        "Computing moments (required for {.pkg CytoTRACEKernel})...",
                        message_type="info",
                        verbose=verbose,
                    )
                    rep_key = get_rep_key()
                    if "connectivities" not in adata.obsp:
                        sc.pp.neighbors(
                            adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=rep_key
                        )

                    import numpy as np
                    from scipy.sparse import issparse, csr_matrix

                    log_message(
                        "Creating smoothed expression layer (Ms) for RNA-only data...",
                        message_type="info",
                        verbose=verbose,
                    )

                    connectivities = adata.obsp["connectivities"].copy()
                    if issparse(connectivities):
                        connectivities = connectivities.toarray()
                    connectivities += np.eye(connectivities.shape[0])
                    row_sums = connectivities.sum(axis=1, keepdims=True)
                    connectivities = connectivities / row_sums

                    X = adata.X.toarray() if issparse(adata.X) else adata.X
                    Ms = connectivities @ X
                    adata.layers["Ms"] = csr_matrix(Ms)

                    log_message(
                        "Smoothed expression layer (Ms) created successfully",
                        message_type="success",
                        verbose=verbose,
                    )

                ctk = cr.kernels.CytoTRACEKernel(adata)
                ctk.compute_cytotrace()
                compute_transition_matrix(ctk, verbose=verbose)
                main_kernel = ctk

                if "ct_score" in adata.obs:
                    adata.obs["cytotrace_score"] = adata.obs.pop("ct_score")
                    adata.obs["cytotrace_pseudotime"] = 1 - adata.obs["cytotrace_score"]
                    adata.obs["cellrank_pseudotime"] = adata.obs["cytotrace_pseudotime"]
                if "ct_pseudotime" in adata.obs:
                    del adata.obs["ct_pseudotime"]
                if "ct_num_exp_genes" in adata.obs:
                    adata.obs["cytotrace_num_exp_genes"] = adata.obs.pop(
                        "ct_num_exp_genes"
                    )

                log_message(
                    "{.pkg CytoTRACEKernel} created successfully",
                    message_type="success",
                    verbose=verbose,
                )
            except Exception as e:
                log_message(
                    "{.pkg CytoTRACEKernel} failed: {.val {e}}. Falling back to ConnectivityKernel.",
                    message_type="warning",
                    verbose=verbose,
                )
                use_cytotrace = False

        elif use_velocity:
            if velocity_weight <= 0 and connectivity_weight <= 0:
                log_message(
                    "Both kernel weights are <= 0. Using equal weights (0.5, 0.5).",
                    message_type="warning",
                    verbose=verbose,
                )
                velocity_weight = 0.5
                connectivity_weight = 0.5
            elif velocity_weight <= 0:
                log_message(
                    "{.arg velocity_weight <= 0}. Using {.pkg ConnectivityKernel} only",
                    message_type="info",
                    verbose=verbose,
                )
                use_velocity = False
                connectivity_weight = 1.0
            elif connectivity_weight <= 0:
                log_message(
                    "{.arg connectivity_weight <= 0}. Using {.pkg VelocityKernel} only.",
                    message_type="info",
                    verbose=verbose,
                )
                use_connectivity_kernel = False
                velocity_weight = 1.0
            else:
                total_weight = velocity_weight + connectivity_weight
                if abs(total_weight - 1.0) > 0.01:
                    log_message(
                        "Normalizing kernel weights ({.val {velocity_weight}}, {.val {connectivity_weight}}) to sum to 1.0",
                        message_type="info",
                        verbose=verbose,
                    )
                    velocity_weight = velocity_weight / total_weight
                    connectivity_weight = connectivity_weight / total_weight

        if use_velocity and velocity_weight > 0:
            log_message(
                "Creating {.pkg VelocityKernel} with {.arg model}={.val {mode[-1]}}, {.arg softmax_scale}={.val {softmax_scale}}...",
                message_type="info",
                verbose=verbose,
            )
            try:
                vk = cr.kernels.VelocityKernel(adata, backward=False)
                vk.compute_transition_matrix(
                    model=mode[-1], softmax_scale=softmax_scale
                )
                compute_transition_matrix(vk, verbose=verbose)
                main_kernel = vk
            except Exception as e:
                log_message(
                    "{.pkg VelocityKernel} computation failed: {.val {e}}",
                    message_type="error",
                    verbose=verbose,
                )
                log_message(
                    "Trying with default parameters...",
                    message_type="info",
                    verbose=verbose,
                )
                try:
                    vk = cr.kernels.VelocityKernel(adata, backward=False)
                    vk.compute_transition_matrix(softmax_scale=softmax_scale)
                    compute_transition_matrix(vk, verbose=verbose)
                    main_kernel = vk
                except Exception as e2:
                    log_message(
                        "{.pkg VelocityKernel} still failed: {.val {e2}}. Using {.pkg ConnectivityKernel} only.",
                        message_type="warning",
                        verbose=verbose,
                    )
                    use_velocity = False

        def ensure_neighbors():
            if "connectivities" not in adata.obsp:
                log_message(
                    "Neighbors not found. Computing neighbors...",
                    message_type="info",
                    verbose=verbose,
                )
                rep_key = get_rep_key()
                sc.pp.neighbors(
                    adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=rep_key
                )
                log_message(
                    "Neighbors computed successfully (n_neighbors={.val {n_neighbors}})",
                    message_type="success",
                    verbose=verbose,
                )

        if (
            main_kernel is not None
            and use_connectivity_kernel
            and connectivity_weight > 0
        ):
            try:
                log_message(
                    "Creating {.pkg ConnectivityKernel}...",
                    message_type="info",
                    verbose=verbose,
                )
                ensure_neighbors()
                ck = cr.kernels.ConnectivityKernel(adata)
                compute_transition_matrix(ck, verbose=verbose)
                final_kernel = velocity_weight * main_kernel + connectivity_weight * ck
                log_message(
                    "Combined kernels with weights: main={.val {format(velocity_weight, '.2f')}}, connectivity={.val {format(connectivity_weight, '.2f')}}",
                    message_type="success",
                    verbose=verbose,
                )
            except Exception as e:
                log_message(
                    "Failed to combine kernels: {.val {e}}. Using main kernel only.",
                    message_type="warning",
                    verbose=verbose,
                )
                final_kernel = main_kernel
        elif main_kernel is not None:
            final_kernel = main_kernel
            log_message(
                "Using {.val {type(main_kernel).__name__}} only",
                message_type="info",
                verbose=verbose,
            )
        else:
            log_message(
                "Creating {.pkg ConnectivityKernel} using neighbors by {.pkg compute_transition_matrix}...",
                message_type="info",
                verbose=verbose,
            )
            ensure_neighbors()
            ck = cr.kernels.ConnectivityKernel(adata)
            compute_transition_matrix(ck, verbose=verbose)
            final_kernel = ck
            log_message(
                "Using {.pkg ConnectivityKernel} only (connectivity-based transitions)",
                message_type="warning",
                verbose=verbose,
            )

        try:
            import numpy as np
            from scipy.sparse import diags, issparse, csr_matrix

            tmat = final_kernel.transition_matrix
            matrix_modified = False

            log_message(
                "Validating and fixing transition matrix for {.pkg GPCCA} compatibility...",
                message_type="info",
                verbose=verbose,
            )

            if issparse(tmat):
                if np.any(np.isnan(tmat.data)) or np.any(np.isinf(tmat.data)):
                    log_message(
                        "Cleaning NaN/Inf values...",
                        message_type="warning",
                        verbose=verbose,
                    )
                    tmat.data[np.isnan(tmat.data)] = 0
                    tmat.data[np.isinf(tmat.data)] = 0
                    tmat.eliminate_zeros()
                    matrix_modified = True
            else:
                if np.any(np.isnan(tmat)) or np.any(np.isinf(tmat)):
                    log_message(
                        "Cleaning NaN/Inf values...",
                        message_type="warning",
                        verbose=verbose,
                    )
                    tmat[np.isnan(tmat)] = 0
                    tmat[np.isinf(tmat)] = 0
                    matrix_modified = True

            if issparse(tmat):
                if np.any(tmat.data < 0):
                    log_message(
                        "Clipping {.val {(tmat.data < 0).sum()}} negative values to zero...",
                        message_type="warning",
                        verbose=verbose,
                    )
                    tmat.data = np.maximum(tmat.data, 0)
                    matrix_modified = True
            else:
                if np.any(tmat < 0):
                    log_message(
                        "Clipping negative values to zero...",
                        message_type="warning",
                        verbose=verbose,
                    )
                    tmat = np.maximum(tmat, 0)
                    matrix_modified = True

            if issparse(tmat):
                diag_values = tmat.diagonal()
                min_self_loop = 0.01
                needs_self_loop = diag_values < min_self_loop

                if np.any(needs_self_loop):
                    log_message(
                        "Adding self-loops to {.val {needs_self_loop.sum()}} cells for matrix primitivity...",
                        message_type="info",
                        verbose=verbose,
                    )

                    tmat = tmat.tolil()
                    for i in np.where(needs_self_loop)[0]:
                        tmat[i, i] = min_self_loop
                    tmat = tmat.tocsr()
                    matrix_modified = True

            row_sums = np.array(tmat.sum(axis=1)).flatten()
            zero_rows = row_sums < 1e-10

            if np.any(zero_rows):
                log_message(
                    "Found {.val {zero_rows.sum()}} zero rows. Setting to self-loops...",
                    message_type="warning",
                    verbose=verbose,
                )

                if issparse(tmat):
                    tmat = tmat.tolil()
                    for i in np.where(zero_rows)[0]:
                        tmat[i, i] = 1.0
                    tmat = tmat.tocsr()
                else:
                    for i in np.where(zero_rows)[0]:
                        tmat[i, :] = 0
                        tmat[i, i] = 1.0

                row_sums = np.array(tmat.sum(axis=1)).flatten()
                matrix_modified = True

            if not np.allclose(row_sums, 1.0, rtol=1e-6, atol=1e-8):
                log_message(
                    "Normalizing rows (current range: {.val {format(row_sums.min(), '.8f')}} - {.val {format(row_sums.max(), '.8f')}})...",
                    message_type="info",
                    verbose=verbose,
                )

                row_sums_safe = np.maximum(row_sums, 1e-10)

                if issparse(tmat):
                    tmat = diags(1.0 / row_sums_safe) @ tmat
                else:
                    tmat = tmat / row_sums_safe[:, np.newaxis]

                matrix_modified = True

            final_row_sums = np.array(tmat.sum(axis=1)).flatten()

            if issparse(tmat):
                has_negative = np.any(tmat.data < -1e-10)
            else:
                has_negative = np.any(tmat < -1e-10)

            rows_sum_to_one = np.allclose(final_row_sums, 1.0, rtol=1e-6, atol=1e-8)

            if has_negative:
                log_message(
                    "Warning: Matrix still has negative values after fixing",
                    message_type="warning",
                    verbose=verbose,
                )

            if not rows_sum_to_one:
                log_message(
                    "Warning: Matrix rows still don't sum to 1 (range: {.val {format(final_row_sums.min(), '.8f')}} - {.val {format(final_row_sums.max(), '.8f')}})",
                    message_type="warning",
                    verbose=verbose,
                )

            if matrix_modified:
                final_kernel._transition_matrix = tmat
                log_message(
                    "Matrix fixed and validated (row sums: {.val {format(final_row_sums.min(), '.8f')}} - {.val {format(final_row_sums.max(), '.8f')}})",
                    message_type="success",
                    verbose=verbose,
                )
            else:
                log_message(
                    "Matrix validation passed (row sums: {.val {format(final_row_sums.min(), '.8f')}} - {.val {format(final_row_sums.max(), '.8f')}})",
                    message_type="success",
                    verbose=verbose,
                )

        except Exception as e:
            log_message(
                "Matrix validation encountered error: {.val {e}}. Proceeding with original matrix...",
                message_type="warning",
                verbose=verbose,
            )

        log_message(
            "Creating {.val {estimator_type}} estimator...",
            message_type="info",
            verbose=verbose,
        )

        gpcca_failed = False

        if estimator_type == "GPCCA":
            try:
                estimator = cr.estimators.GPCCA(final_kernel)

                log_message(
                    "Computing eigendecomposition...",
                    message_type="info",
                    verbose=verbose,
                )
                estimator.compute_eigendecomposition()

                if n_macrostates is None:
                    n_cells = adata.n_obs
                    if n_cells < 100:
                        n_states_schur = 5
                    elif n_cells < 500:
                        n_states_schur = 8
                    elif n_cells < 2000:
                        n_states_schur = 10
                    else:
                        n_states_schur = 15

                    log_message(
                        "Auto-determined n_states={.val {n_states_schur}} for Schur (based on {.val {n_cells}} cells)",
                        message_type="info",
                        verbose=verbose,
                    )
                else:
                    n_states_schur = n_macrostates + 2
                    log_message(
                        "Using n_states={.val {n_states_schur}} for Schur (n_macrostates={.val {n_macrostates}} + 2)",
                        message_type="info",
                        verbose=verbose,
                    )

                log_message(
                    "Computing Schur decomposition (n_states={.val {n_states_schur}}, method={.val {schur_method}})...",
                    message_type="info",
                    verbose=verbose,
                )
                try:
                    estimator.compute_schur(n_states_schur, method=schur_method)
                except (ValueError, RuntimeError) as schur_error:
                    if "subspace_angles" in str(
                        schur_error
                    ) or "invariant subspace" in str(schur_error):
                        if schur_method == "brandts":
                            log_message(
                                "Schur decomposition with {.val {schur_method}} method failed: {.val {schur_error}}",
                                message_type="warning",
                                verbose=verbose,
                            )
                            log_message(
                                "Trying with {.val krylov} method instead...",
                                message_type="info",
                                verbose=verbose,
                            )
                            try:
                                estimator.compute_schur(n_states_schur, method="krylov")
                            except Exception as krylov_error:
                                log_message(
                                    "Schur decomposition with {.val krylov} method also failed: {.val {krylov_error}}",
                                    message_type="warning",
                                    verbose=verbose,
                                )
                                raise schur_error
                        else:
                            raise
                    else:
                        raise

                if n_macrostates is None:
                    n_macro = max(2, n_states_schur - 2)
                else:
                    n_macro = n_macrostates

                log_message(
                    "Computing macrostates (n_states={.val {n_macro}}, n_cells={.val {n_cells_terminal}})...",
                    message_type="info",
                    verbose=verbose,
                )
                try:
                    estimator.compute_macrostates(
                        n_states=n_macro, n_cells=n_cells_terminal
                    )
                except ValueError as macro_error:
                    if "Discretizing leads to a cluster" in str(
                        macro_error
                    ) or "0 samples" in str(macro_error):
                        log_message(
                            "Macrostates computation failed: {.val {macro_error}}",
                            message_type="warning",
                            verbose=verbose,
                        )
                        log_message(
                            "Trying with fewer macrostates...",
                            message_type="info",
                            verbose=verbose,
                        )
                        n_macro_reduced = max(2, n_macro - 2)
                        try:
                            estimator.compute_macrostates(
                                n_states=n_macro_reduced, n_cells=n_cells_terminal
                            )
                            log_message(
                                "Macrostates computed successfully with n_states={.val {n_macro_reduced}}",
                                message_type="success",
                                verbose=verbose,
                            )
                        except Exception as e2:
                            log_message(
                                "Reduced macrostates also failed: {.val {e2}}",
                                message_type="warning",
                                verbose=verbose,
                            )
                            log_message(
                                "Automatically switching to {.pkg CFLARE} estimator (more robust for problematic matrices)...",
                                message_type="warning",
                                verbose=verbose,
                            )
                            gpcca_failed = True
                            estimator_type = "CFLARE"
                            raise
                    else:
                        raise

                if not gpcca_failed:
                    log_message(
                        "Setting terminal states (n_cells={.val {n_cells_terminal}})...",
                        message_type="info",
                        verbose=verbose,
                    )
                    estimator.set_terminal_states(n_cells=n_cells_terminal)

                log_message(
                    "{.pkg GPCCA} estimator completed successfully",
                    message_type="success",
                    verbose=verbose,
                )

            except (ValueError, RuntimeError) as e:
                if (
                    "not a transition matrix" in str(e).lower()
                    or "schur" in str(e).lower()
                    or "gpcca" in str(e).lower()
                    or "Discretizing leads to a cluster" in str(e)
                    or "0 samples" in str(e)
                    or "subspace_angles" in str(e).lower()
                    or "invariant subspace" in str(e).lower()
                ):
                    if not gpcca_failed:
                        log_message(
                            "{.pkg GPCCA} estimator failed: {.val {e}}",
                            message_type="warning",
                            verbose=verbose,
                        )
                        log_message(
                            "Automatically switching to {.pkg CFLARE} estimator (more robust for problematic matrices)...",
                            message_type="warning",
                            verbose=verbose,
                        )
                        gpcca_failed = True
                        estimator_type = "CFLARE"
                    else:
                        # Already tried to switch, just raise
                        raise
                else:
                    log_message(
                        "{.pkg GPCCA} failed with unexpected error: {.val {e}}",
                        message_type="error",
                        verbose=verbose,
                    )
                    raise

        if estimator_type == "CFLARE":
            estimator = cr.estimators.CFLARE(final_kernel)

            log_message(
                "Computing eigendecomposition...", message_type="info", verbose=verbose
            )
            estimator.compute_eigendecomposition()

            predict_method = "leiden"
            log_message(
                "Predicting terminal states (use={.val {n_macrostates}}, method={.val {predict_method}})...",
                message_type="info",
                verbose=verbose,
            )

            try:
                estimator.predict(use=n_macrostates, method=predict_method)
            except Exception as e:
                log_message(
                    "Prediction with {.val {predict_method}} failed: {.val {e}}. Trying {.val kmeans}...",
                    message_type="warning",
                    verbose=verbose,
                )
                estimator.predict(use=n_macrostates, method="kmeans")

            if gpcca_failed:
                log_message(
                    "Successfully switched to {.pkg CFLARE} estimator",
                    message_type="success",
                    verbose=verbose,
                )

        log_message(
            "Computing fate probabilities...", message_type="info", verbose=verbose
        )
        estimator.compute_fate_probabilities()

        log_message(
            "Computing lineage drivers for {.arg cluster_key}={.val {group_by}}...",
            message_type="info",
            verbose=verbose,
        )
        try:
            estimator.compute_lineage_drivers(cluster_key=group_by, use_raw=False)
        except RuntimeError as e:
            if "Compute `.fate_probabilities`" in str(e):
                log_message(
                    "Skipping lineage drivers computation (fate probabilities not available)",
                    message_type="warning",
                    verbose=verbose,
                )
            else:
                raise

        log_message(
            "{.pkg CellRank} analysis completed",
            message_type="success",
            verbose=verbose,
        )

        if "cellrank_pseudotime" not in adata.obs:
            if "to_terminal_states" in adata.obsm or "lineages_fwd" in adata.obsm:
                try:
                    fate_probs_key = (
                        "to_terminal_states"
                        if "to_terminal_states" in adata.obsm
                        else "lineages_fwd"
                    )
                    fate_probs = adata.obsm[fate_probs_key]
                    if hasattr(fate_probs, "X"):
                        import numpy as np

                        max_probs = np.array(fate_probs.X.max(axis=1)).flatten()
                    else:
                        max_probs = np.array(fate_probs.max(axis=1)).flatten()
                    adata.obs["cellrank_pseudotime"] = max_probs
                    log_message(
                        "Created {.val cellrank_pseudotime} from fate probabilities",
                        message_type="info",
                        verbose=verbose,
                    )
                except Exception as e:
                    log_message(
                        "Failed to create {.val cellrank_pseudotime} from fate probabilities: {.val {e}}",
                        message_type="warning",
                        verbose=verbose,
                    )

        if save_plot:
            log_message(
                "Generating visualizations...", message_type="info", verbose=verbose
            )

            def get_save_path(name):
                ext = plot_format if plot_format in ["png", "pdf", "svg"] else "png"
                save_path = os.path.abspath(f"{plot_prefix}_{name}.{ext}")
                log_message(
                    "Save path for {.val {name}}: {.val {save_path}}",
                    message_type="info",
                    verbose=verbose,
                )
                return save_path

            try:
                log_message(
                    "Plotting eigenvalue spectrum...",
                    message_type="info",
                    verbose=verbose,
                )
                estimator.plot_spectrum(
                    real_only=True, dpi=plot_dpi, save=get_save_path("spectrum")
                )
                if show_plot:
                    plt.show()
                plt.close()
            except Exception as e:
                log_message(
                    "Failed to plot spectrum: {.val {e}}",
                    message_type="warning",
                    verbose=verbose,
                )

            if estimator_type == "GPCCA":
                try:
                    log_message(
                        "Plotting Schur matrix...",
                        message_type="info",
                        verbose=verbose,
                    )
                    estimator.plot_schur_matrix(
                        dpi=plot_dpi, save=get_save_path("schur_matrix")
                    )
                    if show_plot:
                        plt.show()
                    plt.close()
                except Exception as e:
                    log_message(
                        "Failed to plot Schur matrix: {.val {e}}",
                        message_type="warning",
                        verbose=verbose,
                    )

                try:
                    log_message(
                        "Plotting coarse-grained transition matrix...",
                        message_type="info",
                        verbose=verbose,
                    )
                    estimator.plot_coarse_T(
                        show_initial_dist=True,
                        show_stationary_dist=True,
                        dpi=plot_dpi,
                        save=get_save_path("coarse_T"),
                    )
                    if show_plot:
                        plt.show()
                    plt.close()
                except Exception as e:
                    log_message(
                        "Failed to plot coarse T: {.val {e}}",
                        message_type="warning",
                        verbose=verbose,
                    )

            try:
                log_message(
                    "Plotting all macrostates...",
                    message_type="info",
                    verbose=verbose,
                )
                estimator.plot_macrostates(
                    which="all",
                    basis=basis,
                    dpi=plot_dpi,
                    save=get_save_path("macrostates_all"),
                )
                if show_plot:
                    plt.show()
                plt.close()

                log_message(
                    "Plotting terminal states...",
                    message_type="info",
                    verbose=verbose,
                )
                estimator.plot_macrostates(
                    which="terminal",
                    basis=basis,
                    dpi=plot_dpi,
                    save=get_save_path("macrostates_terminal"),
                )
                if show_plot:
                    plt.show()
                plt.close()
            except Exception as e:
                log_message(
                    "Failed to plot macrostates: {.val {e}}",
                    message_type="warning",
                    verbose=verbose,
                )

            try:
                log_message(
                    "Plotting fate probabilities...",
                    message_type="info",
                    verbose=verbose,
                )
                estimator.plot_fate_probabilities(
                    basis=basis,
                    ncols=3,
                    dpi=plot_dpi,
                    save=get_save_path("fate_probabilities"),
                )
                if show_plot:
                    plt.show()
                plt.close()
            except Exception as e:
                log_message(
                    "Failed to plot fate probabilities: {.val {e}}",
                    message_type="warning",
                    verbose=verbose,
                )

            try:
                fate_probs_key = (
                    "to_terminal_states"
                    if not hasattr(estimator, "backward") or not estimator.backward
                    else "from_terminal_states"
                )
                if fate_probs_key in adata.obsm:
                    lineages = adata.obsm[fate_probs_key].names
                    log_message(
                        "Plotting lineage drivers for {.val {len(lineages)}} lineages...",
                        message_type="info",
                        verbose=verbose,
                    )

                    for lineage in lineages:
                        try:
                            estimator.plot_lineage_drivers(
                                lineage=lineage,
                                n_genes=8,
                                dpi=plot_dpi,
                                save=get_save_path(f"lineage_drivers_{lineage}"),
                            )
                            if show_plot:
                                plt.show()
                            plt.close()
                        except Exception as e:
                            log_message(
                                "Failed to plot lineage drivers for {.val {lineage}}: {.val {e}}",
                                message_type="warning",
                                verbose=verbose,
                            )
            except Exception as e:
                log_message(
                    "Failed to plot lineage drivers: {.val {e}}",
                    message_type="warning",
                    verbose=verbose,
                )

            try:
                log_message(
                    "Plotting aggregate fate probabilities ({.arg mode}={.val paga_pie})...",
                    message_type="info",
                    verbose=verbose,
                )
                cr.pl.aggregate_fate_probabilities(
                    adata,
                    mode="paga_pie",
                    cluster_key=group_by,
                    basis=basis,
                    dpi=plot_dpi,
                    save=get_save_path("aggregate_fates_paga_pie"),
                )
                if show_plot:
                    plt.show()
                plt.close()
            except Exception as e:
                log_message(
                    "Failed to plot aggregate fates: {.val {e}}",
                    message_type="warning",
                    verbose=verbose,
                )

            if "cellrank_pseudotime" in adata.obs:
                try:
                    log_message(
                        "Plotting gene expression trends...",
                        message_type="info",
                        verbose=verbose,
                    )

                    driver_key = "to_terminal_states_lineage_drivers"
                    if driver_key in adata.varm:
                        drivers_df = adata.varm[driver_key]
                        top_genes = drivers_df.nlargest(
                            20, columns=drivers_df.columns[0]
                        ).index.tolist()
                        log_message(
                            "Using top {.val {len(top_genes)}} driver genes",
                            message_type="info",
                            verbose=verbose,
                        )

                        model = cr.models.GAM()
                        cr.pl.gene_trends(
                            adata,
                            model=model,
                            genes=top_genes[:20],
                            time_key="cellrank_pseudotime",
                            n_jobs=n_jobs,
                            dpi=plot_dpi,
                            save=get_save_path("gene_trends"),
                        )
                        if show_plot:
                            plt.show()
                        plt.close()
                    else:
                        log_message(
                            "No lineage drivers found, skipping gene trends",
                            message_type="warning",
                            verbose=verbose,
                        )
                except Exception as e:
                    log_message(
                        "Failed to plot gene trends: {.val {e}}",
                        message_type="warning",
                        verbose=verbose,
                    )

            try:
                log_message(
                    "Plotting kernel projection...",
                    message_type="info",
                    verbose=verbose,
                )
                proj_basis = basis
                log_message(
                    "Using basis: {.val {proj_basis}} for projection",
                    message_type="info",
                    verbose=verbose,
                )
                if final_kernel is not None:
                    log_message(
                        "Final kernel type: {.val {type(final_kernel).__name__}}",
                        message_type="info",
                        verbose=verbose,
                    )

                    plot_palette = palette
                    if plot_palette is None and group_by is not None:
                        groups = (
                            adata.obs[group_by].cat.categories
                            if hasattr(adata.obs[group_by], "cat")
                            else adata.obs[group_by].unique()
                        )
                        plot_palette = dict(
                            zip(groups, plt.cm.tab10(np.linspace(0, 1, len(groups))))
                        )
                        log_message(
                            "Created default palette for {.val {len(groups)}} groups",
                            message_type="info",
                            verbose=verbose,
                        )

                    plot_kwargs = {}
                    if save_plot:
                        ext = (
                            plot_format
                            if plot_format in ["png", "pdf", "svg"]
                            else "png"
                        )
                        save_path = os.path.abspath(
                            f"{plot_prefix}_projection_{proj_basis}.{ext}"
                        )
                        plot_kwargs["save"] = save_path
                        log_message(
                            "Will save projection to: {.val {save_path}}",
                            message_type="info",
                            verbose=verbose,
                        )
                    if plot_dpi:
                        plot_kwargs["dpi"] = plot_dpi
                        log_message(
                            "Using DPI: {.val {plot_dpi}}",
                            message_type="info",
                            verbose=verbose,
                        )

                    if group_by is not None:
                        plot_kwargs["color"] = group_by
                        if plot_palette is not None:
                            plot_kwargs["palette"] = plot_palette
                        plot_kwargs["legend_loc"] = legend_loc
                        log_message(
                            "Adding color information: color={.val {group_by}}, palette={.val {'provided' if palette is not None else 'default'}}",
                            message_type="info",
                            verbose=verbose,
                        )

                    log_message(
                        "Calling {pkg final_kernel.plot_projection()}...",
                        message_type="info",
                        verbose=verbose,
                    )
                    try:
                        final_kernel.plot_projection(
                            basis=proj_basis,
                            recompute=True,
                            **plot_kwargs,
                        )
                        if show_plot:
                            import matplotlib.pyplot as plt

                            plt.show()
                        if save_plot and plot_format == "pdf":
                            save_path = plot_kwargs.get("save")
                            if save_path and save_path.endswith(".pdf"):
                                png_path = save_path.replace(".pdf", ".png")
                                if os.path.exists(png_path) and not os.path.exists(
                                    save_path
                                ):
                                    if os.path.exists(save_path):
                                        try:
                                            os.remove(save_path)
                                        except Exception:
                                            pass
                                    log_message(
                                        "PDF save failed, {.pkg CellRank} saved as PNG: {.val {png_path}}",
                                        message_type="warning",
                                        verbose=verbose,
                                    )
                        log_message(
                            "Projection plot completed successfully",
                            message_type="success",
                            verbose=verbose,
                        )
                    except Exception as e:
                        log_message(
                            "Failed to plot projection: {.val {e}}",
                            message_type="warning",
                            verbose=verbose,
                        )
                else:
                    log_message(
                        "Cannot plot projection: kernel object is None",
                        message_type="warning",
                        verbose=verbose,
                    )
            except Exception as e:
                log_message(
                    "Failed to plot projection: {.val {e}}",
                    message_type="warning",
                    verbose=verbose,
                )
                import traceback

                log_message(
                    "Projection error traceback: {.val {traceback.format_exc()}}",
                    message_type="warning",
                    verbose=verbose,
                )

        log_message("Visualization completed!", message_type="success", verbose=verbose)

        if use_velocity:
            try:
                if "dynamical" not in mode:
                    log_message(
                        "Running recover_dynamics for latent time computation...",
                        message_type="info",
                        verbose=verbose,
                    )
                    try:
                        vkey = mode[-1] if isinstance(mode, list) else mode
                        velocity_graph_key = vkey + "_graph"

                        if velocity_graph_key in adata.uns:
                            adata.uns["velocity_graph"] = adata.uns[velocity_graph_key]
                            if velocity_graph_key + "_neg" in adata.uns:
                                adata.uns["velocity_graph_neg"] = adata.uns[
                                    velocity_graph_key + "_neg"
                                ]
                            else:
                                log_message(
                                    "{.val velocity_graph_neg} not found, computing velocity graph...",
                                    message_type="info",
                                    verbose=verbose,
                                )
                                scv.tl.velocity_graph(
                                    adata,
                                    vkey=vkey,
                                    n_neighbors=n_neighbors,
                                    n_jobs=n_jobs,
                                )
                        elif "velocity_graph" not in adata.uns:
                            log_message(
                                "{.val velocity_graph} not found, computing velocity graph...",
                                message_type="info",
                                verbose=verbose,
                            )
                            scv.tl.velocity_graph(
                                adata, vkey=vkey, n_neighbors=n_neighbors, n_jobs=n_jobs
                            )

                        scv.tl.recover_dynamics(adata, use_raw=False, n_jobs=n_jobs)
                        terminal_key = (
                            "to_terminal_states_probs"
                            if hasattr(estimator, "terminal_states_probabilities")
                            else "terminal_states_probs"
                        )
                        initial_key = (
                            "to_initial_states_probs"
                            if hasattr(estimator, "initial_states_probabilities")
                            else "initial_states_probs"
                        )

                        scv.tl.recover_latent_time(
                            adata,
                            root_key=initial_key,
                            end_key=terminal_key,
                        )
                        if "latent_time" in adata.obs:
                            adata.obs["cellrank_pseudotime"] = adata.obs["latent_time"]
                        log_message(
                            "Latent time computed successfully",
                            message_type="success",
                            verbose=verbose,
                        )
                    except Exception as e:
                        log_message(
                            "{.pkg recover_dynamics} failed ({.val {e}}), skipping latent time computation...",
                            message_type="warning",
                            verbose=verbose,
                        )
                else:
                    terminal_key = (
                        "to_terminal_states_probs"
                        if hasattr(estimator, "terminal_states_probabilities")
                        else "terminal_states_probs"
                    )
                    initial_key = (
                        "to_initial_states_probs"
                        if hasattr(estimator, "initial_states_probabilities")
                        else "initial_states_probs"
                    )

                    scv.tl.recover_latent_time(
                        adata, root_key=initial_key, end_key=terminal_key
                    )
                    if "latent_time" in adata.obs:
                        adata.obs["cellrank_pseudotime"] = adata.obs["latent_time"]
                    log_message(
                        "Latent time computed successfully",
                        message_type="success",
                        verbose=verbose,
                    )
            except Exception as e:
                log_message(
                    "Latent time computation failed: {.val {e}}",
                    message_type="warning",
                    verbose=verbose,
                )
        else:
            log_message(
                "Skipping latent time computation (no velocity data available)",
                message_type="info",
                verbose=verbose,
            )

    finally:
        try:
            figures_dir = os.path.join(os.getcwd(), "figures")
            if os.path.exists(figures_dir) and os.path.isdir(figures_dir):
                if not os.listdir(figures_dir):
                    os.rmdir(figures_dir)
                    log_message(
                        "Removed empty figures directory: {.val {figures_dir}}",
                        message_type="info",
                        verbose=verbose,
                    )
        except Exception:
            pass

        if save_plot:
            os.chdir(prevdir)

    try:
        adata.__dict__["_raw"].__dict__["_var"] = (
            adata.__dict__["_raw"]
            .__dict__["_var"]
            .rename(columns={"_index": "features"})
        )
    except:
        pass

    log_message(
        "Returning adata, estimator, and kernel objects",
        message_type="info",
        verbose=verbose,
    )
    return adata, estimator, final_kernel


def PAGA(
    adata=None,
    h5ad=None,
    group_by=None,
    palette=None,
    linear_reduction=None,
    nonlinear_reduction=None,
    basis=None,
    n_pcs=30,
    n_neighbors=30,
    use_rna_velocity=False,
    vkey="stochastic",
    embedded_with_PAGA=False,
    paga_layout="fr",
    threshold=0.1,
    point_size=20,
    infer_pseudotime=False,
    root_cell=None,
    root_group=None,
    n_dcs=10,
    n_branchings=0,
    min_group_size=0.01,
    n_jobs=1,
    show_plot=True,
    save_plot=False,
    plot_format="png",
    plot_dpi=600,
    plot_prefix="paga",
    dirpath="./paga",
    dpi=300,
    save=False,
    fileprefix="",
    verbose=True,
    legend_loc="on data",
):
    import os
    import platform
    import sys

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["KMP_WARNINGS"] = "0"
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

    is_apple_silicon = platform.system() == "Darwin" and platform.machine() == "arm64"
    if is_apple_silicon:
        configure_apple_silicon_env(
            scanpy_settings=True,
            scanpy_verbosity=True,
            numba_threading=True,
            numba_disable_jit=True,
            configure_numba_runtime=True,
            numba_runtime_optional=True,
            verbose=verbose,
        )

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        log_message(
            "{.pkg matplotlib} setup failed: {.val {e}}",
            message_type="warning",
            verbose=verbose,
        )
        plt = None

    try:
        import scanpy as sc
    except Exception as e:
        err = str(e)
        if "___kmpc_dispatch_deinit" in err or "libomp.dylib" in err:
            log_message(
                "Detected {.pkg OpenMP} runtime conflict while importing {.pkg scanpy}/{.pkg python-igraph} on macOS.",
                message_type="warning",
                verbose=verbose,
            )
            log_message(
                "Please recreate/update environment with {.pkg python-igraph} from pip and ensure {.pkg libomp} is installed in conda env.",
                message_type="warning",
                verbose=verbose,
            )
        log_message("{.pkg scanpy} import failed: {.val {e}}", message_type="error")
        raise

    try:
        import numpy as np
    except Exception as e:
        log_message("{.pkg numpy} import failed: {.val {e}}", message_type="error")
        raise

    import statistics
    from math import hypot
    import warnings

    warnings.simplefilter("ignore", category=UserWarning)
    warnings.simplefilter("ignore", category=FutureWarning)
    warnings.simplefilter("ignore", category=DeprecationWarning)

    prevdir = os.getcwd()

    expanded_path = os.path.expanduser(dirpath)
    os.makedirs(expanded_path, exist_ok=True)
    os.chdir(expanded_path)

    import scanpy as sc

    sc.settings.figdir = "."

    import platform

    if platform.system() == "Windows":
        import sys, multiprocessing, re

        if re.match(pattern=".*pythonw.exe$", string=sys.executable):
            pythonw = sys.executable
        else:
            pythonw = sys.executable.replace("python.exe", "pythonw.exe")
        sys.executable = pythonw
        sys._base_executable = pythonw
        multiprocessing.set_executable(pythonw)

    try:
        if adata is None and h5ad is None:
            log_message("adata or h5ad must be provided", message_type="error")
            return None

        if adata is None:
            adata = sc.read(h5ad)

        if group_by is None:
            log_message("group_by must be provided", message_type="error")
            return None

        if linear_reduction is None and nonlinear_reduction is None:
            log_message(
                "linear_reduction or nonlinear_reduction must be provided at least one",
                message_type="error",
            )
            return None

        if basis is None:
            if nonlinear_reduction is not None:
                basis = nonlinear_reduction
            else:
                basis = "basis"
                adata.obsm["basis"] = adata.obsm[linear_reduction][:, 0:2]

        if point_size is None:
            point_size = min(100000 / adata.shape[0], 20)

        if infer_pseudotime is True and root_cell is None and root_group is None:
            log_message(
                "root_cell or root_group should be provided", message_type="error"
            )
            return None

        if use_rna_velocity is True:
            adata.uns["velocity_graph"] = adata.uns[vkey + "_graph"]

        adata.obs[group_by] = adata.obs[group_by].astype(dtype="category")

        rep_key = linear_reduction
        try:
            obsm_keys = set(adata.obsm_keys())
            if linear_reduction not in obsm_keys:
                log_message(
                    "{.arg linear_reduction} {.val {linear_reduction}} not found in adata.obsm",
                    message_type="error",
                    verbose=verbose,
                )
                return None
            rep_key = linear_reduction

            neighbors_n_pcs = n_pcs
            if rep_key in obsm_keys:
                rep_dims = adata.obsm[rep_key].shape[1]
                if neighbors_n_pcs is not None and isinstance(
                    neighbors_n_pcs, (int, float)
                ):
                    if rep_dims < int(neighbors_n_pcs):
                        neighbors_n_pcs = None
            else:
                neighbors_n_pcs = n_pcs
        except Exception:
            neighbors_n_pcs = n_pcs

        if "X_diffmap" in adata.obsm_keys():
            X_diffmap = adata.obsm["X_diffmap"]
            del adata.obsm["X_diffmap"]
            sc.pp.neighbors(
                adata, n_pcs=neighbors_n_pcs, use_rep=rep_key, n_neighbors=n_neighbors
            )
            adata.obsm["X_diffmap"] = X_diffmap
        else:
            sc.pp.neighbors(
                adata, n_pcs=neighbors_n_pcs, use_rep=rep_key, n_neighbors=n_neighbors
            )

        sc.tl.paga(adata, groups=group_by, use_rna_velocity=use_rna_velocity)

        try:
            if use_rna_velocity is True:
                sc.pl.paga_compare(
                    adata,
                    basis=basis,
                    palette=palette,
                    threshold=threshold,
                    size=point_size,
                    min_edge_width=1,
                    node_size_scale=1,
                    dashed_edges="connectivities",
                    solid_edges="transitions_confidence",
                    transitions="transitions_confidence",
                    title=basis,
                    frameon=False,
                    edges=True,
                    save=False,
                    show=show_plot,
                )
            else:
                sc.pl.paga_compare(
                    adata,
                    basis=basis,
                    palette=palette,
                    threshold=threshold,
                    size=point_size,
                    title=basis,
                    frameon=False,
                    edges=True,
                    save=False,
                    show=show_plot,
                )

            if save:
                plt.savefig(
                    "./" + ".".join(filter(None, [fileprefix, "paga_compare.pdf"])),
                    dpi=dpi,
                    bbox_inches="tight",
                    facecolor="white",
                )
        except Exception as e:
            log_message(
                "{.pkg PAGA} compare plot failed ({.val {e}}), continuing...",
                message_type="warning",
            )
            if plt is not None:
                plt.clf()
                plt.close("all")

        try:
            sc.pl.paga(
                adata,
                threshold=threshold,
                layout=paga_layout,
                title="PAGA layout: " + paga_layout,
                frameon=False,
                save=False,
                show=show_plot,
            )

            if save:
                plt.savefig(
                    "./" + ".".join(filter(None, [fileprefix, "paga_layout.pdf"])),
                    dpi=dpi,
                    bbox_inches="tight",
                    facecolor="white",
                )
        except Exception as e:
            log_message(
                "{.pkg PAGA} layout plot failed ({.val {e}}), continuing...",
                message_type="warning",
            )
            if plt is not None:
                plt.clf()
                plt.close("all")

        try:
            sc.tl.draw_graph(adata, init_pos="paga", layout=paga_layout)
            log_message(
                "draw_graph computation completed successfully", message_type="success"
            )

            sc.pl.draw_graph(
                adata,
                color=group_by,
                palette=palette,
                title="PAGA layout: " + paga_layout,
                layout=paga_layout,
                frameon=False,
                legend_loc=legend_loc,
                show=show_plot,
            )

            if save:
                plt.savefig(
                    "./" + ".".join(filter(None, [fileprefix, "paga_graph.pdf"])),
                    dpi=dpi,
                    bbox_inches="tight",
                    facecolor="white",
                )
        except Exception as e:
            log_message(
                "{.pkg PAGA} draw_graph failed ({.val {e}}), continuing...",
                message_type="warning",
            )
            if plt is not None:
                plt.clf()
                plt.close("all")

        if embedded_with_PAGA is True:
            try:
                umap2d = sc.tl.umap(adata, init_pos="paga", n_components=2, copy=True)
                adata.obsm["PAGAUMAP2D"] = umap2d.obsm["X_umap"]

                sc.pl.paga_compare(
                    adata,
                    basis="PAGAUMAP2D",
                    palette=palette,
                    threshold=threshold,
                    size=point_size,
                    title="PAGA-initialized UMAP",
                    edges=True,
                    save=False,
                    show=show_plot,
                )

                if save:
                    plt.savefig(
                        "./" + ".".join(filter(None, [fileprefix, "paga_umap.pdf"])),
                        dpi=dpi,
                        bbox_inches="tight",
                        facecolor="white",
                    )
            except Exception as e:
                log_message(
                    "{.pkg PAGA} UMAP failed ({.val {e}}), continuing...",
                    message_type="warning",
                )
                if plt is not None:
                    plt.clf()
                    plt.close("all")

        if infer_pseudotime is True:
            if root_group is not None and root_cell is None:
                cell = adata.obs[group_by].index.values[
                    adata.obs[group_by] == root_group
                ]
                root_group_cell = adata.obsm[basis][adata.obs[group_by] == root_group,][
                    :, [0, 1]
                ]
                x = statistics.median(root_group_cell[:, 0])
                y = statistics.median(root_group_cell[:, 1])
                diff = np.array((x - root_group_cell[:, 0], y - root_group_cell[:, 1]))
                dist = []
                for i in range(diff.shape[1]):
                    dist.append(hypot(diff[0, i], diff[1, i]))

                root_cell = cell[dist.index(min(dist))]

            sc.tl.diffmap(adata, n_comps=n_dcs)
            adata.uns["iroot"] = np.flatnonzero(adata.obs_names == root_cell)[0]
            sc.tl.dpt(
                adata,
                n_dcs=n_dcs,
                n_branchings=n_branchings,
                min_group_size=min_group_size,
            )

            try:
                sc.pl.embedding(
                    adata,
                    basis=basis,
                    color="dpt_pseudotime",
                    save=False,
                    show=show_plot,
                )

                if save:
                    plt.savefig(
                        "./"
                        + ".".join(filter(None, [fileprefix, "dpt_pseudotime.pdf"])),
                        dpi=dpi,
                        bbox_inches="tight",
                        facecolor="white",
                    )
            except Exception as e:
                log_message(
                    "DPT pseudotime plot failed ({.val {e}}), continuing...",
                    message_type="warning",
                )
                if plt is not None:
                    plt.clf()
                    plt.close("all")

    finally:
        try:
            figures_dir = os.path.join(os.getcwd(), "figures")
            if os.path.exists(figures_dir) and os.path.isdir(figures_dir):
                if not os.listdir(figures_dir):
                    os.rmdir(figures_dir)
                    log_message(
                        "Removed empty figures directory: {.val {figures_dir}}",
                        message_type="info",
                        verbose=verbose,
                    )
        except Exception:
            pass

        os.chdir(prevdir)

        if is_apple_silicon and plt is not None:
            try:
                plt.clf()
                plt.close("all")
                import gc

                gc.collect()
            except:
                pass

    try:
        adata.__dict__["_raw"].__dict__["_var"] = (
            adata.__dict__["_raw"]
            .__dict__["_var"]
            .rename(columns={"_index": "features"})
        )
    except:
        pass

    return adata


def Palantir(
    adata=None,
    h5ad=None,
    group_by=None,
    palette=None,
    linear_reduction=None,
    nonlinear_reduction=None,
    basis=None,
    n_pcs=30,
    n_neighbors=30,
    dm_n_components=10,
    dm_alpha=0,
    dm_n_eigs=None,
    early_group=None,
    terminal_groups=None,
    early_cell=None,
    terminal_cells=None,
    num_waypoints=1200,
    scale_components=True,
    use_early_cell_as_start=False,
    adjust_early_cell=False,
    adjust_terminal_cells=False,
    max_iterations=25,
    n_jobs=1,
    point_size=20,
    show_plot=True,
    dpi=300,
    save=False,
    dirpath="./",
    fileprefix="",
    verbose=True,
    legend_loc="on data",
):
    import os
    import platform

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["KMP_WARNINGS"] = "0"
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

    is_apple_silicon = platform.system() == "Darwin" and platform.machine() == "arm64"

    if is_apple_silicon:
        configure_apple_silicon_env(
            scanpy_settings=True,
            configure_numba_runtime=True,
            verbose=verbose,
        )

    import matplotlib.pyplot as plt
    import matplotlib
    import statistics
    from math import hypot
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import palantir

    import warnings

    warnings.simplefilter("ignore", category=UserWarning)
    warnings.simplefilter("ignore", category=FutureWarning)
    warnings.simplefilter("ignore", category=DeprecationWarning)

    prevdir = os.getcwd()
    expanded_path = os.path.expanduser(dirpath)
    os.makedirs(expanded_path, exist_ok=True)
    os.chdir(expanded_path)

    if platform.system() == "Windows":
        import sys, multiprocessing, re

        if re.match(pattern=".*pythonw.exe$", string=sys.executable):
            pythonw = sys.executable
        else:
            pythonw = sys.executable.replace("python.exe", "pythonw.exe")
        sys.executable = pythonw
        sys._base_executable = pythonw
        multiprocessing.set_executable(pythonw)

    try:
        if adata is None and h5ad is None:
            log_message("adata or h5ad must be provided", message_type="error")
            exit()

        if adata is None:
            adata = sc.read(h5ad)

        if group_by is None and (
            early_group is not None or terminal_groups is not None
        ):
            log_message("{.arg group_by} must be provided", message_type="error")
            exit()

        if linear_reduction is None and nonlinear_reduction is None:
            log_message(
                "{.arg linear_reduction} or {.arg nonlinear_reduction} must be provided at least one",
                message_type="error",
            )
            exit()

        # if linear_reduction is None:
        #     sc.pp.pca(adata, n_comps=n_pcs)
        #     linear_reduction = "X_pca"

        if basis is None:
            if nonlinear_reduction is not None:
                basis = nonlinear_reduction
            else:
                basis = linear_reduction

        if point_size is None:
            point_size = min(100000 / adata.shape[0], 20)

        if early_group is not None and early_cell is None:
            cell = adata.obs[group_by].index.values[adata.obs[group_by] == early_group]
            early_group_cell = adata.obsm[basis][adata.obs[group_by] == early_group,][
                :, [0, 1]
            ]
            x = statistics.median(early_group_cell[:, 0])
            y = statistics.median(early_group_cell[:, 1])
            diff = np.array((x - early_group_cell[:, 0], y - early_group_cell[:, 1]))
            dist = []
            for i in range(diff.shape[1]):
                dist.append(hypot(diff[0, i], diff[1, i]))

            early_cell = cell[dist.index(min(dist))]

        if early_cell is None:
            log_message("{.arg early_cell} must be provided", message_type="error")
            exit()
        else:
            log_message(
                "{.arg early_cell}: {.val {early_cell}}",
                message_type="info",
                verbose=verbose,
            )

        terminal_cells_dict = dict()
        if terminal_groups is not None and terminal_cells is None:
            for n in range(len(terminal_groups)):
                terminal_group = terminal_groups[n]
                cell = adata.obs[group_by].index.values[
                    adata.obs[group_by] == terminal_group
                ]
                terminal_group_cell = adata.obsm[basis][
                    adata.obs[group_by] == terminal_group,
                ][:, [0, 1]]
                x = statistics.median(terminal_group_cell[:, 0])
                y = statistics.median(terminal_group_cell[:, 1])
                diff = np.array(
                    (x - terminal_group_cell[:, 0], y - terminal_group_cell[:, 1])
                )
                dist = []
                for i in range(diff.shape[1]):
                    dist.append(hypot(diff[0, i], diff[1, i]))
                terminal_cells_dict[cell[dist.index(min(dist))]] = (
                    terminal_group.replace(" ", ".") + "_diff_potential"
                )

            terminal_cells = list(terminal_cells_dict.keys())

        if terminal_cells is None:
            log_message(
                "{.arg terminal_cells}: None", message_type="info", verbose=verbose
            )
        else:
            log_message(
                "{.arg terminal_cells}: {.val {terminal_cells}}",
                message_type="info",
                verbose=verbose,
            )

        pca_projections = pd.DataFrame(
            adata.obsm[linear_reduction][:, :n_pcs], index=adata.obs_names
        )
        log_message("Running diffusion maps", message_type="info", verbose=verbose)
        dm_res = palantir.utils.run_diffusion_maps(
            pca_projections,
            n_components=dm_n_components,
            knn=n_neighbors,
            alpha=dm_alpha,
        )
        ms_data = palantir.utils.determine_multiscale_space(dm_res, n_eigs=dm_n_eigs)
        log_message("Running palantir", message_type="info", verbose=verbose)
        pr_res = palantir.core.run_palantir(
            data=ms_data,
            early_cell=early_cell,
            terminal_states=terminal_cells,
            knn=n_neighbors,
            num_waypoints=num_waypoints,
            scale_components=scale_components,
            use_early_cell_as_start=use_early_cell_as_start,
            max_iterations=max_iterations,
            n_jobs=n_jobs,
        )

        if adjust_early_cell is True or adjust_terminal_cells is True:
            if adjust_early_cell is True:
                early_cell_group = adata.obs[group_by][early_cell]
                cells = adata.obs[group_by].index.values[
                    adata.obs[group_by] == early_cell_group
                ]
                early_cell = pr_res.pseudotime[cells].index.values[
                    pr_res.pseudotime[cells] == min(pr_res.pseudotime[cells])
                ][0]
            if adjust_terminal_cells is True:
                terminal_cells_dict = dict()
                for n in range(len(terminal_cells)):
                    terminal_cell = terminal_cells[n]
                    terminal_cell_group = adata.obs[group_by][terminal_cell]
                    cells = adata.obs[group_by].index.values[
                        adata.obs[group_by] == terminal_cell_group
                    ]
                    terminal_cells_dict[
                        pr_res.pseudotime[cells].index.values[
                            pr_res.pseudotime[cells] == max(pr_res.pseudotime[cells])
                        ][0]
                    ] = terminal_cell_group.replace(" ", ".") + "_diff_potential"
                terminal_cells = list(terminal_cells_dict.keys())

            pr_res = palantir.core.run_palantir(
                data=ms_data,
                early_cell=early_cell,
                terminal_states=terminal_cells,
                knn=n_neighbors,
                num_waypoints=num_waypoints,
                scale_components=scale_components,
                use_early_cell_as_start=use_early_cell_as_start,
                max_iterations=max_iterations,
                n_jobs=n_jobs,
            )

        adata.obsm["palantir_dm"] = dm_res["T"].toarray()
        adata.uns["dm_kernel"] = dm_res["kernel"]
        if len(terminal_cells_dict) > 0:
            pr_res.branch_probs = pr_res.branch_probs.rename(
                columns=terminal_cells_dict
            )
        for term in np.append(
            pr_res.branch_probs.columns.values,
            np.array(["palantir_pseudotime", "palantir_diff_potential"]),
        ):
            if term in adata.obs.columns:
                adata.obs.drop(term, axis=1, inplace=True)
        adata.obs = adata.obs.join(pr_res.pseudotime.to_frame("palantir_pseudotime"))
        adata.obs = adata.obs.join(pr_res.entropy.to_frame("palantir_diff_potential"))
        adata.obs = adata.obs.join(pr_res.branch_probs)

        sc.pl.embedding(
            adata,
            basis=basis,
            color="palantir_pseudotime",
            size=point_size,
            show=show_plot,
        )
        if save:
            plt.savefig(
                "./" + ".".join(filter(None, [fileprefix, "palantir_pseudotime.pdf"])),
                dpi=dpi,
            )

        sc.pl.embedding(
            adata,
            basis=basis,
            color="palantir_diff_potential",
            size=point_size,
            show=show_plot,
        )
        if save:
            plt.savefig(
                "./"
                + ".".join(filter(None, [fileprefix, "palantir_diff_potential.pdf"])),
                dpi=dpi,
            )

        sc.pl.embedding(
            adata,
            basis=basis,
            color=pr_res.branch_probs.columns.values,
            size=point_size,
            show=show_plot,
        )
        if save:
            plt.savefig(
                "./" + ".".join(filter(None, [fileprefix, "palantir_probs.pdf"])),
                dpi=dpi,
            )

        if group_by is not None:
            sc.pl.embedding(
                adata,
                basis=basis,
                color=group_by,
                size=point_size,
                palette=palette,
                legend_loc=legend_loc,
                show=show_plot,
            )
            if save:
                plt.savefig(
                    "./"
                    + ".".join(filter(None, [fileprefix, "palantir_group_by.pdf"])),
                    dpi=dpi,
                )

    finally:
        try:
            figures_dir = os.path.join(os.getcwd(), "figures")
            if os.path.exists(figures_dir) and os.path.isdir(figures_dir):
                if not os.listdir(figures_dir):
                    os.rmdir(figures_dir)
                    log_message(
                        "Removed empty figures directory: {.val {figures_dir}}",
                        message_type="info",
                        verbose=verbose,
                    )
        except Exception:
            pass

        os.chdir(prevdir)

    try:
        adata.__dict__["_raw"].__dict__["_var"] = (
            adata.__dict__["_raw"]
            .__dict__["_var"]
            .rename(columns={"_index": "features"})
        )
    except:
        pass

    return adata


def WOT(
    adata=None,
    h5ad=None,
    group_by=None,
    palette=None,
    time_field="Time",
    growth_iters=3,
    tmap_out="tmaps/tmap_out",
    time_from=None,
    time_to=None,
    get_coupling=False,
    recalculate=False,
    show_plot=True,
    dpi=300,
    save=False,
    dirpath="./",
    fileprefix="",
    verbose=True,
):
    import os
    import platform

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["KMP_WARNINGS"] = "0"
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

    is_apple_silicon = platform.system() == "Darwin" and platform.machine() == "arm64"

    if is_apple_silicon:
        configure_apple_silicon_env(
            scanpy_settings=True,
            configure_numba_runtime=True,
            verbose=verbose,
        )

    import matplotlib.pyplot as plt
    import scanpy as sc
    import numpy as np
    import statistics
    import pandas as pd
    from math import hypot
    import wot

    import warnings

    warnings.simplefilter("ignore", category=UserWarning)
    warnings.simplefilter("ignore", category=FutureWarning)
    warnings.simplefilter("ignore", category=DeprecationWarning)

    import os

    prevdir = os.getcwd()
    expanded_path = os.path.expanduser(dirpath)
    os.makedirs(expanded_path, exist_ok=True)
    os.chdir(expanded_path)

    import platform

    if platform.system() == "Windows":
        import sys
        import re
        import multiprocessing

        if re.match(pattern=".*pythonw.exe$", string=sys.executable):
            pythonw = sys.executable
        else:
            pythonw = sys.executable.replace("python.exe", "pythonw.exe")
        sys.executable = pythonw
        sys._base_executable = pythonw
        multiprocessing.set_executable(pythonw)

    try:
        if adata is None and h5ad is None:
            log_message(
                "{.arg adata} or {.arg h5ad} must be provided", message_type="error"
            )
            exit()

        if adata is None:
            adata = sc.read(h5ad)

        if group_by is None:
            log_message("{.arg group_by} must be provided", message_type="error")
            exit()

        if time_field is None:
            log_message("{.arg time_field} must be provided", message_type="error")
            exit()

        adata.obs[group_by] = adata.obs[group_by].astype(dtype="category")
        if pd.api.types.is_categorical_dtype(adata.obs[time_field]):
            adata.obs["time_field"] = adata.obs[time_field].cat.codes
        elif not pd.api.types.is_numeric_dtype(adata.obs[time_field]):
            try:
                adata.obs["time_field"] = adata.obs[time_field].astype("float")
            except ValueError:
                log_message(
                    "Unable to convert column {.val {time_field}} to float type",
                    message_type="warning",
                )
        else:
            adata.obs["time_field"] = adata.obs[time_field]

        time_dict = dict(zip(adata.obs[time_field], adata.obs["time_field"]))
        if time_from not in time_dict.keys():
            log_message("{.arg time_from} is incorrect", message_type="error")
            exit()

        ot_model = wot.ot.OTModel(
            adata, growth_iters=growth_iters, day_field="time_field"
        )

        if recalculate is True:
            ot_model.compute_all_transport_maps(tmap_out=tmap_out)
            tmap_model = wot.tmap.TransportMapModel.from_directory(tmap_out)
        else:
            try:
                tmap_model = wot.tmap.TransportMapModel.from_directory(tmap_out)
            except (FileNotFoundError, ValueError):
                ot_model.compute_all_transport_maps(tmap_out=tmap_out)
                tmap_model = wot.tmap.TransportMapModel.from_directory(tmap_out)

        cell_sets = {}
        for k, v in zip(adata.obs[group_by], adata.obs_names):
            if k not in cell_sets:
                cell_sets[k] = []
            cell_sets[k].append(v)

        from_populations = tmap_model.population_from_cell_sets(
            cell_sets, at_time=time_dict[time_from]
        )

        trajectory_ds = tmap_model.trajectories(from_populations)
        trajectory_df = pd.DataFrame(
            trajectory_ds.X,
            index=trajectory_ds.obs_names,
            columns=trajectory_ds.var_names,
        )
        adata.uns["trajectory_" + str(time_from)] = trajectory_df.reindex(
            adata.obs_names
        )

        fates_ds = tmap_model.fates(from_populations)
        fates_df = pd.DataFrame(
            fates_ds.X, index=fates_ds.obs_names, columns=fates_ds.var_names
        )
        existing_rows = fates_df.index.tolist()
        new_rows = list(set(adata.obs_names) - set(existing_rows))
        new_df = pd.DataFrame(0, index=new_rows, columns=fates_df.columns)
        fates_df = pd.concat([fates_df, new_df])
        adata.uns["fates_" + str(time_from)] = fates_df.reindex(adata.obs_names)

        if time_to is not None:
            if time_to not in time_dict.keys():
                log_message("{.arg time_to} is incorrect", message_type="error")
                exit()

            to_populations = tmap_model.population_from_cell_sets(
                cell_sets, at_time=time_dict[time_to]
            )
            transition_table = tmap_model.transition_table(
                from_populations, to_populations
            )
            transition_df = pd.DataFrame(
                transition_table.X,
                index=transition_table.obs_names,
                columns=transition_table.var_names,
            )
            adata.uns["transition_" + str(time_from) + "_to_" + str(time_to)] = fates_df
            if get_coupling is True:
                coupling = tmap_model.get_coupling(
                    time_dict[time_from], time_dict[time_to]
                )
                coupling_df = pd.DataFrame(
                    coupling.X, index=coupling.obs_names, columns=coupling.var_names
                )
                adata.uns["coupling_" + str(time_from) + "_to_" + str(time_to)] = (
                    coupling_df
                )

    finally:
        try:
            figures_dir = os.path.join(os.getcwd(), "figures")
            if os.path.exists(figures_dir) and os.path.isdir(figures_dir):
                if not os.listdir(figures_dir):
                    os.rmdir(figures_dir)
                    log_message(
                        "Removed empty figures directory: {.val {figures_dir}}",
                        message_type="info",
                        verbose=verbose,
                    )
        except Exception:
            pass

        os.chdir(prevdir)

    try:
        adata.__dict__["_raw"].__dict__["_var"] = (
            adata.__dict__["_raw"]
            .__dict__["_var"]
            .rename(columns={"_index": "features"})
        )
    except:
        pass

    return adata


def CellTypistModels(on_the_fly=False, verbose=True):
    """
    Get a list of all available CellTypist models.

    Parameters
    ----------
    on_the_fly : bool
        If True, only show downloaded models.
        If False, show all available models (fetch list from server).
    verbose : bool
        Whether to show detailed information.

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing model names and descriptions.
    """
    import os
    import platform

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["KMP_WARNINGS"] = "0"
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

    is_apple_silicon = platform.system() == "Darwin" and platform.machine() == "arm64"
    if is_apple_silicon:
        configure_apple_silicon_env(verbose=verbose)

    try:
        import celltypist
        from celltypist import models
    except ImportError as e:
        log_message(
            "{.pkg celltypist} import failed: {.val {e}}",
            message_type="error",
            verbose=verbose,
        )
        raise

    try:
        models_df = models.models_description(on_the_fly=on_the_fly)
        return models_df
    except Exception as e:
        log_message(
            "{.pkg celltypist} failed to get model list: {.val {e}}",
            message_type="error",
            verbose=verbose,
        )
        raise


def CellTypist(
    adata=None,
    h5ad=None,
    model="Immune_All_Low.pkl",
    mode="best match",
    p_thres=0.5,
    majority_voting=False,
    over_clustering=None,
    min_prop=0,
    use_GPU=False,
    insert_labels=True,
    insert_conf=True,
    insert_conf_by="predicted_labels",
    insert_prob=False,
    insert_decision=False,
    prefix="",
    verbose=True,
):
    """
    Run CellTypist cell type annotation.

    Parameters
    ----------
    adata : AnnData
        AnnData object (required if h5ad is not provided).
    h5ad : str
        Path to h5ad file (optional).
    model : str
        Model name or path. Default is 'Immune_All_Low.pkl'.
        Supports three formats:
        1. Model name (e.g., 'Immune_All_Low.pkl'): automatically searched in ~/.celltypist/data/models/
        2. Full path (contains '/'): use the provided path directly
        3. None: use default model via celltypist.models.get_default_model()
    mode : str
        Prediction mode: 'best match' or 'prob match'. Default is 'best match'.
    p_thres : float
        Probability threshold for 'prob match' mode. Default is 0.5.
    majority_voting : bool
        Whether to use majority voting. Default is False.
    over_clustering : str, list, or None
        Over-clustering result. Can be:
        - String: column name in adata.obs
        - List/array: over-clustering labels
        - None: use heuristic over-clustering
    min_prop : float
        Minimum proportion for majority voting. Default is 0.
    use_GPU : bool
        Whether to use GPU for over-clustering. Default is False.
    insert_labels : bool
        Whether to insert predicted labels into AnnData. Default is True.
    insert_conf : bool
        Whether to insert confidence scores. Default is True.
    insert_conf_by : str
        Which prediction type to base confidence on. Default is 'predicted_labels'.
    insert_prob : bool
        Whether to insert probability matrix. Default is False.
    insert_decision : bool
        Whether to insert decision matrix. Default is False.
    prefix : str
        Prefix for inserted columns. Default is empty string.
    verbose : bool
        Whether to show detailed information. Default is True.

    Returns
    -------
    AnnData
        AnnData object with CellTypist predictions inserted.
    """
    import os
    import platform
    import sys

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["KMP_WARNINGS"] = "0"
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

    is_apple_silicon = platform.system() == "Darwin" and platform.machine() == "arm64"
    if is_apple_silicon:
        configure_apple_silicon_env(
            scanpy_settings=True,
            scanpy_verbosity=True,
            numba_threading=True,
            numba_disable_jit=True,
            configure_numba_runtime=True,
            numba_runtime_optional=True,
            verbose=verbose,
        )

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        log_message(
            "matplotlib setup failed: {.val {e}}",
            message_type="warning",
            verbose=verbose,
        )
        plt = None

    try:
        import celltypist
        from celltypist import models
    except ImportError as e:
        log_message(
            "{.pkg celltypist} import failed: {.val {e}}",
            message_type="error",
            verbose=verbose,
        )
        raise

    import warnings

    warnings.simplefilter("ignore", category=UserWarning)
    warnings.simplefilter("ignore", category=FutureWarning)
    warnings.simplefilter("ignore", category=DeprecationWarning)

    try:
        if adata is None and h5ad is None:
            log_message(
                "{.arg adata} or {.arg h5ad} must be provided",
                message_type="error",
                verbose=verbose,
            )
            return None

        if adata is None:
            import scanpy as sc

            adata = sc.read(h5ad)

        import numpy as np
        import scanpy as sc

        needs_preprocessing = False

        if adata.X is not None:
            sample_size = min(100, adata.n_obs)
            try:
                if hasattr(adata.X, "toarray"):
                    sample_data = adata.X[:sample_size].toarray()
                else:
                    sample_data = adata.X[:sample_size]

                min_val = float(np.min(sample_data))
                max_val = float(np.max(sample_data))

                if min_val < 0 or max_val > 9.22:
                    needs_preprocessing = True
                    log_message(
                        "Data format invalid for CellTypist (min={.val {format(min_val, '.2f')}}, max={.val {format(max_val, '.2f')}}). Will attempt to preprocess...",
                        message_type="warning",
                        verbose=verbose,
                    )
            except Exception as e:
                log_message(
                    "Could not check data format: {.val {e}}. Attempting preprocessing...",
                    message_type="warning",
                    verbose=verbose,
                )
                needs_preprocessing = True

        if needs_preprocessing:
            log_message(
                "Preprocessing data for CellTypist (normalize to 10000 counts per cell and log1p transform)...",
                message_type="info",
                verbose=verbose,
            )

            if adata.raw is not None:
                log_message(
                    "Using .raw data for preprocessing",
                    message_type="info",
                    verbose=verbose,
                )
                adata = adata.raw.to_adata()
            elif adata.X is not None:
                try:
                    if hasattr(adata.X, "toarray"):
                        sample_check = adata.X[: min(10, adata.n_obs)].toarray()
                    else:
                        sample_check = adata.X[: min(10, adata.n_obs)]

                    if np.all(sample_check >= 0) and np.max(sample_check) > 20:
                        adata.raw = adata.copy()
                        log_message(
                            "Detected raw counts in .X, storing in .raw",
                            message_type="info",
                            verbose=verbose,
                        )
                except:
                    pass

            sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
            sc.pp.log1p(adata)

        if model is None:
            model = models.get_default_model()
            log_message(
                "Using default model: {.val {model}}",
                message_type="info",
                verbose=verbose,
            )
        elif "/" not in model:
            try:
                all_models = models.get_all_models()
                if model in all_models:
                    model = models.get_model_path(model)
                    log_message(
                        "Using model: {.val {model}}",
                        message_type="info",
                        verbose=verbose,
                    )
                else:
                    log_message(
                        "Model {.val {model}} not found in downloaded models. Will try to use as path or download if needed.",
                        message_type="warning",
                        verbose=verbose,
                    )
            except Exception as e:
                log_message(
                    "Could not check model list: {.val {e}}. Using model as provided.",
                    message_type="warning",
                    verbose=verbose,
                )

        over_clustering_value = None
        if over_clustering is not None:
            if isinstance(over_clustering, str):
                if over_clustering in adata.obs.columns:
                    over_clustering_value = adata.obs[over_clustering].values
                else:
                    log_message(
                        "{.val {over_clustering}} not found in adata.obs. Will use heuristic over-clustering.",
                        message_type="warning",
                        verbose=verbose,
                    )
            else:
                over_clustering_value = over_clustering

        if majority_voting and over_clustering_value is None:
            # CellTypist may detect graph matrices in `obsp` but still requires
            # `uns["neighbors"]` metadata for Scanpy Leiden over-clustering.
            has_neighbors_uns = isinstance(adata.uns.get("neighbors"), dict)
            has_graph_in_obsp = (
                "connectivities" in adata.obsp and "distances" in adata.obsp
            )

            if not has_neighbors_uns:
                if has_graph_in_obsp:
                    adata.uns["neighbors"] = {
                        "connectivities_key": "connectivities",
                        "distances_key": "distances",
                        "params": {
                            "method": "umap",
                            "metric": "euclidean",
                        },
                    }
                    log_message(
                        "Added missing {.field adata.uns['neighbors']} from existing graph in {.field adata.obsp} for majority voting.",
                        message_type="info",
                        verbose=verbose,
                    )
                else:
                    log_message(
                        "No precomputed neighbors graph found; computing neighbors for CellTypist majority voting...",
                        message_type="info",
                        verbose=verbose,
                    )
                    sc.pp.neighbors(adata, n_neighbors=15)

        log_message("Running CellTypist annotation...", verbose=verbose)

        predictions = celltypist.annotate(
            filename=adata,
            model=model,
            mode=mode,
            p_thres=p_thres,
            majority_voting=majority_voting,
            over_clustering=over_clustering_value,
            min_prop=min_prop,
            use_GPU=use_GPU,
        )

        adata = predictions.to_adata(
            insert_labels=insert_labels,
            insert_conf=insert_conf,
            insert_conf_by=insert_conf_by,
            insert_decision=insert_decision,
            insert_prob=insert_prob,
            prefix=prefix,
        )

        log_message(
            "CellTypist annotation completed successfully",
            message_type="success",
            verbose=verbose,
        )

        return adata

    except Exception as e:
        log_message(
            "CellTypist annotation failed: {.val {e}}",
            message_type="error",
            verbose=verbose,
        )
        raise
