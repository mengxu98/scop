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
    save=False,
    dpi=300,
    dirpath="./",
    fileprefix="",
):
    # Configure OpenMP settings to prevent conflicts
    import os

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["KMP_WARNINGS"] = "0"
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

    import matplotlib

    matplotlib.use("Agg")  # Use non-interactive backend
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
    os.chdir(os.path.expanduser(dirpath))

    try:
        # Input validation
        if adata is None and h5ad is None:
            raise ValueError("Either 'adata' or 'h5ad' must be provided.")

        if adata is None:
            adata = scv.read(h5ad)

        if group_by is None:
            raise ValueError("'group_by' must be provided.")

        if linear_reduction is None and nonlinear_reduction is None:
            raise ValueError(
                "At least one of 'linear_reduction' or 'nonlinear_reduction' must be provided."
            )

        # Setup basis
        if basis is None:
            if nonlinear_reduction is not None:
                # Check if the nonlinear reduction exists in obsm
                if nonlinear_reduction in adata.obsm:
                    basis = nonlinear_reduction
                elif f"X_{nonlinear_reduction}" in adata.obsm:
                    basis = f"X_{nonlinear_reduction}"
                else:
                    print(
                        f"Warning: {nonlinear_reduction} not found in adata.obsm. Available keys: {list(adata.obsm.keys())}"
                    )
                    basis = (
                        linear_reduction
                        if linear_reduction in adata.obsm
                        else f"X_{linear_reduction}"
                    )
            else:
                basis = (
                    linear_reduction
                    if linear_reduction in adata.obsm
                    else f"X_{linear_reduction}"
                )

        # Ensure the basis exists in obsm
        if basis not in adata.obsm:
            print(
                f"Warning: basis '{basis}' not found in adata.obsm. Available keys: {list(adata.obsm.keys())}"
            )
            # Try to find alternative basis
            if linear_reduction in adata.obsm:
                basis = linear_reduction
                print(f"Using {linear_reduction} as basis instead.")
            elif f"X_{linear_reduction}" in adata.obsm:
                basis = f"X_{linear_reduction}"
                print(f"Using X_{linear_reduction} as basis instead.")
            else:
                # Create a 2D basis from linear reduction
                if linear_reduction in adata.obsm:
                    adata.obsm["basis"] = adata.obsm[linear_reduction][:, 0:2]
                    basis = "basis"
                elif f"X_{linear_reduction}" in adata.obsm:
                    adata.obsm["basis"] = adata.obsm[f"X_{linear_reduction}"][:, 0:2]
                    basis = "basis"
                else:
                    raise ValueError(
                        f"Cannot find suitable basis. Available obsm keys: {list(adata.obsm.keys())}"
                    )

        print(f"Using basis: {basis}")
        print(f"Available embeddings in adata.obsm: {list(adata.obsm.keys())}")

        # Ensure group_by is categorical
        adata.obs[group_by] = adata.obs[group_by].astype("category")

        # === PREPROCESSING PHASE ===
        print("=== Starting preprocessing ===")

        # 1. Gene filtering (optional)
        if filter_genes:
            print("Filtering genes...")
            scv.pp.filter_genes(adata, min_counts=min_counts)
            scv.pp.filter_genes(adata, min_counts_u=min_counts_u)

        # 2. Normalization and transformation
        if normalize_per_cell:
            print("Normalizing per cell...")
            scv.pp.normalize_per_cell(adata)

        if log_transform:
            print("Log transforming...")
            sc.pp.log1p(adata)

        # 3. Magic imputation (if requested)
        if magic_impute:
            print("Performing MAGIC imputation...")
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
                print("Warning: magic-impute not installed. Skipping imputation.")

        # 4. Compute neighbors and moments with version compatibility
        print("Computing neighbors and moments...")
        try:
            # Method 1: Try using scVelo's workflow
            scv.pp.moments(
                adata, n_pcs=n_pcs, n_neighbors=n_neighbors, use_rep=linear_reduction
            )
        except Exception as e:
            print(f"Warning: scVelo moments failed ({e}), using manual computation...")
            # Method 2: Manual computation for compatibility
            sc.pp.neighbors(
                adata, n_pcs=n_pcs, n_neighbors=n_neighbors, use_rep=linear_reduction
            )

            # Manual moments calculation
            connectivities = adata.obsp["connectivities"]
            if sparse.issparse(adata.layers["spliced"]):
                Ms = connectivities @ adata.layers["spliced"]
                Mu = connectivities @ adata.layers["unspliced"]
            else:
                Ms = connectivities @ sparse.csr_matrix(adata.layers["spliced"])
                Mu = connectivities @ sparse.csr_matrix(adata.layers["unspliced"])

            adata.layers["Ms"] = Ms
            adata.layers["Mu"] = Mu

        # === VELOCITY ESTIMATION PHASE ===
        print("=== Starting velocity estimation ===")

        # Process each mode
        for m in mode:
            print(f"Processing mode: {m}")

            if m == "dynamical":
                # Dynamical modeling
                print("Performing dynamical modeling...")
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
                    print(
                        "Warning: No genes found for dynamical modeling. Using all genes."
                    )
                    scv.tl.recover_dynamics(adata, use_raw=use_raw, n_jobs=n_jobs)

            # Compute velocity
            print(f"Computing {m} velocity...")
            scv.tl.velocity(adata, mode=m, diff_kinetics=diff_kinetics)

            # Compute velocity graph
            print(f"Computing {m} velocity graph...")
            scv.tl.velocity_graph(adata, vkey=m, n_neighbors=n_neighbors, n_jobs=n_jobs)

            # === DOWNSTREAM ANALYSIS ===
            print(f"=== Downstream analysis for {m} ===")

            # Velocity embedding
            print("Computing velocity embedding...")
            scv.tl.velocity_embedding(adata, basis=basis, vkey=m)

            # Velocity confidence (with error handling)
            if compute_velocity_confidence:
                print("Computing velocity confidence...")
                try:
                    scv.tl.velocity_confidence(adata, vkey=m)
                except Exception as e:
                    print(
                        f"Warning: velocity confidence failed ({e}), using default values..."
                    )
                    n_obs = adata.n_obs
                    adata.obs[m + "_length"] = np.ones(n_obs) * 0.5
                    adata.obs[m + "_confidence"] = np.ones(n_obs) * 0.5

            # Terminal states
            if compute_terminal_states:
                print("Computing terminal states...")
                try:
                    scv.tl.terminal_states(adata, vkey=m)
                    # Rename for consistency
                    for term in ["root_cells", "end_points"]:
                        if term in adata.obs.columns:
                            adata.obs[m + "_" + term] = adata.obs[term]
                            adata.obs.drop(term, axis=1, inplace=True)
                except Exception as e:
                    print(f"Warning: terminal states computation failed: {e}")

            # Pseudotime
            if compute_pseudotime:
                print("Computing velocity pseudotime...")
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
                    print(f"Warning: pseudotime computation failed: {e}")

            # PAGA
            if compute_paga:
                print("Computing PAGA...")
                try:
                    # Ensure neighbors info is available
                    if "neighbors" not in adata.uns:
                        adata.uns["neighbors"] = {}
                    adata.uns["neighbors"]["distances"] = adata.obsp["distances"]
                    adata.uns["neighbors"]["connectivities"] = adata.obsp[
                        "connectivities"
                    ]

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
                    scv.tl.paga(
                        adata,
                        groups=group_by,
                        vkey=m,
                        root_key=root_key,
                        end_key=end_key,
                    )
                except Exception as e:
                    print(f"Warning: PAGA computation failed: {e}")

            # Velocity genes ranking
            if calculate_velocity_genes:
                print("Ranking velocity genes...")
                try:
                    if m != "dynamical":
                        scv.tl.rank_velocity_genes(adata, vkey=m, groupby=group_by)
                        if "spearmans_score" in adata.var.columns:
                            adata.var[m + "_score"] = adata.var["spearmans_score"]
                    else:
                        scv.tl.rank_dynamical_genes(adata, groupby=group_by)
                except Exception as e:
                    print(f"Warning: velocity genes ranking failed: {e}")

            # === VISUALIZATION ===
            if show_plot:
                print(f"Generating plots for {m}...")

                # Setup palette
                groups = (
                    adata.obs[group_by].cat.categories
                    if hasattr(adata.obs[group_by], "cat")
                    else adata.obs[group_by].unique()
                )
                if palette is None:
                    palette = dict(
                        zip(groups, plt.cm.tab10(np.linspace(0, 1, len(groups))))
                    )

                # Velocity stream plot
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
                        legend_loc="right margin",
                        save=f"{fileprefix}_{m}_stream.png" if save else False,
                        show=show_plot,
                    )
                except Exception as e:
                    print(f"Warning: stream plot failed: {e}")

                # Velocity arrow plot
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
                        save=f"{fileprefix}_{m}_arrow.png" if save else False,
                        show=show_plot,
                    )
                except Exception as e:
                    print(f"Warning: arrow plot failed: {e}")

                # Confidence plots
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
                                    save=f"{fileprefix}_{m}_{metric}.png"
                                    if save
                                    else False,
                                    show=show_plot,
                                )
                            except Exception as e:
                                print(f"Warning: {metric} plot failed: {e}")

        print("=== scVelo analysis completed ===")

    except Exception as e:
        print(f"Error in SCVELO analysis: {e}")
        raise
    finally:
        os.chdir(prevdir)

    # Clean up adata for return
    try:
        if hasattr(adata, "_raw") and adata._raw is not None:
            if hasattr(adata._raw, "_var"):
                adata._raw._var = adata._raw._var.rename(columns={"_index": "features"})
    except Exception:
        pass

    return adata


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
    s_genes=None,
    g2m_genes=None,
    calculate_velocity_genes=False,
    denoise=False,
    kinetics=False,
    n_jobs=1,
    show_plot=True,
    dpi=300,
    save=False,
    dirpath="./",
    fileprefix="",
):
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

    import os

    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(dirpath))

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
            print("adata or h5ad must be provided.")
            exit()

        if adata is None:
            adata = scv.read(h5ad)
        # del adata.uns

        if group_by is None:
            print("group_by must be provided.")
            exit()

        if linear_reduction is None and nonlinear_reduction is None:
            print(
                "linear_reduction or nonlinear_reduction must be provided at least one."
            )
            exit()

        if linear_reduction is None:
            sc.pp.pca(adata, n_comps=n_pcs)
            linear_reduction = "X_pca"

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
            print("'fitting_by' must be one of 'deterministic' and 'stochastic'.")
            exit()

        if not all([m in ["deterministic", "stochastic", "dynamical"] for m in mode]):
            print(
                "Invalid mode name! Must be the 'deterministic', 'stochastic' or 'dynamical'."
            )
            exit()

        adata.obs[group_by] = adata.obs[group_by].astype(dtype="category")

        if mode[-1] + "_graph" not in adata.obs.keys():
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
                s_genes=s_genes,
                g2m_genes=g2m_genes,
                show_plot=show_plot,
                dpi=dpi,
                save=save,
                dirpath=dirpath,
                fileprefix=fileprefix,
            )
        adata.layers["velocity"] = adata.layers[mode[-1]]
        cr.tl.terminal_states(adata, cluster_key=group_by)
        cr.pl.terminal_states(adata)
        cr.tl.initial_states(adata, cluster_key=group_by)
        cr.pl.initial_states(adata)
        cr.tl.lineages(adata)
        cr.pl.lineages(adata, same_plot=False)
        cr.pl.lineages(adata, same_plot=True)

        scv.tl.recover_latent_time(
            adata, root_key="initial_states_probs", end_key="terminal_states_probs"
        )
        scv.tl.paga(
            adata,
            groups=group_by,
            root_key="initial_states_probs",
            end_key="terminal_states_probs",
            use_time_prior="velocity_pseudotime",
        )
        cr.pl.cluster_fates(
            adata,
            mode="paga_pie",
            cluster_key=group_by,
            basis=basis,
            legend_kwargs={"loc": "top right out"},
            legend_loc="top left out",
            node_size_scale=5,
            edge_width_scale=1,
            max_edge_width=4,
            title="directed PAGA",
        )
        if show_plot is True:
            plt.show()

        cr.tl.lineage_drivers(adata, cluster_key=group_by)
        cr.pl.lineage_drivers(adata, lineage=adata.obs[group_by].unique()[1], n_genes=4)
        if show_plot is True:
            plt.show()

    finally:
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
    show_plot=True,
    dpi=300,
    save=False,
    dirpath="./",
    fileprefix="",
):
    # Configure OpenMP settings to prevent conflicts
    import os

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["KMP_WARNINGS"] = "0"

    import matplotlib

    matplotlib.use("Agg")  # Use non-interactive backend
    import matplotlib.pyplot as plt
    import scanpy as sc
    import numpy as np
    import statistics
    from math import hypot

    import warnings

    warnings.simplefilter("ignore", category=UserWarning)
    warnings.simplefilter("ignore", category=FutureWarning)
    warnings.simplefilter("ignore", category=DeprecationWarning)

    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(dirpath))

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
            print("adata or h5ad must be provided.")
            exit()

        if adata is None:
            adata = sc.read(h5ad)

        if group_by is None:
            print("group_by must be provided.")
            exit()

        if linear_reduction is None and nonlinear_reduction is None:
            print(
                "linear_reduction or nonlinear_reduction must be provided at least one."
            )
            exit()

        if linear_reduction is None:
            sc.pp.pca(adata, n_comps=n_pcs)
            linear_reduction = "X_pca"

        if basis is None:
            if nonlinear_reduction is not None:
                basis = nonlinear_reduction
            else:
                basis = "basis"
                adata.obsm["basis"] = adata.obsm[linear_reduction][:, 0:2]

        if point_size is None:
            point_size = min(100000 / adata.shape[0], 20)

        if infer_pseudotime is True and root_cell is None and root_group is None:
            print("root_cell or root_group should be provided.")
            exit()

        if use_rna_velocity is True:
            adata.uns["velocity_graph"] = adata.uns[vkey + "_graph"]
        # del adata.uns

        adata.obs[group_by] = adata.obs[group_by].astype(dtype="category")

        if "X_diffmap" in adata.obsm_keys():
            X_diffmap = adata.obsm["X_diffmap"]
            del adata.obsm["X_diffmap"]
            sc.pp.neighbors(
                adata, n_pcs=n_pcs, use_rep=linear_reduction, n_neighbors=n_neighbors
            )
            adata.obsm["X_diffmap"] = X_diffmap
        else:
            sc.pp.neighbors(
                adata, n_pcs=n_pcs, use_rep=linear_reduction, n_neighbors=n_neighbors
            )

        sc.tl.paga(adata, groups=group_by, use_rna_velocity=use_rna_velocity)

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
                show=False,
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
                show=False,
            )
        if show_plot is True:
            plt.show()
        if save:
            plt.savefig(
                ".".join(filter(None, [fileprefix, "paga_compare.png"])), dpi=dpi
            )

        sc.pl.paga(
            adata,
            threshold=threshold,
            layout=paga_layout,
            title="PAGA layout: " + paga_layout,
            frameon=False,
            save=False,
            show=False,
        )
        if show_plot is True:
            plt.show()
        if save:
            plt.savefig(
                ".".join(filter(None, [fileprefix, "paga_layout.png"])), dpi=dpi
            )

        sc.tl.draw_graph(adata, init_pos="paga", layout=paga_layout)
        sc.pl.draw_graph(
            adata,
            color=group_by,
            palette=palette,
            title="PAGA layout: " + paga_layout,
            layout=paga_layout,
            frameon=False,
            legend_loc="on data",
            show=False,
        )
        if show_plot is True:
            plt.show()
        if save:
            plt.savefig(".".join(filter(None, [fileprefix, "paga_graph.png"])), dpi=dpi)

        if embedded_with_PAGA is True:
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
                show=False,
            )
            if show_plot is True:
                plt.show()
            if save:
                plt.savefig(
                    ".".join(filter(None, [fileprefix, "paga_umap.png"])), dpi=dpi
                )

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

            sc.pl.embedding(
                adata, basis=basis, color="dpt_pseudotime", save=False, show=False
            )
            if show_plot is True:
                plt.show()
            if save:
                plt.savefig(
                    ".".join(filter(None, [fileprefix, "dpt_pseudotime.png"])), dpi=dpi
                )

    finally:
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
):
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

    import os

    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(dirpath))

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
            print("adata or h5ad must be provided.")
            exit()

        if adata is None:
            adata = sc.read(h5ad)
        # del adata.uns

        if group_by is None and (
            early_group is not None or terminal_groups is not None
        ):
            print("group_by must be provided.")
            exit()

        if linear_reduction is None and nonlinear_reduction is None:
            print(
                "linear_reduction or nonlinear_reduction must be provided at least one."
            )
            exit()

        if linear_reduction is None:
            sc.pp.pca(adata, n_comps=n_pcs)
            linear_reduction = "X_pca"

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
            print("early_cell must be provided.")
            exit()
        else:
            print("early_cell: ", early_cell)

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
            print("terminal_cells: None")
        else:
            print("terminal_cells: ", terminal_cells)

        # if HVF is None:
        #   sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        # else:
        #   df = pd.DataFrame([False] * adata.X.shape[1],columns=["highly_variable"])
        #   df = df.set_index(adata.var_names)
        #   df.highly_variable.iloc[:n_top_genes] = True
        #   df.loc[df['channel'].isin(['sale','fullprice'])]
        #   df.highly_variable.iloc[df.index.isin(HVF)] = True
        #   if "highly_variable" in adata.var.columns:
        #     adata.var.drop('highly_variable', axis=1, inplace=True)
        #   adata.var=adata.var.join(df)

        # adata.uns['pca']['variance_ratio']
        # pca_projections=n_comps = np.where(np.cumsum(ad.uns['pca']['variance_ratio']) > 0.85)[0][0]
        # pca_projections, _ = palantir.utils.run_pca(adata, use_hvg=True)

        pca_projections = pd.DataFrame(
            adata.obsm[linear_reduction][:, :n_pcs], index=adata.obs_names
        )
        print("running diffusion maps")
        dm_res = palantir.utils.run_diffusion_maps(
            pca_projections,
            n_components=dm_n_components,
            knn=n_neighbors,
            alpha=dm_alpha,
        )
        ms_data = palantir.utils.determine_multiscale_space(dm_res, n_eigs=dm_n_eigs)
        print("running palantir")
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
            adata, basis=basis, color="palantir_pseudotime", size=point_size
        )
        if save:
            plt.savefig(
                ".".join(filter(None, [fileprefix, "palantir_pseudotime.png"])), dpi=dpi
            )

        sc.pl.embedding(
            adata, basis=basis, color="palantir_diff_potential", size=point_size
        )
        if save:
            plt.savefig(
                ".".join(filter(None, [fileprefix, "palantir_diff_potential.png"])),
                dpi=dpi,
            )

        sc.pl.embedding(
            adata,
            basis=basis,
            color=pr_res.branch_probs.columns.values,
            size=point_size,
        )
        if save:
            plt.savefig(
                ".".join(filter(None, [fileprefix, "palantir_probs.png"])), dpi=dpi
            )

    finally:
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
):
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
    os.chdir(os.path.expanduser(dirpath))

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
            print("adata or h5ad must be provided.")
            exit()

        if adata is None:
            adata = sc.read(h5ad)

        if group_by is None:
            print("group_by must be provided.")
            exit()

        if time_field is None:
            print("time_field must be provided.")
            exit()

        adata.obs[group_by] = adata.obs[group_by].astype(dtype="category")
        if pd.api.types.is_categorical_dtype(adata.obs[time_field]):
            adata.obs["time_field"] = adata.obs[time_field].cat.codes
        elif not pd.api.types.is_numeric_dtype(adata.obs[time_field]):
            try:
                adata.obs["time_field"] = adata.obs[time_field].astype("float")
            except ValueError:
                print("Unable to convert column '" + time_field + "' to float type.")
        else:
            adata.obs["time_field"] = adata.obs[time_field]

        time_dict = dict(zip(adata.obs[time_field], adata.obs["time_field"]))
        if time_from not in time_dict.keys():
            print("'time_from' is incorrect")
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

        # obs_list = wot.tmap.trajectory_trends_from_trajectory(trajectory_ds = trajectory_ds, expression_ds = adata)

        if time_to is not None:
            if time_to not in time_dict.keys():
                print("'time_to' is incorrect")
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
