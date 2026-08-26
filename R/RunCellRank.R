#' @title Run CellRank analysis
#'
#' @description
#' CellRank is a toolkit for studying cellular dynamics using Markov state modeling.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams CellDimPlot
#' @inheritParams srt_to_adata
#' @inheritParams RunStandardWorkflow
#' @inheritParams scop-params
#' @param srt A Seurat object.
#' If provided, `adata` will be ignored.
#' @param adata An anndata object.
#' @param basis The basis to use for reduction, e.g., `"UMAP"`.
#' @param n_pcs Number of principal components to use for linear reduction.
#' @param n_neighbors Number of neighbors to use for constructing the KNN graph.
#' @param cores The number of cores to use for `cellrank`.
#' @param legend.position Position of legend in plots.
#' Can be `"on data"`, `"right margin"`, `"bottom right"`, etc.
#' @param show_plot Whether to show the plot.
#' @param save_plot Whether to save plots to files.
#' @param plot_format Format for saved plots: `"png"` (default), `"pdf"`, or `"svg"`.
#' @param plot_dpi Resolution (DPI) for saved plots.
#' @param plot_prefix Prefix for saved plot filenames.
#' @param dirpath The directory to save the plots.
#' @param return_seurat Whether to return a Seurat object instead of an anndata object.
#' @param mode Velocity estimation models to use.
#' Can be `"deterministic"`, `"stochastic"`, or `"dynamical"`.
#' If the corresponding velocity reduction (e.g. `stochastic_umap`) is not
#' found, falls back to the scVelo convention (e.g. `velocity_umap`), which is
#' what an object converted from an AnnData (via [adata_to_srt]) contains.
#' @param fitting_by Method used to fit gene velocities for dynamical modeling.
#' @param magic_impute Flag indicating whether to perform magic imputation.
#' @param knn The number of nearest neighbors for `magic.MAGIC`.
#' @param t Power to which the diffusion operator is powered for `magic.MAGIC`.
#' @param min_shared_counts Minimum number of counts (both unspliced and spliced) required for a gene.
#' @param stream_smooth Multiplication factor for scale in Gaussian kernel around grid point.
#' @param stream_density Controls the closeness of streamlines.
#' When density = 2 (default), the domain is divided into a 60x60 grid,
#' whereas density linearly scales this grid.
#' Each cell in the grid can have, at most, one traversing streamline.
#' @param arrow_size Size of arrows.
#' @param arrow_length Length of arrows.
#' @param arrow_density Amount of velocities to show.
#' @param calculate_velocity_genes Boolean flag indicating whether to calculate velocity genes.
#' @param denoise Boolean flag indicating whether to denoise.
#' @param kinetics Boolean flag indicating whether to estimate RNA kinetics.
#' @param kernel_type Type of kernel to use: `"velocity"` (default, requires spliced/unspliced),
#' `"pseudotime"` (requires pre-computed pseudotime or auto-computes DPT),
#' `"cytotrace"` (auto-computes CytoTRACE score, suitable for RNA-only data),
#' or `"wot"` (uses Waddington-OT transport maps through CellRank's RealTimeKernel).
#' @param time_key Key in metadata for pseudotime. Used when `kernel_type = "pseudotime"`.
#' If the key doesn't exist, DPT pseudotime will be computed automatically.
#' @param time_field Key in metadata for experimental time. Used when `kernel_type = "wot"`.
#' @param growth_iters Number of growth iterations passed to `wot.ot.OTModel`.
#' @param tmap_out Directory used to store or read Waddington-OT transport maps.
#' @param recalculate Whether to recompute Waddington-OT transport maps even when
#' `tmap_out` already exists.
#' @param estimator_type Type of estimator to use: `"GPCCA"` (default) or `"CFLARE"`.
#' GPCCA provides coarse-grained analysis and Schur decomposition.
#' @param use_connectivity_kernel Whether to combine the main kernel with ConnectivityKernel.
#' @param velocity_weight Weight for the VelocityKernel when combining with ConnectivityKernel.
#' @param connectivity_weight Weight for the ConnectivityKernel when combining with VelocityKernel.
#' Weights are automatically normalized to sum to `1.0`.
#' @param softmax_scale Scaling parameter for softmax transformation of velocity kernel.
#' @param n_macrostates Number of macrostates to compute.
#' If `NULL` (default), automatically determined based on eigenvalue spectrum.
#' @param schur_method Method for Schur decomposition: `"krylov"` or `"brandts"`.
#' Only used for GPCCA estimator.
#' @param n_cells_terminal Minimum number of cells required for a state to be considered terminal.
#' @param schur_n_components Number of Schur components. If `NULL`, retain the
#' existing size heuristic.
#' @param terminal_states Optional macrostate names to set as terminal states.
#' Names are validated against the computed macrostates before fate estimation.
#' @param terminal_state_agg Aggregation policy for combined terminal states.
#' @param driver_lineages Optional lineage names for driver-gene computation.
#' @param compute_lineage_drivers Whether to compute and store lineage drivers.
#'
#' @return
#' Returns a Seurat object if `return_seurat = TRUE` or an anndata object with CellRank results stored in `obsm`, `obs`, and `varm` slots.
#' The estimator and kernel objects are stored in `srt@misc$cellrank`.
#'
#' @export
#' @seealso
#' [RunSCVELO], [RunPAGA], [VelocityPlot], [CellDimPlot], [DynamicPlot],
#' [RunCellRankTrends], [RunCellRankEnrichment], [CellRankPlot]
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#' pancreas_sub <- RunCellRank(
#'   srt = pancreas_sub,
#'   group.by = "SubCellType"
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "term_states_fwd",
#'   reduction = "umap",
#'   label = TRUE
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   features = "latent_time",
#'   reduction = "umap"
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   features = c("stochastic_confidence", "stochastic_length"),
#'   reduction = "umap"
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   lineages = "cellrank_pseudotime",
#'   lineages_span = 0.1,
#'   lineages_trim = c(0.05, 0.95)
#' )
#'
#' DynamicPlot(
#'   pancreas_sub,
#'   lineages = "cellrank_pseudotime",
#'   features = c("Arxes1", "Ncoa2"),
#'   group.by = "SubCellType"
#' )
#' }
#' @param backward Whether to compute backward transitions.
#' @param backend Backend for computation: `"python"` (default) or `"cpp"`.
#' The C++ path is an explicitly opted-in approximation and does not reproduce
#' the complete CellRank estimator.
#' @param allow_approximate Whether to allow the approximate C++ path. This must
#' be `TRUE` when `backend = "cpp"`.
#' @param max_dense_gib Maximum estimated GiB allowed for the dense
#' cell-by-cell working matrices used by the C++ path.
#' @param envname Optional Python environment name. `NULL` uses the current
#' SCOP environment selection.
#' @param conda Conda-compatible executable used by [PrepareEnv].
#' @param recompute_neighbors Whether to rebuild Scanpy neighbors before a
#' pseudotime kernel. If `FALSE`, the existing graph must be complete and
#' match the current AnnData cell order.
RunCellRank <- function(
  srt = NULL,
  assay_x = "RNA",
  layer_x = "counts",
  assay_y = c("spliced", "unspliced"),
  layer_y = "counts",
  adata = NULL,
  group.by = NULL,
  cores = 1,
  linear_reduction = NULL,
  nonlinear_reduction = NULL,
  basis = NULL,
  mode = "stochastic",
  fitting_by = "stochastic",
  magic_impute = FALSE,
  knn = 5,
  t = 2,
  min_shared_counts = 30,
  n_pcs = 30,
  n_neighbors = 30,
  stream_smooth = NULL,
  stream_density = 2,
  arrow_size = 5,
  arrow_length = 5,
  arrow_density = 0.5,
  calculate_velocity_genes = FALSE,
  denoise = FALSE,
  kinetics = FALSE,
  kernel_type = c("velocity", "pseudotime", "cytotrace", "wot"),
  time_key = "dpt_pseudotime",
  time_field = "Time",
  growth_iters = 3L,
  tmap_out = "tmaps/tmap_out",
  recalculate = FALSE,
  estimator_type = c("GPCCA", "CFLARE"),
  use_connectivity_kernel = TRUE,
  velocity_weight = 0.8,
  connectivity_weight = 0.2,
  softmax_scale = 4,
  n_macrostates = NULL,
  schur_method = c("krylov", "brandts"),
  schur_n_components = NULL,
  n_cells_terminal = 10,
  terminal_states = NULL,
  terminal_state_agg = c("top_n", "union"),
  driver_lineages = NULL,
  compute_lineage_drivers = TRUE,
  backward = FALSE,
  backend = c("python", "cpp"),
  allow_approximate = FALSE,
  max_dense_gib = 8,
  show_plot = TRUE,
  save_plot = FALSE,
  plot_format = c("pdf", "png", "svg"),
  plot_dpi = 300,
  plot_prefix = "cellrank",
  legend.position = "on data",
  palette = "Chinese",
  palcolor = NULL,
  dirpath = "./cellrank",
  envname = NULL,
  conda = "auto",
  recompute_neighbors = TRUE,
  return_seurat = !is.null(srt),
  verbose = TRUE
) {
  kernel_type <- match.arg(kernel_type)
  backend <- match.arg(backend)
  estimator_type_upper <- toupper(match.arg(estimator_type))
  terminal_state_agg <- match.arg(terminal_state_agg)

  if (identical(backend, "cpp")) {
    assert_cpp_approximation_opt_in(
      allow_approximate,
      "RunCellRank(backend = \"cpp\")"
    )
    unsupported_cpp <- character()
    if (isTRUE(show_plot)) {
      unsupported_cpp <- c(unsupported_cpp, "show_plot")
    }
    if (isTRUE(save_plot)) {
      unsupported_cpp <- c(unsupported_cpp, "save_plot")
    }
    if (!is.null(schur_n_components)) {
      unsupported_cpp <- c(unsupported_cpp, "schur_n_components")
    }
    if (!is.null(terminal_states)) {
      unsupported_cpp <- c(unsupported_cpp, "terminal_states")
    }
    if (!identical(terminal_state_agg, "top_n")) {
      unsupported_cpp <- c(unsupported_cpp, "terminal_state_agg")
    }
    if (!is.null(driver_lineages)) {
      unsupported_cpp <- c(unsupported_cpp, "driver_lineages")
    }
    if (!isTRUE(compute_lineage_drivers)) {
      unsupported_cpp <- c(unsupported_cpp, "compute_lineage_drivers")
    }
    reject_unsupported_cpp_arguments(
      unsupported_cpp,
      "RunCellRank(backend = \"cpp\")"
    )
    if (!is.null(srt)) {
      assert_cpp_dense_budget(
        n_rows = ncol(srt),
        n_cols = ncol(srt),
        copies = 8,
        max_dense_gib = max_dense_gib,
        context = "RunCellRank(backend = \"cpp\")"
      )
    }
    if (identical(kernel_type, "wot")) {
      log_message(
        "{.arg backend = 'cpp'} does not support {.arg kernel_type = 'wot'}; use {.arg backend = 'python'} for Waddington-OT",
        message_type = "error"
      )
    }
    return(run_cellrank_cpp(
      srt = srt, assay_y = assay_y, layer_y = layer_y,
      group.by = group.by, linear_reduction = linear_reduction,
      nonlinear_reduction = nonlinear_reduction,
      n_pcs = n_pcs, n_neighbors = n_neighbors,
      mode = mode, kernel_type = kernel_type,
      velocity_weight = velocity_weight,
      connectivity_weight = connectivity_weight,
      use_connectivity_kernel = use_connectivity_kernel,
      softmax_scale = softmax_scale,
      n_macrostates = n_macrostates,
      n_cells_terminal = n_cells_terminal,
      estimator_type = tolower(estimator_type_upper),
      backward = backward,
      max_dense_gib = max_dense_gib,
      cores = cores, return_seurat = return_seurat,
      verbose = verbose
    ))
  }

  PrepareEnv(
    envname = envname,
    conda = conda,
    modules = c(
      "cellrank",
      if (kernel_type == "wot") "wot",
      if (isTRUE(magic_impute)) "magic"
    )
  )
  check_python("cellrank", envname = envname, conda = conda, verbose = verbose)
  if (kernel_type == "wot") {
    check_python("wot", envname = envname, conda = conda, verbose = verbose)
  }
  if (isTRUE(magic_impute)) {
    check_python("magic-impute", envname = envname, conda = conda, verbose = verbose)
  }
  if (all(is.null(srt), is.null(adata))) {
    log_message(
      "One of {.arg srt} or {.arg adata} must be provided",
      message_type = "error"
    )
  }
  if (is.null(group.by)) {
    log_message(
      "{.arg group.by} must be provided",
      message_type = "error"
    )
  }

  estimator_type <- match.arg(estimator_type)
  schur_method <- match.arg(schur_method)
  plot_format <- match.arg(plot_format)

  if (use_connectivity_kernel) {
    if (velocity_weight <= 0 && connectivity_weight <= 0) {
      velocity_weight <- 0.5
      connectivity_weight <- 0.5
      log_message(
        "Both kernel weights are non-positive; using equal weights",
        message_type = "warning",
        verbose = verbose
      )
    } else if (velocity_weight <= 0) {
      velocity_weight <- 0
      connectivity_weight <- 1
    } else if (connectivity_weight <= 0) {
      velocity_weight <- 1
      connectivity_weight <- 0
    }
    weight_sum <- velocity_weight + connectivity_weight
    if (abs(weight_sum - 1.0) > 0.01) {
      log_message(
        "Kernel weights ({velocity_weight} + {connectivity_weight} = {weight_sum}) do not sum to 1.0. Normalizing...",
        message_type = "warning",
        verbose = verbose
      )
      velocity_weight <- velocity_weight / weight_sum
      connectivity_weight <- connectivity_weight / weight_sum
    }
  }

  if (is.null(linear_reduction)) {
    linear_reduction <- DefaultReduction(srt)
  } else {
    linear_reduction <- DefaultReduction(srt, pattern = linear_reduction)
  }
  if (!linear_reduction %in% names(srt@reductions)) {
    log_message(
      "{.val {linear_reduction}} is not in the srt reduction names",
      message_type = "error"
    )
  }

  if (is.null(nonlinear_reduction)) {
    nonlinear_reduction <- DefaultReduction(srt)
  } else {
    nonlinear_reduction <- DefaultReduction(srt, pattern = nonlinear_reduction)
  }
  if (!nonlinear_reduction %in% names(srt@reductions)) {
    log_message(
      "{.val {nonlinear_reduction}} is not in the srt reduction names",
      message_type = "error"
    )
  }

  if (is.character(mode) && length(mode) == 1) {
    mode <- list(mode)
  }

  args <- mget(names(formals()))
  args <- lapply(
    args, function(x) {
      if (is.numeric(x)) {
        y <- ifelse(grepl("\\.", as.character(x)), as.double(x), as.integer(x))
      } else {
        y <- x
      }
      return(y)
    }
  )
  call_envir <- parent.frame(1)
  args <- lapply(args, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call_envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call_envir)
    } else {
      arg
    }
  })

  if (!is.null(srt)) {
    args[["adata"]] <- srt_to_adata(
      srt = srt,
      assay_x = assay_x,
      layer_x = layer_x,
      assay_y = assay_y,
      layer_y = layer_y,
      prepare_env = FALSE
    )
  }
  group_by_py <- group.by
  if ("group.by" %in% names(args)) {
    args[["group_by"]] <- args[["group.by"]]
    args[["group.by"]] <- NULL
  }
  groups <- py_to_r2(args[["adata"]]$obs)[[group_by_py]]
  args[["legend_loc"]] <- legend.position
  args[["n_jobs"]] <- cores
  args[["dpi"]] <- plot_dpi
  args[["fileprefix"]] <- plot_prefix
  args[["palette"]] <- palette_colors(
    levels(groups) %||% unique(groups),
    palette = palette,
    palcolor = palcolor
  )
  args <- args[
    !names(args) %in%
      c(
        "srt",
        "assay_x",
        "layer_x",
        "assay_y",
        "layer_y",
        "return_seurat",
        "palette",
        "palcolor",
        "legend.position",
        "cores",
        "plot_dpi",
        "plot_prefix",
        "backend",
        "backward",
        "allow_approximate",
        "max_dense_gib",
        "envname",
        "conda"
      )
  ]

  log_message("Running {.pkg CellRank} analysis...", verbose = verbose)
  functions <- scop_python_import("functions", convert = TRUE)
  result <- do.call(functions$CellRank, args)
  log_message(
    "{.pkg CellRank} analysis completed",
    message_type = "success",
    verbose = verbose
  )

  adata <- result[[1]]
  estimator <- result[[2]]
  kernel <- result[[3]]
  payload <- if (length(result) >= 4L) result[[4]] else list()

  payload_frame <- function(x) {
    if (is.null(x) || is.null(x$values)) return(NULL)
    values <- as.matrix(x$values)
    if (!is.null(x$index) && length(x$index) == nrow(values)) {
      rownames(values) <- as.character(x$index)
    }
    if (!is.null(x$columns) && length(x$columns) == ncol(values)) {
      colnames(values) <- as.character(x$columns)
    }
    values
  }
  payload_fate <- payload_frame(payload$fate_probabilities)
  payload_drivers <- payload_frame(payload$lineage_drivers)

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata, prepare_env = FALSE)

    # ── Normalize Python output to match C++ backend structure ──
    # Map Python obs column names → C++ standard names
    py_to_cpp_cols <- c(
      macrostates_fwd = "cellrank_macrostate",
      term_states_fwd = "cellrank_terminal_states"
    )
    missing_py_cols <- character(0)
    for (py_col in names(py_to_cpp_cols)) {
      cpp_col <- py_to_cpp_cols[[py_col]]
      if (py_col %in% colnames(srt_out@meta.data)) {
        vals <- srt_out@meta.data[[py_col]]
        vals_chr <- as.character(vals)
        empty_state <- vals_chr %in% c("transient", "unassigned", "NA", "nan", "")
        if (any(empty_state, na.rm = TRUE)) {
          non_empty <- vals_chr[!empty_state & !is.na(vals_chr)]
          if (length(non_empty) > 0) {
            state_ids <- as.integer(factor(non_empty))
            out_vals <- integer(length(vals_chr))
            out_vals[!empty_state & !is.na(vals_chr)] <- state_ids
            srt_out@meta.data[[cpp_col]] <- out_vals
          } else {
            srt_out@meta.data[[cpp_col]] <- integer(length(vals_chr))
          }
        } else if (is.factor(vals)) {
          srt_out@meta.data[[cpp_col]] <- as.integer(vals)
        } else if (is.numeric(vals) || is.integer(vals)) {
          srt_out@meta.data[[cpp_col]] <- as.integer(vals)
        } else {
          srt_out@meta.data[[cpp_col]] <- as.integer(factor(vals))
        }
      } else {
        missing_py_cols <- c(missing_py_cols, py_col)
      }
    }
    if (length(missing_py_cols) > 0) {
      log_message(
        "{.pkg CellRank} Python output did not include expected obs column{?s}: {.val {paste(missing_py_cols, collapse = ', ')}}",
        message_type = "warning",
        verbose = verbose
      )
    }

    # ── Extract absorption probabilities from adata.obsm ──
    ap <- NULL
    obsm_keys <- tryCatch(
      if (inherits(adata$obsm, "python.builtin.object")) {
        as.character(reticulate::iterate(adata$obsm$keys()))
      } else {
        names(adata$obsm)
      },
      error = function(e) character(0)
    )
    ap_key <- if ("to_terminal_states" %in% obsm_keys) {
      "to_terminal_states"
    } else if ("lineages_fwd" %in% obsm_keys) {
      "lineages_fwd"
    }
    if (is.null(payload_fate) && !is.null(ap_key)) {
      ap_raw <- py_to_r2(adata$obsm[[ap_key]])
      if (inherits(ap_raw, c("matrix", "Matrix", "dgCMatrix", "data.frame"))) {
        ap <- as.matrix(ap_raw)
        rownames(ap) <- colnames(srt_out)
      }
    }
    if (!is.null(payload_fate)) {
      ap <- as.matrix(payload_fate)
      if (is.null(rownames(ap)) && nrow(ap) == ncol(srt_out)) {
        rownames(ap) <- colnames(srt_out)
      }
    }
    if (!is.null(ap) && nrow(ap) == ncol(srt_out)) {
      ap <- ap[colnames(srt_out), , drop = FALSE]
      fate_names <- colnames(ap) %||% paste0("lineage_", seq_len(ncol(ap)))
      fate_meta_names <- paste0(
        "cellrank_fate_",
        make.names(fate_names, unique = TRUE)
      )
      for (i in seq_along(fate_meta_names)) {
        srt_out@meta.data[[fate_meta_names[[i]]]] <- as.numeric(ap[, i])
      }
      colnames(ap) <- fate_names
      if (ncol(ap) >= 2L && !"cellrank_circular" %in% names(srt_out@reductions)) {
        angles <- seq(0, 2 * pi, length.out = ncol(ap) + 1L)[-(ncol(ap) + 1L)]
        ap_norm <- ap / pmax(rowSums(ap), .Machine$double.eps)
        circular_embedding <- cbind(
          x = as.numeric(ap_norm %*% cos(angles)),
          y = as.numeric(ap_norm %*% sin(angles))
        )
        rownames(circular_embedding) <- rownames(ap_norm)
        colnames(circular_embedding) <- c("CRFATE_1", "CRFATE_2")
        srt_out[["cellrank_circular"]] <- SeuratObject::CreateDimReducObject(
          embeddings = circular_embedding,
          assay = SeuratObject::DefaultAssay(srt_out),
          key = "CRFATE_"
        )
      }
    }
    if ("term_states_fwd_probs" %in% colnames(srt_out@meta.data)) {
      srt_out@meta.data[["cellrank_terminal_membership_confidence"]] <-
        as.numeric(srt_out@meta.data[["term_states_fwd_probs"]])
    }
    fate_confidence_source <- "unavailable"
    if (!is.null(ap)) {
      srt_out@meta.data[["cellrank_fate_confidence"]] <-
        cellrank_fate_confidence_from_absorption(ap)
      fate_confidence_source <- if (!is.null(payload_fate)) {
        "normalized payload row maximum"
      } else {
        paste0(ap_key, " row maximum")
      }
    } else if ("term_states_fwd_probs" %in% colnames(srt_out@meta.data)) {
      # Preserve the historical column when fate probabilities are unavailable,
      # but record that this is only a terminal-membership fallback.
      srt_out@meta.data[["cellrank_fate_confidence"]] <-
        as.numeric(srt_out@meta.data[["term_states_fwd_probs"]])
      fate_confidence_source <- "term_states_fwd_probs fallback"
      log_message(
        "CellRank fate probabilities were unavailable; using terminal-state membership confidence as a fallback",
        message_type = "warning",
        verbose = verbose
      )
    }

    # ── Store standard tools$CellRank slot ──
    transition <- if ("cellrank_transition" %in% names(srt_out@graphs)) {
      srt_out@graphs[["cellrank_transition"]]
    } else {
      NULL
    }
    srt_out@tools[["CellRank"]] <- list(
      backend = "python",
      estimator = estimator_type,
      kernel = kernel_type,
      n_macrostates = if (!is.null(srt_out@meta.data[["cellrank_macrostate"]])) {
        length(unique(srt_out@meta.data[["cellrank_macrostate"]]))
      } else {
        0L
      },
      absorption_probabilities = ap,
      fate_probabilities = ap,
      lineage_drivers = payload_drivers,
      transition_matrix = transition,
      transition_key = payload$transition_key %||% "cellrank_transition",
      versions = payload$versions %||% list(),
      states = list(
        macrostates = payload$macrostates %||% character(),
        terminal_states = payload$terminal_states %||% character(),
        initial_states = payload$initial_states %||% character(),
        macrostate_names = payload$macrostate_names %||% character(),
        terminal_state_names = payload$terminal_state_names %||% character(),
        initial_state_names = payload$initial_state_names %||% character()
      ),
      parameters = list(
        kernel_type = kernel_type,
        estimator_type = estimator_type,
        n_macrostates = n_macrostates,
        n_cells_terminal = n_cells_terminal,
        backward = backward,
        fate_confidence_source = fate_confidence_source,
        schur_method = schur_method,
        schur_n_components = schur_n_components,
        terminal_states = terminal_states,
        terminal_state_agg = terminal_state_agg,
        driver_lineages = driver_lineages,
        compute_lineage_drivers = compute_lineage_drivers,
        actual = payload$parameters %||% list()
      )
    )

    # Keep raw Python objects in misc for debugging
    srt_out@misc$cellrank <- list(
      estimator = estimator,
      kernel = kernel
    )

    if (is.null(srt)) {
      return(srt_out)
    } else {
      merged <- srt_append(srt_raw = srt, srt_append = srt_out)

      if (is.null(merged@misc$cellrank)) {
        merged@misc$cellrank <- srt_out@misc$cellrank
      }
      # Merge tools$CellRank from append
      if (!is.null(srt_out@tools[["CellRank"]])) {
        merged@tools[["CellRank"]] <- srt_out@tools[["CellRank"]]
      }

      return(merged)
    }
  } else {
    adata$uns["cellrank_estimator_type"] <- estimator_type
    adata$uns["cellrank_kernel_type"] <- kernel_type
    return(adata)
  }
}

cellrank_fate_confidence_from_absorption <- function(absorption) {
  absorption <- as.matrix(absorption)
  if (!is.numeric(absorption) || length(dim(absorption)) != 2L || !ncol(absorption)) {
    log_message(
      "CellRank absorption probabilities must be a numeric matrix with at least one lineage",
      message_type = "error"
    )
  }
  if (any(!is.finite(absorption))) {
    log_message("CellRank absorption probabilities contain non-finite values", message_type = "error")
  }
  apply(absorption, 1L, max)
}

cellrank_normalize_sparse_transition <- function(
  transition,
  min_self_loop = 0.01
) {
  transition <- Matrix::drop0(Matrix::Matrix(transition, sparse = TRUE))
  if (nrow(transition) != ncol(transition)) {
    log_message("CellRank transition matrix must be square", message_type = "error")
  }
  if (length(transition@x)) {
    transition@x[!is.finite(transition@x) | transition@x < 0] <- 0
    transition <- Matrix::drop0(transition)
  }
  row_sums <- Matrix::rowSums(transition)
  zero_rows <- which(!is.finite(row_sums) | row_sums <= 1e-12)
  if (length(zero_rows)) {
    transition <- transition + Matrix::sparseMatrix(
      i = zero_rows,
      j = zero_rows,
      x = 1,
      dims = dim(transition)
    )
    row_sums <- Matrix::rowSums(transition)
  }
  transition <- Matrix::Diagonal(x = 1 / pmax(row_sums, 1e-12)) %*% transition

  diagonal <- Matrix::diag(transition)
  needs_loop <- which(!is.finite(diagonal) | diagonal < min_self_loop)
  if (length(needs_loop)) {
    transition <- transition + Matrix::sparseMatrix(
      i = needs_loop,
      j = needs_loop,
      x = min_self_loop - pmax(diagonal[needs_loop], 0),
      dims = dim(transition)
    )
  }
  row_sums <- Matrix::rowSums(transition)
  Matrix::drop0(Matrix::Diagonal(x = 1 / pmax(row_sums, 1e-12)) %*% transition)
}

cellrank_density_normalize_connectivities <- function(connectivities) {
  connectivities <- Matrix::drop0(Matrix::Matrix(connectivities, sparse = TRUE))
  if (nrow(connectivities) != ncol(connectivities)) {
    log_message("CellRank connectivity matrix must be square", message_type = "error")
  }
  density <- Matrix::colSums(connectivities)
  density[!is.finite(density) | density <= 1e-15] <- 1
  density_inverse <- Matrix::Diagonal(x = 1 / density)
  cellrank_normalize_sparse_transition(
    density_inverse %*% connectivities %*% density_inverse,
    min_self_loop = 0
  )
}

cellrank_hard_threshold_kernel <- function(
  connectivities,
  pseudotime,
  frac_to_keep = 0.3,
  backward = FALSE
) {
  connectivities <- Matrix::drop0(Matrix::Matrix(connectivities, sparse = TRUE))
  if (nrow(connectivities) != ncol(connectivities) || length(pseudotime) != nrow(connectivities)) {
    log_message(
      "CellRank connectivities and pseudotime dimensions do not agree",
      message_type = "error"
    )
  }
  if (!is.numeric(frac_to_keep) || length(frac_to_keep) != 1L ||
    !is.finite(frac_to_keep) || frac_to_keep < 0 || frac_to_keep > 1) {
    log_message("{.arg frac_to_keep} must be between 0 and 1", message_type = "error")
  }
  pseudotime <- as.numeric(pseudotime)
  if (any(!is.finite(pseudotime))) {
    log_message("CellRank pseudotime contains non-finite values", message_type = "error")
  }
  connectivities <- methods::as(
    methods::as(connectivities, "generalMatrix"),
    "CsparseMatrix"
  )
  cellrank_hard_threshold_kernel_cpp(
    connectivities = connectivities,
    pseudotime = pseudotime,
    frac_to_keep = frac_to_keep,
    backward = isTRUE(backward)
  )
}

run_cellrank_cpp <- function(
  srt, assay_y, layer_y, group.by,
  linear_reduction, nonlinear_reduction,
  n_pcs, n_neighbors, mode, kernel_type,
  velocity_weight, connectivity_weight, use_connectivity_kernel,
  softmax_scale, n_macrostates, n_cells_terminal,
  estimator_type = c("gpcca", "cflare"),
  backward = FALSE,
  max_dense_gib = 8,
  cores, return_seurat, verbose
) {
  estimator_type <- match.arg(estimator_type)
  if (isTRUE(use_connectivity_kernel)) {
    if (velocity_weight <= 0 && connectivity_weight <= 0) {
      velocity_weight <- 0.5
      connectivity_weight <- 0.5
    } else if (velocity_weight <= 0) {
      velocity_weight <- 0
      connectivity_weight <- 1
    } else if (connectivity_weight <= 0) {
      velocity_weight <- 1
      connectivity_weight <- 0
    } else {
      weight_sum <- velocity_weight + connectivity_weight
      velocity_weight <- velocity_weight / weight_sum
      connectivity_weight <- connectivity_weight / weight_sum
    }
  }
  if (is.null(srt)) {
    log_message("{.arg backend = 'cpp'} requires {.arg srt}", message_type = "error")
  }
  if (!isTRUE(return_seurat)) {
    log_message("{.arg backend = 'cpp'} returns a {.cls Seurat} object only", message_type = "error")
  }
  if (is.null(linear_reduction)) {
    linear_reduction <- DefaultReduction(srt)
  } else {
    linear_reduction <- DefaultReduction(srt, pattern = linear_reduction)
  }
  if (is.null(nonlinear_reduction)) {
    nonlinear_reduction <- DefaultReduction(srt)
  } else {
    nonlinear_reduction <- DefaultReduction(srt, pattern = nonlinear_reduction)
  }
  cells <- colnames(srt)
  nonlinear_embedding <- as.matrix(srt@reductions[[nonlinear_reduction]]@cell.embeddings[cells, , drop = FALSE])
  storage.mode(nonlinear_embedding) <- "double"
  n_cells <- nrow(nonlinear_embedding)
  n_dims <- ncol(nonlinear_embedding)
  dense_estimate_gib <- assert_cpp_dense_budget(
    n_rows = n_cells,
    n_cols = n_cells,
    copies = 8,
    max_dense_gib = max_dense_gib,
    context = "RunCellRank(backend = \"cpp\")"
  )

  velocity_reduction <- paste0(mode, "_", nonlinear_reduction)
  velocity_reduction_fallback <- if (!identical(mode, "velocity")) {
    paste0("velocity", "_", nonlinear_reduction)
  } else {
    NULL
  }
  if (identical(kernel_type, "velocity") &&
    !velocity_reduction %in% names(srt@reductions) &&
    !is.null(velocity_reduction_fallback) &&
    velocity_reduction_fallback %in% names(srt@reductions)) {
    log_message(
      "Using {.val {velocity_reduction_fallback}} as the velocity reduction (e.g. an object converted from scVelo)",
      message_type = "warning",
      verbose = verbose
    )
    velocity_reduction <- velocity_reduction_fallback
  }
  pt_key <- paste0(mode, "_pseudotime")
  graph_connectivities <- if ("connectivities" %in% names(srt@graphs)) {
    srt@graphs[["connectivities"]][cells, cells, drop = FALSE]
  } else {
    NULL
  }
  has_pseudotime <- pt_key %in% colnames(srt@meta.data) ||
    "dpt_pseudotime" %in% colnames(srt@meta.data)
  needs_velocity <- identical(kernel_type, "velocity") ||
    (identical(kernel_type, "pseudotime") && !has_pseudotime)
  if (isTRUE(needs_velocity) && !velocity_reduction %in% names(srt@reductions)) {
    log_message(
      "Running {.fn RunSCVELO} cpp backend before {.fn RunCellRank}",
      verbose = verbose
    )
    srt <- run_scanpy_cpp(
      srt = srt, assay_y = assay_y, layer_y = layer_y, group.by = group.by,
      linear_reduction = linear_reduction, nonlinear_reduction = nonlinear_reduction,
      n_pcs = n_pcs, n_neighbors = n_neighbors, mode = mode,
      filter_genes = TRUE, normalize_per_cell = TRUE, log_transform = TRUE,
      compute_terminal_states = FALSE, compute_pseudotime = FALSE,
      compute_velocity_confidence = FALSE,
      cores = cores, return_seurat = TRUE, verbose = verbose
    )
  } else if (isTRUE(needs_velocity)) {
    log_message(
      "Reusing existing {.field {velocity_reduction}} reduction for {.fn RunCellRank}",
      verbose = verbose
    )
  }
  ve <- NULL
  if (isTRUE(needs_velocity) || identical(kernel_type, "velocity")) {
    if (!velocity_reduction %in% names(srt@reductions)) {
      log_message(
        "Velocity reduction {.val {velocity_reduction}} is not available",
        message_type = "error"
      )
    }
    ve <- as.matrix(srt@reductions[[velocity_reduction]]@cell.embeddings)
    storage.mode(ve) <- "double"
  }
  le <- as.matrix(srt@reductions[[linear_reduction]]@cell.embeddings[
    cells, seq_len(min(n_pcs, ncol(srt@reductions[[linear_reduction]]@cell.embeddings))),
    drop = FALSE
  ])
  storage.mode(le) <- "double"
  knn_k <- max(1L, min(as.integer(n_neighbors) - 1L, n_cells - 1L))
  knn <- run_biocneighbors_knn(
    reference = le,
    k = knn_k,
    metric = "euclidean",
    exclude_self = TRUE,
    n_threads = as.integer(cores)
  )
  graph_transition <- NULL
  if (isTRUE(use_connectivity_kernel) && !is.null(graph_connectivities)) {
    graph_transition <- as.matrix(
      cellrank_density_normalize_connectivities(graph_connectivities)
    )
    storage.mode(graph_transition) <- "double"
  }

  combine_with_connectivity <- function(main_transition) {
    if (!isTRUE(use_connectivity_kernel) || connectivity_weight <= 0) {
      return(list(transition = main_transition, combined = FALSE))
    }
    connectivity_transition <- graph_transition
    if (is.null(connectivity_transition)) {
      connectivity_transition <- cellrank_connectivity_kernel_cpp(
        knn_idx = knn[["idx"]],
        knn_dist = knn[["dist"]]
      )
    }
    if (velocity_weight <= 0) {
      return(list(transition = connectivity_transition, combined = TRUE))
    }
    weight_sum <- velocity_weight + connectivity_weight
    list(
      transition =
        (velocity_weight / weight_sum) * main_transition +
          (connectivity_weight / weight_sum) * connectivity_transition,
      combined = TRUE
    )
  }

  # Build transition matrix based on kernel_type
  T_mat <- NULL
  kernel_used <- kernel_type
  pseudotime_source <- NULL

  if (kernel_type == "velocity") {
    main_transition <- cellrank_velocity_kernel_cpp(
      velocity_embedding = ve,
      embedding = nonlinear_embedding,
      knn_idx = knn[["idx"]],
      backward = isTRUE(backward),
      softmax_scale = softmax_scale
    )
    combined <- combine_with_connectivity(main_transition)
    T_mat <- combined$transition
    if (isTRUE(combined$combined)) kernel_used <- "velocity_connectivity_combined"
  } else if (kernel_type == "pseudotime") {
    pt_use <- pt_key
    if (!pt_use %in% colnames(srt@meta.data) && "dpt_pseudotime" %in% colnames(srt@meta.data)) {
      pt_use <- "dpt_pseudotime"
      log_message(
        "Using {.val {pt_use}} as the pseudotime (e.g. an object converted from scanpy)",
        message_type = "warning",
        verbose = verbose
      )
    }
    if (pt_use %in% colnames(srt@meta.data)) {
      pseudotime <- as.numeric(srt@meta.data[[pt_use]])
    } else {
      log_message(
        "Pseudotime not found; computing velocity pseudotime",
        message_type = "warning",
        verbose = verbose
      )
      ts_result <- scanpy_terminal_states_cpp(
        velocity_embedding = ve, embedding = nonlinear_embedding,
        knn_idx = knn[["idx"]], n_neighbors_velo = knn_k, seed = 0L
      )
      vpt_result <- scanpy_pseudotime_cpp(
        velocity_embedding = ve, embedding = nonlinear_embedding,
        knn_idx = knn[["idx"]], root_cells = ts_result[["root_cells"]],
        end_points = ts_result[["end_points"]],
        n_neighbors_velo = knn_k
      )
      pseudotime <- as.numeric(vpt_result[["pseudotime"]])
      pt_use <- "velocity_pseudotime_fallback"
    }
    pseudotime_source <- pt_use
    main_transition <- if (!is.null(graph_connectivities)) {
      as.matrix(cellrank_hard_threshold_kernel(
        connectivities = graph_connectivities,
        pseudotime = pseudotime,
        frac_to_keep = 0.3,
        backward = isTRUE(backward)
      ))
    } else {
      cellrank_pseudotime_kernel_cpp(
        pseudotime = pseudotime,
        knn_idx = knn[["idx"]],
        cell_weights = rep(1, length(pseudotime)),
        backward = isTRUE(backward)
      )
    }
    combined <- combine_with_connectivity(main_transition)
    T_mat <- combined$transition
    if (isTRUE(combined$combined)) kernel_used <- "pseudotime_connectivity_combined"
  } else if (kernel_type == "cytotrace") {
    cytotrace_expr <- as.matrix(GetAssayData5(srt, assay = assay_y[[1L]], layer = layer_y))
    gene_counts <- colSums(cytotrace_expr > 0)
    main_transition <- cellrank_cytotrace_kernel_cpp(
      gene_counts = as.numeric(gene_counts),
      knn_idx = knn[["idx"]],
      backward = isTRUE(backward)
    )
    combined <- combine_with_connectivity(main_transition)
    T_mat <- combined$transition
    if (isTRUE(combined$combined)) kernel_used <- "cytotrace_connectivity_combined"
  } else {
    # Default: velocity kernel (original inline computation)
    T_mat <- matrix(0, n_cells, n_cells)
    for (i in seq_len(n_cells)) {
      vn <- sqrt(sum(ve[i, ]^2))
      if (vn < 1e-10) {
        T_mat[i, i] <- 1.0
        next
      }
      wsum <- 0
      for (k in seq_len(knn_k)) {
        nb <- knn[["idx"]][i, k]
        if (is.na(nb)) next
        nb <- as.integer(nb)
        if (nb < 1 || nb > n_cells || nb == i) next
        delta <- nonlinear_embedding[nb, ] - nonlinear_embedding[i, ]
        dn <- sqrt(sum(delta^2))
        if (dn < 1e-10) next
        cosine <- sum(ve[i, ] * delta) / (vn * dn)
        if (cosine > 0) {
          T_mat[i, nb] <- exp(cosine * softmax_scale)
          wsum <- wsum + T_mat[i, nb]
        }
      }
      if (wsum > 0) T_mat[i, ] <- T_mat[i, ] / wsum else T_mat[i, i] <- 1.0
    }
  }

  storage.mode(T_mat) <- "double"
  rownames(T_mat) <- cells
  colnames(T_mat) <- cells
  n_mac <- if (is.null(n_macrostates)) {
    n_schur <- if (n_cells < 100L) {
      5L
    } else if (n_cells < 500L) {
      8L
    } else if (n_cells < 2000L) {
      10L
    } else {
      15L
    }
    max(2L, n_schur - 2L)
  } else {
    as.integer(n_macrostates)
  }
  n_mac <- min(n_mac, n_cells)

  # Choose estimator
  if (identical(estimator_type, "cflare")) {
    result <- cellrank_cflare_cpp(T_ = T_mat, n_states = n_mac)
  } else {
    result <- cellrank_gpcca_cpp(T_ = T_mat, n_states = n_mac, n_cells_terminal = as.integer(n_cells_terminal))
  }
  srt[["cellrank_terminal_states"]] <- as.integer(result[["terminal_states"]])
  srt[["cellrank_fate_confidence"]] <- as.numeric(result[["fate_confidence"]])
  srt[["cellrank_macrostate"]] <- as.integer(result[["macrostate_assignment"]])
  if (!is.null(ve)) {
    colnames(ve) <- colnames(nonlinear_embedding)
    rownames(ve) <- cells
    srt[["cellrank_velocity"]] <- SeuratObject::CreateDimReducObject(
      embeddings = ve, assay = SeuratObject::DefaultAssay(srt), key = "CRVELO_"
    )
  }

  # Store absorption probabilities if available
  if ("absorption_probabilities" %in% names(result)) {
    ap <- result[["absorption_probabilities"]]
    rownames(ap) <- cells
    srt@tools[["CellRank"]]$absorption_probabilities <- ap
  }
  if ("lineage_assignment" %in% names(result)) {
    srt@tools[["CellRank"]]$lineage_assignment <- as.integer(result[["lineage_assignment"]])
  }
  if ("chi" %in% names(result)) {
    srt@tools[["CellRank"]]$chi <- result[["chi"]]
  }
  if ("coarse_transition" %in% names(result)) {
    srt@tools[["CellRank"]]$coarse_transition <- result[["coarse_transition"]]
  }
  if ("absorption_method" %in% names(result)) {
    srt@tools[["CellRank"]]$absorption_diagnostics <- list(
      method = result[["absorption_method"]],
      iterations = result[["absorption_iterations"]],
      residual = result[["absorption_residual"]],
      converged = result[["absorption_converged"]],
      n_terminal_states = result[["n_terminal_states"]],
      n_terminal_cells = result[["n_terminal_cells"]]
    )
  }

  srt@tools[["CellRank"]] <- c(srt@tools[["CellRank"]], list(
    backend = "cpp",
    implementation = list(
      exact_reference = FALSE,
      scope = "shared-graph transition kernel, approximate macrostates, and cell-level absorption",
      dense_working_set_gib_lower_bound = dense_estimate_gib
    ),
    estimator = paste0(estimator_type, "_approximation"),
    kernel = kernel_used,
    transition_matrix = Matrix::Matrix(T_mat, sparse = TRUE),
    transition_cells = cells,
    n_macrostates = result[["n_macrostates"]],
    eigenvalues = result[["eigenvalues"]],
    schur_vectors = result[["schur_vectors"]],
    stationary_distribution = result[["stationary_distribution"]],
    parameters = list(
      mode = mode, kernel_type = kernel_type, estimator_type = estimator_type,
      softmax_scale = softmax_scale, n_macrostates = n_mac,
      n_cells_terminal = as.integer(n_cells_terminal),
      backward = backward,
      pseudotime_source = pseudotime_source,
      max_dense_gib = max_dense_gib
    )
  ))
  log_message("{.pkg CellRank} cpp backend ({.val {estimator_type}} + {.val {kernel_used}} kernel) completed",
    message_type = "success", verbose = verbose
  )
  srt
}
