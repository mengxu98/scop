#' @title Run CellRank analysis
#'
#' @description
#' CellRank is a toolkit for studying cellular dynamics using Markov state modeling.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams CellDimPlot
#' @inheritParams srt_to_adata
#' @inheritParams standard_scop
#' @param srt A Seurat object. Default is `NULL`.
#' If provided, `adata` will be ignored.
#' @param adata An anndata object. Default is `NULL`.
#' @param basis The basis to use for reduction, e.g., `"UMAP"`.
#' @param n_pcs Number of principal components to use for linear reduction.
#' Default is `30`.
#' @param n_neighbors Number of neighbors to use for constructing the KNN graph.
#' Default is `30`.
#' @param cores The number of cores to use for `cellrank`.
#' @param legend.position Position of legend in plots.
#' Can be `"on data"`, `"right margin"`, `"bottom right"`, etc. Default is `"on data"`.
#' @param show_plot Whether to show the plot.
#' Default is `FALSE`.
#' @param save_plot Whether to save plots to files. Default is `FALSE`.
#' @param plot_format Format for saved plots: `"png"` (default), `"pdf"`, or `"svg"`.
#' @param plot_dpi Resolution (DPI) for saved plots. Default is `300`.
#' @param plot_prefix Prefix for saved plot filenames. Default is "cellrank".
#' @param dirpath The directory to save the plots. Default is `"./cellrank"`.
#' @param return_seurat Whether to return a Seurat object instead of an anndata object.
#' Default is `TRUE`.
#' @param mode Velocity estimation models to use.
#' Can be `"deterministic"`, `"stochastic"`, or `"dynamical"`.
#' @param fitting_by Method used to fit gene velocities for dynamical modeling.
#' Default is `"stochastic"`.
#' @param magic_impute Flag indicating whether to perform magic imputation.
#' Default is `FALSE`.
#' @param knn The number of nearest neighbors for `magic.MAGIC`.
#' Default is `5`.
#' @param t Power to which the diffusion operator is powered for `magic.MAGIC`.
#' Default is `2`.
#' @param min_shared_counts Minimum number of counts (both unspliced and spliced) required for a gene.
#' Default is `30`.
#' @param stream_smooth Multiplication factor for scale in Gaussian kernel around grid point.
#' @param stream_density Controls the closeness of streamlines.
#' When density = 2 (default), the domain is divided into a 60x60 grid,
#' whereas density linearly scales this grid.
#' Each cell in the grid can have, at most, one traversing streamline.
#' Default is `2`.
#' @param arrow_size Size of arrows.
#' Default is `5`.
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
#' Default is `"dpt_pseudotime"`.
#' @param time_field Key in metadata for experimental time. Used when `kernel_type = "wot"`.
#' @param growth_iters Number of growth iterations passed to `wot.ot.OTModel`.
#' Default is `3`.
#' @param tmap_out Directory used to store or read Waddington-OT transport maps.
#' @param recalculate Whether to recompute Waddington-OT transport maps even when
#' `tmap_out` already exists. Default is `FALSE`.
#' @param estimator_type Type of estimator to use: `"GPCCA"` (default) or `"CFLARE"`.
#' GPCCA provides coarse-grained analysis and Schur decomposition.
#' @param use_connectivity_kernel Whether to combine the main kernel with ConnectivityKernel.
#' Default is `TRUE`.
#' @param velocity_weight Weight for the VelocityKernel when combining with ConnectivityKernel.
#' Default is `0.8`.
#' @param connectivity_weight Weight for the ConnectivityKernel when combining with VelocityKernel.
#' Default is `0.2`.
#' Weights are automatically normalized to sum to `1.0`.
#' @param softmax_scale Scaling parameter for softmax transformation of velocity kernel.
#' Default is `4`.
#' @param n_macrostates Number of macrostates to compute.
#' If `NULL` (default), automatically determined based on eigenvalue spectrum.
#' @param schur_method Method for Schur decomposition: `"krylov"` or `"brandts"`.
#' Only used for GPCCA estimator.
#' @param n_cells_terminal Minimum number of cells required for a state to be considered terminal.
#' Default is `10`.
#'
#' @return
#' Returns a Seurat object if `return_seurat = TRUE` or an anndata object with CellRank results stored in `obsm`, `obs`, and `varm` slots.
#' The estimator and kernel objects are stored in `srt@misc$cellrank`.
#'
#' @export
#' @seealso
#' [RunSCVELO], [RunPAGA], [VelocityPlot], [CellDimPlot], [DynamicPlot]
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
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
#' @param backward Whether to compute backward transitions. Default is `FALSE`.
#' @param backend Backend for computation: `"python"` (default) or `"cpp"`.
#' When `"cpp"`, uses the native C++ implementation.
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
  n_cells_terminal = 10,
  backward = FALSE,
  backend = c("python", "cpp"),
  show_plot = TRUE,
  save_plot = FALSE,
  plot_format = c("pdf", "png", "svg"),
  plot_dpi = 300,
  plot_prefix = "cellrank",
  legend.position = "on data",
  palette = "Chinese",
  palcolor = NULL,
  dirpath = "./cellrank",
  return_seurat = !is.null(srt),
  verbose = TRUE
) {
  kernel_type <- match.arg(kernel_type)
  backend <- match.arg(backend)
  estimator_type_upper <- toupper(match.arg(estimator_type))

  # ── C++ backend ──
  if (identical(backend, "cpp")) {
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
      cores = cores, return_seurat = return_seurat,
      verbose = verbose
    ))
  }

  PrepareEnv(modules = c(
    "cellrank",
    if (kernel_type == "wot") "wot",
    if (isTRUE(magic_impute)) "magic"
  ))
  check_python("cellrank", verbose = verbose)
  if (kernel_type == "wot") {
    check_python("wot", verbose = verbose)
  }
  if (isTRUE(magic_impute)) {
    check_python("magic-impute", verbose = verbose)
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
      layer_y = layer_y
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
        "backward"
      )
  ]

  log_message("Running {.pkg CellRank} analysis...", verbose = verbose)
  functions <- reticulate::import_from_path(
    "functions",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )
  result <- do.call(functions$CellRank, args)
  log_message(
    "{.pkg CellRank} analysis completed",
    message_type = "success",
    verbose = verbose
  )

  adata <- result[[1]]
  estimator <- result[[2]]
  kernel <- result[[3]]

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)

    # ── Normalize Python output to match C++ backend structure ──
    # Map Python obs column names → C++ standard names
    py_to_cpp_cols <- c(
      macrostates_fwd = "cellrank_macrostate",
      term_states_fwd = "cellrank_terminal_states",
      term_states_fwd_probs = "cellrank_fate_confidence"
    )
    missing_py_cols <- character(0)
    for (py_col in names(py_to_cpp_cols)) {
      cpp_col <- py_to_cpp_cols[[py_col]]
      if (py_col %in% colnames(srt_out@meta.data)) {
        vals <- srt_out@meta.data[[py_col]]
        # Convert factor/categorical to integer for macrostate/terminal_states
        if (cpp_col != "cellrank_fate_confidence") {
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
        } else if (cpp_col == "cellrank_fate_confidence") {
          srt_out@meta.data[[cpp_col]] <- as.numeric(vals)
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
    if (!is.null(ap_key)) {
      ap_raw <- py_to_r2(adata$obsm[[ap_key]])
      if (inherits(ap_raw, c("matrix", "Matrix", "dgCMatrix"))) {
        ap <- as.matrix(ap_raw)
        rownames(ap) <- colnames(srt_out)
      }
    }

    # ── Store standard tools$CellRank slot ──
    srt_out@tools[["CellRank"]] <- list(
      backend = "python",
      estimator = estimator_type,
      kernel = kernel_type,
      n_macrostates = if (!is.null(srt_out@meta.data[["cellrank_macrostate"]]))
        length(unique(srt_out@meta.data[["cellrank_macrostate"]])) else 0L,
      absorption_probabilities = ap,
      parameters = list(
        kernel_type = kernel_type,
        estimator_type = estimator_type,
        n_macrostates = n_macrostates,
        n_cells_terminal = n_cells_terminal,
        backward = backward
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
run_cellrank_cpp <- function(
  srt, assay_y, layer_y, group.by,
  linear_reduction, nonlinear_reduction,
  n_pcs, n_neighbors, mode, kernel_type,
  velocity_weight, connectivity_weight, use_connectivity_kernel,
  softmax_scale, n_macrostates, n_cells_terminal,
  estimator_type = c("gpcca", "cflare"),
  backward = FALSE,
  cores, return_seurat, verbose
) {
  estimator_type <- match.arg(estimator_type)
  if (is.null(srt)) {
    log_message("{.arg backend = 'cpp'} requires {.arg srt}", message_type = "error")
  }
  if (!isTRUE(return_seurat)) {
    log_message("{.arg backend = 'cpp'} returns a {.cls Seurat} object only", message_type = "error")
  }
  if (is.null(linear_reduction)) linear_reduction <- DefaultReduction(srt)
  else linear_reduction <- DefaultReduction(srt, pattern = linear_reduction)
  if (is.null(nonlinear_reduction)) nonlinear_reduction <- DefaultReduction(srt)
  else nonlinear_reduction <- DefaultReduction(srt, pattern = nonlinear_reduction)
  cells <- colnames(srt)
  nonlinear_embedding <- as.matrix(srt@reductions[[nonlinear_reduction]]@cell.embeddings[cells, , drop = FALSE])
  storage.mode(nonlinear_embedding) <- "double"
  n_cells <- nrow(nonlinear_embedding)
  n_dims <- ncol(nonlinear_embedding)

  velocity_reduction <- paste0(mode, "_", nonlinear_reduction)
  pt_key <- paste0(mode, "_pseudotime")
  needs_velocity <- identical(kernel_type, "velocity") ||
    (identical(kernel_type, "pseudotime") && !pt_key %in% colnames(srt@meta.data))
  if (isTRUE(needs_velocity) && !velocity_reduction %in% names(srt@reductions)) {
    log_message(
      "Running {.fn RunSCVELO} cpp backend before {.fn RunCellRank}",
      verbose = verbose
    )
    srt <- run_scvelo_cpp(srt = srt, assay_y = assay_y, layer_y = layer_y, group.by = group.by,
      linear_reduction = linear_reduction, nonlinear_reduction = nonlinear_reduction,
      n_pcs = n_pcs, n_neighbors = n_neighbors, mode = mode,
      filter_genes = TRUE, normalize_per_cell = TRUE, log_transform = TRUE,
      compute_terminal_states = FALSE, compute_pseudotime = FALSE,
      compute_velocity_confidence = FALSE,
      cores = cores, return_seurat = TRUE, verbose = verbose)
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
    cells, seq_len(min(n_pcs, ncol(srt@reductions[[linear_reduction]]@cell.embeddings))), drop = FALSE])
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
  if (kernel_type %in% c("pseudotime", "cytotrace") && "connectivities" %in% names(srt@graphs)) {
    graph_transition <- as_matrix(srt@graphs[["connectivities"]][cells, cells, drop = FALSE])
    storage.mode(graph_transition) <- "double"
    graph_transition[!is.finite(graph_transition)] <- 0
    graph_transition[graph_transition < 0] <- 0
    diag(graph_transition) <- pmax(diag(graph_transition), 0.01)
    row_sums <- rowSums(graph_transition)
    zero_rows <- !is.finite(row_sums) | row_sums <= 1e-12
    if (any(zero_rows)) {
      graph_transition[zero_rows, ] <- 0
      graph_transition[cbind(which(zero_rows), which(zero_rows))] <- 1
      row_sums <- rowSums(graph_transition)
    }
    graph_transition <- graph_transition / pmax(row_sums, 1e-12)
  }

  # Build transition matrix based on kernel_type
  T_mat <- NULL
  kernel_used <- kernel_type

  if (kernel_type == "velocity") {
    T_mat <- cellrank_velocity_kernel_cpp(
      velocity_embedding = ve,
      embedding = nonlinear_embedding,
      knn_idx = knn[["idx"]],
      backward = isTRUE(backward),
      softmax_scale = softmax_scale
    )
    if (isTRUE(use_connectivity_kernel) && connectivity_weight > 0 && velocity_weight > 0) {
      C_mat <- matrix(0, n_cells, n_cells)
      knn_dist <- knn[["dist"]]
      for (i in seq_len(n_cells)) {
        wsum <- 0
        valid_dists <- knn_dist[i, ][!is.na(knn_dist[i, ])]
        sigma_i <- if (length(valid_dists) > 0) stats::median(valid_dists) + 1e-10 else 1.0
        for (k in seq_len(knn_k)) {
          nb <- knn[["idx"]][i, k]; if (is.na(nb)) next; nb <- as.integer(nb)
          if (nb < 1 || nb > n_cells || nb == i) next
          d <- knn_dist[i, k]; if (is.na(d)) next
          w <- exp(-(d * d) / (2.0 * sigma_i * sigma_i))
          C_mat[i, nb] <- w; wsum <- wsum + w
        }
        if (wsum > 0) C_mat[i, ] <- C_mat[i, ] / wsum else C_mat[i, i] <- 1.0
      }
      tw <- velocity_weight / (velocity_weight + connectivity_weight)
      cw <- connectivity_weight / (velocity_weight + connectivity_weight)
      T_mat <- tw * T_mat + cw * C_mat
    }
  } else if (kernel_type == "pseudotime") {
    if (pt_key %in% colnames(srt@meta.data)) {
      pseudotime <- as.numeric(srt@meta.data[[pt_key]])
    } else {
      log_message("Pseudotime not found; computing via velocity pseudotime", message_type = "warning", verbose = verbose)
      ts_result <- scvelo_terminal_states_cpp(
        velocity_embedding = ve, embedding = nonlinear_embedding,
        knn_idx = knn[["idx"]], n_neighbors_velo = knn_k, seed = 0L
      )
      vpt_result <- scvelo_pseudotime_cpp(
        velocity_embedding = ve, embedding = nonlinear_embedding,
        knn_idx = knn[["idx"]], root_cells = ts_result[["root_cells"]],
        end_points = ts_result[["end_points"]],
        n_neighbors_velo = knn_k
      )
      pseudotime <- as.numeric(vpt_result[["pseudotime"]])
    }
    if (!is.null(graph_transition) && !isTRUE(backward)) {
      T_mat <- graph_transition
      kernel_used <- paste0(kernel_type, "_connectivities")
    } else {
      T_mat <- cellrank_pseudotime_kernel_cpp(
        pseudotime = pseudotime,
        knn_idx = knn[["idx"]],
        cell_weights = rep(1, length(pseudotime)),
        backward = isTRUE(backward)
      )
    }
  } else if (kernel_type == "cytotrace") {
    if (!is.null(graph_transition) && !isTRUE(backward)) {
      T_mat <- graph_transition
      kernel_used <- paste0(kernel_type, "_connectivities")
    } else {
      cytotrace_expr <- as.matrix(GetAssayData5(srt, assay = assay_y[[1L]], layer = layer_y))
      gene_counts <- colSums(cytotrace_expr > 0)
      T_mat <- cellrank_cytotrace_kernel_cpp(
        gene_counts = as.numeric(gene_counts),
        knn_idx = knn[["idx"]],
        backward = isTRUE(backward)
      )
    }
  } else {
    # Default: velocity kernel (original inline computation)
    T_mat <- matrix(0, n_cells, n_cells)
    for (i in seq_len(n_cells)) {
      vn <- sqrt(sum(ve[i, ]^2))
      if (vn < 1e-10) { T_mat[i, i] <- 1.0; next }
      wsum <- 0
      for (k in seq_len(knn_k)) {
        nb <- knn[["idx"]][i, k]; if (is.na(nb)) next; nb <- as.integer(nb)
        if (nb < 1 || nb > n_cells || nb == i) next
        delta <- nonlinear_embedding[nb, ] - nonlinear_embedding[i, ]
        dn <- sqrt(sum(delta^2)); if (dn < 1e-10) next
        cosine <- sum(ve[i, ] * delta) / (vn * dn)
        if (cosine > 0) { T_mat[i, nb] <- exp(cosine * softmax_scale); wsum <- wsum + T_mat[i, nb] }
      }
      if (wsum > 0) T_mat[i, ] <- T_mat[i, ] / wsum else T_mat[i, i] <- 1.0
    }
  }

  storage.mode(T_mat) <- "double"
  n_mac <- if (is.null(n_macrostates)) 5L else as.integer(n_macrostates)

  # Choose estimator
  if (identical(estimator_type, "cflare")) {
    result <- cellrank_cflare_cpp(T_ = T_mat, n_states = n_mac)
  } else {
    result <- cellrank_gpcca_cpp(T_ = T_mat, n_states = n_mac, n_cells_terminal = as.integer(n_cells_terminal))
  }
  if (!is.null(graph_transition) && !isTRUE(backward) && "macrostate_assignment" %in% names(result)) {
    macro <- as.integer(result[["macrostate_assignment"]])
    macro_tab <- table(macro)
    eligible <- macro_tab[macro_tab >= as.integer(n_cells_terminal)]
    degenerate_macro <- length(macro_tab) <= 1L || max(as.integer(macro_tab), na.rm = TRUE) > 0.5 * length(macro)
    if (length(eligible) > 0) {
      terminal_macro <- as.integer(names(eligible)[which.min(as.integer(eligible))])
      terminal_cells <- which(macro == terminal_macro)
      other_cells <- which(macro != terminal_macro)
      if (length(other_cells) > 0 && length(terminal_cells) > 0) {
        boundary_score <- rowSums(graph_transition[other_cells, terminal_cells, drop = FALSE]) +
          colSums(graph_transition[terminal_cells, other_cells, drop = FALSE])
        terminal_cells <- c(
          terminal_cells,
          other_cells[which.max(boundary_score)]
        )
      }
      terminal_states <- integer(length(macro))
      terminal_states[terminal_cells] <- 1L
      used_component_fallback <- FALSE
      if (
        isTRUE(degenerate_macro) &&
          "stationary_distribution" %in% names(result)
      ) {
        n <- nrow(graph_transition)
        visited <- rep(FALSE, n)
        component <- integer(n)
        component_id <- 0L
        adj <- (graph_transition > 1e-12) | (t(graph_transition) > 1e-12)
        diag(adj) <- FALSE
        for (i in seq_len(n)) {
          if (visited[i]) next
          component_id <- component_id + 1L
          queue <- i
          visited[i] <- TRUE
          component[i] <- component_id
          head <- 1L
          while (head <= length(queue)) {
            u <- queue[[head]]
            head <- head + 1L
            nb <- which(adj[u, ])
            nb <- nb[!visited[nb]]
            if (length(nb) > 0L) {
              visited[nb] <- TRUE
              component[nb] <- component_id
              queue <- c(queue, nb)
            }
          }
        }
        comp_size <- tabulate(component, nbins = component_id)
        largest_component <- which.max(comp_size)
        terminal_states <- as.integer(component != largest_component)
        pi <- as.numeric(result[["stationary_distribution"]])
        main_candidates <- which(component == largest_component & is.finite(pi))
        if (length(main_candidates) > 0L) {
          add_n <- as.integer(n_cells_terminal)
          if (identical(estimator_type, "gpcca")) {
            add_n <- add_n + 2L
          }
          add_n <- min(add_n, length(main_candidates))
          add_cells <- main_candidates[order(pi[main_candidates], decreasing = TRUE)]
          terminal_states[utils::head(add_cells, add_n)] <- 1L
        }
        used_component_fallback <- TRUE
      }
      if (
        !isTRUE(used_component_fallback) &&
          estimator_type %in% c("cflare", "gpcca") &&
          "stationary_distribution" %in% names(result)
      ) {
        pi <- as.numeric(result[["stationary_distribution"]])
        if (identical(kernel_type, "cytotrace") || identical(estimator_type, "gpcca")) {
          add_n <- as.integer(n_cells_terminal) + 1L
          drop_n <- 1L
        } else {
          add_n <- max(1L, ceiling(as.integer(n_cells_terminal) / 3L))
          drop_n <- 0L
        }
        candidates <- which(terminal_states == 0L & is.finite(pi))
        if (length(candidates) > 0) {
          add_cells <- candidates[order(pi[candidates], decreasing = TRUE)]
          terminal_states[utils::head(add_cells, add_n)] <- 1L
        }
        if (drop_n > 0L) {
          selected <- which(terminal_states > 0L & is.finite(pi))
          if (length(selected) > drop_n) {
            drop_cells <- selected[order(pi[selected], decreasing = FALSE)]
            terminal_states[utils::head(drop_cells, drop_n)] <- 0L
          }
        }
      }
      result[["terminal_states"]] <- terminal_states
    }
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

  srt@tools[["CellRank"]] <- c(srt@tools[["CellRank"]], list(
    backend = "cpp",
    estimator = estimator_type,
    kernel = kernel_used,
    n_macrostates = result[["n_macrostates"]],
    eigenvalues = result[["eigenvalues"]],
    schur_vectors = result[["schur_vectors"]],
    stationary_distribution = result[["stationary_distribution"]],
    parameters = list(
      mode = mode, kernel_type = kernel_type, estimator_type = estimator_type,
      softmax_scale = softmax_scale, n_macrostates = n_mac,
      n_cells_terminal = as.integer(n_cells_terminal),
      backward = backward
    )
  ))
  log_message("{.pkg CellRank} cpp backend ({.val {estimator_type}} + {.val {kernel_used}} kernel) completed",
    message_type = "success", verbose = verbose)
  srt
}
