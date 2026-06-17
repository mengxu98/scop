#' @title Run scVelo workflow
#'
#' @description
#' scVelo is a scalable toolkit for RNA velocity analysis in single cells.
#' This function runs an enhanced scVelo workflow on a Seurat object with improved error handling,
#' version compatibility, and modular design.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams RunCellRank
#' @param backend Backend used to compute RNA velocity. `"python"` keeps the
#' original scVelo workflow. `"cpp"` uses the package C++ implementation for a
#' stochastic velocity embedding that is compatible with [VelocityPlot].
#' @param filter_genes Whether to filter genes based on minimum counts.
#' @param min_counts Minimum counts for gene filtering.
#' @param min_counts_u Minimum unspliced counts for gene filtering.
#' @param normalize_per_cell Whether to normalize counts per cell.
#' @param log_transform Whether to apply log transformation.
#' @param use_raw Whether to use raw data for dynamical modeling.
#' @param diff_kinetics Whether to use differential kinetics.
#' @param denoise_topn Number of genes with highest likelihood selected to infer velocity directions.
#' @param kinetics_topn Number of genes with highest likelihood selected to infer velocity directions.
#' @param compute_velocity_confidence Whether to compute velocity confidence metrics.
#' @param compute_terminal_states Whether to compute terminal states (root and end points).
#' @param compute_pseudotime Whether to compute velocity pseudotime.
#' @param compute_paga Whether to compute PAGA (Partition-based graph abstraction).
#' @param top_n The number of top features to plot.
#'
#' @seealso
#' [VelocityPlot], [CellDimPlot], [RunPAGA]
#'
#' @export
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunSCVELO(
#'   pancreas_sub,
#'   assay_x = "RNA",
#'   group.by = "SubCellType",
#'   linear_reduction = "PCA",
#'   nonlinear_reduction = "UMAP",
#'   backend = "cpp"
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   c(
#'     "stochastic_length",
#'     "stochastic_confidence"
#'   )
#' )
#'
#' VelocityPlot(
#'   pancreas_sub,
#'   reduction = "UMAP",
#'   plot_type = "stream"
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   pt.size = NA,
#'   velocity = "stochastic"
#' )
RunSCVELO <- function(
  srt = NULL,
  adata = NULL,
  assay_x = "RNA",
  layer_x = "counts",
  assay_y = c("spliced", "unspliced"),
  layer_y = "counts",
  group.by = NULL,
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
  filter_genes = TRUE,
  min_counts = 3,
  min_counts_u = 3,
  normalize_per_cell = TRUE,
  log_transform = TRUE,
  use_raw = FALSE,
  diff_kinetics = FALSE,
  stream_smooth = NULL,
  stream_density = 2,
  arrow_length = 5,
  arrow_size = 5,
  arrow_density = 0.5,
  denoise = FALSE,
  denoise_topn = 3,
  kinetics = FALSE,
  kinetics_topn = 100,
  calculate_velocity_genes = FALSE,
  compute_velocity_confidence = TRUE,
  compute_terminal_states = TRUE,
  compute_pseudotime = TRUE,
  compute_paga = TRUE,
  top_n = 6,
  cores = 1,
  palette = "Chinese",
  palcolor = NULL,
  legend.position = "on data",
  show_plot = TRUE,
  save_plot = FALSE,
  plot_format = c("pdf", "png", "svg"),
  plot_dpi = 300,
  plot_prefix = "scvelo",
  dirpath = "./scvelo",
  backend = c("python", "cpp"),
  return_seurat = !is.null(srt),
  verbose = TRUE
) {
  backend <- match.arg(backend)
  plot_format <- match.arg(plot_format)

  if (identical(backend, "cpp")) {
    return(run_scvelo_cpp(
      srt = srt,
      assay_y = assay_y,
      layer_y = layer_y,
      group.by = group.by,
      linear_reduction = linear_reduction,
      nonlinear_reduction = nonlinear_reduction,
      n_pcs = n_pcs,
      n_neighbors = n_neighbors,
      mode = mode,
      filter_genes = filter_genes,
      min_counts = min_counts,
      min_counts_u = min_counts_u,
      min_shared_counts = min_shared_counts,
      compute_terminal_states = compute_terminal_states,
      compute_pseudotime = compute_pseudotime,
      compute_velocity_confidence = compute_velocity_confidence,
      fitting_by = fitting_by,
      cores = cores,
      return_seurat = return_seurat,
      verbose = verbose
    ))
  }

  PrepareEnv(modules = c(
    "scvelo",
    if (isTRUE(magic_impute)) "magic"
  ))

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
  args <- lapply(args, function(x) {
    if (is.numeric(x)) {
      y <- ifelse(grepl("\\.", as.character(x)), as.double(x), as.integer(x))
    } else {
      y <- x
    }
    return(y)
  })

  call.envir <- parent.frame(1)
  args <- lapply(args, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call.envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call.envir)
    } else {
      arg
    }
  })

  args[["n_jobs"]] <- cores

  args[["legend_loc"]] <- legend.position

  args[["dpi"]] <- plot_dpi
  args[["fileprefix"]] <- plot_prefix

  params <- c(
    "srt",
    "assay_x",
    "layer_x",
    "assay_y",
    "layer_y",
    "return_seurat",
    "palette",
    "palcolor",
    "cores",
    "legend.position",
    "plot_dpi",
    "plot_prefix",
    "backend"
  )
  args <- args[!names(args) %in% params]

  if (!is.null(srt)) {
    old_skip_python_prepare <- getOption("scop_skip_python_prepare", FALSE)
    options(scop_skip_python_prepare = TRUE)
    on.exit(
      options(scop_skip_python_prepare = old_skip_python_prepare),
      add = TRUE
    )
    args[["adata"]] <- srt_to_adata(
      srt = srt,
      assay_x = assay_x,
      layer_x = layer_x,
      assay_y = assay_y,
      layer_y = layer_y,
      verbose = verbose
    )
  }

  if ("group.by" %in% names(args)) {
    args[["group_by"]] <- args[["group.by"]]
    args[["group.by"]] <- NULL
  }
  groups <- py_to_r2(args[["adata"]]$obs)[[group.by]]
  args[["palette"]] <- palette_colors(
    levels(groups) %||% unique(groups),
    palette = palette,
    palcolor = palcolor
  )

  functions <- reticulate::import_from_path(
    "functions",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )

  adata <- do.call(functions$SCVELO, args)

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      srt_out1 <- srt_append(srt_raw = srt, srt_append = srt_out)
      srt_out2 <- srt_append(
        srt_raw = srt_out1,
        srt_append = srt_out,
        pattern = paste0(
          "(velocity)|(distances)|(connectivities)|(Ms)|(Mu)|(",
          paste(mode, collapse = ")|("),
          ")|(paga)"
        ),
        overwrite = TRUE,
        verbose = FALSE
      )
      return(srt_out2)
    }
  } else {
    return(adata)
  }
}

run_scvelo_cpp <- function(
  srt,
  assay_y,
  layer_y,
  group.by,
  linear_reduction,
  nonlinear_reduction,
  n_pcs,
  n_neighbors,
  mode,
  filter_genes = TRUE,
  min_counts = 3,
  min_counts_u = 3,
  min_shared_counts = 30,
  normalize_per_cell = TRUE,
  log_transform = TRUE,
  compute_terminal_states = TRUE,
  compute_pseudotime = TRUE,
  compute_velocity_confidence = FALSE,
  fitting_by = c("stochastic", "deterministic", "em", "nm"),
  cores,
  return_seurat,
  verbose = TRUE
) {
  if (is.null(srt)) {
    log_message(
      "{.arg backend = 'cpp'} requires {.arg srt}",
      message_type = "error"
    )
  }
  if (!isTRUE(return_seurat)) {
    log_message(
      "{.arg backend = 'cpp'} returns a {.cls Seurat} object only",
      message_type = "error"
    )
  }
  if (is.character(mode) && length(mode) == 1L) {
    mode <- list(mode)
  }
  mode_use <- unlist(mode, use.names = FALSE)
  if (!all(mode_use %in% c("stochastic", "deterministic", "dynamical"))) {
    log_message(
      "{.arg backend = 'cpp'} supports {.arg mode = 'stochastic'}, {.arg mode = 'deterministic'}, and {.arg mode = 'dynamical'}",
      message_type = "error"
    )
  }
  if (length(assay_y) < 2L) {
    log_message(
      "{.arg assay_y} must contain spliced and unspliced assay names",
      message_type = "error"
    )
  }
  spliced_assay <- assay_y[[1L]]
  unspliced_assay <- assay_y[[2L]]
  missing_assays <- setdiff(c(spliced_assay, unspliced_assay), names(srt@assays))
  if (length(missing_assays) > 0L) {
    log_message(
      "Missing velocity assays: {.val {missing_assays}}",
      message_type = "error"
    )
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
  if (!linear_reduction %in% names(srt@reductions)) {
    log_message(
      "{.val {linear_reduction}} is not in the srt reduction names",
      message_type = "error"
    )
  }
  if (!nonlinear_reduction %in% names(srt@reductions)) {
    log_message(
      "{.val {nonlinear_reduction}} is not in the srt reduction names",
      message_type = "error"
    )
  }

  spliced <- GetAssayData5(srt, assay = spliced_assay, layer = layer_y)
  unspliced <- GetAssayData5(srt, assay = unspliced_assay, layer = layer_y)
  features <- intersect(rownames(spliced), rownames(unspliced))
  if (length(features) == 0L) {
    log_message(
      "No shared features are available for {.fn RunSCVELO} cpp backend",
      message_type = "error"
    )
  }
  cells <- colnames(srt)
  spliced <- as.matrix(spliced[features, cells, drop = FALSE])
  unspliced <- as.matrix(unspliced[features, cells, drop = FALSE])
  storage.mode(spliced) <- "double"
  storage.mode(unspliced) <- "double"

  # ── scanpy-compatible filtering and normalization ──
  log_message(
    "Running scanpy-compatible preprocessing ({.val {length(features)}} features -> filter + normalize)...",
    verbose = verbose
  )
  initial_spliced_totals <- colSums(spliced)
  initial_unspliced_totals <- colSums(unspliced)
  keep <- scvelo_filter_genes_scanpy_cpp(
    spliced = spliced,
    unspliced = unspliced,
    min_counts = as.integer(min_counts),
    min_counts_u = as.integer(min_counts_u)
  ) > 0L
  if (sum(keep) < 2L) {
    log_message(
      "Too few genes pass filtering for {.fn RunSCVELO} cpp backend",
      message_type = "error"
    )
  }
  features_out <- features[keep]
  normed <- scvelo_normalize_scanpy_cpp(
    spliced = spliced[keep, , drop = FALSE],
    unspliced = unspliced[keep, , drop = FALSE],
    initial_spliced_totals = initial_spliced_totals,
    initial_unspliced_totals = initial_unspliced_totals
  )
  spliced_n <- normed[["spliced_norm"]]
  unspliced_n <- normed[["unspliced_norm"]]
  linear_embedding_all <- as.matrix(
    srt@reductions[[linear_reduction]]@cell.embeddings[cells, , drop = FALSE]
  )
  dims_use <- seq_len(min(as.integer(n_pcs), ncol(linear_embedding_all)))
  linear_embedding <- linear_embedding_all[, dims_use, drop = FALSE]
  storage.mode(linear_embedding) <- "double"
  knn_nonself <- max(1L, min(as.integer(n_neighbors) - 1L, nrow(linear_embedding) - 1L))
  knn <- scvelo_knn_scanpy_cpp(linear_embedding, knn_nonself, TRUE)
  knn_k <- ncol(knn[["idx"]])
  moments <- scvelo_moments_connectivities_cpp(
    spliced = spliced_n,
    unspliced = unspliced_n,
    knn_idx = knn[["idx"]]
  )
  Ms <- moments[["Ms"]]
  Mu <- moments[["Mu"]]

  # Extract UMAP embedding for velocity projection (visualization only)
  nonlinear_embedding <- as.matrix(
    srt@reductions[[nonlinear_reduction]]@cell.embeddings[cells, , drop = FALSE]
  )
  storage.mode(nonlinear_embedding) <- "double"

  # ── Per-mode velocity computation ──
  srt@tools[["SCVELO"]] <- list(
    backend = "cpp",
    group.by = group.by,
    features = features_out,
    parameters = list(
      spliced_assay = spliced_assay,
      unspliced_assay = unspliced_assay,
      layer_y = layer_y,
      linear_reduction = linear_reduction,
      nonlinear_reduction = nonlinear_reduction,
      use_rep = linear_reduction,
      n_pcs = length(dims_use),
      n_neighbors = as.integer(n_neighbors),
      knn_k = ncol(knn[["idx"]]),
      filter_genes = filter_genes,
      normalize_per_cell = normalize_per_cell,
      log_transform = log_transform,
      n_genes_filtered = as.integer(length(features_out))
    )
  )

  for (m in mode_use) {
    log_message(
      "Running {.pkg scVelo} {.val {m}} mode with {.arg backend = 'cpp'} ({.val {nrow(spliced_n)}} features)",
      verbose = verbose
    )

    # Velocity estimation
    if (identical(m, "stochastic")) {
      velocity <- scvelo_stochastic_cpp(
        Ms = Ms, Mu = Mu,
        Mss = moments[["Mss"]],
        Mus = moments[["Mus"]],
        knn_idx = knn[["idx"]],
        embedding = nonlinear_embedding
      )
    } else if (identical(m, "deterministic")) {
      velocity <- scvelo_deterministic_cpp(
        Ms = Ms, Mu = Mu,
        knn_idx = knn[["idx"]],
        embedding = nonlinear_embedding,
        fit_offset = FALSE,
        perc = 95.0
      )
    } else if (identical(m, "dynamical")) {
      # First fit the dynamical model per gene
      n_genes <- nrow(Ms)
      dyn_genes <- if (n_genes > 200) sample.int(n_genes, min(n_genes, 200)) else seq_len(n_genes)
      fitting_by <- match.arg(fitting_by)
      fitting_by <- switch(fitting_by,
        stochastic = "em",
        deterministic = "nm",
        fitting_by
      )
      if (identical(fitting_by, "em")) {
        dyn_fit <- scvelo_dynamical_em_cpp(
          Ms = Ms, Mu = Mu,
          use_genes = as.integer(dyn_genes),
          max_iter_em = 10L,
          conv_tol = 1e-6
        )
      } else {
        dyn_fit <- scvelo_dynamical_nm_cpp(
          Ms = Ms, Mu = Mu,
          use_genes = as.integer(dyn_genes),
          max_iter = 20L
        )
      }
      # Compute velocity from fitted dynamical parameters
      velocity <- scvelo_dynamical_velocity_cpp(
        Ms = Ms, Mu = Mu,
        alpha = dyn_fit[["alpha"]],
        beta  = dyn_fit[["beta"]],
        gamma = dyn_fit[["gamma"]],
        t_    = dyn_fit[["t_"]],
        knn_idx = knn[["idx"]],
        embedding = nonlinear_embedding
      )
    } else {
      log_message("Unknown mode {.val {m}}", message_type = "error")
    }

    velocity_embedding <- velocity[["velocity_embedding"]]
    rownames(velocity_embedding) <- cells
    colnames(velocity_embedding) <- colnames(nonlinear_embedding)
    velocity_reduction <- paste0(m, "_", nonlinear_reduction)
    srt[[velocity_reduction]] <- SeuratObject::CreateDimReducObject(
      embeddings = velocity_embedding,
      assay = spliced_assay,
      key = paste0(gsub("_", "", velocity_reduction), "_")
    )
    conf_key <- paste0(m, "_confidence")
    len_key  <- paste0(m, "_length")
    srt[[conf_key]] <- as.numeric(velocity[["confidence"]])
    srt[[len_key]]  <- as.numeric(velocity[["velocity_length"]])

    srt@tools[["SCVELO"]][[m]] <- list(
      velocity_reduction = velocity_reduction,
      confidence_key = conf_key,
      length_key = len_key
    )
    if (identical(m, "dynamical")) {
      srt@tools[["SCVELO"]][[m]]$alpha <- dyn_fit[["alpha"]]
      srt@tools[["SCVELO"]][[m]]$beta  <- dyn_fit[["beta"]]
      srt@tools[["SCVELO"]][[m]]$gamma <- dyn_fit[["gamma"]]
      srt@tools[["SCVELO"]][[m]]$t_    <- dyn_fit[["t_"]]
      srt@tools[["SCVELO"]][[m]]$loss  <- dyn_fit[["loss"]]
      srt@tools[["SCVELO"]][[m]]$n_fitted <- dyn_fit[["n_fitted"]]
    } else {
      srt@tools[["SCVELO"]][[m]]$gamma <- velocity[["gamma"]]
    }
    if ("r2" %in% names(velocity)) {
      srt@tools[["SCVELO"]][[m]]$r2 <- velocity[["r2"]]
      srt@tools[["SCVELO"]][[m]]$velocity_genes <- velocity[["velocity_genes"]]
    }
    if ("residual" %in% names(velocity)) {
      srt@tools[["SCVELO"]][[m]]$residual <- velocity[["residual"]]
    }
    vg_residual <- if (!is.null(velocity[["residual"]])) velocity[["residual"]] else velocity[["velocity"]]
    graph_gene_idx <- rep(TRUE, nrow(Ms))
    if ("velocity_genes" %in% names(velocity)) {
      graph_gene_idx <- as.logical(velocity[["velocity_genes"]])
      graph_gene_idx[is.na(graph_gene_idx)] <- FALSE
      if (sum(graph_gene_idx) < 2L) {
        graph_gene_idx <- rep(TRUE, nrow(Ms))
      }
    }
    srt@tools[["SCVELO"]][[m]]$n_velocity_graph_genes <- sum(graph_gene_idx)
    vc_main <- scvelo_velocity_confidence_cpp(
      Ms = Ms[graph_gene_idx, , drop = FALSE],
      residual = vg_residual[graph_gene_idx, , drop = FALSE],
      knn_idx = knn[["idx"]]
    )
    srt[[conf_key]] <- as.numeric(vc_main[["confidence"]])
    srt[[len_key]] <- as.numeric(vc_main[["velocity_length"]])
    srt@tools[["SCVELO"]][[m]]$confidence_detail <- vc_main[["confidence"]]
    srt@tools[["SCVELO"]][[m]]$confidence_diff <- vc_main[["confidence_diff"]]
    # Velocity graph (cosine similarity on gene space, sparse format)
    vg <- scvelo_velocity_graph_cpp(
      Ms = Ms[graph_gene_idx, , drop = FALSE],
      Mu = Mu[graph_gene_idx, , drop = FALSE],
      residual = vg_residual[graph_gene_idx, , drop = FALSE],
      knn_idx = knn[["idx"]],
      sqrt_transform = identical(m, "stochastic"),
      n_recurse_neighbors = 1L
    )
    srt@tools[["SCVELO"]][[m]]$velocity_graph <- list(
      rows = vg[["velocity_graph_rows"]],
      cols = vg[["velocity_graph_cols"]],
      vals = vg[["velocity_graph_vals"]],
      neg_rows = vg[["velocity_graph_neg_rows"]],
      neg_cols = vg[["velocity_graph_neg_cols"]],
      neg_vals = vg[["velocity_graph_neg_vals"]]
    )

    # Terminal states
    if (isTRUE(compute_terminal_states)) {
      ts <- scvelo_terminal_states_graph_cpp(
        graph_rows = vg[["velocity_graph_rows"]],
        graph_cols = vg[["velocity_graph_cols"]],
        graph_vals = vg[["velocity_graph_vals"]],
        graph_neg_rows = vg[["velocity_graph_neg_rows"]],
        graph_neg_cols = vg[["velocity_graph_neg_cols"]],
        graph_neg_vals = vg[["velocity_graph_neg_vals"]],
        knn_idx = knn[["idx"]],
        eps = 1e-3
      )
      rc_key <- paste0(m, "_root_cells")
      ep_key <- paste0(m, "_end_points")
      srt[[rc_key]] <- as.numeric(ts[["root_cells"]])
      srt[[ep_key]] <- as.numeric(ts[["end_points"]])
      srt@tools[["SCVELO"]][[m]]$root_cells <- ts[["root_cells"]]
      srt@tools[["SCVELO"]][[m]]$end_points <- ts[["end_points"]]
    }

    # Velocity pseudotime
    if (isTRUE(compute_pseudotime) && isTRUE(compute_terminal_states)) {
      vpt_result <- scvelo_pseudotime_graph_cpp(
        graph_rows = vg[["velocity_graph_rows"]],
        graph_cols = vg[["velocity_graph_cols"]],
        graph_vals = vg[["velocity_graph_vals"]],
        graph_neg_rows = vg[["velocity_graph_neg_rows"]],
        graph_neg_cols = vg[["velocity_graph_neg_cols"]],
        graph_neg_vals = vg[["velocity_graph_neg_vals"]],
        knn_idx = knn[["idx"]],
        root_cells = ts[["root_cells"]],
        end_points = ts[["end_points"]],
        n_dcs = 10L
      )
      vpt <- as.numeric(vpt_result[["pseudotime"]])
      pt_key <- paste0(m, "_pseudotime")
      srt[[pt_key]] <- as.numeric(vpt)
      srt@tools[["SCVELO"]][[m]]$pseudotime <- vpt
      srt@tools[["SCVELO"]][[m]]$pseudotime_raw <- as.numeric(vpt_result[["pseudotime_root"]])
      srt@tools[["SCVELO"]][[m]]$pseudotime_end_inverse <- as.numeric(vpt_result[["pseudotime_end_inverse"]])
      srt@tools[["SCVELO"]][[m]]$root_cell <- vpt_result[["root_cell"]]
      srt@tools[["SCVELO"]][[m]]$end_cell <- vpt_result[["end_cell"]]
      srt@tools[["SCVELO"]][[m]]$diffusion_components <- vpt_result[["diffusion_components"]]
    }

    # Velocity confidence metrics
    if (isTRUE(compute_velocity_confidence)) {
      vc <- scvelo_velocity_confidence_cpp(
        Ms = Ms[graph_gene_idx, , drop = FALSE],
        residual = vg_residual[graph_gene_idx, , drop = FALSE],
        knn_idx = knn[["idx"]]
      )
      conf_detail_key <- paste0(m, "_confidence_detail")
      srt@tools[["SCVELO"]][[m]]$confidence_detail <- vc[["confidence"]]
      srt@tools[["SCVELO"]][[m]]$confidence_diff <- vc[["confidence_diff"]]
    }

    log_message(
      "{.pkg scVelo} {.val {m}} mode completed",
      message_type = "success",
      verbose = verbose
    )
  }

  srt@tools[["SCVELO"]]$mode <- mode_use
  log_message(
    "{.pkg scVelo} cpp backend completed",
    message_type = "success",
    verbose = verbose
  )
  srt
}
