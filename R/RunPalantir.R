#' @title Run Palantir analysis
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams RunCellRank
#' @param backend Backend used to compute Palantir. `"python"` keeps the
#' original Palantir workflow and remains the default. `"cpp"` uses the
#' C++ implementation and stores results in `srt@misc$palantir`. Default is `"cpp"`.
#' @param dm_n_components The number of diffusion components to calculate.
#' @param dm_alpha Normalization parameter for the diffusion operator.
#' @param dm_n_eigs Number of eigen vectors to use.
#' @param early_group Name of the group to start Palantir analysis from.
#' @param early_cell Name of the cell to start Palantir analysis from.
#' @param terminal_groups Character vector specifying terminal groups for Palantir analysis.
#' @param terminal_cells Character vector specifying terminal cells for Palantir analysis.
#' @param num_waypoints Number of waypoints to be included.
#' @param scale_components Should the cell fate probabilities be scaled for each component independently?
#' @param use_early_cell_as_start Should the starting cell for each terminal group be set as early_cell?
#' @param adjust_early_cell Whether to adjust the early cell to the cell with the minimum pseudotime value.
#' @param adjust_terminal_cells Whether to adjust the terminal cells to the cells with the maximum pseudotime value for each terminal group.
#' @param max_iterations Maximum number of iterations for pseudotime convergence.
#' @param point_size The point size for plotting.
#' @param plot_format Format for saved plots: `"pdf"`, `"png"`, or `"svg"`. Default is `"pdf"`.
#' @param plot_prefix Prefix for saved plot filenames. Default is `"palantir"`.
#' @param dirpath The directory to save the plots. Default is `"./"`.
#'
#' @export
#'
#' @seealso [PalantirTrajectoryPlot]
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunPalantir(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   linear_reduction = "PCA",
#'   nonlinear_reduction = "UMAP",
#'   early_group = "Ductal",
#'   terminal_groups = c("Alpha", "Beta", "Delta", "Epsilon")
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   c("palantir_pseudotime", "palantir_diff_potential")
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   grep(
#'     "TerminalState_.*_diff_potential$",
#'     colnames(pancreas_sub@meta.data),
#'     value = TRUE
#'   )
#' )
#'
#' PalantirTrajectoryPlot(
#'   pancreas_sub,
#'   reduction = "UMAP",
#'   pseudotime_interval = c(0, 0.9)
#' )
#'
#' PalantirTrajectoryPlot(
#'   pancreas_sub,
#'   reduction = "UMAP",
#'   cell_color = "branch_selection",
#'   pseudotime_interval = c(0, 0.9)
#' )
RunPalantir <- function(
  srt = NULL,
  assay_x = "RNA",
  layer_x = "counts",
  assay_y = c("spliced", "unspliced"),
  layer_y = "counts",
  adata = NULL,
  group.by = NULL,
  linear_reduction = NULL,
  nonlinear_reduction = NULL,
  basis = NULL,
  n_pcs = 30,
  n_neighbors = 30,
  dm_n_components = 10,
  dm_alpha = 0,
  dm_n_eigs = NULL,
  early_group = NULL,
  early_cell = NULL,
  terminal_cells = NULL,
  terminal_groups = NULL,
  num_waypoints = 1200,
  scale_components = TRUE,
  use_early_cell_as_start = TRUE,
  adjust_early_cell = FALSE,
  adjust_terminal_cells = FALSE,
  max_iterations = 25,
  cores = 1,
  point_size = 20,
  palette = "Chinese",
  palcolor = NULL,
  legend.position = "on data",
  show_plot = FALSE,
  save_plot = FALSE,
  plot_format = c("pdf", "png", "svg"),
  plot_dpi = 300,
  plot_prefix = "palantir",
  dirpath = "./",
  backend = c("cpp", "python"),
  return_seurat = !is.null(srt),
  verbose = TRUE
) {
  backend <- match.arg(backend)
  plot_format <- match.arg(plot_format)

  if (identical(backend, "cpp")) {
    if (is.null(srt)) {
      log_message(
        "{.arg backend = 'cpp'} requires {.arg srt}; use {.arg backend = 'python'} for AnnData input",
        message_type = "error"
      )
    }
    if (!isTRUE(return_seurat)) {
      log_message(
        "{.arg backend = 'cpp'} returns a {.cls Seurat} object only",
        message_type = "error"
      )
    }
    if (is.null(linear_reduction) && is.null(nonlinear_reduction)) {
      log_message(
        "{.arg linear_reduction} or {.arg nonlinear_reduction} must be provided at least one.",
        message_type = "error"
      )
    }
    if (is.null(early_cell) && is.null(early_group)) {
      log_message(
        "{.arg early_cell} or {.arg early_group} must be provided.",
        message_type = "error"
      )
    }
    unsupported_cpp <- character(0)
    if (isTRUE(adjust_early_cell)) {
      unsupported_cpp <- c(unsupported_cpp, "adjust_early_cell")
    }
    if (isTRUE(adjust_terminal_cells)) {
      unsupported_cpp <- c(unsupported_cpp, "adjust_terminal_cells")
    }
    if (isTRUE(show_plot)) {
      unsupported_cpp <- c(unsupported_cpp, "show_plot")
    }
    if (isTRUE(save_plot)) {
      unsupported_cpp <- c(unsupported_cpp, "save_plot")
    }
    if (length(unsupported_cpp) > 0L) {
      unsupported_text <- paste0(unsupported_cpp, collapse = ", ")
      log_message(
        "{.arg backend = 'cpp'} currently does not support {.arg {unsupported_text}}; use {.arg backend = 'python'} for those features",
        message_type = "warning",
        verbose = verbose
      )
    }
    return(run_palantir_cpp(
      srt = srt,
      group.by = group.by,
      linear_reduction = linear_reduction,
      nonlinear_reduction = nonlinear_reduction,
      n_pcs = n_pcs,
      n_neighbors = n_neighbors,
      dm_n_components = dm_n_components,
      dm_alpha = dm_alpha,
      early_group = early_group,
      early_cell = early_cell,
      terminal_cells = terminal_cells,
      terminal_groups = terminal_groups,
      num_waypoints = num_waypoints,
      scale_components = scale_components,
      use_early_cell_as_start = use_early_cell_as_start,
      max_iterations = max_iterations,
      cores = cores,
      verbose = verbose
    ))
  }

  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1",
    KMP_WARNINGS = "0",
    KMP_DUPLICATE_LIB_OK = "TRUE"
  )
  PrepareEnv(modules = "palantir")
  check_python("palantir", verbose = verbose)
  if (all(is.null(srt), is.null(adata))) {
    log_message(
      "{.arg srt} or {.arg adata} must be provided.",
      message_type = "error"
    )
  }
  if (
    is.null(group.by) && any(!is.null(early_group), !is.null(terminal_groups))
  ) {
    log_message(
      "{.arg group.by} must be provided when {.arg early_group} or {.arg terminal_groups} provided.",
      message_type = "error"
    )
  }
  if (is.null(linear_reduction) && is.null(nonlinear_reduction)) {
    log_message(
      "{.arg linear_reduction} or {.arg nonlinear_reduction} must be provided at least one.",
      message_type = "error"
    )
  }
  if (is.null(early_cell) && is.null(early_group)) {
    log_message(
      "{.arg early_cell} or {.arg early_group} must be provided.",
      message_type = "error"
    )
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

  args[["legend_loc"]] <- legend.position
  args[["n_jobs"]] <- as.integer(cores)
  args[["save"]] <- save_plot
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
    "save_plot",
    "plot_dpi",
    "plot_prefix",
    "legend.position",
    "cores",
    "backend"
  )
  args <- args[!names(args) %in% params]

  if (!is.null(srt)) {
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
      nonlinear_reduction <- DefaultReduction(
        srt,
        pattern = nonlinear_reduction
      )
    }
    if (!nonlinear_reduction %in% names(srt@reductions)) {
      log_message(
        "{.val {nonlinear_reduction}} is not in the srt reduction names",
        message_type = "error"
      )
    }

    args[["adata"]] <- srt_to_adata(
      srt = srt,
      assay_x = assay_x,
      layer_x = layer_x,
      assay_y = assay_y,
      layer_y = layer_y
    )

    if (!is.null(linear_reduction)) {
      args[["linear_reduction"]] <- linear_reduction
    }
    if (!is.null(nonlinear_reduction)) {
      args[["nonlinear_reduction"]] <- nonlinear_reduction
    }
    if (is.null(basis)) {
      if (!is.null(nonlinear_reduction)) {
        args[["basis"]] <- nonlinear_reduction
      } else if (!is.null(linear_reduction)) {
        args[["basis"]] <- linear_reduction
      }
    } else {
      args[["basis"]] <- basis
    }
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
  adata <- do.call(functions$Palantir, args)

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      srt_out1 <- srt_append(srt_raw = srt, srt_append = srt_out)
      srt_out2 <- srt_append(
        srt_raw = srt_out1,
        srt_append = srt_out,
        pattern = "(palantir)|(dm_kernel)|(_diff_potential)",
        overwrite = TRUE,
        verbose = FALSE
      )
      return(srt_out2)
    }
  } else {
    return(adata)
  }
}

run_palantir_cpp <- function(
  srt,
  group.by = NULL,
  linear_reduction = NULL,
  nonlinear_reduction = NULL,
  n_pcs = 30,
  n_neighbors = 30,
  dm_n_components = 10,
  dm_alpha = 0,
  early_group = NULL,
  early_cell = NULL,
  terminal_cells = NULL,
  terminal_groups = NULL,
  num_waypoints = 1200,
  scale_components = TRUE,
  use_early_cell_as_start = TRUE,
  max_iterations = 25,
  cores = 1,
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a Seurat object", message_type = "error")
  }
  if (!is.null(linear_reduction)) {
    linear_reduction <- DefaultReduction(srt, pattern = linear_reduction)
  }
  if (!is.null(nonlinear_reduction)) {
    nonlinear_reduction <- DefaultReduction(srt, pattern = nonlinear_reduction)
  }
  if (
    !is.null(linear_reduction) && linear_reduction %in% names(srt@reductions)
  ) {
    analysis_reduction <- linear_reduction
  } else if (
    !is.null(nonlinear_reduction) &&
      nonlinear_reduction %in% names(srt@reductions)
  ) {
    analysis_reduction <- nonlinear_reduction
  } else {
    log_message("No valid reduction found", message_type = "error")
  }
  basis_reduction <- if (
    !is.null(nonlinear_reduction) &&
      nonlinear_reduction %in% names(srt@reductions)
  ) {
    nonlinear_reduction
  } else {
    analysis_reduction
  }

  embedding_all <- as.matrix(
    srt@reductions[[analysis_reduction]]@cell.embeddings
  )
  dims_use <- seq_len(min(as.integer(n_pcs), ncol(embedding_all)))
  embedding <- embedding_all[, dims_use, drop = FALSE]
  storage.mode(embedding) <- "double"
  n_cells <- nrow(embedding)
  cell_names <- rownames(embedding)

  basis_embedding <- as.matrix(
    srt@reductions[[basis_reduction]]@cell.embeddings
  )
  basis_dims <- seq_len(min(2L, ncol(basis_embedding)))
  basis_embedding <- basis_embedding[, basis_dims, drop = FALSE]
  storage.mode(basis_embedding) <- "double"

  group_medoid_cell <- function(group_value) {
    if (is.null(group.by) || !group.by %in% colnames(srt@meta.data)) {
      log_message(
        "{.arg group.by} must be provided when group labels are used.",
        message_type = "error"
      )
    }
    idx <- which(srt@meta.data[[group.by]] == group_value)
    if (length(idx) == 0L) {
      log_message(
        "Group {.val {group_value}} is not present in {.arg group.by}",
        message_type = "error"
      )
    }
    coords <- basis_embedding[idx, , drop = FALSE]
    center <- apply(coords, 2, stats::median)
    idx[which.min(rowSums((sweep(coords, 2, center, "-"))^2))]
  }

  if (is.null(early_cell)) {
    if (!is.null(early_group)) {
      early_cell <- group_medoid_cell(early_group)
    } else {
      early_cell <- 1L
    }
  }
  if (is.character(early_cell)) {
    early_cell <- match(early_cell, cell_names)
  }
  early_cell <- as.integer(early_cell)
  if (is.na(early_cell) || early_cell < 1L || early_cell > n_cells) {
    log_message("Invalid {.arg early_cell}", message_type = "error")
  }

  terminal_state_cells <- integer(0)
  if (!is.null(terminal_groups)) {
    terminal_state_cells <- c(
      terminal_state_cells,
      vapply(terminal_groups, group_medoid_cell, integer(1))
    )
  }
  if (!is.null(terminal_cells)) {
    terminal_from_cells <- if (is.character(terminal_cells)) {
      match(terminal_cells, cell_names)
    } else {
      as.integer(terminal_cells)
    }
    if (
      any(
        is.na(terminal_from_cells) |
          terminal_from_cells < 1L |
          terminal_from_cells > n_cells
      )
    ) {
      log_message("Invalid {.arg terminal_cells}", message_type = "error")
    }
    terminal_state_cells <- c(terminal_state_cells, terminal_from_cells)
  }
  terminal_state_cells <- unique(setdiff(
    as.integer(terminal_state_cells),
    early_cell
  ))

  knn_k <- max(1L, min(as.integer(n_neighbors), n_cells - 1L))
  kernel_knn_k <- max(1L, min(as.integer(n_neighbors) - 1L, n_cells - 1L))
  log_message("Computing {.pkg Palantir} cpp KNN graph...", verbose = verbose)
  knn_res <- run_cpp_knn(
    reference = embedding,
    query = embedding,
    k = kernel_knn_k,
    metric = "euclidean",
    exclude_self = TRUE,
    n_threads = as.integer(cores)
  )

  kernel <- palantir_compute_kernel_cpp(
    data = embedding,
    knn_idx = knn_res[["idx"]],
    knn_dist = knn_res[["dist"]],
    knn = knn_k,
    alpha = dm_alpha
  )
  norm <- palantir_normalize_kernel_cpp(
    kernel_i = as.integer(kernel[["i"]]),
    kernel_j = as.integer(kernel[["j"]]),
    kernel_x = as.numeric(kernel[["x"]]),
    n = as.integer(kernel[["n"]])
  )

  n_eigs <- min(as.integer(dm_n_components), n_cells - 2L)
  T_dm <- Matrix::sparseMatrix(
    i = as.integer(norm[["T_i"]]),
    j = as.integer(norm[["T_j"]]),
    x = as.numeric(norm[["T_x"]]),
    dims = c(kernel[["n"]], kernel[["n"]])
  )
  eig <- RSpectra::eigs(
    T_dm,
    k = n_eigs,
    which = "LM",
    opts = list(
      tol = 1e-4,
      maxitr = 1000,
      initvec = palantir_numpy_random_sample_cpp(nrow(T_dm), seed = 0L)
    )
  )
  eig_order <- order(Re(eig$values), decreasing = TRUE)
  eigenvalues <- Re(eig$values)[eig_order]
  eigenvectors <- Re(eig$vectors)[, eig_order, drop = FALSE]
  eigenvectors <- apply(eigenvectors, 2, function(x) {
    norm_x <- sqrt(sum(x^2))
    if (norm_x > 0) x / norm_x else x
  })
  ms_data <- palantir_multiscale_space_cpp(eigenvectors, eigenvalues)

  dm_boundaries <- unique(c(
    apply(ms_data, 2, which.max),
    apply(ms_data, 2, which.min)
  ))
  dists_to_ec <- as.matrix(stats::dist(rbind(
    ms_data[early_cell, , drop = FALSE],
    ms_data[dm_boundaries, , drop = FALSE]
  )))
  start_cell <- dm_boundaries[which.min(dists_to_ec[1, -1])]
  if (isTRUE(use_early_cell_as_start)) {
    start_cell <- early_cell
  }

  if (isTRUE(scale_components)) {
    ms_min <- apply(ms_data, 2, min)
    ms_max <- apply(ms_data, 2, max)
    ms_range <- ms_max - ms_min
    ms_range[ms_range < 1e-10] <- 1
    ms_data_scaled <- sweep(sweep(ms_data, 2, ms_min, "-"), 2, ms_range, "/")
  } else {
    ms_data_scaled <- ms_data
  }

  n_wp_target <- min(as.integer(num_waypoints), n_cells)
  wp <- palantir_maxmin_waypoints_cpp(
    ms_data = ms_data_scaled,
    num_waypoints = as.integer(n_wp_target),
    seed = 20L
  )
  wp <- unique(as.integer(c(
    start_cell,
    terminal_state_cells,
    setdiff(c(wp, dm_boundaries), c(start_cell, terminal_state_cells))
  )))
  n_wp <- length(wp)

  pt_res <- palantir_pseudotime_cpp(
    ms_data = ms_data_scaled,
    start_cell = as.integer(start_cell - 1L),
    waypoints = wp,
    knn = as.integer(knn_k),
    max_iterations = as.integer(max_iterations),
    n_jobs = as.integer(cores)
  )
  pseudotime <- as.numeric(pt_res[["pseudotime"]])
  pseudotime_raw <- pseudotime
  names(pseudotime) <- rownames(embedding)
  names(pseudotime_raw) <- rownames(embedding)

  wp_ms <- ms_data_scaled[wp, , drop = FALSE]
  wp_pseudotime <- pseudotime[wp]
  mc <- palantir_markov_chain_cpp(
    wp_data = wp_ms,
    knn = as.integer(min(knn_k, n_wp - 1L)),
    pseudotime = wp_pseudotime
  )

  T_mat <- Matrix::sparseMatrix(
    i = as.integer(mc[["T_i"]]),
    j = as.integer(mc[["T_j"]]),
    x = as.numeric(mc[["T_x"]]),
    dims = c(mc[["n"]], mc[["n"]])
  )
  T_t <- Matrix::t(T_mat)
  n_T <- nrow(T_t)
  if (n_T <= 2L) {
    eig_T <- eigen(as.matrix(T_t))
    lead_idx <- which.max(Re(eig_T$values))
    ranks <- abs(Re(eig_T$vectors[, lead_idx]))
  } else {
    eig_T <- RSpectra::eigs(T_t, k = 1, which = "LM")
    ranks <- abs(Re(eig_T$vectors[, 1]))
  }
  names(ranks) <- wp

  if (length(terminal_state_cells) > 0L) {
    terminal_states <- intersect(terminal_state_cells, wp)
  } else {
    med_r <- stats::median(ranks)
    mad_r <- stats::median(abs(ranks - med_r))
    cutoff <- med_r + 3.09 * mad_r
    high_rank <- wp[ranks > cutoff]
    if (length(high_rank) == 0L) {
      high_rank <- wp[order(ranks, decreasing = TRUE)][seq_len(min(
        3L,
        length(wp)
      ))]
    }
    wp_boundaries <- intersect(wp, dm_boundaries)
    if (length(wp_boundaries) > 0L) {
      high_rank_idx <- match(high_rank, wp)
      wp_boundary_idx <- match(wp_boundaries, wp)
      dist_to_boundary <- as.matrix(stats::dist(rbind(
        wp_ms[high_rank_idx, , drop = FALSE],
        wp_ms[wp_boundary_idx, , drop = FALSE]
      )))
      n_hr <- length(high_rank)
      n_bnd <- length(wp_boundaries)
      if (n_hr > 0L && n_bnd > 0L) {
        dist_part <- dist_to_boundary[
          seq_len(n_hr),
          (n_hr + 1L):(n_hr + n_bnd),
          drop = FALSE
        ]
        nearest <- apply(dist_part, 1, which.min)
        terminal_states <- unique(wp_boundaries[nearest])
      } else {
        terminal_states <- high_rank
      }
    } else {
      terminal_states <- high_rank
    }
  }
  terminal_states <- setdiff(terminal_states, start_cell)
  n_terminal <- length(terminal_states)

  if (n_terminal > 0L && n_terminal < n_wp) {
    ts_indices <- match(terminal_states, wp)
    for (ts in ts_indices) {
      T_mat[ts, ] <- 0
      T_mat[ts, ts] <- 1
    }
    trans_set <- setdiff(seq_len(n_wp), ts_indices)
    n_trans <- length(trans_set)
    if (n_trans > 0L) {
      Q <- T_mat[trans_set, trans_set, drop = FALSE]
      R <- T_mat[trans_set, ts_indices, drop = FALSE]
      I_Q <- Matrix::Diagonal(n_trans) - Q
      B <- tryCatch(
        as.matrix(Matrix::solve(I_Q, R)),
        error = function(e) {
          I_dense <- as.matrix(I_Q)
          R_dense <- as.matrix(R)
          tryCatch(
            qr.solve(I_dense, R_dense),
            error = function(e2) {
              sv <- svd(I_dense)
              keep <- sv$d > max(dim(I_dense)) * max(sv$d) * .Machine$double.eps
              if (!any(keep)) {
                matrix(0, nrow(I_dense), ncol(R_dense))
              } else {
                sv$v[, keep, drop = FALSE] %*%
                  (diag(1 / sv$d[keep], nrow = sum(keep)) %*%
                    (t(sv$u[, keep, drop = FALSE]) %*% R_dense))
              }
            }
          )
        }
      )
      B[B < 0] <- 0
      branch_probs <- matrix(0, n_wp, n_terminal)
      for (i in seq_len(n_terminal)) {
        branch_probs[ts_indices[i], i] <- 1
      }
      for (i in seq_len(n_trans)) {
        branch_probs[trans_set[i], ] <- B[i, ]
      }
    } else {
      branch_probs <- matrix(0, n_wp, n_terminal)
      for (i in seq_len(n_terminal)) {
        branch_probs[ts_indices[i], i] <- 1
      }
    }
  } else {
    branch_probs <- matrix(0, n_wp, 1)
  }

  rs <- rowSums(branch_probs)
  zero_branch_rows <- which(rs < 1e-10)
  if (length(zero_branch_rows) > 0L && ncol(branch_probs) > 0L) {
    branch_probs[zero_branch_rows, ] <- 1 / ncol(branch_probs)
    rs <- rowSums(branch_probs)
  }
  rs[rs < 1e-10] <- 1
  branch_probs <- branch_probs / rs
  rownames(branch_probs) <- cell_names[wp]
  colnames(branch_probs) <- if (length(terminal_groups) == ncol(branch_probs)) {
    paste0("TerminalState_", terminal_groups)
  } else if (length(terminal_cells) == ncol(branch_probs)) {
    paste0("TerminalState_", terminal_cells)
  } else {
    paste0("TerminalState_", seq_len(ncol(branch_probs)))
  }

  W_mat <- as.matrix(pt_res[["W"]])
  if (nrow(W_mat) == nrow(branch_probs)) {
    branch_probs_all <- t(W_mat) %*% branch_probs
  } else {
    log_message(
      "Skipping {.pkg Palantir} branch probability projection because waypoint dimensions do not match",
      message_type = "warning",
      verbose = verbose
    )
    branch_probs_all <- matrix(0, n_cells, ncol(branch_probs))
  }
  branch_probs_all_rs <- rowSums(branch_probs_all)
  zero_projected_rows <- which(branch_probs_all_rs < 1e-10)
  if (length(zero_projected_rows) > 0L && ncol(branch_probs_all) > 0L) {
    branch_probs_all[zero_projected_rows, ] <- 1 / ncol(branch_probs_all)
    branch_probs_all_rs <- rowSums(branch_probs_all)
  }
  branch_probs_all_rs[branch_probs_all_rs < 1e-10] <- 1
  branch_probs_all <- branch_probs_all / branch_probs_all_rs
  rownames(branch_probs_all) <- cell_names
  colnames(branch_probs_all) <- colnames(branch_probs)

  ent <- apply(branch_probs_all, 1, function(p) {
    p <- p[p > 0]
    if (length(p) < 2L) {
      return(0)
    }
    -sum(p * log(p))
  })
  srt[["palantir_pseudotime"]] <- pseudotime[colnames(srt)]
  srt[["palantir_diff_potential"]] <- ent[colnames(srt)]
  branch_meta <- as.data.frame(branch_probs_all[colnames(srt), , drop = FALSE])
  branch_meta_names <- make.unique(paste0(
    colnames(branch_meta),
    "_diff_potential"
  ))
  colnames(branch_meta) <- branch_meta_names
  srt <- Seurat::AddMetaData(srt, metadata = branch_meta)

  srt@misc[["palantir"]] <- list(
    pseudotime = pseudotime,
    pseudotime_raw = pseudotime_raw,
    entropy = ent,
    branch_probs = branch_probs_all,
    waypoints = wp,
    terminal_states = terminal_states,
    start_cell = start_cell,
    diffusion_components = ms_data_scaled,
    backend = "cpp",
    parameters = list(
      reduction = analysis_reduction,
      basis_reduction = basis_reduction,
      n_pcs = length(dims_use),
      n_neighbors = as.integer(n_neighbors),
      knn_k = knn_k,
      kernel_knn_k = kernel_knn_k,
      dm_n_components = as.integer(dm_n_components),
      dm_alpha = dm_alpha,
      num_waypoints = n_wp,
      scale_components = scale_components,
      use_early_cell_as_start = use_early_cell_as_start,
      max_iterations = max_iterations
    )
  )

  log_message(
    "{.pkg Palantir} cpp backend completed",
    message_type = "success",
    verbose = verbose
  )
  srt
}
