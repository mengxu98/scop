#' @title Run PAGA analysis
#'
#' @description
#' PAGA is a graph-based method used to infer cellular trajectories.
#' This function runs the PAGA analysis on a Seurat object.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams RunCellRank
#' @param backend Backend used to compute PAGA. `"python"` keeps the original
#' scanpy workflow and remains the default. `"cpp"` uses the native C++
#' implementation for the standard connectivity graph and tree, plus an
#' approximate R igraph layout stored in `paga$pos`.
#' @param use_rna_velocity Whether to use RNA velocity for PAGA analysis.
#' Default is `FALSE`.
#' @param vkey The name of the RNA velocity data to use if `use_rna_velocity` is `TRUE`.
#' Default is `"stochastic"`.
#' @param embedded_with_PAGA Whether to embed data using PAGA layout.
#' Default is `FALSE`.
#' @param paga_layout The layout for plotting PAGA graph.
#' See \href{https://scanpy.readthedocs.io/en/stable/tutorials/plotting/advanced.html#paga}{layout} param in `scanpy.pl.paga` function.
#' @param threshold The threshold for plotting PAGA graph.
#' Edges for weights below this threshold will not be drawn.
#' @param point_size The point size for plotting.
#' @param infer_pseudotime Whether to infer pseudotime.
#' @param root_group The group to use as the root for pseudotime inference.
#' @param root_cell The cell to use as the root for pseudotime inference.
#' @param n_dcs The number of diffusion components to use for pseudotime inference.
#' @param n_branchings Number of branchings to detect.
#' @param min_group_size The minimum size of a group (as a fraction of the total number of cells) to consider it as a potential branching point.
#'
#' @seealso
#' [PAGAPlot], [CellDimPlot], [RunSCVELO]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunPAGA(
#'   pancreas_sub,
#'   assay_x = "RNA",
#'   group.by = "SubCellType",
#'   linear_reduction = "PCA",
#'   nonlinear_reduction = "UMAP",
#'   backend = "cpp"
#' )
#' PAGAPlot(pancreas_sub, reduction = "UMAP")
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   paga = pancreas_sub@misc$paga
#' )
RunPAGA <- function(
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
  n_pcs = 30,
  n_neighbors = 30,
  use_rna_velocity = FALSE,
  vkey = "stochastic",
  embedded_with_PAGA = FALSE,
  paga_layout = "fr",
  threshold = 0.1,
  point_size = 20,
  infer_pseudotime = FALSE,
  root_group = NULL,
  root_cell = NULL,
  n_dcs = 10,
  n_branchings = 0,
  min_group_size = 0.01,
  palette = "Chinese",
  palcolor = NULL,
  legend.position = "on data",
  cores = 1,
  show_plot = FALSE,
  save_plot = FALSE,
  plot_format = c("pdf", "png", "svg"),
  plot_dpi = 300,
  plot_prefix = "paga",
  dirpath = "./paga",
  backend = c("python", "cpp"),
  return_seurat = !is.null(srt),
  verbose = TRUE
) {
  backend <- match.arg(backend)
  plot_format <- match.arg(plot_format)

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
  if (!is.null(srt)) {
    if (is.null(linear_reduction)) {
      linear_reduction <- DefaultReduction(srt)
    } else {
      linear_reduction <- DefaultReduction(srt, pattern = linear_reduction)
    }

    if (is.null(nonlinear_reduction)) {
      nonlinear_reduction <- DefaultReduction(srt)
    } else {
      nonlinear_reduction <- DefaultReduction(
        srt,
        pattern = nonlinear_reduction
      )
    }
  }

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
    unsupported_cpp <- character(0)
    if (isTRUE(embedded_with_PAGA)) {
      unsupported_cpp <- c(unsupported_cpp, "embedded_with_PAGA")
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
        message_type = "warning"
      )
    }
    srt <- run_paga_cpp(
      srt = srt,
      group.by = group.by,
      linear_reduction = linear_reduction,
      nonlinear_reduction = nonlinear_reduction,
      n_pcs = n_pcs,
      n_neighbors = n_neighbors,
      paga_layout = paga_layout,
      threshold = threshold,
      cores = cores,
      infer_pseudotime = infer_pseudotime,
      root_group = root_group,
      root_cell = root_cell,
      n_dcs = n_dcs,
      n_branchings = n_branchings,
      min_group_size = min_group_size,
      use_rna_velocity = use_rna_velocity,
      vkey = vkey,
      verbose = verbose
    )
    return(srt)
  }

  PrepareEnv(modules = "scanpy")

  args <- mget(names(formals()))
  args <- lapply(
    args, function(x) {
      if (is.numeric(x)) {
        y <- ifelse(
          grepl(
            "\\.",
            as.character(x)
          ),
          as.double(x),
          as.integer(x)
        )
      } else {
        y <- x
      }
      y
    }
  )
  call_envir <- parent.frame(1)
  args <- lapply(
    args, function(arg) {
      if (is.symbol(arg)) {
        eval(arg, envir = call_envir)
      } else if (is.call(arg)) {
        eval(arg, envir = call_envir)
      } else {
        arg
      }
    }
  )

  args[["legend_loc"]] <- legend.position
  args[["n_jobs"]] <- cores
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
    args[["adata"]] <- srt_to_adata(
      srt = srt,
      assay_x = assay_x,
      layer_x = layer_x,
      assay_y = assay_y,
      layer_y = layer_y
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
  log_message("Running {.pkg PAGA} analysis...", verbose = verbose)
  adata <- do.call(functions$PAGA, args)
  log_message(
    "{.pkg PAGA} analysis completed",
    message_type = "success",
    verbose = verbose
  )

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      srt_out1 <- srt_append(
        srt_raw = srt,
        srt_append = srt_out
      )
      srt_out2 <- srt_append(
        srt_raw = srt_out1,
        srt_append = srt_out,
        pattern = "(paga)|(distances)|(connectivities)|(draw_graph)",
        overwrite = TRUE,
        verbose = FALSE
      )
      return(srt_out2)
    }
  } else {
    return(adata)
  }
}

run_paga_cpp <- function(
  srt,
  group.by,
  linear_reduction,
  nonlinear_reduction = NULL,
  n_pcs,
  n_neighbors,
  paga_layout,
  threshold,
  cores,
  infer_pseudotime = FALSE,
  root_group = NULL,
  root_cell = NULL,
  n_dcs = 10,
  n_branchings = 0,
  min_group_size = 0.01,
  use_rna_velocity = FALSE,
  vkey = "stochastic",
  knn_override = NULL,
  verbose = TRUE
) {
  if (is.null(linear_reduction) || !linear_reduction %in% names(srt@reductions)) {
    log_message(
      "{.arg linear_reduction} must identify an existing reduction for {.arg backend = 'cpp'}",
      message_type = "error"
    )
  }
  if (!group.by %in% colnames(srt@meta.data)) {
    log_message(
      "{.arg group.by} {.val {group.by}} is not present in {.arg srt}",
      message_type = "error"
    )
  }

  groups <- srt@meta.data[[group.by]]
  if (!is.factor(groups)) {
    groups <- factor(groups)
    srt@meta.data[[group.by]] <- groups
  }
  if (nlevels(groups) < 2L) {
    log_message(
      "{.arg group.by} must contain at least two groups for {.fn RunPAGA}",
      message_type = "error"
    )
  }

  if (!is.null(knn_override)) {
    log_message(
      "Using externally-provided KNN graph for parity comparison",
      verbose = verbose
    )
    knn <- knn_override
    knn_k <- ncol(knn[["idx"]])
    dims_use <- NULL
  } else {
    embedding <- srt@reductions[[linear_reduction]]@cell.embeddings
    if (nrow(embedding) != ncol(srt)) {
      log_message(
        "{.arg linear_reduction} embeddings must contain one row per cell",
        message_type = "error"
      )
    }
    dims_use <- seq_len(min(as.integer(n_pcs), ncol(embedding)))
    embedding <- as.matrix(embedding[, dims_use, drop = FALSE])
    storage.mode(embedding) <- "double"

    knn_k <- max(1L, min(as.integer(n_neighbors) - 1L, nrow(embedding) - 1L))
    log_message(
      "Running {.pkg PAGA} with {.arg backend = 'cpp'} using {.val {knn_k}} neighbors",
      verbose = verbose
    )
    knn <- run_cpp_knn(
      reference = embedding,
      query = embedding,
      k = knn_k,
      metric = "euclidean",
      exclude_self = TRUE,
      n_threads = as.integer(cores)
    )
  }
  paga <- paga_connectivities_cpp(
    knn_idx = knn[["idx"]],
    groups = as.integer(groups),
    n_groups = nlevels(groups)
  )
  group_levels <- levels(groups)
  for (nm in c("connectivities", "connectivities_tree")) {
    dimnames(paga[[nm]]) <- list(NULL, NULL)
  }
  names(paga[["group_sizes"]]) <- group_levels
  pos <- paga_layout_igraph(
    connectivities = paga[["connectivities"]],
    connectivities_tree = paga[["connectivities_tree"]],
    layout = paga_layout,
    threshold = threshold
  )

  srt@misc[["paga"]] <- list(
    connectivities = paga[["connectivities"]],
    connectivities_tree = paga[["connectivities_tree"]],
    groups = group.by,
    pos = pos,
    backend = "cpp",
    parameters = list(
      linear_reduction = linear_reduction,
      n_pcs = length(dims_use),
      n_neighbors = as.integer(n_neighbors),
      knn_k = knn_k,
      paga_layout = paga_layout,
      threshold = threshold,
      layout_engine = "igraph_r_approximate",
      connectivity_engine = "cpp_exact",
      tree_engine = "cpp_exact",
      layout_seed = 0L
    )
  )
  srt@misc[[paste0(group.by, "_sizes")]] <- as.numeric(paga[["group_sizes"]])
  names(srt@misc[[paste0(group.by, "_sizes")]]) <- group_levels

  # Diffusion pseudotime
  if (isTRUE(infer_pseudotime)) {
    root_grp <- if (!is.null(root_group)) {
      match(root_group, group_levels)
    } else if (!is.null(root_cell)) {
      match(as.character(srt@meta.data[[root_cell, group.by]]), group_levels)
    } else {
      1L
    }
    if (anyNA(root_grp) || length(root_grp) == 0L) root_grp <- 1L

    dpt <- paga_diffusion_pseudotime_cpp(
      connectivities = paga[["connectivities"]],
      root_group = as.integer(root_grp),
      n_dcs = as.integer(n_dcs),
      n_branchings = as.integer(n_branchings),
      group_sizes = as.numeric(paga[["group_sizes"]]),
      min_group_size = min_group_size
    )
    srt@misc[["paga"]]$pseudotime <- dpt$pseudotime
    srt@misc[["paga"]]$diffusion_components <- dpt$diffusion_components
    srt@misc[["paga"]]$diffusion_eigenvalues <- dpt$diffusion_eigenvalues
    srt@misc[["paga"]]$root_group <- dpt$root_group
    srt@misc[["paga"]]$parameters$n_dcs <- n_dcs
    srt@misc[["paga"]]$parameters$infer_pseudotime <- TRUE
    srt@misc[["paga"]]$parameters$n_branchings <- n_branchings
    names(dpt$pseudotime) <- group_levels
    log_message(
      "{.pkg PAGA} diffusion pseudotime computed (root: {.val {group_levels[root_grp]}})",
      message_type = "success",
      verbose = verbose
    )
  }

  # Velocity-based PAGA transitions
  if (isTRUE(use_rna_velocity)) {
    vel_candidates <- unique(c(
      if (!is.null(nonlinear_reduction)) paste0(vkey, "_", nonlinear_reduction),
      paste0(vkey, "_velocity_embedding"),
      vkey
    ))
    vel_key <- vel_candidates[vel_candidates %in% names(srt@reductions)][1]
    if (!vel_key %in% names(srt@reductions)) {
      log_message(
        "Velocity embedding for {.val {vkey}} not found; run {.fn RunSCVELO} first",
        message_type = "warning",
        verbose = verbose
      )
    } else {
      vel_emb <- srt@reductions[[vel_key]]@cell.embeddings
      storage.mode(vel_emb) <- "double"
      vel_trans <- paga_velocity_transitions_cpp(
        velocity_embedding = vel_emb,
        knn_idx = knn[["idx"]],
        groups = as.integer(groups),
        n_groups = nlevels(groups)
      )
      srt@misc[["paga"]]$velocity_transitions <- vel_trans[["transitions_confidence"]]
      srt@misc[["paga"]]$velocity_transitions_tree <- vel_trans[["transitions_confidence_tree"]]
      srt@misc[["paga"]]$velocity_group_sizes <- vel_trans[["group_sizes"]]
      names(vel_trans[["group_sizes"]]) <- group_levels
      log_message(
        "{.pkg PAGA} velocity transitions computed",
        message_type = "success",
        verbose = verbose
      )
    }
  }

  # PAGA root cell selection
  if (isTRUE(infer_pseudotime) && !is.null(nonlinear_reduction) && nonlinear_reduction %in% names(srt@reductions)) {
    root_grp_idx <- if (!is.null(root_group)) match(root_group, group_levels) else 1L
    if (anyNA(root_grp_idx) || length(root_grp_idx) == 0L) root_grp_idx <- 1L
    emb_r <- srt@reductions[[nonlinear_reduction]]@cell.embeddings
    storage.mode(emb_r) <- "double"
    root_cells <- paga_root_cell_cpp(
      embedding = emb_r,
      groups = as.integer(groups),
      root_group = as.integer(root_grp_idx[1])
    )
    srt@misc[["paga"]]$root_cells <- root_cells
    log_message(
      "PAGA root cell candidates: {.val {head(root_cells, 5)}}",
      message_type = "success",
      verbose = verbose
    )
  }

  log_message(
    "{.pkg PAGA} cpp backend completed",
    message_type = "success",
    verbose = verbose
  )
  srt
}

paga_layout_igraph <- function(
  connectivities,
  connectivities_tree = NULL,
  layout = "fr",
  threshold = 0.1
) {
  adj <- as.matrix(connectivities)
  if (!is.null(threshold) && is.finite(threshold) && threshold > 0) {
    adj[adj < threshold] <- 0
  }
  diag(adj) <- 0
  n_groups <- nrow(adj)
  if (n_groups == 1L) {
    pos <- matrix(c(0.5, 0.5), ncol = 2)
    dimnames(pos) <- NULL
    return(pos)
  }

  layout <- layout %||% "fr"
  seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (seed_exists) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit(
    {
      if (seed_exists) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
  set.seed(0)

  graph <- igraph::graph_from_adjacency_matrix(
    adj,
    mode = "max",
    weighted = TRUE,
    diag = FALSE
  )
  weights <- igraph::E(graph)$weight
  pos <- switch(layout,
    fr = {
      init <- matrix(stats::runif(n_groups * 2L), ncol = 2L)
      igraph::layout_with_fr(
        graph = graph,
        coords = init,
        weights = weights,
        niter = 500
      )
    },
    circle = igraph::layout_in_circle(graph),
    kk = igraph::layout_with_kk(graph, weights = weights),
    rt = paga_layout_tree(connectivities_tree, circular = FALSE),
    rt_circular = paga_layout_tree(connectivities_tree, circular = TRUE),
    igraph::layout_with_fr(
      graph = graph,
      coords = matrix(stats::runif(n_groups * 2L), ncol = 2L),
      weights = weights,
      niter = 500
    )
  )
  pos <- as.matrix(pos)
  pos[, 2] <- -pos[, 2]
  dimnames(pos) <- NULL
  pos
}

paga_layout_tree <- function(connectivities_tree, circular = FALSE) {
  if (is.null(connectivities_tree)) {
    log_message(
      "{.arg connectivities_tree} is required for tree PAGA layouts",
      message_type = "error"
    )
  }
  tree <- as.matrix(connectivities_tree)
  tree <- pmax(tree, t(tree))
  graph <- igraph::graph_from_adjacency_matrix(
    tree,
    mode = "max",
    weighted = TRUE,
    diag = FALSE
  )
  igraph::layout_as_tree(
    graph = graph,
    root = 1L,
    circular = circular
  )
}
