#' @title Run PAGA analysis
#'
#' @description
#' PAGA is a graph-based method used to infer cellular trajectories.
#' This function runs the PAGA analysis on a Seurat object.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams RunCellRank
#' @param backend Backend used to compute PAGA. `"cpp"` uses the native C++
#' implementation for the standard connectivity graph. `"python"` keeps the
#' original scanpy workflow for RNA-velocity transitions, PAGA-initialized
#' embeddings, plotting side effects, or DPT pseudotime.
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
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunPAGA(
#'   pancreas_sub,
#'   assay_x = "RNA",
#'   group.by = "SubCellType",
#'   linear_reduction = "PCA",
#'   nonlinear_reduction = "UMAP"
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "draw_graph_fr"
#' )
#'
#' PAGAPlot(pancreas_sub, reduction = "UMAP")
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   paga = pancreas_sub@misc$paga
#' )
#'
#' pancreas_sub <- RunPAGA(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   linear_reduction = "PCA",
#'   nonlinear_reduction = "UMAP",
#'   embedded_with_PAGA = TRUE,
#'   infer_pseudotime = TRUE,
#'   root_group = "Ductal"
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   features = "dpt_pseudotime",
#'   reduction = "PAGAUMAP2D"
#' )
#'
#' PAGAPlot(pancreas_sub, reduction = "PAGAUMAP2D")
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "PAGAUMAP2D",
#'   paga = pancreas_sub@misc$paga
#' )
#' }
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
  backend = c("cpp", "python"),
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
    unsupported_cpp <- c(
      if (isTRUE(use_rna_velocity)) "use_rna_velocity",
      if (isTRUE(embedded_with_PAGA)) "embedded_with_PAGA",
      if (isTRUE(infer_pseudotime)) "infer_pseudotime",
      if (isTRUE(show_plot)) "show_plot",
      if (isTRUE(save_plot)) "save_plot"
    )
    if (length(unsupported_cpp) > 0L) {
      unsupported_text <- paste0(unsupported_cpp, collapse = ", ")
      log_message(
        "{.arg backend = 'cpp'} currently supports the standard PAGA connectivity graph only; use {.arg backend = 'python'} for {.arg {unsupported_text}}",
        message_type = "error"
      )
    }
    srt <- run_paga_cpp(
      srt = srt,
      group.by = group.by,
      linear_reduction = linear_reduction,
      n_pcs = n_pcs,
      n_neighbors = n_neighbors,
      cores = cores,
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
  n_pcs,
  n_neighbors,
  cores,
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
  paga <- paga_connectivities_cpp(
    knn_idx = knn[["idx"]],
    groups = as.integer(groups),
    n_groups = nlevels(groups)
  )
  group_levels <- levels(groups)
  for (nm in c(
    "connectivities",
    "connectivities_tree",
    "expected_n_edges_random",
    "directed_edges"
  )) {
    dimnames(paga[[nm]]) <- list(group_levels, group_levels)
  }
  names(paga[["group_sizes"]]) <- group_levels

  srt@misc[["paga"]] <- list(
    connectivities = paga[["connectivities"]],
    connectivities_tree = paga[["connectivities_tree"]],
    groups = group.by,
    backend = "cpp",
    parameters = list(
      linear_reduction = linear_reduction,
      n_pcs = length(dims_use),
      n_neighbors = as.integer(n_neighbors),
      knn_k = knn_k
    )
  )
  srt@misc[[paste0(group.by, "_sizes")]] <- as.numeric(paga[["group_sizes"]])
  names(srt@misc[[paste0(group.by, "_sizes")]]) <- group_levels
  log_message(
    "{.pkg PAGA} cpp backend completed",
    message_type = "success",
    verbose = verbose
  )
  srt
}
