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
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunSCVELO(
#'   pancreas_sub,
#'   assay_x = "RNA",
#'   group.by = "SubCellType",
#'   linear_reduction = "PCA",
#'   nonlinear_reduction = "UMAP"
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   c(
#'     "stochastic_length",
#'     "stochastic_confidence",
#'     "stochastic_pseudotime"
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
#' }
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
  if (!identical(mode_use, "stochastic")) {
    log_message(
      "{.arg backend = 'cpp'} currently supports {.arg mode = 'stochastic'} only",
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
  loadings <- srt@reductions[[linear_reduction]]@feature.loadings
  if (nrow(loadings) > 0L) {
    features <- intersect(rownames(loadings), features)
  }
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

  linear_embedding <- srt@reductions[[linear_reduction]]@cell.embeddings
  dims_use <- seq_len(min(as.integer(n_pcs), ncol(linear_embedding)))
  linear_embedding <- as.matrix(linear_embedding[cells, dims_use, drop = FALSE])
  storage.mode(linear_embedding) <- "double"

  nonlinear_embedding <- as.matrix(
    srt@reductions[[nonlinear_reduction]]@cell.embeddings[cells, , drop = FALSE]
  )
  storage.mode(nonlinear_embedding) <- "double"
  knn_k <- max(1L, min(as.integer(n_neighbors) - 1L, nrow(linear_embedding) - 1L))
  log_message(
    "Running {.pkg scVelo} stochastic embedding with {.arg backend = 'cpp'} using {.val {length(features)}} features",
    verbose = verbose
  )
  knn <- run_cpp_knn(
    reference = linear_embedding,
    query = linear_embedding,
    k = knn_k,
    metric = "euclidean",
    exclude_self = TRUE,
    n_threads = as.integer(cores)
  )
  velocity <- scvelo_stochastic_embedding_cpp(
    spliced = spliced,
    unspliced = unspliced,
    knn_idx = knn[["idx"]],
    embedding = nonlinear_embedding
  )
  velocity_embedding <- velocity[["velocity_embedding"]]
  rownames(velocity_embedding) <- cells
  colnames(velocity_embedding) <- colnames(nonlinear_embedding)
  velocity_reduction <- paste0("stochastic_", nonlinear_reduction)
  srt[[velocity_reduction]] <- SeuratObject::CreateDimReducObject(
    embeddings = velocity_embedding,
    assay = spliced_assay,
    key = paste0(gsub("_", "", velocity_reduction), "_")
  )
  srt$stochastic_confidence <- as.numeric(velocity[["confidence"]])
  srt$stochastic_length <- as.numeric(velocity[["velocity_length"]])
  srt@tools[["SCVELO"]] <- list(
    backend = "cpp",
    mode = "stochastic",
    velocity_reduction = velocity_reduction,
    group.by = group.by,
    features = features,
    gamma = velocity[["gamma"]],
    parameters = list(
      spliced_assay = spliced_assay,
      unspliced_assay = unspliced_assay,
      layer_y = layer_y,
      linear_reduction = linear_reduction,
      nonlinear_reduction = nonlinear_reduction,
      n_pcs = length(dims_use),
      n_neighbors = as.integer(n_neighbors),
      knn_k = knn_k
    )
  )
  log_message(
    "{.pkg scVelo} cpp stochastic embedding completed",
    message_type = "success",
    verbose = verbose
  )
  srt
}
