#' @title Run scVelo workflow
#'
#' @description
#' scVelo is a scalable toolkit for RNA velocity analysis in single cells.
#' This function runs an enhanced scVelo workflow on a Seurat object with improved error handling,
#' version compatibility, and modular design.
#'
#' @md
#' @inheritParams RunPAGA
#' @param mode Velocity estimation models to use.
#' Can be a vector containing `"deterministic"`, `"stochastic"`, and/or `"dynamical"`.
#' @param fitting_by Method used to fit gene velocities for dynamical modeling, e.g., "stochastic".
#' @param magic_impute Flag indicating whether to perform magic imputation.
#' @param knn The number of nearest neighbors for `magic.MAGIC`.
#' @param t power to which the diffusion operator is powered for `magic.MAGIC`.
#' @param min_shared_counts Minimum number of counts (both unspliced and spliced) required for a gene.
#' @param n_pcs Number of principal components (PCs) used for velocity estimation.
#' @param n_neighbors Number of nearest neighbors used for velocity estimation.
#' @param filter_genes Whether to filter genes based on minimum counts.
#' @param min_counts Minimum counts for gene filtering.
#' @param min_counts_u Minimum unspliced counts for gene filtering.
#' @param normalize_per_cell Whether to normalize counts per cell.
#' @param log_transform Whether to apply log transformation.
#' @param use_raw Whether to use raw data for dynamical modeling.
#' @param diff_kinetics Whether to use differential kinetics.
#' @param stream_smooth Multiplication factor for scale in Gaussian kernel around grid point.
#' @param stream_density Controls the closeness of streamlines.
#' When density = 2 (default), the domain is divided into a 60x60 grid,
#' whereas density linearly scales this grid.
#' Each cell in the grid can have, at most, one traversing streamline.
#' @param arrow_length Length of arrows.
#' @param arrow_size Size of arrows.
#' @param arrow_density Amount of velocities to show.
#' @param denoise Boolean flag indicating whether to denoise.
#' @param denoise_topn Number of genes with highest likelihood selected to infer velocity directions.
#' @param kinetics Boolean flag indicating whether to estimate RNA kinetics.
#' @param kinetics_topn Number of genes with highest likelihood selected to infer velocity directions.
#' @param calculate_velocity_genes Boolean flag indicating whether to calculate velocity genes.
#' @param compute_velocity_confidence Whether to compute velocity confidence metrics.
#' @param compute_terminal_states Whether to compute terminal states (root and end points).
#' @param compute_pseudotime Whether to compute velocity pseudotime.
#' @param compute_paga Whether to compute PAGA (Partition-based graph abstraction).
#' @param top_n The number of top features to plot.
#' @param n_jobs The number of parallel jobs to run.
#'
#' @seealso
#' [srt_to_adata], [VelocityPlot], [CellDimPlot], [RunPAGA]
#'
#' @export
#' @examples
#' \dontrun{
#' PrepareEnv()
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunSCVELO(
#'   pancreas_sub,
#'   assay_x = "RNA",
#'   group_by = "SubCellType",
#'   linear_reduction = "PCA",
#'   nonlinear_reduction = "UMAP"
#' )
#' head(pancreas_sub[[]])
#' names(pancreas_sub@assays)
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   c("stochastic_length", "stochastic_confidence")
#' )
#' FeatureDimPlot(
#'   pancreas_sub,
#'   "stochastic_pseudotime"
#' )
#' VelocityPlot(
#'   pancreas_sub,
#'   reduction = "UMAP",
#'   plot_type = "stream"
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   pt.size = NA,
#'   velocity = "stochastic"
#' )
#'
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunSCVELO(
#'   pancreas_sub,
#'   assay_x = "RNA",
#'   group_by = "SubCellType",
#'   linear_reduction = "PCA",
#'   nonlinear_reduction = "UMAP",
#'   mode = c("deterministic", "stochastic"),
#'   filter_genes = TRUE,
#'   min_counts = 5,
#'   compute_velocity_confidence = TRUE,
#'   compute_terminal_states = TRUE,
#'   compute_pseudotime = TRUE,
#'   compute_paga = TRUE
#' )
#' }
RunSCVELO <- function(
    srt = NULL,
    assay_x = "RNA",
    layer_x = "counts",
    assay_y = c("spliced", "unspliced"),
    layer_y = "counts",
    adata = NULL,
    group_by = NULL,
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
    n_jobs = 1,
    palette = "Paired",
    palcolor = NULL,
    show_plot = TRUE,
    save = FALSE,
    dpi = 300,
    dirpath = "./",
    fileprefix = "",
    return_seurat = !is.null(srt)) {
  check_python("scvelo")
  if (isTRUE(magic_impute)) {
    check_python("magic-impute")
  }

  if (all(is.null(srt), is.null(adata))) {
    log_message(
      "One of {.arg srt} or {.arg adata} must be provided",
      message_type = "error"
    )
  }
  if (is.null(group_by)) {
    log_message(
      "{.arg group_by} must be provided",
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
        "palcolor"
      )
  ]

  if (!is.null(srt)) {
    args[["adata"]] <- srt_to_adata(
      srt = srt,
      assay_x = assay_x,
      layer_x = layer_x,
      assay_y = assay_y,
      layer_y = layer_y
    )
  }

  groups <- py_to_r2(args[["adata"]]$obs)[[group_by]]
  args[["palette"]] <- palette_colors(
    levels(groups) %||% unique(groups),
    palette = palette,
    palcolor = palcolor
  )

  log_message("Running {.pkg scVelo} analysis...")
  scop_analysis <- reticulate::import_from_path(
    "scop_analysis",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )

  adata <- do.call(scop_analysis$SCVELO, args)
  log_message(
    "{.pkg scVelo} analysis completed",
    message_type = "success"
  )
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
