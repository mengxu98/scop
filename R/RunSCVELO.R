#' Run scVelo workflow
#'
#' scVelo is a scalable toolkit for RNA velocity analysis in single cells. This function runs scVelo workflow on a Seurat object.
#'
#' @inheritParams RunPAGA
#' @param mode Velocity estimation model to use, e.g., "stochastic".
#' @param fitting_by Method used to fit gene velocities, e.g., "stochastic".
#' @param magic_impute Flag indicating whether to perform magic imputation.
#' @param knn The number of nearest neighbors for magic.MAGIC.
#' @param t power to which the diffusion operator is powered for magic.MAGIC.
#' @param min_shared_counts Minimum number of counts (both unspliced and spliced) required for a gene.
#' @param n_pcs Number of principal components (PCs) used for velocity estimation.
#' @param n_neighbors Number of nearest neighbors used for velocity estimation.
#' @param stream_smooth Multiplication factor for scale in Gaussian kernel around grid point.
#' @param stream_density Controls the closeness of streamlines. When density = 2 (default), the domain is divided into a 60x60 grid, whereas density linearly scales this grid. Each cell in the grid can have, at most, one traversing streamline.
#' @param arrow_length Length of arrows.
#' @param arrow_size Size of arrows.
#' @param arrow_density Amount of velocities to show.
#' @param denoise Boolean flag indicating whether to denoise.
#' @param denoise_topn Number of genes with highest likelihood selected to infer velocity directions.
#' @param kinetics Boolean flag indicating whether to estimate RNA kinetics.
#' @param kinetics_topn Number of genes with highest likelihood selected to infer velocity directions.
#' @param calculate_velocity_genes Boolean flag indicating whether to calculate velocity genes.
#' @param top_n The number of top features to plot.
#' @param n_jobs The number of parallel jobs to run.
#'
#' @seealso \code{\link{srt_to_adata}} \code{\link{VelocityPlot}} \code{\link{CellDimPlot}} \code{\link{RunPAGA}}
#'
#' @examples
#' \dontrun{
#' data("pancreas_sub")
#' pancreas_sub <- RunSCVELO(
#'   srt = pancreas_sub,
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
#' pancreas_sub <- standard_scop(
#'   pancreas_sub,
#'   normalization_method = "SCT",
#'   nonlinear_reduction = "tsne"
#' )
#' pancreas_sub <- RunSCVELO(
#'   srt = pancreas_sub,
#'   assay_x = "SCT",
#'   group_by = "SubCellType",
#'   linear_reduction = "Standardpca",
#'   nonlinear_reduction = "StandardTSNE2D"
#' )
#' }
#' @export
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
    stop("One of 'srt', 'adata' must be provided.")
  }
  if (is.null(group_by)) {
    stop("'group_by' must be provided.")
  }
  if (is.null(linear_reduction) && is.null(nonlinear_reduction)) {
    stop(
      "'linear_reduction' or 'nonlinear_reduction' must be provided at least one."
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
  groups <- py_to_r_auto(args[["adata"]]$obs)[[group_by]]
  args[["palette"]] <- palette_scop(
    levels(groups) %||% unique(groups),
    palette = palette,
    palcolor = palcolor
  )

  scop_analysis <- reticulate::import_from_path(
    "scop_analysis",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )
  adata <- do.call(scop_analysis$SCVELO, args)

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      srt_out1 <- SrtAppend(srt_raw = srt, srt_append = srt_out)
      srt_out2 <- SrtAppend(
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
