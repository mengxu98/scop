#' @title Run CellRank analysis with kernel-estimator architecture
#'
#' @description
#' CellRank is a powerful toolkit for studying cellular dynamics using Markov state modeling.
#' This function implements the modern kernel-estimator architecture recommended by CellRank,
#' which provides more flexibility and advanced features compared to the legacy API.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt A Seurat object. Default is `NULL`.
#' If provided, `adata` will be ignored.
#' @param adata An anndata object. Default is `NULL`.
#' @param assay_x Assay to convert in the anndata object.
#' @param layer_x Layer name for `assay_x` in the Seurat object.
#' @param assay_y Assay to convert in the anndata object.
#' @param layer_y Layer names for the `assay_y` in the Seurat object.
#' @param group_by Variable to use for grouping cells in the Seurat object.
#' @param linear_reduction Linear reduction method to use, e.g., `"PCA"`.
#' @param nonlinear_reduction Non-linear reduction method to use, e.g., `"UMAP"`.
#' @param basis The basis to use for reduction, e.g., `"UMAP"`.
#' @param n_pcs Number of principal components to use for linear reduction.
#' Default is `30`.
#' @param n_neighbors Number of neighbors to use for constructing the KNN graph.
#' Default is `30`.
#' @param cores The number of cores to use for parallelization with \link[foreach:foreach]{foreach::foreach}.
#' Default is `1`.
#' @param palette The palette to use for coloring cells.
#' @param palcolor A vector of colors to use as the palette.
#' @param legend.position Position of legend in plots. Can be `"on data"`, `"right margin"`, `"bottom right"`, etc. Default is `"on data"`.
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
#' @param t power to which the diffusion operator is powered for `magic.MAGIC`.
#' Default is `2`.
#' @param min_shared_counts Minimum number of counts (both unspliced and spliced) required for a gene.
#' Default is `30`.
#' @param stream_smooth Multiplication factor for scale in Gaussian kernel around grid point.
#' @param stream_density Controls the closeness of streamlines.
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
#' or `"cytotrace"` (auto-computes CytoTRACE score, suitable for RNA-only data).
#' @param time_key Key in metadata for pseudotime. Used when `kernel_type = "pseudotime"`.
#' If the key doesn't exist, DPT pseudotime will be computed automatically.
#' Default is `"dpt_pseudotime"`.
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
#' Lineage pseudotime columns are automatically added to metadata with prefix `"Lineage_"`
#' (e.g., `Lineage_Alpha`, `Lineage_Beta`), compatible with [LineagePlot], [RunDynamicFeatures],
#' and [RunDynamicEnrichment]. These are computed by combining base pseudotime
#' (`cytotrace_pseudotime`, `dpt_pseudotime`, or `latent_time`) with fate probabilities.
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
#'   group_by = "SubCellType",
#'   cores = 6
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
RunCellRank <- function(
    srt = NULL,
    assay_x = "RNA",
    layer_x = "counts",
    assay_y = c("spliced", "unspliced"),
    layer_y = "counts",
    adata = NULL,
    group_by = NULL,
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
    kernel_type = c("velocity", "pseudotime", "cytotrace"),
    time_key = "dpt_pseudotime",
    estimator_type = c("GPCCA", "CFLARE"),
    use_connectivity_kernel = TRUE,
    velocity_weight = 0.8,
    connectivity_weight = 0.2,
    softmax_scale = 4,
    n_macrostates = NULL,
    schur_method = c("krylov", "brandts"),
    n_cells_terminal = 10,
    show_plot = TRUE,
    save_plot = FALSE,
    plot_format = c("pdf", "png", "svg"),
    plot_dpi = 300,
    plot_prefix = "cellrank",
    legend.position = "on data",
    palette = "Paired",
    palcolor = NULL,
    dirpath = "./cellrank",
    return_seurat = !is.null(srt),
    verbose = TRUE) {
  PrepareEnv()
  check_python("cellrank", verbose = verbose)
  if (isTRUE(magic_impute)) {
    check_python("magic-impute", verbose = verbose)
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

  kernel_type <- match.arg(kernel_type)
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
  groups <- py_to_r2(args[["adata"]]$obs)[[group_by]]
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
        "plot_prefix"
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

      return(merged)
    }
  } else {
    adata$uns["cellrank_estimator"] <- estimator
    adata$uns["cellrank_kernel"] <- kernel
    return(adata)
  }
}
