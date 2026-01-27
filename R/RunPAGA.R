#' @title Run PAGA analysis
#'
#' @description
#' PAGA is a graph-based method used to infer cellular trajectories.
#' This function runs the PAGA analysis on a Seurat object.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams RunCellRank
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
#' [srt_to_adata], [PAGAPlot], [CellDimPlot], [RunSCVELO]
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
    palette = "Paired",
    palcolor = NULL,
    legend.position = "on data",
    cores = 1,
    show_plot = FALSE,
    save_plot = FALSE,
    plot_format = c("pdf", "png", "svg"),
    plot_dpi = 300,
    plot_prefix = "paga",
    dirpath = "./paga",
    return_seurat = !is.null(srt),
    verbose = TRUE) {
  PrepareEnv()

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
    "cores"
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
