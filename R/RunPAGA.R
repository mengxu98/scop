#' @title Run PAGA analysis
#'
#' @description
#' PAGA is a graph-based method used to infer cellular trajectories.
#' This function runs the PAGA analysis on a Seurat object.
#'
#' @md
#' @param srt A Seurat object.
#' @param assay_x Assay to convert in the anndata object.
#' @param layer_x Layer name for `assay_x` in the Seurat object.
#' @param assay_y Assay to convert in the anndata object.
#' @param layer_y Layer names for the `assay_y` in the Seurat object.
#' @param adata An anndata object.
#' @param group_by Variable to use for grouping cells in the Seurat object.
#' @param linear_reduction Linear reduction method to use, e.g., `"PCA"`.
#' @param nonlinear_reduction Non-linear reduction method to use, e.g., `"UMAP"`.
#' @param basis The basis to use for reduction, e.g., `"UMAP"`.
#' @param n_pcs Number of principal components to use for linear reduction.
#' Default is `30`.
#' @param n_neighbors Number of neighbors to use for constructing the KNN graph.
#' Default is `30`.
#' @param use_rna_velocity Whether to use RNA velocity for PAGA analysis.
#' Default is `FALSE`.
#' @param vkey The name of the RNA velocity data to use if `use_rna_velocity` is `TRUE`.
#' Default is `"stochastic"`.
#' @param embedded_with_PAGA Whether to embed data using PAGA layout.
#' Default is `FALSE`.
#' @param paga_layout The layout for plotting PAGA graph.
#' See \href{https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.paga.html}{layout} param in `scanpy.pl.paga` function.
#' @param threshold The threshold for plotting PAGA graph.
#' Edges for weights below this threshold will not be drawn.
#' @param point_size The point size for plotting.
#' @param infer_pseudotime Whether to infer pseudotime.
#' @param root_group The group to use as the root for pseudotime inference.
#' @param root_cell The cell to use as the root for pseudotime inference.
#' @param n_dcs The number of diffusion components to use for pseudotime inference.
#' @param n_branchings Number of branchings to detect.
#' @param min_group_size The minimum size of a group (as a fraction of the total number of cells) to consider it as a potential branching point.
#' @param palette The palette to use for coloring cells.
#' @param palcolor A vector of colors to use as the palette.
#' @param show_plot Whether to show the plot.
#' @param dpi The DPI (dots per inch) for saving the plot.
#' @param save Whether to save the plots.
#' @param dirpath The directory to save the plots.
#' @param fileprefix The file prefix to use for the plots.
#' @param return_seurat Whether to return a Seurat object instead of an anndata object.
#' Default is `TRUE`.
#'
#' @seealso
#' [srt_to_adata], [PAGAPlot], [CellDimPlot], [RunSCVELO]
#'
#' @export
#'
#' @examples
#' PrepareEnv()
#' data(pancreas_sub)
#' pancreas_sub <- RunPAGA(
#'   srt = pancreas_sub,
#'   assay_x = "RNA",
#'   group_by = "SubCellType",
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
#'   srt = pancreas_sub,
#'   group_by = "SubCellType",
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
RunPAGA <- function(
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
    show_plot = TRUE,
    save = FALSE,
    dpi = 300,
    dirpath = "./",
    fileprefix = "",
    return_seurat = !is.null(srt)) {
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
  params <- c(
    "srt",
    "assay_x",
    "layer_x",
    "assay_y",
    "layer_y",
    "return_seurat",
    "palette",
    "palcolor"
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
  groups <- py_to_r2(args[["adata"]]$obs)[[group_by]]
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
  log_message("Running {.pkg PAGA} analysis...")
  adata <- do.call(scop_analysis$PAGA, args)
  log_message(
    "{.pkg PAGA} analysis completed",
    message_type = "success"
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
