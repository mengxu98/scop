#' Run PAGA analysis
#'
#' PAGA is a graph-based method used to infer cellular trajectories.
#' This function runs the PAGA analysis on a Seurat object.
#'
#' @param srt A Seurat object.
#' @param assay_x Assay to convert as the main data matrix (X) in the anndata object.
#' @param layer_x Layer name for assay_x in the Seurat object.
#' @param assay_y Assays to convert as layers in the anndata object.
#' @param layer_y Layer names for the assay_y in the Seurat object.
#' @param adata An anndata object.
#' @param group_by Variable to use for grouping cells in the Seurat object.
#' @param linear_reduction Linear reduction method to use, e.g., "PCA".
#' @param nonlinear_reduction Non-linear reduction method to use, e.g., "UMAP".
#' @param basis The basis to use for reduction, e.g., "UMAP".
#' @param n_pcs Number of principal components to use for linear reduction. Default is 30.
#' @param n_neighbors Number of neighbors to use for constructing the KNN graph. Default is 30.
#' @param use_rna_velocity Whether to use RNA velocity for PAGA analysis. Default is FALSE.
#' @param vkey The name of the RNA velocity data to use if \code{use_rna_velocity} is TRUE. Default is "stochastic".
#' @param embedded_with_PAGA Whether to embed data using PAGA layout. Default is FALSE.
#' @param paga_layout The layout for plotting PAGA graph. See See \href{https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.paga.html}{layout} param in scanpy.pl.paga function.
#' @param threshold The threshold for plotting PAGA graph. Edges for weights below this threshold will not draw.
#' @param point_size The point size for plotting.
#' @param infer_pseudotime Whether to infer pseudotime.
#' @param root_group The group to use as the root for pseudotime inference.
#' @param root_cell The cell to use as the root for pseudotime inference.
#' @param n_dcs TThe number of diffusion components to use for pseudotime inference.
#' @param n_branchings Number of branchings to detect.
#' @param min_group_size The minimum size of a group (as a fraction of the total number of cells) to consider it as a potential branching point.
#' @param palette The palette to use for coloring cells.
#' @param palcolor A vector of colors to use as the palette.
#' @param show_plot Whether to show the PAGA plot.
#' @param dpi The DPI (dots per inch) for saving the PAGA plot.
#' @param save Whether to save the PAGA plots.
#' @param dirpath The directory to save the PAGA plots.
#' @param fileprefix The file prefix to use for the PAGA plots.
#' @param return_seurat Whether to return a Seurat object instead of an anndata object. Default is TRUE.
#'
#' @seealso \code{\link{srt_to_adata}} \code{\link{PAGAPlot}} \code{\link{CellDimPlot}} \code{\link{RunSCVELO}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("pancreas_sub")
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
#' PAGAPlot(pancreas_sub, reduction = "UMAP")
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
#' head(pancreas_sub[[]])
#' names(pancreas_sub@reductions)
#' FeatureDimPlot(
#'   pancreas_sub,
#'   features = "dpt_pseudotime",
#'   reduction = "PAGAUMAP2D"
#' )
#' PAGAPlot(pancreas_sub, reduction = "PAGAUMAP2D")
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "PAGAUMAP2D",
#'   paga = pancreas_sub@misc$paga
#' )
#' }
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
  check_python("scanpy")
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
  adata <- do.call(scop_analysis$PAGA, args)

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

#' PAGA plot
#'
#' This function generates a PAGA plot based on the given Seurat object and PAGA result.
#'
#' @param srt A Seurat object containing a PAGA result.
#' @param paga The PAGA result from the Seurat object. Defaults to \code{srt\$misc\$paga}.
#' @param type The type of plot to generate. Possible values are "connectivities" (default) and "connectivities_tree".
#' @param reduction The type of reduction to use for the plot. Defaults to the default reduction in the Seurat object.
#' @param dims The dimensions of the reduction to use for the plot. Defaults to the first two dimensions.
#' @param cells The cells to include in the plot. Defaults to all cells.
#' @param show_transition A logical value indicating whether to display transitions between different cell states. Defaults to \code{FALSE}.
#' @param node_palette The color palette to use for node coloring. Defaults to "Paired".
#' @param node_palcolor A vector of colors to use for node coloring. Defaults to \code{NULL}.
#' @param node_size The size of the nodes in the plot. Defaults to 4.
#' @param node_alpha The transparency of the nodes in the plot. Defaults to 1.
#' @param node_highlight The group(s) to highlight in the plot. Defaults to \code{NULL}.
#' @param node_highlight_color The color to use for highlighting the nodes. Defaults to "red".
#' @param label A logical value indicating whether to display labels for the nodes. Defaults to \code{FALSE}.
#' @param label.size The size of the labels. Defaults to 3.5.
#' @param label.fg The color of the label text. Defaults to "white".
#' @param label.bg The background color of the labels. Defaults to "black".
#' @param label.bg.r The transparency of the label background color. Defaults to 0.1.
#' @param label_insitu A logical value indicating whether to use in-situ labeling for the nodes. Defaults to \code{FALSE}.
#' @param label_repel A logical value indicating whether to use repel mode for labeling nodes. Defaults to \code{FALSE}.
#' @param label_repulsion The repulsion factor for repel mode. Defaults to 20.
#' @param label_point_size The size of the points in the labels. Defaults to 1.
#' @param label_point_color The color of the points in the labels. Defaults to "black".
#' @param label_segment_color The color of the lines connecting the points in the labels. Defaults to "black".
#' @param edge_threshold The threshold for including edges in the plot. Defaults to 0.01.
#' @param edge_line The type of line to use for the edges. Possible values are "straight" and "curved". Defaults to "straight".
#' @param edge_line_curvature The curvature factor for curved edges. Defaults to 0.3.
#' @param edge_line_angle The angle for curved edges. Defaults to 90.
#' @param edge_size The size range for the edges. Defaults to c(0.2, 1).
#' @param edge_color The color of the edges. Defaults to "grey40".
#' @param edge_alpha The transparency of the edges. Defaults to 0.5.
#' @param edge_shorten The amount to shorten the edges. Defaults to 0.
#' @param edge_offset The offset for curved edges. Defaults to 0.
#' @param edge_highlight The group(s) to highlight in the plot. Defaults to \code{NULL}.
#' @param edge_highlight_color The color to use for highlighting the edges. Defaults to "red".
#' @param transition_threshold The threshold for including transitions in the plot. Defaults to 0.01.
#' @param transition_line The type of line to use for the transitions. Possible values are "straight" and "curved". Defaults to "straight".
#' @param transition_line_curvature The curvature factor for curved transitions. Defaults to 0.3.
#' @param transition_line_angle The angle for curved transitions. Defaults to 90.
#' @param transition_size The size range for the transitions. Defaults to c(0.2, 1).
#' @param transition_color The color of the transitions. Defaults to "black".
#' @param transition_alpha The transparency of the transitions. Defaults to 1.
#' @param transition_arrow_type The type of arrow to use for the transitions. Possible values are "closed", "open", and "triangle". Defaults to "closed".
#' @param transition_arrow_angle The angle of the arrow for transitions. Defaults to 20.
#' @param transition_arrow_length The length of the arrow for transitions. Defaults to unit(0.02, "npc").
#' @param transition_shorten The amount to shorten the transitions. Defaults to 0.05.
#' @param transition_offset The offset for curved transitions. Defaults to 0.
#' @param transition_highlight The group(s) to highlight in the plot. Defaults to \code{NULL}.
#' @param transition_highlight_color The color to use for highlighting the transitions. Defaults to "red".
#' @param aspect.ratio The aspect ratio of the plot. Defaults to 1.
#' @param title The title of the plot. Defaults to "PAGA".
#' @param subtitle The subtitle of the plot. Defaults to \code{NULL}.
#' @param xlab The label for the x-axis. Defaults to \code{NULL}.
#' @param ylab The label for the y-axis. Defaults to \code{NULL}.
#' @param legend.position The position of the legend. Possible values are "right", "left", "bottom", and "top". Defaults to "right".
#' @param legend.direction The direction of the legend. Possible values are "vertical" and "horizontal". Defaults to "vertical".
#' @param theme_use The name of the theme to use for the plot. Defaults to "theme_scop".
#' @param theme_args A list of arguments to pass to the theme function. Defaults to an empty list.
#' @param return_layer A logical value indicating whether to return the plot as a ggplot2 layer. Defaults to \code{FALSE}.
#'
#' @seealso \code{\link{RunPAGA}} \code{\link{CellDimPlot}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("pancreas_sub")
#' pancreas_sub <- RunPAGA(
#'   srt = pancreas_sub,
#'   group_by = "SubCellType",
#'   linear_reduction = "PCA",
#'   nonlinear_reduction = "UMAP",
#'   return_seurat = TRUE
#' )
#' PAGAPlot(pancreas_sub)
#' PAGAPlot(
#'   pancreas_sub,
#'   type = "connectivities_tree"
#' )
#' PAGAPlot(
#'   pancreas_sub,
#'   reduction = "PCA"
#' )
#' PAGAPlot(
#'   pancreas_sub,
#'   reduction = "PAGAUMAP2D"
#' )
#' PAGAPlot(
#'   pancreas_sub,
#'   edge_shorten = 0.05
#' )
#' PAGAPlot(
#'   pancreas_sub,
#'   label = TRUE
#' )
#' PAGAPlot(
#'   pancreas_sub,
#'   label = TRUE,
#'   label_insitu = TRUE
#' )
#' PAGAPlot(
#'   pancreas_sub,
#'   label = TRUE,
#'   label_insitu = TRUE,
#'   label_repel = TRUE
#' )
#' PAGAPlot(
#'   pancreas_sub,
#'   edge_line = "curved"
#' )
#' PAGAPlot(
#'   pancreas_sub,
#'   node_size = "GroupSize"
#' )
#' PAGAPlot(
#'   pancreas_sub,
#'   node_highlight = "Ductal"
#' )
#' PAGAPlot(
#'   pancreas_sub,
#'   edge_highlight = paste(
#'     "Pre-endocrine",
#'     levels(pancreas_sub$SubCellType),
#'     sep = "-"
#'   )
#' )
#'
#' pancreas_sub <- RunSCVELO(
#'   srt = pancreas_sub,
#'   group_by = "SubCellType",
#'   linear_reduction = "PCA",
#'   nonlinear_reduction = "UMAP",
#'   return_seurat = TRUE
#' )
#' PAGAPlot(
#'   pancreas_sub,
#'   show_transition = TRUE
#' )
#' PAGAPlot(
#'   pancreas_sub,
#'   show_transition = TRUE,
#'   transition_offset = 0.02
#' )
#' }
PAGAPlot <- function(
    srt,
    paga = srt@misc$paga,
    type = "connectivities",
    reduction = NULL,
    dims = c(1, 2),
    cells = NULL,
    show_transition = FALSE,
    node_palette = "Paired",
    node_palcolor = NULL,
    node_size = 4,
    node_alpha = 1,
    node_highlight = NULL,
    node_highlight_color = "red",
    label = FALSE,
    label.size = 3.5,
    label.fg = "white",
    label.bg = "black",
    label.bg.r = 0.1,
    label_insitu = FALSE,
    label_repel = FALSE,
    label_repulsion = 20,
    label_point_size = 1,
    label_point_color = "black",
    label_segment_color = "black",
    edge_threshold = 0.01,
    edge_line = c("straight", "curved"),
    edge_line_curvature = 0.3,
    edge_line_angle = 90,
    edge_size = c(0.2, 1),
    edge_color = "grey40",
    edge_alpha = 0.5,
    edge_shorten = 0,
    edge_offset = 0,
    edge_highlight = NULL,
    edge_highlight_color = "red",
    transition_threshold = 0.01,
    transition_line = c("straight", "curved"),
    transition_line_curvature = 0.3,
    transition_line_angle = 90,
    transition_size = c(0.2, 1),
    transition_color = "black",
    transition_alpha = 1,
    transition_arrow_type = "closed",
    transition_arrow_angle = 20,
    transition_arrow_length = grid::unit(0.02, "npc"),
    transition_shorten = 0.05,
    transition_offset = 0,
    transition_highlight = NULL,
    transition_highlight_color = "red",
    aspect.ratio = 1,
    title = "PAGA",
    subtitle = NULL,
    xlab = NULL,
    ylab = NULL,
    legend.position = "right",
    legend.direction = "vertical",
    theme_use = "theme_scop",
    theme_args = list(),
    return_layer = FALSE) {
  if (is.null(paga)) {
    stop("Cannot find the paga result.")
  }
  if (type == "connectivities_tree") {
    use_triangular <- "both"
    edge_threshold <- 0
  } else {
    use_triangular <- "upper"
  }
  connectivities <- paga[[type]]
  transition <- paga[["transitions_confidence"]]
  groups <- paga[["groups"]]
  if (!is.factor(srt@meta.data[[groups]])) {
    srt@meta.data[[groups]] <- factor(srt@meta.data[[groups]])
  }
  if (nlevels(srt@meta.data[[groups]]) != nrow(connectivities)) {
    stop("nlevels in ", groups, " is not identical with the group in paga")
  }
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  if (!reduction %in% names(srt@reductions)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }
  reduction_key <- srt@reductions[[reduction]]@key
  dat_dim <- as.data.frame(srt@reductions[[reduction]]@cell.embeddings)
  colnames(dat_dim) <- paste0(reduction_key, seq_len(ncol(dat_dim)))
  rownames(dat_dim) <- rownames(dat_dim) %||% colnames(srt@assays[[1]])
  dat_dim <- dat_dim[, paste0(reduction_key, dims)]
  dat_dim[[groups]] <- srt@meta.data[rownames(dat_dim), groups]
  dat <- stats::aggregate(
    dat_dim[, paste0(reduction_key, dims)],
    by = list(dat_dim[[groups]]),
    FUN = stats::median
  )
  colnames(dat)[1] <- groups
  rownames(dat) <- dat[[groups]]
  dat[["GroupSize"]] <- as.numeric(table(dat_dim[[groups]])[rownames(dat)])
  colnames(connectivities) <- rownames(connectivities) <- rownames(dat)
  if (!is.null(transition)) {
    colnames(transition) <- rownames(transition) <- rownames(dat)
  }
  if (!isTRUE(show_transition)) {
    transition <- NULL
  } else if (isTRUE(show_transition) && is.null(transition)) {
    warning(
      "transitions_confidence need to be calculated first.",
      immediate. = TRUE
    )
  }
  xlab <- xlab %||% paste0(reduction_key, dims[1])
  ylab <- ylab %||% paste0(reduction_key, dims[2])

  if (!is.null(cells)) {
    dat <- dat[intersect(rownames(dat), cells), , drop = FALSE]
    connectivities <- connectivities[rownames(dat), rownames(dat)]
  }

  out <- GraphPlot(
    node = dat,
    edge = Matrix::as.matrix(connectivities),
    node_coord = paste0(reduction_key, dims),
    node_group = groups,
    node_palette = node_palette,
    node_palcolor = node_palcolor,
    node_size = node_size,
    node_alpha = node_alpha,
    node_highlight = node_highlight,
    node_highlight_color = node_highlight_color,
    label = label,
    label.size = label.size,
    label.fg = label.fg,
    label.bg = label.bg,
    label.bg.r = label.bg.r,
    label_insitu = label_insitu,
    label_repel = label_repel,
    label_repulsion = label_repulsion,
    label_point_size = label_point_size,
    label_point_color = label_point_color,
    label_segment_color = label_segment_color,
    edge_threshold = edge_threshold,
    use_triangular = use_triangular,
    edge_line = edge_line,
    edge_line_curvature = edge_line_curvature,
    edge_line_angle = edge_line_angle,
    edge_size = edge_size,
    edge_color = edge_color,
    edge_alpha = edge_alpha,
    edge_shorten = edge_shorten,
    edge_offset = edge_offset,
    edge_highlight = edge_highlight,
    edge_highlight_color = edge_highlight_color,
    transition = transition,
    transition_threshold = transition_threshold,
    transition_line = transition_line,
    transition_line_curvature = transition_line_curvature,
    transition_line_angle = transition_line_angle,
    transition_color = transition_color,
    transition_size = transition_size,
    transition_alpha = transition_alpha,
    transition_arrow_type = transition_arrow_type,
    transition_arrow_angle = transition_arrow_angle,
    transition_arrow_length = transition_arrow_length,
    transition_shorten = transition_shorten,
    transition_offset = transition_offset,
    transition_highlight = transition_highlight,
    transition_highlight_color = transition_highlight_color,
    aspect.ratio = aspect.ratio,
    title = title,
    subtitle = subtitle,
    xlab = xlab,
    ylab = ylab,
    legend.position = legend.position,
    legend.direction = legend.direction,
    theme_use = theme_use,
    theme_args = theme_args,
    return_layer = return_layer
  )
  return(out)
}
