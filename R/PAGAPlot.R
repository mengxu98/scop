#' @title PAGA plot
#'
#' @description
#' This function generates a PAGA plot based on the given Seurat object and PAGA result.
#'
#' @md
#' @inheritParams CellDimPlot
#' @inheritParams GraphPlot
#' @param paga The PAGA result from the Seurat object.
#' Default is `srt@misc$paga`.
#' @param type The type of plot to generate.
#' Possible values are `"connectivities"` (default) and `"connectivities_tree"`.
#' @param show_transition Whether to display transitions between different cell states.
#' Default is `FALSE`.
#' @param title The text for the title.
#' Default is `"PAGA"`.
#'
#' @seealso
#' [RunPAGA], [CellDimPlot], [GraphPlot]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunPAGA(
#'   pancreas_sub,
#'   group_by = "SubCellType",
#'   linear_reduction = "PCA",
#'   nonlinear_reduction = "UMAP",
#'   return_seurat = TRUE
#' )
#'
#' PAGAPlot(pancreas_sub)
#'
#' PAGAPlot(
#'   pancreas_sub,
#'   type = "connectivities_tree"
#' )
#'
#' PAGAPlot(
#'   pancreas_sub,
#'   reduction = "PCA"
#' )
#'
#' PAGAPlot(
#'   pancreas_sub,
#'   reduction = "PAGAUMAP2D"
#' )
#'
#' PAGAPlot(
#'   pancreas_sub,
#'   edge_shorten = 0.05
#' )
#'
#' PAGAPlot(
#'   pancreas_sub,
#'   label = TRUE
#' )
#'
#' PAGAPlot(
#'   pancreas_sub,
#'   label = TRUE,
#'   label_insitu = TRUE
#' )
#'
#' PAGAPlot(
#'   pancreas_sub,
#'   label = TRUE,
#'   label_insitu = TRUE,
#'   label_repel = TRUE
#' )
#'
#' PAGAPlot(
#'   pancreas_sub,
#'   edge_line = "curved"
#' )
#'
#' PAGAPlot(
#'   pancreas_sub,
#'   node_size = "GroupSize"
#' )
#'
#' PAGAPlot(
#'   pancreas_sub,
#'   node_highlight = "Ductal"
#' )
#'
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
#'   pancreas_sub,
#'   group_by = "SubCellType",
#'   linear_reduction = "PCA",
#'   nonlinear_reduction = "UMAP",
#'   return_seurat = TRUE
#' )
#'
#' PAGAPlot(
#'   pancreas_sub,
#'   show_transition = TRUE
#' )
#'
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
    log_message(
      "Cannot find the paga result.",
      message_type = "error"
    )
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
    log_message(
      "nlevels in ", groups, " is not identical with the group in paga",
      message_type = "error"
    )
  }
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  if (!reduction %in% names(srt@reductions)) {
    log_message(
      paste0(reduction, " is not in the srt reduction names."),
      message_type = "error"
    )
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
  if (isFALSE(show_transition)) {
    transition <- NULL
  } else if (isTRUE(show_transition) && is.null(transition)) {
    log_message(
      "transitions_confidence need to be calculated first.",
      message_type = "warning"
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
    edge = as_matrix(connectivities),
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
