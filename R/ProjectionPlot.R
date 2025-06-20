#' Projection Plot
#'
#' This function generates a projection plot, which can be used to compare two groups of cells in a dimensionality reduction space.
#'
#' @param srt_query An object of class Seurat storing the query group cells.
#' @param srt_ref An object of class Seurat storing the reference group cells.
#' @param query_group The grouping variable for the query group cells.
#' @param ref_group The grouping variable for the reference group cells.
#' @param query_reduction The name of the reduction in the query group cells.
#' @param ref_reduction The name of the reduction in the reference group cells.
#' @param query_param A list of parameters for customizing the query group plot. Available parameters: palette (color palette for groups) and cells.highlight (whether to highlight cells).
#' @param ref_param A list of parameters for customizing the reference group plot. Available parameters: palette (color palette for groups) and cells.highlight (whether to highlight cells).
#' @param xlim The x-axis limits for the plot. If not provided, the limits will be calculated based on the data.
#' @param ylim The y-axis limits for the plot. If not provided, the limits will be calculated based on the data.
#' @param pt.size The size of the points in the plot.
#' @param stroke.highlight The size of the stroke highlight for cells.
#'
#' @export
#'
#' @examples
#' data("panc8_sub")
#' srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
#' srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
#' srt_ref <- integration_scop(
#'   srt_ref,
#'   batch = "tech",
#'   integration_method = "Seurat"
#' )
#' CellDimPlot(
#'   srt_ref,
#'   group.by = c("celltype", "tech")
#' )
#'
#' # Projection
#' srt_query <- RunKNNMap(
#'   srt_query = srt_query,
#'   srt_ref = srt_ref,
#'   ref_umap = "SeuratUMAP2D"
#' )
#' ProjectionPlot(
#'   srt_query = srt_query,
#'   srt_ref = srt_ref,
#'   query_group = "celltype",
#'   ref_group = "celltype"
#' )
ProjectionPlot <- function(
    srt_query,
    srt_ref,
    query_group = NULL,
    ref_group = NULL,
    query_reduction = "ref.embeddings",
    ref_reduction = srt_query[[query_reduction]]@misc[["reduction.model"]] %||%
      NULL,
    query_param = list(palette = "Set1", cells.highlight = TRUE),
    ref_param = list(palette = "Paired"),
    xlim = NULL,
    ylim = NULL,
    pt.size = 0.8,
    stroke.highlight = 0.5) {
  if (is.null(ref_reduction)) {
    log_message(
      "Please specify the ref_reduction.",
      message_type = "error"
    )
  }
  query_param[["show_stat"]] <- FALSE
  ref_param[["show_stat"]] <- FALSE

  if (is.null(xlim)) {
    ref_xlim <- range(srt_ref[[ref_reduction]]@cell.embeddings[, 1])
    query_xlim <- range(srt_query[[query_reduction]]@cell.embeddings[, 1])
    xlim <- range(c(ref_xlim, query_xlim))
  }
  if (is.null(ylim)) {
    ref_ylim <- range(srt_ref[[ref_reduction]]@cell.embeddings[, 2])
    query_ylim <- range(srt_query[[query_reduction]]@cell.embeddings[, 2])
    ylim <- range(c(ref_ylim, query_ylim))
  }

  p1 <- do.call(
    CellDimPlot,
    args = c(
      srt = srt_ref,
      reduction = ref_reduction,
      group.by = ref_group,
      ref_param
    )
  ) +
    guides(
      color = guide_legend(
        title = paste0("Ref: ", ref_group),
        override.aes = list(size = 4)
      )
    )
  p1legend <- get_legend(
    p1 + theme(legend.position = "bottom")
  )

  p2 <- do.call(
    CellDimPlot,
    args = c(
      srt = srt_query,
      reduction = query_reduction,
      group.by = query_group,
      query_param
    )
  ) +
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(limits = ylim)
  p2data <- ggplot_build(p2)$data[[1]]
  color <- p2data$colour
  names(color) <- p2$data$group.by
  p2 <- p2 +
    guides(
      color = guide_legend(
        title = paste0("Query: ", query_group),
        override.aes = list(
          size = 4,
          shape = 21,
          color = "black",
          fill = stats::na.omit(color[levels(p2$data$group.by)])
        )
      )
    )
  p2legend <- get_legend(
    p2 + theme(legend.position = "bottom")
  )

  if (!is.null(p1legend) && !is.null(p2legend)) {
    legend <- cbind(p1legend, p2legend)
  } else {
    legend <- p1legend %||% p2legend
  }

  p3 <- p1 +
    ggnewscale::new_scale_fill() +
    ggnewscale::new_scale_color() +
    geom_point(
      data = p2data,
      aes(x = x, y = y),
      color = "black",
      size = pt.size + stroke.highlight
    ) +
    geom_point(
      data = p2data,
      aes(x = x, y = y, color = colour),
      size = pt.size
    ) +
    scale_color_identity() +
    facet_null() +
    theme(legend.position = "none")
  if (is.null(legend)) {
    return(p3)
  } else {
    gtable <- as_grob(p3)
    gtable <- add_grob(gtable, legend, "right")
    p <- patchwork::wrap_plots(gtable)
    return(p)
  }
}
