#' LineagePlot
#'
#' Generate a lineage plot based on the pseudotime.
#'
#' @param srt An object of class Seurat.
#' @param lineages A character vector that specifies the lineages to be included. Typically, use the pseudotime of cells.
#' @param reduction An optional string specifying the dimensionality reduction method to use.
#' @param dims A numeric vector of length 2 specifying the dimensions to plot.
#' @param cells An optional character vector specifying the cells to include in the plot.
#' @param trim A numeric vector of length 2 specifying the quantile range of lineages to include in the plot.
#' @param span A numeric value specifying the span of the loess smoother.
#' @param palette A character string specifying the color palette to use for the lineages.
#' @param palcolor An optional string specifying the color for the palette.
#' @param lineages_arrow An arrow object specifying the arrow for lineages.
#' @param linewidth A numeric value specifying the linewidth for the lineages.
#' @param line_bg A character string specifying the color for the background lines.
#' @param line_bg_stroke A numeric value specifying the stroke width for the background lines.
#' @param whiskers A logical value indicating whether to include whiskers in the plot.
#' @param whiskers_linewidth A numeric value specifying the linewidth for the whiskers.
#' @param whiskers_alpha A numeric value specifying the transparency for the whiskers.
#' @param aspect.ratio A numeric value specifying the aspect ratio of the plot.
#' @param title An optional character string specifying the plot title.
#' @param subtitle An optional character string specifying the plot subtitle.
#' @param xlab An optional character string specifying the x-axis label.
#' @param ylab An optional character string specifying the y-axis label.
#' @param legend.position A character string specifying the position of the legend.
#' @param legend.direction A character string specifying the direction of the legend.
#' @param theme_use A character string specifying the theme to use for the plot.
#' @param theme_args A list of additional arguments to pass to the theme function.
#' @param return_layer A logical value indicating whether to return the plot as a layer.
#' @param seed An optional integer specifying the random seed for reproducibility.
#'
#' @importFrom Seurat Key Embeddings
#' @importFrom ggplot2 aes geom_path geom_segment labs
#' @importFrom grid arrow unit
#' @importFrom stats loess quantile
#' @export
#'
#' @seealso \code{\link{RunSlingshot}} \code{\link{CellDimPlot}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunSlingshot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   show_plot = FALSE
#' )
#' LineagePlot(
#'   pancreas_sub,
#'   lineages = paste0("Lineage", 1:3)
#' )
#' LineagePlot(
#'   pancreas_sub,
#'   lineages = paste0("Lineage", 1:3),
#'   whiskers = TRUE
#' )
LineagePlot <- function(
    srt,
    lineages,
    reduction = NULL,
    dims = c(1, 2),
    cells = NULL,
    trim = c(0.01, 0.99),
    span = 0.75,
    palette = "Dark2",
    palcolor = NULL,
    lineages_arrow = arrow(length = unit(0.1, "inches")),
    linewidth = 1,
    line_bg = "white",
    line_bg_stroke = 0.5,
    whiskers = FALSE,
    whiskers_linewidth = 0.5,
    whiskers_alpha = 0.5,
    aspect.ratio = 1,
    title = NULL,
    subtitle = NULL,
    xlab = NULL,
    ylab = NULL,
    legend.position = "right",
    legend.direction = "vertical",
    theme_use = "theme_scop",
    theme_args = list(),
    return_layer = FALSE,
    seed = 11) {
  set.seed(seed)

  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  if (!reduction %in% names(srt@reductions)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }

  reduction_key <- srt@reductions[[reduction]]@key
  dat_dim <- srt@reductions[[reduction]]@cell.embeddings
  colnames(dat_dim) <- paste0(reduction_key, seq_len(ncol(dat_dim)))
  rownames(dat_dim) <- rownames(dat_dim) %||% colnames(srt@assays[[1]])
  dat_lineages <- srt@meta.data[, unique(lineages), drop = FALSE]
  dat <- cbind(dat_dim, dat_lineages[row.names(dat_dim), , drop = FALSE])
  dat[, "cell"] <- rownames(dat)
  if (!is.null(cells)) {
    dat <- dat[intersect(rownames(dat), cells), , drop = FALSE]
  }

  xlab <- xlab %||% paste0(reduction_key, dims[1])
  ylab <- ylab %||% paste0(reduction_key, dims[2])
  if (identical(theme_use, "theme_blank")) {
    theme_args[["xlab"]] <- xlab
    theme_args[["ylab"]] <- ylab
  }

  colors <- palette_scop(lineages, palette = palette, palcolor = palcolor)
  axes <- paste0(reduction_key, dims)
  fitted_list <- lapply(lineages, function(l) {
    trim_pass <- dat[[l]] > quantile(dat[[l]], trim[1], na.rm = TRUE) &
      dat[[l]] < quantile(dat[[l]], trim[2], na.rm = TRUE)
    na_pass <- !is.na(dat[[l]])
    index <- which(trim_pass & na_pass)
    index <- index[order(dat[index, l])]
    dat_sub <- dat[index, , drop = FALSE]
    # if (is.null(weights)) {
    weights_used <- rep(1, nrow(dat_sub))
    # } else {
    # weights_used <- dat_sub[[weights]]
    # }
    fitted <- lapply(axes, function(x) {
      loess(
        formula(paste(x, l, sep = "~")),
        weights = weights_used,
        data = dat_sub,
        span = span,
        degree = 2
      )$fitted
    })
    names(fitted) <- axes
    fitted[["index"]] <- index
    return(fitted)
  })
  names(fitted_list) <- lineages

  curve_layer <- lapply(lineages, function(l) {
    dat_smooth <- as.data.frame(fitted_list[[l]])
    colnames(dat_smooth) <- c(
      paste0("Axis_", 1:(ncol(dat_smooth) - 1)),
      "index"
    )
    dat_smooth[, "Lineages"] <- factor(l, levels = lineages)
    dat_smooth <- unique(na.omit(dat_smooth))
    curve <- list()
    if (isTRUE(whiskers)) {
      dat_smooth[, "raw_Axis_1"] <- dat[dat_smooth[, "index"], axes[1]]
      dat_smooth[, "raw_Axis_2"] <- dat[dat_smooth[, "index"], axes[2]]
      curve <- c(
        curve,
        geom_segment(
          data = dat_smooth,
          mapping = aes(
            x = Axis_1,
            y = Axis_2,
            xend = raw_Axis_1,
            yend = raw_Axis_2,
            color = Lineages
          ),
          linewidth = whiskers_linewidth,
          alpha = whiskers_alpha,
          show.legend = TRUE,
          inherit.aes = FALSE
        )
      )
    }
    curve <- c(
      curve,
      geom_path(
        data = dat_smooth,
        mapping = aes(x = Axis_1, y = Axis_2),
        color = line_bg,
        linewidth = linewidth + line_bg_stroke,
        arrow = lineages_arrow,
        show.legend = TRUE,
        inherit.aes = FALSE
      ),
      geom_path(
        data = dat_smooth,
        mapping = aes(x = Axis_1, y = Axis_2, color = Lineages),
        linewidth = linewidth,
        arrow = lineages_arrow,
        show.legend = TRUE,
        inherit.aes = FALSE
      )
    )
    return(curve)
  })
  curve_layer <- c(
    unlist(curve_layer),
    list(scale_color_manual(values = colors))
  )

  lab_layer <- list(labs(
    title = title,
    subtitle = subtitle,
    x = xlab,
    y = ylab
  ))
  theme_layer <- list(
    do.call(theme_use, theme_args) +
      theme(
        aspect.ratio = aspect.ratio,
        legend.position = legend.position,
        legend.direction = legend.direction
      )
  )

  if (isTRUE(return_layer)) {
    return(list(
      curve_layer = curve_layer,
      lab_layer = lab_layer,
      theme_layer = theme_layer
    ))
  } else {
    return(
      ggplot() +
        curve_layer +
        lab_layer +
        theme_layer
    )
  }
}
