#' Plot tAge transcriptomic aging-clock predictions
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param object A `Seurat` object returned by [RuntAge()] or a data frame with
#' tAge prediction columns.
#' @param score_col Prediction column to plot. If `NULL`, the first column
#' ending in `"_tAge"` is used.
#' @param group.by Column used to group predictions on the y axis.
#' @param split.by An optional column to split the plot into facets.
#' @param tool_name Name of the Seurat tool entry containing tAge results.
#' @param plot_type Plot type. Currently `"box"` draws a boxplot with jittered
#' pseudobulk samples.
#' @param palette,palcolor Palette forwarded to `thisplot::palette_colors()`.
#' @param alpha Overall point and box alpha. `pt.alpha` and `box_alpha`
#' take precedence when set.
#' @param point_size,point_alpha Deprecated alias(es) for `pt.size`/`pt.alpha`; supply exactly one of the two. It
#' will be removed in scop 1.0.0.
#' @param pt.size,pt.alpha Jittered point size and alpha.
#' @param box_alpha Boxplot fill alpha.
#' @param flip Whether to flip coordinates so group labels run along the y
#' axis.
#' @param title,subtitle,xlab,ylab Plot labels.
#' @param theme_use Theme function name passed to the shared
#' expression-statistics renderer. Default is `"theme_scop"`.
#' @param theme_args Additional arguments passed to `theme_use`.
#' @param legend.position,legend.direction Legend position and direction.
#' @param grid_major Whether to show major panel grid lines.
#' @param grid_major_colour,grid_major_linetype,grid_major_linewidth
#' Appearance of major panel grid lines.
#' @param aspect.ratio Aspect ratio of the plot (`y / x`).
#' @param seed RNG seed for jitter reproducibility.
#' @param ... Arguments passed to other functions.
#' @return A `ggplot` object.
#'
#' @seealso [RuntAge], [FeatureStatPlot]
#'
#' @export
#'
#' @examples
#' predictions <- data.frame(
#'   group = rep(c("Young", "Aged"), each = 4),
#'   Mouse_tAge = c(3.1, 4.2, 5.0, 4.6, 18.2, 20.1, 21.3, 19.6)
#' )
#' tAgePlot(predictions, group.by = "group")
#'
#' \dontrun{
#' pancreas_sub <- RuntAge(pancreas_sub, group.by = "CellType")
#' tAgePlot(pancreas_sub, group.by = "CellType")
#' }
tAgePlot <- function(
  object,
  score_col = NULL,
  group.by = NULL,
  split.by = NULL,
  tool_name = "tAge",
  plot_type = c("box"),
  palette = "Chinese",
  palcolor = NULL,
  alpha = 1,
  pt.size = 2.2,
  point_size = NULL,
  pt.alpha = 0.85,
  point_alpha = NULL,
  box_alpha = 0.15,
  flip = TRUE,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  legend.position = "none",
  legend.direction = "vertical",
  grid_major = TRUE,
  grid_major_colour = "grey80",
  grid_major_linetype = 2,
  grid_major_linewidth = 0.3,
  aspect.ratio = NULL,
  seed = 11,
  verbose = TRUE,
  ...
) {
  plot_type <- match.arg(plot_type)
  if (!is.null(point_alpha)) {
    .Deprecated(msg = paste0(
      "`point_alpha` is deprecated; use `pt.alpha` instead. ",
      "It will be removed in scop 1.0.0."
    ))
    pt.alpha <- point_alpha
  }
  point_alpha <- pt.alpha
  if (!is.null(point_size)) {
    .Deprecated(msg = paste0(
      "`point_size` is deprecated; use `pt.size` instead. ",
      "It will be removed in scop 1.0.0."
    ))
    pt.size <- point_size
  }
  point_size <- pt.size
  plot_df <- tage_plot_data(
    object = object,
    tool_name = tool_name,
    verbose = verbose
  )
  score_cols <- grep("_tAge$", colnames(plot_df), value = TRUE)
  score_col <- score_col %||% if (length(score_cols) > 0L) score_cols[[1]] else NULL
  if (is.null(score_col) || !score_col %in% colnames(plot_df)) {
    log_message(
      "{.arg score_col} must identify a tAge prediction column",
      message_type = "error"
    )
  }
  group.by <- group.by %||% infer_tage_plot_group(plot_df)
  if (is.null(group.by) || !group.by %in% colnames(plot_df)) {
    log_message(
      "{.arg group.by} must identify a column in the tAge prediction table",
      message_type = "error"
    )
  }

  plot_df[[group.by]] <- factor(
    plot_df[[group.by]],
    levels = names(sort(tapply(plot_df[[score_col]], plot_df[[group.by]], stats::median, na.rm = TRUE)))
  )
  plot_df$cells <- rownames(plot_df)
  rownames(plot_df) <- paste0("tAge_", plot_df$cells)
  plot <- ExpressionStatPlot(
    meta.data = plot_df, stat.by = score_col, group.by = group.by,
    plot_type = "box", palette = palette, palcolor = palcolor,
    add_point = TRUE, pt.size = point_size,
    pt.alpha = point_alpha %||% alpha, jitter.width = 0.16, jitter.height = 0,
    box_width = 0.55, alpha = box_alpha %||% alpha,
    title = title, subtitle = subtitle, xlab = xlab %||% "",
    ylab = ylab %||% score_col, legend.title = group.by,
    theme_use = theme_use, theme_args = theme_args,
    legend.position = legend.position, legend.direction = legend.direction,
    grid_major = grid_major, grid_major_colour = grid_major_colour,
    grid_major_linetype = grid_major_linetype,
    grid_major_linewidth = grid_major_linewidth,
    aspect.ratio = aspect.ratio, seed = seed
  )[[1L]] + ggplot2::facet_null()
  background_layers <- vapply(plot$layers, function(layer) inherits(layer$geom, "GeomRect"), logical(1))
  plot$layers <- plot$layers[!background_layers]
  for (index in seq_along(plot$layers)) {
    if (inherits(plot$layers[[index]]$geom, "GeomBoxplot")) {
      plot$layers[[index]]$aes_params$alpha <- box_alpha %||% alpha
    }
    if (inherits(plot$layers[[index]]$position, "PositionJitterdodge")) {
      plot$layers[[index]]$position$seed <- seed
    }
  }
  if (!is.null(split.by)) {
    if (!split.by %in% colnames(plot_df)) {
      log_message(
        "{.arg split.by} must identify a column in the tAge prediction table",
        message_type = "error"
      )
    }
    plot$data[[split.by]] <- plot_df[[split.by]][match(plot$data$cells, plot_df$cells)]
    plot <- plot + ggplot2::facet_wrap(stats::reformulate(split.by))
  }
  if (isTRUE(flip)) {
    plot <- plot + ggplot2::coord_flip()
  }
  plot
}

tage_plot_data <- function(object, tool_name = "tAge", verbose = TRUE) {
  if (inherits(object, "Seurat")) {
    result <- object@tools[[tool_name]]
    if (is.null(result) || is.null(result$predictions)) {
      log_message(
        "{.arg object} does not contain tAge predictions in {.code object@tools[[{tool_name}]]}",
        message_type = "error"
      )
    }
    plot_df <- as.data.frame(result$predictions, check.names = FALSE)
    plot_df$sample_id <- rownames(plot_df)
    meta <- result$pseudobulk_metadata
    if (!is.null(meta)) {
      meta <- as.data.frame(meta, check.names = FALSE)
      meta$sample_id <- rownames(meta)
      add_cols <- setdiff(colnames(meta), colnames(plot_df))
      plot_df <- merge(
        plot_df,
        meta[, c("sample_id", add_cols), drop = FALSE],
        by = "sample_id",
        all.x = TRUE,
        sort = FALSE
      )
    }
    return(plot_df)
  }
  if (!is.data.frame(object)) {
    log_message(
      "{.arg object} must be a {.cls Seurat} object or data frame",
      message_type = "error"
    )
  }
  object
}

infer_tage_plot_group <- function(plot_df) {
  candidates <- setdiff(
    colnames(plot_df),
    c(
      "sample_id",
      grep("_tAge$", colnames(plot_df), value = TRUE),
      grep("_tAge_std$", colnames(plot_df), value = TRUE),
      "cumulative_coverage",
      "n_cells"
    )
  )
  candidates <- candidates[
    vapply(plot_df[candidates], function(x) is.character(x) || is.factor(x), logical(1))
  ]
  candidates[[1]] %||% NULL
}
