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
#' @param alpha Overall point and box alpha. `point_alpha` and `box_alpha`
#' take precedence when set.
#' @param point_size,point_alpha Jittered point size and alpha.
#' @param box_alpha Boxplot fill alpha.
#' @param flip Whether to flip coordinates so group labels run along the y
#' axis. Default is `TRUE`.
#' @param title,subtitle,xlab,ylab Plot labels.
#' @param theme_use Theme function name. Default is `"theme_scop"`, which
#' maps to `thisplot::theme_this()`, matching the style used by
#' `thisplot::StatPlot()`.
#' @param theme_args Additional arguments passed to `theme_use`.
#' @param legend.position,legend.direction Legend position and direction.
#' Default is `"none"`.
#' @param grid_major Whether to show major panel grid lines.
#' @param grid_major_colour,grid_major_linetype,grid_major_linewidth
#' Appearance of major panel grid lines.
#' @param aspect.ratio Aspect ratio of the plot (`y / x`).
#' @param seed RNG seed for jitter reproducibility.
#' @param ... Reserved for future use.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @seealso [RuntAge]
#'
#' @examples
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
  point_size = 2.2,
  point_alpha = 0.85,
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

  group_levels <- unique(as.character(plot_df[[group.by]]))
  colors <- thisplot::palette_colors(
    group_levels,
    palette = palette,
    palcolor = palcolor,
    type = "discrete"
  )
  plot_df[[group.by]] <- factor(
    plot_df[[group.by]],
    levels = names(sort(tapply(plot_df[[score_col]], plot_df[[group.by]], median, na.rm = TRUE)))
  )

  if (identical(theme_use, "theme_scop")) {
    theme_use <- "theme_this"
  }
  theme_fun <- get_namespace_fun("thisplot", theme_use)

  point_alpha <- point_alpha %||% alpha
  box_alpha  <- box_alpha  %||% alpha

  plot <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data[[group.by]],
      y = .data[[score_col]],
      color = .data[[group.by]]
    )
  ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = .data[[group.by]]),
      outlier.shape = NA,
      width = 0.55,
      alpha = box_alpha
    ) +
    ggplot2::geom_jitter(
      ggplot2::position_jitter(width = 0.16, height = 0, seed = seed),
      size = point_size,
      alpha = point_alpha
    ) +
    ggplot2::scale_color_manual(
      values = colors,
      guide = if (identical(legend.position, "none")) "none" else ggplot2::guide_legend()
    ) +
    ggplot2::scale_fill_manual(
      values = colors,
      guide = if (identical(legend.position, "none")) "none" else ggplot2::guide_legend()
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab %||% score_col,
      color = group.by,
      fill = group.by
    ) +
    do.call(theme_fun, theme_args) +
    ggplot2::theme(
      legend.position = legend.position,
      legend.direction = legend.direction
    )

  if (isTRUE(flip)) {
    plot <- plot + ggplot2::coord_flip()
  }

  if (!is.null(aspect.ratio)) {
    plot <- plot + ggplot2::theme(aspect.ratio = aspect.ratio)
  }

  if (!is.null(split.by)) {
    if (!split.by %in% colnames(plot_df)) {
      log_message(
        "{.arg split.by} must identify a column in the tAge prediction table",
        message_type = "error"
      )
    }
    plot <- plot + ggplot2::facet_wrap(stats::as.formula(paste("~", split.by)))
  }

  add_major_grid_theme(
    plot = plot,
    grid_major = grid_major,
    grid_major_colour = grid_major_colour,
    grid_major_linetype = grid_major_linetype,
    grid_major_linewidth = grid_major_linewidth
  )
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
