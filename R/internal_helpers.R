scop_annotate_quadrants <- function(
  plot,
  x,
  y,
  cutoffs,
  line_color = "grey30",
  line_type = "solid",
  line_width = 0.5,
  label_size = 3,
  group = NULL
) {
  cutoffs_x <- NULL
  cutoffs_y <- NULL
  if (is.list(cutoffs)) {
    cutoffs_x <- cutoffs[[1]]
    cutoffs_y <- cutoffs[[2]]
  } else if (is.numeric(cutoffs)) {
    if (length(cutoffs) == 1) {
      cutoffs_x <- cutoffs
      cutoffs_y <- cutoffs
    } else {
      cutoffs_x <- cutoffs[1]
      cutoffs_y <- cutoffs[2]
    }
  }

  if (!is.null(cutoffs_x)) {
    plot <- plot +
      ggplot2::geom_vline(
        xintercept = cutoffs_x,
        linetype = line_type,
        color = line_color,
        linewidth = line_width
      )
  }
  if (!is.null(cutoffs_y)) {
    plot <- plot +
      ggplot2::geom_hline(
        yintercept = cutoffs_y,
        linetype = line_type,
        color = line_color,
        linewidth = line_width
      )
  }

  plot_data <- plot$data
  plot_limits <- ggplot2::ggplot_build(plot)$layout$panel_params[[1]]

  x_breaks <- c(plot_limits$x.range[1], cutoffs_x, plot_limits$x.range[2])
  x_breaks <- unique(sort(x_breaks))
  y_breaks <- c(plot_limits$y.range[1], cutoffs_y, plot_limits$y.range[2])
  y_breaks <- unique(sort(y_breaks))

  x_breaks <- x_breaks[
    x_breaks >= plot_limits$x.range[1] & x_breaks <= plot_limits$x.range[2]
  ]
  y_breaks <- y_breaks[
    y_breaks >= plot_limits$y.range[1] & y_breaks <= plot_limits$y.range[2]
  ]

  x_pos <- (utils::head(x_breaks, -1) + utils::tail(x_breaks, -1)) / 2
  y_pos <- (utils::head(y_breaks, -1) + utils::tail(y_breaks, -1)) / 2

  plot_data$x_cat <- as.integer(cut(
    plot_data[[x]],
    breaks = x_breaks,
    include.lowest = TRUE
  ))
  plot_data$y_cat <- as.integer(cut(
    plot_data[[y]],
    breaks = y_breaks,
    include.lowest = TRUE
  ))

  if (is.null(group)) {
    counts <- as.data.frame(
      table(plot_data[, c("x_cat", "y_cat")]),
      stringsAsFactors = FALSE
    )
    counts$value <- counts$Freq / sum(counts$Freq) * 100
  } else {
    counts <- as.data.frame(
      table(plot_data[, c(group, "x_cat", "y_cat")]),
      stringsAsFactors = FALSE
    )
    counts_list <- split(counts, counts[[group]])
    counts <- do.call(
      rbind,
      lapply(counts_list, function(df) {
        df$value <- df$Freq / sum(df$Freq) * 100
        df
      })
    )
  }

  counts$x_cat <- as.integer(counts$x_cat)
  counts$y_cat <- as.integer(counts$y_cat)

  annot_df <- counts[counts$Freq > 0, , drop = FALSE]
  if (nrow(annot_df) == 0) {
    return(plot)
  }

  annot_df$x <- x_pos[annot_df$x_cat]
  annot_df$y <- y_pos[annot_df$y_cat]
  annot_df$label <- paste0(round(annot_df$value, 1), "%")

  plot +
    ggplot2::geom_text(
      data = annot_df,
      ggplot2::aes(x = x, y = y, label = label),
      size = label_size
    )
}

scop_clip_symmetric_range <- function(data, value_col = "avg_log2FC") {
  values <- data[[value_col]][is.finite(data[[value_col]])]
  if (length(values) == 0) {
    return(list(data = data, limits = c(-1, 1)))
  }

  upper <- stats::quantile(values, c(0.99, 1))
  lower <- stats::quantile(values, c(0.01, 0))
  upper <- ifelse(upper[1] > 0, upper[1], upper[2])
  lower <- ifelse(lower[1] < 0, lower[1], lower[2])

  if (upper > 0 && lower < 0) {
    value_range <- min(abs(c(upper, lower)), na.rm = TRUE)
    upper <- value_range
    lower <- -value_range
  }

  data[data[[value_col]] > upper, value_col] <- upper
  data[data[[value_col]] < lower, value_col] <- lower

  list(data = data, limits = c(lower, upper))
}

scop_jitter_highlighted_points <- function(
  data,
  jitter_width = 0.2,
  jitter_height = 0.2,
  seed = 11
) {
  data[, "x_plot"] <- data[, "x"]
  data[, "y_plot"] <- data[, "y"]

  border_idx <- which(
    data[, "border"] & is.finite(data[, "x"]) & is.finite(data[, "y"])
  )
  if (length(border_idx) == 0) {
    return(data)
  }

  idx <- seq_along(border_idx)
  x_offset <- ((((idx * 0.61803398875) + (seed * 0.01)) %% 1) - 0.5) *
    2 *
    jitter_width
  y_offset <- ((((idx * 0.41421356237) + (seed * 0.01)) %% 1) - 0.5) *
    2 *
    jitter_height

  data[border_idx, "x_plot"] <- data[border_idx, "x"] + x_offset
  data[border_idx, "y_plot"] <- data[border_idx, "y"] + y_offset
  data
}

scop_collapse_sparse_rows <- function(matrix, group) {
  if (length(group) != nrow(matrix)) {
    cli::cli_abort(
      "{.arg group} length must match the number of rows of {.arg matrix}"
    )
  }

  group <- as.character(group)
  keep <- !is.na(group) & nzchar(group)
  matrix <- matrix[keep, , drop = FALSE]
  group <- group[keep]

  levels_use <- unique(group)
  matrix_summary <- Matrix::summary(matrix)
  i_new <- match(group[matrix_summary$i], levels_use)

  Matrix::sparseMatrix(
    i = i_new,
    j = matrix_summary$j,
    x = matrix_summary$x,
    dims = c(length(levels_use), ncol(matrix)),
    dimnames = list(levels_use, colnames(matrix))
  )
}
