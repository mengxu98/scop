#' @title Plot deconvolution results
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param object A `SummarizedExperiment` object containing deconvolution
#' results in `metadata(object)[["Deconvolution"]]`.
#' @param res A deconvolution result data frame. When provided, `object` is
#' ignored.
#' @param plot_type Plot type. One of `"bar"`, `"heatmap"`, or `"box"`.
#' @param sample_order Optional sample order.
#' @param cell_type_order Optional cell-type order.
#' @param palette Palette used for discrete cell-type colors in `"bar"` and
#' `"box"` modes.
#' @param palcolor Optional custom colors for `palette`.
#' @param heatmap_palette Palette used for continuous proportions in
#' `"heatmap"` mode.
#' @param heatmap_palcolor Optional custom colors for `heatmap_palette`.
#' @param sample_annotation Character vector of `colData(object)` columns used
#' as top annotations in `"heatmap"` mode.
#' @param sample_annotation_palette Palette(s) used for `sample_annotation`.
#' @param sample_annotation_palcolor Optional custom colors for
#' `sample_annotation_palette`. Use a list when multiple annotations are
#' provided.
#' @param sample_split Optional `colData(object)` column used to split heatmap
#' columns for grouped comparison.
#' @param cluster_rows,cluster_columns Whether to cluster rows or columns in
#' `"heatmap"` mode.
#' @param show_row_names,show_column_names Whether to show row or column names in
#' `"heatmap"` mode.
#' @param theme_use Theme function name. Default is `"theme_scop"`.
#' @param theme_args Additional theme arguments passed to `theme_use`.
#' @param legend.position Legend position. Default is `"right"`.
#' @param legend.direction Legend direction. Default is `"vertical"`.
#' @param grid_major Whether to show major grid lines for `"bar"` and `"box"`
#' plots. Default is `FALSE`.
#' @param ... Reserved for future use.
#'
#' @return A `ggplot` object.
#' For `plot_type = "heatmap"`, returns a `ComplexHeatmap::Heatmap` object.
#'
#' @seealso [RunDeconvolution]
#'
#' @export
#'
#' @examples
#' data(islet_bulk)
#' data(panc8_sub)
#' islet_bulk <- RunDeconvolution(
#'   islet_bulk,
#'   reference = panc8_sub,
#'   method = "MuSiC",
#'   group.by = "celltype"
#' )
#' DeconvolutionPlot(islet_bulk, plot_type = "bar")
#'
#' DeconvolutionPlot(islet_bulk, plot_type = "box")
#'
#' ht <- DeconvolutionPlot(
#'   islet_bulk,
#'   plot_type = "heatmap",
#'   sample_annotation = "condition",
#'   sample_split = "condition"
#' )
#' ComplexHeatmap::draw(ht)
DeconvolutionPlot <- function(
  object = NULL,
  res = NULL,
  plot_type = c("bar", "heatmap", "box"),
  sample_order = NULL,
  cell_type_order = NULL,
  palette = "Chinese",
  palcolor = NULL,
  heatmap_palette = "YlGnBu",
  heatmap_palcolor = NULL,
  sample_annotation = NULL,
  sample_annotation_palette = "Chinese",
  sample_annotation_palcolor = NULL,
  sample_split = NULL,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  theme_use = "theme_scop",
  theme_args = list(),
  legend.position = "right",
  legend.direction = "vertical",
  grid_major = FALSE,
  verbose = TRUE,
  ...
) {
  plot_type <- match.arg(plot_type)
  df <- resolve_deconvolution_result(object = object, res = res)
  if (nrow(df) == 0) {
    log_message(
      "No deconvolution result is available for plotting.",
      message_type = "error"
    )
  }

  df$sample <- bulk_match_levels(df$sample, sample_order)
  df$cell_type <- bulk_match_levels(df$cell_type, cell_type_order)
  fill_cols <- palette_colors(
    levels(df$cell_type),
    palette = palette,
    palcolor = palcolor
  )
  theme_obj <- bulk_plot_theme(theme_use = theme_use, theme_args = theme_args)

  if (identical(plot_type, "bar")) {
    plot <- ggplot2::ggplot(df, ggplot2::aes(
      x = sample,
      y = proportion,
      fill = cell_type
    )) +
      ggplot2::geom_col(
        width = 0.78,
        color = "black",
        linewidth = 0.2
      ) +
      ggplot2::scale_fill_manual(
        values = fill_cols,
        drop = FALSE
      ) +
      ggplot2::labs(
        x = NULL,
        y = "Proportion",
        fill = "Cell type"
      ) +
      theme_obj +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1,
          vjust = 1
        ),
        legend.position = legend.position,
        legend.direction = legend.direction,
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )
    return(add_major_grid_theme(plot = plot, grid_major = grid_major))
  }

  if (identical(plot_type, "heatmap")) {
    return(
      deconvolution_heatmap(
        object = object,
        df = df,
        sample_annotation = sample_annotation,
        sample_annotation_palette = sample_annotation_palette,
        sample_annotation_palcolor = sample_annotation_palcolor,
        sample_split = sample_split,
        heatmap_palette = heatmap_palette,
        heatmap_palcolor = heatmap_palcolor,
        cluster_rows = cluster_rows,
        cluster_columns = cluster_columns,
        show_row_names = show_row_names,
        show_column_names = show_column_names
      )
    )
  }

  plot <- ggplot2::ggplot(df, ggplot2::aes(
    x = cell_type,
    y = proportion,
    fill = cell_type
  )) +
    ggplot2::geom_boxplot(
      outlier.shape = NA,
      width = 0.68,
      linewidth = 0.35,
      alpha = 0.92
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = cell_type),
      position = ggplot2::position_jitter(width = 0.14, height = 0),
      size = 1.3,
      alpha = 0.75,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = fill_cols,
      drop = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = fill_cols,
      drop = FALSE
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Proportion",
      fill = "Cell type"
    ) +
    theme_obj +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  add_major_grid_theme(plot = plot, grid_major = grid_major)
}

bulk_plot_theme <- function(theme_use = "theme_scop", theme_args = list()) {
  if (identical(theme_use, "theme_scop")) {
    theme_use <- "theme_this"
  }
  theme_fun <- tryCatch(get(theme_use, mode = "function"), error = function(e) NULL)
  if (is.null(theme_fun)) {
    return(ggplot2::theme_bw())
  }
  do.call(theme_fun, theme_args)
}

bulk_match_levels <- function(x, levels_use = NULL) {
  x <- as.character(x)
  if (is.null(levels_use)) {
    return(factor(x, levels = unique(x)))
  }
  levels_use <- unique(as.character(levels_use))
  factor(x, levels = levels_use[levels_use %in% x])
}

deconvolution_heatmap <- function(
  object,
  df,
  sample_annotation = NULL,
  sample_annotation_palette = "Chinese",
  sample_annotation_palcolor = NULL,
  sample_split = NULL,
  heatmap_palette = "YlGnBu",
  heatmap_palcolor = NULL,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE
) {
  mat <- as.matrix(stats::xtabs(proportion ~ cell_type + sample, data = df))
  mat <- mat[
    levels(df$cell_type)[levels(df$cell_type) %in% rownames(mat)],
    levels(df$sample)[levels(df$sample) %in% colnames(mat)],
    drop = FALSE
  ]

  fill_limits <- range(mat, na.rm = TRUE)
  if (!all(is.finite(fill_limits))) {
    fill_limits <- c(0, 1)
  }
  if (diff(fill_limits) <= 0) {
    fill_limits <- c(fill_limits[1], fill_limits[1] + 1e-8)
  }
  col_fun <- circlize::colorRamp2(
    breaks = seq(fill_limits[1], fill_limits[2], length.out = 100),
    colors = palette_colors(
      palette = heatmap_palette,
      palcolor = heatmap_palcolor
    )
  )

  column_meta <- resolve_deconvolution_column_metadata(
    object = object,
    sample_names = colnames(mat),
    sample_annotation = sample_annotation,
    sample_annotation_palette = sample_annotation_palette,
    sample_annotation_palcolor = sample_annotation_palcolor,
    sample_split = sample_split
  )

  ComplexHeatmap::Heatmap(
    matrix = mat,
    col = col_fun,
    name = "Proportion",
    top_annotation = column_meta$top_annotation,
    column_split = column_meta$column_split,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    row_names_side = "left",
    column_names_side = "bottom",
    column_names_rot = 45,
    rect_gp = grid::gpar(col = "white", lwd = 0.8),
    border = TRUE,
    heatmap_legend_param = list(
      title = "Proportion",
      direction = "vertical",
      border = "black"
    )
  )
}

resolve_deconvolution_column_metadata <- function(
  object,
  sample_names,
  sample_annotation = NULL,
  sample_annotation_palette = "Chinese",
  sample_annotation_palcolor = NULL,
  sample_split = NULL
) {
  if (is.null(object) || !inherits(object, "SummarizedExperiment")) {
    if (!is.null(sample_annotation) || !is.null(sample_split)) {
      log_message(
        "{.arg sample_annotation} and {.arg sample_split} require {.arg object} to be a {.cls SummarizedExperiment}.",
        message_type = "warning"
      )
    }
    return(list(top_annotation = NULL, column_split = NULL))
  }

  meta <- as.data.frame(SummarizedExperiment::colData(object))
  meta$sample <- rownames(meta)
  idx <- match(sample_names, meta$sample)
  if (any(is.na(idx))) {
    log_message(
      "Some deconvolution samples were not found in {.fn colData(object)}.",
      message_type = "warning"
    )
  }
  meta <- meta[idx, , drop = FALSE]
  rownames(meta) <- sample_names

  top_annotation <- NULL
  if (!is.null(sample_annotation)) {
    sample_annotation <- intersect(sample_annotation, colnames(meta))
    if (length(sample_annotation) == 0) {
      log_message(
        "No requested {.arg sample_annotation} columns were found in {.fn colData(object)}.",
        message_type = "warning"
      )
    } else {
      palette_use <- rep_len(sample_annotation_palette, length(sample_annotation))
      if (is.null(sample_annotation_palcolor)) {
        palcolor_use <- rep(list(NULL), length(sample_annotation))
      } else if (is.list(sample_annotation_palcolor)) {
        palcolor_use <- rep_len(sample_annotation_palcolor, length(sample_annotation))
      } else {
        palcolor_use <- rep(list(sample_annotation_palcolor), length(sample_annotation))
      }

      anno_df <- meta[, sample_annotation, drop = FALSE]
      col_list <- list()
      for (i in seq_along(sample_annotation)) {
        nm <- sample_annotation[i]
        values <- anno_df[[nm]]
        palette_i <- palette_use[i]
        palcolor_i <- palcolor_use[[i]]
        if (is.logical(values)) {
          values <- factor(values, levels = c(TRUE, FALSE))
          anno_df[[nm]] <- values
        } else if (!is.numeric(values)) {
          values <- factor(values, levels = unique(values))
          anno_df[[nm]] <- values
        }

        if (is.numeric(values)) {
          range_i <- range(values, na.rm = TRUE)
          if (!all(is.finite(range_i))) {
            next
          }
          if (diff(range_i) <= 0) {
            range_i <- c(range_i[1], range_i[1] + 1e-8)
          }
          col_list[[nm]] <- circlize::colorRamp2(
            breaks = seq(range_i[1], range_i[2], length.out = 100),
            colors = palette_colors(
              palette = palette_i,
              palcolor = palcolor_i
            )
          )
        } else {
          col_list[[nm]] <- stats::setNames(
            palette_colors(
              levels(anno_df[[nm]]),
              palette = palette_i,
              palcolor = palcolor_i
            ),
            levels(anno_df[[nm]])
          )
        }
      }

      top_annotation <- ComplexHeatmap::HeatmapAnnotation(
        df = anno_df,
        col = col_list,
        which = "column",
        border = TRUE,
        show_annotation_name = TRUE,
        annotation_name_side = "left",
        simple_anno_size = grid::unit(4, "mm")
      )
    }
  }

  column_split <- NULL
  if (!is.null(sample_split)) {
    if (!sample_split %in% colnames(meta)) {
      log_message(
        "{.arg sample_split} was not found in {.fn colData(object)}.",
        message_type = "warning"
      )
    } else {
      split_values <- meta[[sample_split]]
      if (!is.factor(split_values)) {
        split_values <- factor(split_values, levels = unique(split_values))
      }
      column_split <- split_values
    }
  }

  list(
    top_annotation = top_annotation,
    column_split = column_split
  )
}
