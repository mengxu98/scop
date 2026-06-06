#' @title Immune abundance plots
#'
#' @description
#' Visualize immune abundance or score matrices from deconvolution results,
#' spatial metadata, GSVA-like scores, or user-provided matrices.
#'
#' @md
#' @param object Optional `SummarizedExperiment`, `Seurat`, deconvolution bundle,
#' or matrix-like object.
#' @param immune.data Optional abundance matrix with samples/spots/cells in rows
#' and immune cell types or signatures in columns.
#' @param immune.cols Metadata columns to extract from a `Seurat` object.
#' If `NULL`, columns matching `RCTD_prop_*`, `*_prop_*`, or `*_frac_*` are used.
#' @param group.by Optional grouping column for `bar` and `box` plots.
#' @param group.data Optional named vector or data frame containing sample groups
#' when groups cannot be inferred from `object`.
#' @param plot_type Plot type: `"heatmap"`, `"bar"`, `"box"`, or `"cor"`.
#' @param sample_order Optional sample order.
#' @param cell_type_order Optional immune cell type or signature order.
#' @param scale Scaling for heatmap values: `"none"`, `"row"`, or `"column"`.
#' @param cor_method Correlation method for `plot_type = "cor"`.
#' @param bar_position Bar position for `plot_type = "bar"`.
#' `"stack"` shows raw abundance; `"fill"` rescales each sample to 1.
#' @param show_sample_names Whether to show sample names in heatmap and bar
#' plots.
#' @param show_cor_label Whether to print correlation values in `cor` plots.
#' @param cor_label_size Text size for correlation labels.
#' @param add_stat Whether to add group comparison labels to `box` plots.
#' Requires `ggpubr`.
#' @param comparisons Optional group comparisons passed to
#' `ggpubr::stat_compare_means()`.
#' @param pairwise_method Statistical test for box plot comparisons.
#' @param sig_label Significance label type for `ggpubr::stat_compare_means()`.
#' @param sig_labelsize Significance label size.
#' @param box_width,box_alpha Box plot width and alpha.
#' @param pt.size,pt.alpha,jitter.width Point size, alpha, and jitter width.
#' @param grid_major Whether to show major y grid lines for bar and box plots.
#' @param palette Discrete palette for grouped plots.
#' @param palcolor Optional custom discrete colors.
#' @param heatmap_palette Continuous palette name for heatmaps.
#' @param heatmap_palcolor Optional custom continuous colors.
#' @param title,subtitle Plot title and subtitle.
#' @param legend.position Legend position.
#' @param legend.direction Legend direction.
#' @param theme_use Theme function name.
#' @param theme_args Additional theme arguments.
#'
#' @return A `ggplot` object.
#'
#' @export
ImmuneAbundancePlot <- function(
  object = NULL,
  immune.data = NULL,
  immune.cols = NULL,
  group.by = NULL,
  group.data = NULL,
  plot_type = c("heatmap", "bar", "box", "cor"),
  sample_order = NULL,
  cell_type_order = NULL,
  scale = c("none", "row", "column"),
  cor_method = c("spearman", "pearson", "kendall"),
  bar_position = c("stack", "fill"),
  show_sample_names = FALSE,
  show_cor_label = FALSE,
  cor_label_size = 3,
  add_stat = FALSE,
  comparisons = NULL,
  pairwise_method = "wilcox.test",
  sig_label = c("p.signif", "p.format"),
  sig_labelsize = 3.5,
  box_width = 0.64,
  box_alpha = 0.92,
  pt.size = 1.35,
  pt.alpha = 0.72,
  jitter.width = 0.12,
  grid_major = FALSE,
  palette = "Chinese",
  palcolor = NULL,
  heatmap_palette = "YlGnBu",
  heatmap_palcolor = NULL,
  title = NULL,
  subtitle = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list()
) {
  plot_type <- match.arg(plot_type)
  scale <- match.arg(scale)
  cor_method <- match.arg(cor_method)
  bar_position <- match.arg(bar_position)
  sig_label <- match.arg(sig_label)
  mat <- resolve_immune_abundance(
    object = object,
    immune.data = immune.data,
    immune.cols = immune.cols
  )
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    log_message(
      "No immune abundance matrix is available for plotting.",
      message_type = "error"
    )
  }
  mat <- immune_order_matrix(
    mat = mat,
    sample_order = sample_order,
    cell_type_order = cell_type_order
  )

  theme_obj <- immune_plot_theme(theme_use = theme_use, theme_args = theme_args)
  if (identical(plot_type, "cor")) {
    cor_mat <- stats::cor(mat,
      method = cor_method,
      use = "pairwise.complete.obs"
    )
    cor_mat[!is.finite(cor_mat)] <- 0
    df <- immune_matrix_long(cor_mat, "cell_type_1", "cell_type_2", "cor")
    fill_cols <- immune_continuous_colors(
      palette = "RdBu",
      palcolor = heatmap_palcolor,
      reverse = TRUE
    )
    plot <- ggplot2::ggplot(df, ggplot2::aes(
        x = cell_type_2,
        y = cell_type_1,
        fill = cor
      )) +
        ggplot2::geom_tile(color = "white", linewidth = 0.5) +
        ggplot2::scale_fill_gradientn(
          colors = fill_cols,
          limits = c(-1, 1),
          na.value = "grey90",
          name = "Correlation"
        ) +
        ggplot2::coord_equal() +
        ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
        theme_obj +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid = ggplot2::element_blank(),
          legend.position = legend.position,
          legend.direction = legend.direction
        )
    if (isTRUE(show_cor_label)) {
      plot <- plot + ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.2f", cor)),
        size = cor_label_size,
        color = "grey15"
      )
    }
    return(plot)
  }

  if (identical(plot_type, "heatmap")) {
    mat_plot <- immune_scale_matrix(mat, scale = scale)
    df <- immune_matrix_long(mat_plot, "sample", "cell_type", "value")
    fill_cols <- immune_continuous_colors(
      palette = heatmap_palette,
      palcolor = heatmap_palcolor
    )
    return(
      ggplot2::ggplot(df, ggplot2::aes(
        x = cell_type,
        y = sample,
        fill = value
      )) +
        ggplot2::geom_tile(color = "white", linewidth = 0.35) +
        ggplot2::scale_fill_gradientn(
          colors = fill_cols,
          na.value = "grey90",
          name = "Abundance"
        ) +
        ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
        theme_obj +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.y = if (isTRUE(show_sample_names)) {
            ggplot2::element_text()
          } else {
            ggplot2::element_blank()
          },
          axis.ticks.y = if (isTRUE(show_sample_names)) {
            ggplot2::element_line()
          } else {
            ggplot2::element_blank()
          },
          panel.grid = ggplot2::element_blank(),
          legend.position = legend.position,
          legend.direction = legend.direction
        )
    )
  }

  df <- immune_matrix_long(mat, "sample", "cell_type", "abundance")
  group_df <- resolve_immune_groups(
    object = object,
    samples = rownames(mat),
    group.by = group.by,
    group.data = group.data
  )
  if (!is.null(group_df)) {
    df$group <- group_df$group[match(df$sample, group_df$sample)]
  } else {
    df$group <- "All"
  }
  df$group <- factor(df$group, levels = unique(df$group))

  fill_cols <- palette_colors(
    if (identical(plot_type, "bar")) unique(df$cell_type) else levels(df$group),
    palette = palette,
    palcolor = palcolor
  )

  if (identical(plot_type, "bar")) {
    ylab_use <- if (identical(bar_position, "fill")) {
      "Estimated proportion"
    } else {
      "Abundance"
    }
    return(
      ggplot2::ggplot(df, ggplot2::aes(
        x = sample,
        y = abundance,
        fill = cell_type
      )) +
        ggplot2::geom_col(
          width = 0.78,
          alpha = 0.9,
          color = "black",
          linewidth = 0.16,
          position = bar_position
        ) +
        ggplot2::facet_grid(. ~ group, scales = "free_x", space = "free_x") +
        ggplot2::scale_fill_manual(values = fill_cols, drop = FALSE) +
        ggplot2::labs(
          title = title,
          subtitle = subtitle,
          x = NULL,
          y = ylab_use,
          fill = "Cell type"
        ) +
        theme_obj +
        ggplot2::theme(
          axis.text.x = if (isTRUE(show_sample_names)) {
            ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
          } else {
            ggplot2::element_blank()
          },
          axis.ticks.x = if (isTRUE(show_sample_names)) {
            ggplot2::element_line()
          } else {
            ggplot2::element_blank()
          },
          legend.position = legend.position,
          legend.direction = legend.direction,
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.major.y = if (isTRUE(grid_major)) {
            ggplot2::element_line(colour = "grey86", linewidth = 0.25)
          } else {
            ggplot2::element_blank()
          },
          panel.grid.minor = ggplot2::element_blank(),
          strip.background = ggplot2::element_rect(
            fill = "grey94",
            colour = "grey55",
            linewidth = 0.35
          )
        )
    )
  }

  plot <- ggplot2::ggplot(df, ggplot2::aes(
    x = cell_type,
    y = abundance,
    fill = group
  )) +
    ggplot2::geom_boxplot(
      outlier.shape = NA,
      width = box_width,
      linewidth = 0.38,
      alpha = box_alpha,
      color = "grey20",
      position = ggplot2::position_dodge2(width = 0.78, preserve = "single")
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = group),
      position = ggplot2::position_jitterdodge(
        jitter.width = jitter.width,
        jitter.height = 0,
        dodge.width = 0.78
      ),
      size = pt.size,
      alpha = pt.alpha,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = fill_cols, drop = FALSE) +
    ggplot2::scale_color_manual(values = fill_cols, drop = FALSE) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = "Abundance",
      fill = group.by %||% "Group"
    ) +
    theme_obj +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = legend.position,
      legend.direction = legend.direction,
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = if (isTRUE(grid_major)) {
        ggplot2::element_line(colour = "grey86", linewidth = 0.25)
      } else {
        ggplot2::element_blank()
      },
      panel.grid.minor = ggplot2::element_blank()
    )
  if (isTRUE(add_stat)) {
    check_r("ggpubr", verbose = FALSE)
    stat_args <- list(
      mapping = ggplot2::aes(
        x = cell_type,
        y = abundance,
        group = group
      ),
      method = pairwise_method,
      label = sig_label,
      size = sig_labelsize,
      hide.ns = TRUE
    )
    if (!is.null(comparisons)) {
      stat_args$comparisons <- comparisons
    }
    plot <- plot + do.call(ggpubr::stat_compare_means, stat_args)
  }
  plot
}

immune_order_matrix <- function(
  mat,
  sample_order = NULL,
  cell_type_order = NULL
) {
  if (!is.null(sample_order)) {
    sample_order <- as.character(sample_order)
    sample_order <- c(
      intersect(sample_order, rownames(mat)),
      setdiff(rownames(mat), sample_order)
    )
    if (length(sample_order) == nrow(mat)) {
      mat <- mat[sample_order, , drop = FALSE]
    }
  }
  if (!is.null(cell_type_order)) {
    cell_type_order <- as.character(cell_type_order)
    cell_type_order <- c(
      intersect(cell_type_order, colnames(mat)),
      setdiff(colnames(mat), cell_type_order)
    )
    if (length(cell_type_order) == ncol(mat)) {
      mat <- mat[, cell_type_order, drop = FALSE]
    }
  }
  mat
}

resolve_immune_abundance <- function(
  object = NULL,
  immune.data = NULL,
  immune.cols = NULL
) {
  if (!is.null(immune.data)) {
    return(immune_numeric_matrix(immune.data, row_label = "immune.data"))
  }
  if (inherits(object, c("matrix", "data.frame", "Matrix"))) {
    return(immune_numeric_matrix(object, row_label = "object"))
  }
  if (is.list(object)) {
    mat <- object$details$proportion_matrix %||% NULL
    if (is.null(mat) && !is.null(object$results)) {
      mat <- deconv_mat(object$results)
    }
    if (!is.null(mat)) {
      return(immune_numeric_matrix(mat, row_label = "object"))
    }
  }
  if (methods::is(object, "SummarizedExperiment")) {
    bundle <- S4Vectors::metadata(object)[["Deconvolution"]]
    if (is.null(bundle)) {
      log_message(
        "Cannot find deconvolution results in {.code metadata(object)[['Deconvolution']]}",
        message_type = "error"
      )
    }
    mat <- bundle$details$proportion_matrix %||% deconv_mat(bundle$results)
    return(immune_numeric_matrix(mat, row_label = "Deconvolution"))
  }
  if (inherits(object, "Seurat")) {
    meta <- object@meta.data
    if (is.null(immune.cols)) {
      immune.cols <- grep(
        "(^RCTD_prop_)|(_prop_)|(_frac_)",
        colnames(meta),
        value = TRUE
      )
      immune.cols <- immune.cols[
        vapply(meta[, immune.cols, drop = FALSE], is.numeric, logical(1))
      ]
    }
    if (length(immune.cols) == 0L) {
      log_message(
        "No immune abundance metadata columns were found. Provide {.arg immune.cols} or {.arg immune.data}.",
        message_type = "error"
      )
    }
    missing_cols <- setdiff(immune.cols, colnames(meta))
    if (length(missing_cols) > 0L) {
      log_message(
        "Some {.arg immune.cols} are missing from metadata: {.val {missing_cols}}",
        message_type = "error"
      )
    }
    mat <- as.matrix(meta[, immune.cols, drop = FALSE])
    colnames(mat) <- sub("^.*?_prop_", "", colnames(mat))
    colnames(mat) <- sub("^.*?_frac_", "", colnames(mat))
    return(immune_numeric_matrix(mat, row_label = "metadata"))
  }
  log_message(
    "Provide {.arg immune.data}, a deconvolution result, a {.cls SummarizedExperiment}, or a {.cls Seurat} object.",
    message_type = "error"
  )
}

immune_numeric_matrix <- function(x, row_label = "matrix") {
  if (is.null(x)) {
    log_message(
      "{.arg {row_label}} is empty.",
      message_type = "error"
    )
  }
  mat <- as.matrix(x)
  dim_names <- dimnames(mat)
  mat <- suppressWarnings(matrix(
    as.numeric(mat),
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dim_names
  ))
  if (is.null(rownames(mat))) {
    log_message(
      "{.arg {row_label}} must have rownames for sample alignment.",
      message_type = "error"
    )
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- paste0("cell_type_", seq_len(ncol(mat)))
  }
  mat[!is.finite(mat)] <- NA_real_
  mat
}

immune_matrix_long <- function(mat, row_name, col_name, value_name) {
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c(row_name, col_name, value_name)
  df[[row_name]] <- factor(df[[row_name]], levels = rownames(mat))
  df[[col_name]] <- factor(df[[col_name]], levels = colnames(mat))
  df
}

immune_scale_matrix <- function(mat, scale = "none") {
  if (identical(scale, "none")) {
    return(mat)
  }
  out <- switch(
    scale,
    row = t(scale(t(mat))),
    column = scale(mat),
    mat
  )
  out <- as.matrix(out)
  out[!is.finite(out)] <- 0
  out
}

resolve_immune_groups <- function(
  object,
  samples,
  group.by = NULL,
  group.data = NULL
) {
  if (!is.null(group.data)) {
    if (is.data.frame(group.data)) {
      sample_col <- intersect(c("sample", "Sample", "id", "ID"), colnames(group.data))[1]
      group_col <- group.by %||%
        setdiff(colnames(group.data), sample_col)[1]
      if (is.na(sample_col) || is.na(group_col)) {
        log_message(
          "{.arg group.data} must contain sample and group columns.",
          message_type = "error"
        )
      }
      return(data.frame(
        sample = as.character(group.data[[sample_col]]),
        group = as.character(group.data[[group_col]]),
        stringsAsFactors = FALSE
      ))
    }
    if (!is.null(names(group.data))) {
      return(data.frame(
        sample = names(group.data),
        group = as.character(group.data),
        stringsAsFactors = FALSE
      ))
    }
  }
  if (is.null(group.by)) {
    return(NULL)
  }
  if (methods::is(object, "SummarizedExperiment")) {
    meta <- as.data.frame(SummarizedExperiment::colData(object))
    meta$sample <- rownames(meta)
    if (!group.by %in% colnames(meta)) {
      log_message(
        "{.arg group.by} is not in {.fn colData(object)}.",
        message_type = "error"
      )
    }
    return(data.frame(
      sample = meta$sample,
      group = as.character(meta[[group.by]]),
      stringsAsFactors = FALSE
    ))
  }
  if (inherits(object, "Seurat")) {
    meta <- object@meta.data
    if (!group.by %in% colnames(meta)) {
      log_message(
        "{.arg group.by} is not in Seurat metadata.",
        message_type = "error"
      )
    }
    return(data.frame(
      sample = rownames(meta),
      group = as.character(meta[[group.by]]),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(sample = samples, group = "All", stringsAsFactors = FALSE)
}

immune_plot_theme <- function(theme_use = "theme_scop", theme_args = list()) {
  if (identical(theme_use, "theme_scop")) {
    theme_use <- "theme_this"
  }
  theme_fun <- tryCatch(get(theme_use, mode = "function"), error = function(e) NULL)
  if (is.null(theme_fun)) {
    return(ggplot2::theme_bw())
  }
  do.call(theme_fun, theme_args)
}

immune_continuous_colors <- function(
  palette = "YlGnBu",
  palcolor = NULL,
  reverse = FALSE
) {
  cols <- palette_colors(
    palette = palette,
    palcolor = palcolor
  )
  if (isTRUE(reverse)) {
    cols <- rev(cols)
  }
  cols
}
