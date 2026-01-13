#' @title Transcript-averaged cell scoring (TACS)
#'
#' @description TACS is a method for plotting a FACS-like plot for two features based on sc-RNA-seq data.
#' For each of two query features, 100 features with similar expression patterns are selected and ranked by their Pearson correlation with the query.
#' In a process akin to compensation, the intersection of the feature lists is removed from each list.
#' The log normalized expression of the resulting features are then averaged within each cell, and the resulting quantities are plotted.
#' This function is based on a simple scheme: choose features similar to the ones specified and average them to reduce the noise.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams CellDimPlot
#' @inheritParams standard_scop
#' @inheritParams FeatureDimPlot
#' @param ref_srt A Seurat object.
#' If your dataset is perturbed in a way that would substantially alter feature-feature correlations,
#' for example if different time points are present or certain cell types are mostly depleted,
#' you can feed in a reference srt, and TACS will choose axes based on the reference data.
#' Default is `NULL`.
#' @param assay Which assay to use. Default is `"RNA"`.
#' @param feature1 Horizontal axis on plot mimics this feature.
#' Character, usually length 1 but possibly longer.
#' @param feature2 Vertical axis on plot mimics this feature. Character, usually length 1 but possibly longer.
#' @param features_predetermined If `FALSE`, plot the sum of many features similar to feature1 instead of feature1 alone (same for feature2).
#' See [GetSimilarFeatures].
#' If `TRUE`, plot the sum of only the features given.
#' @param num_features_add Each axis shows a simple sum of similar features. This is how many (before removing overlap).
#' @param aggregator How to combine correlations when finding similar features.
#' Options: `"sum"` (default), `"min"` (for "and"-like filter), `"max"`, or `"mean"`.
#' @param cutoffs If given, divide plot into four quadrants and annotate with percentages.
#' Can be a numeric vector of length 1 or 2, or a list of two numeric vectors for x and y axes respectively.
#' @param density If `TRUE`, plot contours instead of points.
#' @param bins Number of bins for density plot. Default is `20`.
#' @param h Bandwidth for density plot. Default is `NULL`.
#' @param remove_outliers If `TRUE`, remove outliers from the plot. Default is `FALSE`.
#' @param suffix The suffix of the axis labels. Default is `" expression level"`.
#' @param include_all If `TRUE`, include a panel with all cells. Default is `FALSE`.
#' @param all_color The color of the all cells panel. Default is `"grey20"`.
#' @param quadrants_line_color The color of the quadrants lines. Default is `"grey30"`.
#' @param quadrants_line_type The type of the quadrants lines. Default is `"solid"`.
#' @param quadrants_line_width The width of the quadrants lines. Default is `0.3`.
#' @param quadrants_label_size The size of the quadrants labels. Default is `3`.
#' @param density_alpha The alpha of the density plot. Default is `NULL`.
#' @param ... Additional parameters passed to [ggplot2::stat_density2d].
#'
#' @export
#'
#' @references
#' [Kernfeld et al. paper](https://doi.org/10.1016/j.immuni.2018.04.015),
#' [Github](https://github.com/maehrlab/thymusatlastools2/blob/f8b51ad684d56b2eeda780787eb9ad4ff3003eef/R/data_handling_seurat.R#L271)
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' TACSPlot(
#'   pancreas_sub,
#'   feature1 = "H3f3b",
#'   feature2 = "Eif1",
#'   group.by = "CellType"
#' )
#'
#' TACSPlot(
#'   pancreas_sub,
#'   feature1 = "H3f3b",
#'   feature2 = "Eif1",
#'   group.by = "CellType",
#'   density = TRUE,
#'   include_all = TRUE,
#'   cutoffs = c(3, 2.5)
#' )
#'
#' TACSPlot(
#'   pancreas_sub,
#'   feature1 = "H3f3b",
#'   feature2 = "Eif1",
#'   group.by = "CellType",
#'   density = TRUE,
#'   cutoffs = list(x = c(2, 3), y = c(2.5))
#' )
#'
#' TACSPlot(
#'   pancreas_sub,
#'   feature1 = "H3f3b",
#'   feature2 = "Eif1",
#'   group.by = "SubCellType",
#'   density = TRUE
#' )
TACSPlot <- function(
    srt,
    ref_srt = NULL,
    assay = "RNA",
    layer = "data",
    group.by = NULL,
    feature1,
    feature2,
    cutoffs = NULL,
    density = FALSE,
    palette = "Paired",
    num_features_add = 100,
    features_predetermined = FALSE,
    aggregator = "sum",
    remove_outliers = FALSE,
    aspect.ratio = 1,
    title = NULL,
    subtitle = NULL,
    xlab = NULL,
    ylab = NULL,
    suffix = " expression level",
    legend.position = "right",
    legend.direction = "vertical",
    theme_use = "theme_scop",
    theme_args = list(),
    include_all = FALSE,
    all_color = "grey20",
    quadrants_line_color = "grey30",
    quadrants_line_type = "solid",
    quadrants_line_width = 0.3,
    quadrants_label_size = 3,
    density_alpha = NULL,
    bins = 20,
    h = NULL,
    nrow = NULL,
    ncol = NULL,
    verbose = TRUE,
    ...) {
  if (is.null(ref_srt)) {
    ref_srt <- srt
  }

  if (features_predetermined) {
    feature1_similar <- feature1
    feature2_similar <- feature2
  } else {
    feature1_similar <- GetSimilarFeatures(
      ref_srt, feature1, num_features_add,
      aggregator = aggregator
    )
    feature1_similar <- c(feature1_similar, feature1)
    feature2_similar <- GetSimilarFeatures(
      ref_srt, feature2, num_features_add,
      aggregator = aggregator
    )
    feature2_similar <- c(feature2_similar, feature2)
    shared <- intersect(feature1_similar, feature2_similar)
    feature1_similar <- setdiff(feature1_similar, shared)
    feature2_similar <- setdiff(feature2_similar, shared)
  }

  feature1_score <- rowMeans(
    FetchDataZero(
      srt,
      feature1_similar,
      assay = assay,
      layer = layer,
      verbose = verbose,
      ...
    )
  )
  feature2_score <- rowMeans(
    FetchDataZero(
      srt,
      feature2_similar,
      assay = assay,
      layer = layer,
      verbose = verbose,
      ...
    )
  )
  feature1_suffix <- paste0(feature1[1], suffix)
  feature2_suffix <- paste0(feature2[1], suffix)

  srt <- AddMetaData(
    srt, feature1_score,
    col.name = feature1_suffix
  )
  srt <- AddMetaData(
    srt, feature2_score,
    col.name = feature2_suffix
  )
  plot_df <- FetchData(
    srt,
    c(feature1_suffix, feature2_suffix, group.by),
    assay = assay,
    layer = layer,
    ...
  )
  facet_levels <- FetchData(srt, group.by)[[1]] |>
    factor() |>
    levels()
  colors <- palette_colors(
    facet_levels,
    palette = palette
  )
  if (include_all) {
    plot_df_all <- plot_df
    plot_df_all[[group.by]] <- "All"
    plot_df <- rbind(plot_df, plot_df_all)
    facet_levels <- c(rep("All", include_all), facet_levels)
    colors <- c("All" = all_color, colors)
  }

  if (!is.null(group.by)) {
    plot_df[[group.by]] <- factor(
      plot_df[[group.by]],
      levels = facet_levels,
      ordered = TRUE
    )
    plot_df[[group.by]] <- droplevels(plot_df[[group.by]])
  }

  if (remove_outliers) {
    q1_g1 <- stats::quantile(
      plot_df[[feature1_suffix]], 0.25,
      na.rm = TRUE
    )
    q3_g1 <- stats::quantile(
      plot_df[[feature1_suffix]], 0.75,
      na.rm = TRUE
    )
    iqr_g1 <- q3_g1 - q1_g1
    upper_bound_g1 <- q3_g1 + 1.5 * iqr_g1

    q1_g2 <- stats::quantile(
      plot_df[[feature2_suffix]], 0.25,
      na.rm = TRUE
    )
    q3_g2 <- stats::quantile(
      plot_df[[feature2_suffix]], 0.75,
      na.rm = TRUE
    )
    iqr_g2 <- q3_g2 - q1_g2
    upper_bound_g2 <- q3_g2 + 1.5 * iqr_g2

    original_rows <- nrow(plot_df)
    plot_df <- plot_df[
      plot_df[[feature1_suffix]] <= upper_bound_g1 &
        plot_df[[feature2_suffix]] <= upper_bound_g2,
    ]
    removed_count <- original_rows - nrow(plot_df)
    if (removed_count > 0) {
      message(paste("Removed", removed_count, "outlier cells."))
    }
  }

  if (density) {
    plot_df[[feature1_suffix]] <- jitter(plot_df[[feature1_suffix]])
    plot_df[[feature2_suffix]] <- jitter(plot_df[[feature2_suffix]])
  }

  p <- ggplot(plot_df)
  if (density) {
    p <- p + stat_density2d(
      aes(
        x = .data[[feature1_suffix]],
        y = .data[[feature2_suffix]],
        alpha = density_alpha,
        color = .data[[group.by]]
      ),
      bins = bins,
      h = h
    ) +
      scale_alpha_continuous(range = c(0.4, 1)) +
      scale_color_manual(values = colors)
  } else {
    p <- p + geom_point(
      aes(
        x = .data[[feature1_suffix]],
        y = .data[[feature2_suffix]],
        color = .data[[group.by]]
      )
    ) +
      scale_color_manual(values = colors)
  }
  p <- p + expand_limits(y = 0, x = 0)
  if (!is.null(group.by)) {
    p <- p + facet_wrap(
      stats::as.formula(paste0("~", group.by)),
      nrow = nrow,
      ncol = ncol
    )
  }

  xlab_use <- xlab %||% feature1_suffix
  ylab_use <- ylab %||% feature2_suffix

  p <- p +
    labs(title = title, subtitle = subtitle, x = xlab_use, y = ylab_use) +
    do.call(theme_use, theme_args) +
    theme(
      aspect.ratio = aspect.ratio,
      legend.position = legend.position,
      legend.direction = legend.direction
    ) +
    guides(
      color = guide_legend(
        title = group.by,
        title.hjust = 0,
        order = 1,
        override.aes = list(size = 4, alpha = 1)
      )
    )

  if (!is.null(cutoffs)) {
    p <- add_quadrants(
      p,
      feature1_suffix = feature1_suffix,
      feature2_suffix = feature2_suffix,
      cutoffs = cutoffs,
      quadrants_line_color = quadrants_line_color,
      quadrants_line_type = quadrants_line_type,
      quadrants_line_width = quadrants_line_width,
      quadrants_label_size = quadrants_label_size,
      group.by = group.by
    )
  }

  return(p)
}

add_quadrants <- function(
    p,
    feature1_suffix,
    feature2_suffix,
    cutoffs,
    quadrants_line_color = "grey30",
    quadrants_line_type = "solid",
    quadrants_line_width = 0.5,
    quadrants_label_size = 3,
    group.by = NULL) {
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
    p <- p + geom_vline(
      xintercept = cutoffs_x,
      linetype = quadrants_line_type,
      color = quadrants_line_color,
      linewidth = quadrants_line_width
    )
  }
  if (!is.null(cutoffs_y)) {
    p <- p + geom_hline(
      yintercept = cutoffs_y,
      linetype = quadrants_line_type,
      color = quadrants_line_color,
      linewidth = quadrants_line_width
    )
  }

  plot_data <- p$data
  plot_limits <- ggplot_build(p)$layout$panel_params[[1]]

  x_breaks <- c(plot_limits$x.range[1], cutoffs_x, plot_limits$x.range[2])
  x_breaks <- unique(sort(x_breaks))
  y_breaks <- c(plot_limits$y.range[1], cutoffs_y, plot_limits$y.range[2])
  y_breaks <- unique(sort(y_breaks))

  x_breaks <- x_breaks[x_breaks >= plot_limits$x.range[1] & x_breaks <= plot_limits$x.range[2]]
  y_breaks <- y_breaks[y_breaks >= plot_limits$y.range[1] & y_breaks <= plot_limits$y.range[2]]
  x_breaks <- unique(sort(x_breaks))
  y_breaks <- unique(sort(y_breaks))

  x_pos <- (utils::head(x_breaks, -1) + utils::tail(x_breaks, -1)) / 2
  y_pos <- (utils::head(y_breaks, -1) + utils::tail(y_breaks, -1)) / 2

  plot_data$x_cat <- as.integer(
    cut(
      plot_data[[feature1_suffix]],
      breaks = x_breaks,
      include.lowest = TRUE
    )
  )
  plot_data$y_cat <- as.integer(
    cut(
      plot_data[[feature2_suffix]],
      breaks = y_breaks,
      include.lowest = TRUE
    )
  )

  if (is.null(group.by)) {
    counts <- as.data.frame(
      table(
        plot_data[, c("x_cat", "y_cat")]
      ),
      stringsAsFactors = FALSE
    )
    counts$value <- percentify(counts$Freq)
  } else {
    counts <- as.data.frame(
      table(
        plot_data[, c(group.by, "x_cat", "y_cat")]
      ),
      stringsAsFactors = FALSE
    )
    counts_list <- split(counts, counts[[group.by]])

    percentages_list <- lapply(
      counts_list, function(df) {
        df$value <- percentify(df$Freq)
        return(df)
      }
    )
    counts <- do.call(rbind, percentages_list)
  }

  counts$x_cat <- as.integer(counts$x_cat)
  counts$y_cat <- as.integer(counts$y_cat)

  annot_df <- counts[counts$Freq > 0, ]
  if (nrow(annot_df) == 0) {
    return(p)
  }
  annot_df$x <- x_pos[annot_df$x_cat]
  annot_df$y <- y_pos[annot_df$y_cat]
  annot_df$label <- paste0(round(annot_df$value, 1), "%")

  p <- p + geom_text(
    data = annot_df,
    aes(x = .data$x, y = .data$y, label = .data$label),
    size = quadrants_label_size
  )

  return(p)
}

percentify <- function(x) {
  return(100 * round(div_by_sum(x), 3))
}

div_by_sum <- function(x) {
  if (sum(x) == 0) 0 * x else x / sum(x)
}

#' @title Find features with expression patterns similar to provided features
#'
#' @md
#' @inheritParams TACSPlot
#' @inheritParams FetchDataZero
#' @param n An integer; number of results to return.
#' @param features_use A character vector of features eligible to be returned.
#' @param anticorr Whether to allow negatively correlated features.
#' Default is `FALSE`.
#'
#' @return character vector.
#' @export
GetSimilarFeatures <- function(
    srt,
    features,
    n,
    features_use = rownames(srt),
    anticorr = FALSE,
    aggregator = "sum",
    assay = "RNA",
    layer = "data") {
  if (!all(features %in% rownames(srt))) {
    log_message(
      "Some of your features have no data available.",
      message_type = "warning"
    )
  }
  features <- intersect(features, rownames(srt))
  if (!all(features_use %in% rownames(srt))) {
    log_message(
      "Some of your genes have no data available.",
      message_type = "warning"
    )
  }
  features_use <- intersect(features_use, rownames(srt))
  data_use <- GetAssayData5(
    srt,
    assay = assay,
    layer = layer
  )[features_use, , drop = FALSE] |> Matrix::Matrix(sparse = TRUE)

  gene_averages <- Matrix::rowMeans(data_use)
  squares <- Matrix::rowSums(data_use^2)
  gene_vars <- squares / ncol(data_use) - gene_averages^2
  gene_sds <- sqrt(gene_vars)

  query <- as_matrix(Seurat::FetchData(srt, features))
  covariances_all <- data_use %*% query - gene_averages %o% colSums(query)
  correlations_all <- Matrix::Diagonal(x = 1 / gene_sds) %*% covariances_all

  if (ncol(correlations_all) == 1) {
    correlation <- as.vector(correlations_all)
  } else {
    correlation <- switch(
      EXPR = aggregator,
      "sum" = Matrix::rowSums(correlations_all),
      "min" = apply(as_matrix(correlations_all), 1, min),
      "max" = apply(as_matrix(correlations_all), 1, max),
      "mean" = Matrix::rowMeans(correlations_all),
      log_message(
        "{.arg aggregator} must be one of: 'sum', 'min', 'max', 'mean'",
        message_type = "error"
      )
    )
  }
  correlation <- correlation[setdiff(names(correlation), features)]
  if (anticorr) {
    similar_genes <- names(
      sort(abs(correlation), decreasing = TRUE)[seq_len(n)]
    )
  } else {
    similar_genes <- names(
      sort(correlation, decreasing = TRUE)[seq_len(n)]
    )
  }
  return(similar_genes)
}

#' @title FetchData but with zeroes for unavailable genes
#'
#' @md
#' @inheritParams TACSPlot
#' @param features A character vector of feature names.
#' @param ... Other arguments to pass to [Seurat::FetchData].
#'
#' @export
FetchDataZero <- function(
    srt,
    features,
    assay = "RNA",
    layer = "data",
    verbose = TRUE,
    ...) {
  features <- features[stats::complete.cases(features)]
  avail <- intersect(features, rownames(srt))
  unavail <- setdiff(features, rownames(srt))
  log_message(
    "Some features are not available. Returning zeroes",
    message_type = "warning",
    verbose = verbose && length(unavail) > 0
  )

  to_return <- Seurat::FetchData(
    srt,
    avail,
    assay = assay,
    layer = layer,
    ...
  )
  pad <- as.data.frame(
    matrix(0,
      nrow = nrow(to_return),
      ncol = length(unavail),
      dimnames = list(
        rownames(to_return),
        unavail
      )
    )
  )
  to_return <- cbind(to_return, pad)
  to_return <- to_return[, features, drop = FALSE]
  return(to_return)
}
