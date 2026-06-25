#' @title ESTIMATE score plots
#'
#' @description
#' Visualize ESTIMATE stromal, immune, combined ESTIMATE, and tumor-purity
#' scores from a `RunESTIMATE()` result.
#'
#' @md
#' @param object Optional `RunESTIMATE()` bundle, `SummarizedExperiment`, or
#' `Seurat` object containing ESTIMATE results.
#' @param score.data Optional score matrix or data frame with samples in rows.
#' @param group.by Optional grouping column for grouped violin and box plots.
#' @param group.data Optional named vector or data frame containing sample
#' groups.
#' @param plot_type Plot type.
#' @param scores ESTIMATE score columns to plot.
#' @param add_stat Whether to add group comparison labels to violin or box
#' plots. Requires `ggpubr`.
#' @param ... Additional plotting arguments.
#'
#' @return A `ggplot` object.
#'
#' @export
EstimateScorePlot <- function(
  object = NULL,
  score.data = NULL,
  group.by = NULL,
  group.data = NULL,
  plot_type = c("violin", "box", "heatmap", "cor"),
  scores = c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity"),
  add_stat = TRUE,
  ...
) {
  plot_type <- match.arg(plot_type)
  resolved <- resolve_estimate_scores(object = object, score.data = score.data)
  mat <- estimate_select_scores(resolved$matrix, scores = scores)
  dots <- list(...)
  theme_obj <- immune_plot_theme(
    theme_use = dots$theme_use %||% "theme_scop",
    theme_args = dots$theme_args %||% list()
  )

  if (identical(plot_type, "cor")) {
    if (ncol(mat) < 2L) {
      log_message("Need at least two ESTIMATE scores for a correlation plot.", message_type = "error")
    }
    cor_method <- dots$cor_method %||% "spearman"
    cor_mat <- stats::cor(mat, method = cor_method, use = "pairwise.complete.obs")
    cor_mat[!is.finite(cor_mat)] <- 0
    df <- estimate_matrix_long(cor_mat, "score_type_1", "score_type_2", "cor")
    return(
      ggplot2::ggplot(df, ggplot2::aes(score_type_2, score_type_1, fill = cor)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.5) +
        ggplot2::scale_fill_gradientn(
          colors = immune_continuous_colors(
            palette = "RdBu",
            palcolor = dots$heatmap_palcolor %||% NULL,
            reverse = TRUE
          ),
          limits = c(-1, 1),
          na.value = "grey90",
          name = "Correlation"
        ) +
        ggplot2::coord_equal() +
        ggplot2::labs(
          title = dots$title %||% NULL,
          subtitle = dots$subtitle %||% NULL,
          x = NULL,
          y = NULL
        ) +
        theme_obj +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid = ggplot2::element_blank(),
          legend.position = dots$legend.position %||% "right"
        )
    )
  }

  if (identical(plot_type, "heatmap")) {
    scale <- dots$scale %||% "none"
    scale <- match.arg(scale, c("none", "row", "column"))
    mat_plot <- immune_scale_matrix(mat, scale = scale)
    df <- estimate_matrix_long(mat_plot, "sample", "score_type", "score")
    return(
      ggplot2::ggplot(df, ggplot2::aes(score_type, sample, fill = score)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.35) +
        ggplot2::scale_fill_gradientn(
          colors = immune_continuous_colors(
            palette = dots$heatmap_palette %||% "YlGnBu",
            palcolor = dots$heatmap_palcolor %||% NULL
          ),
          na.value = "grey90",
          name = "Score"
        ) +
        ggplot2::labs(
          title = dots$title %||% NULL,
          subtitle = dots$subtitle %||% NULL,
          x = NULL,
          y = NULL
        ) +
        theme_obj +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.y = if (isTRUE(dots$show_sample_names %||% FALSE)) {
            ggplot2::element_text()
          } else {
            ggplot2::element_blank()
          },
          axis.ticks.y = if (isTRUE(dots$show_sample_names %||% FALSE)) {
            ggplot2::element_line()
          } else {
            ggplot2::element_blank()
          },
          panel.grid = ggplot2::element_blank(),
          legend.position = dots$legend.position %||% "right"
        )
    )
  }

  df <- estimate_matrix_long(mat, "sample", "score_type", "score")
  group_df <- resolve_estimate_groups(
    object = object,
    bundle = resolved$bundle,
    samples = rownames(mat),
    group.by = group.by,
    group.data = group.data
  )
  df$group <- group_df$group[match(df$sample, group_df$sample)]
  df$group[is.na(df$group)] <- "All"
  df$group <- factor(df$group, levels = unique(df$group))
  fill_cols <- palette_colors(
    levels(df$group),
    palette = dots$palette %||% "Chinese",
    palcolor = dots$palcolor %||% NULL
  )

  plot <- ggplot2::ggplot(df, ggplot2::aes(score_type, score, fill = group))
  if (identical(plot_type, "violin")) {
    plot <- plot +
      ggplot2::geom_violin(
        trim = FALSE,
        scale = "width",
        linewidth = 0.35,
        alpha = dots$violin_alpha %||% 0.82,
        position = ggplot2::position_dodge(width = 0.78)
      ) +
      ggplot2::geom_boxplot(
        width = dots$box_width %||% 0.16,
        outlier.shape = NA,
        color = "white",
        linewidth = 0.3,
        position = ggplot2::position_dodge(width = 0.78)
      )
  } else {
    plot <- plot +
      ggplot2::geom_boxplot(
        outlier.shape = NA,
        width = dots$box_width %||% 0.64,
        linewidth = 0.38,
        alpha = dots$box_alpha %||% 0.92,
        color = "grey20",
        position = ggplot2::position_dodge2(width = 0.78, preserve = "single")
      )
  }
  plot <- plot +
    ggplot2::geom_point(
      ggplot2::aes(color = group),
      position = ggplot2::position_jitterdodge(
        jitter.width = dots$jitter.width %||% 0.12,
        jitter.height = 0,
        dodge.width = 0.78
      ),
      size = dots$pt.size %||% 1.35,
      alpha = dots$pt.alpha %||% 0.72,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = fill_cols, drop = FALSE) +
    ggplot2::scale_color_manual(values = fill_cols, drop = FALSE) +
    ggplot2::labs(
      title = dots$title %||% NULL,
      subtitle = dots$subtitle %||% NULL,
      x = NULL,
      y = "ESTIMATE score",
      fill = group.by %||% "Group"
    ) +
    theme_obj +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = dots$legend.position %||% "right",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (isTRUE(add_stat) && length(levels(df$group)) > 1L) {
    check_r("ggpubr", verbose = FALSE)
    plot <- plot + ggpubr::stat_compare_means(
      ggplot2::aes(group = group),
      method = dots$pairwise_method %||% "wilcox.test",
      label = dots$sig_label %||% "p.signif",
      size = dots$sig_labelsize %||% 3.5,
      hide.ns = TRUE
    )
  }
  plot
}

#' @title Gene and ESTIMATE score relationship plots
#'
#' @description
#' Compare ESTIMATE scores between target-gene high and low expression groups,
#' or draw continuous gene-expression correlations with ESTIMATE scores.
#'
#' @inheritParams EstimateScorePlot
#' @param gene.data Optional gene expression matrix with genes in rows and
#' samples in columns.
#' @param features Target genes to plot.
#' @param assay Assay used for `Seurat` or `SummarizedExperiment` expression.
#' @param layer Assay layer used for `Seurat` expression.
#' @param split Split method for high/low groups in violin plots.
#' @param quantile Quantile cutoff used when `split = "quantile"`.
#' @param cor_method Correlation method for scatter plots.
#'
#' @return A `ggplot` object.
#'
#' @export
EstimateGenePlot <- function(
  object = NULL,
  score.data = NULL,
  gene.data = NULL,
  features,
  assay = NULL,
  layer = "data",
  plot_type = c("violin", "scatter"),
  split = c("median", "quantile"),
  quantile = 0.5,
  cor_method = c("spearman", "pearson", "kendall"),
  add_stat = TRUE,
  ...
) {
  plot_type <- match.arg(plot_type)
  split <- match.arg(split)
  cor_method <- match.arg(cor_method)
  if (missing(features) || length(features) == 0L) {
    log_message("{.arg features} must contain at least one target gene.", message_type = "error")
  }
  if (length(quantile) != 1L || !is.finite(quantile) || quantile <= 0 || quantile >= 1) {
    log_message("{.arg quantile} must be a number between 0 and 1.", message_type = "error")
  }

  resolved <- resolve_estimate_scores(object = object, score.data = score.data)
  score_mat <- estimate_select_scores(resolved$matrix, scores = estimate_score_columns())
  gene_mat <- resolve_estimate_gene_matrix(
    object = object,
    gene.data = gene.data,
    features = features,
    assay = assay,
    layer = layer,
    bundle = resolved$bundle
  )

  common_samples <- intersect(rownames(score_mat), rownames(gene_mat))
  if (length(common_samples) < 3L) {
    log_message(
      "Need at least three matched samples between gene expression and ESTIMATE scores.",
      message_type = "error"
    )
  }
  score_mat <- score_mat[common_samples, , drop = FALSE]
  gene_mat <- gene_mat[common_samples, , drop = FALSE]

  dots <- list(...)
  theme_obj <- immune_plot_theme(
    theme_use = dots$theme_use %||% "theme_scop",
    theme_args = dots$theme_args %||% list()
  )

  if (identical(plot_type, "scatter")) {
    df <- merge(
      estimate_matrix_long(score_mat, "sample", "score_type", "score"),
      estimate_matrix_long(gene_mat, "sample", "feature", "expression"),
      by = "sample"
    )
    label_df <- estimate_gene_cor_labels(
      df = df,
      method = cor_method
    )
    return(
      ggplot2::ggplot(df, ggplot2::aes(expression, score)) +
        ggplot2::geom_point(
          color = dots$point_color %||% "#2B6CB0",
          size = dots$pt.size %||% 1.6,
          alpha = dots$pt.alpha %||% 0.75
        ) +
        ggplot2::geom_smooth(
          method = "lm",
          se = FALSE,
          linewidth = 0.5,
          color = dots$smooth_color %||% "#D1495B"
        ) +
        ggplot2::geom_text(
          data = label_df,
          ggplot2::aes(label = label, x = x, y = y),
          inherit.aes = FALSE,
          hjust = 0,
          vjust = 1,
          size = dots$label_size %||% 3
        ) +
        ggplot2::facet_grid(feature ~ score_type, scales = "free") +
        ggplot2::labs(
          title = dots$title %||% NULL,
          subtitle = dots$subtitle %||% NULL,
          x = "Expression",
          y = "ESTIMATE score"
        ) +
        theme_obj
    )
  }

  group_df <- estimate_gene_split_groups(
    gene_mat = gene_mat,
    split = split,
    quantile = quantile
  )
  df <- merge(
    estimate_matrix_long(score_mat, "sample", "score_type", "score"),
    group_df,
    by = "sample"
  )
  df$group <- factor(df$group, levels = c("Low", "High"))
  fill_cols <- palette_colors(
    levels(df$group),
    palette = dots$palette %||% "Chinese",
    palcolor = dots$palcolor %||% c("#2B6CB0", "#D1495B")
  )
  plot <- ggplot2::ggplot(df, ggplot2::aes(score_type, score, fill = group)) +
    ggplot2::geom_violin(
      trim = FALSE,
      scale = "width",
      linewidth = 0.35,
      alpha = dots$violin_alpha %||% 0.82,
      position = ggplot2::position_dodge(width = 0.78)
    ) +
    ggplot2::geom_boxplot(
      width = dots$box_width %||% 0.16,
      outlier.shape = NA,
      color = "white",
      linewidth = 0.3,
      position = ggplot2::position_dodge(width = 0.78)
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = group),
      position = ggplot2::position_jitterdodge(
        jitter.width = dots$jitter.width %||% 0.12,
        jitter.height = 0,
        dodge.width = 0.78
      ),
      size = dots$pt.size %||% 1.35,
      alpha = dots$pt.alpha %||% 0.72,
      show.legend = FALSE
    ) +
    ggplot2::facet_wrap(~feature, scales = "free_x") +
    ggplot2::scale_fill_manual(values = fill_cols, drop = FALSE) +
    ggplot2::scale_color_manual(values = fill_cols, drop = FALSE) +
    ggplot2::labs(
      title = dots$title %||% NULL,
      subtitle = dots$subtitle %||% NULL,
      x = NULL,
      y = "ESTIMATE score",
      fill = "Expression"
    ) +
    theme_obj +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  if (isTRUE(add_stat)) {
    check_r("ggpubr", verbose = FALSE)
    plot <- plot + ggpubr::stat_compare_means(
      ggplot2::aes(group = group),
      method = dots$pairwise_method %||% "wilcox.test",
      label = dots$sig_label %||% "p.signif",
      size = dots$sig_labelsize %||% 3.5,
      hide.ns = TRUE
    )
  }
  plot
}

resolve_estimate_scores <- function(object = NULL, score.data = NULL) {
  bundle <- NULL
  if (!is.null(score.data)) {
    mat <- score.data
  } else if (is.list(object) && !is.null(object$scores)) {
    bundle <- object
    mat <- object$scores
  } else if (methods::is(object, "SummarizedExperiment")) {
    bundle <- S4Vectors::metadata(object)[["ESTIMATE"]]
    if (is.null(bundle)) {
      log_message(
        "Cannot find ESTIMATE results in {.code metadata(object)[['ESTIMATE']]}",
        message_type = "error"
      )
    }
    mat <- bundle$scores
  } else if (inherits(object, "Seurat")) {
    bundle <- object@tools$ESTIMATE
    if (is.null(bundle)) {
      log_message(
        "Cannot find ESTIMATE results in {.code object@tools$ESTIMATE}",
        message_type = "error"
      )
    }
    mat <- bundle$scores
  } else if (inherits(object, c("matrix", "data.frame", "Matrix"))) {
    mat <- object
  } else {
    log_message(
      "Provide {.arg score.data}, a {.fn RunESTIMATE} result, a {.cls SummarizedExperiment}, or a {.cls Seurat} object.",
      message_type = "error"
    )
  }
  mat <- estimate_numeric_matrix(mat, row_label = "ESTIMATE scores")
  list(matrix = mat, bundle = bundle)
}

estimate_select_scores <- function(mat, scores) {
  scores <- intersect(scores, colnames(mat))
  if (length(scores) == 0L) {
    log_message("No requested ESTIMATE score columns are available.", message_type = "error")
  }
  mat[, scores, drop = FALSE]
}

estimate_numeric_matrix <- function(x, row_label = "matrix") {
  mat <- as.matrix(x)
  dim_names <- dimnames(mat)
  mat <- suppressWarnings(matrix(
    as.numeric(mat),
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dim_names
  ))
  if (is.null(rownames(mat))) {
    log_message("{.arg {row_label}} must have rownames for sample alignment.", message_type = "error")
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- paste0("score_", seq_len(ncol(mat)))
  }
  mat[!is.finite(mat)] <- NA_real_
  mat
}

estimate_matrix_long <- function(mat, row_name, col_name, value_name) {
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c(row_name, col_name, value_name)
  df[[row_name]] <- factor(df[[row_name]], levels = rownames(mat))
  df[[col_name]] <- factor(df[[col_name]], levels = colnames(mat))
  df
}

resolve_estimate_groups <- function(
  object,
  bundle = NULL,
  samples,
  group.by = NULL,
  group.data = NULL
) {
  if (!is.null(group.data)) {
    if (is.data.frame(group.data)) {
      sample_col <- intersect(c("sample", "Sample", "id", "ID"), colnames(group.data))[1]
      group_col <- group.by %||% setdiff(colnames(group.data), sample_col)[1]
      if (is.na(sample_col) || is.na(group_col)) {
        log_message("{.arg group.data} must contain sample and group columns.", message_type = "error")
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
  sample_meta <- bundle$details$sample_metadata %||% NULL
  if (!is.null(sample_meta)) {
    group_col <- group.by %||% intersect(c("group", "condition"), colnames(sample_meta))[1]
    if (!is.na(group_col) && group_col %in% colnames(sample_meta)) {
      return(data.frame(
        sample = as.character(sample_meta$sample %||% rownames(sample_meta)),
        group = as.character(sample_meta[[group_col]]),
        stringsAsFactors = FALSE
      ))
    }
  }
  if (!is.null(group.by) && methods::is(object, "SummarizedExperiment")) {
    meta <- as.data.frame(SummarizedExperiment::colData(object))
    if (!group.by %in% colnames(meta)) {
      log_message("{.arg group.by} is not in {.fn colData(object)}.", message_type = "error")
    }
    return(data.frame(
      sample = rownames(meta),
      group = as.character(meta[[group.by]]),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(sample = samples, group = "All", stringsAsFactors = FALSE)
}

resolve_estimate_gene_matrix <- function(
  object = NULL,
  gene.data = NULL,
  features,
  assay = NULL,
  layer = "data",
  bundle = NULL
) {
  if (!is.null(gene.data)) {
    mat <- as.matrix(gene.data)
  } else if (methods::is(object, "SummarizedExperiment")) {
    assay_use <- assay %||% SummarizedExperiment::assayNames(object)[1]
    mat <- SummarizedExperiment::assay(object, assay_use)
  } else if (inherits(object, "Seurat")) {
    params <- bundle$parameters %||% list()
    pseudo <- estimate_seurat_pseudobulk(
      object = object,
      assay = assay %||% params$assay,
      layer = layer %||% params$layer %||% "data",
      sample.by = params$sample.by,
      group.by = params$group.by,
      aggregate_fun = params$aggregate_fun %||% "mean",
      features = features
    )
    mat <- pseudo$matrix
  } else if (inherits(object, c("matrix", "data.frame", "Matrix"))) {
    mat <- as.matrix(object)
  } else {
    log_message(
      "Provide {.arg gene.data}, an expression matrix, a {.cls SummarizedExperiment}, or a {.cls Seurat} object.",
      message_type = "error"
    )
  }
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    log_message("{.arg gene.data} must have gene rownames and sample colnames.", message_type = "error")
  }
  feature_key <- toupper(rownames(mat))
  features_key <- toupper(features)
  keep <- feature_key %in% features_key
  if (!any(keep)) {
    log_message("No target genes are available in expression data.", message_type = "error")
  }
  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- feature_key[keep]
  mat <- estimate_prepare_matrix(mat)
  mat <- t(mat)
  colnames(mat) <- rownames(estimate_prepare_matrix(t(mat)))
  mat <- mat[, intersect(features_key, colnames(mat)), drop = FALSE]
  if (ncol(mat) == 0L) {
    log_message("No variable target genes remain after preprocessing.", message_type = "error")
  }
  mat
}

estimate_gene_split_groups <- function(gene_mat, split = "median", quantile = 0.5) {
  out <- lapply(colnames(gene_mat), function(feature) {
    x <- gene_mat[, feature]
    cutoff <- if (identical(split, "median")) {
      stats::median(x, na.rm = TRUE)
    } else {
      stats::quantile(x, probs = quantile, na.rm = TRUE, names = FALSE)
    }
    data.frame(
      sample = rownames(gene_mat),
      feature = feature,
      expression = as.numeric(x),
      group = ifelse(x > cutoff, "High", "Low"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

estimate_gene_cor_labels <- function(df, method = "spearman") {
  splits <- split(df, interaction(df$feature, df$score_type, drop = TRUE))
  out <- lapply(splits, function(x) {
    keep <- is.finite(x$expression) & is.finite(x$score)
    x <- x[keep, , drop = FALSE]
    label <- "rho = NA\np = NA"
    if (nrow(x) >= 3L && stats::var(x$expression) > 0 && stats::var(x$score) > 0) {
      ct <- suppressWarnings(stats::cor.test(x$expression, x$score, method = method))
      label <- paste0(
        "r = ", sprintf("%.2f", unname(ct$estimate)),
        "\np = ", format.pval(ct$p.value, digits = 2, eps = 1e-3)
      )
    }
    data.frame(
      feature = x$feature[1],
      score_type = x$score_type[1],
      x = min(x$expression, na.rm = TRUE),
      y = max(x$score, na.rm = TRUE),
      label = label,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}
