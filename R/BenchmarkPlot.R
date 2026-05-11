#' @title Plot benchmark metrics
#'
#' @description
#' Visualize benchmark results stored in a `Seurat` object or a summary
#' `data.frame`. Per-cell metrics such as LISI are shown as feature plots and/or
#' boxplots. Summary metrics such as integration or mapping benchmark scores are
#' shown as barplots or a `funkyheatmap`.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams CellDimPlot
#' @param srt A `Seurat` object.
#' @param data Optional summary benchmark `data.frame` containing at least
#' `metric` and `value`, and optionally `method` and `workflow`.
#' @param features Metadata columns containing per-cell benchmark scores.
#' Default is `NULL`.
#' @param metrics One or more summary metric names to visualize. Default is
#' `NULL`, which uses all available summary metrics.
#' @param tool_name Tool entries created by benchmark-related workflows.
#' This can be a character vector. For per-cell metrics, benchmark columns are
#' resolved from tool entries that contain `colnames`; for summary metrics,
#' entries containing `summary` or `metrics$summary` are used.
#' @param reduction Dimensional reduction used for per-cell feature plots.
#' Default is `NULL`, which uses the reduction stored in `tool_name` when
#' available, otherwise [DefaultReduction()].
#' @param plot_type Plot type. One of `"auto"`, `"feature"`, `"boxplot"`,
#' `"bar"`, or `"funkyheatmap"`.
#' @param plot_boxplot Whether to add the summary boxplot when per-cell metrics
#' are shown. Default is `TRUE`.
#' @param boxplot_jitter Whether to overlay jittered points on the boxplot.
#' Default is `FALSE`.
#'
#' @return
#' A ggplot, patchwork plot, or `funkyheatmap` object depending on the selected
#' mode. If `combine = FALSE` in per-cell mode, a named list of plots is
#' returned.
#'
#' @examples
#' if (requireNamespace("labeling", quietly = TRUE)) {
#'   metrics_df <- data.frame(
#'     method = c("Raw", "Raw", "Harmony", "Harmony"),
#'     metric = c("batch_ASW_mixing", "celltype_ASW", "batch_ASW_mixing", "celltype_ASW"),
#'     value = c(0.42, 0.71, 0.68, 0.66)
#'   )
#'   BenchmarkPlot(
#'     data = metrics_df,
#'     plot_type = "bar"
#'   )
#'
#'   data("pbmcmultiome_sub", package = "scop")
#'   pbmcmultiome_sub[["MethodA_batch_LISI"]] <-
#'     seq_len(ncol(pbmcmultiome_sub)) / ncol(pbmcmultiome_sub)
#'   pbmcmultiome_sub[["MethodB_batch_LISI"]] <-
#'     rev(pbmcmultiome_sub[["MethodA_batch_LISI", drop = TRUE]])
#'   BenchmarkPlot(
#'     pbmcmultiome_sub,
#'     features = c("MethodA_batch_LISI", "MethodB_batch_LISI"),
#'     plot_type = "boxplot"
#'   )
#' }
#'
#' @export
BenchmarkPlot <- function(
  srt = NULL,
  data = NULL,
  features = NULL,
  metrics = NULL,
  tool_name = NULL,
  reduction = NULL,
  plot_type = c("auto", "feature", "boxplot", "bar", "funkyheatmap"),
  plot_boxplot = TRUE,
  boxplot_jitter = FALSE,
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  pt.size = NULL,
  pt.alpha = 1,
  palette = "Chinese",
  palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  verbose = TRUE,
  ...
) {
  plot_type <- match.arg(plot_type)

  if (is.null(data)) {
    if (!inherits(srt, "Seurat")) {
      log_message(
        "{.arg srt} must be a {.cls Seurat} when {.arg data} is not provided",
        message_type = "error"
      )
    }
  }

  summary_mode_requested <- !is.null(data) ||
    !is.null(metrics) ||
    plot_type %in% c("bar", "funkyheatmap")

  if (isTRUE(summary_mode_requested)) {
    summary_df <- benchmark_collect_summary_data(
      srt = srt,
      data = data,
      tool_name = tool_name,
      metrics = metrics
    )
    if (identical(plot_type, "auto")) {
      plot_type <- if (
        requireNamespace("funkyheatmap", quietly = TRUE) &&
          length(unique(summary_df$method)) > 1 &&
          length(unique(summary_df$metric)) > 1
      ) {
        "funkyheatmap"
      } else {
        "bar"
      }
    }
    if (identical(plot_type, "funkyheatmap")) {
      return(benchmark_summary_funkyheatmap(
        summary_df = summary_df,
        metrics = metrics,
        verbose = verbose
      ))
    }
    return(benchmark_summary_barplot(
      summary_df = summary_df,
      palette = palette,
      palcolor = palcolor,
      theme_use = theme_use,
      theme_args = theme_args
    ))
  }

  benchmark_feature_plot(
    srt = srt,
    features = features,
    tool_name = tool_name,
    reduction = reduction,
    plot_type = plot_type,
    plot_boxplot = plot_boxplot,
    boxplot_jitter = boxplot_jitter,
    combine = combine,
    nrow = nrow,
    ncol = ncol,
    byrow = byrow,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    palette = palette,
    palcolor = palcolor,
    theme_use = theme_use,
    theme_args = theme_args,
    verbose = verbose,
    ...
  )
}

benchmark_feature_plot <- function(
  srt,
  features = NULL,
  tool_name = NULL,
  reduction = NULL,
  plot_type = "auto",
  plot_boxplot = TRUE,
  boxplot_jitter = FALSE,
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  pt.size = NULL,
  pt.alpha = 1,
  palette = "Chinese",
  palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  verbose = TRUE,
  ...
) {
  tool_res <- NULL
  if (!is.null(tool_name)) {
    tool_name <- as.character(tool_name)
    tool_name_use <- tool_name[[1]]
    if (!tool_name_use %in% names(srt@tools)) {
      log_message(
        "Tool entry {.val {tool_name_use}} not found in {.cls Seurat}",
        message_type = "error"
      )
    }
    tool_res <- srt@tools[[tool_name_use]]
  }

  features_missing <- is.null(features)
  if (features_missing) {
    if (!is.null(tool_res) && "colnames" %in% names(tool_res)) {
      features <- tool_res[["colnames"]]
    } else {
      features <- grep("_LISI$", colnames(srt@meta.data), value = TRUE)
    }
  }
  features <- unique(features)
  if (length(features) == 0) {
    log_message(
      "No per-cell benchmark columns found. Please provide {.arg features} or a valid {.arg tool_name}.",
      message_type = "error"
    )
  }
  if (!all(features %in% colnames(srt@meta.data))) {
    missing_cols <- setdiff(features, colnames(srt@meta.data))
    log_message(
      "The following benchmark columns are missing: {.val {missing_cols}}",
      message_type = "error"
    )
  }

  if (
    identical(plot_type, "auto") &&
      length(features) < 2
  ) {
    plot_boxplot <- FALSE
  }

  if (
    isTRUE(plot_boxplot) &&
      isTRUE(features_missing) &&
      !is.null(tool_res) &&
      length(tool_res[["label_colnames"]] %||% character()) > 1
  ) {
    log_message(
      "Skip default benchmark boxplot because the selected tool contains multiple label types. Please provide comparable {.arg features}.",
      message_type = "warning",
      verbose = verbose
    )
    plot_boxplot <- FALSE
  }

  if (identical(plot_type, "boxplot")) {
    return(benchmark_feature_boxplot(
      srt = srt,
      features = features,
      palette = palette,
      palcolor = palcolor,
      boxplot_jitter = boxplot_jitter,
      theme_use = theme_use,
      theme_args = theme_args,
      verbose = verbose
    ))
  }

  tool_reduction <- tool_res[["reduction"]] %||% tool_res[["reductions"]] %||% NULL
  if (length(tool_reduction) > 1) {
    tool_reduction <- NULL
  }
  reduction <- reduction %||% tool_reduction %||% DefaultReduction(srt)
  if (!reduction %in% SeuratObject::Reductions(srt)) {
    log_message(
      "Reduction {.val {reduction}} not found in {.cls Seurat}",
      message_type = "error"
    )
  }

  feature_title_map <- benchmark_clean_feature_labels(features)
  feature_mean_df <- data.frame(
    feature = features,
    title = unname(feature_title_map[features]),
    mean_score = vapply(
      features,
      function(feature) mean(srt@meta.data[[feature]], na.rm = TRUE),
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
  feature_mean_df <- feature_mean_df[order(feature_mean_df$mean_score, decreasing = FALSE), , drop = FALSE]
  feature_plot_features <- feature_mean_df$feature[!duplicated(feature_mean_df$title)]

  plots <- list()
  for (feature in feature_plot_features) {
    plots[[feature]] <- FeatureDimPlot(
      srt,
      features = feature,
      reduction = reduction,
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      combine = FALSE,
      title = feature_title_map[[feature]],
      subtitle = NULL,
      show_stat = FALSE,
      ...
    )[[1]] +
      ggplot2::labs(title = feature_title_map[[feature]], subtitle = NULL) +
      ggplot2::theme(
        plot.subtitle = ggplot2::element_blank(),
        strip.text = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank()
      )
  }

  if (!identical(plot_type, "feature") && isTRUE(plot_boxplot)) {
    plots[["benchmark_boxplot"]] <- benchmark_feature_boxplot(
      srt = srt,
      features = features,
      palette = palette,
      palcolor = palcolor,
      boxplot_jitter = boxplot_jitter,
      theme_use = theme_use,
      theme_args = theme_args,
      verbose = verbose
    )
  }

  if (!isTRUE(combine)) {
    return(plots)
  }
  if (length(plots) == 1) {
    return(plots[[1]])
  }

  do.call(
    patchwork::wrap_plots,
    list(
      plotlist = plots,
      nrow = nrow,
      ncol = ncol,
      byrow = byrow
    )
  )
}

benchmark_collect_summary_data <- function(
  srt = NULL,
  data = NULL,
  tool_name = NULL,
  metrics = NULL
) {
  if (!is.null(data)) {
    summary_df <- as.data.frame(data, stringsAsFactors = FALSE)
  } else {
    tool_name <- tool_name %||% names(srt@tools)
    tool_name <- unique(as.character(tool_name))
    if (length(tool_name) == 0) {
      log_message(
        "No benchmark tool entries were available",
        message_type = "error"
      )
    }
    summary_list <- lapply(tool_name, function(tool_nm) {
      if (!tool_nm %in% names(srt@tools)) {
        return(NULL)
      }
      tool_res <- srt@tools[[tool_nm]]
      summary_df <- NULL
      if (
        is.list(tool_res) &&
          "summary" %in% names(tool_res) &&
          is.data.frame(tool_res[["summary"]])
      ) {
        summary_df <- as.data.frame(tool_res[["summary"]], stringsAsFactors = FALSE)
      }
      if (
        is.null(summary_df) &&
          is.list(tool_res) &&
          "metrics" %in% names(tool_res) &&
          is.list(tool_res[["metrics"]]) &&
          is.data.frame(tool_res[["metrics"]][["summary"]])
      ) {
        summary_df <- as.data.frame(tool_res[["metrics"]][["summary"]], stringsAsFactors = FALSE)
      }
      if (is.null(summary_df)) {
        return(NULL)
      }
      summary_df$tool_name <- tool_nm
      if (!"method" %in% colnames(summary_df)) {
        summary_df$method <- tool_res[["integration_method"]] %||%
          tool_res[["method"]] %||%
          tool_res[["label_method"]] %||%
          gsub("_metrics$", "", tool_nm)
      }
      if (!"workflow" %in% colnames(summary_df)) {
        summary_df$workflow <- if (grepl("LabelTransfer", tool_nm, ignore.case = TRUE)) {
          "LabelTransfer"
        } else if (grepl("ReferenceMapping", tool_nm, ignore.case = TRUE)) {
          "ReferenceMapping"
        } else if (grepl("integration", tool_nm, ignore.case = TRUE) || grepl("_metrics$", tool_nm)) {
          "Integration"
        } else {
          tool_nm
        }
      }
      summary_df
    })
    summary_list <- Filter(Negate(is.null), summary_list)
    if (length(summary_list) == 0) {
      log_message(
        "No summary benchmark data were found in {.arg srt@tools}",
        message_type = "error"
      )
    }
    summary_df <- do.call(rbind, summary_list)
  }

  if (!all(c("metric", "value") %in% colnames(summary_df))) {
    log_message(
      "{.arg data} must contain at least {.val metric} and {.val value} columns",
      message_type = "error"
    )
  }
  if (!"method" %in% colnames(summary_df)) {
    summary_df$method <- "Method1"
  }
  if (!"workflow" %in% colnames(summary_df)) {
    summary_df$workflow <- "Benchmark"
  }
  summary_df$metric <- as.character(summary_df$metric)
  summary_df$method <- as.character(summary_df$method)
  summary_df$workflow <- as.character(summary_df$workflow)
  summary_df$value <- as.numeric(summary_df$value)
  summary_df <- summary_df[!is.na(summary_df$value), , drop = FALSE]

  if (!is.null(metrics)) {
    metrics <- unique(as.character(metrics))
    summary_df <- summary_df[summary_df$metric %in% metrics, , drop = FALSE]
  }
  if (nrow(summary_df) == 0) {
    log_message(
      "No valid summary benchmark rows were available after filtering",
      message_type = "error"
    )
  }

  metric_levels <- metrics %||% unique(summary_df$metric)
  metric_levels <- metric_levels[metric_levels %in% summary_df$metric]
  summary_df$metric <- factor(summary_df$metric, levels = metric_levels)

  method_df <- stats::aggregate(
    value ~ method,
    data = summary_df,
    FUN = mean
  )
  method_df <- method_df[order(method_df$value, decreasing = TRUE), , drop = FALSE]
  summary_df$method <- factor(summary_df$method, levels = method_df$method)
  summary_df$workflow <- factor(summary_df$workflow, levels = unique(summary_df$workflow))
  summary_df
}

benchmark_summary_barplot <- function(
  summary_df,
  palette = "Chinese",
  palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list()
) {
  fill_cols <- palette_colors(
    levels(summary_df$method),
    palette = palette,
    palcolor = palcolor
  )

  p <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(
      x = .data[["method"]],
      y = .data[["value"]],
      fill = .data[["method"]]
    )
  ) +
    ggplot2::geom_col(
      width = 0.7,
      color = "black",
      linewidth = 0.2
    ) +
    ggplot2::scale_fill_manual(values = fill_cols, drop = FALSE) +
    ggplot2::labs(
      x = "",
      y = "Score"
    ) +
    do.call(theme_use, theme_args) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
    )

  if (length(levels(summary_df$workflow)) > 1) {
    p <- p + ggplot2::facet_grid(metric ~ workflow, scales = "free_y")
  } else {
    p <- p + ggplot2::facet_wrap(~ metric, scales = "free_y")
  }
  p
}

benchmark_summary_funkyheatmap <- function(
  summary_df,
  metrics = NULL,
  verbose = TRUE
) {
  check_r("funkyheatmap", verbose = FALSE)
  if (!requireNamespace("funkyheatmap", quietly = TRUE)) {
    log_message(
      "{.pkg funkyheatmap} is required for {.arg plot_type = 'funkyheatmap'}",
      message_type = "error",
      verbose = verbose
    )
  }

  metrics_use <- metrics %||% levels(summary_df$metric)
  metrics_use <- metrics_use[metrics_use %in% as.character(summary_df$metric)]
  wide_df <- stats::reshape(
    as.data.frame(summary_df[, c("method", "metric", "value")], stringsAsFactors = FALSE),
    idvar = "method",
    timevar = "metric",
    direction = "wide"
  )
  colnames(wide_df) <- sub("^value\\.", "", colnames(wide_df))
  keep_cols <- c("method", metrics_use)
  keep_cols <- keep_cols[keep_cols %in% colnames(wide_df)]
  wide_df <- wide_df[, keep_cols, drop = FALSE]

  column_info <- data.frame(
    id = c("method", metrics_use),
    group = c("", rep("Benchmark", length(metrics_use))),
    name = c("", metrics_use),
    geom = c("text", rep("funkyrect", length(metrics_use))),
    palette = c(NA_character_, rep("score_palette", length(metrics_use))),
    stringsAsFactors = FALSE
  )
  column_info$options <- I(c(
    list(list(hjust = 0, width = 4)),
    rep(list(list(width = 1.8, legend = TRUE)), length(metrics_use))
  ))
  column_info$options <- as.list(column_info$options)

  wide_df <- as.data.frame(wide_df, stringsAsFactors = FALSE)

  tryCatch(
    {
      funkyheatmap::funky_heatmap(
        data = wide_df,
        column_info = column_info,
        palettes = list(
          score_palette = grDevices::colorRampPalette(
            c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B")
          )(100)
        ),
        position_args = funkyheatmap::position_arguments(expand_xmax = 4)
      )
    },
    error = function(error) {
      log_message(
        "Failed to render {.pkg funkyheatmap}: {.val {conditionMessage(error)}}",
        message_type = "warning",
        verbose = verbose
      )
      benchmark_summary_barplot(summary_df = summary_df)
    }
  )
}

benchmark_clean_feature_labels <- function(features) {
  clean_method_label <- function(x) {
    x <- gsub("_+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    raw_index <- grepl(
      "pca(UMAP2D|TSNE2D|DM2D|PHATE2D|PACMAP2D|TRIMAP2D|LARGEVIS2D|FR2D)?$",
      x,
      ignore.case = FALSE
    )
    x <- sub("UMAP2D$", "", x, ignore.case = TRUE)
    x <- sub("TSNE2D$", "", x, ignore.case = TRUE)
    x <- sub("DM2D$", "", x, ignore.case = TRUE)
    x <- sub("PHATE2D$", "", x, ignore.case = TRUE)
    x <- sub("PACMAP2D$", "", x, ignore.case = TRUE)
    x <- sub("TRIMAP2D$", "", x, ignore.case = TRUE)
    x <- sub("LARGEVIS2D$", "", x, ignore.case = TRUE)
    x <- sub("FR2D$", "", x, ignore.case = TRUE)
    x <- sub("(?<![A-Za-z])pca$", "", x, perl = TRUE)
    x <- sub("(?<![A-Za-z])PCA$", "", x, perl = TRUE)
    x <- sub("lsi$", "", x, ignore.case = FALSE)
    x <- sub("LSI$", "", x, ignore.case = FALSE)
    x <- gsub("^_+|_+$", "", x)
    x[raw_index | grepl("^pca$", x, ignore.case = TRUE)] <- "Raw"
    x[nchar(x) == 0] <- NA_character_
    x
  }

  common_suffix <- function(x) {
    if (length(x) == 0) {
      return("")
    }
    rev_split <- lapply(x, function(val) {
      strsplit(
        paste(rev(strsplit(val, "")[[1]]), collapse = ""),
        ""
      )[[1]]
    })
    min_len <- min(vapply(rev_split, length, integer(1)))
    chars <- character(0)
    for (i in seq_len(min_len)) {
      current <- vapply(rev_split, `[`, character(1), i)
      if (length(unique(current)) != 1) {
        break
      }
      chars <- c(chars, current[1])
    }
    if (length(chars) == 0) {
      return("")
    }
    paste(rev(chars), collapse = "")
  }

  feature_suffix <- common_suffix(features)
  feature_labels <- features
  if (nzchar(feature_suffix)) {
    feature_labels <- sub(
      paste0(gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", feature_suffix), "$"),
      "",
      features
    )
  }
  feature_labels <- gsub("_+$", "", feature_labels)
  feature_labels <- clean_method_label(feature_labels)
  feature_labels[is.na(feature_labels)] <- features[is.na(feature_labels)]
  feature_labels[!nzchar(feature_labels)] <- features[!nzchar(feature_labels)]
  stats::setNames(feature_labels, features)
}

benchmark_feature_boxplot <- function(
  srt,
  features,
  palette = "Chinese",
  palcolor = NULL,
  boxplot_jitter = FALSE,
  theme_use = "theme_scop",
  theme_args = list(),
  verbose = TRUE
) {
  plot_list <- lapply(features, function(feature) {
    data.frame(
      feature = feature,
      score = srt@meta.data[[feature]],
      stringsAsFactors = FALSE
    )
  })
  plot_df <- do.call(rbind, plot_list)
  plot_df <- plot_df[!is.na(plot_df$score), , drop = FALSE]

  if (nrow(plot_df) == 0) {
    log_message(
      "No valid observations available for benchmark boxplot",
      message_type = "warning",
      verbose = verbose
    )
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  label_map <- benchmark_clean_feature_labels(features)
  plot_df$feature_label <- label_map[plot_df$feature]
  mean_df <- stats::aggregate(
    score ~ feature_label,
    data = plot_df,
    FUN = mean
  )
  mean_df <- mean_df[order(mean_df$score, decreasing = FALSE), , drop = FALSE]
  plot_df$feature_label <- factor(
    plot_df$feature_label,
    levels = mean_df$feature_label
  )
  fill_cols <- palette_colors(
    levels(plot_df$feature_label),
    palette = palette,
    palcolor = palcolor
  )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data[["feature_label"]],
      y = .data[["score"]],
      fill = .data[["feature_label"]]
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.5,
      outlier.shape = NA,
      linewidth = 0.4,
      alpha = 0.9
    )

  if (isTRUE(boxplot_jitter)) {
    p <- p +
      ggplot2::geom_jitter(
        width = 0.12,
        size = 0.25,
        alpha = 0.12,
        color = "grey50"
      )
  }

  if (
    length(levels(plot_df$feature_label)) >= 2 &&
      "Raw" %in% levels(plot_df$feature_label)
  ) {
    check_r("ggpubr", verbose = FALSE)
    comparisons <- lapply(
      setdiff(levels(plot_df$feature_label), "Raw"),
      function(method) c("Raw", method)
    )
    y_top <- max(plot_df$score, na.rm = TRUE) * 1.2
    p <- p +
      ggpubr::stat_compare_means(
        comparisons = comparisons,
        method = "wilcox.test",
        label = "p.signif",
        step.increase = 0.08,
        tip.length = 0.03,
        size = 3
      ) +
      ggplot2::coord_cartesian(ylim = c(0, y_top))
  }

  p +
    ggplot2::scale_fill_manual(values = fill_cols, drop = FALSE) +
    ggplot2::labs(
      x = "",
      y = "Score"
    ) +
    do.call(theme_use, theme_args) +
    ggplot2::theme(
      aspect.ratio = 1,
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "none",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}
