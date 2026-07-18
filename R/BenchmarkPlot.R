#' @title Plot benchmark metrics
#'
#' @description
#' Visualize benchmark results stored in a `Seurat` object, a summary
#' `data.frame`, or a `scop_benchmark` object created by [RunBenchmark()].
#' Spatial benchmark results default to a publication-oriented overview that
#' pairs clustering quality with runtime and peak-memory efficiency. Per-cell
#' metrics such as LISI remain available as feature plots and boxplots.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams CellDimPlot
#' @param srt A `Seurat` object.
#' @param data Optional `scop_benchmark` object or summary benchmark
#' `data.frame` containing at least `metric` and `value`, and optionally
#' `method`, `workflow`, and `direction`.
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
#' @param plot_type Plot type. `"overview"`, `"quality"`, `"efficiency"`, and
#' `"heatmap"` consume `scop_benchmark` results. Existing `"feature"`,
#' `"boxplot"`, `"bar"`, and `"funkyheatmap"` modes remain supported.
#' @param sort_by Method ordering for spatial benchmark plots. `"quality"`
#' sorts by the mean selected quality metric; other choices sort by method,
#' runtime, or peak memory.
#' @param show_values Whether to print raw metric values on quality and heatmap
#' panels.
#' @param show_status Whether the overview should add a status strip for
#' failed, unavailable, or timed-out methods.
#' @param resource_scale Resource-axis transformation. `"auto"` independently
#' uses log10 for runtime or memory when positive values span at least tenfold.
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
#' metrics_df <- data.frame(
#'   method = c("Raw", "Raw", "Harmony", "Harmony"),
#'   metric = c("batch_ASW_mixing", "celltype_ASW", "batch_ASW_mixing", "celltype_ASW"),
#'   value = c(0.42, 0.71, 0.68, 0.66)
#' )
#' BenchmarkPlot(
#'   data = metrics_df,
#'   plot_type = "bar"
#' )
#'
#' data("pbmcmultiome_sub", package = "scop")
#' pbmcmultiome_sub[["MethodA_batch_LISI"]] <-
#'   seq_len(ncol(pbmcmultiome_sub)) / ncol(pbmcmultiome_sub)
#' pbmcmultiome_sub[["MethodB_batch_LISI"]] <-
#'   rev(pbmcmultiome_sub[["MethodA_batch_LISI", drop = TRUE]])
#' BenchmarkPlot(
#'   pbmcmultiome_sub,
#'   features = c("MethodA_batch_LISI", "MethodB_batch_LISI"),
#'   plot_type = "boxplot"
#' )
#'
#' @export
BenchmarkPlot <- function(
  srt = NULL,
  data = NULL,
  features = NULL,
  metrics = NULL,
  tool_name = NULL,
  reduction = NULL,
  plot_type = c(
    "auto", "overview", "quality", "efficiency", "heatmap",
    "feature", "boxplot", "bar", "funkyheatmap"
  ),
  sort_by = c("quality", "method", "runtime", "memory"),
  show_values = TRUE,
  show_status = TRUE,
  resource_scale = c("auto", "linear", "log10"),
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
  sort_by <- match.arg(sort_by)
  resource_scale <- match.arg(resource_scale)
  benchmark_assert_flag(show_values, "show_values")
  benchmark_assert_flag(show_status, "show_status")

  benchmark_result <- if (inherits(data, "scop_benchmark")) data else NULL
  if (!is.null(benchmark_result)) {
    if (identical(plot_type, "auto")) plot_type <- "overview"
    if (plot_type %in% c("overview", "quality", "efficiency", "heatmap")) {
      return(benchmark_result_plot(
        result = benchmark_result,
        plot_type = plot_type,
        metrics = metrics,
        sort_by = sort_by,
        show_values = show_values,
        show_status = show_status,
        resource_scale = resource_scale,
        palette = palette,
        palcolor = palcolor,
        theme_use = theme_use,
        theme_args = theme_args
      ))
    }
    if (plot_type %in% c("feature", "boxplot")) {
      log_message(
        "{.arg plot_type = '{plot_type}'} requires per-cell metadata, not a {.cls scop_benchmark}",
        message_type = "error"
      )
    }
    data <- benchmark_result$metrics
  } else if (plot_type %in% c("overview", "quality", "efficiency", "heatmap")) {
    log_message(
      "{.arg plot_type = '{plot_type}'} requires a {.cls scop_benchmark} object",
      message_type = "error"
    )
  }

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
        isTRUE(check_r("funkyheatmap", verbose = FALSE)) &&
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

benchmark_result_plot <- function(
  result,
  plot_type,
  metrics,
  sort_by,
  show_values,
  show_status,
  resource_scale,
  palette,
  palcolor,
  theme_use,
  theme_args
) {
  quality_metrics <- benchmark_plot_quality_metrics(result, metrics)
  method_levels <- benchmark_plot_method_levels(
    result = result,
    metrics = quality_metrics,
    sort_by = sort_by
  )
  if (identical(plot_type, "quality")) {
    return(benchmark_quality_plot(
      result, quality_metrics, method_levels, show_values,
      palette, palcolor, theme_use, theme_args
    ))
  }
  if (identical(plot_type, "efficiency")) {
    return(benchmark_efficiency_plot(
      result, quality_metrics, method_levels, resource_scale,
      theme_use, theme_args
    ))
  }
  if (identical(plot_type, "heatmap")) {
    return(benchmark_metric_heatmap(
      result, quality_metrics, method_levels, show_values,
      theme_use, theme_args
    ))
  }

  quality <- benchmark_quality_plot(
    result, quality_metrics, method_levels, show_values,
    palette, palcolor, theme_use, theme_args
  )
  efficiency <- benchmark_efficiency_plot(
    result, quality_metrics, method_levels, resource_scale,
    theme_use, theme_args
  )
  overview <- patchwork::wrap_plots(
    quality, efficiency,
    nrow = 1,
    widths = c(1.2, 1)
  )
  failed <- result$summary$status != "success"
  if (isTRUE(show_status) && any(failed)) {
    status <- benchmark_status_plot(result$summary[failed, , drop = FALSE])
    overview <- patchwork::wrap_plots(
      overview, status,
      ncol = 1,
      heights = c(4, max(0.7, 0.42 * sum(failed)))
    )
  }
  overview + patchwork::plot_annotation(
    title = "Spatial clustering benchmark",
    subtitle = "Agreement: higher is better  |  Runtime and peak memory: lower is better",
    tag_levels = "a",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      plot.subtitle = ggplot2::element_text(color = "#5F6368", size = 7.5),
      plot.tag = ggplot2::element_text(face = "bold", size = 8)
    )
  )
}

benchmark_plot_quality_metrics <- function(result, metrics) {
  available <- c("ARI", "NMI", "purity")
  selected <- benchmark_resolve_metrics(metrics %||% available)
  selected <- selected[selected %in% available]
  if (length(selected) == 0L) {
    log_message("No quality metrics are available for benchmark plotting", message_type = "error")
  }
  selected
}

benchmark_plot_method_levels <- function(result, metrics, sort_by) {
  summary <- result$summary
  quality <- rowMeans(summary[, metrics, drop = FALSE], na.rm = TRUE)
  quality[!is.finite(quality)] <- -Inf
  order_index <- switch(sort_by,
    quality = order(quality, decreasing = TRUE, na.last = TRUE),
    method = order(tolower(summary$method)),
    runtime = order(summary$runtime_s, na.last = TRUE),
    memory = order(summary$peak_memory_mb, na.last = TRUE)
  )
  summary$method[order_index]
}

benchmark_quality_plot <- function(
  result,
  metrics,
  method_levels,
  show_values,
  palette,
  palcolor,
  theme_use,
  theme_args
) {
  rows <- lapply(metrics, function(metric) {
    data.frame(
      method = result$summary$method,
      metric = metric,
      value = result$summary[[metric]],
      status = result$summary$status,
      stringsAsFactors = FALSE
    )
  })
  data <- do.call(rbind, rows)
  data$metric <- factor(data$metric, levels = metrics)
  method_breaks <- seq_along(rev(method_levels))
  names(method_breaks) <- rev(method_levels)
  method_labels <- benchmark_plot_method_labels(result, names(method_breaks))
  metric_offsets <- stats::setNames(
    seq(-0.22, 0.22, length.out = length(metrics)),
    metrics
  )
  data$method_position <- unname(method_breaks[data$method]) +
    unname(metric_offsets[as.character(data$metric)])
  values <- data$value[is.finite(data$value)]
  lower <- if (length(values) > 0L && min(values) < 0) min(-0.05, min(values) - 0.04) else 0
  metric_colors <- palette_colors(
    metrics,
    palette = palette,
    palcolor = palcolor
  )
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data[["value"]],
      y = .data[["method_position"]],
      color = .data[["metric"]]
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.35,
      color = "#D9DDE3"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        xend = .data[["value"]],
        yend = .data[["method_position"]]
      ),
      linewidth = 0.8,
      alpha = 0.22,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      size = 2.8,
      stroke = 0.25,
      na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(values = metric_colors, drop = FALSE) +
    ggplot2::scale_x_continuous(
      limits = c(lower, 1.12),
      breaks = scales::breaks_width(0.25),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      breaks = method_breaks,
      labels = benchmark_wrap_method(method_labels),
      expand = ggplot2::expansion(add = 0.55)
    ) +
    ggplot2::labs(
      title = "Clustering agreement",
      subtitle = paste0("n = ", max(result$summary$n_evaluated, na.rm = TRUE), " aligned spots"),
      x = "Agreement score (higher is better)",
      y = NULL,
      color = NULL
    ) +
    do.call(theme_use, theme_args) +
    benchmark_publication_theme() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "#EEF0F3", linewidth = 0.3),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "top",
      legend.justification = "left",
      legend.box.margin = ggplot2::margin(0, 0, 1, 0),
      plot.title = ggplot2::element_text(face = "bold", size = 9),
      plot.subtitle = ggplot2::element_text(color = "#6B7280", size = 7)
    )
  if (isTRUE(show_values)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = ifelse(is.finite(.data[["value"]]), sprintf("%.2f", .data[["value"]]), "")),
      hjust = -0.35,
      size = 2.4,
      show.legend = FALSE,
      na.rm = TRUE
    )
  }
  p + ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(5, 12, 5, 5))
}

benchmark_efficiency_plot <- function(
  result,
  metrics,
  method_levels,
  resource_scale,
  theme_use,
  theme_args
) {
  data <- result$summary
  data$quality <- rowMeans(data[, metrics, drop = FALSE], na.rm = TRUE)
  data$quality[!is.finite(data$quality)] <- NA_real_
  data <- data[
    data$status == "success" &
      is.finite(data$runtime_s) & data$runtime_s > 0 &
      is.finite(data$peak_memory_mb) & data$peak_memory_mb > 0,
    ,
    drop = FALSE
  ]
  if (nrow(data) == 0L) {
    return(benchmark_empty_panel(
      "Computational efficiency",
      "No successful run has both\nruntime and memory measurements"
    ))
  }
  memory_in_gb <- max(data$peak_memory_mb, na.rm = TRUE) >= 1024
  if (memory_in_gb) data$memory_display <- data$peak_memory_mb / 1024 else data$memory_display <- data$peak_memory_mb
  x_log <- benchmark_resource_log(data$runtime_s, resource_scale)
  y_log <- benchmark_resource_log(data$memory_display, resource_scale)
  data$method <- factor(data$method, levels = method_levels)
  data$method_label <- benchmark_wrap_method(
    benchmark_plot_method_labels(result, as.character(data$method)),
    width = 18L
  )
  data$resource_label <- paste0(
    data$method_label,
    "\n",
    benchmark_format_metric(data$runtime_s, "runtime_s"),
    " | ",
    benchmark_format_metric(data$peak_memory_mb, "peak_memory_mb")
  )
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data[["runtime_s"]],
      y = .data[["memory_display"]],
      fill = .data[["quality"]]
    )
  ) +
    ggplot2::geom_point(
      shape = 21,
      size = 4,
      stroke = 0.45,
      color = "#263238"
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = .data[["resource_label"]]),
      size = 2.45,
      color = "#263238",
      lineheight = 0.95,
      box.padding = 0.45,
      point.padding = 0.25,
      min.segment.length = 0,
      segment.color = "#AAB2BD",
      segment.size = 0.25,
      show.legend = FALSE,
      max.overlaps = Inf,
      seed = 1
    ) +
    ggplot2::scale_fill_gradientn(
      colors = c("#DDE7F2", "#7FA9D1", "#174A7E"),
      limits = c(0, 1),
      oob = scales::squish,
      na.value = "#D6D8DC"
    ) +
    ggplot2::labs(
      title = "Computational efficiency",
      subtitle = "Lower-left is better; labels show resource use",
      x = if (x_log) "Runtime (s, log10)" else "Runtime (s)",
      y = paste0("Peak memory (", if (memory_in_gb) "GB" else "MB", if (y_log) ", log10" else "", ")"),
      fill = "Mean quality"
    ) +
    do.call(theme_use, theme_args) +
    benchmark_publication_theme() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "left",
      plot.title = ggplot2::element_text(face = "bold", size = 9),
      plot.subtitle = ggplot2::element_text(color = "#6B7280", size = 7)
    )
  if (x_log) {
    p <- p + ggplot2::scale_x_log10(
      labels = scales::label_number(),
      expand = ggplot2::expansion(mult = c(0.08, 0.18))
    )
  } else {
    p <- p + ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.08, 0.18))
    )
  }
  if (y_log) {
    p <- p + ggplot2::scale_y_log10(
      labels = scales::label_number(),
      expand = ggplot2::expansion(mult = c(0.10, 0.18))
    )
  } else {
    p <- p + ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.10, 0.18))
    )
  }
  p + ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(5, 8, 5, 5))
}

benchmark_wrap_method <- function(x, width = 22L) {
  vapply(
    x,
    function(label) paste(strwrap(label, width = width), collapse = "\n"),
    character(1),
    USE.NAMES = FALSE
  )
}

benchmark_plot_method_labels <- function(result, methods) {
  if (!"tier" %in% colnames(result$summary)) return(methods)
  tiers <- stats::setNames(as.character(result$summary$tier), result$summary$method)
  ifelse(tiers[methods] == "legacy", paste0(methods, " [legacy]"), methods)
}

benchmark_publication_theme <- function() {
  ggplot2::theme(
    text = ggplot2::element_text(color = "#263238"),
    axis.text = ggplot2::element_text(size = 7, color = "#263238"),
    axis.title = ggplot2::element_text(size = 7.5, color = "#263238"),
    axis.line = ggplot2::element_line(linewidth = 0.35, color = "#343A40"),
    panel.border = ggplot2::element_blank(),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    axis.ticks = ggplot2::element_line(linewidth = 0.3, color = "#343A40"),
    axis.ticks.length = grid::unit(1.3, "mm"),
    legend.text = ggplot2::element_text(size = 7),
    legend.title = ggplot2::element_text(size = 7),
    legend.key.height = grid::unit(3.2, "mm"),
    legend.key.width = grid::unit(4.2, "mm"),
    legend.background = ggplot2::element_blank(),
    legend.box.background = ggplot2::element_blank()
  )
}

benchmark_resource_log <- function(values, resource_scale) {
  if (identical(resource_scale, "log10")) return(TRUE)
  if (identical(resource_scale, "linear")) return(FALSE)
  values <- values[is.finite(values) & values > 0]
  length(values) >= 2L && max(values) / min(values) >= 10
}

benchmark_metric_heatmap <- function(
  result,
  quality_metrics,
  method_levels,
  show_values,
  theme_use,
  theme_args
) {
  metrics <- c(quality_metrics, "runtime_s", "peak_memory_mb")
  directions <- stats::setNames(
    c(rep("higher", length(quality_metrics)), "lower", "lower"),
    metrics
  )
  rows <- lapply(metrics, function(metric) {
    values <- result$summary[[metric]]
    score <- benchmark_normalize_metric(values, directions[[metric]])
    data.frame(
      method = result$summary$method,
      metric = metric,
      value = values,
      score = score,
      label = benchmark_format_metric(values, metric),
      stringsAsFactors = FALSE
    )
  })
  data <- do.call(rbind, rows)
  data$method <- factor(data$method, levels = rev(method_levels))
  labels <- c(
    ARI = "ARI", NMI = "NMI", purity = "Purity",
    runtime_s = "Runtime", peak_memory_mb = "Peak memory"
  )
  data$metric <- factor(data$metric, levels = metrics, labels = labels[metrics])
  p <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["metric"]], y = .data[["method"]], fill = .data[["score"]]
  )) +
    ggplot2::geom_tile(color = "white", linewidth = 0.8) +
    ggplot2::scale_fill_gradientn(
      colors = c("#F4F7FA", "#B5CBE1", "#4F83B3", "#123B63"),
      limits = c(0, 1),
      na.value = "#F1F2F4"
    ) +
    ggplot2::scale_y_discrete(labels = function(x) {
      benchmark_wrap_method(benchmark_plot_method_labels(result, x))
    }) +
    ggplot2::labs(
      title = "Direction-aware benchmark matrix",
      subtitle = "Color is normalized within each metric; labels show raw values",
      x = NULL,
      y = NULL,
      fill = "Relative\nperformance"
    ) +
    do.call(theme_use, theme_args) +
    benchmark_publication_theme() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1),
      plot.title = ggplot2::element_text(face = "bold", size = 9),
      plot.subtitle = ggplot2::element_text(color = "#6B7280", size = 7)
    )
  if (isTRUE(show_values)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(
        label = .data[["label"]],
        color = ifelse(is.finite(.data[["score"]]) & .data[["score"]] > 0.68, "light", "dark")
      ),
      size = 2.5,
      show.legend = FALSE
    ) +
      ggplot2::scale_color_manual(values = c(light = "white", dark = "#263238"))
  }
  p
}

benchmark_normalize_metric <- function(values, direction) {
  finite <- is.finite(values)
  out <- rep(NA_real_, length(values))
  if (!any(finite)) return(out)
  limits <- range(values[finite])
  if (diff(limits) == 0) {
    out[finite] <- 0.5
  } else {
    out[finite] <- (values[finite] - limits[[1]]) / diff(limits)
  }
  if (identical(direction, "lower")) out[finite] <- 1 - out[finite]
  out
}

benchmark_format_metric <- function(values, metric) {
  if (metric %in% c("ARI", "NMI", "purity")) {
    return(ifelse(is.finite(values), sprintf("%.2f", values), "NA"))
  }
  if (identical(metric, "runtime_s")) {
    return(ifelse(is.finite(values), paste0(scales::number(values, accuracy = 0.1), " s"), "NA"))
  }
  ifelse(
    is.finite(values),
    ifelse(
      values >= 1024,
      paste0(scales::number(values / 1024, accuracy = 0.01), " GB"),
      paste0(scales::number(values, accuracy = 1), " MB")
    ),
    "NA"
  )
}

benchmark_status_plot <- function(data) {
  status_colors <- c(
    failed = "#C95D63",
    unavailable = "#D89A45",
    timeout = "#8B6FAE"
  )
  data$method <- factor(data$method, levels = rev(data$method))
  data$label <- paste0(
    data$method, "  |  ", data$status,
    ifelse(nzchar(data$error), paste0("  -  ", benchmark_trim_text(data$error, 72L)), "")
  )
  ggplot2::ggplot(data, ggplot2::aes(
    x = 1, y = .data[["method"]], fill = .data[["status"]]
  )) +
    ggplot2::geom_tile(width = 1, height = 0.78, color = "white", linewidth = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data[["label"]]),
      x = 0.53,
      hjust = 0,
      color = "white",
      size = 2.2
    ) +
    ggplot2::scale_fill_manual(values = status_colors, drop = FALSE) +
    ggplot2::coord_cartesian(xlim = c(0.5, 1.5), clip = "off") +
    ggplot2::labs(title = "Incomplete runs") +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", size = 8, hjust = 0),
      plot.margin = ggplot2::margin(1, 4, 1, 4)
    )
}

benchmark_trim_text <- function(x, width) {
  x <- gsub("[\r\n]+", " ", x)
  ifelse(nchar(x) > width, paste0(substr(x, 1L, width - 3L), "..."), x)
}

benchmark_empty_panel <- function(title, subtitle) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0.5, y = 0.5,
      label = subtitle,
      color = "#6B7280",
      size = 3,
      lineheight = 1.1
    ) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 9))
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
    p <- p + ggplot2::facet_wrap(~metric, scales = "free_y")
  }
  p
}

benchmark_summary_funkyheatmap <- function(
  summary_df,
  metrics = NULL,
  verbose = TRUE
) {
  check_r("funkyheatmap", verbose = FALSE)
  check_r("funkyheatmap", verbose = FALSE)

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
    return(ggplot2::ggplot() +
      ggplot2::theme_void())
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
