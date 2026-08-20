#' @title Plot integration benchmark results
#'
#' @description
#' Visualize an [RunIntegrationBenchmark()] result. `"overview"` shows scIB-style
#' bio / batch / overall scores, `"heatmap"` shows scaled metrics, `"scatter"`
#' plots the bio-vs-batch trade-off, `"umap"` draws batch and cell-type
#' embeddings, and `"lisi"` maps iLISI/cLISI on UMAP. `"auto"` returns all of
#' these plots.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams CellDimPlot
#' @param result An `integration_benchmark_result` from
#'   [RunIntegrationBenchmark()].
#' @param plot_type Plot type, or `"auto"` for a named list of all views.
#' @param metrics Optional metric names used in the heatmap.
#'
#' @return A ggplot, patchwork, or a named list of plots when
#' `plot_type = "auto"`.
#'
#' @seealso [RunIntegrationBenchmark], [LISIPlot], [BenchmarkPlot]
#'
#' @export
#'
#' @examples
#' metrics <- data.frame(
#'   method = rep(c("Uncorrected", "Harmony"), each = 2),
#'   metric = rep(c("iLISI", "cLISI"), 2),
#'   category = rep(c("batch", "bio"), 2),
#'   value = c(1.2, 1.1, 2.4, 1.2),
#'   scaled = c(0.1, 0.9, 0.7, 0.85),
#'   direction = "higher",
#'   stringsAsFactors = FALSE
#' )
#' fake <- structure(
#'   list(
#'     summary = data.frame(
#'       method = c("Uncorrected", "Harmony"),
#'       bio = c(0.80, 0.82),
#'       batch = c(0.20, 0.70),
#'       overall = c(0.56, 0.77),
#'       status = "success",
#'       stringsAsFactors = FALSE
#'     ),
#'     metrics = metrics,
#'     runs = data.frame(
#'       method = c("Uncorrected", "Harmony"),
#'       status = "success",
#'       umap = NA_character_,
#'       stringsAsFactors = FALSE
#'     ),
#'     srt = NULL,
#'     batch = "tech",
#'     celltype = "celltype"
#'   ),
#'   class = c("integration_benchmark_result", "list")
#' )
#' IntegrationBenchmarkPlot(fake, plot_type = "overview")
#' IntegrationBenchmarkPlot(fake, plot_type = "heatmap")
#' IntegrationBenchmarkPlot(fake, plot_type = "scatter")
IntegrationBenchmarkPlot <- function(
  result,
  plot_type = c("auto", "overview", "heatmap", "scatter", "umap", "lisi"),
  metrics = NULL,
  palette = "Chinese",
  palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  pt.size = NULL,
  pt.alpha = 1,
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  verbose = TRUE,
  ...
) {
  if (!inherits(result, "integration_benchmark_result")) {
    log_message(
      "{.arg result} must be an {.cls integration_benchmark_result}",
      message_type = "error"
    )
  }
  plot_type <- match.arg(
    plot_type,
    choices = c("auto", "overview", "heatmap", "scatter", "umap", "lisi"),
    several.ok = TRUE
  )
  if ("auto" %in% plot_type) {
    plot_type <- c("overview", "heatmap", "scatter", "umap", "lisi")
  }
  plots <- list()
  if ("overview" %in% plot_type) {
    plots$overview <- integration_benchmark_overview(
      result,
      palette = palette,
      palcolor = palcolor,
      theme_use = theme_use,
      theme_args = theme_args
    )
  }
  if ("heatmap" %in% plot_type) {
    plots$heatmap <- integration_benchmark_heatmap(
      result,
      metrics = metrics,
      theme_use = theme_use,
      theme_args = theme_args
    )
  }
  if ("scatter" %in% plot_type) {
    plots$scatter <- integration_benchmark_scatter(
      result,
      palette = palette,
      palcolor = palcolor,
      theme_use = theme_use,
      theme_args = theme_args
    )
  }
  if ("umap" %in% plot_type) {
    plots$umap <- integration_benchmark_umap(
      result,
      palette = palette,
      palcolor = palcolor,
      theme_use = theme_use,
      theme_args = theme_args,
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      verbose = verbose,
      ...
    )
  }
  if ("lisi" %in% plot_type) {
    plots$lisi <- integration_benchmark_lisi(
      result,
      theme_use = theme_use,
      theme_args = theme_args,
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      verbose = verbose,
      ...
    )
  }
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) == 0L) {
    log_message(
      "No integration benchmark plots could be drawn",
      message_type = "error"
    )
  }
  if (length(plots) == 1L && isTRUE(combine)) {
    return(plots[[1L]])
  }
  plots
}

integration_benchmark_success_methods <- function(result) {
  ok <- result$runs$method[result$runs$status == "success"]
  if ("overall" %in% colnames(result$summary)) {
    ord <- result$summary$method[order(result$summary$overall, decreasing = TRUE)]
    ok <- intersect(as.character(ord), as.character(ok))
  }
  unique(as.character(ok))
}

integration_benchmark_overview <- function(
  result,
  palette,
  palcolor,
  theme_use,
  theme_args
) {
  df <- result$summary
  df <- df[df$status %in% c("success", NA, ""), , drop = FALSE]
  if (nrow(df) == 0L) {
    df <- result$summary
  }
  long <- rbind(
    data.frame(method = df$method, score = "Bio conservation", value = df$bio),
    data.frame(method = df$method, score = "Batch correction", value = df$batch),
    data.frame(method = df$method, score = "Overall", value = df$overall)
  )
  long$score <- factor(
    long$score,
    levels = c("Bio conservation", "Batch correction", "Overall")
  )
  method_levels <- as.character(df$method[order(df$overall, decreasing = TRUE)])
  long$method <- factor(long$method, levels = method_levels)
  cols <- palette_colors(
    levels(long$score),
    palette = palette,
    palcolor = palcolor
  )
  ggplot2::ggplot(
    long,
    ggplot2::aes(
      x = .data[["method"]],
      y = .data[["value"]],
      fill = .data[["score"]]
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.8),
      width = 0.72,
      color = "black",
      linewidth = 0.2
    ) +
    ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(
      title = "Integration scores",
      x = NULL,
      y = "Scaled score (higher is better)",
      fill = NULL
    ) +
    apply_plot_theme(theme_use, theme_args) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      legend.position = "top"
    )
}

integration_benchmark_heatmap <- function(
  result,
  metrics,
  theme_use,
  theme_args
) {
  df <- result$metrics
  df <- df[df$method %in% integration_benchmark_success_methods(result), , drop = FALSE]
  keep <- integration_metric_meta()$metric
  if (!is.null(metrics)) {
    keep <- intersect(as.character(metrics), unique(df$metric))
  }
  df <- df[df$metric %in% keep, , drop = FALSE]
  if (nrow(df) == 0L) {
    log_message(
      "No scaled integration metrics were available for the heatmap",
      message_type = "error"
    )
  }
  df$method <- factor(df$method, levels = integration_benchmark_success_methods(result))
  df$metric <- factor(df$metric, levels = keep[keep %in% df$metric])
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[["metric"]],
      y = .data[["method"]],
      fill = .data[["scaled"]]
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", .data[["scaled"]])),
      size = 3
    ) +
    ggplot2::scale_fill_gradientn(
      colours = palette_colors(palette = "RdBu"),
      limits = c(0, 1),
      name = "Scaled"
    ) +
    ggplot2::labs(
      title = "Scaled integration metrics",
      x = NULL,
      y = NULL
    ) +
    apply_plot_theme(theme_use, theme_args) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )
}

integration_benchmark_scatter <- function(
  result,
  palette,
  palcolor,
  theme_use,
  theme_args
) {
  df <- result$summary
  df <- df[df$method %in% integration_benchmark_success_methods(result), , drop = FALSE]
  cols <- palette_colors(
    as.character(df$method),
    palette = palette,
    palcolor = palcolor
  )
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[["batch"]],
      y = .data[["bio"]],
      color = .data[["method"]],
      label = .data[["method"]]
    )
  ) +
    ggplot2::geom_point(size = 3.5) +
    ggrepel::geom_text_repel(min.segment.length = 0.2, seed = 11, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = cols, drop = FALSE) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      title = "Bio conservation vs batch correction",
      x = "Batch correction (higher is better)",
      y = "Bio conservation (higher is better)",
      color = NULL
    ) +
    apply_plot_theme(theme_use, theme_args) +
    ggplot2::theme(legend.position = "none")
}

integration_benchmark_umap <- function(
  result,
  palette,
  palcolor,
  theme_use,
  theme_args,
  pt.size,
  pt.alpha,
  verbose,
  ...
) {
  srt <- result$srt
  if (!inherits(srt, "Seurat")) {
    log_message(
      "UMAP plots need {.code result$srt} from {.fn RunIntegrationBenchmark}",
      message_type = "warning",
      verbose = verbose
    )
    return(NULL)
  }
  methods <- integration_benchmark_success_methods(result)
  group_cols <- unique(c(result$batch, result$celltype))
  group_cols <- group_cols[!is.na(group_cols) & group_cols %in% colnames(srt@meta.data)]
  if (length(group_cols) == 0L) {
    return(NULL)
  }
  plots <- list()
  for (method in methods) {
    umap <- result$runs$umap[result$runs$method == method]
    umap <- umap[[1L]]
    if (is.na(umap) || !umap %in% SeuratObject::Reductions(srt)) {
      next
    }
    for (group in group_cols) {
      p <- CellDimPlot(
        srt,
        group.by = group,
        reduction = umap,
        title = paste(method, group),
        palette = palette,
        palcolor = palcolor,
        theme_use = theme_use,
        theme_args = theme_args,
        pt.size = pt.size,
        pt.alpha = pt.alpha,
        combine = FALSE,
        ...
      )
      if (is.list(p)) {
        p <- p[[1L]]
      }
      plots[[paste(method, group, sep = "__")]] <- p
    }
  }
  if (length(plots) == 0L) {
    return(NULL)
  }
  patchwork::wrap_plots(
    plots,
    ncol = length(group_cols)
  )
}

integration_benchmark_lisi <- function(
  result,
  theme_use,
  theme_args,
  pt.size,
  pt.alpha,
  verbose,
  ...
) {
  srt <- result$srt
  if (!inherits(srt, "Seurat")) {
    log_message(
      "LISI plots need {.code result$srt} from {.fn RunIntegrationBenchmark}",
      message_type = "warning",
      verbose = verbose
    )
    return(NULL)
  }
  methods <- integration_benchmark_success_methods(result)
  plots <- list()
  for (method in methods) {
    umap <- result$runs$umap[result$runs$method == method]
    umap <- umap[[1L]]
    if (is.na(umap) || !umap %in% SeuratObject::Reductions(srt)) {
      umap <- NULL
    }
    features <- grep(
      paste0("^", make.names(method), "_.*_LISI$"),
      colnames(srt@meta.data),
      value = TRUE
    )
    if (length(features) == 0L) {
      next
    }
    p <- FeatureDimPlot(
      srt,
      features = features,
      reduction = umap,
      theme_use = theme_use,
      theme_args = theme_args,
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      combine = TRUE,
      ncol = length(features),
      ...
    )
    plots[[method]] <- p
  }
  if (length(plots) == 0L) {
    return(NULL)
  }
  if (length(plots) == 1L) {
    return(plots[[1L]])
  }
  patchwork::wrap_plots(plots, ncol = 1)
}
