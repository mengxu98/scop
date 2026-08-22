#' @title Plot integration benchmark results
#'
#' @description
#' Visualize a `Seurat` object returned by [RunIntegrationBenchmark()].
#' `"box"` compares per-cell LISI scores across methods, `"heatmap"` shows
#' scaled summary metrics, `"scatter"` plots the bio-vs-batch trade-off, and
#' `"umap"` draws batch and cell-type embeddings. `"auto"` returns all of these
#' plots.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams CellDimPlot
#' @param srt A `Seurat` object from [RunIntegrationBenchmark()]. `"box"` can
#'   also plot any metadata columns ending in `"_LISI"`.
#' @param plot_type Plot type, or `"auto"` for a named list of all views.
#' @param tool_name `srt@tools` entry created by [RunIntegrationBenchmark()].
#' @param metrics Optional metric names used in the heatmap.
#' @param boxplot_jitter Whether to overlay jittered points on the LISI boxes.
#'
#' @return A ggplot, patchwork, or a named list of plots when
#' `plot_type = "auto"`.
#'
#' @seealso [RunIntegrationBenchmark], [RunLISI]
#'
#' @export
#'
#' @examples
#' srt <- SeuratObject::CreateSeuratObject(
#'   counts = Matrix::Matrix(
#'     matrix(1:20, nrow = 2, dimnames = list(c("g1", "g2"), paste0("c", 1:10))),
#'     sparse = TRUE
#'   )
#' )
#' set.seed(1)
#' srt$Uncorrected_tech_LISI <- stats::runif(10, 1.0, 1.4)
#' srt$Harmony_tech_LISI <- stats::runif(10, 2.0, 2.8)
#' srt$Uncorrected_celltype_LISI <- stats::runif(10, 1.0, 1.2)
#' srt$Harmony_celltype_LISI <- stats::runif(10, 1.0, 1.3)
#' srt@tools$IntegrationBenchmark <- list(
#'   summary = data.frame(
#'     method = c("Uncorrected", "Harmony"),
#'     bio = c(0.80, 0.82),
#'     batch = c(0.20, 0.70),
#'     overall = c(0.56, 0.77),
#'     status = "success",
#'     stringsAsFactors = FALSE
#'   ),
#'   metrics = data.frame(
#'     method = rep(c("Uncorrected", "Harmony"), each = 2),
#'     metric = rep(c("iLISI", "cLISI"), 2),
#'     category = rep(c("batch", "bio"), 2),
#'     value = c(1.2, 1.1, 2.4, 1.2),
#'     scaled = c(0.1, 0.9, 0.7, 0.85),
#'     direction = "higher",
#'     stringsAsFactors = FALSE
#'   ),
#'   runs = data.frame(
#'     method = c("Uncorrected", "Harmony"),
#'     status = "success",
#'     umap = NA_character_,
#'     stringsAsFactors = FALSE
#'   ),
#'   batch = "tech",
#'   celltype = "celltype"
#' )
#' thisplot::print_colored_table(
#'   srt@tools$IntegrationBenchmark$summary,
#'   by = "row",
#'   palette = "Chinese"
#' )
#' IntegrationBenchmarkPlot(srt, plot_type = "box")
#' IntegrationBenchmarkPlot(srt, plot_type = "heatmap")
#' IntegrationBenchmarkPlot(srt, plot_type = "scatter")
IntegrationBenchmarkPlot <- function(
  srt,
  plot_type = c("auto", "box", "heatmap", "scatter", "umap"),
  tool_name = "IntegrationBenchmark",
  metrics = NULL,
  palette = "Chinese",
  palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  pt.size = NULL,
  pt.alpha = 1,
  boxplot_jitter = FALSE,
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  verbose = TRUE,
  ...
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  plot_type <- match.arg(
    plot_type,
    choices = c("auto", "box", "heatmap", "scatter", "umap"),
    several.ok = TRUE
  )
  if ("auto" %in% plot_type) {
    plot_type <- c("box", "heatmap", "scatter", "umap")
  }
  info <- srt@tools[[tool_name]]
  plots <- list()
  requested <- plot_type
  if ("box" %in% plot_type) {
    plots$box <- integration_benchmark_box(
      srt,
      info = info,
      palette = palette,
      palcolor = palcolor,
      theme_use = theme_use,
      theme_args = theme_args,
      boxplot_jitter = boxplot_jitter,
      verbose = verbose
    )
  }
  if ("heatmap" %in% plot_type) {
    if (!is.null(info$metrics)) {
      plots$heatmap <- integration_benchmark_heatmap(
        info,
        metrics = metrics,
        theme_use = theme_use,
        theme_args = theme_args
      )
    } else if (identical(requested, "heatmap")) {
      integration_benchmark_require_info(info, "heatmap")
    }
  }
  if ("scatter" %in% plot_type) {
    if (!is.null(info$summary)) {
      plots$scatter <- integration_benchmark_scatter(
        info,
        palette = palette,
        palcolor = palcolor,
        theme_use = theme_use,
        theme_args = theme_args
      )
    } else if (identical(requested, "scatter")) {
      integration_benchmark_require_info(info, "scatter")
    }
  }
  if ("umap" %in% plot_type) {
    plots$umap <- integration_benchmark_umap(
      srt,
      info = info,
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

integration_benchmark_success_methods <- function(info) {
  if (is.null(info) || is.null(info$runs)) {
    return(character())
  }
  ok <- info$runs$method[info$runs$status == "success"]
  if (!is.null(info$summary) && "overall" %in% colnames(info$summary)) {
    ord <- info$summary$method[order(info$summary$overall, decreasing = TRUE)]
    ok <- intersect(as.character(ord), as.character(ok))
  }
  unique(as.character(ok))
}

integration_benchmark_require_info <- function(info, plot_type) {
  if (is.null(info) || is.null(info$summary) || is.null(info$metrics)) {
    log_message(
      "{.arg plot_type} {.val {plot_type}} needs {.code srt@tools$IntegrationBenchmark} from {.fn RunIntegrationBenchmark}",
      message_type = "error"
    )
  }
  info
}

parse_lisi_feature <- function(feature) {
  hit <- regexec("^(.*)_([^_]+)_LISI$", feature)
  parts <- regmatches(feature, hit)[[1L]]
  if (length(parts) != 3L) {
    return(NULL)
  }
  list(method = parts[[2L]], label = parts[[3L]])
}

integration_lisi_long <- function(srt, info = NULL) {
  methods <- integration_benchmark_success_methods(info)
  labels <- unique(c(info$batch, info$celltype))
  labels <- labels[!is.na(labels) & nzchar(as.character(labels))]
  cols <- grep("_LISI$", colnames(srt@meta.data), value = TRUE)
  cols <- cols[!grepl("_LISI_mean$", cols)]
  rows <- list()
  for (col in cols) {
    parsed <- parse_lisi_feature(col)
    if (is.null(parsed)) {
      next
    }
    if (length(methods) > 0L && !parsed$method %in% methods) {
      next
    }
    if (length(labels) > 0L && !parsed$label %in% labels) {
      next
    }
    metric <- parsed$label
    if (!is.null(info$batch) && identical(parsed$label, info$batch)) {
      metric <- "iLISI (batch mixing)"
    } else if (!is.null(info$celltype) && identical(parsed$label, info$celltype)) {
      metric <- "cLISI (cell-type purity)"
    } else if (identical(parsed$label, "tech") || grepl("batch", parsed$label, ignore.case = TRUE)) {
      metric <- paste0("iLISI (", parsed$label, ")")
    } else {
      metric <- paste0("cLISI (", parsed$label, ")")
    }
    values <- as.numeric(srt[[col, drop = TRUE]])
    keep <- is.finite(values)
    if (!any(keep)) {
      next
    }
    rows[[length(rows) + 1L]] <- data.frame(
      method = parsed$method,
      metric = metric,
      value = values[keep],
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0L) {
    return(NULL)
  }
  out <- do.call(rbind, rows)
  method_levels <- integration_benchmark_success_methods(info)
  if (length(method_levels) == 0L) {
    method_levels <- unique(out$method)
  }
  out$method <- factor(out$method, levels = method_levels)
  out$metric <- factor(out$metric, levels = unique(out$metric))
  out
}

integration_benchmark_box <- function(
  srt,
  info,
  palette,
  palcolor,
  theme_use,
  theme_args,
  boxplot_jitter,
  verbose
) {
  df <- integration_lisi_long(srt, info)
  if (is.null(df) || nrow(df) == 0L) {
    log_message(
      "No per-cell LISI columns were found for the box plot",
      message_type = "warning",
      verbose = verbose
    )
    return(NULL)
  }
  cols <- palette_colors(
    levels(df$method),
    palette = palette,
    palcolor = palcolor
  )
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[["method"]],
      y = .data[["value"]],
      fill = .data[["method"]]
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.55,
      outlier.size = 0.4,
      outlier.alpha = 0.35,
      linewidth = 0.35,
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
  p +
    ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
    ggplot2::facet_wrap(~metric, scales = "free_y") +
    ggplot2::labs(
      title = "LISI by integration method",
      x = NULL,
      y = "LISI",
      fill = NULL
    ) +
    apply_plot_theme(theme_use, theme_args) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA)
    )
}

integration_benchmark_heatmap <- function(
  info,
  metrics,
  theme_use,
  theme_args
) {
  info <- integration_benchmark_require_info(info, "heatmap")
  df <- info$metrics
  methods <- integration_benchmark_success_methods(info)
  df <- df[df$method %in% methods, , drop = FALSE]
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
  df$method <- factor(df$method, levels = methods)
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
  info,
  palette,
  palcolor,
  theme_use,
  theme_args
) {
  info <- integration_benchmark_require_info(info, "scatter")
  df <- info$summary
  methods <- integration_benchmark_success_methods(info)
  df <- df[df$method %in% methods, , drop = FALSE]
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
  srt,
  info,
  palette,
  palcolor,
  theme_use,
  theme_args,
  pt.size,
  pt.alpha,
  verbose,
  ...
) {
  if (is.null(info) || is.null(info$runs)) {
    log_message(
      "UMAP plots need {.code srt@tools$IntegrationBenchmark} from {.fn RunIntegrationBenchmark}",
      message_type = "warning",
      verbose = verbose
    )
    return(NULL)
  }
  methods <- integration_benchmark_success_methods(info)
  group_cols <- unique(c(info$batch, info$celltype))
  group_cols <- group_cols[!is.na(group_cols) & group_cols %in% colnames(srt@meta.data)]
  if (length(group_cols) == 0L) {
    return(NULL)
  }
  plots <- list()
  for (method in methods) {
    umap <- info$runs$umap[info$runs$method == method]
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
