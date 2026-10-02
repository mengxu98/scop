#' @title Plot IREA response or polarization results
#'
#' @description Display enrichment effects and reference-cell significance.
#' @param object An `irea_result`, a Seurat object containing
#'   `object@tools$IREA`, or a named list of results for dotplots and heatmaps.
#' @param plot_type One of `"compass"`, `"radar"`, `"dotplot"`, or `"heatmap"`.
#' @param padjustCutoff Significance threshold used for compass colouring and the
#'   radar display. Radar scores are zero when no positive effect meets this
#'   threshold; otherwise all positive effects are divided by their maximum.
#' @param palette,palcolor Palette name and optional custom colours passed to
#'   [thisplot::palette_colors()]. For compass plots, colours map to Positive,
#'   Negative and Not significant, in that order.
#' @param theme_use Theme name, function or ggplot theme, as in other package plots.
#' @return An editable ggplot object; plotting does not change result statistics.
#' @md
#' @inherit RunIREA references
#' @examples
#' \dontrun{
#' reference <- PrepareDB(
#'   db = "IREA_Macrophage", species = "Homo_sapiens"
#' )[["Homo_sapiens"]][["IREA_Macrophage"]]
#' result <- RunIREA(c("ISG15", "IFIT3", "BST2"), reference = reference)
#' IREAPlot(result, palette = "Chinese")
#' }
#' @export
IREAPlot <- function(object, plot_type = c("compass", "radar", "dotplot", "heatmap"),
                     padjustCutoff = 0.05, palette = "Chinese", palcolor = NULL,
                     theme_use = "theme_scop") {
  result <- object
  plot_theme <- apply_plot_theme(theme_use)
  plot_type <- match.arg(plot_type)
  if (inherits(result, "Seurat")) result <- result@tools$IREA
  multi <- is.list(result) && !inherits(result, "irea_result")
  if (multi) {
    if (plot_type %in% c("compass", "radar")) {
      log_message("Compass and radar plots require one IREA result.", message_type = "error")
    }
    if (!length(result) || !all(vapply(result, inherits, logical(1), "irea_result"))) {
      log_message("Every list element must be an IREA result.", message_type = "error")
    }
    if (length(unique(vapply(result, function(x) x$parameters$analysis, character(1)))) != 1L) {
      log_message("Comparison results must use the same analysis.", message_type = "error")
    }
    for (field in c("method", "mode", "species")) {
      values <- vapply(result, function(x) {
        value <- x$parameters[[field]]
        if (is.null(value)) NA_character_ else as.character(value)
      }, character(1))
      if (anyNA(values) || length(unique(values)) != 1L) {
        log_message("Comparison results must use the same ", field, ".", message_type = "error")
      }
    }
    groups <- names(result)
    if (is.null(groups) || any(!nzchar(groups))) groups <- paste0("Result ", seq_along(result))
    dat <- do.call(rbind, Map(function(x, name) transform(x$table, group = name), result, groups))
  } else {
    if (!inherits(result, "irea_result")) log_message("result must be an IREA result.", message_type = "error")
    dat <- result$table
    dat$group <- "Result"
  }
  if (!is.numeric(padjustCutoff) || length(padjustCutoff) != 1L || is.na(padjustCutoff) ||
    padjustCutoff <= 0 || padjustCutoff > 1) {
    log_message("Invalid padjustCutoff.", message_type = "error")
  }
  dat$significant <- !is.na(dat$fdr) & dat$fdr < padjustCutoff
  caption <- "Experimental IREA; numerical equivalence to the portal is not established."
  if (plot_type == "compass") {
    if (result$parameters$analysis != "cytokine_response") {
      log_message("Compass plot requires a cytokine response result.", message_type = "error")
    }
    dat$direction <- ifelse(!dat$significant, "Not significant",
      ifelse(dat$effect >= 0, "Positive", "Negative")
    )
    dat$term <- factor(dat$term, levels = dat$term[order(dat$effect)])
    return(ggplot2::ggplot(dat, ggplot2::aes(x = .data$term, y = abs(.data$effect), fill = .data$direction)) +
      ggplot2::geom_col(width = 1) +
      ggplot2::coord_polar() +
      ggplot2::scale_fill_manual(values = thisplot::palette_colors(c("Positive", "Negative", "Not significant"),
        palette = palette, palcolor = palcolor
      )) +
      plot_theme +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(size = 6),
        panel.border = ggplot2::element_blank()
      ) +
      ggplot2::labs(x = NULL, y = "Absolute enrichment effect", fill = "Response", caption = caption))
  }
  if (plot_type == "radar") {
    if (result$parameters$analysis != "cell_polarization") {
      log_message("Radar plot requires a cell polarization result.", message_type = "error")
    }
    dat$radar_score <- irea_radar_score(dat$effect, dat$fdr, padjustCutoff)
    colors <- thisplot::palette_colors(c("Response", "Grid"), palette = palette, palcolor = palcolor)
    theta <- pi / 2 - 2 * pi * (seq_len(nrow(dat)) - 1) / nrow(dat)
    dat$x <- dat$radar_score * cos(theta)
    dat$y <- dat$radar_score * sin(theta)
    axes <- data.frame(x = cos(theta), y = sin(theta), term = dat$term)
    circles <- do.call(rbind, lapply(c(0.25, 0.5, 0.75, 1), function(radius) {
      t <- seq(0, 2 * pi, length.out = 181)
      data.frame(x = radius * cos(t), y = radius * sin(t), radius = factor(radius))
    }))
    return(ggplot2::ggplot() +
      ggplot2::geom_path(data = circles, ggplot2::aes(x = .data$x, y = .data$y, group = .data$radius), color = colors[["Grid"]], alpha = 0.3) +
      ggplot2::geom_segment(data = axes, ggplot2::aes(x = 0, y = 0, xend = .data$x, yend = .data$y), color = colors[["Grid"]], alpha = 0.3) +
      ggplot2::geom_polygon(
        data = dat, ggplot2::aes(x = .data$x, y = .data$y), fill = colors[["Response"]],
        alpha = 0.25, color = colors[["Response"]]
      ) +
      ggplot2::geom_point(data = dat, ggplot2::aes(x = .data$x, y = .data$y), color = colors[["Response"]]) +
      ggplot2::geom_text(data = axes, ggplot2::aes(x = 1.13 * .data$x, y = 1.13 * .data$y, label = .data$term), size = 3) +
      ggplot2::coord_fixed(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25)) +
      plot_theme +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      ) +
      ggplot2::labs(title = "Cell polarization", subtitle = "Normalized score (0 to 1)", caption = caption))
  }
  dat$term <- factor(dat$term, levels = unique(dat$term[order(dat$effect)]))
  if (plot_type == "dotplot") {
    return(ggplot2::ggplot(dat, ggplot2::aes(x = .data$group, y = .data$term)) +
      ggplot2::geom_point(ggplot2::aes(
        size = -log10(pmax(.data$fdr, .Machine$double.xmin)),
        color = .data$effect
      )) +
      ggplot2::scale_color_gradientn(colors = thisplot::palette_colors(
        type = "continuous", palette = palette, palcolor = palcolor
      )) +
      plot_theme +
      ggplot2::labs(x = NULL, y = NULL, size = "-log10(FDR)", color = "Effect", caption = caption))
  }
  ggplot2::ggplot(dat, ggplot2::aes(x = .data$group, y = .data$term, fill = .data$effect)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colors = thisplot::palette_colors(
      type = "continuous", palette = palette, palcolor = palcolor
    )) +
    plot_theme +
    ggplot2::labs(x = NULL, y = NULL, fill = "Enrichment effect", caption = caption)
}
