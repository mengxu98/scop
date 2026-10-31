#' @title Plot stored CellOracle results
#'
#' @description
#' Display official CellOracle displacements, grid flow, predicted expression
#' changes or a cluster-specific regulatory network without rerunning Python.
#' @md
#' @param object A Seurat object returned by [RunCellOracle()].
#' @param tool_name Result name in `object@tools`.
#' @param plot_type One of `"quiver"`, `"flow"`, `"effect"`, or `"network"`.
#' @param randomized Show the official randomized displacement reference in
#' quiver/flow plots. This is not an independent biological replicate.
#' @param arrow_scale Positive multiplier of stored displacements in embedding
#' units. This is a display parameter, not matplotlib's inverse `scale` argument.
#' @param pt.size Background cell point size.
#' @param pt.alpha Background cell opacity.
#' @param top_n Number of genes in effect plots or edges in network plots.
#' @param cluster Cluster name required for network plots when multiple networks
#' are stored. Networks are never silently pooled across clusters.
#' @param theme_use Theme function or function name.
#' @param theme_args Named list of arguments passed to the theme function.
#' @return A ggplot object. Flow plots apply the stored upstream `mass_filter`
#' (`TRUE` means excluded), without smoothing or resampling. Effect bars display
#' mean absolute `delta_X` across cells, not fold changes or p-values. Network
#' arrow widths reflect absolute upstream coefficients, not activation/repression.
#' @seealso [RunCellOracle()]
#' @export
CellOraclePlot <- function(
  object, tool_name = "CellOracle",
  plot_type = c("quiver", "flow", "effect", "network"),
  randomized = FALSE, arrow_scale = 1, pt.size = 0.5, pt.alpha = 0.4,
  top_n = 20L, cluster = NULL, theme_use = "theme_scop", theme_args = list()
) {
  plot_type <- match.arg(plot_type)
  if (!inherits(object, "Seurat") || !tool_name %in% names(object@tools)) {
    log_message("Cannot find stored CellOracle results", message_type = "error")
  }
  result <- object@tools[[tool_name]]
  if (!is.numeric(arrow_scale) || length(arrow_scale) != 1L ||
      !is.finite(arrow_scale) || arrow_scale <= 0 ||
      !is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) ||
      top_n < 1 || top_n != as.integer(top_n)) {
    log_message("{.arg arrow_scale} must be positive and {.arg top_n} a positive integer", message_type = "error")
  }
  if (!is.logical(randomized) || length(randomized) != 1L || is.na(randomized)) {
    log_message("{.arg randomized} must be TRUE or FALSE", message_type = "error")
  }
  if (plot_type == "network") {
    networks <- result$filtered_links
    cluster <- cluster %||% if (length(networks) == 1L) names(networks) else NULL
    if (length(cluster) != 1L || !cluster %in% names(networks)) {
      log_message("Select one {.arg cluster} from the saved CellOracle networks", message_type = "error")
    }
    edges <- networks[[cluster]]
    if (!nrow(edges)) log_message("The selected network has no filtered edges", message_type = "error")
    edges <- utils::head(edges[order(-edges$coef_abs), , drop = FALSE], top_n)
    genes <- unique(c(as.character(edges$source), as.character(edges$target)))
    theta <- 2 * pi * seq_along(genes) / length(genes) + pi / 4
    nodes <- data.frame(x = cos(theta), y = sin(theta), gene = genes, row.names = genes)
    edges$x <- nodes[edges$source, "x"]
    edges$y <- nodes[edges$source, "y"]
    edges$xend <- nodes[edges$target, "x"]
    edges$yend <- nodes[edges$target, "y"]
    # Keep directed edges separately: a net-flow matrix would cancel reciprocal regulation.
    plot <- ggplot2::ggplot(edges, ggplot2::aes(x = .data$x, y = .data$y)) +
      ggplot2::geom_curve(ggplot2::aes(xend = .data$xend, yend = .data$yend,
        linewidth = .data$coef_abs), curvature = 0.15,
        arrow = grid::arrow(length = grid::unit(0.1, "inches"))) +
      ggplot2::geom_point(data = nodes, size = 3, color = "#336699") +
      ggplot2::geom_text(data = nodes, ggplot2::aes(label = .data$gene), nudge_y = 0.1) +
      ggplot2::scale_linewidth_continuous(range = c(0.3, 1.2), name = "Absolute coefficient") +
      ggplot2::coord_equal(clip = "off") +
      ggplot2::labs(title = paste("CellOracle network:", cluster), x = NULL, y = NULL)
    return(plot + do.call(match.fun(theme_use), theme_args) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()))
  }
  if (plot_type == "effect") {
    data <- result$gene_effect
    data <- utils::head(data[order(-data$mean_abs_delta), , drop = FALSE], top_n)
    data$gene <- factor(data$gene, levels = rev(data$gene))
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$mean_abs_delta, y = .data$gene)) +
      ggplot2::geom_col(fill = "#336699") +
      ggplot2::labs(x = "Mean absolute predicted expression change", y = NULL)
  } else {
    cells <- result$cell_shift
    vectors <- if (plot_type == "quiver") cells else result$grid_flow[!result$grid_flow$filtered, , drop = FALSE]
    vectors$xend <- vectors$x + arrow_scale * vectors[[if (randomized) "dx_random" else "dx"]]
    vectors$yend <- vectors$y + arrow_scale * vectors[[if (randomized) "dy_random" else "dy"]]
    plot <- ggplot2::ggplot(cells, ggplot2::aes(x = .data$x, y = .data$y)) +
      ggplot2::geom_point(color = "grey65", size = pt.size, alpha = pt.alpha) +
      ggplot2::geom_segment(
        data = vectors, ggplot2::aes(xend = .data$xend, yend = .data$yend),
        arrow = grid::arrow(length = grid::unit(0.06, "inches")), linewidth = 0.25
      ) +
      ggplot2::coord_equal() +
      ggplot2::labs(x = "Embedding 1", y = "Embedding 2",
        title = if (randomized) "CellOracle randomized reference" else "CellOracle predicted shift")
  }
  plot + do.call(match.fun(theme_use), theme_args)
}
