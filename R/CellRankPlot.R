#' @title Plot stored CellRank results
#'
#' @description
#' Plot CellRank outputs without rerunning a Python backend.
#'
#' @param srt A Seurat object returned by [RunCellRank].
#' @param plot_type One of `"fate"`, `"states"`, `"circular"`, `"drivers"`,
#' `"trends"`, `"clusters"`, `"projection"`, or `"random_walks"`.
#' @param lineage Lineage used for driver/trend plots.
#' @param reduction Reduction used for cell-space plots.
#' @param group.by Group column used for state/random-walk plots.
#' @param top_n Number of driver genes to display.
#' @param n_sims Number of random walks.
#' @param max_iter Maximum length of each random walk.
#' @param seed Random seed.
#' @param ... Arguments passed to the underlying SCOP plotting function.
#'
#' @return A ggplot object or a SCOP plot object.
#' @export
CellRankPlot <- function(
  srt,
  plot_type = c("fate", "states", "circular", "drivers", "trends", "clusters", "projection", "random_walks"),
  lineage = NULL,
  reduction = NULL,
  group.by = NULL,
  top_n = 20L,
  n_sims = 100L,
  max_iter = 500L,
  seed = 0L,
  ...
) {
  plot_type <- match.arg(plot_type)
  if (!inherits(srt, "Seurat")) log_message("{.arg srt} must be a Seurat object", message_type = "error")
  if (is.null(srt@tools$CellRank)) log_message("CellRank results are missing", message_type = "error")
  fate <- tryCatch(cellrank_fate_matrix(srt), error = function(e) NULL)

  if (plot_type == "fate") {
    reduction <- reduction %||% DefaultReduction(srt)
    if (is.null(fate)) log_message("CellRank fate probabilities are missing", message_type = "error")
    cols <- paste0("cellrank_fate_", make.names(colnames(fate), unique = TRUE))
    cols <- intersect(cols, colnames(srt@meta.data))
    if (!length(cols)) log_message("No CellRank fate metadata columns are available", message_type = "error")
    return(FeatureDimPlot(srt, features = cols, reduction = reduction, ...))
  }
  if (plot_type == "states") {
    reduction <- reduction %||% DefaultReduction(srt)
    return(CellDimPlot(srt, group.by = group.by %||% "cellrank_terminal_states", reduction = reduction, ...))
  }
  if (plot_type == "projection") {
    reduction <- reduction %||% DefaultReduction(srt)
    time_key <- if ("palantir_pseudotime" %in% colnames(srt@meta.data)) {
      "palantir_pseudotime"
    } else {
      "cellrank_pseudotime"
    }
    return(PseudotimeProjectionPlot(
      srt, reduction = reduction, time_key = time_key,
      plot_type = "stream", group.by = group.by, ...
    ))
  }
  if (plot_type %in% c("drivers", "trends", "clusters") && is.null(lineage)) {
    log_message("{.arg lineage} is required for this plot type", message_type = "error")
  }

  if (plot_type == "circular") {
    if (is.null(fate) || ncol(fate) < 2L) log_message("At least two fate lineages are required", message_type = "error")
    fate <- as.matrix(fate)
    fate <- fate / pmax(rowSums(fate), .Machine$double.eps)
    angles <- seq(0, 2 * pi, length.out = ncol(fate) + 1L)[-(ncol(fate) + 1L)]
    xy <- cbind(
      x = as.numeric(fate %*% cos(angles)),
      y = as.numeric(fate %*% sin(angles))
    )
    xy <- data.frame(xy, cell = rownames(fate), confidence = apply(fate, 1L, max))
    vertices <- data.frame(x = cos(angles), y = sin(angles), label = colnames(fate))
    return(
      ggplot2::ggplot(xy, ggplot2::aes(x, y, color = confidence)) +
        ggplot2::geom_point(alpha = 0.7) +
        ggplot2::geom_point(data = vertices, ggplot2::aes(x, y), inherit.aes = FALSE, size = 3) +
        ggplot2::geom_text(data = vertices, ggplot2::aes(x, y, label = label), inherit.aes = FALSE, vjust = -1) +
        ggplot2::coord_equal(
          xlim = c(-1.25, 1.25),
          ylim = c(-1.25, 1.25),
          clip = "off"
        ) +
        ggplot2::scale_color_viridis_c(
          name = "Fate confidence",
          limits = c(0, 1),
          oob = scales::squish
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
          plot.background = ggplot2::element_rect(fill = "white", colour = NA),
          panel.background = ggplot2::element_rect(fill = "white", colour = NA),
          legend.background = ggplot2::element_rect(fill = "white", colour = NA),
          plot.margin = ggplot2::margin(12, 36, 12, 36)
        )
    )
  }

  if (plot_type == "drivers") {
    drivers <- srt@tools$CellRank$lineage_drivers
    corr_col <- paste0(lineage, "_corr")
    if (is.null(drivers) || !corr_col %in% colnames(drivers)) log_message("Stored lineage drivers are missing", message_type = "error")
    corr_values <- as.numeric(drivers[, corr_col])
    keep <- is.finite(corr_values)
    tab <- data.frame(
      gene = rownames(drivers)[keep],
      correlation = corr_values[keep],
      stringsAsFactors = FALSE
    )
    tab <- head(tab[order(tab$correlation, decreasing = TRUE), , drop = FALSE], as.integer(top_n))
    tab$gene <- factor(tab$gene, levels = rev(tab$gene))
    return(ggplot2::ggplot(tab, ggplot2::aes(gene, correlation)) + ggplot2::geom_col() + ggplot2::coord_flip() + ggplot2::theme_bw() + ggplot2::labs(x = NULL, y = "Correlation"))
  }

  trend <- srt@tools$CellRank$trends[[lineage]]
  if (is.null(trend)) log_message("Run {.fn RunCellRankTrends} before plotting trends", message_type = "error")
  mat <- as.matrix(trend$trend_matrix)
  if (plot_type == "trends") {
    genes <- intersect(trend$heatmap_genes %||% rownames(mat), rownames(mat))
    long <- reshape2::melt(mat[genes, , drop = FALSE], varnames = c("gene", "time"), value.name = "value")
    return(ggplot2::ggplot(long, ggplot2::aes(time, gene, fill = value)) + ggplot2::geom_raster() + ggplot2::scale_fill_viridis_c() + ggplot2::theme_minimal() + ggplot2::labs(x = "Pseudotime", y = NULL))
  }
  clusters <- trend$cluster_table
  if (plot_type == "clusters") {
    selected <- clusters$gene
    long <- reshape2::melt(mat[selected, , drop = FALSE], varnames = c("gene", "time"), value.name = "value")
    long$cluster <- clusters$cluster[match(long$gene, clusters$gene)]
    return(ggplot2::ggplot(long, ggplot2::aes(time, value, group = gene, color = cluster)) + ggplot2::geom_line(alpha = 0.25) + ggplot2::stat_summary(ggplot2::aes(group = cluster), fun = mean, geom = "line", linewidth = 1.2) + ggplot2::theme_bw() + ggplot2::labs(x = "Pseudotime", y = "Normalized trend"))
  }

  transition <- srt@tools$CellRank$transition_matrix %||% srt@graphs[["cellrank_transition"]]
  if (is.null(transition)) log_message("Stored CellRank transition matrix is missing", message_type = "error")
  reduction <- reduction %||% DefaultReduction(srt)
  coords <- SeuratObject::Embeddings(srt[[reduction]])[, 1:2, drop = FALSE]
  set.seed(seed)
  starts <- seq_len(nrow(coords))
  if (!is.null(group.by) && group.by %in% colnames(srt@meta.data)) starts <- which(as.character(srt@meta.data[[group.by]]) == unique(as.character(srt@meta.data[[group.by]]))[1L])
  paths <- vector("list", as.integer(n_sims))
  for (sim in seq_len(as.integer(n_sims))) {
    current <- sample(starts, 1L)
    visited <- integer(as.integer(max_iter))
    visited[1L] <- current
    for (step in 2:as.integer(max_iter)) {
      probs <- as.numeric(transition[current, ])
      if (!any(probs > 0)) { visited <- visited[seq_len(step - 1L)]; break }
      current <- sample.int(nrow(coords), 1L, prob = probs / sum(probs))
      visited[step] <- current
    }
    paths[[sim]] <- data.frame(sim = sim, step = seq_along(visited), coords[visited, , drop = FALSE])
  }
  path_df <- do.call(rbind, paths)
  names(path_df)[3:4] <- c("x", "y")
  ggplot2::ggplot(path_df, ggplot2::aes(x, y, group = sim, color = sim)) + ggplot2::geom_path(alpha = 0.35) + ggplot2::coord_equal() + ggplot2::theme_void() + ggplot2::guides(color = "none")
}
