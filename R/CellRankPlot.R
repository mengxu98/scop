#' @title Plot stored CellRank results
#'
#' @description
#' Plot CellRank outputs without rerunning a Python backend.
#'
#' @param srt A Seurat object returned by [RunCellRank].
#' @param plot_type One of `"fate"`, `"states"`, `"circular"`, `"drivers"`,
#' `"trends"`, `"clusters"`, `"enrichment"`, `"projection"`, or
#' `"random_walks"`.
#' @param lineage Lineage used for driver/trend plots.
#' @param database Enrichment database used when `plot_type = "enrichment"`.
#' If `NULL`, the first stored non-empty database is used.
#' @param palette,palcolor Discrete SCOP palette for states, modules, and cell
#' groups.
#' @param feature_palette,feature_palcolor Continuous SCOP palette for fate,
#' trends, pseudotime, and enrichment strength.
#' @param theme_use,theme_args SCOP plot theme and arguments passed to it.
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
  plot_type = c("fate", "states", "circular", "drivers", "trends", "clusters", "enrichment", "projection", "random_walks"),
  lineage = NULL,
  database = NULL,
  reduction = NULL,
  group.by = NULL,
  top_n = 20L,
  n_sims = 100L,
  max_iter = 500L,
  seed = 0L,
  palette = "Chinese",
  palcolor = NULL,
  feature_palette = "Spectral",
  feature_palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  ...
) {
  plot_type <- match.arg(plot_type)
  if (!inherits(srt, "Seurat")) log_message("{.arg srt} must be a Seurat object", message_type = "error")
  if (is.null(srt@tools$CellRank)) log_message("CellRank results are missing", message_type = "error")
  fate <- tryCatch(cellrank_fate_matrix(srt), error = function(e) NULL)
  theme_fun <- resolve_plot_theme_use(theme_use)
  theme_layer <- do.call(theme_fun, theme_args)
  continuous_colors <- unname(palette_colors(
    n = 11L,
    palette = feature_palette,
    palcolor = feature_palcolor
  ))

  if (plot_type == "fate") {
    reduction <- reduction %||% DefaultReduction(srt)
    if (is.null(fate)) log_message("CellRank fate probabilities are missing", message_type = "error")
    cols <- paste0("cellrank_fate_", make.names(colnames(fate), unique = TRUE))
    cols <- intersect(cols, colnames(srt@meta.data))
    if (!length(cols)) log_message("No CellRank fate metadata columns are available", message_type = "error")
    return(FeatureDimPlot(
      srt,
      features = cols,
      reduction = reduction,
      palette = feature_palette,
      palcolor = feature_palcolor,
      theme_use = theme_use,
      theme_args = theme_args,
      ...
    ))
  }
  if (plot_type == "states") {
    reduction <- reduction %||% DefaultReduction(srt)
    return(CellDimPlot(
      srt,
      group.by = group.by %||% "cellrank_terminal_states",
      reduction = reduction,
      palette = palette,
      palcolor = palcolor,
      theme_use = theme_use,
      theme_args = theme_args,
      ...
    ))
  }
  if (plot_type == "projection") {
    reduction <- reduction %||% DefaultReduction(srt)
    transition <- srt@tools$CellRank$transition_matrix %||%
      srt@graphs[["cellrank_transition"]]
    if (is.null(transition)) {
      log_message("Stored CellRank transition matrix is missing", message_type = "error")
    }
    cells <- colnames(srt)
    coords <- SeuratObject::Embeddings(srt[[reduction]])[, 1:2, drop = FALSE]
    if (!is.null(rownames(transition)) && !is.null(colnames(transition))) {
      missing_cells <- union(
        setdiff(cells, rownames(transition)),
        setdiff(cells, colnames(transition))
      )
      if (length(missing_cells)) {
        log_message("Stored CellRank transition matrix is not aligned to current cells", message_type = "error")
      }
      transition <- transition[cells, cells, drop = FALSE]
    } else if (!identical(dim(transition), c(length(cells), length(cells)))) {
      log_message("Stored CellRank transition matrix has incompatible dimensions", message_type = "error")
    }
    next_xy <- as.matrix(transition %*% coords)
    delta <- next_xy - coords
    magnitude <- sqrt(rowSums(delta ^ 2))
    scale_ref <- as.numeric(stats::quantile(magnitude[magnitude > 0], 0.95, na.rm = TRUE))
    if (!is.finite(scale_ref) || scale_ref <= 0) scale_ref <- 1
    delta <- delta / scale_ref * 0.8
    n_arrows <- min(nrow(coords), 250L)
    arrow_idx <- unique(as.integer(round(seq(1, nrow(coords), length.out = n_arrows))))
    arrow_df <- data.frame(
      x = coords[arrow_idx, 1],
      y = coords[arrow_idx, 2],
      xend = coords[arrow_idx, 1] + delta[arrow_idx, 1],
      yend = coords[arrow_idx, 2] + delta[arrow_idx, 2]
    )
    cell_df <- data.frame(x = coords[, 1], y = coords[, 2])
    if (!is.null(group.by) && group.by %in% colnames(srt@meta.data)) {
      cell_df$group <- srt@meta.data[[group.by]]
      group_colors <- palette_colors(
        x = levels(factor(cell_df$group)),
        palette = palette,
        palcolor = palcolor
      )
      p <- ggplot2::ggplot(cell_df, ggplot2::aes(x, y, color = group)) +
        ggplot2::geom_point(size = 0.8, alpha = 0.75) +
        ggplot2::scale_color_manual(values = group_colors)
    } else {
      time_key <- if ("palantir_pseudotime" %in% colnames(srt@meta.data)) {
        "palantir_pseudotime"
      } else {
        "cellrank_pseudotime"
      }
      cell_df$value <- srt@meta.data[[time_key]]
      p <- ggplot2::ggplot(cell_df, ggplot2::aes(x, y, color = value)) +
        ggplot2::geom_point(size = 0.8, alpha = 0.75) +
        ggplot2::scale_color_gradientn(colors = continuous_colors, name = time_key)
    }
    return(
      p +
        ggplot2::geom_segment(
          data = arrow_df,
          ggplot2::aes(x, y, xend = xend, yend = yend),
          inherit.aes = FALSE,
          color = "black",
          linewidth = 0.18,
          alpha = 0.55,
          arrow = grid::arrow(length = grid::unit(1.3, "mm"), type = "closed")
        ) +
        ggplot2::coord_equal() +
        theme_layer +
        ggplot2::labs(x = colnames(coords)[1], y = colnames(coords)[2], color = group.by)
    )
  }
  if (plot_type %in% c("drivers", "trends", "clusters", "enrichment") && is.null(lineage)) {
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
        ggplot2::scale_color_gradientn(
          colors = continuous_colors,
          name = "Fate confidence",
          limits = c(0, 1),
          oob = scales::squish
        ) +
        theme_layer +
        ggplot2::theme(
          axis.title = ggplot2::element_blank(),
          axis.text = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank(),
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
    bar_color <- unname(palette_colors(n = 1L, palette = palette, palcolor = palcolor))
    return(ggplot2::ggplot(tab, ggplot2::aes(gene, correlation)) + ggplot2::geom_col(fill = bar_color) + ggplot2::coord_flip() + theme_layer + ggplot2::labs(x = NULL, y = "Correlation"))
  }

  if (plot_type == "enrichment") {
    enrichment <- srt@tools$CellRank$enrichment[[lineage]]$databases
    if (is.null(enrichment) || !length(enrichment)) {
      log_message("Run {.fn RunCellRankEnrichment} before plotting enrichment", message_type = "error")
    }
    if (is.null(database)) {
      nonempty <- names(enrichment)[vapply(enrichment, nrow, integer(1)) > 0L]
      database <- if (length(nonempty)) nonempty[[1L]] else NULL
    }
    if (is.null(database) || !database %in% names(enrichment) || !nrow(enrichment[[database]])) {
      log_message("Stored enrichment table is empty for the requested database", message_type = "error")
    }
    tab <- as.data.frame(enrichment[[database]], stringsAsFactors = FALSE)
    required <- c("Description", "p.adjust", "Count", "cluster")
    if (!all(required %in% names(tab))) {
      log_message("Stored enrichment table is missing plotting columns", message_type = "error")
    }
    tab <- do.call(rbind, lapply(split(tab, tab$cluster), function(x) {
      head(x[order(x$p.adjust, -x$Count), , drop = FALSE], as.integer(top_n))
    }))
    tab$cluster <- factor(tab$cluster, levels = unique(as.character(tab$cluster)))
    tab$Description <- factor(tab$Description, levels = rev(unique(as.character(tab$Description))))
    tab$minus_log10_q <- -log10(pmax(as.numeric(tab$p.adjust), .Machine$double.xmin))
    return(
      ggplot2::ggplot(
        tab,
        ggplot2::aes(cluster, Description, size = Count, color = minus_log10_q)
      ) +
        ggplot2::geom_point(alpha = 0.9) +
        ggplot2::scale_color_gradientn(colors = continuous_colors, name = expression(-log[10](adjusted~italic(P)))) +
        ggplot2::labs(x = "Trend module", y = NULL, title = database) +
        theme_layer
    )
  }

  if (plot_type %in% c("trends", "clusters")) {
    trend <- srt@tools$CellRank$trends[[lineage]]
    if (is.null(trend)) log_message("Run {.fn RunCellRankTrends} before plotting trends", message_type = "error")
    mat <- as.matrix(trend$trend_matrix)
    if (plot_type == "trends") {
      genes <- intersect(trend$heatmap_genes %||% rownames(mat), rownames(mat))
      long <- reshape2::melt(mat[genes, , drop = FALSE], varnames = c("gene", "time"), value.name = "value")
      return(ggplot2::ggplot(long, ggplot2::aes(time, gene, fill = value)) + ggplot2::geom_raster() + ggplot2::scale_fill_gradientn(colors = continuous_colors) + theme_layer + ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()) + ggplot2::labs(x = "Pseudotime", y = NULL))
    }
    clusters <- trend$cluster_table
    selected <- clusters$gene
    long <- reshape2::melt(mat[selected, , drop = FALSE], varnames = c("gene", "time"), value.name = "value")
    long$cluster <- clusters$cluster[match(long$gene, clusters$gene)]
    cluster_colors <- palette_colors(
      x = levels(factor(long$cluster)),
      palette = palette,
      palcolor = palcolor
    )
    return(ggplot2::ggplot(long, ggplot2::aes(time, value, group = gene, color = cluster)) + ggplot2::geom_line(alpha = 0.25) + ggplot2::stat_summary(ggplot2::aes(group = cluster), fun = mean, geom = "line", linewidth = 1.2) + ggplot2::scale_color_manual(values = cluster_colors) + ggplot2::facet_wrap(~cluster, scales = "free_y") + theme_layer + ggplot2::theme(legend.position = "none") + ggplot2::labs(x = "Pseudotime", y = "Normalized trend"))
  }

  transition <- srt@tools$CellRank$transition_matrix %||% srt@graphs[["cellrank_transition"]]
  if (is.null(transition)) log_message("Stored CellRank transition matrix is missing", message_type = "error")
  reduction <- reduction %||% DefaultReduction(srt)
  coords <- SeuratObject::Embeddings(srt[[reduction]])[, 1:2, drop = FALSE]
  cells <- rownames(coords)
  if (!is.null(rownames(transition)) && !is.null(colnames(transition))) {
    if (length(setdiff(cells, rownames(transition))) ||
        length(setdiff(cells, colnames(transition)))) {
      log_message("Stored CellRank transition matrix is not aligned to current cells", message_type = "error")
    }
    transition <- transition[cells, cells, drop = FALSE]
  } else if (!identical(dim(transition), c(length(cells), length(cells)))) {
    log_message("Stored CellRank transition matrix has incompatible dimensions", message_type = "error")
  }
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
  walk_palette <- unname(palette_colors(n = 2L, palette = palette, palcolor = palcolor))
  walk_color <- walk_palette[[min(2L, length(walk_palette))]]
  ggplot2::ggplot(path_df, ggplot2::aes(x, y, group = sim)) + ggplot2::geom_path(alpha = 0.35, color = walk_color, linewidth = 0.35) + ggplot2::coord_equal() + theme_layer + ggplot2::labs(x = colnames(coords)[1], y = colnames(coords)[2])
}
