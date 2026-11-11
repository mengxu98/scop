#' @title Plot CellRank results
#'
#' @param srt Deprecated alias for `object`; supply exactly one of the two. It
#' will be removed in scop 1.0.0.
#' @param object A Seurat object returned by [RunCellRank].
#' @param plot_type One of `"fate"`, `"states"`, `"circular"`, `"drivers"`,
#' `"projection"`, `"random_walks"`, or `"flow"`.
#' @param lineage Lineage used for driver plots.
#' @param palette,palcolor Discrete SCOP palette for states and cell
#' groups.
#' @param feature_palette,feature_palcolor Continuous SCOP palette for fate,
#' confidence, and pseudotime.
#' @param theme_use,theme_args SCOP plot theme and arguments passed to it.
#' @param reduction Reduction used for cell-space plots.
#' @param group.by Group column used for state/random-walk plots.
#' @param top_n Number of driver genes to display.
#' @param n_sims Number of random walks.
#' @param max_iter Maximum length of each random walk.
#' @param seed Random seed. Random-walk plotting restores the caller's RNG state,
#'   including when plotting fails.
#' @param start_cells Optional cell names used as random-walk starts. For a
#' moscot result, omitted starts default to cells in the earliest time point.
#' @param cluster Starting group for `plot_type = "flow"`.
#' @param clusters Optional groups to display in a real-time flow plot.
#' @param ... Arguments passed to other functions.
#' @return A ggplot or patchwork object.
#'
#' @seealso [DynamicHeatmap] and [DynamicPlot] for feature dynamics;
#' [EnrichmentPlot] for enrichment results.
#'
#' @export
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#' pancreas_sub <- RunCellRank(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   kernel_type = "pseudotime",
#'   time_key = "palantir_pseudotime",
#'   schur_method = "brandts",
#'   backend = "cpp"
#' )
#' CellRankPlot(
#'   pancreas_sub,
#'   plot_type = "fate",
#'   reduction = "UMAP",
#'   xlab = "UMAP_1",
#'   ylab = "UMAP_2"
#' )
#' CellRankPlot(
#'   pancreas_sub,
#'   plot_type = "circular"
#' )
CellRankPlot <- function(
  object,
  plot_type = c(
    "fate",
    "states",
    "circular",
    "drivers",
    "projection",
    "random_walks",
    "flow"
  ),
  lineage = NULL,
  reduction = NULL,
  group.by = NULL,
  top_n = 20L,
  n_sims = 100L,
  max_iter = 500L,
  seed = 0L,
  start_cells = NULL,
  cluster = NULL,
  clusters = NULL,
  palette = "Chinese",
  palcolor = NULL,
  feature_palette = "Spectral",
  feature_palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  srt = NULL,
  ...
) {
  srt <- resolve_deprecated_srt(object, srt, missing(object))
  plot_type <- match.arg(plot_type)
  plot_args <- list(...)
  draw_plot <- function(fun, args) {
    args <- utils::modifyList(args, plot_args, keep.null = TRUE)
    if (plot_type %in% c("projection", "random_walks")) args$combine <- FALSE
    do.call(fun, args)
  }
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a Seurat object", message_type = "error")
  }
  if (is.null(srt@tools$CellRank)) {
    log_message("CellRank results are missing", message_type = "error")
  }
  if (!is.null(reduction)) {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  fate <- tryCatch(cellrank_fate_matrix(srt), error = function(e) NULL)
  theme_layer <- apply_plot_theme(theme_use, theme_args)
  continuous_colors <- unname(palette_colors(
    n = 11L,
    palette = feature_palette,
    palcolor = feature_palcolor
  ))

  if (plot_type == "fate") {
    reduction <- reduction %||% DefaultReduction(srt)
    if (is.null(fate)) {
      log_message(
        "CellRank fate probabilities are missing",
        message_type = "error"
      )
    }
    cols <- paste0("cellrank_fate_", make.names(colnames(fate), unique = TRUE))
    cols <- intersect(cols, colnames(srt@meta.data))
    if (!length(cols)) {
      log_message(
        "No CellRank fate metadata columns are available",
        message_type = "error"
      )
    }
    lineage_labels <- colnames(fate)[match(cols, paste0("cellrank_fate_", make.names(colnames(fate), unique = TRUE)))]
    fate_plots <- draw_plot(FeatureDimPlot, list(
      object = srt,
      features = cols,
      reduction = reduction,
      palette = feature_palette,
      palcolor = feature_palcolor,
      theme_use = theme_use,
      theme_args = theme_args,
      show_stat = FALSE,
      legend.title = "Fate probability",
      lower_cutoff = 0,
      upper_cutoff = 1,
      combine = FALSE
    ))
    if (!is.list(fate_plots) || inherits(fate_plots, "ggplot")) {
      return(fate_plots)
    }
    for (index in seq_along(fate_plots)) {
      fate_plots[[index]] <- fate_plots[[index]] +
        ggplot2::labs(title = plot_args$title %||% lineage_labels[[index]]) +
        ggplot2::theme(strip.text.x = ggplot2::element_blank())
    }
    if (identical(plot_args$combine, FALSE)) {
      return(fate_plots)
    }
    return(patchwork::wrap_plots(
      fate_plots,
      ncol = plot_args$ncol %||% min(3L, length(cols)),
      nrow = plot_args$nrow, byrow = plot_args$byrow %||% TRUE
    ))
  }
  if (plot_type == "states") {
    reduction <- reduction %||% DefaultReduction(srt)
    return(draw_plot(CellDimPlot, list(
      object = srt,
      group.by = group.by %||% "cellrank_terminal_states",
      reduction = reduction,
      palette = palette,
      palcolor = palcolor,
      theme_use = theme_use,
      theme_args = theme_args,
      show_stat = FALSE
    )))
  }
  if (plot_type == "projection") {
    reduction <- reduction %||% DefaultReduction(srt)
    transition <- srt@tools$CellRank$transition_matrix %||%
      srt@graphs[["cellrank_transition"]]
    if (is.null(transition)) {
      log_message(
        "Stored CellRank transition matrix is missing",
        message_type = "error"
      )
    }
    cells <- colnames(srt)
    coords <- SeuratObject::Embeddings(srt[[reduction]])[cells, 1:2, drop = FALSE]
    if (!is.null(rownames(transition)) && !is.null(colnames(transition))) {
      missing_cells <- union(
        setdiff(cells, rownames(transition)),
        setdiff(cells, colnames(transition))
      )
      if (length(missing_cells)) {
        log_message(
          "Stored CellRank transition matrix is not aligned to current cells",
          message_type = "error"
        )
      }
      transition <- transition[cells, cells, drop = FALSE]
    } else if (!identical(dim(transition), c(length(cells), length(cells)))) {
      log_message(
        "Stored CellRank transition matrix has incompatible dimensions",
        message_type = "error"
      )
    }
    next_xy <- as.matrix(transition %*% coords)
    delta <- next_xy - coords
    magnitude <- sqrt(rowSums(delta^2))
    scale_ref <- as.numeric(stats::quantile(
      magnitude[magnitude > 0],
      0.95,
      na.rm = TRUE
    ))
    if (!is.finite(scale_ref) || scale_ref <= 0) {
      scale_ref <- 1
    }
    delta <- delta / scale_ref * 0.8
    n_arrows <- min(nrow(coords), 250L)
    arrow_idx <- unique(as.integer(round(seq(
      1,
      nrow(coords),
      length.out = n_arrows
    ))))
    arrow_df <- data.frame(
      x = coords[arrow_idx, 1],
      y = coords[arrow_idx, 2],
      xend = coords[arrow_idx, 1] + delta[arrow_idx, 1],
      yend = coords[arrow_idx, 2] + delta[arrow_idx, 2]
    )
    if (!is.null(group.by) && group.by %in% colnames(srt@meta.data)) {
      p <- draw_plot(CellDimPlot, list(
        object = srt, group.by = group.by, reduction = reduction,
        palette = palette, palcolor = palcolor,
        theme_use = theme_use, theme_args = theme_args,
        show_stat = FALSE, combine = FALSE
      ))[[1L]]
    } else {
      time_key <- if ("palantir_pseudotime" %in% colnames(srt@meta.data)) {
        "palantir_pseudotime"
      } else {
        "cellrank_pseudotime"
      }
      if (!time_key %in% colnames(srt@meta.data)) {
        time_key <- "cellrank_fate_confidence"
      }
      p <- draw_plot(FeatureDimPlot, list(
        object = srt, features = time_key, reduction = reduction,
        palette = feature_palette, palcolor = feature_palcolor,
        theme_use = theme_use, theme_args = theme_args,
        show_stat = FALSE, combine = FALSE
      ))[[1L]]
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
        )
    )
  }
  if (
    plot_type == "drivers" &&
      is.null(lineage)
  ) {
    log_message(
      "{.arg lineage} is required for this plot type",
      message_type = "error"
    )
  }

  if (plot_type == "circular") {
    if (is.null(fate) || ncol(fate) < 2L) {
      log_message(
        "At least two fate lineages are required",
        message_type = "error"
      )
    }
    fate <- as.matrix(fate)
    fate <- fate / pmax(rowSums(fate), .Machine$double.eps)
    angles <- seq(0, 2 * pi, length.out = ncol(fate) + 1L)[-(ncol(fate) + 1L)]
    xy <- cbind(
      x = as.numeric(fate %*% cos(angles)),
      y = as.numeric(fate %*% sin(angles))
    )
    xy <- data.frame(
      xy,
      cell = rownames(fate),
      confidence = apply(fate, 1L, max)
    )
    vertices <- data.frame(
      x = cos(angles),
      y = sin(angles),
      label = colnames(fate)
    )
    circular_plot <- draw_plot(thisplot::DimDataPlot, list(
      data = xy, x = "x", y = "y",
      palette = feature_palette, palcolor = feature_palcolor,
      pt.size = 0.8, pt.alpha = 0, legend.title = "Fate confidence",
      theme_use = theme_use, theme_args = theme_args
    ))
    circular_plot$layers <- list()
    return(
      circular_plot +
        ggplot2::geom_point(
          data = xy, ggplot2::aes(x, y, color = confidence),
          inherit.aes = FALSE, size = plot_args$pt.size %||% 0.8,
          alpha = plot_args$pt.alpha %||% 0.7
        ) +
        ggplot2::geom_point(
          data = vertices,
          ggplot2::aes(x, y),
          inherit.aes = FALSE,
          size = 3
        ) +
        ggplot2::geom_text(
          data = vertices,
          ggplot2::aes(x, y, label = label),
          inherit.aes = FALSE,
          vjust = -1
        ) +
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
          legend.background = ggplot2::element_rect(
            fill = "white",
            colour = NA
          ),
          plot.margin = ggplot2::margin(12, 36, 12, 36)
        )
    )
  }

  if (plot_type == "drivers") {
    drivers <- srt@tools$CellRank$lineage_drivers
    corr_col <- paste0(lineage, "_corr")
    if (is.null(drivers) || !corr_col %in% colnames(drivers)) {
      log_message("Stored lineage drivers are missing", message_type = "error")
    }
    corr_values <- as.numeric(drivers[, corr_col])
    keep <- is.finite(corr_values)
    tab <- data.frame(
      gene = rownames(drivers)[keep],
      correlation = corr_values[keep],
      stringsAsFactors = FALSE
    )
    tab <- head(
      tab[order(tab$correlation, decreasing = TRUE), , drop = FALSE],
      as.integer(top_n)
    )
    return(draw_plot(thisplot::StatPlot, list(
      meta.data = tab, stat.by = "gene", value.by = "correlation",
      stat_type = "value", plot_type = "bar", flip = TRUE,
      rank.by = "value", palette = palette, palcolor = palcolor,
      theme_use = theme_use, theme_args = theme_args,
      xlab = NULL, ylab = "Correlation", legend.position = "none"
    )))
  }

  if (plot_type == "flow") {
    if (is.null(group.by) || !group.by %in% colnames(srt@meta.data)) {
      log_message(
        "{.arg group.by} is required for a real-time flow plot",
        message_type = "error"
      )
    }
    temporal <- srt@tools$CellRank$temporal
    if (is.null(temporal) || is.null(temporal$time_key)) {
      log_message(
        "Stored real-time metadata are missing; run RunCellRank(kernel_type = 'moscot') first",
        message_type = "error"
      )
    }
    time_key <- as.character(temporal$time_key)
    if (!time_key %in% colnames(srt@meta.data)) {
      log_message(
        "Stored real-time column {.val {time_key}} is missing from Seurat metadata",
        message_type = "error"
      )
    }
    transition <- srt@tools$CellRank$transition_matrix %||%
      srt@graphs[["cellrank_transition"]]
    if (is.null(transition)) {
      log_message(
        "Stored CellRank transition matrix is missing",
        message_type = "error"
      )
    }
    cells <- colnames(srt)
    if (!is.null(rownames(transition)) && !is.null(colnames(transition))) {
      transition <- transition[cells, cells, drop = FALSE]
    }
    time_numeric <- suppressWarnings(as.numeric(as.character(srt@meta.data[[
      time_key
    ]])))
    if (any(!is.finite(time_numeric))) {
      log_message(
        "Stored real-time values are not numeric",
        message_type = "error"
      )
    }
    time_levels <- as.numeric(
      temporal$time_values %||% sort(unique(time_numeric))
    )
    time_levels <- sort(unique(time_levels))
    groups <- as.character(srt@meta.data[[group.by]])
    earliest <- time_levels[[1L]]
    if (is.null(cluster)) {
      cluster <- unique(groups[time_numeric == earliest])[1L]
    }
    starts <- which(time_numeric == earliest & groups == as.character(cluster))
    if (!length(starts)) {
      log_message(
        "No cells found for the requested initial cluster and earliest time",
        message_type = "error"
      )
    }
    keep_groups <- as.character(clusters %||% sort(unique(groups)))
    mass <- numeric(length(cells))
    mass[starts] <- 1 / length(starts)
    rows <- vector("list", length(time_levels))
    for (i in seq_along(time_levels)) {
      if (i > 1L) {
        mass <- as.numeric(mass %*% transition)
        mass[time_numeric != time_levels[[i]]] <- 0
        total <- sum(mass)
        if (is.finite(total) && total > 0) mass <- mass / total
      }
      tab <- stats::aggregate(
        mass,
        by = list(cluster = groups, time = time_numeric),
        FUN = sum
      )
      tab <- tab[
        tab$time == time_levels[[i]] & tab$cluster %in% keep_groups, ,
        drop = FALSE
      ]
      names(tab)[3L] <- "probability"
      rows[[i]] <- tab
    }
    flow <- do.call(rbind, rows)
    if (is.null(flow) || !nrow(flow)) {
      log_message("No real-time flow could be computed", message_type = "error")
    }
    complete_flow <- expand.grid(
      cluster = keep_groups,
      time = time_levels,
      stringsAsFactors = FALSE
    )
    flow <- merge(
      complete_flow,
      flow,
      by = c("cluster", "time"),
      all.x = TRUE,
      sort = FALSE
    )
    flow$probability[is.na(flow$probability)] <- 0
    flow$cluster <- factor(flow$cluster, levels = keep_groups)
    flow$time <- as.numeric(flow$time)
    return(
      ggplot2::ggplot(flow, ggplot2::aes(time, .data$probability, fill = cluster)) +
        ggplot2::geom_area(alpha = 0.85, colour = "white", linewidth = 0.2) +
        ggplot2::scale_fill_manual(values = palette_colors(
          x = keep_groups, palette = palette, palcolor = palcolor
        ), drop = FALSE) +
        ggplot2::scale_x_continuous(breaks = time_levels) +
        theme_layer +
        ggplot2::theme(aspect.ratio = plot_args$aspect.ratio %||% 0.7) +
        ggplot2::labs(
          x = plot_args$xlab %||% "Experimental time",
          y = plot_args$ylab %||% "Propagated probability mass", fill = group.by,
          title = plot_args$title %||% paste0("Flow from ", cluster, " at time ", earliest)
        )
    )
  }

  transition <- srt@tools$CellRank$transition_matrix %||%
    srt@graphs[["cellrank_transition"]]
  if (is.null(transition)) {
    log_message(
      "Stored CellRank transition matrix is missing",
      message_type = "error"
    )
  }
  reduction <- reduction %||% DefaultReduction(srt)
  coords <- SeuratObject::Embeddings(srt[[reduction]])[, 1:2, drop = FALSE]
  cells <- rownames(coords)
  if (!is.null(rownames(transition)) && !is.null(colnames(transition))) {
    if (
      length(setdiff(cells, rownames(transition))) ||
        length(setdiff(cells, colnames(transition)))
    ) {
      log_message(
        "Stored CellRank transition matrix is not aligned to current cells",
        message_type = "error"
      )
    }
    transition <- transition[cells, cells, drop = FALSE]
  } else if (!identical(dim(transition), c(length(cells), length(cells)))) {
    log_message(
      "Stored CellRank transition matrix has incompatible dimensions",
      message_type = "error"
    )
  }
  check_r("withr", verbose = FALSE)
  withr::local_preserve_seed()
  set.seed(seed)
  starts <- seq_len(nrow(coords))
  if (!is.null(start_cells)) {
    start_cells <- intersect(as.character(start_cells), cells)
    if (!length(start_cells)) {
      log_message(
        "{.arg start_cells} does not match stored CellRank cells",
        message_type = "error"
      )
    }
    starts <- match(start_cells, cells)
  } else if (
    !is.null(srt@tools$CellRank$temporal$time_key) &&
      srt@tools$CellRank$temporal$time_key %in% colnames(srt@meta.data)
  ) {
    time_values <- suppressWarnings(as.numeric(as.character(srt@meta.data[[
      srt@tools$CellRank$temporal$time_key
    ]])))
    if (any(is.finite(time_values))) {
      starts <- which(time_values == min(time_values, na.rm = TRUE))
    }
  } else if (!is.null(group.by) && group.by %in% colnames(srt@meta.data)) {
    starts <- which(
      as.character(srt@meta.data[[group.by]]) ==
        unique(as.character(srt@meta.data[[group.by]]))[1L]
    )
  }
  paths <- vector("list", as.integer(n_sims))
  for (sim in seq_len(as.integer(n_sims))) {
    current <- sample(starts, 1L)
    visited <- integer(as.integer(max_iter))
    visited[1L] <- current
    for (step in 2:as.integer(max_iter)) {
      probs <- as.numeric(transition[current, ])
      if (!any(probs > 0)) {
        visited <- visited[seq_len(step - 1L)]
        break
      }
      current <- sample.int(nrow(coords), 1L, prob = probs / sum(probs))
      visited[step] <- current
    }
    paths[[sim]] <- data.frame(
      sim = sim,
      step = seq_along(visited),
      coords[visited, , drop = FALSE]
    )
  }
  path_df <- do.call(rbind, paths)
  names(path_df)[3:4] <- c("x", "y")
  walk_palette <- unname(palette_colors(
    n = 2L,
    palette = palette,
    palcolor = palcolor
  ))
  walk_color <- walk_palette[[min(2L, length(walk_palette))]]
  if (is.null(group.by)) {
    group.by <- "cellrank_plot_cells"
    srt[[group.by]] <- factor(rep("Cells", ncol(srt)))
  }
  background <- draw_plot(CellDimPlot, list(
    object = srt, group.by = group.by,
    reduction = reduction, palette = palette, palcolor = palcolor,
    theme_use = theme_use, theme_args = theme_args,
    show_stat = FALSE, pt.alpha = 0.35, combine = FALSE
  ))[[1L]]
  background + ggplot2::geom_path(
    data = path_df, ggplot2::aes(x, y, group = sim), inherit.aes = FALSE,
    alpha = 0.45, color = walk_color, linewidth = 0.35
  )
}
