#' @title Plot Giotto backend results
#'
#' @description
#' Plot standalone Giotto backend results with scop plotting conventions. The
#' input `Seurat` object, when supplied, is copied internally for plotting and
#' is not modified.
#'
#' @md
#' @param x A result returned by one of the `RunGiotto*()` functions.
#' @param ... Arguments passed to S3 methods.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' cluster_result <- list(
#'   clusters = data.frame(
#'     cluster = paste0("cluster_", (seq_len(ncol(spatial)) - 1) %% 3 + 1),
#'     row.names = colnames(spatial)
#'   ),
#'   parameters = list(
#'     cluster_colname = "Giotto_cluster",
#'     coord.cols = c("x", "y"),
#'     k = 8,
#'     resolution = 0.4
#'   )
#' )
#' class(cluster_result) <- c("giotto2_cluster", "giotto2_result", "list")
#' GiottoPlot(
#'   cluster_result,
#'   srt = spatial,
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#'
#' proximity <- list(
#'   enrichment = data.frame(
#'     group_1 = c("Tumor", "Tumor", "Stroma", "Immune"),
#'     group_2 = c("Stroma", "Immune", "Immune", "Tumor"),
#'     enrichment = c(1.6, -0.7, 0.9, -1.2),
#'     type_int = c("enriched", "depleted", "enriched", "depleted")
#'   ),
#'   parameters = list(network_method = "Delaunay", number_of_simulations = 100)
#' )
#' class(proximity) <- c("giotto2_cell_proximity", "giotto2_result", "list")
#' GiottoPlot(proximity)
#'
#' spatial_genes <- list(
#'   results = data.frame(
#'     feat_ID = c("COL1A1", "KRT19", "MS4A1", "PECAM1"),
#'     spatGeneRank = c(41.2, 32.8, 18.4, 11.9)
#'   ),
#'   top_features = c("COL1A1", "KRT19", "MS4A1"),
#'   parameters = list(assay = "Spatial", layer = "data")
#' )
#' class(spatial_genes) <- c("giotto2_spatial_genes", "giotto2_result", "list")
#' GiottoPlot(spatial_genes, plot_type = "ranking", top_n = 4)
#'
#' modules <- list(
#'   module_tables = list(
#'     result.cor_DT = expand.grid(
#'       feat_ID = c("COL1A1", "KRT19", "MS4A1"),
#'       variable = c("COL1A1", "KRT19", "MS4A1")
#'     )
#'   ),
#'   features = c("COL1A1", "KRT19", "MS4A1")
#' )
#' modules$module_tables$result.cor_DT$spat_cor <- c(
#'   1, 0.35, -0.20,
#'   0.35, 1, 0.15,
#'   -0.20, 0.15, 1
#' )
#' class(modules) <- c("giotto2_spatial_modules", "giotto2_result", "list")
#' GiottoPlot(modules, top_n = 3)
#'
#' @export
GiottoPlot <- function(x, ...) {
  UseMethod("GiottoPlot")
}

#' @export
GiottoPlot.default <- function(x, ...) {
  log_message(
    "{.arg x} must be a standalone result returned by a {.fn RunGiotto*} function",
    message_type = "error"
  )
}

#' @rdname GiottoPlot
#' @param srt Original `Seurat` object used to create the Giotto result. Required
#' for spatial spot plots.
#' @param image Name of the Seurat spatial image. Required when multiple images
#' are present; a single image is selected automatically when `NULL`.
#' @param coord.cols Metadata coordinate columns used when no image is available.
#' @param overlay_image Whether to draw the spatial image beneath spots.
#' @param crop Whether to crop spatial panels to plotted spots.
#' @param pt.size Point size for spatial plots.
#' @param pt.alpha Point alpha for spatial plots.
#' @param stroke Point border width for discrete spatial plots.
#' @param palette Discrete palette used for groups.
#' @param feature_palette Continuous palette used for spatial expression plots.
#' @param bg_color Point border color for discrete spatial plots.
#' @param legend.position Legend position.
#' @param theme_use Theme function name used by scop plots.
#' @param theme_args Additional arguments passed to `theme_use`.
#' @param title,subtitle Plot title and subtitle. If `NULL`, sensible defaults
#' are used.
#' @export
GiottoPlot.giotto2_cluster <- function(
  x,
  srt,
  image = x$parameters$image %||% NULL,
  coord.cols = x$parameters$coord.cols %||% c("x", "y"),
  overlay_image = TRUE,
  crop = TRUE,
  pt.size = NULL,
  pt.alpha = 0.95,
  stroke = 0.08,
  palette = "Chinese",
  feature_palette = "Spectral",
  bg_color = "grey25",
  legend.position = "right",
  theme_use = "theme_blank",
  theme_args = list(),
  title = "Giotto Leiden clusters",
  subtitle = NULL,
  ...
) {
  giotto_plot_require_srt(srt)
  giotto_plot_require_table(x$clusters, "clusters")
  if (!"cluster" %in% colnames(x$clusters)) {
    log_message("Giotto cluster result must contain a {.field cluster} column", message_type = "error")
  }
  cells <- intersect(rownames(x$clusters), colnames(srt))
  if (length(cells) == 0L) {
    log_message("No Giotto cluster cells match {.arg srt}", message_type = "error")
  }

  cluster_colname <- x$parameters$cluster_colname %||% "Giotto_cluster"
  plot_srt <- srt
  plot_srt@meta.data[[cluster_colname]] <- NA_character_
  plot_srt@meta.data[cells, cluster_colname] <- as.character(x$clusters[cells, "cluster"])
  plot_srt@meta.data[[cluster_colname]] <- factor(
    plot_srt@meta.data[[cluster_colname]],
    levels = sort(unique(stats::na.omit(plot_srt@meta.data[[cluster_colname]])))
  )

  subtitle <- subtitle %||% giotto_cluster_subtitle(x, length(cells))
  p <- SpatialSpotPlot(
    plot_srt,
    group.by = cluster_colname,
    image = image,
    overlay_image = overlay_image,
    crop = crop,
    coord.cols = coord.cols,
    flip.y = TRUE,
    cells = cells,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    stroke = stroke,
    palette = palette,
    bg_color = bg_color,
    legend.position = legend.position,
    legend.title = cluster_colname,
    theme_use = theme_use,
    theme_args = theme_args,
    ...
  )
  giotto_plot_titles(p, title = title, subtitle = subtitle)
}

#' @export
plot.giotto2_cluster <- function(x, y = NULL, ...) {
  GiottoPlot(x, ...)
}

#' @rdname GiottoPlot
#' @param heatmap_palette Continuous palette used for heatmaps.
#' @param heatmap_palcolor Optional custom colors used to create
#' `heatmap_palette`.
#' @export
GiottoPlot.giotto2_cell_proximity <- function(
  x,
  heatmap_palette = "RdBu",
  heatmap_palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  title = "Giotto cell proximity enrichment",
  subtitle = NULL,
  ...
) {
  giotto_plot_require_table(x$enrichment, "enrichment")
  dat <- x$enrichment
  score_col <- giotto_plot_pick_numeric_col(
    dat,
    c("enrichm", "enrichment", "PI_value", "observed_over_expected", "zscore", "score")
  )
  if (!all(c("group_1", "group_2") %in% colnames(dat))) {
    dat <- giotto_add_interaction_columns(dat)
  }
  if (!all(c("group_1", "group_2") %in% colnames(dat))) {
    log_message("Giotto proximity result must contain pairwise group columns", message_type = "error")
  }
  dat[[".pair"]] <- paste(dat[["group_1"]], dat[["group_2"]], sep = "--")
  dat[[".type"]] <- if ("type_int" %in% colnames(dat)) as.character(dat[["type_int"]]) else "pair"
  dat <- dat[order(dat[[".type"]], dat[["group_1"]], dat[["group_2"]]), , drop = FALSE]
  dat[[".type"]] <- factor(dat[[".type"]], levels = unique(dat[[".type"]]))
  dat[[".pair"]] <- factor(dat[[".pair"]], levels = unique(dat[[".pair"]]))

  subtitle <- subtitle %||% giotto_proximity_subtitle(x)
  cols <- giotto_plot_diverging_colors(heatmap_palette, heatmap_palcolor)
  ggplot2::ggplot(dat, ggplot2::aes(x = .data$.pair, y = .data$.type, fill = .data[[score_col]])) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::scale_fill_gradient2(
      low = cols[[1L]],
      mid = cols[[2L]],
      high = cols[[3L]],
      midpoint = 0,
      na.value = "grey90"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = NULL,
      fill = giotto_plot_pretty_label(score_col)
    ) +
    giotto_plot_theme(theme_use, theme_args) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.text.y = ggplot2::element_text(size = 9),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(color = "grey35", size = 10),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9),
      legend.position = "right"
    )
}

#' @export
plot.giotto2_cell_proximity <- function(x, y = NULL, ...) {
  GiottoPlot(x, ...)
}

#' @rdname GiottoPlot
#' @param plot_type Plot type for spatial gene results. `"ranking"` plots the
#' feature-level table; `"feature"` plots expression of one feature on spatial
#' coordinates.
#' @param feature Feature to draw for `plot_type = "feature"`. If `NULL`, the
#' top Giotto feature is used.
#' @param top_n Number of rows shown in ranking plots.
#' @param assay Assay used for spatial feature expression plots.
#' @param layer Assay layer used for spatial feature expression plots.
#' @export
GiottoPlot.giotto2_spatial_genes <- function(
  x,
  srt = NULL,
  plot_type = c("ranking", "feature"),
  feature = NULL,
  top_n = 20,
  assay = x$parameters$assay %||% NULL,
  layer = x$parameters$layer %||% "data",
  image = x$parameters$image %||% NULL,
  coord.cols = x$parameters$coord.cols %||% c("x", "y"),
  overlay_image = TRUE,
  crop = TRUE,
  pt.size = NULL,
  pt.alpha = 0.95,
  palette = "Chinese",
  feature_palette = "Spectral",
  heatmap_palette = "RdBu",
  heatmap_palcolor = NULL,
  legend.position = "right",
  theme_use = "theme_scop",
  theme_args = list(),
  title = NULL,
  subtitle = NULL,
  ...
) {
  plot_type <- match.arg(plot_type)
  if (identical(plot_type, "feature")) {
    giotto_plot_require_srt(srt)
    feature <- feature %||% x$top_features[[1L]]
    if (is.null(feature) || is.na(feature) || !feature %in% rownames(srt)) {
      log_message("Requested Giotto feature is not present in {.arg srt}", message_type = "error")
    }
    p <- SpatialSpotPlot(
      srt,
      features = feature,
      assay = assay,
      layer = layer,
      image = image,
      overlay_image = overlay_image,
      crop = crop,
      coord.cols = coord.cols,
      flip.y = TRUE,
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      palette = feature_palette,
      legend.position = legend.position,
      legend.title = feature,
      theme_use = "theme_blank",
      theme_args = list(),
      ...
    )
    return(giotto_plot_titles(
      p,
      title = title %||% paste0("Top Giotto spatial gene: ", feature),
      subtitle = subtitle %||% "Normalized expression on tissue coordinates"
    ))
  }

  giotto_plot_require_table(x$results, "results")
  dat <- utils::head(x$results, as.integer(top_n))
  feature_col <- giotto_plot_pick_col(dat, c("feats", "feat_ID", "gene", "feature", "features"))
  score_col <- giotto_plot_pick_numeric_col(dat, c("score", "spatGeneRank", "estimate", "statistic", "p.value", "adj.p.value"))
  dat[[".feature"]] <- factor(as.character(dat[[feature_col]]), levels = rev(unique(as.character(dat[[feature_col]]))))
  dat[[".score"]] <- dat[[score_col]]
  if (grepl("p[._]?value|p\\.adj|adj", score_col, ignore.case = TRUE)) {
    dat[[".score"]] <- -log10(pmax(dat[[".score"]], .Machine$double.xmin))
  }
  fill_col <- palette_colors(palette = palette, n = 1L)[[1L]]

  ggplot2::ggplot(dat, ggplot2::aes(x = .data$.score, y = .data$.feature)) +
    ggplot2::geom_col(fill = fill_col, width = 0.72) +
    ggplot2::geom_point(color = "white", size = 1.6) +
    ggplot2::labs(
      title = title %||% "Giotto binSpect spatial genes",
      subtitle = subtitle %||% paste0("Top ", nrow(dat), " genes ranked by Giotto ", giotto_plot_pretty_label(score_col)),
      x = giotto_plot_pretty_label(score_col),
      y = NULL
    ) +
    giotto_plot_theme(theme_use, theme_args) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(color = "grey35", size = 10),
      axis.text.y = ggplot2::element_text(size = 9)
    )
}

#' @export
plot.giotto2_spatial_genes <- function(x, y = NULL, ...) {
  GiottoPlot(x, ...)
}

#' @rdname GiottoPlot
#' @param features Features used for spatial co-expression heatmaps. If `NULL`,
#' top features from the Giotto result are used.
#' @export
GiottoPlot.giotto2_spatial_modules <- function(
  x,
  features = NULL,
  top_n = 20,
  heatmap_palette = "RdBu",
  heatmap_palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  title = "Giotto spatial co-expression",
  subtitle = "Spatial correlation among top features",
  ...
) {
  cor_table <- giotto_plot_module_cor_table(x)
  required <- c("feat_ID", "variable", "spat_cor")
  if (!all(required %in% colnames(cor_table))) {
    log_message(
      "Giotto spatial module plot requires columns {.val {required}} in the extracted correlation table",
      message_type = "error"
    )
  }
  features <- features %||% giotto_plot_module_features(cor_table, x, top_n)
  features <- intersect(unique(features), unique(c(as.character(cor_table$feat_ID), as.character(cor_table$variable))))
  if (length(features) == 0L) {
    log_message("No requested module features are available for plotting", message_type = "error")
  }
  dat <- cor_table[
    as.character(cor_table$feat_ID) %in% features &
      as.character(cor_table$variable) %in% features,
    ,
    drop = FALSE
  ]
  dat[["feat_ID"]] <- factor(as.character(dat[["feat_ID"]]), levels = features)
  dat[["variable"]] <- factor(as.character(dat[["variable"]]), levels = rev(features))
  cols <- giotto_plot_diverging_colors(heatmap_palette, heatmap_palcolor)

  ggplot2::ggplot(dat, ggplot2::aes(x = .data$feat_ID, y = .data$variable, fill = .data$spat_cor)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.12) +
    ggplot2::scale_fill_gradient2(
      low = cols[[1L]],
      mid = cols[[2L]],
      high = cols[[3L]],
      midpoint = 0,
      limits = c(-1, 1),
      na.value = "grey90"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = NULL,
      fill = "Spatial\ncorrelation"
    ) +
    giotto_plot_theme(theme_use, theme_args) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.text.y = ggplot2::element_text(size = 8),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(color = "grey35", size = 10),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9),
      legend.position = "right"
    )
}

#' @export
plot.giotto2_spatial_modules <- function(x, y = NULL, ...) {
  GiottoPlot(x, ...)
}

giotto_plot_require_srt <- function(srt) {
  if (missing(srt) || is.null(srt) || !inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be the original {.cls Seurat} object used for the Giotto result", message_type = "error")
  }
  invisible(TRUE)
}

giotto_plot_require_table <- function(x, name) {
  if (is.null(x) || !inherits(x, "data.frame") || nrow(x) == 0L) {
    log_message("Giotto result does not contain a non-empty {.field {name}} table", message_type = "error")
  }
  invisible(TRUE)
}

giotto_plot_pick_col <- function(dat, candidates) {
  hit <- candidates[candidates %in% colnames(dat)][1L]
  if (is.na(hit)) {
    log_message("None of the expected columns were found: {.val {candidates}}", message_type = "error")
  }
  hit
}

giotto_plot_pick_numeric_col <- function(dat, candidates) {
  numeric_cols <- colnames(dat)[vapply(dat, is.numeric, logical(1))]
  hit <- intersect(candidates, numeric_cols)[1L]
  if (is.na(hit)) {
    hit <- numeric_cols[1L]
  }
  if (is.na(hit)) {
    log_message("Giotto result does not contain a numeric column for plotting", message_type = "error")
  }
  hit
}

giotto_plot_theme <- function(theme_use = "theme_scop", theme_args = list()) {
  if (is.null(theme_use)) {
    return(ggplot2::theme_minimal())
  }
  if (inherits(theme_use, "theme")) {
    return(theme_use)
  }
  if (identical(theme_use, "theme_scop")) {
    theme_use <- "theme_this"
  }
  theme_fun <- if (is.character(theme_use)) {
    get(theme_use, mode = "function", inherits = TRUE)
  } else {
    theme_use
  }
  do.call(theme_fun, theme_args)
}

giotto_plot_titles <- function(p, title = NULL, subtitle = NULL) {
  p +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey35")
    )
}

giotto_plot_diverging_colors <- function(palette = "RdBu", palcolor = NULL) {
  cols <- palette_colors(palette = palette, palcolor = palcolor, n = 9L)
  if (length(cols) >= 9L) {
    cols <- cols[c(1L, 6L, 9L)]
  } else if (length(cols) < 3L) {
    cols <- c("#2166ac", "white", "#b2182b")
  } else {
    cols <- cols[c(1L, ceiling(length(cols) / 2), length(cols))]
  }
  cols
}

giotto_plot_pretty_label <- function(x) {
  out <- gsub("[._]+", " ", x)
  substr(out, 1L, 1L) <- toupper(substr(out, 1L, 1L))
  out
}

giotto_cluster_subtitle <- function(x, n_cells) {
  parts <- c(
    paste0(n_cells, " spots"),
    if (!is.null(x$parameters$k)) paste0("k=", x$parameters$k) else NULL,
    if (!is.null(x$parameters$resolution)) paste0("resolution=", x$parameters$resolution) else NULL
  )
  paste(parts, collapse = " | ")
}

giotto_proximity_subtitle <- function(x) {
  parts <- c(
    x$parameters$network_method %||% NULL,
    if (!is.null(x$parameters$number_of_simulations)) {
      paste0(x$parameters$number_of_simulations, " simulations")
    } else {
      NULL
    }
  )
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0L) {
    return(NULL)
  }
  paste(parts, collapse = ", ")
}

giotto_plot_module_cor_table <- function(x) {
  if (!is.null(x$module_tables) && "result.cor_DT" %in% names(x$module_tables)) {
    return(x$module_tables[["result.cor_DT"]])
  }
  if (!is.null(x$module_tables)) {
    for (tab in x$module_tables) {
      if (inherits(tab, "data.frame") && all(c("feat_ID", "variable", "spat_cor") %in% colnames(tab))) {
        return(tab)
      }
    }
  }
  log_message("Giotto spatial module result does not contain a spatial correlation table", message_type = "error")
}

giotto_plot_module_features <- function(cor_table, x, top_n) {
  if (!is.null(x$features)) {
    return(utils::head(x$features, as.integer(top_n)))
  }
  counts <- sort(table(as.character(cor_table$feat_ID)), decreasing = TRUE)
  utils::head(names(counts), as.integer(top_n))
}
