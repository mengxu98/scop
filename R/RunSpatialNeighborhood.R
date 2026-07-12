#' @title Run spatial neighborhood statistics
#'
#' @description
#' Build a scop-style spatial neighborhood result bundle and optionally dispatch
#' to a supported backend for colocalization or local-effect statistics.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt A `Seurat` object.
#' @param group.by Metadata column containing spatial cell or spot labels.
#' @param method Neighborhood calculation. `"observed"` returns native KNN or
#' radius summaries. `"spicyR"` runs differential neighborhood statistics and
#' requires `split.by`.
#' @param assay Assay used when `features` are requested.
#' @param layer Assay layer used when `features` are requested.
#' @param coord.cols Metadata coordinate columns used when no Seurat image
#' coordinates are available.
#' @param image Name of the Seurat spatial image. Required when multiple images
#' are present; a single image is selected automatically when `NULL`.
#' @param coordinate_space Coordinate system used to build neighbor distances.
#' @param sample.by Metadata column identifying images or samples. If `NULL`,
#' all spots are treated as one sample.
#' @param split.by Optional metadata column identifying conditions for
#' differential neighborhood statistics.
#' @param subject.by Optional metadata column identifying subjects for backends
#' that support paired or repeated designs.
#' @param radius Optional spatial radius used for scop-native neighborhood
#' summaries.
#' @param k Optional number of nearest neighbors used for scop-native
#' neighborhood summaries. When both `radius` and `k` are `NULL`, `k = 6` is
#' used.
#' @param features Optional features to extract into the backend input table.
#' @param from,to Optional cell or spot label filters.
#' @param tool_name Name used to store results in `srt@tools`.
#' @param store_results Whether to store results in `srt@tools`.
#' @param ... Additional arguments passed to the selected backend.
#'
#' @return A `Seurat` object with results stored in `srt@tools[[tool_name]]`.
#' @concept spatial-producer
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' spatial <- RunSpatialNeighborhood(
#'   spatial,
#'   group.by = "coda_label",
#'   coord.cols = c("x", "y"),
#'   k = 4,
#'   verbose = FALSE
#' )
#'
#' SpatialNeighborhoodPlot(spatial, plot_type = "heatmap")
#' SpatialNeighborhoodPlot(spatial, plot_type = "network", top_n = 12)
#' SpatialNeighborhoodPlot(spatial, plot_type = "stat", top_n = 12)
#' SpatialNeighborhoodPlot(
#'   spatial,
#'   plot_type = "spatial",
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
RunSpatialNeighborhood <- function(
  srt,
  group.by,
  method = c("observed", "spicyR"),
  assay = NULL,
  layer = "data",
  coord.cols = c("col", "row"),
  image = NULL,
  sample.by = NULL,
  split.by = NULL,
  subject.by = NULL,
  radius = NULL,
  k = NULL,
  features = NULL,
  from = NULL,
  to = NULL,
  tool_name = "SpatialNeighborhood",
  store_results = TRUE,
  verbose = TRUE,
  coordinate_space = c("legacy_display", "raw"),
  ...
) {
  method <- match.arg(method)
  coordinate_space <- match.arg(coordinate_space)
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (identical(method, "spicyR") && is.null(split.by)) {
    log_message(
      "{.arg split.by} is required when {.arg method = 'spicyR'}; the backend was not run",
      message_type = "error"
    )
  }
  spatial_neighborhood_assert_string(group.by, "group.by")
  spatial_neighborhood_assert_string(tool_name, "tool_name")

  input <- spatial_neighborhood_input(
    srt = srt,
    group.by = group.by,
    assay = assay,
    layer = layer,
    coord.cols = coord.cols,
    image = image,
    coordinate_space = coordinate_space,
    sample.by = sample.by,
    split.by = split.by,
    subject.by = subject.by,
    features = features
  )
  observed <- spatial_neighborhood_observed_pairs(
    cells = input$cells,
    radius = radius,
    k = k
  )

  backend <- switch(method,
    observed = list(raw = NULL, table = NULL),
    spicyR = spatial_neighborhood_run_spicyr(
      cells = input$cells,
      group.by = group.by,
      sample.by = sample.by,
      split.by = split.by,
      subject.by = subject.by,
      verbose = verbose,
      ...
    )
  )

  pair_table <- spatial_neighborhood_standardize_pair_table(
    backend = backend,
    observed = observed$pair_table,
    method = method
  )
  pair_table <- spatial_neighborhood_filter_pairs(
    pair_table,
    from = from,
    to = to
  )
  observed$pair_table <- spatial_neighborhood_filter_pairs(
    observed$pair_table,
    from = from,
    to = to
  )
  observed$edge_table <- spatial_neighborhood_filter_edges(
    observed$edge_table,
    from = from,
    to = to
  )
  long_table <- spatial_neighborhood_long_table(
    pair_table = pair_table,
    observed = observed$edge_table,
    method = method
  )

  method_bundle <- list(
    method = method,
    pair_table = pair_table,
    long_table = long_table,
    observed_table = observed$pair_table,
    edge_table = observed$edge_table,
    raw = backend$raw,
    input = input$cells,
    summary = scop_spatial_neighborhood_summary(pair_table, observed$edge_table),
    parameters = list(
      method = method,
      group.by = group.by,
      assay = input$assay,
      layer = layer,
      coord.cols = coord.cols,
      image = image,
      coordinate_space = coordinate_space,
      sample.by = sample.by,
      split.by = split.by,
      subject.by = subject.by,
      radius = radius,
      k = k %||% if (is.null(radius)) 6L else NULL,
      features = features,
      from = from,
      to = to
    )
  )

  if (isTRUE(store_results)) {
    old <- srt@tools[[tool_name]] %||% list()
    methods_store <- old$methods %||% list()
    methods_store[[method]] <- method_bundle
    srt@tools[[tool_name]] <- list(
      method = "SpatialNeighborhood",
      active_method = method,
      methods = methods_store,
      pair_table = pair_table,
      long_table = long_table,
      summary = method_bundle$summary,
      parameters = method_bundle$parameters
    )
    srt@tools[[tool_name]] <- spatial_result_build(
      bundle = srt@tools[[tool_name]],
      method = "SpatialNeighborhood",
      result_type = "neighborhood",
      provenance = list(
        producer = "RunSpatialNeighborhood",
        backend_id = if (identical(method, "observed")) "core" else "spicyr"
      )
    )
  }

  log_message(
    "Spatial neighborhood analysis completed ({.val {method}})",
    message_type = "success",
    verbose = verbose
  )
  srt
}

#' @title Spatial neighborhood plot
#'
#' @description
#' Visualize results produced by [RunSpatialNeighborhood()] using scop spatial,
#' statistical, and network plotting conventions.
#'
#' @md
#' @inheritParams RunSpatialNeighborhood
#' @inheritParams SpatialSpotPlot
#' @param method Stored neighborhood method to plot. If `NULL`, the active
#' method is used.
#' @param plot_type Plot type. One of `"heatmap"`, `"network"`, `"stat"`, or
#' `"spatial"`.
#' @param comparison,condition Optional filters for stored result tables.
#' @param value Column used as the plotted effect value.
#' @param FDR_threshold FDR cutoff used to mark significant pairs.
#' @param top_n Number of pairs to show for network and statistic plots.
#' @param layout Network layout.
#' @param edge_size Network edge size range.
#' @param cols.enriched,cols.depleted,cols.ns Colors for direction categories.
#' @param pair Pair to visualize for `plot_type = "spatial"`, either
#' `"from|to"` or a length-2 character vector.
#' @param seed Random seed used by layouts and jittered plot layers.
#'
#' @return A `ggplot`, `patchwork`, or list of `ggplot` objects.
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' spatial <- RunSpatialNeighborhood(
#'   spatial,
#'   group.by = "coda_label",
#'   coord.cols = c("x", "y"),
#'   k = 4,
#'   verbose = FALSE
#' )
#'
#' SpatialNeighborhoodPlot(spatial, plot_type = "heatmap")
#' SpatialNeighborhoodPlot(spatial, plot_type = "network", top_n = 12)
#' SpatialNeighborhoodPlot(spatial, plot_type = "stat", top_n = 12)
#' SpatialNeighborhoodPlot(
#'   spatial,
#'   plot_type = "spatial",
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
SpatialNeighborhoodPlot <- function(
  srt,
  method = NULL,
  plot_type = c("heatmap", "network", "stat", "spatial"),
  comparison = NULL,
  condition = NULL,
  value = c("estimate", "fraction", "count"),
  FDR_threshold = 0.05,
  top_n = 30,
  layout = c("fr", "nicely", "kk", "circle", "mds"),
  edge_size = c(0.4, 2),
  cols.enriched = "#d7301f",
  cols.depleted = "#2b8cbe",
  cols.ns = "grey75",
  pair = NULL,
  image = NULL,
  overlay_image = TRUE,
  coord.cols = c("col", "row"),
  split.by = NULL,
  palette = "RdBu",
  palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  seed = 11,
  verbose = TRUE,
  ...
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  plot_type <- match.arg(plot_type)
  value <- match.arg(value)
  layout <- match.arg(layout)
  bundle <- spatial_neighborhood_get_bundle(srt = srt, method = method)

  if (identical(plot_type, "spatial")) {
    return(spatial_neighborhood_spatial_plot(
      srt = srt,
      bundle = bundle,
      pair = pair,
      image = image,
      overlay_image = overlay_image,
      coord.cols = coord.cols,
      split.by = split.by,
      palette = palette,
      palcolor = palcolor,
      legend.position = legend.position,
      legend.direction = legend.direction,
      legend.title = legend.title,
      theme_use = theme_use,
      theme_args = theme_args,
      combine = combine,
      nrow = nrow,
      ncol = ncol,
      byrow = byrow,
      verbose = verbose,
      ...
    ))
  }

  df <- spatial_neighborhood_filter_pair_table(
    bundle$pair_table,
    comparison = comparison,
    condition = condition
  )
  if (nrow(df) == 0L) {
    return(scop_spatial_empty_plot(
      "No spatial neighborhood records remain after filtering",
      title = legend.title %||% value,
      theme_use = theme_use,
      theme_args = theme_args
    ))
  }
  df <- spatial_neighborhood_prepare_plot_table(
    df = df,
    value = value,
    FDR_threshold = FDR_threshold
  )

  switch(plot_type,
    heatmap = spatial_neighborhood_heatmap_plot(
      df = df,
      value = value,
      palette = palette,
      palcolor = palcolor,
      legend.position = legend.position,
      legend.direction = legend.direction,
      legend.title = legend.title %||% value,
      theme_use = theme_use,
      theme_args = theme_args
    ),
    network = spatial_neighborhood_network_plot(
      df = df,
      value = value,
      FDR_threshold = FDR_threshold,
      top_n = top_n,
      layout = layout,
      edge_size = edge_size,
      cols.enriched = cols.enriched,
      cols.depleted = cols.depleted,
      cols.ns = cols.ns,
      legend.position = legend.position,
      legend.direction = legend.direction,
      legend.title = legend.title,
      theme_use = theme_use,
      theme_args = theme_args,
      seed = seed
    ),
    stat = spatial_neighborhood_stat_plot(
      df = df,
      value = value,
      top_n = top_n,
      cols.enriched = cols.enriched,
      cols.depleted = cols.depleted,
      cols.ns = cols.ns,
      legend.position = legend.position,
      legend.direction = legend.direction,
      legend.title = legend.title,
      theme_use = theme_use,
      theme_args = theme_args
    )
  )
}

spatial_neighborhood_input <- function(
  srt,
  group.by,
  assay = NULL,
  layer = "data",
  coord.cols = c("col", "row"),
  image = NULL,
  coordinate_space = c("legacy_display", "raw"),
  sample.by = NULL,
  split.by = NULL,
  subject.by = NULL,
  features = NULL
) {
  meta <- srt@meta.data
  needed <- c(group.by, sample.by, split.by, subject.by)
  missing <- setdiff(needed, colnames(meta))
  if (length(missing) > 0L) {
    log_message(
      "Missing metadata column{?s}: {.val {missing}}",
      message_type = "error"
    )
  }

  coords <- spatial_analysis_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    coordinate_space = coordinate_space
  )$data
  cells <- intersect(rownames(coords), rownames(meta))
  if (length(cells) == 0L) {
    log_message("No cells or spots match spatial coordinates", message_type = "error")
  }

  out <- data.frame(
    cell = cells,
    x = suppressWarnings(as.numeric(coords[cells, "x", drop = TRUE])),
    y = suppressWarnings(as.numeric(coords[cells, "y", drop = TRUE])),
    group = as.character(meta[cells, group.by, drop = TRUE]),
    sample = if (is.null(sample.by)) "sample1" else as.character(meta[cells, sample.by, drop = TRUE]),
    condition = if (is.null(split.by)) "all" else as.character(meta[cells, split.by, drop = TRUE]),
    subject = if (is.null(subject.by)) {
      if (is.null(sample.by)) "sample1" else as.character(meta[cells, sample.by, drop = TRUE])
    } else {
      as.character(meta[cells, subject.by, drop = TRUE])
    },
    stringsAsFactors = FALSE
  )
  keep <- is.finite(out$x) & is.finite(out$y) &
    !is.na(out$group) & nzchar(out$group) &
    !is.na(out$sample) & nzchar(out$sample) &
    !is.na(out$condition) & nzchar(out$condition)
  out <- out[keep, , drop = FALSE]
  if (nrow(out) < 2L) {
    log_message(
      "At least two cells or spots with finite coordinates and labels are required",
      message_type = "error"
    )
  }
  rownames(out) <- out$cell

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!is.null(features)) {
    expr <- GetAssayData5(srt, assay = assay, layer = layer)
    missing_features <- setdiff(features, rownames(expr))
    if (length(missing_features) > 0L) {
      log_message(
        "Feature{?s} not found in assay {.val {assay}}: {.val {missing_features}}",
        message_type = "error"
      )
    }
    expr <- as.matrix(expr[features, out$cell, drop = FALSE])
    for (feature in features) {
      out[[feature]] <- as.numeric(expr[feature, out$cell, drop = TRUE])
    }
  }

  list(cells = out, assay = assay)
}

spatial_neighborhood_filter_pairs <- function(df, from = NULL, to = NULL) {
  if (!is.null(from)) {
    df <- df[df$from %in% from, , drop = FALSE]
  }
  if (!is.null(to)) {
    df <- df[df$to %in% to, , drop = FALSE]
  }
  if (nrow(df) == 0L) {
    log_message(
      "No spatial neighborhood pairs remain after {.arg from}/{.arg to} filtering",
      message_type = "error"
    )
  }
  df
}

spatial_neighborhood_filter_edges <- function(df, from = NULL, to = NULL) {
  if (!is.null(from)) {
    df <- df[df$from %in% from, , drop = FALSE]
  }
  if (!is.null(to)) {
    df <- df[df$to %in% to, , drop = FALSE]
  }
  df
}

spatial_neighborhood_observed_pairs <- function(cells, radius = NULL, k = NULL) {
  if (is.null(k) && is.null(radius)) {
    k <- 6L
  }
  if (!is.null(k)) {
    k <- as.integer(k[1L])
    if (is.na(k) || k < 1L) {
      log_message("{.arg k} must be a positive integer", message_type = "error")
    }
  }
  if (!is.null(radius)) {
    radius <- suppressWarnings(as.numeric(radius[1L]))
    if (is.na(radius) || radius <= 0) {
      log_message("{.arg radius} must be a positive number", message_type = "error")
    }
  }

  edge_list <- lapply(split(cells, cells$sample), function(sample_cells) {
    if (nrow(sample_cells) < 2L) {
      return(NULL)
    }
    coords <- sample_cells[, c("x", "y"), drop = FALSE]
    coords$cell_id <- sample_cells$cell
    if (!is.null(radius)) {
      graph <- spatial_graph_compute(
        coords = coords,
        method = "radius",
        radius = radius,
        directed = TRUE,
        weight = "binary"
      )
      if (!is.null(k) && nrow(graph$edges) > 0L) {
        keep <- unlist(lapply(
          split(seq_len(nrow(graph$edges)), graph$edges$from),
          function(i) utils::head(i, k)
        ), use.names = FALSE)
        graph$edges <- graph$edges[keep, , drop = FALSE]
      }
    } else {
      graph <- spatial_graph_compute(
        coords = coords,
        method = "knn",
        k = min(k, nrow(coords) - 1L),
        directed = TRUE,
        weight = "binary"
      )
    }
    if (nrow(graph$edges) == 0L) return(NULL)
    from_idx <- graph$edges$from
    to_idx <- graph$edges$to
    data.frame(
      cell = sample_cells$cell[from_idx],
      neighbor = sample_cells$cell[to_idx],
      from = sample_cells$group[from_idx],
      to = sample_cells$group[to_idx],
      sample = sample_cells$sample[from_idx],
      condition = sample_cells$condition[from_idx],
      subject = sample_cells$subject[from_idx],
      distance = graph$edges$distance,
      stringsAsFactors = FALSE
    )
  })
  edge_table <- do.call(rbind, edge_list)
  if (is.null(edge_table) || nrow(edge_table) == 0L) {
    log_message(
      "No spatial neighbor pairs were found with the selected {.arg radius}/{.arg k}",
      message_type = "error"
    )
  }
  rownames(edge_table) <- NULL

  count_df <- stats::aggregate(
    distance ~ sample + condition + subject + from + to,
    data = edge_table,
    FUN = length
  )
  colnames(count_df)[colnames(count_df) == "distance"] <- "count"
  total_df <- stats::aggregate(
    count ~ sample + condition,
    data = count_df,
    FUN = sum
  )
  colnames(total_df)[colnames(total_df) == "count"] <- "total"
  count_df <- merge(count_df, total_df, by = c("sample", "condition"), all.x = TRUE)
  count_df$fraction <- count_df$count / pmax(count_df$total, 1)
  count_df$method <- "observed"
  count_df$comparison <- count_df$condition
  count_df$estimate <- count_df$fraction
  count_df$statistic <- NA_real_
  count_df$pval <- NA_real_
  count_df$FDR <- NA_real_
  count_df$direction <- "observed"
  count_df <- count_df[, c(
    "method", "comparison", "condition", "from", "to", "estimate",
    "statistic", "pval", "FDR", "direction", "sample", "subject",
    "count", "total", "fraction"
  ), drop = FALSE]

  list(pair_table = count_df, edge_table = edge_table)
}

spatial_neighborhood_run_spicyr <- function(
  cells,
  group.by,
  sample.by = NULL,
  split.by = NULL,
  subject.by = NULL,
  verbose = TRUE,
  ...
) {
  check_r("spicyR", verbose = FALSE)
  spicy <- get_namespace_fun("spicyR", "spicy")
  backend_cells <- cells
  backend_cells[[".scop_cell_type"]] <- backend_cells$group
  backend_cells[[".scop_image_id"]] <- backend_cells$sample
  backend_cells[[".scop_condition"]] <- backend_cells$condition
  backend_cells[[".scop_subject"]] <- ifelse(
    is.na(backend_cells$subject),
    backend_cells$sample,
    backend_cells$subject
  )

  args <- list(
    cells = backend_cells,
    condition = ".scop_condition",
    subject = ".scop_subject",
    cellType = ".scop_cell_type",
    imageID = ".scop_image_id",
    spatialCoords = c("x", "y")
  )
  dots <- list(...)
  args <- utils::modifyList(args, dots)
  raw <- invoke_fun(spicy, args)
  list(raw = raw, table = spatial_neighborhood_as_data_frame(raw))
}

spatial_neighborhood_standardize_pair_table <- function(backend, observed, method) {
  raw_df <- backend$table
  if (is.null(raw_df) || nrow(raw_df) == 0L) {
    out <- observed
    out$method <- method
    return(out)
  }

  df <- as.data.frame(raw_df, stringsAsFactors = FALSE, check.names = FALSE)
  from_col <- spatial_neighborhood_first_col(df, c("from", "cellType1", "cell_type1", "source", "sender"))
  to_col <- spatial_neighborhood_first_col(df, c("to", "cellType2", "cell_type2", "target", "receiver"))
  pair_col <- spatial_neighborhood_first_col(df, c("pair", "testPair", "cellTypePair", "contrast"))
  if ((is.null(from_col) || is.null(to_col)) && !is.null(pair_col)) {
    parts <- strsplit(as.character(df[[pair_col]]), "[|:~_]", perl = TRUE)
    df$from <- vapply(parts, function(x) x[1L] %||% NA_character_, character(1))
    df$to <- vapply(parts, function(x) x[2L] %||% NA_character_, character(1))
    from_col <- "from"
    to_col <- "to"
  }
  if (is.null(from_col) || is.null(to_col)) {
    obs_summary <- stats::aggregate(
      cbind(count, total, fraction) ~ from + to,
      data = observed,
      FUN = mean
    )
    obs_summary$method <- method
    obs_summary$comparison <- "spicyR"
    obs_summary$condition <- NA_character_
    obs_summary$estimate <- obs_summary$fraction
    obs_summary$statistic <- NA_real_
    obs_summary$pval <- NA_real_
    obs_summary$FDR <- NA_real_
    obs_summary$direction <- "observed"
    obs_summary$sample <- NA_character_
    obs_summary$subject <- NA_character_
    return(obs_summary[, colnames(observed), drop = FALSE])
  }

  estimate_col <- spatial_neighborhood_first_col(df, c("estimate", "coefficient", "coef", "logFC", "effect"))
  statistic_col <- spatial_neighborhood_first_col(df, c("statistic", "t", "t.value", "z", "z.value"))
  pval_col <- spatial_neighborhood_first_col(df, c("pval", "p.value", "p_value", "PValue", "P.Value"))
  fdr_col <- spatial_neighborhood_first_col(df, c("FDR", "fdr", "adj.P.Val", "q.value", "qval"))
  comparison_col <- spatial_neighborhood_first_col(df, c("comparison", "contrast", "condition"))

  out <- data.frame(
    method = method,
    comparison = if (is.null(comparison_col)) method else as.character(df[[comparison_col]]),
    condition = if (is.null(comparison_col)) NA_character_ else as.character(df[[comparison_col]]),
    from = as.character(df[[from_col]]),
    to = as.character(df[[to_col]]),
    estimate = spatial_neighborhood_numeric_col(df, estimate_col),
    statistic = spatial_neighborhood_numeric_col(df, statistic_col),
    pval = spatial_neighborhood_numeric_col(df, pval_col),
    FDR = spatial_neighborhood_numeric_col(df, fdr_col),
    direction = NA_character_,
    sample = NA_character_,
    subject = NA_character_,
    count = NA_real_,
    total = NA_real_,
    fraction = NA_real_,
    stringsAsFactors = FALSE
  )
  if (all(is.na(out$FDR)) && any(!is.na(out$pval))) {
    out$FDR <- stats::p.adjust(out$pval, method = "fdr")
  }
  out$direction <- ifelse(
    !is.na(out$FDR) & out$FDR <= 0.05 & out$estimate > 0,
    "enriched",
    ifelse(
      !is.na(out$FDR) & out$FDR <= 0.05 & out$estimate < 0,
      "depleted",
      "ns"
    )
  )
  out
}

spatial_neighborhood_long_table <- function(pair_table, observed, method) {
  pair_table$pair <- paste(pair_table$from, pair_table$to, sep = "|")
  pair_table$record_type <- "pair"
  pair_table
}

spatial_neighborhood_get_bundle <- function(srt, method = NULL, tool_name = "SpatialNeighborhood") {
  store <- srt@tools[[tool_name]]
  if (is.null(store)) {
    log_message(
      "Cannot find SpatialNeighborhood results. Run {.fn RunSpatialNeighborhood} first",
      message_type = "error"
    )
  }
  method <- method %||% store$active_method
  if (is.null(method) || !method %in% names(store$methods)) {
    log_message(
      "SpatialNeighborhood method {.val {method}} was not found",
      message_type = "error"
    )
  }
  store$methods[[method]]
}

spatial_neighborhood_filter_pair_table <- function(df, comparison = NULL, condition = NULL) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (!is.null(comparison) && "comparison" %in% colnames(df)) {
    df <- df[df$comparison %in% comparison, , drop = FALSE]
  }
  if (!is.null(condition) && "condition" %in% colnames(df)) {
    df <- df[df$condition %in% condition, , drop = FALSE]
  }
  df
}

spatial_neighborhood_prepare_plot_table <- function(df, value, FDR_threshold) {
  if (!value %in% colnames(df)) {
    log_message(
      "{.arg value} {.val {value}} is not available in SpatialNeighborhood results",
      message_type = "error"
    )
  }
  df[[value]] <- suppressWarnings(as.numeric(df[[value]]))
  df$direction <- ifelse(
    !is.na(df$FDR) & df$FDR <= FDR_threshold & df[[value]] > 0,
    "enriched",
    ifelse(
      !is.na(df$FDR) & df$FDR <= FDR_threshold & df[[value]] < 0,
      "depleted",
      ifelse(is.na(df$direction) | !nzchar(df$direction), "ns", df$direction)
    )
  )
  df$pair <- paste(df$from, df$to, sep = "|")
  df
}

spatial_neighborhood_heatmap_plot <- function(
  df,
  value,
  palette = "RdBu",
  palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = NULL,
  theme_use = "theme_scop",
  theme_args = list()
) {
  cols <- palette_colors(type = "continuous", palette = palette, palcolor = palcolor)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$to, y = .data$from, fill = .data[[value]])) +
    do.call(ggplot2::geom_tile, c(list(color = "white"), spatial_neighborhood_linewidth_arg(0.2))) +
    ggplot2::scale_fill_gradientn(colors = cols, na.value = "grey90") +
    ggplot2::labs(x = "To", y = "From", fill = legend.title %||% value) +
    spatial_neighborhood_theme(theme_use, theme_args) +
    ggplot2::theme(
      legend.position = legend.position,
      legend.direction = legend.direction,
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

spatial_neighborhood_network_plot <- function(
  df,
  value,
  FDR_threshold,
  top_n,
  layout,
  edge_size,
  cols.enriched,
  cols.depleted,
  cols.ns,
  legend.position,
  legend.direction,
  legend.title,
  theme_use,
  theme_args,
  seed
) {
  df <- df[order(abs(df[[value]]), decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  if (!is.null(top_n) && is.finite(top_n) && nrow(df) > top_n) {
    df <- df[seq_len(top_n), , drop = FALSE]
  }
  edges <- data.frame(
    from = df$from,
    to = df$to,
    weight = df[[value]],
    direction = factor(df$direction, levels = c("enriched", "depleted", "ns", "observed")),
    stringsAsFactors = FALSE
  )
  pdata <- spatial_neighborhood_network_plot_data(
    edges = edges,
    layout = layout,
    seed = seed
  )
  edge_plot <- pdata$edge_plot
  node_plot <- pdata$nodes
  edge_plot$abs_weight <- abs(suppressWarnings(as.numeric(edge_plot$weight)))
  if (all(!is.finite(edge_plot$abs_weight))) {
    edge_plot$abs_weight <- 1
  }
  edge_aes <- ggplot2::aes(
    x = .data$x,
    y = .data$y,
    xend = .data$x_end,
    yend = .data$y_end,
    color = .data$direction
  )
  edge_aes[[spatial_neighborhood_linewidth_name()]] <- rlang::expr(.data$abs_weight)
  edge_scale <- if (identical(spatial_neighborhood_linewidth_name(), "linewidth")) {
    ggplot2::scale_linewidth_continuous(range = edge_size, guide = "none")
  } else {
    ggplot2::scale_size_continuous(range = edge_size, guide = "none")
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_curve(
      data = edge_plot,
      edge_aes,
      curvature = 0.18,
      alpha = 0.7,
      arrow = grid::arrow(type = "closed", length = grid::unit(0.018, "npc")),
      lineend = "round"
    ) +
    ggplot2::geom_point(
      data = node_plot,
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 21,
      size = 5,
      fill = "white",
      color = "grey20",
      stroke = 0.3
    ) +
    ggplot2::geom_text(
      data = node_plot,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$name),
      size = 3
    ) +
    ggplot2::scale_color_manual(
      values = c(
        enriched = cols.enriched,
        depleted = cols.depleted,
        ns = cols.ns,
        observed = cols.ns
      ),
      drop = FALSE
    ) +
    edge_scale +
    ggplot2::coord_equal() +
    ggplot2::labs(x = NULL, y = NULL, color = legend.title %||% "Direction") +
    spatial_neighborhood_theme(theme_use, theme_args) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = legend.position,
      legend.direction = legend.direction
    )
  p
}

spatial_neighborhood_stat_plot <- function(
  df,
  value,
  top_n,
  cols.enriched,
  cols.depleted,
  cols.ns,
  legend.position,
  legend.direction,
  legend.title,
  theme_use,
  theme_args
) {
  df <- df[order(abs(df[[value]]), decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  if (!is.null(top_n) && is.finite(top_n) && nrow(df) > top_n) {
    df <- df[seq_len(top_n), , drop = FALSE]
  }
  df$pair <- factor(df$pair, levels = rev(unique(df$pair)))
  ggplot2::ggplot(df, ggplot2::aes(x = .data$pair, y = .data[[value]], fill = .data$direction)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = c(
        enriched = cols.enriched,
        depleted = cols.depleted,
        ns = cols.ns,
        observed = cols.ns
      ),
      drop = FALSE
    ) +
    ggplot2::labs(x = NULL, y = value, fill = legend.title %||% "Direction") +
    spatial_neighborhood_theme(theme_use, theme_args) +
    ggplot2::theme(
      legend.position = legend.position,
      legend.direction = legend.direction
    )
}

spatial_neighborhood_spatial_plot <- function(
  srt,
  bundle,
  pair = NULL,
  image = NULL,
  overlay_image = TRUE,
  coord.cols = c("col", "row"),
  split.by = NULL,
  palette = "RdBu",
  palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  verbose = TRUE,
  ...
) {
  edges <- bundle$edge_table
  if (is.null(edges) || nrow(edges) == 0L) {
    log_message(
      "Spatial edge table is not available for spatial neighborhood plotting",
      message_type = "error"
    )
  }
  pair_use <- spatial_neighborhood_resolve_pair(pair, edges)
  hit <- edges$from == pair_use[1L] & edges$to == pair_use[2L]
  score <- tabulate(match(edges$cell[hit], colnames(srt)), nbins = ncol(srt))
  names(score) <- colnames(srt)
  theme_use <- theme_use %||% "theme_scop"
  SpatialSpotPlot(
    srt = srt,
    values = score,
    image = image,
    overlay_image = overlay_image,
    coord.cols = coord.cols,
    split.by = split.by,
    palette = palette,
    palcolor = palcolor,
    legend.position = legend.position,
    legend.direction = legend.direction,
    legend.title = legend.title %||% paste(pair_use, collapse = " -> "),
    theme_use = theme_use,
    theme_args = theme_args,
    combine = combine,
    nrow = nrow,
    ncol = ncol,
    byrow = byrow,
    verbose = verbose,
    ...
  )
}

spatial_neighborhood_resolve_pair <- function(pair, edges) {
  if (is.null(pair)) {
    counts <- sort(table(paste(edges$from, edges$to, sep = "|")), decreasing = TRUE)
    pair <- names(counts)[1L]
  }
  if (length(pair) == 1L) {
    pair <- strsplit(as.character(pair), "\\|", perl = TRUE)[[1]]
  }
  if (length(pair) != 2L) {
    log_message(
      "{.arg pair} must be a length-2 character vector or a single {.val from|to} string",
      message_type = "error"
    )
  }
  as.character(pair)
}

spatial_neighborhood_network_plot_data <- function(edges, layout = "fr", seed = 11) {
  check_r("igraph", verbose = FALSE)
  edges <- edges[!is.na(edges$from) & nzchar(edges$from) & !is.na(edges$to) & nzchar(edges$to), , drop = FALSE]
  if (nrow(edges) == 0L) {
    log_message("No spatial neighborhood network edges are available", message_type = "error")
  }
  nodes <- data.frame(
    name = unique(c(edges$from, edges$to)),
    stringsAsFactors = FALSE
  )
  graph <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
  xy <- spatial_neighborhood_igraph_layout(graph, layout = layout, seed = seed)
  nodes <- igraph::as_data_frame(graph, what = "vertices")
  nodes$x <- xy[match(nodes$name, rownames(xy)), 1]
  nodes$y <- xy[match(nodes$name, rownames(xy)), 2]
  edges_out <- igraph::as_data_frame(graph, what = "edges")
  from_idx <- match(edges_out$from, nodes$name)
  to_idx <- match(edges_out$to, nodes$name)
  keep <- !is.na(from_idx) & !is.na(to_idx)
  edge_plot <- edges_out[keep, , drop = FALSE]
  from_idx <- from_idx[keep]
  to_idx <- to_idx[keep]
  edge_plot$x <- nodes$x[from_idx]
  edge_plot$y <- nodes$y[from_idx]
  edge_plot$x_end <- nodes$x[to_idx]
  edge_plot$y_end <- nodes$y[to_idx]
  list(nodes = nodes, edges = edges_out, edge_plot = edge_plot)
}

spatial_neighborhood_igraph_layout <- function(graph, layout = "fr", seed = 11) {
  n_vertices <- igraph::vcount(graph)
  if (n_vertices == 1L) {
    xy <- matrix(c(0, 0), ncol = 2)
    rownames(xy) <- igraph::V(graph)$name
    colnames(xy) <- c("x", "y")
    return(xy)
  }
  weights <- NULL
  if ("weight" %in% igraph::edge_attr_names(graph)) {
    weights <- abs(suppressWarnings(as.numeric(igraph::E(graph)$weight)))
    weights[!is.finite(weights) | weights <= 0] <- 1
  }
  set.seed(seed)
  xy <- tryCatch(
    switch(layout,
      fr = igraph::layout_with_fr(graph, weights = weights),
      nicely = igraph::layout_nicely(graph),
      kk = igraph::layout_with_kk(graph, weights = weights),
      circle = igraph::layout_in_circle(graph),
      mds = igraph::layout_with_mds(graph),
      igraph::layout_nicely(graph)
    ),
    error = function(e) igraph::layout_nicely(graph)
  )
  xy <- as.matrix(xy)
  if (ncol(xy) < 2L) {
    xy <- cbind(xy[, 1], 0)
  }
  xy <- xy[, seq_len(2), drop = FALSE]
  rownames(xy) <- igraph::V(graph)$name
  colnames(xy) <- c("x", "y")
  xy
}

spatial_neighborhood_as_data_frame <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.data.frame(x)) {
    return(x)
  }
  if (is.list(x)) {
    data_items <- x[vapply(x, is.data.frame, logical(1))]
    if (length(data_items) > 0L) {
      return(data_items[[1L]])
    }
  }
  tryCatch(as.data.frame(x), error = function(e) NULL)
}

spatial_neighborhood_first_col <- function(df, candidates) {
  hit <- candidates[tolower(candidates) %in% tolower(colnames(df))]
  if (length(hit) == 0L) {
    return(NULL)
  }
  colnames(df)[match(tolower(hit[1L]), tolower(colnames(df)))]
}

spatial_neighborhood_numeric_col <- function(df, col) {
  if (is.null(col)) {
    return(rep(NA_real_, nrow(df)))
  }
  suppressWarnings(as.numeric(df[[col]]))
}

spatial_neighborhood_theme <- function(theme_use = "theme_scop", theme_args = list()) {
  if (is.null(theme_use)) {
    return(NULL)
  }
  tryCatch(
    {
      if (is.function(theme_use)) {
        do.call(theme_use, theme_args)
      } else {
        do.call(theme_use, theme_args)
      }
    },
    error = function(e) ggplot2::theme_bw()
  )
}

spatial_neighborhood_linewidth_name <- function() {
  if (utils::packageVersion("ggplot2") >= "3.4.0") {
    "linewidth"
  } else {
    "size"
  }
}

spatial_neighborhood_linewidth_arg <- function(value) {
  stats::setNames(list(value), spatial_neighborhood_linewidth_name())
}

spatial_neighborhood_assert_string <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    log_message(
      "{.arg {arg}} must be a non-empty character string",
      message_type = "error"
    )
  }
}
