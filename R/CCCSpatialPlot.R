# Shared spatial CCC payloads and slice-network rendering ---------------------

.ccc_spatial_plot_methods <- c("SpatialCellChat", "SpaTalk", "COMMOT")

ccc_spatial_plot_coordinates <- function(coordinates, source = list()) {
  coordinates <- as.data.frame(coordinates, stringsAsFactors = FALSE)
  required <- c("cell_id", "x", "y")
  if (!all(required %in% colnames(coordinates))) {
    log_message(
      "Spatial CCC coordinates must contain {.field cell_id}, {.field x}, and {.field y}",
      message_type = "error"
    )
  }
  if (nrow(coordinates) == 0L) {
    log_message("Spatial CCC coordinates cannot be empty", message_type = "error")
  }
  coordinates$cell_id <- as.character(coordinates$cell_id)
  if (anyNA(coordinates$cell_id) || any(!nzchar(coordinates$cell_id)) || anyDuplicated(coordinates$cell_id)) {
    log_message("Spatial CCC coordinate IDs must be unique and non-empty", message_type = "error")
  }
  for (column in c("x", "y")) {
    coordinates[[column]] <- suppressWarnings(as.numeric(coordinates[[column]]))
    if (any(!is.finite(coordinates[[column]]))) {
      log_message("Spatial CCC coordinates must be finite", message_type = "error")
    }
  }
  if (!"x_raw" %in% colnames(coordinates)) coordinates$x_raw <- coordinates$x
  if (!"y_raw" %in% colnames(coordinates)) coordinates$y_raw <- coordinates$y
  if (!"x_display" %in% colnames(coordinates)) coordinates$x_display <- coordinates$x
  if (!"y_display" %in% colnames(coordinates)) coordinates$y_display <- coordinates$y
  for (column in c("x_raw", "y_raw", "x_display", "y_display")) {
    coordinates[[column]] <- suppressWarnings(as.numeric(coordinates[[column]]))
    if (any(!is.finite(coordinates[[column]]))) {
      log_message("Spatial CCC display coordinates must be finite", message_type = "error")
    }
  }
  coordinates$image <- if ("image" %in% colnames(coordinates)) {
    as.character(coordinates$image)
  } else {
    rep(NA_character_, nrow(coordinates))
  }
  coordinates <- coordinates[, c(
    "cell_id", "x", "y", "x_raw", "y_raw", "x_display", "y_display", "image"
  ), drop = FALSE]
  rownames(coordinates) <- coordinates$cell_id
  source <- source %||% list()
  source$coordinate_space <- source$coordinate_space %||% "raw"
  source$plot_coordinate_space <- source$plot_coordinate_space %||% "display"
  source$coordinate_contract_version <- .spatial_coordinate_contract_version
  list(coordinates = coordinates, source = source)
}

ccc_spatial_plot_payload_validate <- function(payload, producer = "spatial CCC producer") {
  if (!is.list(payload)) {
    log_message("{.val {producer}} spatial plot payload must be a list", message_type = "error")
  }
  coords <- payload$coordinates
  checked <- ccc_spatial_plot_coordinates(coords, payload$source %||% list())
  labels <- payload$labels %||% NULL
  composition <- payload$composition %||% NULL
  if (!xor(!is.null(labels), !is.null(composition))) {
    log_message(
      "Spatial CCC payload must contain exactly one of {.field labels} or {.field composition}",
      message_type = "error"
    )
  }
  ids <- checked$coordinates$cell_id
  group_levels <- as.character(payload$group_levels %||% character(0))
  if (anyNA(group_levels) || any(!nzchar(group_levels)) || anyDuplicated(group_levels)) {
    log_message("Spatial CCC payload group levels must be unique and non-empty", message_type = "error")
  }
  if (!is.null(labels)) {
    if (is.null(names(labels))) {
      if (length(labels) != length(ids)) {
        log_message("Spatial CCC labels do not align with coordinate rows", message_type = "error")
      }
      names(labels) <- ids
    }
    labels <- as.character(labels[ids])
    if (anyNA(labels) || any(!nzchar(labels))) {
      log_message("Spatial CCC labels must be present for every coordinate", message_type = "error")
    }
    observed <- unique(labels)
    if (length(group_levels) == 0L) group_levels <- observed
    if (!all(observed %in% group_levels)) {
      log_message("Spatial CCC labels are absent from {.field group_levels}", message_type = "error")
    }
  }
  if (!is.null(composition)) {
    composition <- as.matrix(composition)
    storage.mode(composition) <- "double"
    if (is.null(rownames(composition)) || is.null(colnames(composition)) ||
        anyNA(rownames(composition)) || anyNA(colnames(composition)) ||
        anyDuplicated(rownames(composition)) || anyDuplicated(colnames(composition)) ||
        !identical(as.character(rownames(composition)), ids)) {
      log_message("Spatial CCC composition must align one-to-one with coordinate IDs", message_type = "error")
    }
    if (any(!is.finite(composition)) || any(composition < 0) || any(rowSums(composition) <= 0)) {
      log_message("Spatial CCC composition must contain finite non-negative positive rows", message_type = "error")
    }
    totals <- as.numeric(rowSums(composition))
    if (any(abs(totals - 1) > 1e-6)) {
      log_message("Spatial CCC composition rows must sum to one", message_type = "error")
    }
    composition <- composition / rowSums(composition)
    if (length(group_levels) == 0L) group_levels <- colnames(composition)
    if (!identical(as.character(colnames(composition)), group_levels)) {
      log_message("Spatial CCC composition columns must match {.field group_levels}", message_type = "error")
    }
  }
  if (length(group_levels) == 0L) {
    log_message("Spatial CCC payload has no cell-type groups", message_type = "error")
  }
  payload$coordinates <- checked$coordinates
  payload$source <- checked$source
  payload$group_levels <- group_levels
  payload$analysis_level <- as.character(payload$analysis_level %||% "unknown")[[1L]]
  if (!is.null(labels)) {
    payload$labels <- stats::setNames(labels, ids)
    payload$composition <- NULL
  } else {
    payload$composition <- composition
    payload$labels <- NULL
  }
  spatial_tag_coordinate_contract(payload)
}

ccc_spatial_plot_payload_from_labels <- function(
  coordinates,
  labels,
  source = list(),
  analysis_level = "cell",
  group_levels = NULL
) {
  checked <- ccc_spatial_plot_coordinates(coordinates, source)
  label_names <- names(labels)
  labels <- as.character(labels)
  if (!is.null(label_names)) names(labels) <- label_names
  ids <- checked$coordinates$cell_id
  if (!is.null(names(labels))) labels <- labels[ids]
  if (length(labels) != length(ids)) {
    log_message("Spatial CCC labels do not align with coordinate rows", message_type = "error")
  }
  if (is.null(group_levels)) group_levels <- unique(labels)
  ccc_spatial_plot_payload_validate(list(
    coordinates = checked$coordinates,
    labels = stats::setNames(labels, ids),
    composition = NULL,
    group_levels = as.character(group_levels),
    analysis_level = analysis_level,
    source = checked$source
  ))
}

ccc_spatial_plot_payload_from_composition <- function(
  coordinates,
  composition,
  source = list(),
  analysis_level = "composition",
  group_levels = NULL
) {
  checked <- ccc_spatial_plot_coordinates(coordinates, source)
  composition <- as.matrix(composition)
  storage.mode(composition) <- "double"
  ids <- checked$coordinates$cell_id
  if (is.null(rownames(composition))) {
    if (nrow(composition) != length(ids)) {
      log_message("Spatial CCC composition rows do not align with coordinate rows", message_type = "error")
    }
    rownames(composition) <- ids
  }
  if (is.null(colnames(composition))) {
    log_message("Spatial CCC composition must have cell-type column names", message_type = "error")
  }
  if (!identical(as.character(rownames(composition)), ids)) {
    missing <- setdiff(ids, rownames(composition))
    if (length(missing) > 0L) {
      log_message("Spatial CCC composition is missing selected coordinate IDs", message_type = "error")
    }
    composition <- composition[ids, , drop = FALSE]
  }
  if (is.null(group_levels)) group_levels <- colnames(composition)
  if (length(group_levels) != ncol(composition)) {
    log_message("Spatial CCC composition group levels do not match matrix columns", message_type = "error")
  }
  totals <- rowSums(composition)
  if (any(!is.finite(totals)) || any(totals <= 0)) {
    log_message("Spatial CCC composition rows must have positive finite totals", message_type = "error")
  }
  composition <- composition / totals
  colnames(composition) <- as.character(group_levels)
  ccc_spatial_plot_payload_validate(list(
    coordinates = checked$coordinates,
    labels = NULL,
    composition = composition,
    group_levels = as.character(group_levels),
    analysis_level = analysis_level,
    source = checked$source
  ))
}

ccc_spatial_plot_payload_from_legacy <- function(result, producer = "spatial producer") {
  if (is.list(result$spatial_plot)) {
    return(ccc_spatial_plot_payload_validate(result$spatial_plot, producer))
  }
  coordinates <- result$coordinates %||% NULL
  if (!is.data.frame(coordinates) || !"label" %in% colnames(coordinates)) {
    log_message(
      "Stored {producer} result lacks a spatial plot payload; rerun the producer with the current SCOP version",
      message_type = "error"
    )
  }
  ccc_spatial_plot_payload_from_labels(
    coordinates = coordinates,
    labels = coordinates$label,
    source = result$source %||% list(),
    analysis_level = result$diagnostics$interpretation %||% "unknown"
  )
}

ccc_spatial_payload_from_raw <- function(coordinates, source = list(), transform = NULL) {
  coordinates <- as.data.frame(coordinates, stringsAsFactors = FALSE)
  if (!all(c("x", "y") %in% colnames(coordinates))) {
    log_message("Raw spatial coordinates must contain x and y", message_type = "error")
  }
  display <- coordinates
  transform <- transform %||% attr(coordinates, "spatial_transform", exact = TRUE)
  if (is.list(transform)) {
    display <- tryCatch(
      spatial_coords_to_display(coordinates, transform),
      error = function(e) coordinates
    )
  }
  coordinates$x_raw <- as.numeric(coordinates$x)
  coordinates$y_raw <- as.numeric(coordinates$y)
  coordinates$x_display <- as.numeric(display$x)
  coordinates$y_display <- as.numeric(display$y)
  ccc_spatial_plot_coordinates(coordinates, source)
}

ccc_spatial_stored <- function(object, method = NULL, condition = NULL, sample = NULL) {
  method <- detect_method(object, method)
  if (identical(method, "CCC") || !method %in% .ccc_spatial_plot_methods) {
    log_message(
      "Spatial CCC plotting requires one of {.val {paste(.ccc_spatial_plot_methods, collapse = ', ')}}; {.val {method}} is not supported",
      message_type = "error"
    )
  }
  bundle <- object@tools[[method]]
  if (is.null(bundle)) log_message("{.val {method}} results are absent", message_type = "error")
  result.name <- condition %||% bundle$active_result
  if (is.null(result.name) || !nzchar(as.character(result.name))) {
    log_message("A stored spatial CCC result name is required", message_type = "error")
  }
  result <- bundle$results[[as.character(result.name)]]
  if (is.null(result)) log_message("Unknown spatial CCC result {.val {result.name}}", message_type = "error")
  if (identical(method, "SpatialCellChat")) {
    samples <- names(result)
    if (is.null(sample)) {
      if (length(samples) != 1L) {
        log_message("Multiple spatial samples are stored; select {.arg sample}: {.val {samples}}", message_type = "error")
      }
      sample <- samples[[1L]]
    }
    if (!sample %in% samples) log_message("Unknown spatial sample {.val {sample}}", message_type = "error")
    entry <- result[[sample]]
    spatial_require_coordinate_contract(entry, "RunSpatialCellChat()")
    payload <- ccc_spatial_plot_payload_from_legacy(entry, "RunSpatialCellChat()")
    table <- entry$interactions %||% entry$long_table %||% data.frame()
    return(list(
      method = method, result.name = as.character(result.name), sample = as.character(sample),
      bundle = bundle, result = entry, payload = payload, table = table
    ))
  }
  if (!is.null(sample)) {
    log_message("{.arg sample} is only supported for multi-sample SpatialCellChat results", message_type = "error")
  }
  spatial_require_coordinate_contract(result, paste0("Run", method, "()"))
  payload <- ccc_spatial_plot_payload_from_legacy(result, paste0("Run", method, "()"))
  table <- result$long_table %||% result$lr_table %||% data.frame()
  list(
    method = method, result.name = as.character(result.name), sample = NULL,
    bundle = bundle, result = result, payload = payload, table = table
  )
}

ccc_spatial_network_data <- function(
  stored,
  sender.use = NULL,
  receiver.use = NULL,
  ligand.use = NULL,
  receptor.use = NULL,
  interaction.use = NULL,
  signaling = NULL,
  pairLR.use = NULL,
  thresh = 0.05,
  value = "sum",
  edge_value = "sum",
  edge_threshold = 0,
  top_n = 20
) {
  df <- standardize_long_df(stored$table)
  if (nrow(df) == 0L) {
    log_message("The stored spatial CCC result has no communication rows", message_type = "error")
  }
  df <- filter_long_df(
    df = df,
    sender.use = sender.use,
    receiver.use = receiver.use,
    ligand.use = ligand.use,
    receptor.use = receptor.use,
    interaction.use = interaction.use,
    signaling = signaling,
    pairLR.use = pairLR.use
  )
  if (nrow(df) == 0L) {
    log_message("No spatial CCC records remain after the requested filters", message_type = "error")
  }
  df <- ccc_assign_plot_score(df, value = value)
  df <- ccc_mark_significance(df, thresh = thresh)
  df <- prepare_plot_df(df)
  pair_df <- pair_plot_df(df)
  if (is.null(pair_df) || nrow(pair_df) == 0L) {
    log_message("No sender-receiver spatial CCC pairs remain after filtering", message_type = "error")
  }
  edge_value <- as.character(edge_value)[[1L]]
  if (!edge_value %in% c("sum", "mean", "max", "count")) {
    log_message("{.arg edge_value} is not available for spatial CCC edges", message_type = "error")
  }
  pair_df$spatial_weight <- suppressWarnings(as.numeric(pair_df[[edge_value]]))
  pair_df$spatial_weight[!is.finite(pair_df$spatial_weight)] <- 0
  if (!is.numeric(edge_threshold) || length(edge_threshold) != 1L ||
      !is.finite(edge_threshold) || edge_threshold < 0) {
    log_message("{.arg edge_threshold} must be one non-negative finite value", message_type = "error")
  }
  pair_df <- pair_df[pair_df$spatial_weight >= edge_threshold, , drop = FALSE]
  if (nrow(pair_df) == 0L) {
    log_message("No spatial CCC edges remain after {.arg edge_threshold}", message_type = "error")
  }
  if (is.numeric(top_n) && length(top_n) == 1L && is.finite(top_n) && top_n > 0) {
    pair_df <- top_pairs(pair_df, top_n = as.integer(top_n), value_col = "spatial_weight")
  }
  pair_df
}

ccc_spatial_membership <- function(payload) {
  ids <- payload$coordinates$cell_id
  groups <- payload$group_levels
  if (!is.null(payload$labels)) {
    membership <- matrix(0, nrow = length(ids), ncol = length(groups), dimnames = list(ids, groups))
    membership[cbind(seq_along(ids), match(payload$labels[ids], groups))] <- 1
    return(membership)
  }
  payload$composition[ids, groups, drop = FALSE]
}

ccc_spatial_render_data <- function(payload, pair_df) {
  coords <- payload$coordinates
  membership <- ccc_spatial_membership(payload)
  groups <- payload$group_levels
  coords$group <- if (!is.null(payload$labels)) payload$labels[coords$cell_id] else NA_character_
  node_rows <- lapply(seq_along(groups), function(i) {
    weight <- membership[, i]
    total <- sum(weight, na.rm = TRUE)
    if (!is.finite(total) || total <= 0) return(NULL)
    data.frame(
      group = groups[[i]],
      x = sum(coords$x_display * weight, na.rm = TRUE) / total,
      y = sum(coords$y_display * weight, na.rm = TRUE) / total,
      abundance = total,
      stringsAsFactors = FALSE
    )
  })
  nodes <- do.call(rbind, Filter(Negate(is.null), node_rows))
  rownames(nodes) <- NULL
  lookup <- stats::setNames(seq_len(nrow(nodes)), nodes$group)
  edges <- pair_df[
    as.character(pair_df$sender) %in% nodes$group &
      as.character(pair_df$receiver) %in% nodes$group, , drop = FALSE
  ]
  if (nrow(edges) == 0L) {
    log_message("Spatial CCC communication groups do not match the stored spatial payload", message_type = "error")
  }
  edges$sender <- as.character(edges$sender)
  edges$receiver <- as.character(edges$receiver)
  edges$weight <- suppressWarnings(as.numeric(edges$spatial_weight %||% edges$sum))
  edges$weight[!is.finite(edges$weight)] <- 0
  edges$x <- nodes$x[lookup[edges$sender]]
  edges$y <- nodes$y[lookup[edges$sender]]
  edges$xend <- nodes$x[lookup[edges$receiver]]
  edges$yend <- nodes$y[lookup[edges$receiver]]
  list(coords = coords, nodes = nodes, edges = edges, membership = membership)
}

ccc_spatial_auto_radius <- function(coords) {
  span <- diff(range(c(coords$x_display, coords$y_display), finite = TRUE))
  if (!is.finite(span) || span <= 0) span <- 1
  sqrt(span^2 / max(1, nrow(coords))) * 0.45
}

ccc_spatial_network_plot <- function(
  stored,
  pair_df,
  cell_palette = "Chinese",
  cell_palcolor = NULL,
  link_palette = "Chinese",
  link_palcolor = NULL,
  edge_size = c(0.5, 1.8),
  edge_color = NULL,
  edge_alpha = 0.6,
  edge_line = "curved",
  edge_curvature = 0.2,
  directed = TRUE,
  arrow_type = "closed",
  arrow_angle = 20,
  arrow_length = grid::unit(0.02, "npc"),
  node_size = 5,
  node_alpha = 0.9,
  spot_size = 1.2,
  spot_alpha = 0.35,
  composition_display = c("pie", "dominant", "none"),
  composition_radius = NULL,
  label = TRUE,
  label.size = 3.2,
  title = NULL,
  subtitle = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = NULL,
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list()
) {
  composition_display <- match.arg(composition_display)
  render <- ccc_spatial_render_data(stored$payload, pair_df)
  coords <- render$coords
  nodes <- render$nodes
  edges <- render$edges
  groups <- stored$payload$group_levels
  cell_cols <- palette_colors(groups, palette = cell_palette, palcolor = cell_palcolor, NA_keep = TRUE)
  link_cols <- palette_colors(groups, palette = link_palette, palcolor = link_palcolor, NA_keep = TRUE)
  if (length(edge_size) != 2L || any(!is.finite(edge_size)) || any(edge_size <= 0)) {
    log_message("{.arg edge_size} must contain two positive finite values", message_type = "error")
  }
  if (length(node_size) != 1L || !is.finite(node_size) || node_size <= 0) {
    log_message("{.arg node_size} must be one positive finite value for spatial plots", message_type = "error")
  }
  if (length(spot_size) != 1L || !is.finite(spot_size) || spot_size <= 0) {
    log_message("{.arg spot_size} must be one positive finite value", message_type = "error")
  }
  p <- ggplot2::ggplot()
  if (!is.null(stored$payload$composition) && identical(composition_display, "pie")) {
    radius <- composition_radius %||% ccc_spatial_auto_radius(coords)
    if (length(radius) != 1L || !is.finite(radius) || radius <= 0) {
      log_message("{.arg composition_radius} must be one positive finite display-coordinate value", message_type = "error")
    }
    pie_rows <- list()
    row_i <- 0L
    for (i in seq_len(nrow(coords))) {
      values <- stored$payload$composition[i, groups, drop = TRUE]
      values <- values[values > 0]
      if (length(values) == 0L) next
      values <- values / sum(values)
      ends <- cumsum(values) * 2 * pi
      starts <- c(0, head(ends, -1L))
      for (j in seq_along(values)) {
        row_i <- row_i + 1L
        pie_rows[[row_i]] <- data.frame(
          x0 = coords$x_display[[i]], y0 = coords$y_display[[i]],
          r0 = 0, r = radius, start = starts[[j]], end = ends[[j]],
          group = names(values)[[j]], stringsAsFactors = FALSE
        )
      }
    }
    pies <- do.call(rbind, pie_rows)
    p <- p + ggforce::geom_arc_bar(
      data = pies,
      ggplot2::aes(
        x0 = .data[["x0"]], y0 = .data[["y0"]], r0 = .data[["r0"]],
        r = .data[["r"]], start = .data[["start"]], end = .data[["end"]],
        fill = .data[["group"]]
      ),
      color = NA, alpha = spot_alpha, inherit.aes = FALSE, show.legend = FALSE
    )
  } else if (!is.null(stored$payload$composition) && identical(composition_display, "dominant")) {
    coords$group <- groups[max.col(stored$payload$composition[coords$cell_id, groups, drop = FALSE], ties.method = "first")]
    p <- p + ggplot2::geom_point(
      data = coords,
      ggplot2::aes(x = .data[["x_display"]], y = .data[["y_display"]], color = .data[["group"]]),
      size = spot_size, alpha = spot_alpha, show.legend = FALSE
    )
  } else if (is.null(stored$payload$composition)) {
    p <- p + ggplot2::geom_point(
      data = coords,
      ggplot2::aes(x = .data[["x_display"]], y = .data[["y_display"]], color = .data[["group"]]),
      size = spot_size, alpha = spot_alpha, show.legend = FALSE
    )
  } else if (identical(composition_display, "none")) {
    p <- p + ggplot2::geom_point(
      data = coords,
      ggplot2::aes(x = .data[["x_display"]], y = .data[["y_display"]]),
      color = "grey85", size = spot_size, alpha = spot_alpha,
      show.legend = FALSE, inherit.aes = FALSE
    )
  }
  p <- p +
    ggplot2::scale_color_manual(values = cell_cols, drop = FALSE, guide = "none") +
    ggnewscale::new_scale_color()
  edge_data <- edges[edges$sender != edges$receiver, , drop = FALSE]
  self_data <- edges[edges$sender == edges$receiver, , drop = FALSE]
  edge_arrow <- if (isTRUE(directed)) grid::arrow(type = arrow_type, angle = arrow_angle, length = arrow_length) else NULL
  if (nrow(edge_data) > 0L) {
    edge_geom <- if (identical(edge_line, "straight")) ggplot2::geom_segment else ggplot2::geom_curve
    edge_args <- list(
      data = edge_data,
      mapping = ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]], xend = .data[["xend"]],
        yend = .data[["yend"]], linewidth = .data[["weight"]], color = .data[["sender"]]
      ),
      alpha = edge_alpha, arrow = edge_arrow, show.legend = FALSE, inherit.aes = FALSE
    )
    if (!is.null(edge_color)) {
      edge_args$mapping <- ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]], xend = .data[["xend"]],
        yend = .data[["yend"]], linewidth = .data[["weight"]]
      )
      edge_args$color <- edge_color
    }
    if (identical(edge_line, "curved")) edge_args$curvature <- edge_curvature
    p <- p + do.call(edge_geom, edge_args)
  }
  if (nrow(self_data) > 0L) {
    span_x <- diff(range(c(coords$x_display, coords$y_display), finite = TRUE))
    loop <- max(0.05, if (is.finite(span_x)) span_x * 0.04 else 0.1)
    self_data$xend <- self_data$x + loop
    self_data$yend <- self_data$y + loop
    self_mapping <- if (is.null(edge_color)) {
      ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]], xend = .data[["xend"]],
        yend = .data[["yend"]], linewidth = .data[["weight"]], color = .data[["sender"]]
      )
    } else {
      ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]], xend = .data[["xend"]],
        yend = .data[["yend"]], linewidth = .data[["weight"]]
      )
    }
    self_args <- list(
      data = self_data,
      mapping = self_mapping,
      curvature = 1.3, alpha = edge_alpha, arrow = edge_arrow, show.legend = FALSE, inherit.aes = FALSE
    )
    if (!is.null(edge_color)) self_args$color <- edge_color
    p <- p + do.call(ggplot2::geom_curve, self_args)
  }
  p <- p +
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]], fill = .data[["group"]],
        size = .data[["abundance"]]
      ),
      shape = 21, color = "grey20", stroke = 0.7, alpha = node_alpha, show.legend = TRUE
    ) +
    ggplot2::scale_color_manual(values = link_cols, drop = FALSE, guide = "none") +
    ggplot2::scale_fill_manual(
      values = cell_cols, drop = FALSE,
      name = legend.title %||% "Cell type",
      guide = ggplot2::guide_legend(override.aes = list(shape = 21, size = 3, color = "grey20"))
    ) +
    ggplot2::scale_size_continuous(range = c(node_size * 0.7, node_size * 1.8), guide = "none") +
    ggplot2::scale_linewidth_continuous(range = edge_size, guide = "none") +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = title %||% paste(stored$method, stored$sample %||% "", "spatial communication"),
      subtitle = subtitle,
      x = NULL, y = NULL
    ) +
    apply_plot_theme(
      theme_use = theme_use,
      theme_args = theme_args,
      allow_null = TRUE
    ) +
    ggplot2::theme(
      legend.position = legend.position,
      legend.direction = legend.direction,
      text = ggplot2::element_text(size = font.size),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
  if (isTRUE(label)) {
    p <- p + ggrepel::geom_text_repel(
      data = nodes,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]], label = .data[["group"]]),
      seed = 1, size = label.size, color = "grey15", box.padding = 0.35,
      point.padding = 0.25, min.segment.length = 0, max.overlaps = Inf,
      inherit.aes = FALSE, show.legend = FALSE
    )
  }
  p
}
