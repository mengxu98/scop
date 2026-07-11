#' @title Build a native spatial network
#'
#' @description
#' Build a k-nearest-neighbor or radius spatial network from raw Seurat spatial
#' coordinates. Results are stored as named graphs in
#' `srt@tools$SpatialNetwork`.
#'
#' @param srt A `Seurat` object.
#' @param method Network method, either `"knn"` or `"radius"`.
#' @param image Seurat image name. A single image is selected automatically;
#'   multi-image objects require an explicit value.
#' @param coord.cols Metadata columns used when the object has no image.
#' @param k Number of neighbors for `method = "knn"`.
#' @param radius Positive distance threshold for `method = "radius"`, expressed
#'   in the raw coordinate units.
#' @param graph.name Optional graph name. If `NULL`, a deterministic name is
#'   generated from the image, method, and method parameter.
#' @param overwrite Whether an existing graph with the same name may be
#'   replaced.
#' @inheritParams thisutils::log_message
#'
#' @return The input `Seurat` object with a `SpatialNetwork` result in
#'   `srt@tools`.
#'
#' @examples
#' counts <- matrix(
#'   c(3, 1, 0, 2, 0, 4, 1, 0, 2, 1, 3, 0),
#'   nrow = 3,
#'   dimnames = list(paste0("gene", 1:3), paste0("spot", 1:4))
#' )
#' srt <- SeuratObject::CreateSeuratObject(counts)
#' srt$col <- c(0, 1, 0, 1)
#' srt$row <- c(0, 0, 1, 1)
#' srt <- RunSpatialNetwork(srt, k = 2, verbose = FALSE)
#' SpatialNetworkPlot(srt)
#'
#' @export
RunSpatialNetwork <- function(
  srt,
  method = c("knn", "radius"),
  image = NULL,
  coord.cols = c("col", "row"),
  k = 6,
  radius = NULL,
  graph.name = NULL,
  overwrite = FALSE,
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  if (!is.null(image) && (!is.character(image) || length(image) != 1L || is.na(image) || !nzchar(image))) {
    log_message("{.arg image} must be one non-empty image name", message_type = "error")
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    log_message("{.arg overwrite} must be TRUE or FALSE", message_type = "error")
  }
  method <- match.arg(method)
  coord_result <- spatial_coords_raw(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    image_policy = "strict"
  )
  selected_image <- coord_result$source$image
  coords <- coord_result$data
  cells <- intersect(colnames(srt), coords$cell_id)
  coords <- coords[match(cells, coords$cell_id), , drop = FALSE]
  rownames(coords) <- coords$cell_id
  if (length(cells) < 2L) {
    log_message("At least two cells or spots with finite spatial coordinates are required", message_type = "error")
  }

  if (identical(method, "knn")) {
    k_numeric <- suppressWarnings(as.numeric(k))
    if (
      length(k_numeric) != 1L || is.na(k_numeric) || !is.finite(k_numeric) ||
        k_numeric != floor(k_numeric)
    ) {
      log_message("{.arg k} must be one positive integer", message_type = "error")
    }
    k <- as.integer(k_numeric)
    if (is.na(k) || k < 1L || k >= nrow(coords)) {
      log_message("{.arg k} must be between 1 and the number of nodes minus one", message_type = "error")
    }
  } else {
    if (is.null(radius) || length(radius) != 1L) {
      log_message("{.arg radius} must be one positive finite number", message_type = "error")
    }
    radius <- suppressWarnings(as.numeric(radius))
    if (is.na(radius) || !is.finite(radius) || radius <= 0) {
      log_message("{.arg radius} must be one positive finite number", message_type = "error")
    }
  }
  graph <- spatial_graph_compute(
    coords = coords,
    method = method,
    k = k,
    radius = radius,
    directed = FALSE,
    weight = "binary"
  )
  nodes <- graph$nodes[, c("cell_id", "x", "y"), drop = FALSE]
  nodes$image <- if (is.null(selected_image)) NA_character_ else selected_image
  edges <- graph$edges
  edges$from <- graph$nodes$cell_id[edges$from]
  edges$to <- graph$nodes$cell_id[edges$to]
  parameters <- graph$parameters[c("method", "k", "radius")]

  sanitize_name <- function(x) {
    x <- gsub("[^A-Za-z0-9]+", "_", x)
    gsub("^_+|_+$", "", x)
  }
  if (is.null(graph.name)) {
    prefix <- if (is.null(selected_image)) "" else paste0(sanitize_name(selected_image), "_")
    if (identical(method, "knn")) {
      graph.name <- paste0(prefix, "knn_k", k)
    } else {
      radius_token <- format(radius, scientific = FALSE, trim = TRUE, digits = 15)
      radius_token <- gsub("-", "m", gsub("\\.", "p", radius_token))
      graph.name <- paste0(prefix, "radius_r", radius_token)
    }
  }
  if (!is.character(graph.name) || length(graph.name) != 1L || is.na(graph.name) || !nzchar(graph.name)) {
    log_message("{.arg graph.name} must be one non-empty string", message_type = "error")
  }

  store <- srt@tools[["SpatialNetwork"]] %||% list(
    method = "SpatialNetwork",
    active_graph = NULL,
    graphs = list()
  )
  store$graphs <- store$graphs %||% list()
  if (!is.null(store$graphs[[graph.name]]) && !isTRUE(overwrite)) {
    log_message(
      "Spatial graph {.val {graph.name}} already exists; set {.arg overwrite = TRUE} to replace it",
      message_type = "error"
    )
  }
  store$graphs[[graph.name]] <- list(
    nodes = nodes,
    edges = edges,
    parameters = c(parameters, list(weight = "binary")),
    source = c(coord_result$source, list(transform = coord_result$transform))
  )
  store$method <- "SpatialNetwork"
  store$active_graph <- graph.name
  srt@tools[["SpatialNetwork"]] <- store
  log_message(
    "Stored spatial graph {.val {graph.name}} with {.val {nrow(nodes)}} nodes and {.val {nrow(edges)}} edges",
    message_type = "success",
    verbose = verbose
  )
  srt
}

#' @title Plot a native spatial network
#'
#' @description
#' Plot a graph produced by [RunSpatialNetwork()]. The graph can be read from a
#' Seurat object or supplied directly as `srt@tools$SpatialNetwork`.
#'
#' @param object Optional `Seurat` object containing the graph and metadata.
#' @param res Optional plain result list from `object@tools$SpatialNetwork`.
#' @param graph.name Stored graph name. The active graph is used when `NULL`.
#' @param group.by Node column or Seurat metadata column used for coloring.
#' @param edge.color,edge.linewidth Edge appearance.
#' @param pt.size,pt.alpha Node appearance.
#' @param palette,palcolor Palette name or explicit colors.
#' @param raster Whether to rasterize only the node layer.
#' @param raster.dpi Node rasterization resolution.
#' @param theme_use,theme_args scop theme and its arguments.
#'
#' @return A `ggplot` object.
#'
#' @export
SpatialNetworkPlot <- function(
  object = NULL,
  res = NULL,
  graph.name = NULL,
  group.by = NULL,
  edge.color = "grey80",
  edge.linewidth = 0.2,
  pt.size = NULL,
  pt.alpha = 1,
  palette = "Paired",
  palcolor = NULL,
  raster = FALSE,
  raster.dpi = 300,
  theme_use = "theme_scop",
  theme_args = list()
) {
  if (!is.null(object) && !inherits(object, "Seurat")) {
    log_message("{.arg object} must be a {.cls Seurat} object", message_type = "error")
  }
  if (!is.null(graph.name) && (!is.character(graph.name) || length(graph.name) != 1L || is.na(graph.name) || !nzchar(graph.name))) {
    log_message("{.arg graph.name} must be one non-empty string", message_type = "error")
  }
  if (!is.null(group.by) && (!is.character(group.by) || length(group.by) != 1L || is.na(group.by) || !nzchar(group.by))) {
    log_message("{.arg group.by} must be one non-empty column name", message_type = "error")
  }
  if (is.null(res)) {
    if (is.null(object)) {
      log_message("Provide either {.arg object} or {.arg res}", message_type = "error")
    }
    res <- object@tools[["SpatialNetwork"]]
  }
  if (!is.list(res) || is.null(res$graphs)) {
    log_message("{.arg res} must be a SpatialNetwork result list", message_type = "error")
  }
  graph.name <- graph.name %||% res$active_graph
  if (is.null(graph.name) || is.null(res$graphs[[graph.name]])) {
    log_message("Spatial graph {.val {graph.name}} was not found", message_type = "error")
  }
  graph <- res$graphs[[graph.name]]
  nodes <- as.data.frame(graph$nodes, stringsAsFactors = FALSE)
  edges <- as.data.frame(graph$edges, stringsAsFactors = FALSE)
  if (!all(c("cell_id", "x", "y") %in% colnames(nodes))) {
    log_message("SpatialNetwork nodes must contain cell_id, x, and y", message_type = "error")
  }

  image_info <- NULL
  plot_space <- "raw"
  if (is.null(object) && !is.null(graph$source$transform)) {
    nodes <- spatial_coords_to_display(nodes, graph$source$transform)
    plot_space <- "display"
  }
  if (!is.null(object)) {
    keep <- nodes$cell_id %in% colnames(object)
    if (any(!keep)) {
      log_message(
        "Drop {.val {sum(!keep)}} graph nodes that are absent from {.arg object}",
        message_type = "warning"
      )
    }
    nodes <- nodes[keep, , drop = FALSE]
    source_image <- graph$source$image %||% NULL
    object_images <- tryCatch(SeuratObject::Images(object), error = function(e) character())
    if (!is.null(source_image) && source_image %in% object_images) {
      display <- spatial_dim_coords(object, image = source_image, overlay_image = TRUE)
      matched <- match(nodes$cell_id, rownames(display$data))
      valid <- !is.na(matched)
      nodes <- nodes[valid, , drop = FALSE]
      matched <- matched[valid]
      nodes$x <- display$data$x[matched]
      nodes$y <- display$data$y[matched]
      image_info <- display$image
      plot_space <- "display"
    }
  }
  if (nrow(nodes) == 0L) {
    log_message("No graph nodes remain for plotting", message_type = "error")
  }

  if (!is.null(group.by)) {
    if (group.by %in% colnames(nodes)) {
      nodes$.group <- nodes[[group.by]]
    } else if (!is.null(object) && group.by %in% colnames(object@meta.data)) {
      nodes$.group <- object@meta.data[nodes$cell_id, group.by]
    } else {
      log_message(
        "{.arg group.by} {.val {group.by}} was not found in graph nodes or Seurat metadata",
        message_type = "error"
      )
    }
  }

  from_idx <- match(edges$from, nodes$cell_id)
  to_idx <- match(edges$to, nodes$cell_id)
  keep_edges <- !is.na(from_idx) & !is.na(to_idx)
  edge_dat <- edges[keep_edges, , drop = FALSE]
  edge_dat$x <- nodes$x[from_idx[keep_edges]]
  edge_dat$y <- nodes$y[from_idx[keep_edges]]
  edge_dat$xend <- nodes$x[to_idx[keep_edges]]
  edge_dat$yend <- nodes$y[to_idx[keep_edges]]

  p <- ggplot2::ggplot()
  if (!is.null(image_info)) {
    p <- p + ggplot2::annotation_raster(
      spatial_dim_raster(image_info$image, 1),
      xmin = 0,
      xmax = image_info$width,
      ymin = 0,
      ymax = image_info$height
    )
  }
  p <- p + ggplot2::geom_segment(
    data = edge_dat,
    ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
    color = edge.color,
    linewidth = edge.linewidth
  )
  pt.size <- pt.size %||% 1
  if (is.null(group.by)) {
    node_layer <- ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = "grey20",
      size = pt.size,
      alpha = pt.alpha
    )
  } else if (is.numeric(nodes$.group)) {
    node_layer <- ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = .data$x, y = .data$y, color = .data$.group),
      size = pt.size,
      alpha = pt.alpha
    )
  } else {
    node_layer <- ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = .data$x, y = .data$y, fill = .data$.group),
      shape = 21,
      color = "grey20",
      stroke = 0.1,
      size = pt.size,
      alpha = pt.alpha
    )
  }
  if (isTRUE(raster)) {
    raster.dpi <- suppressWarnings(as.numeric(raster.dpi[[1L]]))
    if (is.na(raster.dpi) || !is.finite(raster.dpi) || raster.dpi <= 0) {
      log_message("{.arg raster.dpi} must be one positive number", message_type = "error")
    }
    check_r("ggrastr", verbose = FALSE)
    rasterise <- get_namespace_fun("ggrastr", "rasterise")
    node_layer <- rasterise(node_layer, dpi = raster.dpi)
  }
  p <- p + node_layer
  if (!is.null(group.by) && is.numeric(nodes$.group)) {
    p <- p + ggplot2::scale_color_gradientn(
      colors = palette_colors(type = "continuous", palette = palette, palcolor = palcolor),
      na.value = "grey80"
    ) + ggplot2::labs(color = group.by)
  } else if (!is.null(group.by)) {
    lvls <- levels(factor(nodes$.group))
    p <- p + ggplot2::scale_fill_manual(
      values = palette_colors(lvls, palette = palette, palcolor = palcolor),
      na.value = "grey80"
    ) + ggplot2::labs(fill = group.by)
  }
  x_range <- range(nodes$x, na.rm = TRUE)
  y_range <- range(nodes$y, na.rm = TRUE)
  x_pad <- max(diff(x_range) * 0.04, .Machine$double.eps)
  y_pad <- max(diff(y_range) * 0.04, .Machine$double.eps)
  p +
    ggplot2::coord_equal(
      xlim = x_range + c(-x_pad, x_pad),
      ylim = y_range + c(-y_pad, y_pad)
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      subtitle = if (identical(plot_space, "raw")) "coordinate space: raw" else NULL
    ) +
    spatial_dim_drop_coord(do.call(theme_use, theme_args))
}
