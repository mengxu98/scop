# Pure coordinate and graph primitives for spatial analyses.

spatial_coords_raw <- function(
  srt,
  image = NULL,
  coord.cols = c("col", "row"),
  image_policy = "strict"
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  image_policy <- match.arg(image_policy, c("strict", "legacy_first"))
  images <- tryCatch(SeuratObject::Images(srt), error = function(e) character())
  selected_image <- image
  if (length(images) == 0L) {
    if (!is.null(selected_image)) {
      log_message("{.arg image} was supplied but {.arg srt} has no spatial images", message_type = "error")
    }
    coords <- scop_spatial_metadata_coords(srt, coord.cols = coord.cols)
    cells <- rownames(coords)
    transform <- list(
      scale = 1,
      y_flip = FALSE,
      image_width = NA_real_,
      image_height = NA_real_,
      raw_x_col = coord.cols[[1L]],
      raw_y_col = coord.cols[[2L]]
    )
  } else {
    if (is.null(selected_image)) {
      if (length(images) > 1L && identical(image_policy, "strict")) {
        log_message(
          "Multiple spatial images are available; select one with {.arg image}: {.val {images}}",
          message_type = "error"
        )
      }
      selected_image <- images[[1L]]
    }
    if (!selected_image %in% images) {
      log_message(
        "{.arg image} {.val {selected_image}} is not present in {.cls Seurat}; available images: {.val {images}}",
        message_type = "error"
      )
    }
    raw <- as.data.frame(SeuratObject::GetTissueCoordinates(srt[[selected_image]]))
    cell_col <- if ("cell" %in% colnames(raw)) "cell" else NULL
    cells <- if (is.null(cell_col)) rownames(raw) else as.character(raw[[cell_col]])
    if (is.null(cells) || anyNA(cells) || any(!nzchar(cells)) || anyDuplicated(cells)) {
      log_message("Spatial image coordinates must have unique, non-missing cell or spot identifiers", message_type = "error")
    }
    x_col <- spatial_dim_pick_col(raw, c("x", "pxl_col_in_fullres", "imagecol"))
    y_col <- spatial_dim_pick_col(raw, c("y", "pxl_row_in_fullres", "imagerow"))
    coords <- data.frame(
      x = suppressWarnings(as.numeric(raw[[x_col]])),
      y = suppressWarnings(as.numeric(raw[[y_col]])),
      row.names = cells,
      stringsAsFactors = FALSE
    )
    scale <- tryCatch(spatial_dim_image_scale(srt[[selected_image]]), error = function(e) 1)
    image_dims <- tryCatch(dim(srt[[selected_image]]@image), error = function(e) NULL)
    transform <- list(
      scale = if (length(scale) == 1L && is.finite(scale) && scale > 0) scale else 1,
      y_flip = !is.null(image_dims) && length(image_dims) >= 2L,
      image_width = if (!is.null(image_dims)) image_dims[[2L]] else NA_real_,
      image_height = if (!is.null(image_dims)) image_dims[[1L]] else NA_real_,
      raw_x_col = x_col,
      raw_y_col = y_col
    )
  }
  if (anyDuplicated(cells)) {
    log_message("Spatial coordinates contain duplicated cell or spot identifiers", message_type = "error")
  }
  finite <- is.finite(coords$x) & is.finite(coords$y)
  if (any(!finite)) {
    log_message(
      "Drop {.val {sum(!finite)}} cell or spot coordinate{?s} with non-finite x or y values",
      message_type = "warning"
    )
  }
  coords <- coords[finite, , drop = FALSE]
  coords$cell_id <- rownames(coords)
  coords <- coords[, c("cell_id", "x", "y"), drop = FALSE]
  object_cells <- colnames(srt)
  keep <- object_cells[object_cells %in% coords$cell_id]
  coords <- coords[match(keep, coords$cell_id), , drop = FALSE]
  rownames(coords) <- coords$cell_id
  coords$image <- if (is.null(selected_image)) NA_character_ else selected_image
  coords <- coords[, c("cell_id", "x", "y", "image"), drop = FALSE]
  result <- list(
    data = coords,
    transform = transform,
    source = list(
      image = selected_image,
      coord.cols = if (is.null(selected_image)) coord.cols[1:2] else NULL,
      image_policy = image_policy,
      coordinate_space = "raw"
    )
  )
  attr(result$data, "spatial_source") <- result$source
  attr(result$data, "spatial_transform") <- result$transform
  result
}

spatial_coords_to_display <- function(raw, transform) {
  out <- as.data.frame(raw, stringsAsFactors = FALSE)
  if (!all(c("x", "y") %in% colnames(out))) {
    log_message("{.arg raw} must contain x and y columns", message_type = "error")
  }
  scale <- transform$scale %||% 1
  out$x <- as.numeric(out$x) * scale
  out$y <- as.numeric(out$y) * scale
  if (isTRUE(transform$y_flip)) {
    height <- transform$image_height
    if (length(height) != 1L || !is.finite(height)) {
      log_message("Display conversion requires a finite image height", message_type = "error")
    }
    out$y <- height - out$y
  }
  out
}

spatial_coords_to_raw <- function(display, transform) {
  out <- as.data.frame(display, stringsAsFactors = FALSE)
  if (!all(c("x", "y") %in% colnames(out))) {
    log_message("{.arg display} must contain x and y columns", message_type = "error")
  }
  scale <- transform$scale %||% 1
  if (length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    log_message("Raw conversion requires a positive finite scale", message_type = "error")
  }
  if (isTRUE(transform$y_flip)) {
    height <- transform$image_height
    if (length(height) != 1L || !is.finite(height)) {
      log_message("Raw conversion requires a finite image height", message_type = "error")
    }
    out$y <- height - as.numeric(out$y)
  }
  out$x <- as.numeric(out$x) / scale
  out$y <- as.numeric(out$y) / scale
  out
}

spatial_analysis_coords <- function(
  srt,
  image = NULL,
  coord.cols = c("col", "row"),
  coordinate_space = c("legacy_display", "raw"),
  image_policy = "strict"
) {
  coordinate_space <- match.arg(coordinate_space)
  if (identical(coordinate_space, "legacy_display")) {
    result <- list(
      data = spatial_dim_coords(
        srt = srt,
        image = image,
        coord.cols = coord.cols,
        overlay_image = FALSE
      )$data,
      transform = NULL,
      source = list(image = image, coord.cols = coord.cols[1:2], coordinate_space = coordinate_space)
    )
    attr(result$data, "spatial_source") <- result$source
    attr(result$data, "spatial_transform") <- result$transform
    return(result)
  }
  spatial_coords_raw(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    image_policy = image_policy
  )
}

#' @title Read spatial coordinates with an explicit coordinate contract
#'
#' @description
#' Return raw analysis coordinates or display coordinates together with their
#' source and reversible transform. This function does not modify the object.
#'
#' @param object A `Seurat` object.
#' @param image Optional Seurat image name.
#' @param coord.cols Metadata columns used when no image is available.
#' @param space Coordinate space to return.
#' @param image_policy Multi-image selection policy. The default requires an
#'   explicit image when more than one image is available.
#'
#' @return A plain list with `data`, `source`, and `transform` entries.
#' @export
SpatialCoordinates <- function(
  object,
  image = NULL,
  coord.cols = c("col", "row"),
  space = c("raw", "display"),
  image_policy = c("strict", "legacy_first")
) {
  space <- match.arg(space)
  image_policy <- match.arg(image_policy)
  result <- spatial_coords_raw(
    srt = object,
    image = image,
    coord.cols = coord.cols,
    image_policy = image_policy
  )
  if (identical(space, "display")) {
    result$data <- spatial_coords_to_display(result$data, result$transform)
    result$source$coordinate_space <- "display"
  }
  result
}

spatial_graph_weights <- function(distance, method, sigma = NULL) {
  switch(
    method,
    binary = rep(1, length(distance)),
    inverse_distance = 1 / (1 + distance),
    gaussian = {
      if (is.null(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
        log_message("{.arg sigma} must be positive and finite for Gaussian weights", message_type = "error")
      }
      exp(-(distance^2) / (2 * sigma^2))
    }
  )
}

spatial_graph_compute <- function(
  coords,
  method = c("knn", "radius"),
  k = 6,
  radius = NULL,
  directed = FALSE,
  weight = c("binary", "inverse_distance", "gaussian"),
  sigma = NULL
) {
  method <- match.arg(method)
  weight <- match.arg(weight)
  coords <- as.data.frame(coords, stringsAsFactors = FALSE)
  if (!all(c("x", "y") %in% colnames(coords))) {
    log_message("{.arg coords} must contain x and y columns", message_type = "error")
  }
  xy <- as.matrix(coords[, c("x", "y"), drop = FALSE])
  if (nrow(xy) < 2L || any(!is.finite(xy))) {
    log_message("At least two rows with finite spatial coordinates are required", message_type = "error")
  }
  ids <- if ("cell_id" %in% colnames(coords)) as.character(coords$cell_id) else rownames(coords)
  if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    ids <- as.character(seq_len(nrow(coords)))
  }
  check_r("BiocNeighbors", verbose = FALSE)
  kmknn_param <- get_namespace_fun("BiocNeighbors", "KmknnParam")
  if (identical(method, "knn")) {
    if (is.null(k) || length(k) != 1L) {
      log_message("{.arg k} must be one positive integer", message_type = "error")
    }
    k <- suppressWarnings(as.integer(k[[1L]]))
    if (is.na(k) || k < 1L || k >= nrow(xy)) {
      log_message("{.arg k} must be between 1 and the number of nodes minus one", message_type = "error")
    }
    find_knn <- get_namespace_fun("BiocNeighbors", "findKNN")
    nearest <- find_knn(
      xy,
      k = k,
      get.index = TRUE,
      get.distance = TRUE,
      BNPARAM = kmknn_param()
    )
    from <- rep(seq_len(nrow(xy)), each = k)
    to <- as.vector(t(nearest$index))
    distance <- as.vector(t(nearest$distance))
  } else {
    if (is.null(radius) || length(radius) != 1L) {
      log_message("{.arg radius} must be one positive finite number", message_type = "error")
    }
    radius <- suppressWarnings(as.numeric(radius[[1L]]))
    if (is.na(radius) || !is.finite(radius) || radius <= 0) {
      log_message("{.arg radius} must be one positive finite number", message_type = "error")
    }
    find_neighbors <- get_namespace_fun("BiocNeighbors", "findNeighbors")
    nearest <- find_neighbors(
      xy,
      threshold = radius,
      get.index = TRUE,
      get.distance = TRUE,
      BNPARAM = kmknn_param()
    )
    from <- rep(seq_len(nrow(xy)), lengths(nearest$index))
    to <- unlist(nearest$index, use.names = FALSE)
    distance <- unlist(nearest$distance, use.names = FALSE)
  }
  edges <- data.frame(
    from = as.integer(from),
    to = as.integer(to),
    distance = as.numeric(distance),
    stringsAsFactors = FALSE
  )
  edges <- edges[
    edges$from != edges$to & is.finite(edges$distance),
    ,
    drop = FALSE
  ]
  if (!isTRUE(directed) && nrow(edges) > 0L) {
    left <- pmin(edges$from, edges$to)
    right <- pmax(edges$from, edges$to)
    ord <- order(left, right, edges$distance, ids[edges$to])
    edges <- edges[ord, , drop = FALSE]
    key <- paste(pmin(edges$from, edges$to), pmax(edges$from, edges$to), sep = "\r")
    edges <- edges[!duplicated(key), , drop = FALSE]
    canonical_from <- pmin(edges$from, edges$to)
    canonical_to <- pmax(edges$from, edges$to)
    edges$from <- canonical_from
    edges$to <- canonical_to
  } else if (nrow(edges) > 0L) {
    edges <- edges[order(edges$from, edges$distance, ids[edges$to]), , drop = FALSE]
  }
  edges$weight <- spatial_graph_weights(edges$distance, weight, sigma = sigma)
  rownames(edges) <- NULL
  list(
    nodes = data.frame(
      node_id = seq_len(nrow(coords)),
      cell_id = ids,
      x = as.numeric(coords$x),
      y = as.numeric(coords$y),
      stringsAsFactors = FALSE
    ),
    edges = edges,
    parameters = list(
      method = method,
      k = if (identical(method, "knn")) k else NULL,
      radius = if (identical(method, "radius")) radius else NULL,
      directed = isTRUE(directed),
      weight = weight,
      sigma = sigma
    )
  )
}

spatial_boundary_validate <- function(boundaries, image = NULL) {
  if (!is.data.frame(boundaries) || nrow(boundaries) == 0L) {
    log_message("{.arg boundaries} must be a non-empty data frame", message_type = "error")
  }
  pick_col <- function(candidates, required = FALSE) {
    index <- match(tolower(candidates), tolower(colnames(boundaries)))
    index <- index[!is.na(index)]
    if (length(index) == 0L) {
      if (isTRUE(required)) {
        log_message(
          "Boundary data do not contain a required column among {.val {candidates}}",
          message_type = "error"
        )
      }
      return(NULL)
    }
    colnames(boundaries)[index[[1L]]]
  }
  cell_col <- pick_col(c("cell_id", "cell", "barcode", "id"), required = TRUE)
  x_col <- pick_col(c("x", "imagecol", "pxl_col_in_fullres"), required = TRUE)
  y_col <- pick_col(c("y", "imagerow", "pxl_row_in_fullres"), required = TRUE)
  polygon_col <- pick_col(c("polygon_id", "polygon", "poly_id"))
  ring_col <- pick_col(c("ring_id", "ring", "subgroup"))
  order_col <- pick_col(c("vertex_order", "order", "index", "vertex"))
  image_col <- pick_col(c("image", "image_id", "slice"))
  out <- boundaries
  out$cell_id <- as.character(boundaries[[cell_col]])
  out$x <- suppressWarnings(as.numeric(boundaries[[x_col]]))
  out$y <- suppressWarnings(as.numeric(boundaries[[y_col]]))
  out$polygon_id <- if (is.null(polygon_col)) out$cell_id else as.character(boundaries[[polygon_col]])
  out$ring_id <- if (is.null(ring_col)) "1" else as.character(boundaries[[ring_col]])
  out$vertex_order <- if (is.null(order_col)) seq_len(nrow(out)) else suppressWarnings(as.numeric(boundaries[[order_col]]))
  out$image <- if (!is.null(image_col)) {
    as.character(boundaries[[image_col]])
  } else if (is.null(image)) {
    NA_character_
  } else {
    image
  }
  valid <- !is.na(out$cell_id) & nzchar(out$cell_id) &
    !is.na(out$polygon_id) & nzchar(out$polygon_id) &
    !is.na(out$ring_id) & nzchar(out$ring_id) &
    is.finite(out$x) & is.finite(out$y) & is.finite(out$vertex_order)
  if (any(!valid)) {
    log_message(
      "Drop {.val {sum(!valid)}} invalid boundary vertices",
      message_type = "warning"
    )
  }
  out <- out[valid, , drop = FALSE]
  if (nrow(out) == 0L) {
    log_message("No valid segmentation boundaries remain", message_type = "error")
  }
  ring_key <- interaction(out$cell_id, out$polygon_id, out$ring_id, drop = TRUE, lex.order = TRUE)
  valid_ring <- vapply(
    split(seq_len(nrow(out)), ring_key),
    function(i) nrow(unique(out[i, c("x", "y"), drop = FALSE])) >= 3L,
    logical(1)
  )
  if (any(!valid_ring)) {
    log_message(
      "Each segmentation ring must contain at least three distinct vertices; cell or spot centers are not polygons",
      message_type = "error"
    )
  }
  out$.polygon_group <- interaction(out$cell_id, out$polygon_id, drop = TRUE, lex.order = TRUE)
  out <- out[order(out$.polygon_group, out$ring_id, out$vertex_order), , drop = FALSE]
  rownames(out) <- NULL
  out
}
