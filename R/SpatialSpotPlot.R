#' @title Spatial spot plot
#'
#' @md
#' @inheritParams CellDimPlot
#' @param group.by Metadata columns to color spots by.
#' @param features Features to color spots by. When provided, expression values
#' are read from `assay` and `layer`.
#' @param assay Assay used for `features`. If `NULL`, the default assay is used.
#' @param layer Assay layer used for `features`.
#' @param values Optional vector, matrix, or data.frame with spot-level values.
#' Row names or vector names must match spatial spot names.
#' @param plot_type Plot type. `"point"` keeps the default spot plot behavior.
#' `"pie"` draws spot-level pies from numeric metadata columns supplied to
#' `group.by` or from a numeric matrix/data.frame supplied to `values`.
#' @param plot.data Optional long-format data.frame for plotting repeated
#' spatial points, such as cell-to-spot assignments.
#' @param spot.by Column in `plot.data` containing spot names.
#' @param color.by Column in `plot.data` used to color repeated spatial points.
#' @param geom Geometry used for `plot.data`: `"point"` or `"jitter"`.
#' @param image Name of the Seurat spatial image. If `NULL`, the first image is
#' used when present.
#' @param overlay_image Whether to draw the spatial image beneath spots.
#' @param image.alpha Transparency of the spatial image.
#' @param crop Whether to crop the panel to plotted spots.
#' @param coord.cols Metadata coordinate columns used when no image is available.
#' @param flip.y Whether to reverse the y axis for metadata coordinates.
#' @param pt.size Point size.
#' @param pie.radius,pie.radius.scale Radius controls for `plot_type = "pie"`.
#' If `pie.radius` is `NULL`, the radius is estimated from spot spacing and
#' multiplied by `pie.radius.scale`.
#' @param pt.alpha Point alpha.
#' @param stroke Point border width.
#' @param bg_color Point border color.
#' @param jitter_width,jitter_height Jitter size used when `geom = "jitter"`.
#'
#' @return A `ggplot`, `patchwork`, or list of `ggplot` objects.
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' SpatialSpotPlot(
#'   visium_human_pancreas_sub,
#'   group.by = "coda_label"
#' )
#'
#' SpatialSpotPlot(
#'   visium_human_pancreas_sub,
#'   features = rownames(visium_human_pancreas_sub)[1:2],
#'   layer = "counts"
#' )
SpatialSpotPlot <- function(
  srt,
  group.by = NULL,
  features = NULL,
  assay = NULL,
  layer = "data",
  values = NULL,
  plot_type = c("point", "pie"),
  plot.data = NULL,
  spot.by = NULL,
  color.by = NULL,
  geom = c("point", "jitter"),
  image = NULL,
  overlay_image = TRUE,
  image.alpha = 1,
  crop = TRUE,
  coord.cols = c("col", "row"),
  flip.y = TRUE,
  split.by = NULL,
  cells = NULL,
  show_na = FALSE,
  pt.size = NULL,
  pie.radius = NULL,
  pie.radius.scale = 0.45,
  pt.alpha = 0.9,
  stroke = 0.1,
  jitter_width = 0.25,
  jitter_height = 0.25,
  palette = ifelse(is.null(features), "Chinese", "Spectral"),
  palcolor = NULL,
  bg_color = "grey20",
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  plot_type <- match.arg(plot_type)
  geom <- match.arg(geom)
  if (!is.null(plot.data)) {
    if (!identical(plot_type, "point")) {
      log_message(
        "{.arg plot.data} is only supported when {.arg plot_type = 'point'}",
        message_type = "error"
      )
    }
    return(spatial_dim_long_plot(
      srt = srt,
      plot.data = plot.data,
      spot.by = spot.by,
      color.by = color.by,
      geom = geom,
      image = image,
      overlay_image = overlay_image,
      image.alpha = image.alpha,
      crop = crop,
      coord.cols = coord.cols,
      flip.y = flip.y,
      split.by = split.by,
      cells = cells,
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      stroke = stroke,
      jitter_width = jitter_width,
      jitter_height = jitter_height,
      palette = palette,
      palcolor = palcolor,
      bg_color = bg_color,
      legend.position = legend.position,
      legend.direction = legend.direction,
      legend.title = legend.title,
      theme_use = theme_use,
      theme_args = theme_args,
      combine = combine,
      nrow = nrow,
      ncol = ncol,
      byrow = byrow
    ))
  }
  if (identical(plot_type, "pie")) {
    return(spatial_dim_pie_plot(
      srt = srt,
      group.by = group.by,
      values = values,
      image = image,
      overlay_image = overlay_image,
      image.alpha = image.alpha,
      crop = crop,
      coord.cols = coord.cols,
      flip.y = flip.y,
      split.by = split.by,
      cells = cells,
      pie.radius = pie.radius,
      pie.radius.scale = pie.radius.scale,
      pt.alpha = pt.alpha,
      palette = palette,
      palcolor = palcolor,
      bg_color = bg_color,
      legend.position = legend.position,
      legend.direction = legend.direction,
      legend.title = legend.title,
      theme_use = theme_use,
      theme_args = theme_args
    ))
  }
  if (is.null(group.by) && is.null(features) && is.null(values)) {
    log_message(
      "One of {.arg group.by}, {.arg features}, {.arg values}, or {.arg plot.data} must be provided",
      message_type = "error"
    )
  }

  coords <- spatial_dim_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    overlay_image = overlay_image
  )
  dat <- coords$data
  if (!is.null(cells)) {
    dat <- dat[intersect(rownames(dat), cells), , drop = FALSE]
  }
  if (nrow(dat) == 0L) {
    log_message("No spots are available for plotting", message_type = "error")
  }
  if (is.null(pt.size)) {
    pt.size <- min(3000 / nrow(dat), 2)
  }

  if (is.null(split.by)) {
    dat[[".split"]] <- factor("All")
    split.by <- ".split"
  } else {
    if (!split.by %in% colnames(srt@meta.data)) {
      log_message(
        "{.arg split.by} {.val {split.by}} is not in srt meta.data",
        message_type = "error"
      )
    }
    dat[[split.by]] <- srt@meta.data[rownames(dat), split.by, drop = TRUE]
  }

  value_items <- spatial_dim_value_items(values)
  plot_items <- c(group.by, features, value_items)
  value_type <- c(
    rep("metadata", length(group.by)),
    rep("feature", length(features)),
    rep("values", length(value_items))
  )
  names(value_type) <- plot_items
  plots <- lapply(plot_items, function(item) {
    plot_dat <- dat
    plot_dat[[".value"]] <- spatial_dim_values(
      srt = srt,
      item = item,
      type = value_type[[item]],
      assay = assay,
      layer = layer,
      values = values,
      cells = rownames(plot_dat),
      show_na = show_na
    )
    if (isFALSE(show_na)) {
      plot_dat <- plot_dat[!is.na(plot_dat[[".value"]]), , drop = FALSE]
    }
    spatial_dim_single_plot(
      plot_dat = plot_dat,
      value_col = ".value",
      value_name = item,
      split.by = split.by,
      image_info = coords$image,
      overlay_image = overlay_image,
      image.alpha = image.alpha,
      crop = crop,
      flip.y = flip.y && isFALSE(coords$uses_image),
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      stroke = stroke,
      palette = palette,
      palcolor = palcolor,
      bg_color = bg_color,
      legend.position = legend.position,
      legend.direction = legend.direction,
      legend.title = legend.title %||% item,
      theme_use = theme_use,
      theme_args = theme_args
    )
  })
  names(plots) <- plot_items

  if (isFALSE(combine)) {
    return(plots)
  }
  if (length(plots) == 1L) {
    return(plots[[1L]])
  }
  patchwork::wrap_plots(plots, nrow = nrow, ncol = ncol, byrow = byrow)
}

spatial_dim_long_plot <- function(
  srt,
  plot.data,
  spot.by,
  color.by,
  geom = c("point", "jitter"),
  image = NULL,
  overlay_image = TRUE,
  image.alpha = 1,
  crop = TRUE,
  coord.cols = c("col", "row"),
  flip.y = TRUE,
  split.by = NULL,
  cells = NULL,
  pt.size = NULL,
  pt.alpha = 0.9,
  stroke = 0.1,
  jitter_width = 0.25,
  jitter_height = 0.25,
  palette = "Chinese",
  palcolor = NULL,
  bg_color = "grey20",
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = NULL,
  theme_use = "theme_blank",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE
) {
  geom <- match.arg(geom)
  df <- as.data.frame(plot.data, check.names = FALSE)
  if (is.null(spot.by) || !spot.by %in% colnames(df)) {
    log_message(
      "{.arg spot.by} must be a column in {.arg plot.data}",
      message_type = "error"
    )
  }
  if (is.null(color.by) || !color.by %in% colnames(df)) {
    log_message(
      "{.arg color.by} must be a column in {.arg plot.data}",
      message_type = "error"
    )
  }
  if (!is.null(cells)) {
    df <- df[df[[spot.by]] %in% cells, , drop = FALSE]
  }
  coords <- spatial_dim_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    overlay_image = overlay_image
  )
  keep <- df[[spot.by]] %in% rownames(coords$data)
  df <- df[keep, , drop = FALSE]
  if (nrow(df) == 0L) {
    log_message(
      "No rows in {.arg plot.data} match spatial spots",
      message_type = "error"
    )
  }
  df$x <- coords$data[df[[spot.by]], "x"]
  df$y <- coords$data[df[[spot.by]], "y"]
  if (!is.null(split.by)) {
    if (!split.by %in% colnames(srt@meta.data)) {
      log_message(
        "{.arg split.by} {.val {split.by}} is not in srt meta.data",
        message_type = "error"
      )
    }
    df[[split.by]] <- srt@meta.data[df[[spot.by]], split.by, drop = TRUE]
  }
  if (is.null(pt.size)) {
    pt.size <- min(3000 / nrow(df), 2)
  }

  theme_obj <- spatial_dim_drop_coord(do.call(theme_use, theme_args))
  values <- df[[color.by]]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y))
  if (isTRUE(overlay_image) && !is.null(coords$image)) {
    p <- p +
      ggplot2::annotation_raster(
        spatial_dim_raster(coords$image$image, image.alpha),
        xmin = 0,
        xmax = coords$image$width,
        ymin = 0,
        ymax = coords$image$height
      )
  }

  if (is.numeric(values)) {
    cols <- palette_colors(
      type = "continuous",
      palette = palette,
      palcolor = palcolor
    )
    point_layer <- if (geom == "jitter") {
      ggplot2::geom_jitter(
        ggplot2::aes(color = .data[[color.by]]),
        width = jitter_width,
        height = jitter_height,
        size = pt.size,
        alpha = pt.alpha
      )
    } else {
      ggplot2::geom_point(
        ggplot2::aes(color = .data[[color.by]]),
        size = pt.size,
        alpha = pt.alpha
      )
    }
    p <- p +
      point_layer +
      ggplot2::scale_color_gradientn(colors = cols, na.value = "grey80") +
      ggplot2::labs(x = NULL, y = NULL, color = legend.title %||% color.by)
  } else {
    values <- as.character(values)
    df[[color.by]] <- factor(values, levels = unique(values))
    cols <- palette_colors(
      levels(df[[color.by]]),
      palette = palette,
      palcolor = palcolor
    )
    point_layer <- if (geom == "jitter") {
      ggplot2::geom_jitter(
        ggplot2::aes(fill = .data[[color.by]]),
        shape = 21,
        color = bg_color,
        stroke = stroke,
        width = jitter_width,
        height = jitter_height,
        size = pt.size,
        alpha = pt.alpha
      )
    } else {
      ggplot2::geom_point(
        ggplot2::aes(fill = .data[[color.by]]),
        shape = 21,
        color = bg_color,
        stroke = stroke,
        size = pt.size,
        alpha = pt.alpha
      )
    }
    p <- p +
      point_layer +
      ggplot2::scale_fill_manual(values = cols, na.value = "grey80") +
      ggplot2::labs(x = NULL, y = NULL, fill = legend.title %||% color.by)
  }

  p <- p +
    ggplot2::theme(
      legend.position = legend.position,
      legend.direction = legend.direction
    ) +
    theme_obj

  if (isTRUE(flip.y) && isFALSE(coords$uses_image)) {
    p <- p + ggplot2::scale_y_reverse()
  }
  if (!is.null(split.by)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", split.by)))
  }
  if (isTRUE(crop)) {
    xpad <- diff(range(df$x, na.rm = TRUE)) * 0.04
    ypad <- diff(range(df$y, na.rm = TRUE)) * 0.04
    p <- p +
      ggplot2::coord_equal(
        xlim = range(df$x, na.rm = TRUE) + c(-xpad, xpad),
        ylim = range(df$y, na.rm = TRUE) + c(-ypad, ypad)
      )
  } else {
    p <- p + ggplot2::coord_equal()
  }
  p
}

spatial_dim_pie_plot <- function(
  srt,
  group.by = NULL,
  values = NULL,
  image = NULL,
  overlay_image = TRUE,
  image.alpha = 1,
  crop = TRUE,
  coord.cols = c("col", "row"),
  flip.y = TRUE,
  split.by = NULL,
  cells = NULL,
  pie.radius = NULL,
  pie.radius.scale = 0.45,
  pt.alpha = 0.9,
  palette = "Chinese",
  palcolor = NULL,
  bg_color = "grey20",
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = NULL,
  theme_use = "theme_blank",
  theme_args = list()
) {
  check_r("scatterpie", verbose = FALSE)
  coords <- spatial_dim_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    overlay_image = overlay_image
  )
  dat <- coords$data
  if (!is.null(cells)) {
    dat <- dat[intersect(rownames(dat), cells), , drop = FALSE]
  }
  if (nrow(dat) == 0L) {
    log_message("No spots are available for plotting", message_type = "error")
  }

  mat <- spatial_dim_pie_values(
    srt = srt,
    group.by = group.by,
    values = values,
    cells = rownames(dat)
  )
  mat[!is.finite(mat) | mat < 0] <- 0
  keep <- rowSums(mat, na.rm = TRUE) > 0
  dat <- dat[keep, , drop = FALSE]
  mat <- mat[keep, , drop = FALSE]
  if (nrow(dat) == 0L) {
    log_message(
      "No spots with positive pie values are available for plotting",
      message_type = "error"
    )
  }
  mat <- sweep(mat, 1, rowSums(mat), "/")
  plot_dat <- cbind(dat, as.data.frame(mat, check.names = FALSE))

  if (!is.null(split.by)) {
    if (!split.by %in% colnames(srt@meta.data)) {
      log_message(
        "{.arg split.by} {.val {split.by}} is not in srt meta.data",
        message_type = "error"
      )
    }
    plot_dat[[split.by]] <- srt@meta.data[rownames(plot_dat), split.by, drop = TRUE]
  }

  plot_dat[[".radius"]] <- spatial_dim_pie_radius(
    coords = plot_dat[, c("x", "y"), drop = FALSE],
    radius = pie.radius,
    scale = pie.radius.scale
  )
  cols <- palette_colors(
    colnames(mat),
    palette = palette,
    palcolor = palcolor
  )
  theme_obj <- spatial_dim_drop_coord(do.call(theme_use, theme_args))
  p <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = x, y = y))
  if (isTRUE(overlay_image) && !is.null(coords$image)) {
    p <- p +
      ggplot2::annotation_raster(
        spatial_dim_raster(coords$image$image, image.alpha),
        xmin = 0,
        xmax = coords$image$width,
        ymin = 0,
        ymax = coords$image$height
      )
  }
  p <- p +
    scatterpie::geom_scatterpie(
      ggplot2::aes(x = x, y = y, r = .data[[".radius"]]),
      cols = colnames(mat),
      color = bg_color,
      alpha = pt.alpha
    ) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(x = NULL, y = NULL, fill = legend.title %||% "Proportion") +
    ggplot2::theme(
      legend.position = legend.position,
      legend.direction = legend.direction
    ) +
    theme_obj

  if (isTRUE(flip.y) && isFALSE(coords$uses_image)) {
    p <- p + ggplot2::scale_y_reverse()
  }
  if (!is.null(split.by)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", split.by)))
  }
  if (isTRUE(crop)) {
    radius_max <- max(plot_dat[[".radius"]], na.rm = TRUE)
    xpad <- max(diff(range(plot_dat$x, na.rm = TRUE)) * 0.04, radius_max)
    ypad <- max(diff(range(plot_dat$y, na.rm = TRUE)) * 0.04, radius_max)
    p <- p +
      ggplot2::coord_equal(
        xlim = range(plot_dat$x, na.rm = TRUE) + c(-xpad, xpad),
        ylim = range(plot_dat$y, na.rm = TRUE) + c(-ypad, ypad)
      )
  } else {
    p <- p + ggplot2::coord_equal()
  }
  p
}

spatial_dim_pie_values <- function(srt, group.by = NULL, values = NULL, cells) {
  if (!is.null(values)) {
    if (is.atomic(values) && is.null(dim(values))) {
      log_message(
        "{.arg values} must be a numeric matrix/data.frame for {.arg plot_type = 'pie'}",
        message_type = "error"
      )
    }
    mat <- as.data.frame(values, check.names = FALSE)
    if (is.null(rownames(mat)) || !any(cells %in% rownames(mat))) {
      log_message(
        "{.arg values} row names must match spatial spot names",
        message_type = "error"
      )
    }
    mat <- mat[cells, , drop = FALSE]
  } else {
    if (is.null(group.by) || length(group.by) == 0L) {
      log_message(
        "{.arg group.by} or {.arg values} must be provided for {.arg plot_type = 'pie'}",
        message_type = "error"
      )
    }
    missing_cols <- setdiff(group.by, colnames(srt@meta.data))
    if (length(missing_cols) > 0L) {
      log_message(
        "{.arg group.by} columns are not in srt meta.data: {.val {missing_cols}}",
        message_type = "error"
      )
    }
    mat <- srt@meta.data[cells, group.by, drop = FALSE]
  }
  if (ncol(mat) == 0L) {
    log_message(
      "Pie plotting requires at least one numeric column",
      message_type = "error"
    )
  }
  is_numeric <- vapply(mat, is.numeric, logical(1))
  if (!all(is_numeric)) {
    log_message(
      "All pie columns must be numeric",
      message_type = "error"
    )
  }
  as.matrix(mat)
}

spatial_dim_pie_radius <- function(coords, radius = NULL, scale = 0.45) {
  if (
    length(scale) != 1L ||
      !is.numeric(scale) ||
      is.na(scale) ||
      scale <= 0
  ) {
    log_message(
      "{.arg pie.radius.scale} must be a single positive number",
      message_type = "error"
    )
  }
  if (!is.null(radius)) {
    if (!is.numeric(radius)) {
      log_message(
        "{.arg pie.radius} must be numeric",
        message_type = "error"
      )
    }
    if (length(radius) == 1L) {
      radius <- rep(radius, nrow(coords))
      spatial_dim_validate_pie_radius(radius)
      return(radius)
    }
    if (!is.null(names(radius))) {
      radius <- as.numeric(radius[rownames(coords)])
      spatial_dim_validate_pie_radius(radius)
      return(radius)
    }
    if (length(radius) != nrow(coords)) {
      log_message(
        "{.arg pie.radius} must have length 1 or one value per spot",
        message_type = "error"
      )
    }
    spatial_dim_validate_pie_radius(radius)
    return(radius)
  }
  x <- sort(unique(coords$x[is.finite(coords$x)]))
  y <- sort(unique(coords$y[is.finite(coords$y)]))
  dx <- diff(x)
  dy <- diff(y)
  steps <- c(dx[dx > 0], dy[dy > 0])
  step <- suppressWarnings(stats::median(steps, na.rm = TRUE))
  if (!is.finite(step) || step <= 0) {
    xrange <- diff(range(coords$x, na.rm = TRUE))
    yrange <- diff(range(coords$y, na.rm = TRUE))
    step <- max(xrange, yrange, 1, na.rm = TRUE) / sqrt(max(nrow(coords), 1))
  }
  rep(step * scale, nrow(coords))
}

spatial_dim_validate_pie_radius <- function(radius) {
  if (any(!is.finite(radius) | radius <= 0)) {
    log_message(
      "{.arg pie.radius} values must be finite positive numbers",
      message_type = "error"
    )
  }
}

spatial_dim_coords <- function(
  srt,
  image = NULL,
  coord.cols = c("col", "row"),
  overlay_image = TRUE
) {
  images <- tryCatch(SeuratObject::Images(srt), error = function(e) character())
  image_info <- NULL
  if (length(images) > 0L) {
    image <- image %||% images[1L]
    if (!image %in% images) {
      log_message(
        "{.arg image} {.val {image}} is not present in {.cls Seurat}",
        message_type = "error"
      )
    }
    coords <- as.data.frame(SeuratObject::GetTissueCoordinates(srt[[image]]))
    cell_col <- if ("cell" %in% colnames(coords)) "cell" else NULL
    cells <- if (is.null(cell_col)) rownames(coords) else coords[[cell_col]]
    rownames(coords) <- cells
    x_col <- spatial_dim_pick_col(
      coords,
      c("x", "pxl_col_in_fullres", "imagecol")
    )
    y_col <- spatial_dim_pick_col(
      coords,
      c("y", "pxl_row_in_fullres", "imagerow")
    )
    scale <- spatial_dim_image_scale(srt[[image]])
    image_array <- srt[[image]]@image
    image_height <- dim(image_array)[1L]
    image_width <- dim(image_array)[2L]
    out <- data.frame(
      x = coords[[x_col]] * scale,
      y = image_height - coords[[y_col]] * scale,
      row.names = cells,
      stringsAsFactors = FALSE
    )
    image_info <- list(
      image = image_array,
      width = image_width,
      height = image_height
    )
    return(list(data = out, image = image_info, uses_image = TRUE))
  }

  if (!all(coord.cols %in% colnames(srt@meta.data))) {
    log_message(
      "Spatial coordinates were not found. Provide a Seurat image or metadata columns {.val {coord.cols}}.",
      message_type = "error"
    )
  }
  out <- data.frame(
    x = srt@meta.data[[coord.cols[1L]]],
    y = srt@meta.data[[coord.cols[2L]]],
    row.names = rownames(srt@meta.data),
    stringsAsFactors = FALSE
  )
  list(data = out, image = image_info, uses_image = FALSE)
}

spatial_dim_value_items <- function(values) {
  if (is.null(values)) {
    return(character())
  }
  if (is.atomic(values) && is.null(dim(values))) {
    return("value")
  }
  values <- as.data.frame(values, check.names = FALSE)
  if (ncol(values) == 0L) {
    log_message(
      "{.arg values} must contain at least one column",
      message_type = "error"
    )
  }
  colnames(values)
}

spatial_dim_values_from_input <- function(
  values,
  item,
  cells,
  show_na = FALSE
) {
  if (is.atomic(values) && is.null(dim(values))) {
    if (is.null(names(values))) {
      log_message(
        "{.arg values} vector must be named with spatial spot names",
        message_type = "error"
      )
    }
    out <- values[cells]
    names(out) <- cells
  } else {
    values <- as.data.frame(values, check.names = FALSE)
    if (is.null(rownames(values)) || !any(cells %in% rownames(values))) {
      log_message(
        "{.arg values} matrix/data.frame row names must match spatial spot names",
        message_type = "error"
      )
    }
    if (!item %in% colnames(values)) {
      log_message(
        "{.arg values} does not contain column {.val {item}}",
        message_type = "error"
      )
    }
    out <- values[cells, item, drop = TRUE]
  }
  if (!is.numeric(out)) {
    out <- as.character(out)
    if (isTRUE(show_na)) {
      out[is.na(out)] <- "NA"
    }
    out <- factor(out, levels = unique(out))
  }
  out
}

spatial_dim_values <- function(
  srt,
  item,
  type = c("metadata", "feature", "values"),
  assay = NULL,
  layer = "data",
  values = NULL,
  cells,
  show_na = FALSE
) {
  type <- match.arg(type)
  if (type == "metadata") {
    if (!item %in% colnames(srt@meta.data)) {
      log_message(
        "{.arg group.by} {.val {item}} is not in srt meta.data",
        message_type = "error"
      )
    }
    values <- srt@meta.data[cells, item, drop = TRUE]
    if (!is.numeric(values)) {
      values <- as.character(values)
      if (isTRUE(show_na)) {
        values[is.na(values)] <- "NA"
      }
      values <- factor(values, levels = unique(values))
    }
    return(values)
  }

  if (type == "values") {
    return(spatial_dim_values_from_input(values, item, cells, show_na))
  }

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  expr <- GetAssayData5(srt, assay = assay, layer = layer)
  if (!item %in% rownames(expr)) {
    log_message(
      "{.arg features} {.val {item}} is not in assay {.val {assay}}",
      message_type = "error"
    )
  }
  as.numeric(expr[item, cells, drop = TRUE])
}

spatial_dim_single_plot <- function(
  plot_dat,
  value_col,
  value_name,
  split.by,
  image_info = NULL,
  overlay_image = TRUE,
  image.alpha = 1,
  crop = TRUE,
  flip.y = FALSE,
  pt.size = 1,
  pt.alpha = 0.9,
  stroke = 0.1,
  palette = "Chinese",
  palcolor = NULL,
  bg_color = "grey20",
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = value_name,
  theme_use = "theme_blank",
  theme_args = list()
) {
  theme_obj <- spatial_dim_drop_coord(do.call(theme_use, theme_args))
  p <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = .data$x, y = .data$y))
  if (isTRUE(overlay_image) && !is.null(image_info)) {
    p <- p +
      ggplot2::annotation_raster(
        spatial_dim_raster(image_info$image, image.alpha),
        xmin = 0,
        xmax = image_info$width,
        ymin = 0,
        ymax = image_info$height
      )
  }

  values <- plot_dat[[value_col]]
  if (is.numeric(values)) {
    cols <- palette_colors(
      type = "continuous",
      palette = palette,
      palcolor = palcolor
    )
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(color = .data[[value_col]]),
        size = pt.size,
        alpha = pt.alpha
      ) +
      ggplot2::scale_color_gradientn(colors = cols, na.value = "grey80")
  } else {
    lvls <- levels(factor(values))
    cols <- palette_colors(lvls, palette = palette, palcolor = palcolor)
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(fill = .data[[value_col]]),
        shape = 21,
        color = bg_color,
        stroke = stroke,
        size = pt.size,
        alpha = pt.alpha
      ) +
      ggplot2::scale_fill_manual(values = cols, na.value = "grey80")
  }

  legend_labs <- if (is.numeric(values)) {
    ggplot2::labs(x = NULL, y = NULL, color = legend.title)
  } else {
    ggplot2::labs(x = NULL, y = NULL, fill = legend.title)
  }
  p <- p +
    legend_labs +
    ggplot2::theme(
      legend.position = legend.position,
      legend.direction = legend.direction
    ) +
    theme_obj

  if (isTRUE(flip.y)) {
    p <- p + ggplot2::scale_y_reverse()
  }
  if (isTRUE(crop)) {
    xpad <- diff(range(plot_dat$x, na.rm = TRUE)) * 0.04
    ypad <- diff(range(plot_dat$y, na.rm = TRUE)) * 0.04
    p <- p +
      ggplot2::coord_equal(
        xlim = range(plot_dat$x, na.rm = TRUE) + c(-xpad, xpad),
        ylim = range(plot_dat$y, na.rm = TRUE) + c(-ypad, ypad)
      )
  } else {
    p <- p + ggplot2::coord_equal()
  }
  if (!identical(split.by, ".split")) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", split.by)))
  }
  p
}

spatial_dim_image_scale <- function(image) {
  sf <- image@scale.factors
  if (!is.null(sf$lowres)) {
    return(sf$lowres)
  }
  if (!is.null(sf$hires)) {
    return(sf$hires)
  }
  1
}

spatial_dim_raster <- function(image, alpha = 1) {
  alpha <- max(min(alpha, 1), 0)
  if (alpha < 1) {
    image <- image * alpha + (1 - alpha)
  }
  as.raster(image)
}

spatial_dim_pick_col <- function(x, candidates) {
  nm <- colnames(x)
  hit <- candidates[tolower(candidates) %in% tolower(nm)][1L]
  if (is.na(hit)) {
    log_message(
      "Unable to resolve spatial coordinate columns from {.val {colnames(x)}}",
      message_type = "error"
    )
  }
  nm[match(tolower(hit), tolower(nm))]
}

spatial_dim_drop_coord <- function(x) {
  if (!is.list(x) || inherits(x, "theme")) {
    return(x)
  }
  if (
    all(vapply(
      x,
      function(item) is.list(item) && length(item) == 1L,
      logical(1)
    ))
  ) {
    x <- unlist(x, recursive = FALSE)
  }
  Filter(function(item) !inherits(item, "Coord"), x)
}
