#' @title Spatial spot plot
#'
#' @md
#' @inheritParams CellDimPlot
#' @inheritParams scop-params
#' @param srt A `Seurat` object. The same object may be supplied as
#' `object =` for consistency with spatial plotting APIs.
#' @param object Optional alias for `srt`. Supply exactly one of `srt` or
#' `object`.
#' @param group.by Metadata columns to color spots by.
#' @param features Features to color spots by (from `assay`/`layer`).
#' @param values Spot-level values (named vector, matrix, or data.frame).
#' @param plot_type `"point"` or `"pie"`. Pie uses numeric `group.by` columns,
#' `values`, or `"<prefix>_prop_*"`/`"<prefix>_frac_*"` when `group.by` is
#' `"<prefix>_dominant_type"`.
#' @param plot.data,spot.by,color.by Long-format data for repeated spatial points
#' (e.g. cell-to-spot assignments).
#' @param geom `"point"` or `"jitter"` for `plot.data`.
#' @param overlay_image,image.alpha Draw the spatial image beneath spots.
#' @param image.scale Image scale factor matching the raster stored in the
#'   selected image. Use `"hires"` for a hires raster; do not modify Seurat
#'   scale-factor slots.
#' @param crop Crop the panel to plotted spots.
#' @param flip.y Reverse the y axis for metadata coordinates.
#' @param show_axes Keep axis text, ticks, and grid when `theme_use` is
#' `"theme_spatial"`.
#' @param theme_use Theme name or function.
#' @param pie.radius,pie.radius.scale Pie radius. `NULL` estimates from spot
#' spacing, then multiplies by `pie.radius.scale`.
#' @param stroke Point border width.
#' @param bg_color Point border color.
#' @param jitter_width,jitter_height Jitter when `geom = "jitter"`.
#'
#' @return A `ggplot`, `patchwork`, or list of `ggplot` objects.
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' SpatialSpotPlot(
#'   object = visium_human_pancreas_sub,
#'   group.by = "coda_label"
#' )
#'
#' SpatialSpotPlot(
#'   object = visium_human_pancreas_sub,
#'   features = rownames(visium_human_pancreas_sub)[1:2],
#'   layer = "counts"
#' )
SpatialSpotPlot <- function(
  srt = NULL,
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
  show_axes = FALSE,
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
  palette = "Spectral",
  palcolor = NULL,
  bg_color = "grey20",
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = NULL,
  theme_use = "theme_spatial",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  verbose = TRUE,
  object = NULL,
  image.scale = c("lowres", "hires")
) {
  srt <- spatial_resolve_srt(srt = srt, object = object)
  plot_type <- match.arg(plot_type)
  geom <- match.arg(geom)
  image.scale <- match.arg(image.scale)
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
      image.scale = image.scale,
      overlay_image = overlay_image,
      image.alpha = image.alpha,
      crop = crop,
      coord.cols = coord.cols,
      flip.y = flip.y,
      show_axes = show_axes,
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
      image.scale = image.scale,
      overlay_image = overlay_image,
      image.alpha = image.alpha,
      crop = crop,
      coord.cols = coord.cols,
      flip.y = flip.y,
      show_axes = show_axes,
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
    image.scale = image.scale,
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
      show_axes = show_axes,
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
  image.scale = c("lowres", "hires"),
  overlay_image = TRUE,
  image.alpha = 1,
  crop = TRUE,
  coord.cols = c("col", "row"),
  flip.y = TRUE,
  show_axes = FALSE,
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
  theme_use = "theme_spatial",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE
) {
  geom <- match.arg(geom)
  image.scale <- match.arg(image.scale)
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
    image.scale = image.scale,
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

  theme_args$show_axes <- show_axes
  theme_obj <- apply_plot_theme(
    theme_use = theme_use,
    theme_args = theme_args
  )
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
      spatial_dim_continuous_scale(values, aesthetic = "color", colors = cols) +
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
    limits <- spatial_crop_limits(df$x, df$y)
    p <- p +
      ggplot2::coord_equal(
        xlim = limits$xlim,
        ylim = limits$ylim
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
  image.scale = c("lowres", "hires"),
  overlay_image = TRUE,
  image.alpha = 1,
  crop = TRUE,
  coord.cols = c("col", "row"),
  flip.y = TRUE,
  show_axes = FALSE,
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
  theme_use = "theme_spatial",
  theme_args = list()
) {
  check_r("scatterpie", verbose = FALSE)
  image.scale <- match.arg(image.scale)
  coords <- spatial_dim_coords(
    srt = srt,
    image = image,
    image.scale = image.scale,
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
    return(
      spatial_empty_plot(
        "No spots with positive pie values",
        title = legend.title %||% "Proportion",
        theme_use = theme_use,
        theme_args = theme_args
      )
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
  theme_args$show_axes <- show_axes
  theme_obj <- apply_plot_theme(
    theme_use = theme_use,
    theme_args = theme_args
  )
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
    limits <- spatial_crop_limits(plot_dat$x, plot_dat$y, min_pad = radius_max)
    p <- p +
      ggplot2::coord_equal(
        xlim = limits$xlim,
        ylim = limits$ylim
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
    mat <- spatial_dim_pie_metadata(
      srt = srt,
      group.by = group.by,
      cells = cells
    )
  }
  if (ncol(mat) == 0L) {
    log_message(
      "Pie plotting requires at least one numeric column",
      message_type = "error"
    )
  }
  is_numeric <- vapply(mat, is.numeric, logical(1))
  if (!all(is_numeric)) {
    non_numeric <- colnames(mat)[!is_numeric]
    log_message(
      paste(
        "All pie columns must be numeric. Non-numeric column(s):",
        "{.val {non_numeric}}.",
        "For RCTD pies, use {.arg group.by} with {.val RCTD_prop_*}",
        "columns, or use {.val RCTD_dominant_type} when matching",
        "{.val RCTD_prop_*} columns are present."
      ),
      message_type = "error"
    )
  }
  as.matrix(mat)
}

spatial_dim_pie_metadata <- function(srt, group.by, cells) {
  mat <- srt@meta.data[cells, group.by, drop = FALSE]
  is_numeric <- vapply(mat, is.numeric, logical(1))
  if (all(is_numeric)) {
    return(mat)
  }
  inferred_cols <- spatial_dim_infer_pie_columns(srt, group.by)
  if (length(inferred_cols) > 0L) {
    return(srt@meta.data[cells, inferred_cols, drop = FALSE])
  }
  mat
}

spatial_dim_infer_pie_columns <- function(srt, group.by) {
  if (length(group.by) != 1L || !grepl("_dominant_type$", group.by)) {
    return(character())
  }
  prefix <- sub("_dominant_type$", "", group.by)
  meta_cols <- colnames(srt@meta.data)
  candidates <- meta_cols[
    startsWith(meta_cols, paste0(prefix, "_prop_")) |
      startsWith(meta_cols, paste0(prefix, "_frac_"))
  ]
  if (length(candidates) == 0L) {
    return(character())
  }
  is_numeric <- vapply(
    srt@meta.data[, candidates, drop = FALSE],
    is.numeric,
    logical(1)
  )
  candidates[is_numeric]
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
  image.scale = c("lowres", "hires"),
  overlay_image = TRUE,
  image_policy = "strict"
) {
  image.scale <- match.arg(image.scale)
  raw <- spatial_coords_raw(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    image.scale = image.scale,
    require_scale = TRUE,
    image_policy = image_policy
  )
  selected_image <- raw$source$image
  out <- spatial_coords_to_display(raw$data, raw$transform)
  out <- out[, c("x", "y"), drop = FALSE]
  attr(out, "spatial_source") <- utils::modifyList(
    raw$source,
    list(coordinate_space = "legacy_display")
  )
  attr(out, "spatial_transform") <- raw$transform
  image_info <- NULL
  if (!is.null(selected_image)) {
    image_array <- tryCatch(srt[[selected_image]]@image, error = function(e) NULL)
    if (!is.null(image_array)) {
      image_height <- dim(image_array)[1L]
      image_width <- dim(image_array)[2L]
      image_info <- list(
        image = image_array,
        width = image_width,
        height = image_height
      )
      return(list(
        data = out,
        image = image_info,
        uses_image = TRUE,
        source = attr(out, "spatial_source"),
        transform = raw$transform
      ))
    }
  }
  list(
    data = out,
    image = image_info,
    uses_image = FALSE,
    source = attr(out, "spatial_source"),
    transform = raw$transform
  )
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
  theme_use = "theme_spatial",
  theme_args = list(),
  show_axes = FALSE
) {
  theme_args$show_axes <- show_axes
  theme_obj <- apply_plot_theme(
    theme_use = theme_use,
    theme_args = theme_args
  )
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
  if (nrow(plot_dat) == 0L || all(is.na(values))) {
    return(spatial_empty_plot(
      "No values available for plotting",
      title = value_name,
      theme_use = theme_use,
      theme_args = theme_args
    ))
  }
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
      spatial_dim_continuous_scale(values, aesthetic = "color", colors = cols)
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
    limits <- spatial_crop_limits(plot_dat$x, plot_dat$y)
    p <- p +
      ggplot2::coord_equal(
        xlim = limits$xlim,
        ylim = limits$ylim
      )
  } else {
    p <- p + ggplot2::coord_equal()
  }
  if (!identical(split.by, ".split")) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", split.by)))
  }
  p
}

spatial_dim_continuous_scale <- function(values, aesthetic = c("color", "fill"), colors) {
  aesthetic <- match.arg(aesthetic)
  finite <- values[is.finite(values)]
  limits <- NULL
  breaks <- ggplot2::waiver()
  if (length(unique(finite)) == 1L) {
    center <- unique(finite)
    span <- max(abs(center) * 0.1, 0.5)
    limits <- center + c(-span, span)
    breaks <- center
  }
  if (identical(aesthetic, "color")) {
    return(ggplot2::scale_color_gradientn(
      colors = colors,
      na.value = "grey80",
      limits = limits,
      breaks = breaks
    ))
  }
  ggplot2::scale_fill_gradientn(
    colors = colors,
    na.value = "grey80",
    limits = limits,
    breaks = breaks
  )
}

spatial_dim_image_scale <- function(
  image,
  image.scale = c("lowres", "hires")
) {
  spatial_image_scale_info(image, image.scale = image.scale)$scale
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
