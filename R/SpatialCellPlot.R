#' @title Plot spatial cell boundaries
#'
#' @description
#' Plot cell segmentation polygons from a boundary table or Seurat spatial image.
#'
#' @param object Optional `Seurat` object used to extract boundaries or values.
#' @param res Optional result list containing a `boundaries` data frame.
#' @param boundaries Optional boundary data frame.
#' @param cells Optional cell or spot identifiers to retain.
#' @param image Seurat image name. Multi-image objects require an explicit name.
#' @param boundary Boundary name in the image. `NULL` selects segmentation when
#' available and otherwise plots centroids. Explicit names must exist.
#' @param xlim,ylim Optional increasing ROI limits in raw coordinate units.
#' Cells are selected by centroids; boundary vertices are clipped at the limits.
#' @param image.scale,overlay_image,image.alpha Raster scale, visibility, and opacity.
#' @param flip.y Reverse the y axis. `NULL` uses the official FOV orientation;
#' use `FALSE` for the historical FOV layout. Visium rasters include the flip.
#' @param pt.size Centroid point size.
#' @inheritParams SpatialSpotPlot
#' @param crop Whether to crop the plot to the selected boundaries.
#' @param group.by Boundary column or Seurat metadata column used for filling.
#' @param features Features to display. Multiple features return a patchwork.
#' @param assay Assay used when `features` are fetched from `object`.
#' @param layer Layer used when `features` are fetched from `object`.
#' @param palette,palcolor Palette name or explicit colors. The default uses
#' `"Paired"` for categories and sequential `"YlGnBu"` for numeric values.
#' Named colors are matched to category names.
#' @param fill.alpha Polygon fill opacity.
#' @param boundary.color,boundary.linewidth Boundary appearance.
#' @param theme_use,theme_args Theme used to style the plot. Default is
#' `"theme_spatial"`.
#' @param ... Arguments passed to other functions.
#' @param combine Return combined panels, or a named list when `FALSE`.
#' @param nrow,ncol,byrow Layout controls for multiple feature panels.
#' @param legend.position,legend.direction,legend.title Legend controls.
#'
#' @details
#' Boundary tables require `cell_id`, `x`, and `y` columns in polygon vertex
#' order, with optional `polygon_id` and `ring_id` columns for multipart cells
#' and holes. Object-derived boundaries use the same display transform as
#' [SpatialSpotPlot()]. ROI filtering does not subset the input object.
#'
#' @return A `ggplot`, patchwork, or named list of plots.
#'
#'
#' @export
#' @examples
#' # Constructed polygons demonstrate plotting, not a segmentation algorithm.
#' boundaries <- data.frame(
#'   cell_id = rep(c("cell1", "cell2"), each = 4),
#'   x = c(0, 2, 2, 0, 3, 5, 5, 3), y = rep(c(0, 0, 1, 1), 2),
#'   cell_type = rep(c("A", "B"), each = 4)
#' )
#' SpatialCellPlot(boundaries = boundaries, group.by = "cell_type")
SpatialCellPlot <- function(
  object = NULL,
  res = NULL,
  boundaries = NULL,
  cells = NULL,
  image = NULL,
  crop = TRUE,
  group.by = NULL,
  features = NULL,
  palette = NULL,
  palcolor = NULL,
  fill.alpha = 0.7,
  boundary.color = "grey30",
  boundary.linewidth = 0.1,
  theme_use = "theme_spatial",
  theme_args = list(),
  assay = NULL,
  layer = "data",
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = NULL,
  boundary = NULL,
  xlim = NULL,
  ylim = NULL,
  image.scale = c("lowres", "hires"),
  overlay_image = TRUE,
  image.alpha = 1,
  flip.y = NULL,
  pt.size = 1,
  molecules = NULL,
  mols.size = 0.1,
  mols.alpha = 1,
  mols.colors = NULL,
  ...
) {
  if (!is.null(object) && !inherits(object, "Seurat")) {
    log_message("{.arg object} must be a {.cls Seurat} object", message_type = "error")
  }
  if (!is.null(image) && (!is.character(image) || length(image) != 1L || is.na(image) || !nzchar(image))) {
    log_message("{.arg image} must be one non-empty image name", message_type = "error")
  }
  if (!is.null(group.by) && (!is.character(group.by) || length(group.by) != 1L || is.na(group.by) || !nzchar(group.by))) {
    log_message("{.arg group.by} must be one non-empty column name", message_type = "error")
  }
  if (is.null(boundaries) && !is.null(res)) {
    boundaries <- if (is.data.frame(res)) res else res$boundaries
  }
  display <- NULL
  roi <- NULL
  if (is.null(boundaries) && !is.null(object)) {
    image <- spatial_image_resolve(
      srt = object,
      image = image,
      image_policy = "strict"
    )$image
    if (is.null(image)) {
      log_message("No Seurat spatial image is available for boundary extraction", message_type = "error")
    }
    info <- spatial_input_info(object, assay = assay, image = image)
    assay <- info$assay
    raw <- SpatialCoordinates(object, image = image)
    cells <- spatial_roi_cells(raw$data, cells, xlim, ylim)
    boundary_name <- boundary %||% spatial_segmentation_name(object[[image]])
    if (!is.null(boundary_name) && !boundary_name %in% SeuratObject::Boundaries(object[[image]])) {
      log_message("Requested boundary is not present in the selected image", message_type = "error")
    }
    if (is.null(boundary_name) || inherits(object[[image]][[boundary_name]], "Centroids")) {
      return(SpatialSpotPlot(object,
        cells = cells, image = image, assay = assay, layer = layer,
        group.by = group.by, features = features,
        values = if (is.null(group.by) && !length(features)) stats::setNames(rep("cell", length(cells)), cells) else NULL,
        palette = palette %||% if (length(features)) "YlGnBu" else "Paired", palcolor = palcolor,
        crop = crop, xlim = xlim, ylim = ylim, image.scale = image.scale,
        overlay_image = overlay_image, image.alpha = image.alpha, flip.y = flip.y,
        pt.size = pt.size, pt.alpha = fill.alpha, stroke = 0,
        show_na = TRUE, theme_use = theme_use, theme_args = theme_args,
        combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,
        legend.position = legend.position, legend.direction = legend.direction, legend.title = legend.title,
        molecules = molecules, mols.size = mols.size, mols.alpha = mols.alpha, mols.colors = mols.colors
      ))
    }
    boundaries <- spatial_segmentation_table(object[[image]][[boundary_name]], cells)
    if (inherits(object[[image]], "VisiumV2") && spatial_image_x_orientation(object, image) == "vertical") {
      boundaries[c("x", "y")] <- boundaries[c("y", "x")]
    }
    display <- spatial_dim_coords(object, image = image, image.scale = image.scale, overlay_image = overlay_image)
    if (is.null(flip.y)) flip.y <- inherits(object[[image]], "FOV") && !inherits(object[[image]], "VisiumV2")
    boundaries <- spatial_coords_to_display(boundaries, display$transform)
    roi <- spatial_coords_to_display(data.frame(
      x = xlim %||% range(raw$data$x),
      y = ylim %||% range(raw$data$y)
    ), display$transform)
  }
  if (is.null(boundaries) || !is.data.frame(boundaries) || nrow(boundaries) == 0L) {
    log_message("Provide real segmentation data through {.arg boundaries}, {.arg res}, or a Seurat image", message_type = "error")
  }

  assay <- assay %||% if (!is.null(object)) SeuratObject::DefaultAssay(object) else NULL
  if (!is.null(features) && !is.null(object) && is.null(assay)) {
    log_message("An {.arg assay} is required when fetching features from {.arg object}", message_type = "error")
  }

  boundaries <- spatial_boundary_validate(boundaries, image = image)
  if (is.null(display) && (!is.null(xlim) || !is.null(ylim))) {
    log_message("ROI selection requires an object with cell centroids; use cells for boundary tables", message_type = "error")
  }
  if (!is.null(cells)) {
    boundaries <- boundaries[boundaries$cell_id %in% cells, , drop = FALSE]
  }
  if (nrow(boundaries) == 0L) {
    log_message("No segmentation boundaries remain after filtering", message_type = "error")
  }
  if (!is.null(group.by) && length(features) > 0L) {
    log_message("Use either {.arg group.by} or {.arg features}, not both", message_type = "error")
  }

  value_tables <- list()
  if (!is.null(group.by)) {
    if (group.by %in% colnames(boundaries)) {
      value_tables[[group.by]] <- boundaries[[group.by]]
    } else if (!is.null(object) && group.by %in% colnames(object@meta.data)) {
      cell_idx <- match(boundaries$cell_id, rownames(object@meta.data))
      if (anyNA(cell_idx)) {
        log_message("Boundary cell IDs do not match the selected Seurat object", message_type = "error")
      }
      value_tables[[group.by]] <- object@meta.data[[group.by]][cell_idx]
    } else {
      log_message("{.arg group.by} {.val {group.by}} was not found", message_type = "error")
    }
  } else if (length(features) > 0L) {
    metadata_features <- if (!is.null(object)) intersect(features, colnames(object[[]])) else character()
    needed_features <- setdiff(unique(features), c(colnames(boundaries), metadata_features))
    expr <- NULL
    if (length(needed_features) > 0L) {
      if (is.null(object)) {
        log_message("A Seurat {.arg object} is required to fetch features", message_type = "error")
      }
      expr <- tryCatch(
        GetAssayData5(
          object,
          assay = assay,
          layer = layer,
          features = needed_features,
          cells = unique(boundaries$cell_id)
        ),
        error = function(e) NULL
      )
      missing_features <- if (is.null(expr)) {
        needed_features
      } else {
        needed_features[!needed_features %in% rownames(expr)]
      }
      if (length(missing_features) > 0L) {
        log_message(
          "Feature {.val {missing_features[[1L]]}} is not available in assay {.val {assay}} layer {.val {layer}}",
          message_type = "error"
        )
      }
      if (!all(boundaries$cell_id %in% colnames(expr))) {
        log_message("Boundary cell IDs do not match the selected Seurat object", message_type = "error")
      }
    }
    for (feature in unique(features)) {
      if (feature %in% colnames(boundaries)) {
        value_tables[[feature]] <- boundaries[[feature]]
      } else if (feature %in% metadata_features) {
        value_tables[[feature]] <- object[[]][boundaries$cell_id, feature]
      } else {
        value_tables[[feature]] <- as.numeric(expr[feature, boundaries$cell_id, drop = TRUE])
      }
    }
  } else {
    value_tables[["Spatial cells"]] <- rep("cell", nrow(boundaries))
  }

  polygon_args <- merge_call_args(list(show.legend = TRUE), list(...))
  plots <- lapply(names(value_tables), function(value_name) {
    dat <- boundaries
    dat$.value <- value_tables[[value_name]]
    mapping <- ggplot2::aes(
      x = .data$x,
      y = .data$y,
      group = .data$.polygon_group,
      subgroup = .data$ring_id,
      fill = .data$.value
    )
    polygon_layer <- do.call(
      ggplot2::geom_polygon,
      c(
        list(
          mapping = mapping,
          data = dat,
          color = boundary.color,
          linewidth = boundary.linewidth,
          alpha = fill.alpha,
          rule = "evenodd"
        ),
        polygon_args
      )
    )
    p <- ggplot2::ggplot()
    if (isTRUE(overlay_image) && !is.null(display$image)) {
      p <- p + ggplot2::annotation_raster(spatial_dim_raster(display$image$image, image.alpha),
        xmin = display$image$xmin %||% 0, xmax = display$image$width,
        ymin = display$image$ymin %||% 0, ymax = display$image$height
      )
    }
    p <- p + polygon_layer
    if (is.numeric(dat$.value)) {
      p <- p + spatial_dim_continuous_scale(
        dat$.value,
        aesthetic = "fill",
        colors = spatial_palette_colors(type = "continuous", palette = palette %||% "YlGnBu", palcolor = palcolor)
      ) + ggplot2::guides(fill = spatial_colorbar_guide(legend.direction))
    } else {
      lvls <- levels(spatial_plot_factor(dat$.value))
      p <- p + ggplot2::scale_fill_manual(
        values = spatial_palette_colors(lvls, palette = palette %||% "Paired", palcolor = palcolor),
        na.value = "grey80", drop = FALSE
      )
    }
    p <- p +
      ggplot2::labs(
        x = NULL, y = NULL, fill = legend.title %||% value_name,
        title = if (length(value_tables) > 1L) value_name else NULL
      ) +
      apply_plot_theme(theme_use = theme_use, theme_args = theme_args) +
      ggplot2::theme(legend.position = legend.position, legend.direction = legend.direction)
    if (isTRUE(flip.y) && !isTRUE(display$uses_image)) p <- p + ggplot2::scale_y_reverse()
    if (isTRUE(crop) || !is.null(xlim) || !is.null(ylim)) {
      limits <- spatial_crop_limits(dat$x, dat$y)
      if (!is.null(xlim)) limits$xlim <- range(roi$x)
      if (!is.null(ylim)) limits$ylim <- range(roi$y)
      p <- p + ggplot2::coord_equal(
        xlim = limits$xlim,
        ylim = limits$ylim
      )
    } else {
      p <- p + ggplot2::coord_equal()
    }
    p
  })
  names(plots) <- names(value_tables)
  if (length(molecules)) {
    layers <- spatial_molecule_layers(
      object, image, molecules, xlim, ylim, image.scale,
      mols.size, mols.alpha, mols.colors
    )
    plots <- lapply(plots, function(p) p + layers)
  }
  combine_plot_list(plots, combine = combine, nrow = nrow, ncol = ncol, byrow = byrow)
}

spatial_plot_add_layers <- function(plots, layers) {
  if (inherits(plots, "patchwork")) {
    return(plots & layers)
  }
  if (inherits(plots, "ggplot")) {
    return(plots + layers)
  }
  lapply(plots, function(p) p + layers)
}

spatial_molecule_layers <- function(object, image, molecules, xlim, ylim, image.scale,
                                    size, alpha, colors) {
  image <- spatial_image_resolve(object, image = image)$image
  if (is.null(image)) log_message("Molecule overlays require a spatial image", message_type = "error")
  fov <- object[[image]]
  sets <- SeuratObject::Molecules(fov)
  available <- unique(unlist(lapply(sets, function(n) SeuratObject::Features(fov[[n]]))))
  if (length(setdiff(molecules, available))) {
    log_message("Requested molecules have not been loaded; use the official molecule loader first", message_type = "error")
  }
  dat <- do.call(rbind, lapply(sets, function(n) {
    use <- intersect(molecules, SeuratObject::Features(fov[[n]]))
    if (!length(use)) {
      return(NULL)
    }
    SeuratObject::GetTissueCoordinates(fov[[n]], features = use)
  }))
  if (inherits(fov, "VisiumV2") && spatial_image_x_orientation(object, image) == "vertical") {
    dat[c("x", "y")] <- dat[c("y", "x")]
  }
  if (!is.null(xlim)) dat <- dat[dat$x >= xlim[1L] & dat$x <= xlim[2L], , drop = FALSE]
  if (!is.null(ylim)) dat <- dat[dat$y >= ylim[1L] & dat$y <= ylim[2L], , drop = FALSE]
  transform <- SpatialCoordinates(object, image = image, image.scale = image.scale)$transform
  dat <- spatial_coords_to_display(dat, transform)
  list(
    ggnewscale::new_scale_color(),
    ggplot2::geom_point(
      data = dat, mapping = ggplot2::aes(x = .data$x, y = .data$y, color = .data$molecule),
      inherit.aes = FALSE, size = size, alpha = alpha
    ),
    ggplot2::scale_color_manual(
      values = spatial_palette_colors(molecules, palette = "Paired", palcolor = colors),
      name = "Molecule", drop = FALSE
    )
  )
}
