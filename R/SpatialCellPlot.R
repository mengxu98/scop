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
#' @param crop Whether to crop the plot to the selected boundaries.
#' @param group.by Boundary column or Seurat metadata column used for filling.
#' @param features Features to display. Multiple features return a patchwork.
#' @param assay Assay used when `features` are fetched from `object`.
#' @param layer Layer used when `features` are fetched from `object`.
#' @param palette,palcolor Palette name or explicit colors.
#' @param fill.alpha Polygon fill opacity.
#' @param boundary.color,boundary.linewidth Boundary appearance.
#' @param theme_use,theme_args Theme used to style the plot. Default is
#' `"theme_spatial"`.
#' @param ... Additional arguments passed to `ggplot2::geom_polygon()`.
#'
#' @details
#' Boundary tables require `cell_id`, `x`, and `y` columns in polygon vertex
#' order. Seurat images must contain segmentation boundaries.
#'
#' @return A `ggplot` or patchwork object.
#'
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' SpatialCellPlot(visium_human_pancreas_sub, group.by = "CellType")
#' @export
SpatialCellPlot <- function(
  object = NULL,
  res = NULL,
  boundaries = NULL,
  cells = NULL,
  image = NULL,
  crop = TRUE,
  group.by = NULL,
  features = NULL,
  palette = "Paired",
  palcolor = NULL,
  fill.alpha = 0.7,
  boundary.color = "grey30",
  boundary.linewidth = 0.1,
  theme_use = "theme_spatial",
  theme_args = list(),
  assay = NULL,
  layer = "data",
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
  if (is.null(boundaries) && !is.null(object)) {
    image <- spatial_image_resolve(
      srt = object,
      image = image,
      image_policy = "strict"
    )$image
    if (is.null(image)) {
      log_message("No Seurat spatial image is available for boundary extraction", message_type = "error")
    }
    boundary_name <- spatial_segmentation_name(object[[image]], required = TRUE)
    boundaries <- tryCatch(
      as.data.frame(SeuratObject::GetTissueCoordinates(object[[image]][[boundary_name]])),
      error = function(e) {
        tryCatch(
          as.data.frame(SeuratObject::GetTissueCoordinates(object[[image]], which = boundary_name)),
          error = function(e2) NULL
        )
      }
    )
  }
  if (is.null(boundaries) || !is.data.frame(boundaries) || nrow(boundaries) == 0L) {
    log_message("Provide real segmentation data through {.arg boundaries}, {.arg res}, or a Seurat image", message_type = "error")
  }

  assay <- assay %||% if (!is.null(object)) SeuratObject::DefaultAssay(object) else NULL
  if (!is.null(features) && !is.null(object) && is.null(assay)) {
    log_message("An {.arg assay} is required when fetching features from {.arg object}", message_type = "error")
  }

  boundaries <- spatial_boundary_validate(boundaries, image = image)
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
      value_tables[[group.by]] <- object@meta.data[[group.by]][cell_idx]
    } else {
      log_message("{.arg group.by} {.val {group.by}} was not found", message_type = "error")
    }
  } else if (length(features) > 0L) {
    needed_features <- setdiff(unique(features), colnames(boundaries))
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
      } else {
        value_tables[[feature]] <- as.numeric(expr[feature, boundaries$cell_id, drop = TRUE])
      }
    }
  } else {
    value_tables[["Spatial cells"]] <- rep("cell", nrow(boundaries))
  }

  polygon_args <- list(...)
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
    p <- ggplot2::ggplot() + polygon_layer
    if (is.numeric(dat$.value)) {
      p <- p + ggplot2::scale_fill_gradientn(
        colors = palette_colors(type = "continuous", palette = palette, palcolor = palcolor),
        na.value = "grey80"
      )
    } else {
      lvls <- levels(factor(dat$.value))
      p <- p + ggplot2::scale_fill_manual(
        values = palette_colors(lvls, palette = palette, palcolor = palcolor),
        na.value = "grey80"
      )
    }
    p <- p +
      ggplot2::labs(x = NULL, y = NULL, fill = value_name) +
      apply_plot_theme(theme_use = theme_use, theme_args = theme_args)
    if (isTRUE(crop)) {
      limits <- spatial_crop_limits(dat$x, dat$y)
      p <- p + ggplot2::coord_equal(
        xlim = limits$xlim,
        ylim = limits$ylim
      )
    } else {
      p <- p + ggplot2::coord_equal()
    }
    p
  })
  if (length(plots) == 1L) {
    plots[[1L]]
  } else {
    check_r("patchwork", verbose = FALSE)
    wrap_plots <- get_namespace_fun("patchwork", "wrap_plots")
    wrap_plots(plots)
  }
}
