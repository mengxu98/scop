#' @title Plot spatial cell boundaries
#'
#' @description
#' Plot real cell segmentation polygons supplied directly, stored in a result
#' object, or extracted from a Seurat spatial image. Spot centers are never
#' converted into synthetic polygons.
#'
#' @param object Optional `Seurat` object used to extract boundaries or values.
#' @param res Optional result list containing a `boundaries` data frame.
#' @param boundaries Optional boundary data frame.
#' @param cells Optional cell or spot identifiers to retain.
#' @param image Seurat image name. Multi-image objects require an explicit name.
#' @param crop Whether to crop the plot to the selected boundaries.
#' @param group.by Boundary column or Seurat metadata column used for filling.
#' @param features Features to display. Multiple features return a patchwork.
#' @param palette,palcolor Palette name or explicit colors.
#' @param fill.alpha Polygon fill opacity.
#' @param boundary.color,boundary.linewidth Boundary appearance.
#' @param theme_use,theme_args scop theme and its arguments.
#' @param ... Additional arguments passed to `ggplot2::geom_polygon()`.
#'
#' @return A `ggplot` or patchwork object.
#'
#' @examples
#' boundaries <- data.frame(
#'   cell_id = rep(c("cell1", "cell2"), each = 4),
#'   polygon_id = rep(c("p1", "p2"), each = 4),
#'   ring_id = 1,
#'   vertex_order = rep(1:4, 2),
#'   x = c(0, 1, 1, 0, 1.2, 2.2, 2.2, 1.2),
#'   y = c(0, 0, 1, 1, 0, 0, 1, 1),
#'   cell_type = rep(c("A", "B"), each = 4)
#' )
#' SpatialCellPlot(boundaries = boundaries, group.by = "cell_type")
#'
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
  theme_use = "theme_scop",
  theme_args = list(),
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
    images <- tryCatch(SeuratObject::Images(object), error = function(e) character())
    if (length(images) == 0L) {
      log_message("No Seurat spatial image is available for boundary extraction", message_type = "error")
    }
    if (is.null(image)) {
      if (length(images) > 1L) {
        log_message(
          "Multiple spatial images are available; select one with {.arg image}: {.val {images}}",
          message_type = "error"
        )
      }
      image <- images[[1L]]
    } else if (!image %in% images) {
      log_message("{.arg image} {.val {image}} is not present in {.cls Seurat}", message_type = "error")
    }
    boundary_names <- tryCatch(SeuratObject::Boundaries(object[[image]]), error = function(e) character())
    if (!"segmentation" %in% boundary_names) {
      log_message("The selected image does not contain segmentation boundaries", message_type = "error")
    }
    boundary_name <- "segmentation"
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
      value_tables[[group.by]] <- object@meta.data[boundaries$cell_id, group.by]
    } else {
      log_message("{.arg group.by} {.val {group.by}} was not found", message_type = "error")
    }
  } else if (length(features) > 0L) {
    for (feature in unique(features)) {
      if (feature %in% colnames(boundaries)) {
        value_tables[[feature]] <- boundaries[[feature]]
      } else {
        if (is.null(object)) {
          log_message("A Seurat {.arg object} is required to fetch feature {.val {feature}}", message_type = "error")
        }
        fetched <- tryCatch(
          SeuratObject::FetchData(object, vars = feature, layer = "data"),
          error = function(e) SeuratObject::FetchData(object, vars = feature)
        )
        value_tables[[feature]] <- fetched[boundaries$cell_id, feature]
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
    layer <- do.call(
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
    p <- ggplot2::ggplot() + layer
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
      spatial_dim_drop_coord(do.call(theme_use, theme_args))
    if (isTRUE(crop)) {
      xr <- range(dat$x, na.rm = TRUE)
      yr <- range(dat$y, na.rm = TRUE)
      xp <- max(diff(xr) * 0.04, .Machine$double.eps)
      yp <- max(diff(yr) * 0.04, .Machine$double.eps)
      p <- p + ggplot2::coord_equal(
        xlim = xr + c(-xp, xp),
        ylim = yr + c(-yp, yp)
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
