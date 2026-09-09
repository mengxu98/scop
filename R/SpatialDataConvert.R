#' @title Convert Seurat to SpatialExperiment
#'
#' @description
#' Create a lightweight `SpatialExperiment` from a spatial Seurat object using
#' one assay layer, metadata, and resolved spatial coordinates.
#'
#' @md
#' @param srt A `Seurat` object.
#' @param assay Assay to export. If `NULL`, the default assay is used.
#' @param layer Assay layer to export.
#' @param coord.cols Metadata coordinate columns. By default, SCOP resolves
#' `x/y` first and then `col/row`.
#' @param image Optional Seurat image name. When present, image-derived
#' coordinates are used.
#' @param coordinate_space Coordinate space exported to `spatialCoords`.
#'   `"legacy_display"` preserves the historical scaled/y-flipped behavior;
#'   `"raw"` (the default) preserves analysis distances. The source and
#'   transform are stored in `metadata(spe)$scop_spatial_coordinates` so that
#'   an explicit display export can be inverted by [spe_to_srt()].
#' @param include_meta Whether to include Seurat metadata as `colData`.
#'
#' @return A `SpatialExperiment`.
#' @export
srt_to_spe <- function(
  srt,
  assay = NULL,
  layer = "counts",
  coord.cols = c("col", "row"),
  image = NULL,
  include_meta = TRUE,
  coordinate_space = c("raw", "legacy_display")
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  check_r(c("SpatialExperiment", "SummarizedExperiment", "S4Vectors"), verbose = FALSE)
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!assay %in% SeuratObject::Assays(srt)) {
    log_message("{.arg assay} {.val {assay}} is not present in {.cls Seurat}", message_type = "error")
  }
  expr <- GetAssayData5(srt, assay = assay, layer = layer)
  coordinate_space <- match.arg(coordinate_space)
  resolved <- spatial_analysis_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    coordinate_space = coordinate_space
  )
  coords <- resolved$data
  cells <- intersect(colnames(expr), rownames(coords))
  if (length(cells) == 0L) {
    log_message("No assay cells match spatial coordinates", message_type = "error")
  }
  expr <- expr[, cells, drop = FALSE]
  coords <- coords[cells, c("x", "y"), drop = FALSE]
  coldata <- if (isTRUE(include_meta)) {
    S4Vectors::DataFrame(srt@meta.data[cells, , drop = FALSE])
  } else {
    S4Vectors::DataFrame(row.names = cells)
  }
  SpatialExperiment::SpatialExperiment(
    assays = list(scop_input = expr),
    colData = coldata,
    spatialCoords = as.matrix(coords),
    metadata = list(scop_spatial_coordinates = list(
      source = resolved$source,
      transform = resolved$transform
    ))
  )
}

#' @title Convert SpatialExperiment to Seurat
#'
#' @description
#' Create a Seurat object from a `SpatialExperiment`, preserving `colData` and
#' spatial coordinates as metadata columns.
#' SCOP display exports are inverted to raw coordinates using their saved
#' transform. External SpatialExperiment coordinates without SCOP provenance
#' are treated as raw coordinates in their supplied units. Original import
#' provenance is retained in `object@misc$spatial_coordinate_import`.
#'
#' @md
#' @param spe A `SpatialExperiment` or `SummarizedExperiment`.
#' @param assay Assay name for the created Seurat assay.
#' @param layer Assay from `spe` to use as counts. If `NULL`, the first assay is
#' used.
#' @param coord.cols Metadata names used for spatial coordinates in Seurat.
#' @param project Project name passed to `Seurat::CreateSeuratObject()`.
#'
#' @return A `Seurat` object.
#' @export
spe_to_srt <- function(
  spe,
  assay = "Spatial",
  layer = NULL,
  coord.cols = c("x", "y"),
  project = "SpatialExperiment"
) {
  check_r(c("SpatialExperiment", "SummarizedExperiment"), verbose = FALSE)
  if (!inherits(spe, "SummarizedExperiment")) {
    log_message("{.arg spe} must be a {.cls SummarizedExperiment} object", message_type = "error")
  }
  assay_names <- SummarizedExperiment::assayNames(spe)
  if (length(assay_names) == 0L) {
    log_message("{.arg spe} must contain at least one assay", message_type = "error")
  }
  layer <- layer %||% assay_names[1L]
  if (!layer %in% assay_names) {
    log_message("{.arg layer} {.val {layer}} is not an assay in {.arg spe}", message_type = "error")
  }
  counts <- SummarizedExperiment::assay(spe, layer)
  if (length(coord.cols) != 2L || anyNA(coord.cols) ||
    any(!nzchar(coord.cols)) || anyDuplicated(coord.cols)) {
    log_message("{.arg coord.cols} must contain two unique non-empty names", message_type = "error")
  }
  meta <- as.data.frame(SummarizedExperiment::colData(spe), optional = TRUE)
  if (nrow(meta) == 0L) {
    meta <- data.frame(row.names = colnames(counts))
  }
  coords <- tryCatch(SpatialExperiment::spatialCoords(spe), error = function(e) NULL)
  provenance <- S4Vectors::metadata(spe)$scop_spatial_coordinates
  if (!is.null(coords) && ncol(coords) >= 2L) {
    if (nrow(coords) != ncol(counts)) {
      log_message("SpatialExperiment coordinates must match all assay columns", message_type = "error")
    }
    if (!is.null(rownames(coords))) {
      if (anyDuplicated(rownames(coords)) || !setequal(rownames(coords), colnames(counts))) {
        log_message("SpatialExperiment coordinate IDs do not match assay columns", message_type = "error")
      }
      coords <- coords[colnames(counts), , drop = FALSE]
    }
    x <- spatial_coordinate_numeric(coords[, 1L])
    y <- spatial_coordinate_numeric(coords[, 2L])
    if (!is.null(provenance)) {
      space <- provenance$source$coordinate_space
      if (!is.character(space) || length(space) != 1L || is.na(space) ||
        !space %in% c("raw", "display", "legacy_display")) {
        log_message("Invalid SCOP SpatialExperiment coordinate provenance", message_type = "error")
      }
      if (space != "raw") {
        transform <- provenance$transform
        scale <- transform$scale
        if (length(scale) != 1L || !is.finite(scale) || scale <= 0) {
          log_message("Display import requires a positive finite saved scale", message_type = "error")
        }
        if (isTRUE(transform$y_flip)) {
          if (length(transform$image_height) != 1L || !is.finite(transform$image_height)) {
            log_message("Display import requires a finite saved image height", message_type = "error")
          }
          y <- transform$image_height - y
        }
        x <- x / scale
        y <- y / scale
      }
    }
    if (any(!is.finite(x) | !is.finite(y))) {
      log_message("SpatialExperiment coordinates must be finite", message_type = "error")
    }
    meta <- meta[colnames(counts), , drop = FALSE]
    meta[[coord.cols[1L]]] <- x
    meta[[coord.cols[2L]]] <- y
  }
  out <- Seurat::CreateSeuratObject(
    counts = counts,
    assay = assay,
    meta.data = meta,
    project = project
  )
  if (!is.null(provenance)) out@misc$spatial_coordinate_import <- provenance
  out
}
