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
#'   `"raw"` preserves analysis distances.
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
  coordinate_space = c("legacy_display", "raw")
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
  coords <- spatial_analysis_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    coordinate_space = coordinate_space
  )$data
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
    spatialCoords = as.matrix(coords)
  )
}

#' @title Convert SpatialExperiment to Seurat
#'
#' @description
#' Create a Seurat object from a `SpatialExperiment`, preserving `colData` and
#' spatial coordinates as metadata columns.
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
  meta <- as.data.frame(SummarizedExperiment::colData(spe), optional = TRUE)
  if (nrow(meta) == 0L) {
    meta <- data.frame(row.names = colnames(counts))
  }
  coords <- tryCatch(SpatialExperiment::spatialCoords(spe), error = function(e) NULL)
  if (!is.null(coords) && ncol(coords) >= 2L) {
    meta[[coord.cols[1L]]] <- as.numeric(coords[, 1L])
    meta[[coord.cols[2L]]] <- as.numeric(coords[, 2L])
  }
  Seurat::CreateSeuratObject(
    counts = counts,
    assay = assay,
    meta.data = meta,
    project = project
  )
}
