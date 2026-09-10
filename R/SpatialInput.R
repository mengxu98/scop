# Inspect the analysis unit and coordinates of spatial data
#
# Resolve one spatial context without changing the object or materializing its
# expression matrix. Import provenance takes precedence over image-class
# evidence. Metadata-only coordinates have unknown physical units unless
# explicitly supplied. Selecting an image never implies that pixels are microns.
#
# @param object A Seurat object.
# @param assay Assay to inspect; NULL uses the default assay.
# @param image Image name. Multiple images require an explicit selection.
# @param coord.cols Coordinate columns for metadata-only input.
# @param data_type Observation type: auto, spot, bin, or cell. Auto returns
#   unknown when no reliable evidence is available.
# @param coordinate_units Optional explicit units: pixel, micron, or unknown.
# @return A list containing the assay, image, observation type, units, selected
#   cell IDs, dimensions, estimated dense matrix bytes, and coordinate source.
spatial_input_info <- function(object, assay = NULL, image = NULL,
                            coord.cols = c("col", "row"), data_type = "auto",
                            coordinate_units = NULL) {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object", call. = FALSE)
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  validate_scalar_string(assay, "assay")
  if (!assay %in% SeuratObject::Assays(object)) stop("assay is not present", call. = FALSE)
  data_type <- match.arg(data_type, c("auto", "spot", "bin", "cell"))
  resolved <- SpatialCoordinates(object, image = image, coord.cols = coord.cols, space = "raw")
  image <- resolved$source$image %||% image
  ids <- resolved$data$cell_id
  if (!length(ids) || anyNA(ids) || anyDuplicated(ids) || any(!nzchar(ids))) {
    stop("Spatial context must contain unique, nonmissing cell IDs", call. = FALSE)
  }
  assay_cells <- colnames(object[[assay]])
  if (!all(ids %in% assay_cells)) {
    stop("Selected image/context contains cells absent from assay; select the matching assay and resolution", call. = FALSE)
  }
  record <- object@misc$scop_spatial_input[[image %||% assay]]
  if (!is.null(record) && (!identical(record$assay, assay) ||
      !all(ids %in% record$cells))) record <- NULL
  inferred <- record$data_type %||% "unknown"
  units <- record$coordinate_units %||% "unknown"
  evidence <- if (is.null(record)) "unknown" else "import"
  if (is.null(record) && !is.null(image)) {
    img <- object[[image]]
    if (inherits(img, c("VisiumV1", "VisiumV2"))) {
      inferred <- if (grepl("[.]\\d+um$", assay)) "bin" else "spot"
      units <- "pixel"
      evidence <- "Visium image"
    } else if (!is.null(spatial_segmentation_name(img))) {
      inferred <- "cell"
      evidence <- "segmentation"
    }
  }
  if (data_type != "auto") {
    if (inferred != "unknown" && inferred != data_type) {
      stop("data_type conflicts with image/import evidence", call. = FALSE)
    }
    inferred <- data_type
    evidence <- paste(evidence, "explicit data_type", sep = "; ")
  }
  if (!is.null(coordinate_units)) {
    coordinate_units <- match.arg(coordinate_units, c("pixel", "micron", "unknown"))
    if (units != "unknown" && coordinate_units != units) {
      stop("coordinate_units conflicts with import/image evidence; transform coordinates explicitly", call. = FALSE)
    }
    units <- coordinate_units
  }
  list(assay = assay, image = image, data_type = inferred, coordinate_units = units,
    evidence = evidence, resolution_um = record$resolution_um %||% NA_real_,
    cells = ids, n_observations = length(ids), n_features = nrow(object[[assay]]),
    estimated_dense_bytes = 8 * as.double(length(ids)) * nrow(object[[assay]]),
    layers = SeuratObject::Layers(object[[assay]]), source = resolved$source)
}

#' Load and validate Visium, Visium HD, or Xenium data
#'
#' Delegate file parsing to Seurat, validate each image against its assay, and
#' record observation type, resolution and coordinate units. This function does
#' not download data, normalize counts, segment cells, or install packages.
#'
#' @param data.dir Existing vendor output directory.
#' @param technology One of visium, visium_hd, or xenium.
#' @param bin.size Positive integer HD bin sizes in microns. Ignored for other
#'   technologies; Seurat stores each resolution in a separate assay/image.
#' @param sample_id Nonempty sample identifier stored in spatial_sample.
#' @param ... Extra arguments to Seurat's Load10X_Spatial or LoadXenium. Xenium
#'   defaults to cell segmentation and does not load molecule coordinates unless
#'   requested; cell segmentation must already exist in the vendor output.
#' @return A validated Seurat object with misc$scop_spatial_input provenance.
#' @seealso SpatialCellPlot
#' @export
ReadSpatialData <- function(data.dir, technology = c("visium", "visium_hd", "xenium"),
                            bin.size = c(8L, 16L), sample_id = "sample1", ...) {
  technology <- match.arg(technology)
  validate_scalar_string(data.dir, "data.dir")
  validate_scalar_string(sample_id, "sample_id")
  if (!dir.exists(data.dir)) stop("data.dir does not exist", call. = FALSE)
  extra <- list(...)
  validate_named_list(extra, "...")
  standard_spatial_fixed_args(extra, c("data.dir", "bin.size"))
  if (technology == "visium_hd" && (!is.numeric(bin.size) || !length(bin.size) ||
      any(!is.finite(bin.size)) || any(bin.size <= 0 | bin.size != as.integer(bin.size)) || anyDuplicated(bin.size))) {
    stop("bin.size must contain unique positive integers", call. = FALSE)
  }
  if (technology == "xenium") {
    args <- merge_call_args(list(data.dir = data.dir, segmentations = "cell",
      molecule.coordinates = FALSE), extra)
    out <- do.call(Seurat::LoadXenium, args)
  } else {
    args <- merge_call_args(list(data.dir = data.dir,
      bin.size = if (technology == "visium_hd") bin.size else NULL), extra)
    out <- do.call(Seurat::Load10X_Spatial, args)
  }
  images <- SeuratObject::Images(out)
  if (!length(images)) stop("Loader returned no spatial images", call. = FALSE)
  records <- stats::setNames(vector("list", length(images)), images)
  for (image in images) {
    assay <- SeuratObject::DefaultAssay(out[[image]])
    info <- spatial_input_info(out, assay = assay, image = image)
    resolution <- if (technology == "visium_hd") {
      if (!grepl("[.]\\d+um$", assay)) stop("HD assay has no explicit bin resolution", call. = FALSE)
      as.numeric(sub(".*[.](\\d+)um$", "\\1", assay))
    } else NA_real_
    records[[image]] <- list(technology = technology, assay = assay, image = image,
      data_type = switch(technology, visium = "spot", visium_hd = "bin", xenium = "cell"),
      coordinate_units = if (technology == "xenium") "micron" else "pixel",
      resolution_um = resolution, cells = info$cells, sample_id = sample_id,
      loader = if (technology == "xenium") "Seurat::LoadXenium" else "Seurat::Load10X_Spatial",
      loader_version = as.character(utils::packageVersion("Seurat")),
      data.dir = normalizePath(data.dir, winslash = "/", mustWork = TRUE))
  }
  if (technology == "visium_hd" && !setequal(vapply(records, `[[`, numeric(1), "resolution_um"), bin.size)) {
    stop("Loader did not return every requested HD resolution", call. = FALSE)
  }
  out$spatial_sample <- sample_id
  out@misc$scop_spatial_input <- records
  out
}

standard_spatial_fixed_args <- function(args, managed) {
  bad <- intersect(names(args), managed)
  if (length(bad)) stop(paste("Arguments managed by the workflow:", paste(bad, collapse = ", ")), call. = FALSE)
  invisible(TRUE)
}
