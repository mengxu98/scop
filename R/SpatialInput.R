#' Inspect the analysis unit and coordinates of spatial data
#'
#' Resolve one spatial context without changing the object or materializing its
#' expression matrix. Import provenance takes precedence over image-class
#' evidence. Metadata-only coordinates have unknown physical units unless
#' explicitly supplied. Selecting an image never implies that pixels are microns.
#'
#' @param object A Seurat object.
#' @param assay Assay to inspect; NULL uses the default assay.
#' @param image Image name. Multiple images require an explicit selection.
#' @param coord.cols Coordinate columns for metadata-only input.
#' @param data_type Observation type: auto, spot, bin, or cell. Auto returns
#'   unknown when no reliable evidence is available.
#' @param coordinate_units Optional explicit units: pixel, micron, or unknown.
#' @return A list containing the assay, image, observation type, units, selected
#'   cell IDs, dimensions, estimated dense matrix bytes, and coordinate source.
#' @export
SpatialDataInfo <- function(object, assay = NULL, image = NULL,
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
#' @seealso SpatialDataInfo, SpatialSegmentationQC, RunSpatialSketch
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
    info <- SpatialDataInfo(out, assay = assay, image = image)
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

#' Inspect segmentation area and transcript counts per cell
#'
#' Return descriptive cell-level QC without filtering or changing segmentation.
#' Missing boundaries are retained as not_evaluated. Polygon rings are validated
#' using the same boundary contract as SpatialCellPlot. Areas use squared raw
#' coordinate units, not necessarily square microns. No universal area cutoff
#' is imposed; thresholds must be chosen for the acquisition and tissue.
#'
#' @param object A Seurat object.
#' @param assay,image,coord.cols See SpatialDataInfo.
#' @param boundaries Optional real boundary table, as accepted by SpatialCellPlot.
#' @param min_counts Minimum transcripts; NULL does not evaluate this criterion.
#' @param area_range Optional finite c(min, max) area thresholds in squared raw units.
#' @return A data.frame of cell IDs, counts, detected genes, area, boundary status,
#'   and QC status. Attributes retain coordinate source and parameters. This is
#'   descriptive QC; it does not detect every segmentation error.
#'   Vendor boundary polygons may simplify the true segmentation masks; area
#'   refers to the supplied polygon representation.
#' @export
SpatialSegmentationQC <- function(object, assay = NULL, image = NULL,
                                  coord.cols = c("x", "y"), boundaries = NULL,
                                  min_counts = NULL, area_range = NULL) {
  info <- SpatialDataInfo(object, assay, image, coord.cols)
  if (!is.null(min_counts) && (length(min_counts) != 1L || !is.numeric(min_counts) ||
      !is.finite(min_counts) || min_counts < 0)) stop("min_counts must be nonnegative", call. = FALSE)
  if (!is.null(area_range) && (!is.numeric(area_range) || length(area_range) != 2L ||
      any(!is.finite(area_range)) || area_range[1] < 0 || area_range[2] <= area_range[1])) {
    stop("area_range must contain increasing nonnegative finite bounds", call. = FALSE)
  }
  counts <- GetAssayData5(object, assay = info$assay, layer = "counts", cells = info$cells)
  if (!identical(colnames(counts), info$cells) || !nrow(counts)) stop("Counts do not match selected cells", call. = FALSE)
  if (is.null(boundaries) && !is.null(info$image)) {
    boundaries <- spatial_segmentation_table(object[[info$image]])
  }
  area <- stats::setNames(rep(NA_real_, length(info$cells)), info$cells)
  if (!is.null(boundaries)) {
    b <- spatial_boundary_validate(boundaries, image = info$image)
    b <- b[b$cell_id %in% info$cells, , drop = FALSE]
    if (!is.null(info$image)) b <- b[is.na(b$image) | b$image == info$image, , drop = FALSE]
    multiple_rings <- vapply(split(b$ring_id, interaction(b$cell_id, b$polygon_id, drop = TRUE)),
      function(x) length(unique(x)) > 1L, logical(1))
    if (any(multiple_rings) && !"hole" %in% names(b)) {
      stop("Multi-ring area requires an explicit logical hole column", call. = FALSE)
    }
    if ("hole" %in% names(b) && (!is.logical(b$hole) || anyNA(b$hole))) {
      stop("hole must be logical and nonmissing", call. = FALSE)
    }
    rings <- split(b, interaction(b$cell_id, b$polygon_id, b$ring_id, drop = TRUE))
    ring_area <- lapply(rings, function(r) {
      if ("hole" %in% names(r) && length(unique(r$hole)) != 1L) stop("hole must be constant within each ring", call. = FALSE)
      next_i <- c(seq_len(nrow(r))[-1L], 1L)
      value <- abs(sum(r$x * r$y[next_i] - r$y * r$x[next_i])) / 2
      data.frame(cell_id = r$cell_id[1], area = if (isTRUE(r$hole[1])) -value else value)
    })
    if (length(ring_area)) {
      sums <- stats::aggregate(area ~ cell_id, do.call(rbind, ring_area), sum)
      area[sums$cell_id] <- sums$area
    }
  }
  out <- data.frame(cell_id = info$cells, counts = as.numeric(Matrix::colSums(counts)),
    detected_genes = as.numeric(Matrix::colSums(counts > 0)), area = unname(area),
    boundary_status = ifelse(is.na(area), "missing", ifelse(area > 0, "available", "invalid")),
    stringsAsFactors = FALSE)
  failed <- !is.na(out$area) & out$area <= 0
  if (!is.null(min_counts)) failed <- failed | out$counts < min_counts
  if (!is.null(area_range)) failed <- failed | (!is.na(out$area) &
    (out$area < area_range[1] | out$area > area_range[2]))
  evaluated <- (!is.null(min_counts) || !is.null(area_range)) & !is.na(out$area)
  out$status <- ifelse(failed, "fail", ifelse(evaluated, "pass", "not_evaluated"))
  attr(out, "source") <- info
  attr(out, "parameters") <- list(min_counts = min_counts, area_range = area_range)
  out
}
