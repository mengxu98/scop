spatial_input_info <- function(object, assay = NULL, image = NULL,
                               coord.cols = c("col", "row"), data_type = "auto",
                               coordinate_units = NULL) {
  if (!inherits(object, "Seurat")) log_message("{.emph object} must be a Seurat object", message_type = "error")
  image <- spatial_image_resolve(object, image = image)$image
  assay <- assay %||% if (is.null(image)) SeuratObject::DefaultAssay(object) else SeuratObject::DefaultAssay(object[[image]])
  validate_scalar_string(assay, "assay")
  if (!assay %in% SeuratObject::Assays(object)) log_message("{.emph assay} is not present", message_type = "error")
  data_type <- match.arg(data_type, c("auto", "spot", "bin", "cell"))
  resolved <- SpatialCoordinates(object, image = image, coord.cols = coord.cols, space = "raw")
  image <- resolved$source$image %||% image
  ids <- resolved$data$cell_id
  if (!length(ids) || anyNA(ids) || anyDuplicated(ids) || any(!nzchar(ids))) {
    log_message("Spatial context must contain unique, nonmissing cell IDs", message_type = "error")
  }
  assay_cells <- colnames(object[[assay]])
  if (is.null(image)) ids <- ids[ids %in% assay_cells]
  if (!length(ids)) log_message("No selected assay cells have spatial coordinates", message_type = "error")
  if (!all(ids %in% assay_cells)) {
    log_message(
      "Selected image/context contains cells absent from assay; select the matching assay and resolution",
      message_type = "error"
    )
  }
  record <- object@misc$scop_spatial_input[[image %||% assay]]
  if (!is.null(record) && (!identical(record$assay, assay) ||
    !all(ids %in% record$cells))) {
    record <- NULL
  }
  if (is.null(record)) {
    candidates <- Filter(
      function(x) identical(x$assay, assay) && all(ids %in% x$cells),
      object@misc$scop_spatial_input
    )
    if (length(candidates) == 1L) record <- candidates[[1L]]
  }
  inferred <- record$data_type %||% "unknown"
  units <- record$coordinate_units %||% "unknown"
  evidence <- if (is.null(record)) "unknown" else "import"
  if (is.null(record) && !is.null(image)) {
    img <- object[[image]]
    if (inherits(img, c("VisiumV1", "VisiumV2"))) {
      inferred <- if (length(spatial_segmentation_names(img))) "cell" else if (grepl("[.]\\d+um$", assay)) "bin" else "spot"
      units <- "pixel"
      evidence <- "Visium image"
    } else if (length(spatial_segmentation_names(img))) {
      inferred <- "cell"
      evidence <- "segmentation"
    }
    chemistry <- object@misc$run_metadata$chemistry_version
    atera_source <- identical(object@misc$atera.molecules$fov, image) ||
      (length(SeuratObject::Images(object)) == 1L && is.character(chemistry) &&
        length(chemistry) == 1L && !is.na(chemistry) && grepl("^Atera ", chemistry))
    if (inherits(img, "FOV") && !inherits(img, "VisiumV2") && atera_source) {
      inferred <- "cell"
      units <- "micron"
      evidence <- "official Atera provenance"
    }
  }
  if (data_type != "auto") {
    if (inferred != "unknown" && inferred != data_type) {
      log_message("{.emph data_type} conflicts with image/import evidence", message_type = "error")
    }
    inferred <- data_type
    evidence <- paste(evidence, "explicit data_type", sep = "; ")
  }
  if (!is.null(coordinate_units)) {
    coordinate_units <- match.arg(coordinate_units, c("pixel", "micron", "unknown"))
    if (units != "unknown" && coordinate_units != units) {
      log_message(
        "{.emph coordinate_units} conflicts with import/image evidence; transform coordinates explicitly",
        message_type = "error"
      )
    }
    units <- coordinate_units
  }
  list(
    assay = assay, image = image, data_type = inferred, coordinate_units = units,
    evidence = evidence, resolution_um = record$resolution_um %||% NA_real_,
    cells = ids, n_observations = length(ids), n_features = nrow(object[[assay]]),
    estimated_dense_bytes = 8 * as.double(length(ids)) * nrow(object[[assay]]),
    layers = SeuratObject::Layers(object[[assay]]), source = resolved$source
  )
}

#' Load and validate Visium, Visium HD, Xenium, or Atera data
#'
#' Delegate file parsing to Seurat, validate each image against its assay, and
#' record observation type, resolution and coordinate units. This function does
#' not download data, normalize counts, segment cells, or install packages.
#'
#' @param data.dir Existing vendor output directory.
#' @param technology One of visium, visium_hd, xenium, or atera.
#' @param bin.size Positive integer HD bin sizes in microns, `"polygons"`, or a
#'   mixture such as `c(8, "polygons")`. Ignored for other
#'   technologies; Seurat stores each resolution in a separate assay/image.
#' @param sample_id Nonempty sample identifier stored in spatial_sample.
#' @param ... Arguments to the official Seurat loader. For HD, use
#'   `segmentation.type = "cell"` or `"nucleus"` and `compact = TRUE`.
#'   Atera retains official defaults, including deferred segmentations and
#'   molecules; use `bpcells.dir` for its counts cache. `genes` selects molecules,
#'   not expression-matrix rows. Additional Atera data can be loaded later with
#'   Seurat's molecule, segmentation, and codeword loaders.
#' @return A validated Seurat object with misc$scop_spatial_input provenance.
#' @seealso SpatialCellPlot
#' @export
ReadSpatialData <- function(data.dir, technology = c("visium", "visium_hd", "xenium", "atera"),
                            bin.size = c(8L, 16L), sample_id = "sample1", ...) {
  technology <- match.arg(technology)
  validate_scalar_string(data.dir, "data.dir")
  validate_scalar_string(sample_id, "sample_id")
  if (!dir.exists(data.dir)) log_message("{.emph data.dir} does not exist", message_type = "error")
  extra <- list(...)
  validate_named_list(extra, "...")
  standard_spatial_fixed_args(extra, c("data.dir", "bin.size"))
  if (technology == "visium_hd") {
    numeric_bins <- suppressWarnings(as.numeric(bin.size[bin.size != "polygons"]))
    if (!(is.numeric(bin.size) || is.character(bin.size)) || !length(bin.size) || anyNA(bin.size) ||
      any(!is.finite(numeric_bins)) || any(numeric_bins <= 0 | numeric_bins %% 1 != 0) ||
      anyDuplicated(numeric_bins) || sum(bin.size == "polygons") > 1L) {
      log_message("{.emph bin.size} must contain unique positive integers and optionally 'polygons'", message_type = "error")
    }
    if ("polygons" %in% bin.size && !"segmentation.type" %in% names(formals(Seurat::Load10X_Spatial))) {
      log_message("HD polygons require a Seurat version supporting Load10X_Spatial(segmentation.type)", message_type = "error")
    }
  }
  loader <- switch(technology,
    atera = "LoadAtera",
    xenium = "LoadXenium",
    "Load10X_Spatial"
  )
  if (!loader %in% getNamespaceExports("Seurat")) {
    log_message(paste0("This Seurat installation does not provide ", loader, "; Atera requires Seurat 6 or later"), message_type = "error")
  }
  if (technology == "atera") {
    out <- do.call(getExportedValue("Seurat", loader), c(list(data.dir = data.dir), extra))
  } else if (technology == "xenium") {
    args <- merge_call_args(list(
      data.dir = data.dir, segmentations = "cell",
      molecule.coordinates = FALSE
    ), extra)
    out <- do.call(Seurat::LoadXenium, args)
  } else {
    args <- merge_call_args(list(
      data.dir = data.dir,
      bin.size = if (technology == "visium_hd") bin.size else NULL
    ), extra)
    out <- tryCatch(do.call(Seurat::Load10X_Spatial, args), error = function(e) {
      if ("polygons" %in% bin.size && grepl("coords_x_orientation", conditionMessage(e), fixed = TRUE)) {
        log_message(paste0("The installed Seurat HD loader and VisiumV2 class are incompatible; use a matching Seurat installation (verified with Seurat 6). Upstream error: ", conditionMessage(e)), message_type = "error")
      }
      stop(e)
    })
  }
  images <- SeuratObject::Images(out)
  if (!length(images)) log_message("Loader returned no spatial images", message_type = "error")
  records <- stats::setNames(vector("list", length(images)), images)
  for (image in images) {
    assay <- SeuratObject::DefaultAssay(out[[image]])
    info <- spatial_input_info(out, assay = assay, image = image)
    segmented <- length(spatial_segmentation_names(out[[image]])) > 0L
    resolution <- if (technology == "visium_hd" && !segmented) {
      if (!grepl("[.]\\d+um$", assay)) log_message("HD assay has no explicit bin resolution", message_type = "error")
      as.numeric(sub(".*[.](\\d+)um$", "\\1", assay))
    } else {
      NA_real_
    }
    records[[image]] <- list(
      technology = technology, assay = assay, image = image,
      data_type = switch(technology,
        visium = "spot",
        visium_hd = if (segmented) "cell" else "bin",
        xenium = "cell",
        atera = "cell"
      ),
      coordinate_units = if (technology %in% c("xenium", "atera")) "micron" else "pixel",
      resolution_um = resolution, cells = info$cells, sample_id = sample_id,
      loader = paste0("Seurat::", loader),
      loader_version = as.character(utils::packageVersion("Seurat")),
      data.dir = normalizePath(data.dir, winslash = "/", mustWork = TRUE)
    )
  }
  if (technology == "visium_hd") {
    returned_bins <- vapply(records, `[[`, numeric(1), "resolution_um")
    returned_cells <- any(vapply(records, function(x) x$data_type == "cell", logical(1)))
    if (!setequal(returned_bins[!is.na(returned_bins)], numeric_bins) ||
      returned_cells != ("polygons" %in% bin.size)) {
      log_message("Loader did not return every requested HD resolution or segmentation", message_type = "error")
    }
  }
  out$spatial_sample <- sample_id
  out@misc$scop_spatial_input <- records
  out
}

standard_spatial_fixed_args <- function(args, managed) {
  bad <- intersect(names(args), managed)
  if (length(bad)) log_message(paste("Arguments managed by the workflow:", paste(bad, collapse = ", ")), message_type = "error")
  invisible(TRUE)
}
