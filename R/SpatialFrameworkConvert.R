#' @title Convert Seurat to a native Giotto object
#'
#' @description
#' Convert one Seurat spatial image into a native Giotto object without running
#' a Giotto workflow or changing the input object.
#'
#' @details
#' This bridge is separate from the legacy [SeuratToScopGiotto()] and
#' [RunGiottoWorkflow()] interface. It returns only the native Giotto object.
#' The scop-controlled converter is used by default to preserve exact
#' single-image selection across Giotto versions; pass `use_official = TRUE`
#' through `...` only when the installed GiottoClass converter is desired.
#'
#' @param srt A `Seurat` object.
#' @param image Seurat image name. Multi-image objects require an explicit name.
#' @param ... Additional arguments passed to [SeuratToScopGiotto()].
#'
#' @return A native Giotto object.
#'
#' @export
srt_to_giotto <- function(srt, image = NULL, ...) {
  image <- spatial_framework_image(srt, image = image)
  srt_use <- spatial_framework_subset_image(srt, image = image)
  extra <- list(...)
  if (!"use_official" %in% names(extra)) {
    extra$use_official <- FALSE
  }
  converted <- do.call(
    SeuratToScopGiotto,
    c(list(srt = srt_use, image = image), extra)
  )
  converted$giotto
}

#' @title Convert Giotto to Seurat
#'
#' @description Convert a native Giotto object into a new Seurat object.
#'
#' @details
#' The bridge selects the GiottoClass v4 or v5 converter according to the
#' installed Seurat major version. Older Giotto installations without these
#' converters should continue using the legacy scop Giotto workflow or be
#' updated before reverse conversion.
#'
#' @param giotto A native Giotto object.
#' @param ... Additional arguments passed to the Giotto converter.
#'
#' @return A new `Seurat` object.
#'
#' @export
giotto_to_srt <- function(giotto, ...) {
  giotto_require(verbose = FALSE)
  converter_name <- spatial_giotto_converter_name()
  converter <- tryCatch(
    get_namespace_fun("GiottoClass", converter_name),
    error = function(e) {
      log_message(
        "Installed Giotto does not provide {.fn {converter_name}} required for this Seurat version",
        message_type = "error"
      )
    }
  )
  giotto_call(converter, c(list(gobject = giotto), list(...)))
}

#' @title Convert Seurat to a native SPATA2 object
#'
#' @description
#' Convert one Seurat spatial image into a native SPATA2 object without running
#' SPATA2 analyses or changing the input object.
#'
#' @param srt A `Seurat` object.
#' @param image Seurat image name. Multi-image objects require an explicit name.
#' @param ... Additional arguments passed to `SPATA2::asSPATA2()`.
#'
#' @return A native SPATA2 object.
#'
#' @export
srt_to_spata2 <- function(srt, image = NULL, ...) {
  image <- spatial_framework_image(srt, image = image)
  srt_use <- spatial_framework_subset_image(srt, image = image)
  check_r("theMILOlab/SPATA2", verbose = FALSE)
  as_spata2 <- get_namespace_fun("SPATA2", "asSPATA2")
  extra <- list(...)
  if ("img_name" %in% names(extra)) {
    log_message("Use {.arg image} instead of passing {.arg img_name} through dots", message_type = "error")
  }
  if (!"sample_name" %in% names(extra)) {
    sample_name <- srt@project.name
    extra$sample_name <- if (is.null(sample_name) || !nzchar(sample_name)) "scop" else sample_name
  }
  args <- c(list(object = srt_use), if (is.null(image)) list() else list(img_name = image), extra)
  spata2_call_converter(as_spata2, args)
}

#' @title Convert SPATA2 to Seurat
#'
#' @description Convert a native SPATA2 object into a new Seurat object.
#'
#' @param spata2 A native SPATA2 object.
#' @param ... Additional arguments passed to `SPATA2::asSeurat()`.
#'
#' @return A new `Seurat` object.
#'
#' @export
spata2_to_srt <- function(spata2, ...) {
  check_r("theMILOlab/SPATA2", verbose = FALSE)
  as_seurat <- get_namespace_fun("SPATA2", "asSeurat")
  extra <- list(...)
  converted <- suppressWarnings(
    tryCatch(
      do.call(as_seurat, c(list(object = spata2), extra)),
      error = function(e) NULL
    )
  )
  if (inherits(converted, "Seurat")) {
    return(converted)
  }
  spata2_to_srt_fallback(spata2)
}

spatial_framework_image <- function(srt, image = NULL) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  if (!is.null(image) && (!is.character(image) || length(image) != 1L || is.na(image) || !nzchar(image))) {
    log_message("{.arg image} must be one non-empty image name", message_type = "error")
  }
  images <- tryCatch(SeuratObject::Images(srt), error = function(e) character())
  if (length(images) == 0L) {
    if (!is.null(image)) {
      log_message("{.arg image} was supplied but {.arg srt} has no spatial images", message_type = "error")
    }
    return(NULL)
  }
  if (is.null(image)) {
    if (length(images) > 1L) {
      log_message(
        "Multiple spatial images are available; select one with {.arg image}: {.val {images}}",
        message_type = "error"
      )
    }
    return(images[[1L]])
  }
  if (!image %in% images) {
    log_message(
      "{.arg image} {.val {image}} is not present in {.cls Seurat}; available images: {.val {images}}",
      message_type = "error"
    )
  }
  image
}

spatial_framework_subset_image <- function(srt, image = NULL) {
  if (is.null(image)) {
    return(srt)
  }
  coords <- as.data.frame(SeuratObject::GetTissueCoordinates(srt[[image]]))
  image_cells <- if ("cell" %in% colnames(coords)) {
    as.character(coords$cell)
  } else {
    rownames(coords)
  }
  image_cells <- intersect(image_cells, colnames(srt))
  if (length(image_cells) == 0L) {
    log_message("The selected image does not contain cells present in {.arg srt}", message_type = "error")
  }
  subset(srt, cells = image_cells)
}

spata2_call_converter <- function(fun, args) {
  # SPATA2 3.1.4 references its exported signatures dataset as a free
  # variable. Supply that binding only for the duration of the official
  # converter call and restore the user's global environment exactly.
  signatures <- tryCatch(
    getExportedValue("SPATA2", "signatures"),
    error = function(e) NULL
  )
  if (is.null(signatures)) {
    return(do.call(fun, args))
  }
  env <- globalenv()
  existed <- exists("signatures", envir = env, inherits = FALSE)
  old <- if (existed) get("signatures", envir = env, inherits = FALSE) else NULL
  assign("signatures", signatures, envir = env)
  on.exit({
    if (existed) {
      assign("signatures", old, envir = env)
    } else if (exists("signatures", envir = env, inherits = FALSE)) {
      rm(list = "signatures", envir = env)
    }
  }, add = TRUE)
  do.call(fun, args)
}

spata2_to_srt_fallback <- function(spata2) {
  get_count_matrix <- get_namespace_fun("SPATA2", "getCountMatrix")
  get_meta_df <- get_namespace_fun("SPATA2", "getMetaDf")
  get_coords_df <- get_namespace_fun("SPATA2", "getCoordsDf")
  get_sample_name <- get_namespace_fun("SPATA2", "getSampleName")

  counts <- get_count_matrix(spata2)
  meta <- as.data.frame(get_meta_df(spata2), stringsAsFactors = FALSE)
  coords <- as.data.frame(get_coords_df(spata2), stringsAsFactors = FALSE)
  sample_name <- as.character(get_sample_name(spata2))
  sample_name <- if (length(sample_name) == 0L || is.na(sample_name[[1L]]) || !nzchar(sample_name[[1L]])) {
    "SPATA2"
  } else {
    sample_name[[1L]]
  }

  first_col <- function(candidates, data) {
    hit <- intersect(candidates, colnames(data))
    if (length(hit) == 0L) NULL else hit[[1L]]
  }
  barcode_col <- first_col(c("barcodes", "barcode", "cell_id", "cell"), meta)
  if (is.null(barcode_col)) {
    if (nrow(meta) != ncol(counts)) {
      log_message("SPATA2 metadata could not be matched to the count matrix", message_type = "error")
    }
    rownames(meta) <- colnames(counts)
  } else {
    rownames(meta) <- as.character(meta[[barcode_col]])
  }
  if (anyNA(rownames(meta)) || any(!nzchar(rownames(meta))) || anyDuplicated(rownames(meta))) {
    log_message("SPATA2 metadata must contain unique, non-missing cell identifiers", message_type = "error")
  }
  cells <- intersect(colnames(counts), rownames(meta))
  if (length(cells) == 0L) {
    log_message("SPATA2 counts and metadata do not share cell identifiers", message_type = "error")
  }
  counts <- counts[, cells, drop = FALSE]
  meta <- meta[cells, , drop = FALSE]
  srt <- SeuratObject::CreateSeuratObject(
    counts = counts,
    project = sample_name,
    meta.data = meta
  )

  coord_barcode_col <- first_col(c("barcodes", "barcode", "cell_id", "cell"), coords)
  x_col <- first_col(c("x_orig", "x"), coords)
  y_col <- first_col(c("y_orig", "y"), coords)
  if (!is.null(coord_barcode_col) && !is.null(x_col) && !is.null(y_col)) {
    spatial <- data.frame(
      x = suppressWarnings(as.numeric(coords[[x_col]])),
      y = suppressWarnings(as.numeric(coords[[y_col]])),
      row.names = as.character(coords[[coord_barcode_col]]),
      stringsAsFactors = FALSE
    )
    spatial <- spatial[intersect(cells, rownames(spatial)), , drop = FALSE]
    spatial <- spatial[is.finite(spatial$x) & is.finite(spatial$y), , drop = FALSE]
    if (nrow(spatial) > 0L) {
      image_name <- gsub("[^A-Za-z0-9]+", "_", sample_name)
      image_name <- gsub("^_+|_+$", "", image_name)
      image_name <- if (nzchar(image_name)) image_name else "SPATA2"
      key <- paste0("spata", gsub("[^A-Za-z0-9]", "", tolower(image_name)), "_")
      srt[[image_name]] <- SeuratObject::CreateFOV(
        spatial,
        type = "centroids",
        assay = SeuratObject::DefaultAssay(srt),
        key = key
      )
    }
  }
  srt
}
