#' @title Convert Seurat to a native Giotto object
#'
#' @description
#' Convert one Seurat spatial image into a native Giotto object without running
#' a Giotto workflow or changing the input object.
#'
#' @param srt A `Seurat` object.
#' @param image Seurat image name. Multi-image objects require an explicit name.
#' @param ... Additional arguments passed to the Giotto converter.
#'
#' @return A native Giotto object.
#'
#' @export
srt_to_giotto <- function(srt, image = NULL, ...) {
  image <- {
    .img_resolved <- spatial_image_resolve(srt = srt, image = image, image_policy = "strict")
    .img_resolved$image
  }
  srt_use <- spatial_framework_subset_image(srt, image = image)
  check_r("drieslab/Giotto", verbose = FALSE)
  converter_name <- if (seurat_major_version() >= 5L) {
    "seuratToGiottoV5"
  } else {
    "seuratToGiottoV4"
  }
  converter <- tryCatch(
    get_namespace_fun("GiottoClass", converter_name),
    error = function(e) {
      log_message(
        "Installed Giotto does not provide {.fn {converter_name}} required for this Seurat version",
        message_type = "error"
      )
    }
  )
  extra <- list(...)
  if ("sobject" %in% names(extra)) {
    log_message(
      "{.arg sobject} is managed by {.fn srt_to_giotto} and cannot be supplied through {.arg ...}",
      message_type = "error"
    )
  }
  args <- list(sobject = srt_use)
  if (!"verbose" %in% names(extra)) {
    args$verbose <- FALSE
  }
  if (!is.null(image) && !"spatial_assay" %in% names(extra)) {
    args$spatial_assay <- tryCatch(
      SeuratObject::DefaultAssay(srt_use[[image]]),
      error = function(e) SeuratObject::DefaultAssay(srt_use)
    )
  }
  do.call(converter, c(args, extra))
}

seurat_major_version <- function() {
  suppressWarnings(as.integer(
    strsplit(as.character(utils::packageVersion("Seurat")), "\\.")[[1L]][1L]
  ))
}

#' @title Convert Giotto to Seurat
#'
#' @description Convert a native Giotto object into a new Seurat object.
#'
#' @details
#' The bridge selects the GiottoClass v4 or v5 converter according to the
#' installed Seurat major version. For round trips produced by
#' [srt_to_giotto()], normalize the Seurat input first (for example with
#' [Seurat::NormalizeData()]); the official GiottoClass converter records an
#' empty normalized layer otherwise, which its reverse converter cannot read
#' back.
#'
#' @param giotto A native Giotto object.
#' @param ... Additional arguments passed to the Giotto converter.
#'
#' @return A new `Seurat` object.
#'
#' @export
giotto_to_srt <- function(giotto, ...) {
  check_r("drieslab/Giotto", verbose = FALSE)
  converter_name <- if (seurat_major_version() >= 5L) {
    "giottoToSeuratV5"
  } else {
    "giottoToSeuratV4"
  }
  converter <- tryCatch(
    get_namespace_fun("GiottoClass", converter_name),
    error = function(e) {
      log_message(
        "Installed Giotto does not provide {.fn {converter_name}} required for this Seurat version",
        message_type = "error"
      )
    }
  )
  do.call(converter, c(list(gobject = giotto), list(...)))
}

spatial_framework_subset_image <- function(srt, image = NULL) {
  if (is.null(image)) {
    return(srt)
  }
  coords <- spatial_coords_raw(
    srt = srt,
    image = image,
    image_policy = "strict"
  )$data
  image_cells <- as.character(coords$cell_id)
  image_cells <- intersect(image_cells, colnames(srt))
  if (length(image_cells) == 0L) {
    log_message("The selected image does not contain cells present in {.arg srt}", message_type = "error")
  }
  subset(srt, cells = image_cells)
}
