#' @title Convert Seurat to a Giotto object
#'
#' @description
#' Convert one Seurat spatial image into a Giotto object without running
#' a Giotto workflow or changing the input object. Requires GiottoClass.
#'
#' @param srt Deprecated alias for `object`; supply exactly one of the two. It
#' will be removed in scop 1.0.0.
#' @param object A `Seurat` object.
#' @param image Seurat image name. Multi-image objects require an explicit name.
#' @param ... Arguments passed to other functions.
#' @return A Giotto object.
#'
#' @export
#' @examples
#' \dontrun{
#' data(visium_human_pancreas_sub)
#' spatial <- Seurat::NormalizeData(
#'   visium_human_pancreas_sub,
#'   assay = "Spatial",
#'   verbose = FALSE
#' )
#' srt_to_giotto(spatial, image = "slice1")
#' }
srt_to_giotto <- function(
  object,
  image = NULL,
  srt = NULL,
  ...
) {
  srt <- resolve_deprecated_srt(object, srt, missing(object))
  image <- {
    .img_resolved <- spatial_image_resolve(srt = srt, image = image, image_policy = "strict")
    .img_resolved$image
  }
  srt_use <- spatial_framework_subset_image(srt, image = image)
  ensure_giotto_bridge()
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
  with_giotto_bridge_options(do.call(converter, c(args, extra)))
}

seurat_major_version <- function() {
  suppressWarnings(as.integer(
    strsplit(as.character(utils::packageVersion("Seurat")), "\\.")[[1L]][1L]
  ))
}

#' @title Convert Giotto to Seurat
#'
#' @description Convert a Giotto object into a new Seurat object.
#' Requires GiottoClass.
#'
#' @md
#' @param giotto A Giotto object. For a round trip, normalize the Seurat object
#'   with [Seurat::NormalizeData()] before converting it to Giotto.
#' @param ... Arguments passed to other functions.
#' @return A new `Seurat` object.
#'
#' @export
#' @examples
#' \dontrun{
#' data(visium_human_pancreas_sub)
#' spatial <- Seurat::NormalizeData(
#'   visium_human_pancreas_sub,
#'   assay = "Spatial",
#'   verbose = FALSE
#' )
#' giotto <- srt_to_giotto(spatial, image = "slice1")
#' giotto_to_srt(giotto)
#' }
giotto_to_srt <- function(giotto, ...) {
  ensure_giotto_bridge()
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
  with_giotto_bridge_options(
    do.call(converter, c(list(gobject = giotto), list(...)))
  )
}

ensure_giotto_bridge <- function() {
  giotto_class <- check_r("GiottoClass", install = FALSE, verbose = FALSE)
  if (isTRUE(giotto_class[["GiottoClass"]])) {
    invisible(TRUE)
  } else {
    check_r("drieslab/Giotto", verbose = FALSE)
  }
}

with_giotto_bridge_options <- function(expr) {
  old <- options(
    giotto.use_conda = FALSE,
    giotto.check_version = FALSE
  )
  on.exit(options(old), add = TRUE)
  force(expr)
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
