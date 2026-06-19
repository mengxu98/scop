#' @title Run semla spatial network construction
#'
#' @description
#' Use the optional `semla` package as a backend to prepare a Staffli-enabled
#' Seurat object and compute spot-level spatial networks. The network is stored
#' in `srt@tools[[tool_name]]` when `store_results = TRUE`.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt A `Seurat` object with spatial image data.
#' @param image_type Image scale used by `semla::UpdateSeuratForSemla()` when
#' the object does not already contain a Staffli object.
#' @param nNeighbors Number of nearest spatial neighbors.
#' @param maxDist Optional maximum neighbor distance.
#' @param minK Minimum number of retained neighbors per spot.
#' @param coords Coordinate system passed to `semla::GetSpatialNetwork()`.
#' @param tool_name Name used to store results in `srt@tools`.
#' @param store_results Whether to store the semla spatial network in
#' `srt@tools`.
#' @param ... Additional arguments passed to semla.
#'
#' @return A `Seurat` object.
#' @export
RunSemlaSpatialNetwork <- function(
  srt,
  image_type = "tissue_lowres",
  nNeighbors = 6,
  maxDist = NULL,
  minK = 0,
  coords = "pixels",
  tool_name = "SemlaSpatialNetwork",
  store_results = TRUE,
  verbose = TRUE,
  ...
) {
  semla_validate_srt(srt)
  semla_require(verbose = verbose)
  image_type <- match.arg(image_type, c("tissue_lowres", "tissue_hires"))
  coords <- match.arg(coords, c("pixels", "array"))
  semla_validate_scalar_string(tool_name, "tool_name")

  srt <- semla_prepare_srt(
    srt = srt,
    image_type = image_type,
    verbose = verbose
  )
  spatial_network <- semla_get_fun("GetSpatialNetwork")(
    srt,
    nNeighbors = nNeighbors,
    maxDist = maxDist,
    minK = minK,
    coords = coords,
    ...
  )

  if (isTRUE(store_results)) {
    srt@tools[[tool_name]] <- list(
      network = spatial_network,
      parameters = list(
        image_type = image_type,
        nNeighbors = nNeighbors,
        maxDist = maxDist,
        minK = minK,
        coords = coords,
        tool_name = tool_name,
        store_results = store_results
      )
    )
  }
  log_message(
    "{.pkg semla} spatial network completed",
    message_type = "success",
    verbose = verbose
  )
  srt
}

#' @title Run semla local G spatial autocorrelation
#'
#' @description
#' Use `semla::RunLocalG()` on a Staffli-enabled Seurat object. Results are
#' written by semla to metadata or to an assay, depending on
#' `store_in_metadata`.
#'
#' @md
#' @inheritParams RunSemlaSpatialNetwork
#' @param features Features passed to `semla::RunLocalG()`.
#' @param alternative Alternative hypothesis passed to semla. Use `NULL` to
#' keep semla's default behavior.
#' @param store_in_metadata Whether semla should store results in metadata.
#' @param assay_name Assay name used by semla when `store_in_metadata = FALSE`.
#'
#' @return A `Seurat` object.
#' @export
RunSemlaLocalG <- function(
  srt,
  features,
  alternative = NULL,
  store_in_metadata = TRUE,
  assay_name = "GiScores",
  image_type = "tissue_lowres",
  verbose = TRUE,
  ...
) {
  semla_validate_srt(srt)
  semla_require(verbose = verbose)
  image_type <- match.arg(image_type, c("tissue_lowres", "tissue_hires"))
  srt <- semla_prepare_srt(
    srt = srt,
    image_type = image_type,
    verbose = verbose
  )
  semla_get_fun("RunLocalG")(
    srt,
    features = features,
    alternative = alternative,
    store_in_metadata = store_in_metadata,
    assay_name = assay_name,
    verbose = verbose,
    ...
  )
}

#' @title Run semla region neighbor detection
#'
#' @description
#' Use `semla::RegionNeighbors()` to identify neighboring spots for selected
#' metadata labels and write the returned columns to Seurat metadata.
#'
#' @md
#' @inheritParams RunSemlaSpatialNetwork
#' @param column_name Metadata column containing labels.
#' @param column_labels Labels to find neighbors for. If `NULL`, semla uses all
#' labels in `column_name`.
#' @param mode Neighbor selection mode passed to semla.
#' @param column_key Prefix for metadata columns returned by semla.
#'
#' @return A `Seurat` object.
#' @export
RunSemlaRegionNeighbors <- function(
  srt,
  column_name,
  column_labels = NULL,
  mode = "outer",
  column_key = NULL,
  image_type = "tissue_lowres",
  verbose = TRUE,
  ...
) {
  semla_validate_srt(srt)
  semla_require(verbose = verbose)
  image_type <- match.arg(image_type, c("tissue_lowres", "tissue_hires"))
  mode <- match.arg(mode, c("outer", "inner", "inner_outer", "all_inner_outer"))
  srt <- semla_prepare_srt(
    srt = srt,
    image_type = image_type,
    verbose = verbose
  )
  semla_get_fun("RegionNeighbors")(
    srt,
    column_name = column_name,
    column_labels = column_labels,
    mode = mode,
    column_key = column_key,
    verbose = verbose,
    ...
  )
}

#' @title Run semla radial distance analysis
#'
#' @description
#' Use `semla::RadialDistance()` to calculate distances from selected spatial
#' regions and write the returned columns to Seurat metadata.
#'
#' @md
#' @inheritParams RunSemlaSpatialNetwork
#' @param column_name Metadata column containing region labels.
#' @param selected_groups Region labels used by semla. If `NULL`, semla uses all
#' labels in `column_name`.
#' @param column_suffix Optional suffix for metadata columns returned by semla.
#'
#' @return A `Seurat` object.
#' @export
RunSemlaRadialDistance <- function(
  srt,
  column_name,
  selected_groups = NULL,
  column_suffix = NULL,
  image_type = "tissue_lowres",
  verbose = TRUE,
  ...
) {
  semla_validate_srt(srt)
  semla_require(verbose = verbose)
  image_type <- match.arg(image_type, c("tissue_lowres", "tissue_hires"))
  srt <- semla_prepare_srt(
    srt = srt,
    image_type = image_type,
    verbose = verbose
  )
  semla_get_fun("RadialDistance")(
    srt,
    column_name = column_name,
    selected_groups = selected_groups,
    column_suffix = column_suffix,
    verbose = verbose,
    ...
  )
}

semla_prepare_srt <- function(
  srt,
  image_type = "tissue_lowres",
  verbose = TRUE
) {
  semla_validate_srt(srt)
  semla_require(verbose = verbose)
  image_type <- match.arg(image_type, c("tissue_lowres", "tissue_hires"))
  if (!is.null(srt@tools[["Staffli"]])) {
    return(srt)
  }
  log_message(
    "Preparing {.cls Seurat} object for {.pkg semla}",
    verbose = verbose
  )
  semla_get_fun("UpdateSeuratForSemla")(
    srt,
    image_type = image_type,
    verbose = verbose
  )
}

semla_require <- function(verbose = TRUE) {
  if (semla_pkg_available()) {
    return(invisible(TRUE))
  }
  log_message(
    paste(
      "The optional R package {.pkg semla} is required.",
      "Install it with {.code thisutils::check_r('spatial-research/semla')}",
      "or {.code pak::pak('spatial-research/semla')}."
    ),
    message_type = "error",
    verbose = verbose
  )
}

semla_pkg_available <- function() {
  requireNamespace("semla", quietly = TRUE)
}

semla_get_fun <- function(fun) {
  semla_require(verbose = FALSE)
  tryCatch(
    getExportedValue("semla", fun),
    error = function(e) {
      log_message(
        "{.pkg semla} does not export {.fn {fun}} in the installed version",
        message_type = "error"
      )
    }
  )
}

semla_validate_srt <- function(srt) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

semla_validate_scalar_string <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    log_message(
      paste0(arg, " must be a non-empty string"),
      message_type = "error"
    )
  }
  invisible(TRUE)
}
