#' @title Run semla spatial network construction
#'
#' @description
#' Use the optional `semla` package as a backend to prepare a Staffli-enabled
#' Seurat object and compute spot-level spatial networks. The network is stored
#' in `srt@tools[[tool_name]]` when `store_results = TRUE`.
#' SCOP provides no dedicated plot for this result; retrieve it with
#' [GetSpatialResult()] and use an existing generic spatial plot when needed.
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
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' spatial@tools$SemlaSpatialNetwork <- list(
#'   network = data.frame(
#'     from = colnames(spatial)[1:6],
#'     to = colnames(spatial)[2:7],
#'     distance = sqrt(diff(spatial$x[1:7])^2 + diff(spatial$y[1:7])^2)
#'   )
#' )
#'
#' head(spatial@tools$SemlaSpatialNetwork$network)
#' SpatialSpotPlot(
#'   spatial,
#'   group.by = "coda_label",
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#'
#' spatial <- RunSemlaSpatialNetwork(
#'   spatial,
#'   nNeighbors = 6,
#'   coords = "pixels",
#'   verbose = FALSE
#' )
#' @concept spatial-producer
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
    srt@tools[[tool_name]] <- spatial_result_build(
      bundle = list(network = spatial_network, cells = colnames(srt)),
      method = "SemlaSpatialNetwork",
      result_type = "neighborhood",
      source = semla_spatial_source(srt, coords = coords),
      provenance = list(producer = "RunSemlaSpatialNetwork", backend_id = "semla"),
      parameters = list(
        image_type = image_type,
        nNeighbors = nNeighbors,
        maxDist = maxDist,
        minK = minK,
        coords = coords,
        tool_name = tool_name,
        store_results = store_results
      ),
      summary = list(n_networks = length(spatial_network))
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
#' SCOP provides no dedicated plot for this result; retrieve its schema record
#' with [GetSpatialResult()] and inspect the recorded output columns or assay.
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
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' spatial <- Seurat::NormalizeData(spatial, assay = "Spatial", verbose = FALSE)
#' features <- rownames(spatial)[1:3]
#' spatial[[paste0(features[1], "_localG")]] <- as.numeric(scale(spatial$x))
#'
#' SpatialSpotPlot(
#'   spatial,
#'   group.by = paste0(features[1], "_localG"),
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#'
#' spatial <- RunSemlaLocalG(
#'   spatial,
#'   features = features,
#'   store_in_metadata = TRUE,
#'   verbose = FALSE
#' )
#' @concept spatial-producer
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
  backend_args <- list(...)
  before_metadata <- colnames(srt@meta.data)
  out <- do.call(
    semla_get_fun("RunLocalG"),
    c(
      list(
        object = srt,
        features = features,
        alternative = alternative,
        store_in_metadata = store_in_metadata,
        assay_name = assay_name,
        verbose = verbose
      ),
      backend_args
    )
  )
  out@tools[["SemlaLocalG"]] <- spatial_result_build(
    bundle = list(
      output_columns = setdiff(colnames(out@meta.data), before_metadata),
      cells = colnames(out)
    ),
    method = "SemlaLocalG", result_type = "neighborhood",
    source = semla_spatial_source(
      out,
      coords = semla_backend_coords(backend_args)
    ),
    provenance = list(producer = "RunSemlaLocalG", backend_id = "semla"),
    parameters = list(
      features = features,
      alternative = alternative,
      store_in_metadata = store_in_metadata,
      assay_name = assay_name,
      image_type = image_type,
      backend_args = backend_args
    )
  )
  out
}

#' @title Run semla region neighbor detection
#'
#' @description
#' Use `semla::RegionNeighbors()` to identify neighboring spots for selected
#' metadata labels and write the returned columns to Seurat metadata.
#' SCOP provides no dedicated plot for this result; retrieve its schema record
#' with [GetSpatialResult()] and inspect the recorded metadata columns.
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
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' spatial$region <- ifelse(
#'   spatial$x > stats::median(spatial$x),
#'   "right",
#'   "left"
#' )
#' spatial$right_border <- spatial$region == "right" &
#'   abs(spatial$x - stats::median(spatial$x)) < stats::sd(spatial$x) * 0.25
#'
#' SpatialSpotPlot(
#'   spatial,
#'   group.by = c("region", "right_border"),
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#'
#' spatial <- RunSemlaRegionNeighbors(
#'   spatial,
#'   column_name = "region",
#'   column_labels = "right",
#'   mode = "outer",
#'   column_key = "right_border",
#'   verbose = FALSE
#' )
#' @concept spatial-producer
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
  backend_args <- list(...)
  before_metadata <- colnames(srt@meta.data)
  out <- do.call(
    semla_get_fun("RegionNeighbors"),
    c(
      list(
        object = srt,
        column_name = column_name,
        column_labels = column_labels,
        mode = mode,
        column_key = column_key,
        verbose = verbose
      ),
      backend_args
    )
  )
  out@tools[["SemlaRegionNeighbors"]] <- spatial_result_build(
    bundle = list(
      output_columns = setdiff(colnames(out@meta.data), before_metadata),
      cells = colnames(out)
    ),
    method = "SemlaRegionNeighbors", result_type = "neighborhood",
    source = semla_spatial_source(
      out,
      coords = semla_backend_coords(backend_args)
    ),
    provenance = list(producer = "RunSemlaRegionNeighbors", backend_id = "semla"),
    parameters = list(
      column_name = column_name,
      column_labels = column_labels,
      mode = mode,
      column_key = column_key,
      image_type = image_type,
      backend_args = backend_args
    )
  )
  out
}

#' @title Run semla radial distance analysis
#'
#' @description
#' Use `semla::RadialDistance()` to calculate distances from selected spatial
#' regions and write the returned columns to Seurat metadata.
#' SCOP provides no dedicated plot for this result; retrieve its schema record
#' with [GetSpatialResult()] and inspect the recorded metadata columns.
#'
#' @md
#' @inheritParams RunSemlaSpatialNetwork
#' @param column_name Metadata column containing region labels.
#' @param selected_groups Region labels used by semla. If `NULL`, semla uses all
#' labels in `column_name`.
#' @param column_suffix Optional suffix for metadata columns returned by semla.
#'
#' @return A `Seurat` object.
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' spatial$region <- ifelse(
#'   spatial$y > stats::median(spatial$y),
#'   "upper",
#'   "lower"
#' )
#' upper_center <- c(
#'   stats::median(spatial$x[spatial$region == "upper"]),
#'   stats::median(spatial$y[spatial$region == "upper"])
#' )
#' spatial$upper_distance <- sqrt(
#'   (spatial$x - upper_center[1])^2 + (spatial$y - upper_center[2])^2
#' )
#'
#' SpatialSpotPlot(
#'   spatial,
#'   group.by = "upper_distance",
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#'
#' spatial <- RunSemlaRadialDistance(
#'   spatial,
#'   column_name = "region",
#'   selected_groups = "upper",
#'   column_suffix = "upper_distance",
#'   verbose = FALSE
#' )
#' @concept spatial-producer
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
  backend_args <- list(...)
  before_metadata <- colnames(srt@meta.data)
  out <- do.call(
    semla_get_fun("RadialDistance"),
    c(
      list(
        object = srt,
        column_name = column_name,
        selected_groups = selected_groups,
        column_suffix = column_suffix,
        verbose = verbose
      ),
      backend_args
    )
  )
  out@tools[["SemlaRadialDistance"]] <- spatial_result_build(
    bundle = list(
      output_columns = setdiff(colnames(out@meta.data), before_metadata),
      cells = colnames(out)
    ),
    method = "SemlaRadialDistance", result_type = "neighborhood",
    source = semla_spatial_source(
      out,
      coords = "pixels",
      output_unit = if (isTRUE(backend_args$convert_to_microns)) {
        "micrometre"
      } else {
        "full_resolution_pixel"
      }
    ),
    provenance = list(producer = "RunSemlaRadialDistance", backend_id = "semla"),
    parameters = list(
      column_name = column_name,
      selected_groups = selected_groups,
      column_suffix = column_suffix,
      image_type = image_type,
      backend_args = backend_args
    )
  )
  out
}

semla_backend_coords <- function(backend_args) {
  coords <- backend_args$coords %||% "pixels"
  if (
    !is.character(coords) || length(coords) != 1L || is.na(coords) ||
      !coords %in% c("pixels", "array")
  ) {
    log_message(
      "{.arg coords} must be one of {.val pixels} or {.val array}",
      message_type = "error"
    )
  }
  coords
}

semla_spatial_source <- function(
  srt,
  coords = c("pixels", "array"),
  output_unit = NULL
) {
  coords <- match.arg(coords)
  images <- tryCatch(SeuratObject::Images(srt), error = function(e) character())
  staffli <- srt@tools[["Staffli"]]
  staffli_metadata <- tryCatch(
    if (isS4(staffli)) {
      methods::slot(staffli, "meta_data")
    } else {
      staffli$meta_data
    },
    error = function(e) NULL
  )
  sample_ids <- tryCatch(
    unique(as.character(staffli_metadata$sampleID)),
    error = function(e) character()
  )
  sample_ids <- sample_ids[!is.na(sample_ids) & nzchar(sample_ids)]
  list(
    image = if (length(images) == 1L) as.character(images) else NA_character_,
    images = images,
    image_policy = "native_multi_image",
    coordinate_space = "raw",
    coord.cols = if (identical(coords, "array")) {
      c("x", "y")
    } else {
      c("pxl_col_in_fullres", "pxl_row_in_fullres")
    },
    unit = if (identical(coords, "array")) {
      "array_index"
    } else {
      "full_resolution_pixel"
    },
    output_unit = output_unit %||% if (identical(coords, "array")) {
      "array_index"
    } else {
      "full_resolution_pixel"
    },
    selection_strategy = if (length(images) > 1L) {
      "all_images_partitioned_by_staffli_sampleID"
    } else if (length(images) == 1L) {
      "single_image_staffli"
    } else {
      "existing_staffli_metadata"
    },
    sample_ids = sample_ids
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
  status <- tryCatch(
    check_r("spatial-research/semla", verbose = FALSE),
    error = function(e) FALSE
  )
  if (isTRUE(unname(unlist(status))[1])) {
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

semla_get_fun <- function(fun) {
  semla_require(verbose = FALSE)
  tryCatch(
    get_namespace_fun("semla", fun),
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
