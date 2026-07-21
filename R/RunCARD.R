#' @title Run CARD spatial deconvolution
#'
#' @description
#' Estimate spot-level cell-type proportions from a spatial `Seurat` object
#' using a single-cell `Seurat` reference and the optional `CARD`/`CARDspa`
#' backend.
#'
#' @md
#' @inheritParams RunRCTD
#' @param sample_varname Optional metadata column in `reference` containing
#' sample labels. When `NULL`, all reference cells are assigned to one sample.
#' @param minCountGene,minCountSpot Filtering parameters passed to
#' `CARD::createCARDObject()` or `CARDspa::createCARDObject()` when supported.
#' @param ct_select Optional cell types to keep in CARD.
#' @param tool_name Name used to store detailed results in `srt@tools`.
#' @param create_card_params Additional parameters passed to
#' `createCARDObject()`.
#' @param card_deconvolution_params Additional parameters passed to
#' `CARD_deconvolution()`.
#'
#' @return A `Seurat` object with CARD proportion columns in metadata and
#' dominant cell type summaries. When `store_results = TRUE`, detailed results
#' are stored in `srt@tools[[tool_name]]`.
#' @concept spatial-producer
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' card_weights <- data.frame(
#'   CARD_prop_Ductal = seq(0.70, 0.20, length.out = ncol(spatial)),
#'   CARD_prop_Endocrine = seq(0.20, 0.70, length.out = ncol(spatial)),
#'   CARD_prop_Stromal = 0.10,
#'   row.names = colnames(spatial)
#' )
#' card_weights <- card_weights / rowSums(card_weights)
#' spatial <- Seurat::AddMetaData(spatial, card_weights)
#' spatial$CARD_dominant_type <- sub(
#'   "^CARD_prop_",
#'   "",
#'   colnames(card_weights)[max.col(card_weights)]
#' )
#' spatial$CARD_max_prop <- apply(card_weights, 1, max)
#'
#' SpatialSpotPlot(
#'   spatial,
#'   group.by = "CARD_dominant_type",
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#'   SpatialSpotPlot(
#'     spatial,
#'     group.by = "CARD_dominant_type",
#'     plot_type = "pie",
#'     overlay_image = FALSE,
#'     coord.cols = c("x", "y")
#'   )
#'
#' data(pancreas_sub)
#' features_use <- head(intersect(rownames(spatial), rownames(pancreas_sub)), 300)
#' spatial <- RunCARD(
#'   spatial,
#'   reference = pancreas_sub,
#'   reference_label = "CellType",
#'   assay = "Spatial",
#'   reference_assay = "RNA",
#'   features = features_use,
#'   verbose = FALSE
#' )
RunCARD <- function(
  srt,
  reference,
  reference_label,
  assay = NULL,
  reference_assay = NULL,
  layer = "counts",
  reference_layer = "counts",
  features = NULL,
  image = NULL,
  coord.cols = c("col", "row"),
  sample_varname = NULL,
  minCountGene = 100,
  minCountSpot = 5,
  ct_select = NULL,
  prefix = "CARD",
  tool_name = "CARD",
  store_results = TRUE,
  round_counts = TRUE,
  create_card_params = list(),
  card_deconvolution_params = list(),
  verbose = TRUE,
  ...,
  coordinate_space = c("raw", "legacy_display")
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (!inherits(reference, "Seurat")) {
    log_message(
      "{.arg reference} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  card_assert_scalar_string(prefix, "prefix")
  card_assert_scalar_string(tool_name, "tool_name")
  card_validate_param_list(create_card_params, "create_card_params")
  card_validate_param_list(card_deconvolution_params, "card_deconvolution_params")
  coordinate_space <- match.arg(coordinate_space)

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  reference_assay <- reference_assay %||% SeuratObject::DefaultAssay(reference)

  labels <- rctd_resolve_reference_labels(reference, reference_label)
  names(labels) <- colnames(reference)
  keep_ref <- !is.na(labels) & nzchar(as.character(labels))
  if (!all(keep_ref)) {
    log_message(
      "Drop {.val {sum(!keep_ref)}} reference cells with missing {.arg reference_label}",
      verbose = verbose
    )
    reference <- reference[, keep_ref]
    labels <- labels[keep_ref]
    names(labels) <- colnames(reference)
  }
  labels <- factor(as.character(labels), levels = unique(as.character(labels)))
  if (length(levels(labels)) == 0L) {
    log_message(
      "{.arg reference_label} must contain at least one non-missing class",
      message_type = "error"
    )
  }

  features_use <- rctd_common_features(
    srt = srt,
    reference = reference,
    assay = assay,
    reference_assay = reference_assay,
    features = features
  )
  if (length(features_use) == 0L) {
    log_message(
      "No shared features are available for {.fn RunCARD}",
      message_type = "error"
    )
  }
  st_counts <- rctd_get_count_matrix(
    srt,
    assay = assay,
    layer = layer,
    features = features_use,
    data_label = "Spatial",
    round_counts = round_counts,
    verbose = verbose
  )
  ref_counts <- rctd_get_count_matrix(
    reference,
    assay = reference_assay,
    layer = reference_layer,
    features = features_use,
    data_label = "Reference",
    round_counts = round_counts,
    verbose = verbose
  )
  keep_features <- Matrix::rowSums(st_counts) > 0 & Matrix::rowSums(ref_counts) > 0
  if (!any(keep_features)) {
    log_message(
      "No shared features have non-zero counts in both spatial and reference data",
      message_type = "error"
    )
  }
  if (sum(keep_features) < length(features_use)) {
    log_message(
      "Use {.val {sum(keep_features)}} shared non-zero features for CARD",
      verbose = verbose
    )
  }
  st_counts <- st_counts[keep_features, , drop = FALSE]
  ref_counts <- ref_counts[keep_features, , drop = FALSE]

  keep_spots <- Matrix::colSums(st_counts) > 0
  keep_ref_cells <- Matrix::colSums(ref_counts) > 0
  st_counts <- st_counts[, keep_spots, drop = FALSE]
  ref_counts <- ref_counts[, keep_ref_cells, drop = FALSE]
  labels <- labels[keep_ref_cells]
  names(labels) <- colnames(ref_counts)
  if (ncol(st_counts) == 0L || ncol(ref_counts) == 0L) {
    log_message(
      "No spatial spots or reference cells remain after zero-count filtering",
      message_type = "error"
    )
  }

  coords <- rctd_get_spatial_coords(
    srt = srt,
    spot_ids = colnames(st_counts),
    image = image,
    coord.cols = coord.cols,
    coordinate_space = coordinate_space
  )
  ref_meta <- card_reference_metadata(
    reference = reference,
    labels = labels,
    sample_varname = sample_varname
  )
  ct_select <- ct_select %||% levels(labels)

  extra_params <- list(...)
  if (length(extra_params) > 0L) {
    card_validate_param_list(extra_params, "...")
    card_deconvolution_params <- c(card_deconvolution_params, extra_params)
  }

  log_message(
    "Run {.pkg CARD} with {.val {nrow(st_counts)}} features, {.val {ncol(st_counts)}} spatial spots, and {.val {ncol(ref_counts)}} reference cells",
    verbose = verbose
  )
  backend <- card_run_backend(
    st_counts = st_counts,
    ref_counts = ref_counts,
    ref_meta = ref_meta,
    coords = coords,
    ct_select = ct_select,
    minCountGene = minCountGene,
    minCountSpot = minCountSpot,
    create_card_params = create_card_params,
    card_deconvolution_params = card_deconvolution_params
  )

  weights <- card_orient_weights(
    weights = backend$weights,
    spot_ids = colnames(st_counts)
  )
  weight_summary <- scop_spatial_finalize_weights(
    weights = weights,
    all_spots = colnames(srt)
  )
  weights <- weight_summary$weights
  srt <- scop_spatial_add_deconv_metadata(
    srt,
    weights = weights,
    prefix = prefix,
    metadata = weight_summary
  )

  if (isTRUE(store_results)) {
    srt@tools[[tool_name]] <- list(
      weights = weights,
      proportions = weight_summary$full_weights,
      cells = colnames(srt),
      coords = coords,
      features = rownames(st_counts),
      reference_metadata = ref_meta,
      backend_package = backend$package,
      object = backend$object,
      summary = scop_spatial_weight_summary(weights),
      parameters = list(
        assay = assay,
        reference_assay = reference_assay,
        layer = layer,
        reference_layer = reference_layer,
        reference_label = reference_label,
        image = image,
        coord.cols = coord.cols,
        coordinate_space = coordinate_space,
        sample_varname = sample_varname,
        minCountGene = minCountGene,
        minCountSpot = minCountSpot,
        ct_select = ct_select,
        prefix = prefix,
        tool_name = tool_name,
        round_counts = round_counts,
        create_card_params = create_card_params,
        card_deconvolution_params = card_deconvolution_params
      )
    )
    srt@tools[[tool_name]] <- spatial_result_build(
      bundle = srt@tools[[tool_name]],
      method = "CARD",
      result_type = "deconvolution",
      source = c(
        attr(coords, "spatial_source") %||% list(),
        list(transform = attr(coords, "spatial_transform"))
      ),
      provenance = list(producer = "RunCARD", backend_id = "card")
    )
  }

  log_message(
    "{.pkg CARD} proportions stored in metadata columns with prefix {.val {prefix}_prop_}",
    verbose = verbose
  )
  srt
}

card_run_backend <- function(
  st_counts,
  ref_counts,
  ref_meta,
  coords,
  ct_select,
  minCountGene,
  minCountSpot,
  create_card_params,
  card_deconvolution_params
) {
  package <- card_resolve_backend_package()
  create_fun <- get_namespace_fun(package, "createCARDObject")
  deconv_fun <- get_namespace_fun(package, "CARD_deconvolution")
  create_args <- c(
    list(
      sc_count = ref_counts,
      sc_meta = ref_meta,
      spatial_count = st_counts,
      spatial_location = coords,
      ct.varname = ".scop_cell_type",
      ct_varname = ".scop_cell_type",
      sample.varname = ".scop_sample",
      sample_varname = ".scop_sample",
      ct.select = ct_select,
      ct_select = ct_select,
      minCountGene = minCountGene,
      minCountSpot = minCountSpot,
      mincountgene = minCountGene,
      mincountspot = minCountSpot
    ),
    create_card_params
  )
  deconv_formals <- names(formals(deconv_fun))
  if ("CARD_object" %in% deconv_formals) {
    card_obj <- do.call(create_fun, card_match_formals(create_fun, create_args))
    deconv_args <- c(
      list(CARD_object = card_obj),
      card_deconvolution_params
    )
  } else {
    deconv_args <- c(create_args, card_deconvolution_params)
  }
  card_obj <- do.call(deconv_fun, card_match_formals(deconv_fun, deconv_args))
  weights <- card_extract_weights(card_obj)
  list(
    package = package,
    object = card_obj,
    weights = weights
  )
}

card_resolve_backend_package <- function() {
  packages <- c("CARD", "CARDspa")
  available <- vapply(
    packages,
    function(pkg) !is.null(tryCatch(asNamespace(pkg), error = function(e) NULL)),
    logical(1)
  )
  if (!any(available)) {
    invisible(lapply(packages, function(pkg) thisutils::check_r(pkg, verbose = FALSE)))
    available <- vapply(
      packages,
      function(pkg) !is.null(tryCatch(asNamespace(pkg), error = function(e) NULL)),
      logical(1)
    )
  }
  if (!any(available)) {
    log_message(
      "Please install {.pkg CARD} or {.pkg CARDspa} before running {.fn RunCARD}",
      message_type = "error"
    )
  }
  packages[which(available)[1L]]
}

card_reference_metadata <- function(reference, labels, sample_varname = NULL) {
  meta <- reference[[]][names(labels), , drop = FALSE]
  meta[[".scop_cell_type"]] <- as.character(labels)
  if (is.null(sample_varname)) {
    meta[[".scop_sample"]] <- "sample1"
  } else {
    if (length(sample_varname) != 1L || !sample_varname %in% colnames(meta)) {
      log_message(
        "{.arg sample_varname} must be a single reference metadata column",
        message_type = "error"
      )
    }
    meta[[".scop_sample"]] <- as.character(meta[[sample_varname]])
  }
  keep <- !is.na(meta[[".scop_sample"]]) & nzchar(meta[[".scop_sample"]])
  if (!all(keep)) {
    log_message(
      "{.arg sample_varname} contains missing values for CARD reference cells",
      message_type = "error"
    )
  }
  meta
}

card_extract_weights <- function(card_obj) {
  weights <- tryCatch(card_obj$Proportion_CARD, error = function(e) NULL)
  if (is.null(weights)) {
    log_message(
      "{.pkg CARD} did not return {.val Proportion_CARD}",
      message_type = "error"
    )
  }
  weights
}

card_orient_weights <- function(weights, spot_ids) {
  weights <- as.matrix(weights)
  rn_match <- if (is.null(rownames(weights))) 0L else sum(rownames(weights) %in% spot_ids)
  cn_match <- if (is.null(colnames(weights))) 0L else sum(colnames(weights) %in% spot_ids)
  if (cn_match > rn_match) {
    weights <- t(weights)
  }
  if (is.null(rownames(weights))) {
    log_message(
      "{.pkg CARD} proportion matrix must contain spatial spot names",
      message_type = "error"
    )
  }
  spots_use <- spot_ids[spot_ids %in% rownames(weights)]
  if (length(spots_use) == 0L) {
    log_message(
      "{.pkg CARD} proportion matrix could not be matched to spatial spot names",
      message_type = "error"
    )
  }
  weights <- weights[spots_use, , drop = FALSE]
  if (is.null(colnames(weights)) || any(!nzchar(colnames(weights)))) {
    colnames(weights) <- paste0("CellType", seq_len(ncol(weights)))
  }
  colnames(weights) <- make.unique(as.character(colnames(weights)), sep = "_")
  weights
}

card_match_formals <- function(fun, args) {
  fmls <- names(formals(fun))
  if ("..." %in% fmls) {
    return(args)
  }
  args[names(args) %in% fmls]
}

card_validate_param_list <- function(x, arg_name) {
  if (!is.list(x)) {
    log_message(
      "{.arg {arg_name}} must be a list",
      message_type = "error"
    )
  }
  if (length(x) == 0L) {
    return(invisible(TRUE))
  }
  nms <- names(x)
  if (is.null(nms) || any(is.na(nms) | !nzchar(nms))) {
    log_message(
      "{.arg {arg_name}} must contain named arguments only",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

card_assert_scalar_string <- function(x, arg) {
  if (is.null(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    log_message(
      "{.arg {arg}} must be a single non-empty string",
      message_type = "error"
    )
  }
}
