#' @title Run RCTD spatial deconvolution
#'
#' @description
#' Estimate spot-level cell type proportions from a spatial `Seurat` object
#' using a single-cell `Seurat` reference and `spacexr` RCTD.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt Spatial `Seurat` object used as the RCTD query.
#' @param reference Reference `Seurat` object containing annotated single cells.
#' @param reference_label Metadata column in `reference` with cell type labels.
#' @param assay Assay used in `srt`. If `NULL`, the default assay is used.
#' @param reference_assay Assay used in `reference`.
#' @param layer,reference_layer Assay layers used for spatial and reference
#' raw counts.
#' @param features Features used for RCTD. If `NULL`, shared features are used.
#' @param image Name of the Seurat spatial image used to recover coordinates
#' when `coord.cols` are not available.
#' @param coord.cols Metadata coordinate columns used when no image coordinate
#' source is requested or available.
#' @param coordinate_space Coordinate space used for distance-sensitive input.
#' The default is raw acquisition coordinates, so backend distances use raw
#' coordinate units. Use `"legacy_display"` explicitly to reproduce the
#' display-scaled coordinates used before scop 0.9.0.
#' @param rctd_mode RCTD mode passed to `spacexr`. `"full"` is the default for
#' Visium spot deconvolution.
#' @param max_cores Number of cores passed to `spacexr`.
#' @param min_cells Minimum number of reference cells required for each cell
#' type. Old `spacexr` RCTD requires at least 25 cells per type.
#' @param prefix Prefix for metadata columns.
#' @param tool_name Name used to store the plain result bundle in `srt@tools`.
#' @param store_results Whether to store detailed RCTD results in `srt@tools`.
#' @param round_counts Whether to round non-integer counts to the nearest
#' integer before passing data to `spacexr`. RCTD requires integer count
#' matrices; this defaults to `TRUE` so bundled example data with scaled
#' non-integer reference counts can run directly.
#' @param create_rctd_params Additional parameters passed to
#' `spacexr::createRctd()` or `spacexr::create.RCTD()`.
#' @param run_rctd_params Additional parameters passed to
#' `spacexr::runRctd()` or `spacexr::run.RCTD()`.
#' @param ... Additional parameters passed to the RCTD run step.
#'
#' @return A `Seurat` object with RCTD proportion columns in metadata and
#' dominant cell type summaries. When `store_results = TRUE`, detailed results
#' are also stored in `srt@tools[[tool_name]]`.
#' @export
#'
#' @examples
#' \dontrun{
#' data(visium_human_pancreas_sub)
#' data(panc8_sub)
#' keep_spots <- unique(round(seq(
#'   1,
#'   ncol(visium_human_pancreas_sub),
#'   length.out = 120
#' )))
#' spatial <- visium_human_pancreas_sub[, keep_spots]
#' reference <- panc8_sub[, panc8_sub@meta.data[["celltype"]] %in%
#'   c("ductal", "alpha", "beta")]
#' reference <- Seurat::FindVariableFeatures(reference, nfeatures = 300, verbose = FALSE)
#' features_use <- head(intersect(
#'   SeuratObject::VariableFeatures(reference),
#'   rownames(spatial)
#' ), 300)
#' spatial <- RunRCTD(
#'   srt = spatial,
#'   reference = reference,
#'   reference_label = "celltype",
#'   assay = "Spatial",
#'   reference_assay = "RNA",
#'   layer = "counts",
#'   reference_layer = "counts",
#'   features = features_use,
#'   coord.cols = c("x", "y"),
#'   rctd_mode = "full",
#'   max_cores = 1,
#'   min_cells = 25,
#'   prefix = "RCTD",
#'   verbose = FALSE
#' )
#' SpatialDeconvolutionPlot(
#'   spatial,
#'   tool_name = "RCTD",
#'   plot_type = "dominant",
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#' }
RunRCTD <- function(
  srt,
  reference,
  reference_label = "celltype",
  assay = NULL,
  reference_assay = NULL,
  layer = "counts",
  reference_layer = "counts",
  features = NULL,
  image = NULL,
  coord.cols = c("x", "y"),
  rctd_mode = c("full", "multi", "doublet"),
  max_cores = 1,
  min_cells = 25,
  prefix = "RCTD",
  store_results = TRUE,
  round_counts = TRUE,
  create_rctd_params = list(),
  run_rctd_params = list(),
  verbose = TRUE,
  ...,
  coordinate_space = c("raw", "legacy_display"),
  tool_name = "RCTD"
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
  validate_scalar_string(prefix, "prefix")
  validate_scalar_string(tool_name, "tool_name")
  validate_scalar_flag(store_results, "store_results")
  if (!is.list(create_rctd_params) || !is.list(run_rctd_params)) {
    log_message(
      "{.arg create_rctd_params} and {.arg run_rctd_params} must be lists",
      message_type = "error"
    )
  }
  validate_named_param_list(create_rctd_params, "create_rctd_params")
  validate_named_param_list(run_rctd_params, "run_rctd_params")
  rctd_mode <- match.arg(rctd_mode)
  coordinate_space <- match.arg(coordinate_space)
  max_cores <- rctd_check_single_positive(max_cores, "max_cores")
  min_cells <- rctd_check_single_positive(min_cells, "min_cells")
  validate_scalar_flag(round_counts, "round_counts")
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  reference_assay <- reference_assay %||% SeuratObject::DefaultAssay(reference)

  labels <- resolve_reference_labels(reference, reference_label)
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
  if (length(labels) == 0L) {
    log_message(
      "{.arg reference_label} must contain at least one non-missing class",
      message_type = "error"
    )
  }
  labels <- rctd_filter_labels_by_min_cells(
    labels = labels,
    min_cells = min_cells,
    verbose = verbose
  )
  dropped_cell_types <- attr(labels, "dropped_cell_types")
  reference <- reference[, names(labels)]
  labels <- factor(as.character(labels), levels = unique(as.character(labels)))
  names(labels) <- colnames(reference)
  label_map <- rctd_backend_label_map(labels)

  features_use <- resolve_common_features(
    srt = srt,
    reference = reference,
    assay = assay,
    reference_assay = reference_assay,
    features = features
  )
  if (length(features_use) == 0L) {
    log_message(
      "No shared features are available for {.fn RunRCTD}",
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
  count_quality <- rctd_sparse_quality_cpp(
    st_counts = st_counts,
    ref_counts = ref_counts
  )
  nonzero_features <- rownames(st_counts)[count_quality$keep_features]
  if (length(nonzero_features) == 0L) {
    log_message(
      "No shared features have non-zero counts in both spatial and reference data",
      message_type = "error"
    )
  }
  if (length(nonzero_features) < length(features_use)) {
    log_message(
      "Use {.val {length(nonzero_features)}} shared non-zero features for RCTD",
      verbose = verbose
    )
  }
  st_counts <- st_counts[nonzero_features, , drop = FALSE]
  ref_counts <- ref_counts[nonzero_features, , drop = FALSE]

  ref_numi <- count_quality$ref_numi
  names(ref_numi) <- colnames(ref_counts)
  keep_ref_numi <- is.finite(ref_numi) & ref_numi > 0
  if (!all(keep_ref_numi)) {
    log_message(
      "Drop {.val {sum(!keep_ref_numi)}} reference cells with zero UMI",
      verbose = verbose
    )
    ref_counts <- ref_counts[, keep_ref_numi, drop = FALSE]
    labels <- labels[keep_ref_numi]
    names(labels) <- colnames(ref_counts)
    labels <- rctd_filter_labels_by_min_cells(
      labels = labels,
      min_cells = min_cells,
      verbose = verbose
    )
    dropped_cell_types <- rbind(
      dropped_cell_types,
      attr(labels, "dropped_cell_types")
    )
    ref_counts <- ref_counts[, names(labels), drop = FALSE]
    ref_numi <- ref_numi[names(labels)]
    labels <- factor(as.character(labels), levels = unique(as.character(labels)))
    names(labels) <- colnames(ref_counts)
    label_map <- rctd_backend_label_map(labels)
  }
  if (ncol(ref_counts) == 0L) {
    log_message(
      "No reference cells remain after filtering zero-UMI cells",
      message_type = "error"
    )
  }

  st_numi <- count_quality$st_numi
  names(st_numi) <- colnames(st_counts)
  keep_spots <- is.finite(st_numi) & st_numi > 0
  if (!all(keep_spots)) {
    log_message(
      "Drop {.val {sum(!keep_spots)}} spatial spots with zero UMI for RCTD",
      verbose = verbose
    )
    st_counts <- st_counts[, keep_spots, drop = FALSE]
    st_numi <- st_numi[keep_spots]
  }
  if (ncol(st_counts) == 0L) {
    log_message(
      "No spatial spots remain after filtering zero-UMI spots",
      message_type = "error"
    )
  }
  coords <- resolve_spatial_spot_coords(
    srt = srt,
    spot_ids = colnames(st_counts),
    image = image,
    coord.cols = coord.cols,
    coordinate_space = coordinate_space
  )

  extra_run_params <- list(...)
  if (length(extra_run_params) > 0L) {
    validate_named_param_list(extra_run_params, "...")
    run_rctd_params <- c(run_rctd_params, extra_run_params)
  }

  log_message(
    "Run {.pkg spacexr} RCTD with {.val {nrow(st_counts)}} features, {.val {ncol(st_counts)}} spatial spots, and {.val {ncol(ref_counts)}} reference cells",
    verbose = verbose
  )
  backend <- rctd_run_spacexr(
    st_counts = st_counts,
    coords = coords,
    st_numi = st_numi,
    ref_counts = ref_counts,
    ref_labels = label_map$labels,
    ref_numi = ref_numi,
    rctd_mode = rctd_mode,
    max_cores = max_cores,
    create_rctd_params = create_rctd_params,
    run_rctd_params = run_rctd_params
  )

  weights <- rctd_orient_weights(
    weights = backend$weights,
    spot_ids = colnames(st_counts),
    label_map = label_map
  )
  weight_summary <- spatial_finalize_weights(
    weights = weights,
    all_spots = colnames(srt)
  )
  weights <- weight_summary$weights
  srt <- spatial_add_deconv_metadata(
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
      metadata = backend$metadata,
      coords = coords,
      features = nonzero_features,
      cell_types = label_map$cell_types,
      cell_type_map = label_map$map,
      dropped_cell_types = dropped_cell_types,
      summary = spatial_weight_summary(weights),
      backend_api = backend$api,
      parameters = list(
        assay = assay,
        reference_assay = reference_assay,
        layer = layer,
        reference_layer = reference_layer,
        reference_label = reference_label,
        image = image,
        coord.cols = coord.cols,
        coordinate_space = coordinate_space,
        coordinate_dependent = TRUE,
        rctd_mode = rctd_mode,
        max_cores = max_cores,
        min_cells = min_cells,
        prefix = prefix,
        tool_name = tool_name,
        round_counts = round_counts,
        create_rctd_params = create_rctd_params,
        run_rctd_params = run_rctd_params
      ),
      object = backend$object
    )
    srt@tools[[tool_name]] <- spatial_tag_coordinate_contract(srt@tools[[tool_name]])
  }

  log_message(
    "{.pkg RCTD} proportions stored in metadata columns with prefix {.val {prefix}_prop_}",
    verbose = verbose
  )
  srt
}


rctd_check_single_positive <- function(x, arg_name) {
  if (
    length(x) != 1L ||
      !is.numeric(x) ||
      is.na(x) ||
      x < 1
  ) {
    log_message(
      "{.arg {arg_name}} must be a single positive number",
      message_type = "error"
    )
  }
  as.integer(x)
}


rctd_filter_labels_by_min_cells <- function(labels, min_cells, verbose = TRUE) {
  label_counts <- table(labels)
  keep_types <- names(label_counts)[label_counts >= min_cells]
  drop_types <- names(label_counts)[label_counts < min_cells]
  dropped <- data.frame(
    cell_type = drop_types,
    n_cells = as.integer(label_counts[drop_types]),
    stringsAsFactors = FALSE
  )
  if (length(drop_types) > 0L) {
    drop_summary <- paste0(drop_types, "=", unname(label_counts[drop_types]), collapse = ", ")
    log_message(
      "Drop reference cell types with fewer than {.val {min_cells}} cells: {.val {drop_summary}}",
      verbose = verbose
    )
  }
  keep <- as.character(labels) %in% keep_types
  labels <- labels[keep]
  attr(labels, "dropped_cell_types") <- dropped
  if (length(unique(as.character(labels))) == 0L) {
    log_message(
      "No reference cell types remain after filtering by {.arg min_cells}",
      message_type = "error"
    )
  }
  labels
}

rctd_backend_label_map <- function(labels) {
  labels_chr <- as.character(labels)
  cell_types <- unique(labels_chr)
  backend <- make.names(cell_types)
  backend <- gsub("/", "_", backend, fixed = TRUE)
  backend[!nzchar(backend)] <- paste0("CellType", which(!nzchar(backend)))
  backend <- make.unique(backend, sep = "_")
  names(backend) <- cell_types
  backend_labels <- factor(
    unname(backend[labels_chr]),
    levels = unname(backend)
  )
  names(backend_labels) <- names(labels)
  list(
    labels = backend_labels,
    cell_types = cell_types,
    backend = unname(backend),
    map = data.frame(
      cell_type = cell_types,
      backend_label = unname(backend),
      stringsAsFactors = FALSE
    )
  )
}

rctd_get_count_matrix <- function(
  srt,
  assay,
  layer,
  features,
  data_label = "Input",
  round_counts = TRUE,
  verbose = TRUE
) {
  mat <- GetAssayData5(srt, assay = assay, layer = layer)
  mat <- mat[features, , drop = FALSE]
  if (!inherits(mat, "Matrix")) {
    mat <- Matrix::Matrix(mat, sparse = TRUE)
  }
  if (!inherits(mat, "dgCMatrix")) {
    mat <- methods::as(mat, "dgCMatrix")
  }
  mat@x[!is.finite(mat@x) | mat@x < 0] <- 0
  non_integer <- abs(mat@x - round(mat@x)) > sqrt(.Machine$double.eps)
  if (any(non_integer)) {
    if (!isTRUE(round_counts)) {
      log_message(
        "{.val {data_label}} counts contain non-integer values, but {.arg round_counts = FALSE}. RCTD requires integer counts.",
        message_type = "error"
      )
    }
    log_message(
      "{.val {data_label}} counts contain {.val {sum(non_integer)}} non-integer values; round to nearest integers for RCTD",
      verbose = verbose
    )
    mat@x <- round(mat@x)
    mat@x[mat@x < 0] <- 0
    mat <- methods::as(mat, "dgCMatrix")
    mat <- Matrix::drop0(mat)
  }
  mat
}


rctd_spacexr_exports <- function() {
  getNamespaceExports("spacexr")
}

rctd_run_spacexr <- function(
  st_counts,
  coords,
  st_numi,
  ref_counts,
  ref_labels,
  ref_numi,
  rctd_mode,
  max_cores,
  create_rctd_params,
  run_rctd_params
) {
  rctd_require_namespaces("spacexr")
  exports <- rctd_spacexr_exports()
  new_api <- c("createRctd", "runRctd")
  legacy_api <- c("SpatialRNA", "Reference", "create.RCTD", "run.RCTD")
  if (all(new_api %in% exports)) {
    return(rctd_run_spacexr_new(
      st_counts = st_counts,
      coords = coords,
      st_numi = st_numi,
      ref_counts = ref_counts,
      ref_labels = ref_labels,
      ref_numi = ref_numi,
      rctd_mode = rctd_mode,
      max_cores = max_cores,
      create_rctd_params = create_rctd_params,
      run_rctd_params = run_rctd_params
    ))
  }
  if (all(legacy_api %in% exports)) {
    return(rctd_run_spacexr_old(
      st_counts = st_counts,
      coords = coords,
      st_numi = st_numi,
      ref_counts = ref_counts,
      ref_labels = ref_labels,
      ref_numi = ref_numi,
      rctd_mode = rctd_mode,
      max_cores = max_cores,
      create_rctd_params = create_rctd_params,
      run_rctd_params = run_rctd_params
    ))
  }

  log_message(
    "{.pkg spacexr} does not expose a supported RCTD API",
    message_type = "error"
  )
}

rctd_run_spacexr_new <- function(
  st_counts,
  coords,
  st_numi,
  ref_counts,
  ref_labels,
  ref_numi,
  rctd_mode,
  max_cores,
  create_rctd_params,
  run_rctd_params
) {
  rctd_require_namespaces(c("SpatialExperiment", "SummarizedExperiment", "S4Vectors"))
  spatial_coldata <- S4Vectors::DataFrame(nUMI = st_numi)
  rownames(spatial_coldata) <- colnames(st_counts)
  reference_coldata <- S4Vectors::DataFrame(cell_type = ref_labels, nUMI = ref_numi)
  rownames(reference_coldata) <- colnames(ref_counts)
  spatial_experiment <- get_namespace_fun("SpatialExperiment", "SpatialExperiment")
  spatial_spe <- spatial_experiment(
    assays = list(counts = st_counts),
    colData = spatial_coldata,
    spatialCoords = as.matrix(coords)
  )
  reference_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = ref_counts),
    colData = reference_coldata
  )
  create_rctd <- get_namespace_fun("spacexr", "createRctd")
  run_rctd <- get_namespace_fun("spacexr", "runRctd")
  create_args <- c(
    list(spatial_spe, reference_se),
    utils::modifyList(list(cell_type_col = "cell_type"), create_rctd_params)
  )
  rctd_data <- do.call(create_rctd, create_args)
  run_args <- c(
    list(rctd_data),
    utils::modifyList(
      list(rctd_mode = rctd_mode, max_cores = max_cores),
      run_rctd_params
    )
  )
  result <- do.call(run_rctd, run_args)
  weights <- SummarizedExperiment::assay(result, "weights")
  metadata <- as.data.frame(SummarizedExperiment::colData(result))
  list(
    api = "createRctd/runRctd",
    weights = weights,
    metadata = metadata,
    object = result
  )
}

rctd_require_namespaces <- function(pkgs) {
  install_specs <- pkgs
  install_specs[install_specs == "spacexr"] <- "dmcable/spacexr"
  available <- unlist(
    check_r(install_specs, verbose = FALSE),
    use.names = TRUE
  )
  if (
    !is.logical(available) || length(available) != length(pkgs) ||
      is.null(names(available)) || anyNA(available) ||
      !setequal(names(available), pkgs)
  ) {
    log_message(
      "Package availability checks returned an invalid result for {.fn RunRCTD}",
      message_type = "error"
    )
  }
  available <- available[match(pkgs, names(available))]
  if (!all(available)) {
    log_message(
      "Please install required package(s) before running {.fn RunRCTD}: {.val {paste(pkgs[!available], collapse = ', ')}}",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

rctd_run_spacexr_old <- function(
  st_counts,
  coords,
  st_numi,
  ref_counts,
  ref_labels,
  ref_numi,
  rctd_mode,
  max_cores,
  create_rctd_params,
  run_rctd_params
) {
  spatial_rna <- get_namespace_fun("spacexr", "SpatialRNA")
  reference_fun <- get_namespace_fun("spacexr", "Reference")
  create_rctd <- get_namespace_fun("spacexr", "create.RCTD")
  run_rctd <- get_namespace_fun("spacexr", "run.RCTD")
  puck <- spatial_rna(coords, st_counts, st_numi)
  reference_obj <- reference_fun(ref_counts, ref_labels, ref_numi)
  create_args <- c(
    list(puck, reference_obj),
    utils::modifyList(list(max_cores = max_cores), create_rctd_params)
  )
  rctd_obj <- do.call(create_rctd, create_args)
  run_args <- c(
    list(rctd_obj),
    utils::modifyList(list(doublet_mode = rctd_mode), run_rctd_params)
  )
  result <- do.call(run_rctd, run_args)
  weights <- result@results$weights
  if (is.null(weights) && is.list(result@results)) {
    first <- result@results[[1L]]
    if (is.list(first) && !is.null(first$all_weights)) {
      cell_types <- result@cell_type_info$renorm[[2]]
      weights <- do.call(rbind, lapply(result@results, function(x) x$all_weights))
      if (is.null(rownames(weights))) {
        rownames(weights) <- colnames(result@spatialRNA@counts)
      }
      colnames(weights) <- cell_types
      weights <- methods::as(Matrix::Matrix(weights, sparse = TRUE), "CsparseMatrix")
    }
  }
  if (is.null(weights)) {
    log_message(
      "spacexr returned no weight matrix in {.val {rctd_mode}} mode; check the installed spacexr version",
      message_type = "error"
    )
  }
  if ("normalize_weights" %in% getNamespaceExports("spacexr")) {
    normalize_weights <- get_namespace_fun("spacexr", "normalize_weights")
    weights <- normalize_weights(weights)
  }
  metadata <- tryCatch(
    as.data.frame(result@results$results_df),
    error = function(e) data.frame(row.names = colnames(st_counts))
  )
  list(
    api = "SpatialRNA/Reference/create.RCTD/run.RCTD",
    weights = weights,
    metadata = metadata,
    object = result
  )
}

rctd_orient_weights <- function(weights, spot_ids, label_map) {
  mat <- as.matrix(weights)
  rn_match <- if (is.null(rownames(mat))) 0L else sum(rownames(mat) %in% spot_ids)
  cn_match <- if (is.null(colnames(mat))) 0L else sum(colnames(mat) %in% spot_ids)
  if (cn_match > rn_match) {
    mat <- t(mat)
  }
  if (is.null(rownames(mat))) {
    log_message(
      "RCTD weights must contain spatial spot names",
      message_type = "error"
    )
  }
  spots_use <- spot_ids[spot_ids %in% rownames(mat)]
  if (length(spots_use) == 0L) {
    log_message(
      "RCTD weights could not be matched to spatial spot names",
      message_type = "error"
    )
  }
  mat <- mat[spots_use, , drop = FALSE]
  if (is.null(colnames(mat)) || any(!nzchar(colnames(mat)))) {
    colnames(mat) <- label_map$backend[seq_len(ncol(mat))]
  }
  display <- label_map$cell_types[match(colnames(mat), label_map$backend)]
  display[is.na(display)] <- colnames(mat)[is.na(display)]
  colnames(mat) <- make.unique(display, sep = "_")
  mat
}
