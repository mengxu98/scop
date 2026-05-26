#' @title Run RCTD spatial deconvolution
#'
#' @description
#' Estimate spot-level cell type proportions from a spatial `Seurat` object
#' using a single-cell `Seurat` reference and `spacexr` RCTD.
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param reference Reference `Seurat` object containing annotated single cells.
#' @param reference_label Metadata column in `reference` with cell type labels.
#' @param reference_assay Assay used in `reference`.
#' @param layer,reference_layer Assay layers used for spatial and reference
#' raw counts.
#' @param features Features used for RCTD. If `NULL`, shared features are used.
#' @param image Name of the Seurat spatial image used to recover coordinates
#' when `coord.cols` are not available.
#' @param coord.cols Metadata coordinate columns used when no image coordinate
#' source is requested or available.
#' @param rctd_mode RCTD mode passed to `spacexr`. `"full"` is the default for
#' Visium spot deconvolution.
#' @param max_cores Number of cores passed to `spacexr`.
#' @param min_cells Minimum number of reference cells required for each cell
#' type. Old `spacexr` RCTD requires at least 25 cells per type.
#' @param prefix Prefix for metadata columns.
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
#' detailed results stored in `srt@tools[["RCTD"]]`.
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' data(panc8_sub)
#'
#' spatial <- RunRCTD(
#'   visium_human_pancreas_sub,
#'   reference = panc8_sub,
#'   reference_label = "celltype",
#'   rctd_mode = "full",
#'   max_cores = 1
#' )
#'
#' rctd_cols <- grep("^RCTD_prop_", colnames(spatial@meta.data), value = TRUE)
#' SpatialSpotPlot(spatial, group.by = rctd_cols[1:4])
#' SpatialSpotPlot(spatial, group.by = rctd_cols, plot_type = "pie")
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
  ...
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
  if (!is.list(create_rctd_params) || !is.list(run_rctd_params)) {
    log_message(
      "{.arg create_rctd_params} and {.arg run_rctd_params} must be lists",
      message_type = "error"
    )
  }
  rctd_validate_named_param_list(create_rctd_params, "create_rctd_params")
  rctd_validate_named_param_list(run_rctd_params, "run_rctd_params")
  rctd_mode <- match.arg(rctd_mode)
  max_cores <- rctd_check_max_cores(max_cores)
  min_cells <- rctd_check_min_cells(min_cells)
  if (length(round_counts) != 1L || !is.logical(round_counts) || is.na(round_counts)) {
    log_message(
      "{.arg round_counts} must be TRUE or FALSE",
      message_type = "error"
    )
  }
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
  if (length(labels) == 0L) {
    log_message(
      "{.arg reference_label} must contain at least one non-missing class",
      message_type = "error"
    )
  }
  labels <- rctd_filter_labels_by_min_cells(
    reference = reference,
    labels = labels,
    min_cells = min_cells,
    verbose = verbose
  )
  dropped_cell_types <- attr(labels, "dropped_cell_types")
  reference <- reference[, names(labels)]
  labels <- factor(as.character(labels), levels = unique(as.character(labels)))
  names(labels) <- colnames(reference)
  label_map <- rctd_backend_label_map(labels)

  features_use <- rctd_common_features(
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
  nonzero_features <- rctd_nonzero_shared_features(
    st_counts = st_counts,
    ref_counts = ref_counts
  )
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

  ref_numi <- Matrix::colSums(ref_counts)
  keep_ref_numi <- is.finite(ref_numi) & ref_numi > 0
  if (!all(keep_ref_numi)) {
    log_message(
      "Drop {.val {sum(!keep_ref_numi)}} reference cells with zero UMI",
      verbose = verbose
    )
    ref_counts <- ref_counts[, keep_ref_numi, drop = FALSE]
    labels <- labels[keep_ref_numi]
    names(labels) <- colnames(ref_counts)
    label_map <- rctd_backend_label_map(labels)
    ref_numi <- ref_numi[keep_ref_numi]
  }
  if (ncol(ref_counts) == 0L) {
    log_message(
      "No reference cells remain after filtering zero-UMI cells",
      message_type = "error"
    )
  }

  st_numi <- Matrix::colSums(st_counts)
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
  coords <- rctd_get_spatial_coords(
    srt = srt,
    spot_ids = colnames(st_counts),
    image = image,
    coord.cols = coord.cols
  )

  extra_run_params <- list(...)
  if (length(extra_run_params) > 0L) {
    rctd_validate_named_param_list(extra_run_params, "...")
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
  weights <- rctd_normalize_weights(weights)
  srt <- rctd_add_metadata(srt, weights = weights, prefix = prefix)

  if (isTRUE(store_results)) {
    srt@tools[["RCTD"]] <- list(
      weights = weights,
      metadata = backend$metadata,
      coords = coords,
      features = nonzero_features,
      cell_types = label_map$cell_types,
      cell_type_map = label_map$map,
      dropped_cell_types = dropped_cell_types,
      backend_api = backend$api,
      parameters = list(
        assay = assay,
        reference_assay = reference_assay,
        layer = layer,
        reference_layer = reference_layer,
        reference_label = reference_label,
        image = image,
        coord.cols = coord.cols,
        rctd_mode = rctd_mode,
        max_cores = max_cores,
        min_cells = min_cells,
        prefix = prefix,
        round_counts = round_counts,
        create_rctd_params = create_rctd_params,
        run_rctd_params = run_rctd_params
      ),
      object = backend$object
    )
  }

  log_message(
    "{.pkg RCTD} proportions stored in metadata columns with prefix {.val {prefix}_prop_}",
    verbose = verbose
  )
  srt
}

rctd_check_max_cores <- function(max_cores) {
  if (
    length(max_cores) != 1L ||
      !is.numeric(max_cores) ||
      is.na(max_cores) ||
      max_cores < 1
  ) {
    log_message(
      "{.arg max_cores} must be a single positive number",
      message_type = "error"
    )
  }
  as.integer(max_cores)
}

rctd_check_min_cells <- function(min_cells) {
  if (
    length(min_cells) != 1L ||
      !is.numeric(min_cells) ||
      is.na(min_cells) ||
      min_cells < 1
  ) {
    log_message(
      "{.arg min_cells} must be a single positive number",
      message_type = "error"
    )
  }
  as.integer(min_cells)
}

rctd_validate_named_param_list <- function(x, arg_name) {
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

rctd_resolve_reference_labels <- function(reference, reference_label) {
  if (
    missing(reference_label) ||
      is.null(reference_label) ||
      length(reference_label) != 1L
  ) {
    log_message(
      "{.arg reference_label} must be a single reference metadata column",
      message_type = "error"
    )
  }
  if (!reference_label %in% colnames(reference[[]])) {
    log_message(
      "{.arg reference_label} {.val {reference_label}} is not present in {.arg reference}",
      message_type = "error"
    )
  }
  reference[[reference_label, drop = TRUE]]
}

rctd_filter_labels_by_min_cells <- function(reference, labels, min_cells, verbose = TRUE) {
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
  names(labels) <- colnames(reference)[keep]
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

rctd_common_features <- function(
  srt,
  reference,
  assay,
  reference_assay,
  features = NULL
) {
  common <- intersect(
    rownames(srt[[assay]]),
    rownames(reference[[reference_assay]])
  )
  if (!is.null(features)) {
    common <- intersect(features, common)
  }
  common
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
    mat <- Matrix::Matrix(as.matrix(mat), sparse = TRUE)
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

rctd_nonzero_shared_features <- function(st_counts, ref_counts) {
  keep <- Matrix::rowSums(st_counts) > 0 & Matrix::rowSums(ref_counts) > 0
  rownames(st_counts)[keep]
}

rctd_get_spatial_coords <- function(
  srt,
  spot_ids,
  image = NULL,
  coord.cols = c("x", "y")
) {
  meta <- srt[[]]
  if (is.null(image) && all(coord.cols %in% colnames(meta))) {
    coords <- meta[spot_ids, coord.cols, drop = FALSE]
    coords <- rctd_validate_coords(coords, spot_ids)
    return(coords)
  }

  image_coords <- rctd_get_image_coords(srt, image = image)
  if (!is.null(image_coords)) {
    common <- intersect(spot_ids, rownames(image_coords))
    if (length(common) == length(spot_ids)) {
      coords <- image_coords[spot_ids, , drop = FALSE]
      coords <- rctd_validate_coords(coords, spot_ids)
      return(coords)
    }
  }

  if (all(coord.cols %in% colnames(meta))) {
    coords <- meta[spot_ids, coord.cols, drop = FALSE]
    coords <- rctd_validate_coords(coords, spot_ids)
    return(coords)
  }

  log_message(
    "Spatial coordinates were not found. Provide a Seurat image or metadata columns {.val {coord.cols}}.",
    message_type = "error"
  )
}

rctd_get_image_coords <- function(srt, image = NULL) {
  images <- tryCatch(SeuratObject::Images(srt), error = function(e) character())
  if (length(images) == 0L) {
    return(NULL)
  }
  image <- image %||% images[1L]
  if (!image %in% images) {
    log_message(
      "{.arg image} {.val {image}} is not present in {.cls Seurat}",
      message_type = "error"
    )
  }
  coords <- tryCatch(
    as.data.frame(SeuratObject::GetTissueCoordinates(srt[[image]])),
    error = function(e) NULL
  )
  if (is.null(coords) || nrow(coords) == 0L) {
    return(NULL)
  }
  if ("cell" %in% colnames(coords)) {
    rownames(coords) <- coords[["cell"]]
  }
  x_col <- rctd_pick_col(coords, c("x", "pxl_col_in_fullres", "imagecol", "col"))
  y_col <- rctd_pick_col(coords, c("y", "pxl_row_in_fullres", "imagerow", "row"))
  if (is.null(x_col) || is.null(y_col)) {
    return(NULL)
  }
  data.frame(
    x = coords[[x_col]],
    y = coords[[y_col]],
    row.names = rownames(coords),
    stringsAsFactors = FALSE
  )
}

rctd_pick_col <- function(x, candidates) {
  nm <- colnames(x)
  hit <- candidates[tolower(candidates) %in% tolower(nm)][1L]
  if (is.na(hit)) {
    return(NULL)
  }
  nm[match(tolower(hit), tolower(nm))]
}

rctd_validate_coords <- function(coords, spot_ids) {
  coords <- as.data.frame(coords, check.names = FALSE)
  if (ncol(coords) < 2L) {
    log_message(
      "Spatial coordinates must contain at least two columns",
      message_type = "error"
    )
  }
  coords <- coords[, seq_len(2L), drop = FALSE]
  colnames(coords) <- c("x", "y")
  coords[] <- lapply(coords, function(x) suppressWarnings(as.numeric(x)))
  if (!identical(rownames(coords), spot_ids)) {
    coords <- coords[spot_ids, , drop = FALSE]
  }
  keep <- stats::complete.cases(coords) &
    is.finite(coords$x) &
    is.finite(coords$y)
  if (!all(keep)) {
    log_message(
      "Spatial coordinates contain missing or non-finite values for {.val {sum(!keep)}} spots",
      message_type = "error"
    )
  }
  coords
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
  check_r("spacexr", verbose = FALSE)
  ns <- getNamespace("spacexr")
  has_new_api <- exists("createRctd", envir = ns, inherits = FALSE) &&
    exists("runRctd", envir = ns, inherits = FALSE)
  has_old_api <- all(vapply(
    c("SpatialRNA", "Reference", "create.RCTD", "run.RCTD"),
    exists,
    logical(1),
    envir = ns,
    inherits = FALSE
  ))

  if (isTRUE(has_new_api)) {
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
  if (isTRUE(has_old_api)) {
    return(rctd_run_spacexr_old(
      ns = ns,
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
  check_r(c("SpatialExperiment", "SummarizedExperiment", "S4Vectors"), verbose = FALSE)
  spatial_coldata <- S4Vectors::DataFrame(nUMI = st_numi)
  rownames(spatial_coldata) <- colnames(st_counts)
  reference_coldata <- S4Vectors::DataFrame(cell_type = ref_labels, nUMI = ref_numi)
  rownames(reference_coldata) <- colnames(ref_counts)
  spatial_spe <- SpatialExperiment::SpatialExperiment(
    assay = st_counts,
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

rctd_run_spacexr_old <- function(
  ns,
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
  spatial_rna <- get("SpatialRNA", envir = ns)
  reference_fun <- get("Reference", envir = ns)
  create_rctd <- get("create.RCTD", envir = ns)
  run_rctd <- get("run.RCTD", envir = ns)
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
  if (exists("normalize_weights", envir = ns, inherits = FALSE)) {
    weights <- get("normalize_weights", envir = ns)(weights)
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

rctd_normalize_weights <- function(weights) {
  weights <- as.matrix(weights)
  weights[!is.finite(weights) | weights < 0] <- 0
  totals <- rowSums(weights)
  keep <- is.finite(totals) & totals > 0
  weights[keep, ] <- sweep(weights[keep, , drop = FALSE], 1, totals[keep], "/")
  weights[!keep, ] <- 0
  weights
}

rctd_add_metadata <- function(srt, weights, prefix = "RCTD") {
  all_spots <- colnames(srt)
  full_weights <- matrix(
    NA_real_,
    nrow = length(all_spots),
    ncol = ncol(weights),
    dimnames = list(all_spots, colnames(weights))
  )
  common <- intersect(all_spots, rownames(weights))
  full_weights[common, ] <- weights[common, , drop = FALSE]
  meta <- as.data.frame(full_weights, check.names = FALSE)
  prop_cols <- paste0(
    prefix,
    "_prop_",
    make.unique(make.names(colnames(full_weights)), sep = "_")
  )
  colnames(meta) <- prop_cols

  dominant <- rep(NA_character_, nrow(full_weights))
  max_prop <- rep(NA_real_, nrow(full_weights))
  row_has_value <- rowSums(!is.na(full_weights)) > 0
  if (any(row_has_value)) {
    mat <- full_weights[row_has_value, , drop = FALSE]
    mat[is.na(mat)] <- 0
    max_idx <- max.col(mat, ties.method = "first")
    row_max <- mat[cbind(seq_len(nrow(mat)), max_idx)]
    dominant[row_has_value] <- ifelse(
      row_max > 0,
      colnames(full_weights)[max_idx],
      NA_character_
    )
    max_prop[row_has_value] <- row_max
  }
  meta[[paste0(prefix, "_dominant_type")]] <- dominant
  meta[[paste0(prefix, "_max_prop")]] <- max_prop
  Seurat::AddMetaData(srt, metadata = meta)
}
