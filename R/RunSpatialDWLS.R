#' @title Run lightweight SpatialDWLS-style deconvolution
#'
#' @description
#' Estimate spot-level cell-type proportions by fitting spatial expression to
#' reference cell-type signatures with non-negative least squares-style
#' coefficients. The output follows SCOP's deconvolution metadata convention so
#' it can be plotted directly with [SpatialSpotPlot()].
#'
#' @md
#' @inheritParams RunCARD
#' @param min_cells Minimum reference cells required per cell type.
#' @param normalize Whether to library-size normalize and `log1p` transform
#' spatial and reference matrices before fitting.
#' @param coordinate_space Coordinate system used for spatial positions.
#'
#' @return A `Seurat` object with `"<prefix>_prop_*"`,
#' `"<prefix>_dominant_type"`, and `"<prefix>_max_prop"` metadata columns.
#' Detailed results are stored in `srt@tools[[tool_name]]`.
#' @export
#'
#' @examples
#' \dontrun{
#' spatial <- RunSpatialDWLS(
#'   spatial,
#'   reference = reference,
#'   reference_label = "celltype",
#'   coord.cols = c("x", "y")
#' )
#' SpatialSpotPlot(spatial, group.by = "SpatialDWLS_dominant_type")
#' SpatialSpotPlot(spatial, group.by = "SpatialDWLS_dominant_type", plot_type = "pie")
#' }
RunSpatialDWLS <- function(
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
  min_cells = 2,
  prefix = "SpatialDWLS",
  tool_name = "SpatialDWLS",
  normalize = TRUE,
  store_results = TRUE,
  verbose = TRUE,
  coordinate_space = c("legacy_display", "raw")
) {
  coordinate_space <- match.arg(coordinate_space)
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  if (!inherits(reference, "Seurat")) {
    log_message("{.arg reference} must be a {.cls Seurat} object", message_type = "error")
  }
  spatial_dwls_assert_string(reference_label, "reference_label")
  spatial_dwls_assert_string(prefix, "prefix")
  spatial_dwls_assert_string(tool_name, "tool_name")
  if (!reference_label %in% colnames(reference@meta.data)) {
    log_message(
      "{.arg reference_label} {.val {reference_label}} is not in reference meta.data",
      message_type = "error"
    )
  }
  min_cells <- as.integer(min_cells[1L])
  if (is.na(min_cells) || min_cells < 1L) {
    log_message("{.arg min_cells} must be a positive integer", message_type = "error")
  }

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  reference_assay <- reference_assay %||% SeuratObject::DefaultAssay(reference)
  st_expr <- GetAssayData5(srt, assay = assay, layer = layer)
  ref_expr <- GetAssayData5(reference, assay = reference_assay, layer = reference_layer)
  shared <- intersect(rownames(st_expr), rownames(ref_expr))
  if (!is.null(features)) {
    shared <- intersect(shared, unique(features))
  }
  if (length(shared) == 0L) {
    log_message("No shared features are available for {.fn RunSpatialDWLS}", message_type = "error")
  }

  labels <- as.character(reference@meta.data[colnames(ref_expr), reference_label, drop = TRUE])
  keep_ref <- !is.na(labels) & nzchar(labels)
  ref_expr <- ref_expr[, keep_ref, drop = FALSE]
  labels <- labels[keep_ref]
  label_counts <- table(labels)
  keep_labels <- names(label_counts)[label_counts >= min_cells]
  if (length(keep_labels) == 0L) {
    log_message(
      "No reference cell type has at least {.val {min_cells}} cell{?s}",
      message_type = "error"
    )
  }
  keep_ref <- labels %in% keep_labels
  ref_expr <- ref_expr[shared, keep_ref, drop = FALSE]
  labels <- labels[keep_ref]
  st_expr <- st_expr[shared, , drop = FALSE]

  coords <- spatial_analysis_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    coordinate_space = coordinate_space
  )$data
  spots <- intersect(colnames(st_expr), rownames(coords))
  if (length(spots) == 0L) {
    log_message("No spatial coordinates match assay spots", message_type = "error")
  }
  st_expr <- st_expr[, spots, drop = FALSE]
  coords <- coords[spots, , drop = FALSE]

  if (isTRUE(normalize)) {
    st_expr <- spatial_dwls_normalize_matrix(st_expr)
    ref_expr <- spatial_dwls_normalize_matrix(ref_expr)
  } else {
    st_expr <- as.matrix(st_expr)
    ref_expr <- as.matrix(ref_expr)
  }

  signatures <- spatial_dwls_signatures(ref_expr, labels)
  keep_features <- Matrix::rowSums(as.matrix(signatures) > 0) > 0 &
    Matrix::rowSums(as.matrix(st_expr) > 0) > 0
  if (!any(keep_features)) {
    log_message("No non-zero shared features remain for {.fn RunSpatialDWLS}", message_type = "error")
  }
  signatures <- signatures[keep_features, , drop = FALSE]
  st_expr <- st_expr[keep_features, , drop = FALSE]

  log_message(
    "Run SpatialDWLS-style fitting with {.val {nrow(st_expr)}} features, {.val {ncol(st_expr)}} spots, and {.val {ncol(signatures)}} cell types",
    message_type = "running",
    verbose = verbose
  )

  weights <- spatial_dwls_fit_weights(signatures = signatures, spatial_expr = st_expr)
  weight_summary <- scop_spatial_finalize_weights(weights = weights, all_spots = colnames(srt))
  weights <- weight_summary$weights
  srt <- scop_spatial_add_deconv_metadata(srt, weights = weights, prefix = prefix, metadata = weight_summary)
  summary <- scop_spatial_weight_summary(weights)

  if (isTRUE(store_results)) {
    srt@tools[[tool_name]] <- list(
      weights = weights,
      signatures = signatures,
      coords = coords,
      features = rownames(st_expr),
      summary = summary,
      parameters = list(
        assay = assay,
        reference_assay = reference_assay,
        layer = layer,
        reference_layer = reference_layer,
        reference_label = reference_label,
        image = image,
        coord.cols = coord.cols,
        coordinate_space = coordinate_space,
        min_cells = min_cells,
        prefix = prefix,
        tool_name = tool_name,
        normalize = normalize
      )
    )
  }

  log_message(
    "SpatialDWLS-style proportions stored in metadata columns with prefix {.val {prefix}_prop_}",
    message_type = "success",
    verbose = verbose
  )
  srt
}

spatial_dwls_normalize_matrix <- function(x) {
  x <- as.matrix(x)
  totals <- colSums(x)
  totals[!is.finite(totals) | totals <= 0] <- 1
  scale_factor <- stats::median(totals[totals > 0], na.rm = TRUE)
  if (!is.finite(scale_factor) || scale_factor <= 0) {
    scale_factor <- 1
  }
  log1p(t(t(x) / totals) * scale_factor)
}

spatial_dwls_signatures <- function(ref_expr, labels) {
  cell_types <- sort(unique(labels))
  sig <- vapply(
    cell_types,
    function(cell_type) {
      Matrix::rowMeans(ref_expr[, labels == cell_type, drop = FALSE])
    },
    numeric(nrow(ref_expr))
  )
  rownames(sig) <- rownames(ref_expr)
  colnames(sig) <- cell_types
  sig
}

spatial_dwls_fit_weights <- function(signatures, spatial_expr) {
  signatures <- as.matrix(signatures)
  spatial_expr <- as.matrix(spatial_expr)
  weights <- matrix(
    0,
    nrow = ncol(spatial_expr),
    ncol = ncol(signatures),
    dimnames = list(colnames(spatial_expr), colnames(signatures))
  )
  qr_sig <- qr(signatures)
  for (spot in colnames(spatial_expr)) {
    coef <- tryCatch(
      qr.coef(qr_sig, spatial_expr[, spot]),
      error = function(e) rep(0, ncol(signatures))
    )
    coef[!is.finite(coef) | coef < 0] <- 0
    weights[spot, ] <- coef
  }
  weights
}

spatial_dwls_assert_string <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    log_message("{.arg {arg}} must be a single non-empty string", message_type = "error")
  }
}
