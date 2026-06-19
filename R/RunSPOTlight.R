#' @title Run SPOTlight spatial deconvolution
#'
#' @description
#' Estimate spot-level cell type proportions from a spatial `Seurat` object
#' using a single-cell `Seurat` reference and the optional `SPOTlight` package.
#'
#' @md
#' @inheritParams RunRCTD
#' @param mgs Optional marker-gene table passed to `SPOTlight`. It must contain
#' columns named by `gene_id`, `group_id`, and `weight_id`. If `NULL`, a simple
#' group-vs-rest marker table is generated from the reference expression matrix.
#' @param marker_top_n Number of automatically generated marker genes retained
#' per reference cell type.
#' @param marker_min_logfc Minimum group-vs-rest log2 fold-change used when
#' generating markers automatically.
#' @param gene_id,group_id,weight_id Column names in `mgs` for gene IDs, cell
#' type labels, and marker weights.
#' @param min_prop Minimum cell-type proportion passed to `SPOTlight`.
#' @param scale Whether `SPOTlight` scales expression internally.
#' @param tool_name Name used to store detailed results in `srt@tools`.
#' @param ... Additional parameters passed to `SPOTlight::SPOTlight()`.
#'
#' @return A `Seurat` object with `SPOTlight` proportion columns in metadata and
#' dominant cell type summaries. When `store_results = TRUE`, detailed results
#' are stored in `srt@tools[[tool_name]]`.
#' @export
#'
#' @examples
#' \dontrun{
#' data(visium_human_pancreas_sub)
#' data(panc8_sub)
#'
#' spatial <- RunSPOTlight(
#'   srt = visium_human_pancreas_sub,
#'   reference = panc8_sub,
#'   reference_label = "celltype",
#'   assay = "Spatial",
#'   reference_assay = "RNA"
#' )
#'
#' SpatialSpotPlot(spatial, group.by = "SPOTlight_dominant_type")
#'
#' spotlight_cols <- grep(
#'   "^SPOTlight_prop_",
#'   colnames(spatial@meta.data),
#'   value = TRUE
#' )
#' SpatialSpotPlot(spatial, group.by = spotlight_cols[1:4])
#'
#' SpatialSpotPlot(
#'   spatial,
#'   group.by = "SPOTlight_dominant_type",
#'   plot_type = "pie"
#' )
#' }
RunSPOTlight <- function(
  srt,
  reference,
  reference_label = "celltype",
  assay = NULL,
  reference_assay = NULL,
  layer = "counts",
  reference_layer = "counts",
  features = NULL,
  mgs = NULL,
  marker_top_n = 100,
  marker_min_logfc = 0,
  gene_id = "gene",
  group_id = "cluster",
  weight_id = "weight",
  min_prop = 0.01,
  scale = TRUE,
  prefix = "SPOTlight",
  tool_name = "SPOTlight",
  store_results = TRUE,
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
  spotlight_assert_scalar_string(prefix, "prefix")
  spotlight_assert_scalar_string(tool_name, "tool_name")

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  reference_assay <- reference_assay %||% SeuratObject::DefaultAssay(reference)

  labels <- spotlight_reference_labels(reference, reference_label)
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
  if (length(levels(labels)) < 1L) {
    log_message(
      "{.arg reference_label} must contain at least one non-missing class",
      message_type = "error"
    )
  }

  features_use <- spotlight_common_features(
    srt = srt,
    reference = reference,
    assay = assay,
    reference_assay = reference_assay,
    features = features
  )
  if (length(features_use) == 0L) {
    log_message(
      "No shared features are available for {.fn RunSPOTlight}",
      message_type = "error"
    )
  }

  st_counts <- spotlight_get_matrix(
    srt,
    assay = assay,
    layer = layer,
    features = features_use,
    data_label = "Spatial"
  )
  ref_counts <- spotlight_get_matrix(
    reference,
    assay = reference_assay,
    layer = reference_layer,
    features = features_use,
    data_label = "Reference"
  )
  keep_features <- Matrix::rowSums(st_counts) > 0 & Matrix::rowSums(ref_counts) > 0
  if (!any(keep_features)) {
    log_message(
      "No shared features have non-zero expression in both spatial and reference data",
      message_type = "error"
    )
  }
  if (sum(keep_features) < length(features_use)) {
    log_message(
      "Use {.val {sum(keep_features)}} shared non-zero features for SPOTlight",
      verbose = verbose
    )
  }
  features_use <- rownames(st_counts)[keep_features]
  st_counts <- st_counts[features_use, , drop = FALSE]
  ref_counts <- ref_counts[features_use, , drop = FALSE]

  mgs <- spotlight_prepare_mgs(
    mgs = mgs,
    ref_counts = ref_counts,
    labels = labels,
    gene_id = gene_id,
    group_id = group_id,
    weight_id = weight_id,
    marker_top_n = marker_top_n,
    marker_min_logfc = marker_min_logfc
  )

  check_r("SPOTlight", verbose = FALSE)

  log_message(
    "Run {.pkg SPOTlight} with {.val {nrow(ref_counts)}} features, {.val {ncol(st_counts)}} spatial spots, and {.val {ncol(ref_counts)}} reference cells",
    verbose = verbose
  )
  result <- spotlight_run_backend(
    ref_counts = ref_counts,
    st_counts = st_counts,
    labels = labels,
    mgs = mgs,
    gene_id = gene_id,
    group_id = group_id,
    weight_id = weight_id,
    min_prop = min_prop,
    scale = scale,
    verbose = verbose,
    ...
  )

  weights <- spotlight_extract_weights(result, spot_ids = colnames(st_counts))
  weight_summary <- rctd_finalize_weights_cpp(
    weights = weights,
    all_spots = colnames(srt)
  )
  weights <- weight_summary$weights
  srt <- rctd_add_metadata(
    srt,
    weights = weights,
    prefix = prefix,
    metadata = weight_summary
  )

  if (isTRUE(store_results)) {
    srt@tools[[tool_name]] <- list(
      weights = weights,
      residual = if (is.list(result)) result$res_ss %||% NULL else NULL,
      marker_genes = mgs,
      features = features_use,
      result = result,
      parameters = list(
        assay = assay,
        reference_assay = reference_assay,
        layer = layer,
        reference_layer = reference_layer,
        reference_label = reference_label,
        marker_top_n = marker_top_n,
        marker_min_logfc = marker_min_logfc,
        gene_id = gene_id,
        group_id = group_id,
        weight_id = weight_id,
        min_prop = min_prop,
        scale = scale,
        prefix = prefix,
        tool_name = tool_name
      )
    )
  }

  log_message(
    "{.pkg SPOTlight} proportions stored in metadata columns with prefix {.val {prefix}_prop_}",
    verbose = verbose
  )
  srt
}

spotlight_reference_labels <- function(reference, reference_label) {
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

spotlight_common_features <- function(
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

spotlight_get_matrix <- function(
  srt,
  assay,
  layer,
  features,
  data_label = "Input"
) {
  mat <- GetAssayData5(srt, assay = assay, layer = layer)
  mat <- mat[features, , drop = FALSE]
  if (!inherits(mat, "Matrix")) {
    mat <- Matrix::Matrix(as.matrix(mat), sparse = TRUE)
  }
  if (!inherits(mat, "dgCMatrix")) {
    mat <- methods::as(mat, "dgCMatrix")
  }
  mat@x[!is.finite(mat@x) | mat@x < 0] <- 0
  if (nrow(mat) == 0L || ncol(mat) == 0L) {
    log_message(
      "{.val {data_label}} matrix is empty after feature filtering",
      message_type = "error"
    )
  }
  Matrix::drop0(mat)
}

spotlight_prepare_mgs <- function(
  mgs,
  ref_counts,
  labels,
  gene_id,
  group_id,
  weight_id,
  marker_top_n,
  marker_min_logfc
) {
  if (is.null(mgs)) {
    return(spotlight_auto_mgs(
      ref_counts = ref_counts,
      labels = labels,
      marker_top_n = marker_top_n,
      marker_min_logfc = marker_min_logfc,
      gene_id = gene_id,
      group_id = group_id,
      weight_id = weight_id
    ))
  }
  mgs <- as.data.frame(mgs, check.names = FALSE)
  required <- c(gene_id, group_id, weight_id)
  missing <- setdiff(required, colnames(mgs))
  if (length(missing) > 0L) {
    log_message(
      "{.arg mgs} is missing required column{?s}: {.val {missing}}",
      message_type = "error"
    )
  }
  mgs[[gene_id]] <- as.character(mgs[[gene_id]])
  mgs[[group_id]] <- as.character(mgs[[group_id]])
  mgs[[weight_id]] <- suppressWarnings(as.numeric(mgs[[weight_id]]))
  keep <- mgs[[gene_id]] %in% rownames(ref_counts) &
    mgs[[group_id]] %in% levels(labels) &
    is.finite(mgs[[weight_id]]) &
    mgs[[weight_id]] > 0
  mgs <- mgs[keep, , drop = FALSE]
  spotlight_validate_mgs_groups(mgs, labels, gene_id, group_id)
  mgs
}

spotlight_auto_mgs <- function(
  ref_counts,
  labels,
  marker_top_n = 100,
  marker_min_logfc = 0,
  gene_id = "gene",
  group_id = "cluster",
  weight_id = "weight"
) {
  if (
    length(marker_top_n) != 1L ||
      !is.numeric(marker_top_n) ||
      is.na(marker_top_n) ||
      marker_top_n < 1
  ) {
    log_message(
      "{.arg marker_top_n} must be a single positive number",
      message_type = "error"
    )
  }
  marker_top_n <- as.integer(marker_top_n)
  marker_min_logfc <- as.numeric(marker_min_logfc)
  if (
    length(marker_min_logfc) != 1L ||
      is.na(marker_min_logfc)
  ) {
    log_message(
      "{.arg marker_min_logfc} must be a single number",
      message_type = "error"
    )
  }
  labels_chr <- as.character(labels)
  groups <- levels(labels)
  pcs <- 1e-6
  markers <- lapply(groups, function(group) {
    in_group <- labels_chr == group
    if (!any(in_group)) {
      return(NULL)
    }
    group_mean <- Matrix::rowMeans(ref_counts[, in_group, drop = FALSE])
    if (all(in_group)) {
      rest_mean <- rep(0, length(group_mean))
    } else {
      rest_mean <- Matrix::rowMeans(ref_counts[, !in_group, drop = FALSE])
    }
    logfc <- log2((group_mean + pcs) / (rest_mean + pcs))
    expressed <- group_mean > 0
    keep <- expressed & is.finite(logfc) & logfc >= marker_min_logfc
    if (!any(keep)) {
      keep <- expressed & is.finite(logfc)
    }
    if (!any(keep)) {
      return(NULL)
    }
    score <- logfc
    score[!is.finite(score)] <- 0
    score <- pmax(score, pcs)
    idx <- which(keep)
    idx <- idx[order(score[idx], group_mean[idx], decreasing = TRUE)]
    idx <- utils::head(idx, marker_top_n)
    data.frame(
      gene = rownames(ref_counts)[idx],
      cluster = group,
      weight = score[idx],
      stringsAsFactors = FALSE
    )
  })
  markers <- markers[!vapply(markers, is.null, logical(1))]
  if (length(markers) == 0L) {
    log_message(
      "Unable to generate marker genes for {.fn RunSPOTlight}",
      message_type = "error"
    )
  }
  mgs <- do.call(rbind, markers)
  rownames(mgs) <- NULL
  colnames(mgs) <- c(gene_id, group_id, weight_id)
  spotlight_validate_mgs_groups(mgs, labels, gene_id, group_id)
  mgs
}

spotlight_validate_mgs_groups <- function(mgs, labels, gene_id, group_id) {
  if (nrow(mgs) == 0L) {
    log_message(
      "{.arg mgs} must contain at least one valid marker gene",
      message_type = "error"
    )
  }
  missing_groups <- setdiff(levels(labels), unique(as.character(mgs[[group_id]])))
  if (length(missing_groups) > 0L) {
    log_message(
      "{.arg mgs} must contain marker genes for every reference cell type. Missing: {.val {missing_groups}}",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

spotlight_run_backend <- function(
  ref_counts,
  st_counts,
  labels,
  mgs,
  gene_id,
  group_id,
  weight_id,
  min_prop,
  scale,
  verbose,
  ...
) {
  spotlight_fun <- get_namespace_fun("SPOTlight", "SPOTlight")
  spotlight_fun(
    x = ref_counts,
    y = st_counts,
    groups = as.character(labels),
    mgs = mgs,
    gene_id = gene_id,
    group_id = group_id,
    weight_id = weight_id,
    min_prop = min_prop,
    scale = scale,
    verbose = verbose,
    ...
  )
}

spotlight_extract_weights <- function(result, spot_ids) {
  weights <- if (is.list(result) && !is.null(result$mat)) {
    result$mat
  } else {
    result
  }
  if (is.data.frame(weights)) {
    weights <- as.matrix(weights)
  }
  if (!is.matrix(weights)) {
    log_message(
      "{.pkg SPOTlight} did not return a valid proportion matrix",
      message_type = "error"
    )
  }
  rn_match <- if (is.null(rownames(weights))) 0L else sum(rownames(weights) %in% spot_ids)
  cn_match <- if (is.null(colnames(weights))) 0L else sum(colnames(weights) %in% spot_ids)
  if (cn_match > rn_match) {
    weights <- t(weights)
  }
  if (is.null(rownames(weights))) {
    log_message(
      "{.pkg SPOTlight} proportion matrix must contain spatial spot names",
      message_type = "error"
    )
  }
  spots_use <- spot_ids[spot_ids %in% rownames(weights)]
  if (length(spots_use) == 0L) {
    log_message(
      "{.pkg SPOTlight} proportion matrix could not be matched to spatial spot names",
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

spotlight_assert_scalar_string <- function(x, arg) {
  if (
    is.null(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !nzchar(x)
  ) {
    log_message(
      "{.arg {arg}} must be a single non-empty string",
      message_type = "error"
    )
  }
}
