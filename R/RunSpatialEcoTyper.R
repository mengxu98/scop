#' @title Run SpatialEcoTyper spatial ecotype discovery
#'
#' @description
#' Run single-sample de novo spatial ecotype discovery with the optional
#' `SpatialEcoTyper` package and write spatial ecotype labels back to a
#' `Seurat` object.
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param celltype.by Metadata column containing cell type annotations.
#' @param x.by,y.by Metadata columns containing single-cell spatial coordinates.
#' @param features Optional feature vector used to subset the expression matrix.
#' @param outprefix Output prefix passed to `SpatialEcoTyper`. Use `NULL` to
#' avoid writing result files to the working directory.
#' @param radius Spatial neighborhood radius, in the same units as `x.by` and
#' `y.by`.
#' @param resolution Louvain clustering resolution used by `SpatialEcoTyper`.
#' @param nfeatures Number of variable features used by `SpatialEcoTyper`.
#' @param min.cts.per.region Minimum number of cell types required in a spatial
#' neighborhood.
#' @param npcs Number of principal components used for similarity networks.
#' @param min.cells Minimum number of cells or spatial meta-cells expressing a
#' feature.
#' @param min.features Minimum number of features detected in a cell or spatial
#' meta-cell.
#' @param iterations Number of similarity network fusion iterations.
#' @param minibatch Number of columns processed per mini-batch in SNF.
#' @param ncores Number of CPU cores used by `SpatialEcoTyper`.
#' @param grid.size Spatial grid size used to discretize coordinates.
#' @param filter.region.by.celltypes Optional cell types used to restrict spatial
#' neighborhoods.
#' @param k Number of spatial nearest neighbors used to construct spatial
#' meta-cells.
#' @param k.sn Number of nearest neighbors used to construct similarity networks.
#' @param dropcell Whether cells without spatial ecotype assignments are removed
#' from the returned `SpatialEcoTyper` metadata.
#' @param prefix Prefix used for the output metadata column.
#' @param tool_name Name used to store detailed results in `srt@tools`.
#' @param store_results Whether to store raw results in `srt@tools`.
#' @param allow_partial Whether to allow missing SE labels for cells absent from
#' the returned `SpatialEcoTyper` metadata. Default is `FALSE` to avoid silent
#' partial annotations.
#' @param ... Additional arguments passed to
#' `SpatialEcoTyper::SpatialEcoTyper()`.
#'
#' @return A `Seurat` object with spatial ecotype labels in metadata and raw
#' results stored in `srt@tools[[tool_name]]` when `store_results = TRUE`.
#' @export
#'
#' @examples
#' \dontrun{
#' srt <- RunSpatialEcoTyper(
#'   srt,
#'   celltype.by = "CellType",
#'   x.by = "X",
#'   y.by = "Y"
#' )
#' }
RunSpatialEcoTyper <- function(
  srt,
  assay = NULL,
  layer = "data",
  celltype.by,
  x.by = "X",
  y.by = "Y",
  features = NULL,
  outprefix = NULL,
  radius = 50,
  resolution = 0.5,
  nfeatures = 300,
  min.cts.per.region = 2,
  npcs = 20,
  min.cells = 5,
  min.features = 10,
  iterations = 10,
  minibatch = 5000,
  ncores = 4,
  grid.size = round(radius * 1.4),
  filter.region.by.celltypes = NULL,
  k = 20,
  k.sn = 50,
  dropcell = FALSE,
  prefix = "SpatialEcoTyper",
  tool_name = "SpatialEcoTyper",
  store_results = TRUE,
  allow_partial = FALSE,
  verbose = TRUE,
  ...
) {
  check_r("digitalcytometry/SpatialEcoTyper", verbose = FALSE)

  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  spatialecotyper_assert_scalar_string(celltype.by, "celltype.by")
  spatialecotyper_assert_scalar_string(x.by, "x.by")
  spatialecotyper_assert_scalar_string(y.by, "y.by")
  spatialecotyper_assert_scalar_string(prefix, "prefix")
  spatialecotyper_assert_scalar_string(tool_name, "tool_name")

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  normdata <- GetAssayData5(srt, assay = assay, layer = layer)
  if (is.null(rownames(normdata)) || is.null(colnames(normdata))) {
    log_message(
      "{.arg layer} must contain feature and cell names",
      message_type = "error"
    )
  }

  if (!is.null(features)) {
    features <- unique(as.character(features))
    missing_features <- setdiff(features, rownames(normdata))
    features <- intersect(features, rownames(normdata))
    if (length(features) == 0L) {
      log_message(
        "No requested {.arg features} are present in the selected assay/layer",
        message_type = "error"
      )
    }
    if (length(missing_features) > 0L) {
      log_message(
        "Ignoring {.val {length(missing_features)}} requested features absent from the selected assay/layer",
        message_type = "warning",
        verbose = verbose
      )
    }
    normdata <- normdata[features, , drop = FALSE]
  }

  meta <- srt[[]]
  needed_cols <- c(celltype.by, x.by, y.by)
  missing_cols <- setdiff(needed_cols, colnames(meta))
  if (length(missing_cols) > 0L) {
    log_message(
      "Missing metadata column{?s}: {.val {missing_cols}}",
      message_type = "error"
    )
  }
  missing_cells <- setdiff(colnames(normdata), rownames(meta))
  if (length(missing_cells) > 0L) {
    log_message(
      "Metadata is missing {.val {length(missing_cells)}} cell{?s} from the selected assay/layer",
      message_type = "error"
    )
  }

  metadata <- spatialecotyper_build_metadata(
    meta = meta[colnames(normdata), , drop = FALSE],
    celltype.by = celltype.by,
    x.by = x.by,
    y.by = y.by
  )

  log_message(
    "Run {.pkg SpatialEcoTyper} on {.val {ncol(normdata)}} cells and {.val {nrow(normdata)}} features",
    verbose = verbose
  )
  spatialecotyper_fun <- get_namespace_fun(
    "SpatialEcoTyper",
    "SpatialEcoTyper"
  )
  result <- spatialecotyper_fun(
    normdata = normdata,
    metadata = metadata,
    outprefix = outprefix,
    radius = radius,
    resolution = resolution,
    nfeatures = nfeatures,
    min.cts.per.region = min.cts.per.region,
    npcs = npcs,
    min.cells = min.cells,
    min.features = min.features,
    iterations = iterations,
    minibatch = minibatch,
    ncores = ncores,
    grid.size = grid.size,
    filter.region.by.celltypes = filter.region.by.celltypes,
    k = k,
    k.sn = k.sn,
    dropcell = dropcell,
    ...
  )

  result_metadata <- spatialecotyper_extract_result_metadata(result)
  se_col <- paste0(prefix, "_SE")
  se_labels <- spatialecotyper_match_se_labels(
    result_metadata = result_metadata,
    cells = colnames(srt),
    allow_partial = allow_partial,
    verbose = verbose
  )
  add_meta <- data.frame(
    se_labels,
    row.names = names(se_labels),
    check.names = FALSE
  )
  colnames(add_meta) <- se_col
  srt <- Seurat::AddMetaData(srt, metadata = add_meta)

  if (isTRUE(store_results)) {
    srt@tools[[tool_name]] <- list(
      result = result,
      metadata = result_metadata,
      labels = se_labels,
      parameters = list(
        assay = assay,
        layer = layer,
        celltype.by = celltype.by,
        x.by = x.by,
        y.by = y.by,
        features = features,
        outprefix = outprefix,
        radius = radius,
        resolution = resolution,
        nfeatures = nfeatures,
        min.cts.per.region = min.cts.per.region,
        npcs = npcs,
        min.cells = min.cells,
        min.features = min.features,
        iterations = iterations,
        minibatch = minibatch,
        ncores = ncores,
        grid.size = grid.size,
        filter.region.by.celltypes = filter.region.by.celltypes,
        k = k,
        k.sn = k.sn,
        dropcell = dropcell,
        prefix = prefix,
        tool_name = tool_name,
        allow_partial = allow_partial
      )
    )
  }

  log_message(
    "{.pkg SpatialEcoTyper} SE labels stored in metadata column {.val {se_col}}",
    message_type = "success",
    verbose = verbose
  )
  srt
}

spatialecotyper_assert_scalar_string <- function(x, arg) {
  if (
    missing(x) ||
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

spatialecotyper_build_metadata <- function(meta, celltype.by, x.by, y.by) {
  x <- meta[[x.by]]
  y <- meta[[y.by]]
  celltype <- meta[[celltype.by]]
  if (is.factor(x)) {
    x <- as.character(x)
  }
  if (is.factor(y)) {
    y <- as.character(y)
  }

  metadata <- data.frame(
    X = suppressWarnings(as.numeric(x)),
    Y = suppressWarnings(as.numeric(y)),
    CellType = as.character(celltype),
    row.names = rownames(meta),
    stringsAsFactors = FALSE
  )
  invalid <- is.na(metadata$X) |
    is.na(metadata$Y) |
    is.na(metadata$CellType) |
    !nzchar(metadata$CellType)
  if (any(invalid)) {
    log_message(
      "SpatialEcoTyper metadata contains {.val {sum(invalid)}} invalid cell{?s}",
      message_type = "error"
    )
  }
  metadata
}

spatialecotyper_extract_result_metadata <- function(result) {
  if (
    !is.list(result) ||
      is.null(result$metadata) ||
      !is.data.frame(result$metadata)
  ) {
    log_message(
      "{.pkg SpatialEcoTyper} did not return a valid {.field metadata} data frame",
      message_type = "error"
    )
  }
  metadata <- result$metadata
  if (!"SE" %in% colnames(metadata)) {
    log_message(
      "{.pkg SpatialEcoTyper} did not return metadata with an {.field SE} column",
      message_type = "error"
    )
  }
  if (is.null(rownames(metadata)) || any(!nzchar(rownames(metadata)))) {
    log_message(
      "{.pkg SpatialEcoTyper} returned metadata without cell row names",
      message_type = "error"
    )
  }
  metadata
}

spatialecotyper_match_se_labels <- function(
  result_metadata,
  cells,
  allow_partial = FALSE,
  verbose = TRUE
) {
  se_labels <- stats::setNames(rep(NA_character_, length(cells)), cells)
  common <- intersect(cells, rownames(result_metadata))
  missing_cells <- setdiff(cells, rownames(result_metadata))
  if (length(missing_cells) > 0L && !isTRUE(allow_partial)) {
    log_message(
      "{.pkg SpatialEcoTyper} returned no SE label for {.val {length(missing_cells)}} cell{?s}. Set {.arg allow_partial = TRUE} to keep missing labels as {.val NA}.",
      message_type = "error"
    )
  }
  if (length(missing_cells) > 0L) {
    log_message(
      "{.pkg SpatialEcoTyper} returned no SE label for {.val {length(missing_cells)}} cell{?s}; storing {.val NA}",
      message_type = "warning",
      verbose = verbose
    )
  }
  se_values <- as.character(result_metadata[common, "SE"])
  se_labels[common] <- se_values
  missing_se <- names(se_labels)[is.na(se_labels) | !nzchar(se_labels)]
  if (length(missing_se) > 0L && !isTRUE(allow_partial)) {
    log_message(
      "{.pkg SpatialEcoTyper} returned missing SE values for {.val {length(missing_se)}} cell{?s}. Set {.arg allow_partial = TRUE} to keep missing labels as {.val NA}.",
      message_type = "error"
    )
  }
  se_labels
}
