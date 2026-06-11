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
#' @param outprefix Output prefix passed to `SpatialEcoTyper`. The default
#' `NULL` avoids writing result files to the working directory.
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
#' @param prefix Prefix used for the output metadata column and `srt@tools` entry.
#' @param store_results Whether to store raw results in `srt@tools`.
#'
#' @return A `Seurat` object with spatial ecotype labels in metadata and raw
#' results stored in `srt@tools[[prefix]]` when `store_results = TRUE`.
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
  store_results = TRUE,
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  if (missing(celltype.by) || length(celltype.by) != 1L || is.na(celltype.by)) {
    log_message("{.arg celltype.by} must be a single metadata column", message_type = "error")
  }
  if (!requireNamespace("SpatialEcoTyper", quietly = TRUE)) {
    log_message("Please install {.pkg SpatialEcoTyper} before running {.fn RunSpatialEcoTyper}", message_type = "error")
  }

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  meta <- srt[[]]
  missing_cols <- setdiff(c(celltype.by, x.by, y.by), colnames(meta))
  if (length(missing_cols) > 0L) {
    log_message("Missing metadata column{?s}: {.val {missing_cols}}", message_type = "error")
  }

  normdata <- GetAssayData5(srt, assay = assay, layer = layer)
  if (!is.null(features)) {
    features <- intersect(features, rownames(normdata))
    if (length(features) == 0L) {
      log_message("No requested {.arg features} are present", message_type = "error")
    }
    normdata <- normdata[features, , drop = FALSE]
  }
  if (!inherits(normdata, "Matrix") && !is.matrix(normdata)) {
    normdata <- as.matrix(normdata)
  }

  meta <- meta[colnames(normdata), , drop = FALSE]
  x <- meta[[x.by]]
  y <- meta[[y.by]]
  if (is.factor(x)) x <- as.character(x)
  if (is.factor(y)) y <- as.character(y)
  metadata <- data.frame(
    X = suppressWarnings(as.numeric(x)),
    Y = suppressWarnings(as.numeric(y)),
    CellType = as.character(meta[[celltype.by]]),
    row.names = rownames(meta),
    stringsAsFactors = FALSE
  )
  invalid <- is.na(metadata$X) | is.na(metadata$Y) |
    is.na(metadata$CellType) | !nzchar(metadata$CellType)
  if (any(invalid)) {
    log_message(
      "SpatialEcoTyper metadata contains {.val {sum(invalid)}} invalid cell{?s}",
      message_type = "error"
    )
  }
  if (!identical(rownames(metadata), colnames(normdata))) {
    log_message(
      "SpatialEcoTyper metadata row names must match expression matrix column names",
      message_type = "error"
    )
  }

  log_message(
    "Run {.pkg SpatialEcoTyper} on {.val {ncol(normdata)}} cells and {.val {nrow(normdata)}} features",
    verbose = verbose
  )
  spatialecotyper_fun <- getExportedValue("SpatialEcoTyper", "SpatialEcoTyper")
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
    dropcell = dropcell
  )
  if (!is.list(result) || is.null(result$metadata) || !"SE" %in% colnames(result$metadata)) {
    log_message("{.pkg SpatialEcoTyper} did not return metadata with an {.val SE} column", message_type = "error")
  }

  se_col <- paste0(prefix, "_SE")
  cells <- colnames(srt)
  se_labels <- stats::setNames(rep(NA_character_, length(cells)), cells)
  common <- intersect(cells, rownames(result$metadata))
  se_labels[common] <- as.character(result$metadata[common, "SE"])
  add_meta <- data.frame(SE = unname(se_labels), row.names = cells, stringsAsFactors = FALSE)
  colnames(add_meta) <- se_col
  srt <- Seurat::AddMetaData(srt, metadata = add_meta)

  if (isTRUE(store_results)) {
    srt@tools[[prefix]] <- list(
      result = result,
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
        prefix = prefix
      )
    )
  }

  log_message(
    "{.pkg SpatialEcoTyper} SE labels stored in metadata column {.val {se_col}}",
    verbose = verbose
  )
  srt
}
