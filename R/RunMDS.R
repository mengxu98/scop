#' @title Run MDS (multi-dimensional scaling)
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams RunUMAP2
#' @param object An object. This can be a Seurat object, an assay object, or a matrix-like object.
#' @param nmds The number of dimensions to be computed.
#' Default is `50`.
#' @param dist.method The distance metric to be used.
#' Currently supported values are `"euclidean"`, `"chisquared"`, `"kullback"`, `"jeffreys"`,
#' `"jensen"`, `"manhattan"`, `"maximum"`, `"canberra"`, `"minkowski"`, and `"hamming"`.
#' Default is `"euclidean"`.
#' @param mds.method The MDS algorithm to be used.
#' Currently supported values are `"cmdscale"`, `"isoMDS"`, and `"sammon"`.
#' Default is `"cmdscale"`.
#' @param rev.mds Whether to perform reverse MDS (i.e., transpose the input matrix) before running the analysis.
#' Default is `FALSE`.
#' @param reduction.name The name of the reduction to be stored in the Seurat object.
#' Default is `"mds"`.
#' @param reduction.key The prefix for the column names of the basis vectors.
#' Default is `"MDS_"`.
#' @param ... Additional arguments to be passed to [stats::cmdscale], [MASS::isoMDS] or [MASS::sammon].
#'
#' @rdname RunMDS
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunMDS(pancreas_sub)
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "mds"
#' )
RunMDS <- function(object, ...) {
  UseMethod(generic = "RunMDS", object = object)
}

#' @rdname RunMDS
#' @method RunMDS Seurat
#' @export
RunMDS.Seurat <- function(
    object,
    assay = NULL,
    layer = "data",
    features = NULL,
    nmds = 50,
    dist.method = "euclidean",
    mds.method = "cmdscale",
    rev.mds = FALSE,
    reduction.name = "mds",
    reduction.key = "MDS_",
    verbose = TRUE,
    seed.use = 11,
    ...) {
  features <- features %||% SeuratObject::VariableFeatures(object = object)
  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  assay_data <- Seurat::GetAssay(
    object = object,
    assay = assay
  )
  reduction_data <- RunMDS(
    object = assay_data,
    assay = assay,
    layer = layer,
    features = features,
    nmds = nmds,
    dist.method = dist.method,
    mds.method = mds.method,
    rev.mds = rev.mds,
    verbose = verbose,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  object[[reduction.name]] <- reduction_data
  object <- Seurat::LogSeuratCommand(object = object)
  return(object)
}

#' @rdname RunMDS
#' @method RunMDS Assay
#' @export
RunMDS.Assay <- function(
    object,
    assay = NULL,
    layer = "data",
    features = NULL,
    nmds = 50,
    dist.method = "euclidean",
    mds.method = "cmdscale",
    rev.mds = FALSE,
    reduction.key = "MDS_",
    verbose = TRUE,
    seed.use = 11,
    ...) {
  features <- features %||% SeuratObject::VariableFeatures(object = object)
  data_use <- GetAssayData5(
    object = object,
    layer = layer
  )
  features_var <- apply(
    X = data_use[features, ],
    MARGIN = 1,
    FUN = stats::var
  )
  features_keep <- features[features_var > 0]
  data_use <- data_use[features_keep, ]
  reduction_data <- RunMDS(
    object = data_use,
    assay = assay,
    layer = layer,
    nmds = nmds,
    dist.method = dist.method,
    mds.method = mds.method,
    rev.mds = rev.mds,
    verbose = verbose,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  return(reduction_data)
}

#' @rdname RunMDS
#' @method RunMDS Assay5
#' @export
RunMDS.Assay5 <- function(
    object,
    assay = NULL,
    layer = "data",
    features = NULL,
    nmds = 50,
    dist.method = "euclidean",
    mds.method = "cmdscale",
    rev.mds = FALSE,
    reduction.key = "MDS_",
    verbose = TRUE,
    seed.use = 11,
    ...) {
  features <- features %||% SeuratObject::VariableFeatures(object = object)
  data_use <- GetAssayData5(
    object = object,
    layer = layer
  )
  features_var <- apply(
    X = data_use[features, ],
    MARGIN = 1,
    FUN = stats::var
  )
  features_keep <- features[features_var > 0]
  data_use <- data_use[features_keep, ]
  reduction_data <- RunMDS(
    object = data_use,
    assay = assay,
    layer = layer,
    nmds = nmds,
    dist.method = dist.method,
    mds.method = mds.method,
    rev.mds = rev.mds,
    verbose = verbose,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  return(reduction_data)
}

#' @rdname RunMDS
#' @method RunMDS default
#' @export
RunMDS.default <- function(
    object,
    assay = NULL,
    layer = "data",
    nmds = 50,
    dist.method = "euclidean",
    mds.method = "cmdscale",
    rev.mds = FALSE,
    reduction.key = "MDS_",
    verbose = TRUE,
    seed.use = 11,
    ...) {
  if (!is.null(seed.use)) {
    set.seed(seed = seed.use)
  }
  if (rev.mds) {
    object <- Matrix::t(object)
  }
  nmds <- min(nmds, nrow(x = object) - 1)
  x <- Matrix::t(
    as_matrix(object)
  )
  cell_dist <- stats::as.dist(
    proxyC::dist(x = x, method = dist.method)
  )
  if (mds.method == "cmdscale") {
    mds_results <- stats::cmdscale(
      cell_dist,
      k = nmds,
      eig = TRUE,
      ...
    )
  }
  if (mds.method == "isoMDS") {
    check_r("MASS", verbose = FALSE)
    mds_results <- MASS::isoMDS(cell_dist, k = nmds, ...)
  }
  if (mds.method == "sammon") {
    check_r("MASS", verbose = FALSE)
    mds_results <- MASS::sammon(cell_dist, k = nmds, ...)
  }
  cell_embeddings <- mds_results$points
  rownames(x = cell_embeddings) <- colnames(x = object)
  colnames(x = cell_embeddings) <- paste0(reduction.key, seq_len(nmds))
  reduction_data <- SeuratObject::CreateDimReducObject(
    embeddings = cell_embeddings,
    assay = assay,
    key = reduction.key,
    misc = list(
      slot = layer,
      mds.results = mds_results
    )
  )
  return(reduction_data)
}
