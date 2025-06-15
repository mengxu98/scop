#' Run MDS (multi-dimensional scaling)
#'
#' @param object An object. This can be a Seurat object, an assay object, or a matrix-like object.
#' @param assay A character string specifying the assay to be used for the analysis. Default is NULL.
#' @param layer A character string specifying the layer to be used for the analysis. Default is "data".
#' @param features A character vector specifying the features to be used for the analysis. Default is NULL, which uses all variable features.
#' @param nmds An integer specifying the number of dimensions to be computed. Default is 50.
#' @param dist.method A character string specifying the distance metric to be used. Currently supported values are "euclidean", "chisquared","kullback", "jeffreys", "jensen", "manhattan", "maximum", "canberra", "minkowski", and "hamming". Default is "euclidean".
#' @param mds.method A character string specifying the MDS algorithm to be used. Currently supported values are "cmdscale", "isoMDS", and "sammon". Default is "cmdscale".
#' @param rev.mds A logical value indicating whether to perform reverse MDS (i.e., transpose the input matrix) before running the analysis. Default is FALSE.
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object. Default is "mds".
#' @param reduction.key A character string specifying the prefix for the column names of the basis vectors. Default is "MDS_".
#' @param verbose A logical value indicating whether to print verbose output. Default is TRUE.
#' @param seed.use An integer specifying the random seed to be used. Default is 11.
#' @param ... Additional arguments to be passed to the stats::\link[stats]{cmdscale}, MASS::\link[MASS]{isoMDS} or MASS::\link[MASS]{sammon} function.
#'
#' @examples
#' pancreas_sub <- RunMDS(object = pancreas_sub)
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "mds"
#' )
#'
#' @rdname RunMDS
#' @export
#'
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
  assay.data <- Seurat::GetAssay(object = object, assay = assay)
  reduction.data <- RunMDS(
    object = assay.data,
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
  object[[reduction.name]] <- reduction.data
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
  data.use <- GetAssayData5(object = object, layer = layer)
  features.var <- apply(
    X = data.use[features, ],
    MARGIN = 1,
    FUN = stats::var
  )
  features.keep <- features[features.var > 0]
  data.use <- data.use[features.keep, ]
  reduction.data <- RunMDS(
    object = data.use,
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
  return(reduction.data)
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
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (rev.mds) {
    object <- Matrix::t(object)
  }
  nmds <- min(nmds, nrow(x = object) - 1)
  x <- Matrix::t(
    Matrix::as.matrix(object)
  )
  cell.dist <- stats::as.dist(
    proxyC::dist(x = x, method = dist.method)
  )
  if (mds.method == "cmdscale") {
    mds.results <- stats::cmdscale(
      cell.dist,
      k = nmds,
      eig = TRUE,
      ...
    )
  }
  if (mds.method == "isoMDS") {
    check_r("MASS")
    mds.results <- MASS::isoMDS(cell.dist, k = nmds, ...)
  }
  if (mds.method == "sammon") {
    check_r("MASS")
    mds.results <- MASS::sammon(cell.dist, k = nmds, ...)
  }
  cell.embeddings <- mds.results$points

  rownames(x = cell.embeddings) <- colnames(x = object)
  colnames(x = cell.embeddings) <- paste0(reduction.key, 1:nmds)
  reduction.data <- SeuratObject::CreateDimReducObject(
    embeddings = cell.embeddings,
    assay = assay,
    key = reduction.key,
    misc = list(slot = layer, mds.results = mds.results)
  )
  return(reduction.data)
}
