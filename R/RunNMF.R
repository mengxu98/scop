#' Run NMF (non-negative matrix factorization)
#'
#' @md
#' @param object An object. This can be a Seurat object, an Assay object, or a matrix-like object.
#' @param assay A character string specifying the assay to be used for the analysis. Default is NULL.
#' @param layer A character string specifying the layer to be used for the analysis. Default is "data".
#' @param features A character vector specifying the features to be used for the analysis. Default is NULL, which uses all variable features.
#' @param nbes An integer specifying the number of basis vectors (components) to be computed. Default is 50.
#' @param nmf.method A character string specifying the NMF algorithm to be used. Currently supported values are "RcppML" and "NMF". Default is "RcppML".
#' @param tol A numeric value specifying the tolerance for convergence (only applicable when nmf.method is "RcppML"). Default is 1e-5.
#' @param maxit An integer specifying the maximum number of iterations for convergence (only applicable when nmf.method is "RcppML"). Default is 100.
#' @param rev.nmf A logical value indicating whether to perform reverse NMF (i.e., transpose the input matrix) before running the analysis. Default is FALSE.
#' @param ndims.print An integer vector specifying the dimensions (number of basis vectors) to print in the output. Default is 1:5.
#' @param nfeatures.print An integer specifying the number of features to print in the output. Default is 30.
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object. Default is "nmf".
#' @param reduction.key A character string specifying the prefix for the column names of the basis vectors. Default is "BE_".
#' @param verbose A logical value indicating whether to print verbose output. Default is TRUE.
#' @param seed.use An integer specifying the random seed to be used. Default is 11.
#' @param ... Additional arguments passed to [RcppML::nmf] or [NMF::nmf] function.
#'
#' @rdname RunNMF
#' @export
#'
#' @examples
#' library(Matrix)
#' pancreas_sub <- RunNMF(object = pancreas_sub)
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "nmf"
#' )
RunNMF <- function(object, ...) {
  UseMethod(generic = "RunNMF", object = object)
}

#' @rdname RunNMF
#' @method RunNMF Seurat
#' @export
RunNMF.Seurat <- function(
    object,
    assay = NULL,
    layer = "data",
    features = NULL,
    nbes = 50,
    nmf.method = "RcppML",
    tol = 1e-5,
    maxit = 100,
    rev.nmf = FALSE,
    ndims.print = 1:5,
    nfeatures.print = 30,
    reduction.name = "nmf",
    reduction.key = "BE_",
    verbose = TRUE,
    seed.use = 11,
    ...) {
  features <- features %||% SeuratObject::VariableFeatures(object = object)
  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  assay_data <- Seurat::GetAssay(object = object, assay = assay)
  reduction_data <- RunNMF(
    object = assay_data,
    assay = assay,
    layer = layer,
    features = features,
    nbes = nbes,
    nmf.method = nmf.method,
    tol = tol,
    maxit = maxit,
    rev.nmf = rev.nmf,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  object[[reduction.name]] <- reduction_data
  object <- Seurat::LogSeuratCommand(object = object)
  return(object)
}

#' @rdname RunNMF
#' @method RunNMF Assay
#' @export
RunNMF.Assay <- function(
    object,
    assay = NULL,
    layer = "data",
    features = NULL,
    nbes = 50,
    nmf.method = "RcppML",
    tol = 1e-5,
    maxit = 100,
    rev.nmf = FALSE,
    ndims.print = 1:5,
    nfeatures.print = 30,
    reduction.key = "BE_",
    verbose = TRUE,
    seed.use = 11,
    ...) {
  features <- features %||% SeuratObject::VariableFeatures(object = object)
  data_use <- GetAssayData5(object = object, layer = layer)
  features_var <- apply(
    X = data_use[features, ],
    MARGIN = 1,
    FUN = stats::var
  )
  features_keep <- features[features_var > 0]
  data_use <- data_use[features_keep, ]
  reduction_data <- RunNMF(
    object = data_use,
    assay = assay,
    layer = layer,
    nbes = nbes,
    nmf.method = nmf.method,
    tol = tol,
    maxit = maxit,
    rev.nmf = rev.nmf,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  return(reduction_data)
}

#' @rdname RunNMF
#' @method RunNMF Assay5
#' @export
RunNMF.Assay5 <- function(
    object,
    assay = NULL,
    layer = "data",
    features = NULL,
    nbes = 50,
    nmf.method = "RcppML",
    tol = 1e-5,
    maxit = 100,
    rev.nmf = FALSE,
    ndims.print = 1:5,
    nfeatures.print = 30,
    reduction.key = "BE_",
    verbose = TRUE,
    seed.use = 11,
    ...) {
  features <- features %||% SeuratObject::VariableFeatures(object = object)
  data_use <- GetAssayData5(object = object, layer = layer)
  features_var <- apply(
    X = data_use[features, ],
    MARGIN = 1,
    FUN = stats::var
  )
  features_keep <- features[features_var > 0]
  data_use <- data_use[features_keep, ]
  reduction_data <- RunNMF(
    object = data_use,
    assay = assay,
    layer = layer,
    nbes = nbes,
    nmf.method = nmf.method,
    tol = tol,
    maxit = maxit,
    rev.nmf = rev.nmf,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  return(reduction_data)
}

#' @rdname RunNMF
#' @method RunNMF default
#' @export
RunNMF.default <- function(
    object,
    assay = NULL,
    layer = "data",
    nbes = 50,
    nmf.method = "RcppML",
    tol = 1e-5,
    maxit = 100,
    rev.nmf = FALSE,
    ndims.print = 1:5,
    nfeatures.print = 30,
    reduction.key = "BE_",
    verbose = TRUE,
    seed.use = 11,
    ...) {
  set.seed(seed = seed.use)

  if (rev.nmf) {
    object <- Matrix::t(x = object)
  }
  nbes <- min(nbes, nrow(x = object) - 1)
  if (nmf.method == "RcppML") {
    check_r("zdebruine/RcppML")
    options("RcppML.verbose" = verbose)
    options("RcppML.threads" = 0)

    nmf_results <- RcppML::nmf(
      Matrix::t(object),
      k = nbes,
      tol = tol,
      maxit = maxit,
      verbose = verbose,
      ...
    )
    cell_embeddings <- nmf_results$w
    feature_loadings <- Matrix::t(nmf_results$h)
  }

  if (nmf.method == "NMF") {
    check_r("NMF")
    nmf_results <- NMF::nmf(
      x = Matrix::as.matrix(
        Matrix::t(object)
      ),
      rank = nbes
    )
    cell_embeddings <- nmf_results@fit@W
    feature_loadings <- Matrix::t(nmf_results@fit@H)
  }

  rownames(x = feature_loadings) <- rownames(x = object)
  colnames(x = feature_loadings) <- paste0(reduction.key, 1:nbes)
  rownames(x = cell_embeddings) <- colnames(x = object)
  colnames(x = cell_embeddings) <- colnames(x = feature_loadings)
  reduction_data <- Seurat::CreateDimReducObject(
    embeddings = cell_embeddings,
    loadings = feature_loadings,
    assay = assay,
    key = reduction.key,
    misc = list(
      slot = layer,
      nmf_results = nmf_results
    )
  )
  if (verbose) {
    msg <- utils::capture.output(
      print(
        x = reduction_data,
        dims = ndims.print,
        nfeatures = nfeatures.print
      )
    )
    log_message(
      paste(msg, collapse = "\n"),
      multiline_indent = FALSE,
      timestamp = FALSE
    )
  }

  return(reduction_data)
}
