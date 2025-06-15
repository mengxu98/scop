#' Run GLMPCA (generalized version of principal components analysis)
#'
#' @md
#' @param object An object. This can be a Seurat object, an assay object, or a matrix-like object.
#' @param assay A character string specifying the assay to be used for the analysis.
#' Default is NULL.
#' @param layer A character string specifying the layer to be used for the analysis.
#' Default is "counts".
#' @param features A character vector specifying the features to be used for the analysis.
#' Default is NULL, which uses all variable features.
#' @param L An integer specifying the number of components to be computed.
#' Default is 5.
#' @param fam A character string specifying the family of the generalized linear model to be used.
#' Currently supported values are "poi", "nb", "nb2", "binom", "mult", and "bern".
#' Default is "poi".
#' @param rev.gmlpca A logical value indicating whether to perform reverse GLMPCA (i.e., transpose the input matrix) before running the analysis. Default is FALSE.
#' @param ndims.print An integer vector specifying the dimensions (number of components) to print in the output.
#' Default is 1:5.
#' @param nfeatures.print An integer specifying the number of features to print in the output.
#' Default is 30.
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object.
#' Default is "glmpca".
#' @param reduction.key A character string specifying the prefix for the column names of the basis vectors.
#' Default is "GLMPC_".
#' @param verbose A logical value indicating whether to print verbose output.
#' Default is TRUE.
#' @param seed.use An integer specifying the random seed to be used.
#' Default is 11.
#' @param ... Additional arguments to be passed to the [glmpca::glmpca] function.
#'
#' @examples
#' pancreas_sub <- RunGLMPCA(object = pancreas_sub)
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "glmpca"
#' )
#'
#' @rdname RunGLMPCA
#' @export
RunGLMPCA <- function(object, ...) {
  UseMethod(generic = "RunGLMPCA", object = object)
}

#' @rdname RunGLMPCA
#' @method RunGLMPCA Seurat
#' @export
RunGLMPCA.Seurat <- function(
  object,
  assay = NULL,
  layer = "counts",
  features = NULL,
  L = 5,
  fam = c("poi", "nb", "nb2", "binom", "mult", "bern"),
  rev.gmlpca = FALSE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.name = "glmpca",
  reduction.key = "GLMPC_",
  verbose = TRUE,
  seed.use = 11,
  ...
) {
  features <- features %||% SeuratObject::VariableFeatures(object = object)
  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  assay.data <- Seurat::GetAssay(object = object, assay = assay)
  reduction.data <- RunGLMPCA(
    object = assay.data,
    assay = assay,
    layer = layer,
    L = L,
    fam = fam,
    rev.gmlpca = rev.gmlpca,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    verbose = verbose,
    ...
  )
  object[[reduction.name]] <- reduction.data
  object <- Seurat::LogSeuratCommand(object = object)
  return(object)
}

#' @rdname RunGLMPCA
#' @method RunGLMPCA Assay
#' @export
RunGLMPCA.Assay <- function(
  object,
  assay = NULL,
  layer = "counts",
  features = NULL,
  L = 5,
  fam = c("poi", "nb", "nb2", "binom", "mult", "bern"),
  rev.gmlpca = FALSE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "GLMPC_",
  verbose = TRUE,
  seed.use = 11,
  ...
) {
  features <- features %||% SeuratObject::VariableFeatures(object = object)
  data.use <- GetAssayData5(object = object, layer = layer)
  features.var <- apply(
    X = data.use[features, ],
    MARGIN = 1,
    FUN = stats::var
  )
  features.keep <- features[features.var > 0]
  data.use <- data.use[features.keep, ]
  reduction.data <- RunGLMPCA(
    object = data.use,
    assay = assay,
    layer = layer,
    L = L,
    fam = fam,
    rev.gmlpca = rev.gmlpca,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    verbose = verbose,
    ...
  )
  return(reduction.data)
}

#' @rdname RunGLMPCA
#' @method RunGLMPCA default
#' @export
RunGLMPCA.default <- function(
  object,
  assay = NULL,
  layer = "counts",
  features = NULL,
  L = 5,
  fam = c("poi", "nb", "nb2", "binom", "mult", "bern"),
  rev.gmlpca = FALSE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "GLMPC_",
  verbose = TRUE,
  seed.use = 11,
  ...
) {
  check_r("glmpca")
  if (inherits(object, "dgCMatrix")) {
    object <- Matrix::as.matrix(object)
  }
  fam <- match.arg(fam)
  glmpca_results <- glmpca::glmpca(Y = object, L = L, fam = fam, ...)
  glmpca_dimnames <- paste0(reduction.key, seq_len(L))
  factors <- Matrix::as.matrix(glmpca_results$factors)
  loadings <- Matrix::as.matrix(glmpca_results$loadings)
  colnames(x = factors) <- glmpca_dimnames
  colnames(x = loadings) <- glmpca_dimnames
  factors_l2_norm <- sqrt(Matrix::colSums(factors^2))
  class(glmpca_results) <- NULL
  glmpca_results$factors <- glmpca_results$loadings <- NULL
  reduction.data <- Seurat::CreateDimReducObject(
    embeddings = factors,
    key = reduction.key,
    loadings = loadings,
    stdev = factors_l2_norm,
    assay = assay,
    global = TRUE,
    misc = list(slot = layer, glmpca.results = glmpca_results)
  )
  if (verbose) {
    msg <- utils::capture.output(
      print(
        x = reduction.data,
        dims = ndims.print,
        nfeatures = nfeatures.print
      )
    )
    message(paste(msg, collapse = "\n"))
  }
  return(reduction.data)
}
