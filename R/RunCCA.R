#' @export
RunCCA.default <- function(
  object1,
  object2,
  standardize = TRUE,
  num.cc = 20,
  seed.use = 42,
  verbose = FALSE,
  ...
) {
  if (!is.null(seed.use)) {
    set.seed(seed.use)
  }
  cells1 <- colnames(object1)
  cells2 <- colnames(object2)
  if (standardize) {
    standardize_fun <- get_namespace_fun("Seurat", "Standardize")
    object1 <- standardize_fun(mat = object1, display_progress = FALSE)
    object2 <- standardize_fun(mat = object2, display_progress = FALSE)
  }
  if (inherits(object1, "sparseMatrix")) {
    object1 <- as.matrix(object1)
  }
  if (inherits(object2, "sparseMatrix")) {
    object2 <- as.matrix(object2)
  }

  cca_mult <- function(L, R) {
    Lm <- if (is.null(dim(L))) matrix(L, nrow = 1) else L
    Rm <- if (is.null(dim(R))) matrix(R, ncol = 1) else R
    matrix_product(Lm, Rm)
  }
  storage.mode(object1) <- "double"
  storage.mode(object2) <- "double"
  A <- cca_crossprod_matrix(object1, object2)
  nc <- min(as.integer(num.cc), nrow(A) - 1L, ncol(A) - 1L)
  ir <- irlba::irlba(A, nv = nc, nu = nc, mult = cca_mult)
  ord <- order(-ir$d)
  U <- ir$u[, ord, drop = FALSE]
  V <- ir$v[, ord, drop = FALSE]
  cca.data <- rbind(U, V)
  colnames(cca.data) <- paste0("CC", seq_len(nc))
  rownames(cca.data) <- c(cells1, cells2)
  cca.data <- apply(cca.data, 2, function(x) {
    if (sign(x[1]) == -1) {
      x <- x * -1
    }
    x
  })
  list(ccv = cca.data, d = ir$d[ord])
}

#' @export
RunCCA.Seurat <- function(
  object1,
  object2,
  assay1 = NULL,
  assay2 = NULL,
  num.cc = 20,
  features = NULL,
  renormalize = FALSE,
  rescale = FALSE,
  compute.gene.loadings = TRUE,
  add.cell.id1 = NULL,
  add.cell.id2 = NULL,
  verbose = TRUE,
  ...
) {
  if (!inherits(object2, "Seurat")) {
    stop("RunCCA.Seurat requires two Seurat objects.", call. = FALSE)
  }
  if (
    isTRUE(renormalize) || isTRUE(rescale) || !isTRUE(compute.gene.loadings)
  ) {
    stop(
      "RunCCA.Seurat supports renormalize = FALSE, rescale = FALSE, and compute.gene.loadings = TRUE.",
      call. = FALSE
    )
  }
  op <- options(
    Seurat.object.assay.version = "v3",
    Seurat.object.assay.calcn = FALSE
  )
  on.exit(options(op), add = TRUE)
  if (is.null(assay1)) {
    assay1 <- SeuratObject::DefaultAssay(object1)
  }
  if (is.null(assay2)) {
    assay2 <- SeuratObject::DefaultAssay(object2)
  }
  if (is.null(features)) {
    features <- union(
      SeuratObject::VariableFeatures(object1),
      SeuratObject::VariableFeatures(object2)
    )
  }
  data.use1 <- SeuratObject::GetAssayData(
    object1,
    assay = assay1,
    layer = "scale.data"
  )
  data.use2 <- SeuratObject::GetAssayData(
    object2,
    assay = assay2,
    layer = "scale.data"
  )
  check_features <- get_namespace_fun("Seurat", "CheckFeatures")
  features <- check_features(
    data.use = data.use1,
    features = features,
    object.name = "object1",
    verbose = FALSE
  )
  features <- check_features(
    data.use = data.use2,
    features = features,
    object.name = "object2",
    verbose = FALSE
  )
  data1 <- data.use1[features, ]
  data2 <- data.use2[features, ]

  cca.results <- RunCCA.default(
    object1 = data1,
    object2 = data2,
    standardize = TRUE,
    num.cc = num.cc,
    verbose = FALSE
  )

  combined.object <- merge(x = object1, y = object2, merge.data = FALSE, ...)
  rownames(cca.results$ccv) <- SeuratObject::Cells(combined.object)
  colnames(data1) <- SeuratObject::Cells(combined.object)[1:ncol(data1)]
  colnames(data2) <- SeuratObject::Cells(combined.object)[
    (ncol(data1) + 1):length(SeuratObject::Cells(combined.object))
  ]
  combined.object[["cca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = cca.results$ccv[colnames(combined.object), ],
    assay = assay1,
    key = "CC_"
  )
  combined.object[["cca"]]@assay.used <- SeuratObject::DefaultAssay(
    combined.object
  )
  combined.object <- SeuratObject::SetAssayData(
    combined.object,
    new.data = cbind(data1, data2),
    layer = "scale.data"
  )
  combined.object
}

#' Run canonical correlation analysis
#'
#' @param object1 First object or matrix.
#' @param object2 Second object or matrix.
#' @param ... Passed to methods.
#'
#' @return CCA results or an object containing CCA results.
#' @export
RunCCA <- function(object1, object2, ...) {
  UseMethod("RunCCA", object1)
}
