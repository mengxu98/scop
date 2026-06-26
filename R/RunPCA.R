#' @export
RunPCA.default <- function(
  object,
  assay = NULL,
  npcs = 50,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "PC_",
  seed.use = 42,
  approx = TRUE,
  feature.names = rownames(object),
  cell.names = colnames(object),
  ...
) {
  if (!is.matrix(object) && !inherits(object, "Matrix")) {
    stop("RunPCA supports Seurat, assay, matrix, and Matrix objects.", call. = FALSE)
  }
  if (isTRUE(rev.pca) || !isTRUE(approx)) {
    stop("RunPCA.default supports rev.pca = FALSE and approx = TRUE.", call. = FALSE)
  }
  if (!is.null(seed.use)) {
    set.seed(seed = seed.use)
  }
  npcs <- min(npcs, nrow(object) - 1)
  nfeatures <- nrow(object)
  n <- ncol(object)

  total_variance <- sum(fast_row_vars(object), na.rm = TRUE)
  obj <- object
  if (!is.matrix(obj)) {
    obj <- as.matrix(obj)
  }
  if (!is.double(obj)) {
    storage.mode(obj) <- "double"
  }
  nv <- tryCatch(
    pca_backend_run(obj, as.integer(npcs), isTRUE(weight.by.var)),
    error = function(e) NULL
  )
  if (!is.null(nv)) {
    feature.loadings <- nv$loadings
    cell.embeddings <- nv$embeddings
    sdev <- as.numeric(nv$sdev)
  } else {
    pca <- irlba::irlba(A = Matrix::t(obj), nv = npcs, ...)
    feature.loadings <- pca$v
    if (isTRUE(weight.by.var)) {
      cell.embeddings <- pca$u %*% diag(pca$d, nrow = length(pca$d))
    } else {
      cell.embeddings <- pca$u
    }
    sdev <- pca$d / sqrt(max(1, n - 1))
  }
  rownames(feature.loadings) <- feature.names
  colnames(feature.loadings) <- paste0(reduction.key, 1:npcs)
  rownames(cell.embeddings) <- cell.names
  colnames(cell.embeddings) <- colnames(feature.loadings)
  methods::new(
    "DimReduc",
    cell.embeddings = cell.embeddings,
    feature.loadings = feature.loadings,
    feature.loadings.projected = matrix(nrow = 0L, ncol = 0L),
    assay.used = assay %||% character(0),
    global = FALSE,
    stdev = as.numeric(sdev),
    jackstraw = methods::new("JackStrawData"),
    misc = list(total.variance = total_variance),
    key = reduction.key
  )
}

#' @export
RunPCA.StdAssay <- function(
  object,
  assay = NULL,
  features = NULL,
  layer = "scale.data",
  npcs = 50,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "PC_",
  seed.use = 42,
  ...
) {
  if (
    (!is.null(layer) && !identical(as.character(layer), "scale.data")) ||
      isFALSE(list(...)[["approx"]])
  ) {
    stop("RunPCA.StdAssay supports layer = 'scale.data' and approx = TRUE.", call. = FALSE)
  }
  if (inherits(object, "Assay5")) {
    layer <- SeuratObject::Layers(object, search = layer)
    data.use <- methods::slot(object, "layers")[[layer]]
    feature.names <- SeuratObject::Features(object, layer = layer)
  } else {
    data.use <- methods::slot(object, "scale.data")
    feature.names <- rownames(data.use)
  }
  cell.names <- colnames(object)
  if (is.null(features)) {
    features <- SeuratObject::VariableFeatures(object)
  }
  features <- features[!is.na(features)]
  idx <- match(features, feature.names, nomatch = 0L)
  idx <- idx[idx > 0L]
  data.use <- data.use[idx, , drop = FALSE]
  feature.names <- feature.names[idx]
  RunPCA.default(
    object = data.use,
    assay = assay,
    npcs = npcs,
    rev.pca = rev.pca,
    weight.by.var = weight.by.var,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    feature.names = feature.names,
    cell.names = cell.names,
    ...
  )
}

#' @export
RunPCA.Seurat <- function(
  object,
  assay = NULL,
  features = NULL,
  npcs = 50,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.name = "pca",
  reduction.key = "PC_",
  seed.use = 42,
  ...
) {
  extra <- list(...)
  run_seurat <- function() {
    seurat_runpca <- get("RunPCA.Seurat", envir = asNamespace("Seurat"))
    seurat_runpca(
      object = object,
      assay = assay,
      features = features,
      npcs = npcs,
      rev.pca = rev.pca,
      weight.by.var = weight.by.var,
      verbose = verbose,
      ndims.print = ndims.print,
      nfeatures.print = nfeatures.print,
      reduction.name = reduction.name,
      reduction.key = reduction.key,
      seed.use = seed.use,
      ...
    )
  }
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  if (
      isTRUE(rev.pca) ||
      isFALSE(extra[["approx"]]) ||
      !assay %in% SeuratObject::Assays(object)
  ) {
    return(run_seurat())
  }
  reduction <- tryCatch(
    RunPCA.StdAssay(
      object = object[[assay]],
      assay = assay,
      features = features,
      npcs = npcs,
      rev.pca = rev.pca,
      weight.by.var = weight.by.var,
      verbose = verbose,
      ndims.print = ndims.print,
      nfeatures.print = nfeatures.print,
      reduction.key = reduction.key,
      seed.use = seed.use,
      ...
    ),
    error = function(e) NULL
  )
  if (!is.null(reduction)) {
    object[[reduction.name]] <- reduction
    return(SeuratObject::LogSeuratCommand(object = object))
  }
  run_seurat()
}

#' Run principal component analysis
#'
#' @param object Object containing expression data.
#' @param ... Passed to methods.
#'
#' @return PCA results or an object containing PCA results.
#' @export
RunPCA <- function(object, ...) {
  UseMethod("RunPCA")
}
