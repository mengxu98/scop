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

  obj <- object
  if (!is.matrix(obj)) {
    obj <- as.matrix(obj)
  }
  if (!is.double(obj)) {
    storage.mode(obj) <- "double"
  }
  nv <- pca_backend_run(obj, as.integer(npcs), isTRUE(weight.by.var))
  feature.loadings <- nv$loadings
  cell.embeddings <- nv$embeddings
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
    stdev = as.numeric(nv$sdev),
    jackstraw = methods::new("JackStrawData"),
    misc = list(),
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
  if (!isTRUE(setequal(features, feature.names))) {
    features <- features[!is.na(features)]
    idx <- match(features, feature.names, nomatch = 0L)
    data.use <- data.use[idx[idx > 0L], , drop = FALSE]
    feature.names <- feature.names[idx[idx > 0L]]
  }
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
  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  if (!assay %in% SeuratObject::Assays(object)) {
    stop(sprintf("Assay '%s' is not present in object.", assay), call. = FALSE)
  }
  reduction.data <- RunPCA.StdAssay(
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
  )
  object[[reduction.name]] <- reduction.data
  SeuratObject::LogSeuratCommand(object = object)
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
