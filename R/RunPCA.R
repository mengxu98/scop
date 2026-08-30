pca_native_auto_supported <- function(object, npcs, extra_args) {
  is.matrix(object) &&
    is.numeric(object) &&
    length(extra_args) == 0L &&
    nrow(object) >= 2L &&
    nrow(object) <= 5000L &&
    ncol(object) >= 8L * nrow(object) &&
    npcs >= 1L &&
    npcs + 1L < min(nrow(object), ncol(object))
}

pca_native_candidate_acceptable <- function(eigvals, npcs, min_gap = 0.005) {
  length(eigvals) >= npcs + 1L &&
    all(is.finite(eigvals[seq_len(npcs + 1L)])) &&
    eigvals[[npcs]] > 0 &&
    (eigvals[[npcs]] - eigvals[[npcs + 1L]]) / eigvals[[npcs]] >= min_gap
}

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
    log_message("RunPCA supports Seurat, assay, matrix, and Matrix objects.", message_type = "error")
  }
  dots <- list(...)
  backend <- dots[["backend"]] %||% "auto"
  backend <- match.arg(backend, c("auto", "irlba", "cpp"))
  dots[["backend"]] <- NULL
  logical_flags <- list(rev.pca, weight.by.var, approx)
  if (
    !all(vapply(
      logical_flags,
      function(x) is.logical(x) && length(x) == 1L && !is.na(x),
      logical(1)
    )) ||
      !is.numeric(npcs) || length(npcs) != 1L || !is.finite(npcs) ||
      npcs < 1 || npcs != as.integer(npcs) ||
      isTRUE(rev.pca) || !isTRUE(approx) ||
      (inherits(object, "Matrix") && !identical(backend, "cpp"))
  ) {
    return(do.call(
      utils::getFromNamespace("RunPCA.default", "Seurat"),
      c(
        list(
          object = object,
          assay = assay,
          npcs = npcs,
          rev.pca = rev.pca,
          weight.by.var = weight.by.var,
          verbose = verbose,
          ndims.print = ndims.print,
          nfeatures.print = nfeatures.print,
          reduction.key = reduction.key,
          seed.use = seed.use,
          approx = approx
        ),
        dots
      )
    ))
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
  nv <- NULL
  native_auto <- identical(backend, "auto") &&
    pca_native_auto_supported(obj, npcs, dots)
  if (identical(backend, "cpp") || native_auto) {
    native_npcs <- if (native_auto) npcs + 1L else npcs
    nv <- tryCatch(
      pca_backend_run(obj, as.integer(native_npcs), isTRUE(weight.by.var)),
      error = function(e) NULL
    )
    if (
      native_auto &&
        !is.null(nv) &&
        !pca_native_candidate_acceptable(nv$eigvals, npcs)
    ) {
      nv <- NULL
    }
  }
  if (!is.null(nv)) {
    keep <- seq_len(npcs)
    feature.loadings <- nv$loadings[, keep, drop = FALSE]
    cell.embeddings <- nv$embeddings[, keep, drop = FALSE]
    sdev <- as.numeric(nv$sdev[keep])
  } else {
    pca <- do.call(
      irlba::irlba,
      c(list(A = Matrix::t(obj), nv = npcs), dots)
    )
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
    assay.used = assay %||% "RNA",
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
    log_message("RunPCA.StdAssay supports layer = 'scale.data' and approx = TRUE.", message_type = "error")
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
  approx <- extra[["approx"]] %||% TRUE
  if (
    !all(vapply(
      list(rev.pca, weight.by.var, approx),
      function(x) is.logical(x) && length(x) == 1L && !is.na(x),
      logical(1)
    )) ||
      isTRUE(rev.pca) ||
      isFALSE(approx) ||
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
