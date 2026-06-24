runumap_equal_scalar <- function(x, y) {
  length(x) == 1L && !is.na(x) && identical(as.character(x), as.character(y))
}

runumap_default_epochs <- function(n_cells, n_epochs) {
  if (!is.null(n_epochs)) {
    return(as.integer(n_epochs))
  }
  if (n_cells <= 10000L) 500L else 200L
}

runumap_supported <- function(
  extra,
  dims,
  features,
  graph,
  nn.name,
  reduction,
  umap.method,
  reduction.model,
  return.model,
  n.neighbors,
  n.components,
  metric,
  learning.rate,
  repulsion.strength,
  negative.sample.rate,
  uwot.sgd,
  uwot.approx_pow,
  metric.kwds,
  angular.rp.forest,
  densmap
) {
  length(extra) == 0L &&
    !is.null(dims) &&
    is.null(features) &&
    is.null(graph) &&
    is.null(nn.name) &&
    length(reduction) == 1L &&
    !is.na(reduction) &&
    runumap_equal_scalar(umap.method, "uwot") &&
    is.null(reduction.model) &&
    identical(return.model, FALSE) &&
    as.integer(n.neighbors) >= 2L &&
    as.integer(n.components) == 2L &&
    runumap_equal_scalar(metric, "cosine") &&
    isTRUE(all.equal(as.numeric(learning.rate), 1)) &&
    isTRUE(all.equal(as.numeric(repulsion.strength), 1)) &&
    as.integer(negative.sample.rate) == 5L &&
    identical(uwot.sgd, FALSE) &&
    identical(uwot.approx_pow, FALSE) &&
    is.null(metric.kwds) &&
    identical(angular.rp.forest, FALSE) &&
    identical(densmap, FALSE)
}

runumap_embedding <- function(
  X,
  n.neighbors,
  local.connectivity,
  set.op.mix.ratio,
  seed.use,
  n.epochs,
  n.components,
  spread,
  min.dist,
  a,
  b,
  repulsion.strength,
  learning.rate,
  negative.sample.rate,
  cores
) {
  epochs <- runumap_default_epochs(nrow(X), n.epochs)
  if (is.null(a) || is.null(b)) {
    find_ab_params <- get("find_ab_params", envir = asNamespace("uwot"))
    ab <- find_ab_params(spread = spread, min_dist = min.dist)
    if (is.null(a)) {
      a <- unname(ab[["a"]])
    }
    if (is.null(b)) {
      b <- unname(ab[["b"]])
    }
  }
  if (!is.null(seed.use)) {
    set.seed(seed.use)
  }
  use_umap2 <- exists("umap2", envir = asNamespace("uwot"), inherits = FALSE)
  umap_fun <- if (use_umap2) {
    get("umap2", envir = asNamespace("uwot"))
  } else {
    get("umap", envir = asNamespace("uwot"))
  }
  args <- list(
    X = X,
    n_threads = if (is.null(cores)) NULL else as.integer(cores),
    n_neighbors = as.integer(n.neighbors),
    n_trees = 50L,
    search_k = 2L * as.integer(n.neighbors) * 50L,
    n_components = as.integer(n.components),
    metric = "cosine",
    n_epochs = as.integer(epochs),
    learning_rate = learning.rate,
    min_dist = min.dist,
    spread = spread,
    set_op_mix_ratio = set.op.mix.ratio,
    local_connectivity = local.connectivity,
    repulsion_strength = repulsion.strength,
    negative_sample_rate = negative.sample.rate,
    a = a,
    b = b,
    fast_sgd = FALSE,
    approx_pow = FALSE,
    verbose = FALSE,
    ret_model = FALSE
  )
  if (use_umap2) {
    args$seed <- seed.use
  }
  do.call(umap_fun, args)
}

runumap_store <- function(
  object,
  reduction.name,
  reduction.key,
  assay,
  emb
) {
  key <- reduction.key %||%
    SeuratObject::Key(object = reduction.name, quiet = TRUE)
  colnames(emb) <- paste0(key, seq_len(ncol(emb)))
  object[[reduction.name]] <- Seurat::CreateDimReducObject(
    embeddings = emb,
    key = key,
    assay = assay,
    global = TRUE
  )
  SeuratObject::LogSeuratCommand(object = object)
}

#' @export
RunUMAP.Seurat <- function(
  object,
  dims = NULL,
  reduction = "pca",
  features = NULL,
  graph = NULL,
  assay = SeuratObject::DefaultAssay(object = object),
  nn.name = NULL,
  slot = "data",
  umap.method = "uwot",
  reduction.model = NULL,
  return.model = FALSE,
  n.neighbors = 30L,
  n.components = 2L,
  metric = "cosine",
  n.epochs = NULL,
  learning.rate = 1,
  min.dist = 0.3,
  spread = 1,
  set.op.mix.ratio = 1,
  local.connectivity = 1L,
  repulsion.strength = 1,
  negative.sample.rate = 5L,
  a = NULL,
  b = NULL,
  uwot.sgd = FALSE,
  uwot.approx_pow = FALSE,
  seed.use = 42L,
  metric.kwds = NULL,
  angular.rp.forest = FALSE,
  densmap = FALSE,
  dens.lambda = 2,
  dens.frac = 0.3,
  dens.var.shift = 0.1,
  verbose = TRUE,
  reduction.name = "umap",
  reduction.key = NULL,
  cores = NULL,
  ...
) {
  extra <- list(...)
  if (!is.null(cores)) {
    cores <- suppressWarnings(as.integer(cores[[1L]]))
    if (is.na(cores) || cores < 1L) {
      stop("cores must be a positive integer.", call. = FALSE)
    }
    cores <- as.integer(min(cores, 16L))
  }
  if (
    !runumap_supported(
      extra = extra,
      dims = dims,
      features = features,
      graph = graph,
      nn.name = nn.name,
      reduction = reduction,
      umap.method = umap.method,
      reduction.model = reduction.model,
      return.model = return.model,
      n.neighbors = n.neighbors,
      n.components = n.components,
      metric = metric,
      learning.rate = learning.rate,
      repulsion.strength = repulsion.strength,
      negative.sample.rate = negative.sample.rate,
      uwot.sgd = uwot.sgd,
      uwot.approx_pow = uwot.approx_pow,
      metric.kwds = metric.kwds,
      angular.rp.forest = angular.rp.forest,
      densmap = densmap
    )
  ) {
    seurat_runumap <- get("RunUMAP.Seurat", envir = asNamespace("Seurat"))
    return(seurat_runumap(
      object = object,
      dims = dims,
      reduction = reduction,
      features = features,
      graph = graph,
      assay = assay,
      nn.name = nn.name,
      slot = slot,
      umap.method = umap.method,
      reduction.model = reduction.model,
      return.model = return.model,
      n.neighbors = n.neighbors,
      n.components = n.components,
      metric = metric,
      n.epochs = n.epochs,
      learning.rate = learning.rate,
      min.dist = min.dist,
      spread = spread,
      set.op.mix.ratio = set.op.mix.ratio,
      local.connectivity = local.connectivity,
      repulsion.strength = repulsion.strength,
      negative.sample.rate = negative.sample.rate,
      a = a,
      b = b,
      uwot.sgd = uwot.sgd,
      uwot.approx_pow = uwot.approx_pow,
      seed.use = seed.use,
      metric.kwds = metric.kwds,
      angular.rp.forest = angular.rp.forest,
      densmap = densmap,
      dens.lambda = dens.lambda,
      dens.frac = dens.frac,
      dens.var.shift = dens.var.shift,
      verbose = verbose,
      reduction.name = reduction.name,
      reduction.key = reduction.key,
      ...
    ))
  }

  emb_src <- SeuratObject::Embeddings(object[[reduction]])
  if (identical(as.integer(dims), seq_len(ncol(emb_src)))) {
    X <- emb_src
  } else {
    X <- emb_src[, dims, drop = FALSE]
  }
  if (ncol(X) < as.integer(n.components)) {
    stop("insufficient dimensions for UMAP", call. = FALSE)
  }
  assay <- SeuratObject::DefaultAssay(object = object[[reduction]])
  emb <- runumap_embedding(
    X = X,
    n.neighbors = n.neighbors,
    local.connectivity = local.connectivity,
    set.op.mix.ratio = set.op.mix.ratio,
    seed.use = seed.use,
    n.epochs = n.epochs,
    n.components = n.components,
    spread = spread,
    min.dist = min.dist,
    a = a,
    b = b,
    repulsion.strength = repulsion.strength,
    learning.rate = learning.rate,
    negative.sample.rate = negative.sample.rate,
    cores = cores
  )
  rownames(emb) <- rownames(X)
  runumap_store(object, reduction.name, reduction.key, assay, emb)
}

#' Run UMAP
#'
#' @param object Object containing reduced-dimensional data.
#' @param ... Passed to methods.
#'
#' @return The input object with a UMAP reduction added.
#' @export
RunUMAP <- function(object, ...) {
  UseMethod("RunUMAP")
}
