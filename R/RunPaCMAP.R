#' @title Run PaCMAP (Pairwise Controlled Manifold Approximation)
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams RunUMAP2
#' @inheritParams RunDM
#' @param n_components The number of PaCMAP components.
#' Default is `2`.
#' @param n.neighbors A number of neighbors considered in the k-Nearest Neighbor graph.
#' Default is `10` for dataset whose sample size is smaller than 10000.
#' For large dataset whose sample size (n) is larger than 10000, the default value is: `10 + 15 * (log10(n) - 4)`.
#' @param MN_ratio The ratio of the ratio of the number of mid-near pairs to the number of neighbors.
#' Default is `0.5`.
#' @param FP_ratio The ratio of the ratio of the number of further pairs to the number of neighbors.
#' Default is `2`.
#' @param distance_method The distance metric to be used.
#' Default is `"euclidean"`.
#' @param lr The learning rate of the Adam optimizer.
#' Default is `1`.
#' @param num_iters The number of iterations for PaCMAP optimization.
#' Default is `450`.
#' @param apply_pca Whether pacmap should apply PCA to the data before constructing the k-Nearest Neighbor graph.
#' Using PCA to preprocess the data can largely accelerate the DR process without losing too much accuracy.
#' Notice that this option does not affect the initialization of the optimization process.
#' Default is `TRUE`.
#' @param init The initialization of the lower dimensional embedding.
#' One of `"pca"` or `"random"`.
#' Default is `"random"`.
#' @param reduction.name The name of the reduction to be stored in the Seurat object.
#' Default is `"pacmap"`.
#' @param reduction.key The prefix for the column names of the PaCMAP embeddings.
#' Default is `"PaCMAP_"`.
#' @param backend PaCMAP backend. `"cpp"` uses a compiled pair sampler and
#' Adam optimizer; `"python"` retains the official pacmap package.
#' @param ... Additional arguments to be passed to pacmap.PaCMAP.
#'
#' @rdname RunPaCMAP
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#' pancreas_sub <- RunPaCMAP(
#'   object = pancreas_sub,
#'   features = SeuratObject::VariableFeatures(pancreas_sub)
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "pacmap"
#' )
#' }
RunPaCMAP <- function(object, ...) {
  UseMethod(generic = "RunPaCMAP", object = object)
}

#' @rdname RunPaCMAP
#' @method RunPaCMAP Seurat
#' @export
RunPaCMAP.Seurat <- function(
  object,
  reduction = "pca",
  dims = NULL,
  features = NULL,
  assay = NULL,
  layer = "data",
  n_components = 2,
  n.neighbors = NULL,
  MN_ratio = 0.5,
  FP_ratio = 2,
  distance_method = "euclidean",
  lr = 1,
  num_iters = 450L,
  apply_pca = TRUE,
  init = "random",
  reduction.name = "pacmap",
  reduction.key = "PaCMAP_",
  verbose = TRUE,
  seed.use = 11L,
  backend = c("cpp", "python"),
  ...
) {
  if (sum(c(is.null(dims), is.null(features))) < 1) {
    log_message(
      "Please specify only one of the following arguments: dims, features, or graph",
      message_type = "error"
    )
  }
  if (!is.null(features)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- as_matrix(
      Matrix::t(
        GetAssayData5(
          object = object,
          layer = layer,
          assay = assay
        )[features, ]
      )
    )
    if (ncol(data.use) < n_components) {
      log_message(
        "Please provide as many or more features than n_components: ",
        length(features),
        " features provided, ",
        n_components,
        " PaCMAP components requested",
        message_type = "error"
      )
    }
  } else if (!is.null(dims)) {
    reduction <- DefaultReduction(object, pattern = reduction)
    data.use <- Embeddings(object[[reduction]])[, dims]
    assay <- DefaultAssay(object = object[[reduction]])
    if (length(dims) < n_components) {
      log_message(
        "Please provide as many or more dims than n_components: ",
        length(dims),
        " dims provided, ",
        n_components,
        " PaCMAP components requested",
        message_type = "error"
      )
    }
  } else {
    log_message(
      "Please specify one of dims or features",
      message_type = "error"
    )
  }
  object[[reduction.name]] <- RunPaCMAP(
    object = data.use,
    assay = assay,
    n_components = n_components,
    n.neighbors = n.neighbors,
    MN_ratio = MN_ratio,
    FP_ratio = FP_ratio,
    distance_method = distance_method,
    lr = lr,
    num_iters = num_iters,
    apply_pca = apply_pca,
    init = init,
    reduction.key = reduction.key,
    verbose = verbose,
    seed.use = seed.use,
    backend = backend,
    ...
  )
  object <- Seurat::LogSeuratCommand(object = object)
  return(object)
}

#' @rdname RunPaCMAP
#' @method RunPaCMAP default
#' @export
RunPaCMAP.default <- function(
  object,
  assay = NULL,
  n_components = 2,
  n.neighbors = NULL,
  MN_ratio = 0.5,
  FP_ratio = 2,
  distance_method = "euclidean",
  lr = 1,
  num_iters = 450L,
  apply_pca = TRUE,
  init = "random",
  reduction.key = "PaCMAP_",
  verbose = TRUE,
  seed.use = 11L,
  backend = c("cpp", "python"),
  ...
) {
  backend <- match.arg(backend)
  if (!is.null(seed.use)) {
    set.seed(seed = seed.use)
  }

  if (identical(backend, "python")) {
    PrepareEnv(modules = "pacmap")
    check_python("pacmap", verbose = verbose)
    pacmap <- reticulate::import("pacmap")

    operator <- pacmap$PaCMAP(
      n_components = as.integer(n_components),
      n_neighbors = n.neighbors,
      MN_ratio = MN_ratio,
      FP_ratio = FP_ratio,
      distance = distance_method,
      lr = lr,
      num_iters = num_iters,
      apply_pca = apply_pca,
      verbose = verbose,
      random_state = as.integer(seed.use)
    )
    embedding <- operator$fit_transform(object, init = init)
  } else {
    dots <- list(...)
    if (length(dots) > 0L) {
      cli::cli_warn(
        "Additional PaCMAP arguments are used only by {.arg backend = 'python'}."
      )
    }
    data <- native_pacmap_prepare_data(object, apply_pca = apply_pca)
    n.neighbors <- n.neighbors %||% if (nrow(data) <= 10000L) {
      10L
    } else {
      as.integer(round(10 + 15 * (log10(nrow(data)) - 4)))
    }
    n.neighbors <- min(as.integer(n.neighbors), nrow(data) - 1L)
    metric_id <- native_manifold_metric_id(distance_method)
    knn <- manifold_exact_knn_cpp(
      data,
      k = n.neighbors + 1L,
      metric = metric_id
    )
    initial <- native_manifold_initialization(
      data,
      n_components = as.integer(n_components),
      init = match.arg(init, c("pca", "random")),
      seed = seed.use %||% 0L
    )
    embedding <- pacmap_optimize_cpp(
      data = data,
      initial = initial,
      knn_index = knn[["idx"]],
      n_mid = max(1L, as.integer(round(n.neighbors * MN_ratio))),
      n_far = max(1L, as.integer(round(n.neighbors * FP_ratio))),
      learning_rate = as.numeric(lr),
      iterations = if (length(num_iters) == 1L) {
        as.integer(200L + num_iters)
      } else {
        as.integer(sum(num_iters))
      },
      seed = as.integer(seed.use %||% 0L),
      metric = metric_id
    )
  }

  colnames(x = embedding) <- paste0(reduction.key, seq_len(ncol(embedding)))
  rownames(x = embedding) <- rownames(object)

  reduction <- SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    key = reduction.key,
    assay = assay,
    global = TRUE
  )
  return(reduction)
}
