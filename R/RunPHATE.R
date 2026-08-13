#' @title Run PHATE (Potential of Heat-diffusion for Affinity-based Trajectory Embedding)
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams RunUMAP2
#' @inheritParams RunDM
#' @param n_components The number of PHATE components.
#' Default is `2`.
#' @param knn A number of nearest neighbors on which to build kernel.
#' Default is `5`.
#' @param decay The sets decay rate of kernel tails.
#' Default is `40`.
#' @param n_landmark A number of landmarks to use in fast PHATE.
#' Default is `2000`.
#' @param t The power to which the diffusion operator is powered.
#' This sets the level of diffusion.
#' If `"auto"`, `t` is selected according to the knee point in the Von Neumann Entropy of the diffusion operator.
#' Default is `"auto"`.
#' @param gamma The informational distance constant between `-1` and `1`.
#' `gamma=1` gives the PHATE log potential, `gamma=0` gives a square root potential.
#' Default is `1`.
#' @param n_pca A number of principal components to use for calculating neighborhoods.
#' For extremely large datasets, using `n_pca < 20` allows neighborhoods to be calculated in roughly `log(n_samples)` time.
#' Default is `100`.
#' @param knn_dist The distance metric for k-nearest neighbors.
#' Recommended values: `"euclidean"`, `"cosine"`, `"precomputed"`.
#' Default is `"euclidean"`.
#' @param knn_max The maximum number of neighbors for which alpha decaying kernel is computed for each point.
#' For very large datasets, setting `knn_max` to a small multiple of `knn` can speed up computation significantly.
#' Default is `NULL`.
#' @param t_max The maximum `t` to test.
#' Default is `100`.
#' @param do_cluster Whether to perform clustering on the PHATE embeddings.
#' Default is `FALSE`.
#' @param backend PHATE backend. `"python"` calls the upstream `phate` Python
#' package and `"cpp"` uses the native C++ helper path.
#' @param mds MDS algorithm passed to PHATE. The native C++ backend currently
#' implements the `"classic"` path.
#' @param mds_solver Metric MDS solver passed to the Python backend.
#' @param n_clusters A number of clusters to be identified.
#' Default is `"auto"`.
#' @param max_clusters The maximum number of clusters to test.
#' Default is `100`.
#' @param reduction.name The name of the reduction to be stored in the Seurat object.
#' Default is `"phate"`.
#' @param reduction.key The prefix for the column names of the PHATE embeddings.
#' Default is `"PHATE_"`.
#' @param ... Additional arguments to be passed to phate.PHATE.
#'
#' @rdname RunPHATE
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#' pancreas_sub <- RunPHATE(
#'   object = pancreas_sub,
#'   features = SeuratObject::VariableFeatures(pancreas_sub)
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "phate"
#' )
#' }
RunPHATE <- function(object, ...) {
  UseMethod(generic = "RunPHATE", object = object)
}
#' @rdname RunPHATE
#' @method RunPHATE Seurat
#' @export
RunPHATE.Seurat <- function(
  object,
  reduction = "pca",
  dims = NULL,
  features = NULL,
  assay = NULL,
  layer = "data",
  n_components = 2,
  knn = 5,
  decay = 40,
  n_landmark = 2000,
  t = "auto",
  gamma = 1,
  n_pca = 100,
  knn_dist = "euclidean",
  knn_max = NULL,
  t_max = 100,
  do_cluster = FALSE,
  backend = c("python", "cpp"),
  mds = "metric",
  mds_solver = "sgd",
  n_clusters = "auto",
  max_clusters = 100,
  reduction.name = "phate",
  reduction.key = "PHATE_",
  verbose = TRUE,
  seed.use = 11,
  ...
) {
  if (sum(c(is.null(dims), is.null(features))) == 2) {
    log_message(
      "Please specify only one of the following arguments: dims, features",
      message_type = "error"
    )
  }
  if (!is.null(features)) {
    assay <- assay %||% DefaultAssay(object = object)
    data_use <- as_matrix(
      Matrix::t(
        GetAssayData5(
          object = object,
          layer = layer,
          assay = assay
        )[features, ]
      )
    )
    if (ncol(data_use) < n_components) {
      log_message(
        "Please provide as many or more features than n_components: ",
        length(features),
        " features provided, ",
        n_components,
        " PHATE components requested",
        message_type = "error"
      )
    }
  } else if (!is.null(dims)) {
    data_use <- Embeddings(object[[reduction]])[, dims]
    assay <- DefaultAssay(object = object[[reduction]])
    if (length(dims) < n_components) {
      log_message(
        "Please provide as many or more dims than n_components: ",
        length(dims),
        " dims provided, ",
        n_components,
        " PHATE components requested",
        message_type = "error"
      )
    }
  } else {
    log_message(
      "Please specify one of dims or features",
      message_type = "error"
    )
  }
  backend <- match.arg(backend)
  if (identical(backend, "cpp") && missing(mds)) {
    mds <- "classic"
  }
  object[[reduction.name]] <- RunPHATE(
    object = data_use,
    assay = assay,
    n_components = n_components,
    knn = knn,
    decay = decay,
    n_landmark = n_landmark,
    t = t,
    gamma = gamma,
    n_pca = n_pca,
    knn_dist = knn_dist,
    knn_max = knn_max,
    t_max = t_max,
    do_cluster = do_cluster,
    backend = backend,
    mds = mds,
    mds_solver = mds_solver,
    n_clusters = n_clusters,
    max_clusters = max_clusters,
    reduction.key = reduction.key,
    verbose = verbose,
    seed.use = seed.use,
    ...
  )
  object <- Seurat::LogSeuratCommand(object = object)
  return(object)
}

#' @rdname RunPHATE
#' @method RunPHATE default
#' @export
RunPHATE.default <- function(
  object,
  assay = NULL,
  n_components = 2,
  knn = 5,
  decay = 40,
  n_landmark = 2000,
  t = "auto",
  gamma = 1,
  n_pca = 100,
  knn_dist = "euclidean",
  knn_max = NULL,
  t_max = 100,
  do_cluster = FALSE,
  backend = c("python", "cpp"),
  mds = "metric",
  mds_solver = "sgd",
  n_clusters = "auto",
  max_clusters = 100,
  reduction.key = "PHATE_",
  verbose = TRUE,
  seed.use = 11,
  ...
) {
  set.seed(seed = seed.use)
  backend <- match.arg(backend)
  if (identical(backend, "cpp") && missing(mds)) {
    mds <- "classic"
  }
  if (identical(backend, "cpp")) {
    return(run_phate_cpp_reduction(
      object = object,
      assay = assay,
      n_components = n_components,
      knn = knn,
      decay = decay,
      n_landmark = n_landmark,
      t = t,
      gamma = gamma,
      n_pca = n_pca,
      knn_dist = knn_dist,
      knn_max = knn_max,
      t_max = t_max,
      do_cluster = do_cluster,
      mds = mds,
      reduction.key = reduction.key
    ))
  }

  PrepareEnv(modules = "phate")
  check_python("phate", verbose = verbose)
  phate <- reticulate::import("phate")

  if (is.numeric(knn_max) && length(knn_max) > 0) {
    knn_max <- as.integer(knn_max)
  } else {
    knn_max <- NULL
  }
  n_pca_arg <- if (is.null(n_pca)) NULL else as.integer(n_pca)
  operator <- do.call(
    phate$PHATE,
    c(
      list(
        n_components = as.integer(n_components),
        knn = as.integer(knn),
        decay = as.integer(decay),
        n_landmark = as.integer(n_landmark),
        t = if (identical(t, "auto")) "auto" else as.integer(t),
        gamma = as.numeric(gamma),
        n_pca = n_pca_arg,
        knn_dist = as.character(knn_dist),
        knn_max = knn_max,
        mds = as.character(mds),
        mds_solver = as.character(mds_solver),
        random_state = as.integer(seed.use),
        verbose = as.integer(verbose)
      ),
      list(...)
    )
  )
  embedding <- operator$fit_transform(object, t_max = as.integer(t_max))
  colnames(embedding) <- paste0(reduction.key, seq_len(ncol(embedding)))
  rownames(embedding) <- rownames(object)

  reduction <- Seurat::CreateDimReducObject(
    embeddings = embedding,
    key = reduction.key,
    assay = assay,
    global = TRUE
  )
  if (isTRUE(do_cluster)) {
    if (is.numeric(n_clusters)) {
      n_clusters <- as.integer(n_clusters)
    }
    if (is.numeric(max_clusters)) {
      max_clusters <- as.integer(max_clusters)
    }
    clusters <- phate$cluster$kmeans(
      operator,
      n_clusters = n_clusters,
      max_clusters = max_clusters,
      random_state = as.integer(seed.use)
    )
    clusters <- clusters + 1
    clusters <- factor(clusters, levels = sort(unique(clusters)))
    names(clusters) <- rownames(embedding)
    SeuratObject::Misc(reduction, slot = "clusters") <- clusters
  }
  return(reduction)
}

run_phate_cpp_reduction <- function(
  object,
  assay = NULL,
  n_components = 2,
  knn = 5,
  decay = 40,
  n_landmark = 2000,
  t = "auto",
  gamma = 1,
  n_pca = 100,
  knn_dist = "euclidean",
  knn_max = NULL,
  t_max = 100,
  do_cluster = FALSE,
  mds = "classic",
  reduction.key = "PHATE_"
) {
  if (!identical(mds, "classic")) {
    log_message(
      "{.fn RunPHATE} cpp backend currently supports {.arg mds = 'classic'} only",
      message_type = "error"
    )
  }
  if (!identical(knn_dist, "euclidean")) {
    log_message(
      "{.fn RunPHATE} cpp backend currently supports {.val euclidean} distance only",
      message_type = "error"
    )
  }
  if (!identical(as.numeric(gamma), 1)) {
    log_message(
      "{.fn RunPHATE} cpp backend currently supports {.arg gamma = 1} only",
      message_type = "error"
    )
  }
  if (isTRUE(do_cluster)) {
    log_message(
      "{.fn RunPHATE} cpp backend currently does not support {.arg do_cluster = TRUE}",
      message_type = "error"
    )
  }

  data_use <- as.matrix(object)
  if (!is.numeric(data_use)) {
    storage.mode(data_use) <- "double"
  }
  if (is.null(rownames(data_use))) {
    rownames(data_use) <- paste0("cell_", seq_len(nrow(data_use)))
  }
  cell_names <- rownames(data_use)
  n_pca_used <- NULL
  if (!is.null(n_pca)) {
    n_pca <- max(1L, min(as.integer(n_pca), ncol(data_use), nrow(data_use) - 1L))
    if (n_pca < ncol(data_use) && nrow(data_use) > 1L) {
      check_r("irlba", verbose = FALSE)
      data_use <- irlba::prcomp_irlba(
        data_use,
        n = n_pca,
        center = TRUE,
        scale. = FALSE
      )$x
      rownames(data_use) <- cell_names
      n_pca_used <- n_pca
      invisible(gc(full = TRUE))
    }
  }
  n_cells <- nrow(data_use)
  if (n_cells < 1) {
    log_message("PHATE input must contain at least one cell", message_type = "error")
  }
  n_components <- as.integer(n_components)
  knn <- min(as.integer(knn), max(1L, n_cells - 1L))
  knn_search <- if (is.null(knn_max)) {
    knn
  } else {
    min(max(knn, as.integer(knn_max)), max(1L, n_cells - 1L))
  }
  t_max <- as.integer(t_max)
  minimum_landmarks <- if (n_cells > 1L) 2L else 1L
  n_landmark <- max(
    minimum_landmarks,
    min(as.integer(n_landmark), n_cells)
  )

  if (n_cells > n_landmark) {
    landmark_index <- sort(sample.int(n_cells, n_landmark, replace = FALSE))
    landmark_reduction <- run_phate_cpp_reduction(
      object = data_use[landmark_index, , drop = FALSE],
      assay = assay,
      n_components = n_components,
      knn = min(knn, n_landmark - 1L),
      decay = decay,
      n_landmark = n_landmark,
      t = t,
      gamma = gamma,
      n_pca = NULL,
      knn_dist = knn_dist,
      knn_max = knn_max,
      t_max = t_max,
      do_cluster = FALSE,
      mds = mds,
      reduction.key = reduction.key
    )
    landmark_embedding <- SeuratObject::Embeddings(landmark_reduction)
    projection_k <- max(1L, min(knn_search, n_landmark))
    check_r("BiocNeighbors", verbose = FALSE)
    projection <- BiocNeighbors::queryKNN(
      X = data_use[landmark_index, , drop = FALSE],
      query = data_use,
      k = projection_k,
      BNPARAM = BiocNeighbors::HnswParam(
        ef.search = max(50L, 3L * projection_k),
        distance = "Euclidean"
      ),
      num.threads = 1L
    )
    bandwidth_column <- min(knn, ncol(projection$distance))
    bandwidth <- pmax(
      projection$distance[, bandwidth_column],
      .Machine$double.eps
    )
    weights <- exp(-(
      projection$distance / bandwidth
    )^as.numeric(decay))
    weights[!is.finite(weights)] <- 0
    weight_totals <- rowSums(weights)
    zero_weight <- weight_totals <= 0
    if (any(zero_weight)) {
      weights[zero_weight, ] <- 0
      weights[cbind(which(zero_weight), 1L)] <- 1
      weight_totals[zero_weight] <- 1
    }
    weights <- weights / weight_totals
    embedding <- matrix(0, n_cells, n_components)
    for (neighbor in seq_len(ncol(projection$index))) {
      embedding <- embedding +
        landmark_embedding[projection$index[, neighbor], , drop = FALSE] *
          weights[, neighbor]
    }
    embedding[landmark_index, ] <- landmark_embedding
    colnames(embedding) <- paste0(reduction.key, seq_len(ncol(embedding)))
    rownames(embedding) <- cell_names
    return(Seurat::CreateDimReducObject(
      embeddings = embedding,
      key = reduction.key,
      assay = assay,
      global = TRUE,
      misc = list(
        backend = "cpp",
        t = SeuratObject::Misc(landmark_reduction, slot = "t"),
        landmark = TRUE,
        n_landmarks = n_landmark,
        landmark_cells = cell_names[landmark_index],
        projection_k = projection_k,
        n_pca = n_pca_used
      )
    ))
  }

  affinity <- phate_graphtools_affinity_data_cpp(
    data = data_use,
    knn = knn,
    decay = as.numeric(decay),
    thresh = 1e-4,
    knn_max = if (is.null(knn_max)) -1L else as.integer(knn_max)
  )
  diffusion_t <- if (identical(t, "auto")) {
    phate_find_optimal_t_cpp(
      rows = affinity$affinity_rows,
      cols = affinity$affinity_cols,
      vals = affinity$affinity_vals,
      n_cells = affinity$n_cells,
      t_max = t_max
    )
  } else {
    as.integer(t)
  }
  log_transition <- phate_diffusion_operator_cpp(
    rows = affinity$affinity_rows,
    cols = affinity$affinity_cols,
    vals = affinity$affinity_vals,
    n_cells = affinity$n_cells,
    t_max = diffusion_t
  )
  distance <- phate_potential_distance_cpp(
    log_transition = log_transition,
    n_landmarks = as.integer(n_landmark)
  )
  embedding <- phate_metric_mds_cpp(
    D = distance,
    n_components = n_components
  )
  colnames(embedding) <- paste0(reduction.key, seq_len(ncol(embedding)))
  rownames(embedding) <- rownames(data_use)

  Seurat::CreateDimReducObject(
    embeddings = embedding,
    key = reduction.key,
    assay = assay,
    global = TRUE,
    misc = list(
      backend = "cpp",
      t = diffusion_t,
      landmark = FALSE,
      n_landmarks = n_cells,
      n_pca = n_pca_used
    )
  )
}
