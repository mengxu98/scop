#' @title Run UMAP (Uniform Manifold Approximation and Projection)
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param object An object. This can be a Seurat object, a matrix-like object, a Neighbor object, or a Graph object.
#' @param reduction The reduction to be used. Default is `"pca"`.
#' @param dims The dimensions to be used. Default is `NULL`.
#' @param features The features to be used. Default is `NULL`.
#' @param neighbor The name of the Neighbor object to be used. Default is `NULL`.
#' @param graph The name of the Graph object to be used. Default is `NULL`.
#' @param assay The assay to be used. Default is `NULL`.
#' @param layer The layer to be used. Default is `"data"`.
#' @param umap.method The UMAP method to be used.
#' Options are `"naive"` and `"uwot"`.
#' Default is `"uwot"`.
#' @param n_threads Num of threads used.
#' @param reduction.model A DimReduc object containing a pre-trained UMAP model.
#' Default is `NULL`.
#' @param return.model Whether to return the UMAP model. Default is `FALSE`.
#' @param n.neighbors A number of nearest neighbors to be used. Default is `30`.
#' @param n.components A number of UMAP components. Default is `2`.
#' @param metric The metric or a function to be used for distance calculations.
#' When using a string, available metrics are: `euclidean`, `manhattan`.
#' Other available generalized metrics are: cosine, pearson, pearson2.
#' Note the triangle inequality may not be satisfied by some generalized metrics, hence knn search may not be optimal.
#' When using metric.function as a function,
#' the signature must be function(matrix, origin, target) and should compute a distance between the origin column and the target columns.
#' Default is `"cosine"`.
#' @param n.epochs A number of iterations performed during layout optimization for UMAP.
#' Default is `200`.
#' @param spread The spread parameter for UMAP, used during automatic estimation of a/b parameters.
#' Default is `1`.
#' @param min.dist The minimum distance between UMAP embeddings, determines how close points appear in the final layout.
#' Default is `0.3`.
#' @param set.op.mix.ratio Interpolate between (fuzzy) union and intersection as the set operation used to combine local fuzzy simplicial sets to obtain a global fuzzy simplicial sets.
#' Both fuzzy set operations use the product t-norm.
#' The value of this parameter should be between `0.0` and `1.0`;
#' a value of `1.0` will use a pure fuzzy union, while `0.0` will use a pure fuzzy intersection.
#' @param local.connectivity The local connectivity, used during construction of fuzzy simplicial set.
#' Default is `1`.
#' @param negative.sample.rate The negative sample rate for UMAP optimization.
#' Determines how many non-neighbor points are used per point and per iteration during layout optimization.
#' Default is `5`.
#' @param a The parameter a for UMAP optimization.
#' Contributes to gradient calculations during layout optimization.
#' When left at NA, a suitable value will be estimated automatically.
#' Default is `NULL`.
#' @param b The parameter b for UMAP optimization. Details see parameter `a`.
#' @param learning.rate The initial value of "learning rate" of layout optimization.
#' Default is `1`.
#' @param repulsion.strength A numeric value determines, together with alpha, the learning rate of layout optimization.
#' Default is `1`.
#' @param reduction.name The name of the reduction to be stored in the Seurat object.
#' Default is `"umap"`.
#' @param reduction.key The prefix for the column names of the UMAP embeddings.
#' Default is `"UMAP_"`.
#' @param seed.use The random seed to be used.
#' Default is `11`.
#' @param ... Additional arguments to be passed to UMAP.
#'
#' @rdname RunUMAP2
#' @export
#'
#' @examples
#' pancreas_sub <- Seurat::FindVariableFeatures(pancreas_sub)
#' pancreas_sub <- RunUMAP2(
#'   object = pancreas_sub,
#'   features = SeuratObject::VariableFeatures(pancreas_sub)
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "umap"
#' )
RunUMAP2 <- function(object, ...) {
  UseMethod(generic = "RunUMAP2", object = object)
}

#' @rdname RunUMAP2
#' @method RunUMAP2 Seurat
#' @export
RunUMAP2.Seurat <- function(
    object,
    reduction = "pca",
    dims = NULL,
    features = NULL,
    neighbor = NULL,
    graph = NULL,
    assay = NULL,
    layer = "data",
    umap.method = "uwot",
    reduction.model = NULL,
    n_threads = NULL,
    return.model = FALSE,
    n.neighbors = 30L,
    n.components = 2L,
    metric = "cosine",
    n.epochs = 200L,
    spread = 1,
    min.dist = 0.3,
    set.op.mix.ratio = 1,
    local.connectivity = 1L,
    negative.sample.rate = 5L,
    a = NULL,
    b = NULL,
    learning.rate = 1,
    repulsion.strength = 1,
    reduction.name = "umap",
    reduction.key = "UMAP_",
    verbose = TRUE,
    seed.use = 11,
    ...) {
  if (
    sum(c(
      is.null(x = dims),
      is.null(x = features),
      is.null(neighbor),
      is.null(x = graph)
    )) ==
      4
  ) {
    log_message(
      "Please specify only one of the following arguments: dims, features, neighbor or graph",
      message_type = "error"
    )
  }
  if (!is.null(x = features)) {
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
    if (ncol(x = data.use) < n.components) {
      log_message(
        "Please provide as many or more features than n.components: ",
        length(x = features),
        " features provided, ",
        n.components,
        " UMAP components requested",
        message_type = "error"
      )
    }
  } else if (!is.null(x = dims)) {
    data.use <- as_matrix(Embeddings(object[[reduction]])[, dims])
    assay <- DefaultAssay(object = object[[reduction]])
    if (length(x = dims) < n.components) {
      log_message(
        "Please provide as many or more dims than n.components: ",
        length(x = dims),
        " dims provided, ",
        n.components,
        " UMAP components requested",
        message_type = "error"
      )
    }
  } else if (!is.null(x = neighbor)) {
    if (!inherits(x = object[[neighbor]], what = "Neighbor")) {
      log_message(
        "Please specify a Neighbor object name, ",
        "instead of the name of a ",
        class(object[[neighbor]]),
        " object",
        message_type = "error"
      )
    }
    data.use <- object[[neighbor]]
  } else if (!is.null(x = graph)) {
    if (!inherits(x = object[[graph]], what = "Graph")) {
      log_message(
        "Please specify a Graph object name, ",
        "instead of the name of a ",
        class(object[[graph]]),
        " object",
        message_type = "error"
      )
    }
    data.use <- object[[graph]]
  } else {
    log_message(
      "Please specify one of dims, features, neighbor, or graph",
      message_type = "error"
    )
  }
  object[[reduction.name]] <- RunUMAP2(
    object = data.use,
    assay = assay,
    umap.method = umap.method,
    reduction.model = reduction.model,
    return.model = return.model,
    n.neighbors = n.neighbors,
    n.components = n.components,
    metric = metric,
    n.epochs = n.epochs,
    spread = spread,
    min.dist = min.dist,
    set.op.mix.ratio = set.op.mix.ratio,
    local.connectivity = local.connectivity,
    negative.sample.rate = negative.sample.rate,
    a = a,
    b = b,
    learning.rate = learning.rate,
    repulsion.strength = repulsion.strength,
    seed.use = seed.use,
    verbose = verbose,
    reduction.key = reduction.key
  )
  object <- Seurat::LogSeuratCommand(object = object)
  return(object)
}

#' @rdname RunUMAP2
#' @method RunUMAP2 default
#' @export
RunUMAP2.default <- function(
    object,
    assay = NULL,
    umap.method = "uwot",
    reduction.model = NULL,
    n_threads = NULL,
    return.model = FALSE,
    n.neighbors = 30L,
    n.components = 2L,
    metric = "cosine",
    n.epochs = 200L,
    spread = 1,
    min.dist = 0.3,
    set.op.mix.ratio = 1,
    local.connectivity = 1L,
    negative.sample.rate = 5L,
    a = NULL,
    b = NULL,
    learning.rate = 1,
    repulsion.strength = 1,
    reduction.key = "UMAP_",
    verbose = TRUE,
    seed.use = 11L,
    ...) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (return.model) {
    log_message(
      "UMAP will return its model",
      verbose = verbose
    )
  }
  if (!is.null(x = reduction.model)) {
    log_message(
      "Running UMAP projection",
      verbose = verbose
    )
    if (
      is.null(x = reduction.model) ||
        !inherits(x = reduction.model, what = "DimReduc")
    ) {
      log_message(
        "If running projection UMAP, please pass a DimReduc object with the model stored to reduction.model",
        message_type = "error"
      )
    }
    model <- SeuratObject::Misc(object = reduction.model, slot = "model")
    if (length(x = model) == 0) {
      log_message(
        "The provided reduction.model does not have a model stored",
        message_type = "error"
      )
    }
    umap.method <- ifelse(
      "layout" %in% names(model),
      "naive-predict",
      "uwot-predict"
    )
  }
  n.epochs <- as.integer(n.epochs)
  n.neighbors <- as.integer(n.neighbors)
  n.components <- as.integer(n.components)
  local.connectivity <- as.integer(local.connectivity)
  negative.sample.rate <- as.integer(negative.sample.rate)

  if (inherits(x = object, what = "Neighbor")) {
    object <- list(
      idx = SeuratObject::Indices(object),
      dist = SeuratObject::Distances(object)
    )
  }

  if (umap.method == "naive") {
    umap.config <- umap::umap.defaults
    umap.config$n_neighbors <- n.neighbors
    umap.config$n_components <- n.components
    umap.config$metric <- metric
    umap.config$n_epochs <- ifelse(is.null(n.epochs), 200, n.epochs)
    umap.config$spread <- spread
    umap.config$min_dist <- min.dist
    umap.config$set_op_mix_ratio <- set.op.mix.ratio
    umap.config$local_connectivity <- local.connectivity
    umap.config$a <- ifelse(is.null(a), NA, a)
    umap.config$b <- ifelse(is.null(b), NA, b)
    umap.config$gamma <- repulsion.strength
    umap.config$alpha <- learning.rate
    umap.config$negative_sample_rate <- negative.sample.rate
    umap.config$random_state <- seed.use
    umap.config$transform_state <- seed.use
    umap.config$verbose <- verbose
    if (is.na(umap.config$a) || is.na(umap.config$b)) {
      umap.config[c("a", "b")] <- get_namespace_fun(
        "umap", "find.ab.params"
      )(
        umap.config$spread,
        umap.config$min_dist
      )
      umap.config$min_dist <- umap::umap.defaults$min_dist
    }

    if (inherits(x = object, what = "dist")) {
      knn <- get_namespace_fun(
        "umap", "knn.from.dist"
      )(d = object, k = n.neighbors)
      out <- umap::umap(
        d = matrix(nrow = nrow(object)),
        config = umap.config,
        knn = knn
      )
      embeddings <- out$layout
      rownames(x = embeddings) <- attr(object, "Labels")
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- SeuratObject::CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      if (return.model) {
        SeuratObject::Misc(reduction, slot = "model") <- out
      }
      return(reduction)
    }
    if (inherits(x = object, what = "list")) {
      knn <- umap::umap.knn(
        indexes = object[["idx"]],
        distances = object[["dist"]]
      )
      out <- umap::umap(
        d = matrix(nrow = nrow(object[["idx"]])),
        config = umap.config,
        knn = knn
      )
      embeddings <- out$layout
      rownames(x = embeddings) <- rownames(object[["idx"]])
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- SeuratObject::CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      if (return.model) {
        SeuratObject::Misc(reduction, slot = "model") <- out
      }
      return(reduction)
    }
    if (inherits(x = object, what = "Graph")) {
      if (!inherits(object, "dgCMatrix")) {
        object <- SeuratObject::as.sparse(object[1:nrow(object), ])
      }
      diag(object) <- 0
      if (ncol(object) > 10000) {
        obs_sample <- sample(1:ncol(object), size = 10000)
      } else {
        obs_sample <- 1:ncol(object)
      }
      if (!isSymmetric(as_matrix(object[obs_sample, obs_sample]))) {
        log_message(
          "Graph must be a symmetric matrix.",
          message_type = "error"
        )
      }

      coo <- matrix(ncol = 3, nrow = length(object@x))
      coo[, 1] <- rep(1:ncol(object), diff(object@p))
      coo[, 2] <- object@i + 1
      coo[, 3] <- object@x
      colnames(coo) <- c("from", "to", "value")
      coo <- coo[order(coo[, 1], coo[, 2]), ]
      graph <- list(
        coo = coo,
        names = rownames(object),
        n.elements = nrow(object)
      )
      class(graph) <- "coo"

      umap.config$init <- "spectral"
      initial <- get_namespace_fun(
        "umap", "make.initial.embedding"
      )(
        V = graph$n.elements,
        config = umap.config,
        g = graph
      )
      embeddings <- get_namespace_fun(
        "umap", "naive.simplicial.set.embedding"
      )(
        g = graph,
        embedding = initial,
        config = umap.config
      )
      embeddings <- get_namespace_fun(
        "umap", "center.embedding"
      )(embeddings)
      rownames(x = embeddings) <- rownames(x = object)
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- SeuratObject::CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      if (return.model) {
        log_message(
          "{.arg return.model} does not support 'Graph' input",
          message_type = "warning"
        )
      }
      return(reduction)
    }
    if (
      inherits(x = object, what = "matrix") ||
        inherits(x = object, what = "Matrix")
    ) {
      out <- umap::umap(d = object, config = umap.config, method = "naive")
      embeddings <- out$layout
      rownames(x = embeddings) <- rownames(object)
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- SeuratObject::CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      if (return.model) {
        SeuratObject::Misc(reduction, slot = "model") <- out
      }
      return(reduction)
    }
  }

  if (umap.method == "uwot") {
    if (inherits(x = object, what = "dist")) {
      embeddings <- uwot::umap(
        X = object,
        n_neighbors = n.neighbors,
        n_threads = n_threads,
        n_components = n.components,
        metric = metric,
        n_epochs = n.epochs,
        learning_rate = learning.rate,
        min_dist = min.dist,
        spread = spread,
        set_op_mix_ratio = set.op.mix.ratio,
        local_connectivity = local.connectivity,
        repulsion_strength = repulsion.strength,
        negative_sample_rate = negative.sample.rate,
        a = a,
        b = b,
        ret_model = FALSE
      )
      rownames(x = embeddings) <- attr(object, "Labels")
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- SeuratObject::CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      if (return.model) {
        log_message(
          "{.arg return.model} does not support 'dist' input",
          message_type = "warning"
        )
      }
      return(reduction)
    }
    if (inherits(x = object, what = "list")) {
      out <- uwot::umap(
        X = NULL,
        nn_method = object,
        n_threads = n_threads,
        n_components = n.components,
        metric = metric,
        n_epochs = n.epochs,
        learning_rate = learning.rate,
        min_dist = min.dist,
        spread = spread,
        set_op_mix_ratio = set.op.mix.ratio,
        local_connectivity = local.connectivity,
        repulsion_strength = repulsion.strength,
        negative_sample_rate = negative.sample.rate,
        a = a,
        b = b,
        ret_model = return.model
      )
      if (return.model) {
        embeddings <- out$embedding
      } else {
        embeddings <- out
      }
      rownames(x = embeddings) <- row.names(object[["idx"]])
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- SeuratObject::CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      if (return.model) {
        out$nn_index <- NULL
        SeuratObject::Misc(reduction, slot = "model") <- out
      }
      return(reduction)
    }
    if (inherits(x = object, what = "Graph")) {
      if (!inherits(object, "dgCMatrix")) {
        object <- SeuratObject::as.sparse(object[1:nrow(object), ])
      }
      diag(object) <- 0
      if (ncol(object) > 10000) {
        obs_sample <- sample(1:ncol(object), size = 10000)
      } else {
        obs_sample <- 1:ncol(object)
      }
      if (!isSymmetric(as_matrix(object[obs_sample, obs_sample]))) {
        log_message(
          "Graph must be a symmetric matrix.",
          message_type = "error"
        )
      }
      val <- split(object@x, rep(1:ncol(object), diff(object@p)))
      pos <- split(object@i + 1, rep(1:ncol(object), diff(object@p)))
      idx <- Matrix::t(mapply(
        function(x, y) {
          out <- y[utils::head(order(x, decreasing = TRUE), n.neighbors)]
          length(out) <- n.neighbors
          return(out)
        },
        x = val,
        y = pos
      ))
      connectivity <- Matrix::t(mapply(
        function(x, y) {
          out <- y[utils::head(order(x, decreasing = TRUE), n.neighbors)]
          length(out) <- n.neighbors
          out[is.na(out)] <- 0
          return(out)
        },
        x = val,
        y = val
      ))
      idx[is.na(idx)] <- sample(
        1:nrow(object),
        size = sum(is.na(idx)),
        replace = TRUE
      )
      nn <- list(
        idx = idx,
        dist = max(connectivity) -
          connectivity +
          min(diff(range(connectivity)), 1) / 1e50
      )
      out <- uwot::umap(
        X = NULL,
        nn_method = nn,
        n_threads = n_threads,
        n_components = n.components,
        metric = metric,
        n_epochs = n.epochs,
        learning_rate = learning.rate,
        min_dist = min.dist,
        spread = spread,
        set_op_mix_ratio = set.op.mix.ratio,
        local_connectivity = local.connectivity,
        repulsion_strength = repulsion.strength,
        negative_sample_rate = negative.sample.rate,
        a = a,
        b = b,
        ret_model = return.model
      )
      if (return.model) {
        embeddings <- out$embedding
      } else {
        embeddings <- out
      }

      rownames(x = embeddings) <- row.names(object)
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- SeuratObject::CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      if (return.model) {
        out$nn_index <- NULL
        SeuratObject::Misc(reduction, slot = "model") <- out
      }
      return(reduction)
    }
    if (
      inherits(x = object, what = "matrix") ||
        inherits(x = object, what = "Matrix")
    ) {
      out <- uwot::umap(
        X = object,
        n_neighbors = n.neighbors,
        n_threads = n_threads,
        n_components = n.components,
        metric = metric,
        n_epochs = n.epochs,
        learning_rate = learning.rate,
        min_dist = min.dist,
        spread = spread,
        set_op_mix_ratio = set.op.mix.ratio,
        local_connectivity = local.connectivity,
        repulsion_strength = repulsion.strength,
        negative_sample_rate = negative.sample.rate,
        a = a,
        b = b,
        ret_model = return.model
      )
      if (return.model) {
        embeddings <- out$embedding
      } else {
        embeddings <- out
      }
      rownames(x = embeddings) <- row.names(object)
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- SeuratObject::CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      if (return.model) {
        out$nn_index <- NULL
        SeuratObject::Misc(reduction, slot = "model") <- out
      }
      return(reduction)
    }
  }

  if (umap.method == "naive-predict") {
    if (
      inherits(x = object, what = "matrix") ||
        inherits(x = object, what = "Matrix")
    ) {
      class(model) <- "umap"
      if (any(!colnames(model$data) %in% colnames(object))) {
        log_message(
          "query data must contain the same features with the model:\n",
          paste(utils::head(colnames(model$data), 10), collapse = ","),
          " ......",
          message_type = "error"
        )
      }
      embeddings <- stats::predict(model, object[, colnames(model$data)])
      rownames(x = embeddings) <- row.names(object)
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- SeuratObject::CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      return(reduction)
    } else {
      log_message(
        "naive umap model only support 'matrix' input.",
        message_type = "error"
      )
    }
  }

  if (umap.method == "uwot-predict") {
    if (inherits(x = object, what = "list")) {
      if (ncol(object[["idx"]]) != model$n_neighbors) {
        log_message(
          "Number of neighbors between query and reference is not equal to the number of neighbros within reference",
          message_type = "warning"
        )
        model$n_neighbors <- ncol(object[["idx"]])
      }
      if (
        is.null(model$num_precomputed_nns) || model$num_precomputed_nns == 0
      ) {
        model$num_precomputed_nns <- 1
      }
      embeddings <- uwot::umap_transform(
        X = NULL,
        nn_method = object,
        model = model,
        n_epochs = n.epochs,
        n_threads = n_threads,
        verbose = FALSE
      )
      rownames(x = embeddings) <- row.names(object[["idx"]])
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- SeuratObject::CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      return(reduction)
    }
    if (inherits(x = object, what = "Graph")) {
      match_k <- Matrix::t(as_matrix(apply(
        object,
        2,
        function(x) order(x, decreasing = TRUE)[1:n.neighbors]
      )))
      match_k_connectivity <- Matrix::t(as_matrix(apply(
        object,
        2,
        function(x) x[order(x, decreasing = TRUE)[1:n.neighbors]]
      )))
      object <- list(
        idx = match_k,
        dist = max(match_k_connectivity) - match_k_connectivity
      )
      if (
        is.null(model$num_precomputed_nns) || model$num_precomputed_nns == 0
      ) {
        model$num_precomputed_nns <- 1
      }
      embeddings <- uwot::umap_transform(
        X = NULL,
        nn_method = object,
        model = model,
        n_epochs = n.epochs,
        n_threads = n_threads,
        verbose = FALSE
      )
      rownames(x = embeddings) <- row.names(object[["idx"]])
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- SeuratObject::CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      return(reduction)
    }
    if (
      inherits(x = object, what = "matrix") ||
        inherits(x = object, what = "Matrix")
    ) {
      embeddings <- uwot::umap_transform(
        X = object,
        model = model,
        n_epochs = n.epochs,
        n_threads = n_threads,
        verbose = FALSE
      )
      rownames(x = embeddings) <- row.names(object)
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- SeuratObject::CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      return(reduction)
    }
  }
}
