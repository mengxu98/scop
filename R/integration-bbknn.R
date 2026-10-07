#' @title BBKNN integration function
#'
#' @inheritParams RunIntegration
#' @param bbknn_params A list of parameters for the bbknn.matrix.bbknn function, default is an empty list.
#' @param backend BBKNN graph backend. `"cpp"` uses the compiled cross-batch
#' KNN graph implementation; `"python"` retains the official bbknn package.
#'
#' @export
BBKNN_integrate <- function(
  srt_merge = NULL,
  batch = NULL,
  append = TRUE,
  srt_list = NULL,
  assay = NULL,
  do_normalization = NULL,
  normalization_method = "LogNormalize",
  do_HVF_finding = TRUE,
  HVF_source = "separate",
  HVF_method = "vst",
  nHVF = 2000,
  HVF_min_intersection = 1,
  HVF = NULL,
  do_scaling = TRUE,
  vars_to_regress = NULL,
  regression_model = "linear",
  scale_within_batch = FALSE,
  linear_reduction = "pca",
  linear_reduction_dims = 50,
  linear_reduction_dims_use = NULL,
  linear_reduction_params = list(),
  force_linear_reduction = FALSE,
  nonlinear_reduction = "umap",
  nonlinear_reduction_dims = c(2, 3),
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  bbknn_params = list(),
  verbose = TRUE,
  seed = 11,
  backend = c("cpp", "python")
) {
  backend <- match.arg(backend)
  if (identical(backend, "python")) {
    PrepareEnv(modules = "bbknn")
    check_python("bbknn")
    bbknn <- reticulate::import("bbknn")
  }

  if (length(linear_reduction) > 1) {
    log_message(
      "Only the first method in the {.arg linear_reduction} will be used",
      message_type = "warning"
    )
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c(
    "pca",
    "svd",
    "ica",
    "nmf",
    "mds",
    "glmpca"
  )
  if (!is.null(srt_merge)) {
    reduc_test <- c(reduc_test, SeuratObject::Reductions(srt_merge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    log_message(
      "{.arg linear_reduction} must be one of {.val {reduc_test}}",
      message_type = "error"
    )
  }
  if (
    !is.null(linear_reduction_dims_use) &&
      max(linear_reduction_dims_use) > linear_reduction_dims
  ) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  validate_nonlinear_reductions(
    nonlinear_reduction,
    allowed = c("umap", "umap-naive", "fr")
  )
  cluster_algorithm_index <- resolve_cluster_algorithm_index(cluster_algorithm)

  set.seed(seed)

  validate_integration_input_cells(srt_list, srt_merge)
  srt_merge_raw <- srt_merge
  if (!is.null(srt_list)) {
    checked <- CheckDataList(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srt_merge <- Reduce(merge, srt_list)
    SeuratObject::VariableFeatures(srt_merge) <- HVF
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    checked <- CheckDataMerge(
      srt_merge = srt_merge,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      verbose = verbose,
      seed = seed
    )
    srt_merge <- checked[["srt_merge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    log_message(
      "{.arg normalization_method} is {.pkg TFIDF}. Use {.pkg lsi} workflow..."
    )
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }
  scale_features <- rownames(
    GetAssayData5(
      srt_merge,
      layer = "scale.data",
      assay = SeuratObject::DefaultAssay(srt_merge)
    )
  )
  if (
    isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))
  ) {
    log_message("Perform {.fn Seurat::ScaleData}")
    assay_merge <- SeuratObject::DefaultAssay(srt_merge)
    if (inherits(srt_merge[[assay_merge]], "Assay5")) {
      srt_merge[[assay_merge]] <- SeuratObject::JoinLayers(
        srt_merge[[assay_merge]]
      )
    }
    srt_merge <- ScaleData(
      object = srt_merge,
      split.by = if (isTRUE(scale_within_batch)) batch else NULL,
      assay = SeuratObject::DefaultAssay(srt_merge),
      features = HVF,
      vars.to.regress = vars_to_regress,
      model.use = regression_model,
      verbose = FALSE
    )
  }

  log_message(
    "Perform {.val {linear_reduction}} linear dimension reduction"
  )
  srt_merge <- RunDimsReduction(
    srt_merge,
    prefix = "BBKNN",
    features = HVF,
    assay = SeuratObject::DefaultAssay(srt_merge),
    linear_reduction = linear_reduction,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = force_linear_reduction,
    verbose = verbose,
    seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- resolve_linear_dims_use(
      srt = srt_merge,
      reduction = paste0("BBKNN", linear_reduction),
      normalization_method = normalization_method,
      reduction_method = linear_reduction
    )
  }

  log_message(
    "Perform {.pkg BBKNN} integration with the {.val {backend}} backend",
    verbose = verbose
  )
  log_message(
    "Using {.val {paste0('BBKNN', linear_reduction)}} ({.val {min(linear_reduction_dims_use)}}:{.val {max(linear_reduction_dims_use)}}) as input",
    verbose = verbose
  )
  emb <- Embeddings(srt_merge, reduction = paste0("BBKNN", linear_reduction))[,
    linear_reduction_dims_use,
    drop = FALSE
  ]
  if (identical(backend, "python")) {
    params <- list(
      pca = emb,
      batch_list = srt_merge[[batch, drop = TRUE]]
    )
    for (nm in names(bbknn_params)) {
      params[[nm]] <- bbknn_params[[nm]]
    }
    bem <- invoke_fun(bbknn$matrix$bbknn, params)
  } else {
    bem <- bbknn_native_matrix(
      embedding = emb,
      batches = srt_merge[[batch, drop = TRUE]],
      params = bbknn_params
    )
  }
  n.neighbors <- bem[[3]]$n_neighbors
  srt_integrated <- srt_merge

  bbknn_graph <- SeuratObject::as.sparse(
    bem[[2]][1:nrow(bem[[2]]), , drop = FALSE]
  )
  rownames(bbknn_graph) <- colnames(bbknn_graph) <- rownames(emb)
  bbknn_graph <- SeuratObject::as.Graph(bbknn_graph)
  bbknn_graph@assay.used <- SeuratObject::DefaultAssay(srt_integrated)
  srt_integrated@graphs[["BBKNN"]] <- bbknn_graph

  bbknn_dist <- Matrix::t(
    SeuratObject::as.sparse(
      bem[[1]][1:nrow(bem[[1]]), , drop = FALSE]
    )
  )
  rownames(bbknn_dist) <- colnames(bbknn_dist) <- rownames(emb)
  bbknn_dist <- SeuratObject::as.Graph(bbknn_dist)
  bbknn_dist@assay.used <- SeuratObject::DefaultAssay(srt_integrated)
  srt_integrated@graphs[["BBKNN_dist"]] <- bbknn_dist

  val <- split(
    bbknn_dist@x,
    rep(seq_len(ncol(bbknn_dist)), diff(bbknn_dist@p))
  )
  pos <- split(
    bbknn_dist@i + 1,
    rep(seq_len(ncol(bbknn_dist)), diff(bbknn_dist@p))
  )
  idx <- Matrix::t(
    mapply(
      function(x, y) {
        out <- y[utils::head(order(x, decreasing = FALSE), n.neighbors - 1)]
        length(out) <- n.neighbors - 1
        out
      },
      x = val,
      y = pos
    )
  )
  idx[is.na(idx)] <- sample(
    seq_len(nrow(idx)),
    size = sum(is.na(idx)),
    replace = TRUE
  )
  idx <- cbind(seq_len(nrow(idx)), idx)
  dist <- Matrix::t(mapply(
    function(x, y) {
      out <- y[utils::head(order(x, decreasing = FALSE), n.neighbors - 1)]
      length(out) <- n.neighbors - 1
      out[is.na(out)] <- 0
      out
    },
    x = val,
    y = val
  ))
  dist <- cbind(0, dist)
  srt_integrated[["BBKNN_neighbors"]] <- methods::new(
    Class = "Neighbor",
    nn.idx = idx,
    nn.dist = dist,
    alg.info = list(),
    cell.names = rownames(emb)
  )
  nonlinear_reduction_params[["n.neighbors"]] <- n.neighbors

  srt_integrated <- find_neighbors_and_clusters(
    srt = srt_integrated,
    reduction = NULL,
    dims_use = NULL,
    graph_prefix = "BBKNN_",
    graph_snn = "BBKNN",
    cluster_colname = "BBKNNclusters",
    HVF = HVF,
    neighbor_metric = "euclidean",
    neighbor_k = 20L,
    cluster_algorithm = cluster_algorithm,
    cluster_algorithm_index = cluster_algorithm_index,
    cluster_resolution = cluster_resolution,
    run_find_neighbors = FALSE,
    verbose = verbose
  )

  srt_integrated <- run_nonlinear_reduction(
    srt = srt_integrated,
    prefix = "BBKNN",
    graph_use = "BBKNN",
    neighbor_use = "BBKNN_neighbors",
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed,
    verbose = verbose
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[[
    "BBKNN_HVF"
  ]] <- HVF
  srt_integrated@misc[["BBKNN_backend"]] <- backend
  srt_integrated@misc[["BBKNN_parameters"]] <- bem[[3]]

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_integrated,
      pattern = paste0(assay, "|BBKNN|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_integrated)
  }
}


bbknn_seurat_annoy_cross <- function(reference, query, k, n_trees, metric) {
  seurat_metric <- switch(metric,
    angular = "cosine",
    euclidean = "euclidean",
    manhattan = "manhattan",
    metric
  )
  if (is.null(rownames(reference))) {
    rownames(reference) <- paste0("r", seq_len(nrow(reference)))
  }
  if (is.null(rownames(query))) {
    rownames(query) <- paste0("q", seq_len(nrow(query)))
  }
  result <- get_namespace_fun("Seurat", "AnnoyNN")(
    data = reference,
    query = query,
    metric = seurat_metric,
    n.trees = n_trees,
    k = k,
    include.distance = TRUE
  )
  list(idx = result$nn.idx, distance = result$nn.dists)
}


bbknn_native_matrix <- function(embedding, batches, params = list()) {
  embedding <- as.matrix(embedding)
  batches <- as.factor(batches)
  if (anyNA(batches)) {
    log_message("{.arg batch} contains missing values", message_type = "error")
  }
  neighbors_within_batch <- as.integer(
    params[["neighbors_within_batch"]] %||% 3L
  )
  n_pcs <- min(
    ncol(embedding),
    as.integer(params[["n_pcs"]] %||% 50L)
  )
  if (n_pcs < 1L) {
    log_message("{.arg n_pcs} must be positive", message_type = "error")
  }
  embedding <- embedding[, seq_len(n_pcs), drop = FALSE]
  computation <- tolower(params[["computation"]] %||% "annoy")
  if (computation %in% c("ckdtree", "kdtree", "faiss")) {
    computation <- "exact"
  }
  if (!computation %in% c("annoy", "exact")) {
    log_message(
      "Native BBKNN supports {.val annoy} and {.val exact} computation",
      message_type = "error"
    )
  }
  metric <- tolower(params[["metric"]] %||% "euclidean")
  supported_metrics <- if (identical(computation, "annoy")) {
    c("euclidean", "angular", "manhattan")
  } else {
    c("euclidean", "cosine", "angular")
  }
  if (!metric %in% supported_metrics) {
    log_message(
      "Unsupported native BBKNN metric {.val {metric}} for {.val {computation}}",
      message_type = "error"
    )
  }
  annoy_n_trees <- as.integer(params[["annoy_n_trees"]] %||% 10L)
  if (annoy_n_trees < 1L) {
    log_message(
      "{.arg annoy_n_trees} must be positive",
      message_type = "error"
    )
  }
  cores <- as.integer(
    params[["cores"]] %||%
      thisutils::detect_cores(max_threads = 8L)
  )
  batch_indices <- split(seq_len(nrow(embedding)), batches, drop = TRUE)
  if (any(lengths(batch_indices) < neighbors_within_batch)) {
    log_message(
      "Each batch must contain at least {.val {neighbors_within_batch}} cells",
      message_type = "error"
    )
  }

  neighbor_blocks <- lapply(batch_indices, function(reference_index) {
    result <- if (identical(computation, "annoy")) {
      bbknn_seurat_annoy_cross(
        reference = embedding[reference_index, , drop = FALSE],
        query = embedding,
        k = neighbors_within_batch,
        n_trees = annoy_n_trees,
        metric = metric
      )
    } else {
      exact_metric <- if (identical(metric, "angular")) "cosine" else metric
      cross_knn_f32(
        reference = embedding[reference_index, , drop = FALSE],
        query = embedding,
        k = neighbors_within_batch,
        metric = exact_metric,
        cores = cores
      )
    }
    list(
      idx = matrix(
        reference_index[result[["idx"]]],
        nrow = nrow(embedding),
        ncol = neighbors_within_batch
      ),
      dist = result[["distance"]]
    )
  })
  idx <- do.call(cbind, lapply(neighbor_blocks, `[[`, "idx"))
  dist <- do.call(cbind, lapply(neighbor_blocks, `[[`, "dist"))
  fuzzy <- bbknn_fuzzy_membership_cpp(
    index = idx,
    distance = dist,
    local_connectivity = as.numeric(params[["local_connectivity"]] %||% 1),
    bandwidth = as.numeric(params[["bandwidth"]] %||% 1)
  )
  idx <- fuzzy[["idx"]]
  dist <- fuzzy[["distance"]]
  membership <- fuzzy[["membership"]]
  n_neighbors <- ncol(idx)
  rows <- rep(seq_len(nrow(idx)), times = ncol(idx))
  cols <- as.integer(idx)
  distance_values <- as.numeric(dist)
  nonself <- rows != cols & is.finite(distance_values)
  rows <- rows[nonself]
  cols <- cols[nonself]
  distance_values <- distance_values[nonself]

  connectivity_values <- as.numeric(membership)[nonself]
  distance_graph <- Matrix::sparseMatrix(
    i = rows,
    j = cols,
    x = distance_values,
    dims = c(nrow(embedding), nrow(embedding))
  )
  connectivity <- Matrix::sparseMatrix(
    i = rows,
    j = cols,
    x = connectivity_values,
    dims = c(nrow(embedding), nrow(embedding))
  )
  connectivity_transpose <- Matrix::t(connectivity)
  connectivity_product <- connectivity * connectivity_transpose
  set_op_mix_ratio <- as.numeric(params[["set_op_mix_ratio"]] %||% 1)
  if (!is.finite(set_op_mix_ratio) ||
    set_op_mix_ratio < 0 ||
    set_op_mix_ratio > 1) {
    log_message(
      "{.arg set_op_mix_ratio} must be between 0 and 1",
      message_type = "error"
    )
  }
  connectivity_union <-
    connectivity + connectivity_transpose - connectivity_product
  connectivity <- set_op_mix_ratio * connectivity_union +
    (1 - set_op_mix_ratio) * connectivity_product
  connectivity <- Matrix::drop0(connectivity)
  trim <- params[["trim"]]
  if (is.null(trim)) {
    trim <- 10L * n_neighbors
  }
  trim <- as.integer(trim)
  if (trim < 0L) {
    log_message("{.arg trim} must be non-negative", message_type = "error")
  }
  if (trim > 0L) {
    connectivity <- methods::as(connectivity, "dgCMatrix")
    connectivity_summary <- list(
      i = connectivity@i + 1L,
      j = rep.int(
        seq_len(ncol(connectivity)),
        diff(connectivity@p)
      ),
      x = connectivity@x
    )
    thresholds <- numeric(nrow(connectivity))
    split_values <- split(
      connectivity_summary$x,
      factor(connectivity_summary$i, levels = seq_len(nrow(connectivity)))
    )
    thresholds <- vapply(split_values, function(values) {
      if (length(values) <= trim) {
        0
      } else {
        sort(values, partial = length(values) - trim + 1L)[
          length(values) - trim + 1L
        ]
      }
    }, numeric(1))
    keep <- connectivity_summary$x >= thresholds[connectivity_summary$i] &
      connectivity_summary$x >= thresholds[connectivity_summary$j]
    connectivity <- Matrix::sparseMatrix(
      i = connectivity_summary$i[keep],
      j = connectivity_summary$j[keep],
      x = connectivity_summary$x[keep],
      dims = dim(connectivity)
    )
  }
  rownames(connectivity) <- colnames(connectivity) <- rownames(embedding)
  rownames(distance_graph) <- colnames(distance_graph) <- rownames(embedding)

  list(
    distance_graph,
    connectivity,
    list(
      n_neighbors = n_neighbors,
      neighbors_within_batch = neighbors_within_batch,
      metric = metric,
      n_pcs = n_pcs,
      trim = trim,
      computation = computation,
      annoy_n_trees = if (identical(computation, "annoy")) {
        annoy_n_trees
      } else {
        NULL
      },
      backend = if (identical(computation, "annoy")) {
        "seurat_annoy"
      } else {
        paste0("cpp_", computation)
      }
    )
  )
}
