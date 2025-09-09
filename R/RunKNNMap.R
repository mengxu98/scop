#' Single-cell reference mapping with KNN method
#'
#' This function performs single-cell reference mapping using the K-nearest neighbor (KNN) method. It takes two single-cell datasets as input: srt_query and srt_ref. The function maps cells from the srt_query dataset to the srt_ref dataset based on their similarity or distance.
#'
#' @param srt_query An object of class Seurat storing the query cells.
#' @param srt_ref An object of class Seurat storing the reference cells.
#' @param query_assay A character string specifying the assay name for the query cells. If not provided, the default assay for the query object will be used.
#' @param ref_assay A character string specifying the assay name for the reference cells. If not provided, the default assay for the reference object will be used.
#' @param ref_umap A character string specifying the name of the UMAP reduction in the reference object. If not provided, the first UMAP reduction found in the reference object will be used.
#' @param ref_group A character string specifying a metadata column name in the reference object to use for grouping.
#' @param features A character vector specifying the features to be used for calculating the distance metric. If not provided, the function will use the variable features calculated by the Seurat package.
#' @param nfeatures A integer specifying the number of highly variable features to be calculated if \code{features} is not provided.
#' @param query_reduction A character string specifying the name of a dimensionality reduction in the query object to use for calculating the distance metric.
#' @param ref_reduction A character string specifying the name of a dimensionality reduction in the reference object to use for calculating the distance metric.
#' @param query_dims A numeric vector specifying the dimension indices from the query reduction to be used for calculating the distance metric.
#' @param ref_dims A numeric vector specifying the dimension indices from the reference reduction to be used for calculating the distance metric.
#' @param projection_method A character string specifying the projection method to use. Options are "model" and "knn". If "model" is selected, the function will try to use a pre-trained UMAP model in the reference object for projection. If "knn" is selected, the function will directly find the nearest neighbors using the distance metric.
#' @param nn_method A character string specifying the nearest neighbor search method to use. Options are "raw", "annoy", and "rann". If "raw" is selected, the function will use the brute-force method to find the nearest neighbors. If "annoy" is selected, the function will use the Annoy library for approximate nearest neighbor search. If "rann" is selected, the function will use the RANN library for approximate nearest neighbor search. If not provided, the function will choose the search method based on the size of the query and reference datasets.
#' @param k An integer specifying the number of nearest neighbors to find for each cell in the query object.
#' @param distance_metric A character string specifying the distance metric to use for calculating the pairwise distances between cells. Options include: "pearson", "spearman", "cosine", "correlation", "jaccard", "ejaccard", "dice", "edice", "hamman", "simple matching", and "faith". Additional distance metrics can also be used, such as "euclidean", "manhattan", "hamming", etc.
#' @param vote_fun A character string specifying the function to be used for aggregating the nearest neighbors in the reference object. Options are "mean", "median", "sum", "min", "max", "sd", "var", etc. If not provided, the default is "mean".
#'
#' @export
#'
#' @examples
#' data(panc8_sub)
#' srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
#' srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
#' srt_ref <- integration_scop(
#'   srt_ref,
#'   batch = "tech",
#'   integration_method = "Seurat"
#' )
#' CellDimPlot(
#'   srt_ref,
#'   group.by = c("celltype", "tech")
#' )
#'
#' # Set the number of threads for RcppParallel
#' # details see: ?RcppParallel::setThreadOptions
#' # if (requireNamespace("RcppParallel", quietly = TRUE)) {
#' #   RcppParallel::setThreadOptions()
#' # }
#' # Projection
#' srt_query <- RunKNNMap(
#'   srt_query = srt_query,
#'   srt_ref = srt_ref,
#'   ref_umap = "SeuratUMAP2D"
#' )
#' ProjectionPlot(
#'   srt_query = srt_query,
#'   srt_ref = srt_ref,
#'   query_group = "celltype", ref_group = "celltype"
#' )
RunKNNMap <- function(
    srt_query,
    srt_ref,
    query_assay = NULL,
    ref_assay = NULL,
    ref_umap = NULL,
    ref_group = NULL,
    features = NULL,
    nfeatures = 2000,
    query_reduction = NULL,
    ref_reduction = NULL,
    query_dims = 1:30,
    ref_dims = 1:30,
    projection_method = c("model", "knn"),
    nn_method = NULL,
    k = 30,
    distance_metric = "cosine",
    vote_fun = "mean") {
  query_assay <- query_assay %||% SeuratObject::DefaultAssay(srt_query)
  ref_assay <- ref_assay %||% SeuratObject::DefaultAssay(srt_ref)
  if (!is.null(ref_group)) {
    if (length(ref_group) == ncol(srt_ref)) {
      srt_ref[["ref_group"]] <- ref_group
    } else if (length(ref_group) == 1) {
      if (!ref_group %in% colnames(srt_ref@meta.data)) {
        log_message(
          "ref_group must be one of the column names in the meta.data",
          message_type = "error"
        )
      } else {
        srt_ref[["ref_group"]] <- srt_ref[[ref_group]]
      }
    } else {
      log_message(
        "Length of ref_group must be one or length of srt_ref.",
        message_type = "error"
      )
    }
    ref_group <- "ref_group"
  }
  if (is.null(ref_umap)) {
    ref_umap <- sort(
      SeuratObject::Reductions(srt_ref)[grep(
        "umap",
        SeuratObject::Reductions(srt_ref),
        ignore.case = TRUE
      )]
    )[1]
    if (length(ref_umap) == 0) {
      log_message(
        "Cannot find UMAP reduction in the srt_ref",
        message_type = "error"
      )
    } else {
      log_message("Set ref_umap to ", ref_umap)
    }
  }
  projection_method <- match.arg(projection_method)
  if (projection_method == "model" && !"model" %in% names(srt_ref[[ref_umap]]@misc)) {
    log_message("No UMAP model detected. Set the projection_method to 'knn'")
    projection_method <- "knn"
  }
  distance_metrics <- c("euclidean", "cosine", "manhattan", "hamming")
  if (projection_method == "model" && !distance_metric %in% distance_metrics) {
    log_message(
      "{.arg distance_metric} must be one of {.val {distance_metrics}} when projection_method = 'model'",
      message_type = "error"
    )
  }
  simil_method <- c(
    "pearson",
    "spearman",
    "cosine",
    "correlation",
    "jaccard",
    "ejaccard",
    "dice",
    "edice",
    "hamman",
    "simple matching",
    "faith"
  )
  dist_method <- c(
    "euclidean",
    "chisquared",
    "kullback",
    "manhattan",
    "maximum",
    "canberra",
    "minkowski",
    "hamming"
  )
  if (!distance_metric %in% c(simil_method, dist_method)) {
    log_message(
      "{.arg distance_metric} must be one of {.val {c(simil_method, dist_method)}}",
      message_type = "error"
    )
  }
  if (projection_method == "model") {
    model <- srt_ref[[ref_umap]]@misc$model
    if ("layout" %in% names(model)) {
      if (k != model$config$n_neighbors) {
        k <- model$config$n_neighbors
        log_message("Set k to {.val k} which is used in the umap model")
      }
    } else if ("embedding" %in% names(model)) {
      if (k != model$n_neighbors) {
        k <- model$n_neighbors
        log_message("Set k to {.val k} which is used in the umap model")
      }
    }
  }

  if (!is.null(query_reduction) && !is.null(ref_reduction)) {
    log_message("Use the reduction to calculate distance metric")
    if (!is.null(query_dims) && !is.null(ref_dims) && length(query_dims) == length(ref_dims)) {
      query <- Seurat::Embeddings(srt_query, reduction = query_reduction)
      query <- query[, query_dims]
      ref <- Seurat::Embeddings(srt_ref, reduction = ref_reduction)
      ref <- ref[, ref_dims]
    } else {
      log_message(
        "{.arg query_dims} and {.arg ref_dims} must be provided with the same length",
        message_type = "error"
      )
    }
  } else {
    log_message("Use the features to calculate distance metric")
    status_query <- CheckDataType(
      data = GetAssayData5(
        srt_query,
        layer = "data",
        assay = query_assay,
        verbose = FALSE
      )
    )
    status_ref <- CheckDataType(
      data = GetAssayData5(
        srt_ref,
        layer = "data",
        assay = ref_assay,
        verbose = FALSE
      )
    )
    if (status_ref != status_query) {
      log_message(
        "Data type is different between {.arg srt_query} and {.arg srt_ref}",
        message_type = "warning"
      )
    }
    if (any(status_query == "unknown", status_ref == "unknown")) {
      log_message(
        "Data type is unknown in {.arg srt_query} or {.arg srt_ref}",
        message_type = "warning"
      )
    }
    if (length(features) == 0) {
      if (length(SeuratObject::VariableFeatures(srt_ref, assay = ref_assay)) == 0) {
        srt_ref <- Seurat::FindVariableFeatures(
          srt_ref,
          nfeatures = nfeatures,
          assay = ref_assay
        )
      }
      if (length(SeuratObject::VariableFeatures(srt_query, assay = query_assay)) == 0) {
        srt_query <- Seurat::FindVariableFeatures(
          srt_query,
          nfeatures = nfeatures,
          assay = query_assay
        )
      }
      features <- intersect(
        SeuratObject::VariableFeatures(srt_query, assay = query_assay),
        SeuratObject::VariableFeatures(srt_ref, assay = ref_assay)
      )
    }
    features_common <- Reduce(
      intersect,
      list(
        features,
        rownames(srt_query[[query_assay]]),
        rownames(srt_ref[[ref_assay]])
      )
    )
    log_message("Use ", length(features_common), " features to calculate distance.")
    query <- Matrix::t(
      GetAssayData5(
        srt_query,
        layer = "data",
        assay = query_assay
      )[
        features_common,
      ]
    )
    ref <- Matrix::t(
      GetAssayData5(
        srt_ref,
        layer = "data",
        assay = ref_assay
      )[
        features_common,
      ]
    )
  }

  if (
    projection_method == "model" &&
      "layout" %in% names(model) &&
      is.null(ref_group)
  ) {
    srt_query[["ref.embeddings"]] <- RunUMAP2(
      object = query,
      reduction.model = srt_ref[[ref_umap]],
      assay = query_assay
    )
    srt_query[["ref.embeddings"]]@misc[["reduction.model"]] <- ref_umap
    return(srt_query)
  }

  if (is.null(nn_method)) {
    if (as.numeric(nrow(query)) * as.numeric(nrow(ref)) >= 1e8) {
      nn_method <- "annoy"
    } else {
      nn_method <- "raw"
    }
  }
  log_message("Use '", nn_method, "' method to find neighbors.")
  if (!nn_method %in% c("raw", "annoy", "rann")) {
    log_message(
      "nn_method must be one of raw, rann and annoy",
      message_type = "error"
    )
  }
  if (
    nn_method == "annoy" &&
      !distance_metric %in% c("euclidean", "cosine", "manhattan", "hamming")
  ) {
    log_message(
      "distance_metric must be one of euclidean, cosine, manhattan, and hamming when nn_method='annoy'",
      message_type = "error"
    )
  }

  if (nn_method %in% c("annoy", "rann")) {
    query.neighbor <- Seurat::FindNeighbors(
      query = query,
      object = ref,
      k.param = k,
      nn.method = nn_method,
      annoy.metric = distance_metric,
      return.neighbor = TRUE
    )
    match_k <- query.neighbor@nn.idx
    rownames(match_k) <- rownames(query)
    match_k_cell <- apply(match_k, c(1, 2), function(x) rownames(ref)[x])
    knn_cells <- c(match_k_cell)
    match_k_distance <- query.neighbor@nn.dist
    rownames(match_k_distance) <- rownames(query)
    refumap_all <- srt_ref[[ref_umap]]@cell.embeddings[knn_cells, ]
    group <- rep(query.neighbor@cell.names, ncol(query.neighbor@nn.idx))
  } else {
    if (distance_metric %in% c(simil_method, "pearson", "spearman")) {
      if (distance_metric %in% c("pearson", "spearman")) {
        if (distance_metric == "spearman") {
          ref <- Matrix::t(apply(ref, 1, rank))
          query <- Matrix::t(apply(query, 1, rank))
        }
        distance_metric <- "correlation"
      }
      d <- 1 -
        proxyC::simil(
          x = SeuratObject::as.sparse(ref),
          y = SeuratObject::as.sparse(query),
          method = distance_metric,
          use_nan = TRUE
        )
    } else if (distance_metric %in% dist_method) {
      d <- proxyC::dist(
        x = SeuratObject::as.sparse(ref),
        y = SeuratObject::as.sparse(query),
        method = distance_metric,
        use_nan = TRUE
      )
    }
    if (k == 1) {
      match_k <- as_matrix(apply(
        d,
        2,
        function(x) order(x, decreasing = FALSE)[1]
      ))
      match_k_cell <- as_matrix(apply(
        d,
        2,
        function(x) names(x)[order(x, decreasing = FALSE)[1]]
      ))
      match_k_distance <- as_matrix(apply(
        d,
        2,
        function(x) x[order(x, decreasing = FALSE)[1]]
      ))
    } else {
      match_k <- Matrix::t(
        as_matrix(apply(
          d,
          2,
          function(x) order(x, decreasing = FALSE)[1:k]
        ))
      )
      match_k_cell <- Matrix::t(
        as_matrix(apply(
          d,
          2,
          function(x) names(x)[order(x, decreasing = FALSE)[1:k]]
        ))
      )
      match_k_distance <- Matrix::t(
        as_matrix(apply(
          d,
          2,
          function(x) x[order(x, decreasing = FALSE)[1:k]]
        ))
      )
    }
    knn_cells <- match_k_cell
    refumap_all <- srt_ref[[ref_umap]]@cell.embeddings[
      knn_cells, ,
      drop = FALSE
    ]
    group <- rep(colnames(d), k)
  }

  if (projection_method == "model") {
    if ("layout" %in% names(model)) {
      srt_query[["ref.embeddings"]] <- RunUMAP2(
        object = query,
        reduction.model = srt_ref[[ref_umap]],
        assay = query_assay
      )
      srt_query[["ref.embeddings"]]@misc[["reduction.model"]] <- ref_umap
    } else if ("embedding" %in% names(model)) {
      neighborlist <- list(idx = match_k, dist = match_k_distance)
      srt_query[["ref.embeddings"]] <- RunUMAP2(
        object = neighborlist,
        reduction.model = srt_ref[[ref_umap]],
        assay = query_assay
      )
      srt_query[["ref.embeddings"]]@misc[["reduction.model"]] <- ref_umap
    }
  } else {
    refumap <- stats::aggregate(refumap_all, by = list(group), FUN = vote_fun)
    rownames(refumap) <- refumap[, 1]
    refumap[, 1] <- NULL
    colnames(refumap) <- paste0("Dim_", seq_len(ncol(refumap)))
    refumap <- as_matrix(refumap)
    srt_query[["ref.embeddings"]] <- SeuratObject::CreateDimReducObject(
      embeddings = refumap,
      key = "Dim_",
      assay = query_assay,
      misc = list(reduction.model = ref_umap)
    )
  }

  if (!is.null(ref_group)) {
    log_message("Predicting cell types based on ref_group.") ## slow
    level <- as.character(unique(srt_ref[["ref_group", drop = TRUE]]))
    if (k == 1) {
      match_best <- srt_ref[["ref_group", drop = TRUE]][match_k_cell[, 1]]
      names(match_best) <- names(match_k_cell[, 1])
    } else {
      rn <- rownames(match_k_cell)
      match_k_cell <- matrix(
        srt_ref[["ref_group", drop = TRUE]][match_k_cell],
        nrow = nrow(match_k_cell),
        ncol = ncol(match_k_cell)
      )
      rownames(match_k_cell) <- rn
      match_freq <- apply(match_k_cell, 1, table)
      if (!inherits(match_freq, "list")) {
        match_freq <- as.list(stats::setNames(object = rep(k, nrow(match_k_cell)), rn))
        match_freq <- lapply(
          stats::setNames(names(match_freq), names(match_freq)),
          function(x) stats::setNames(k, match_k_cell[x, 1])
        )
      }
      match_prob <- lapply(
        match_freq, function(x) {
          x[level[!level %in% names(x)]] <- 0
          x <- x / sum(x)
          return(x)
        }
      ) %>%
        dplyr::bind_rows()
      match_prob <- as_matrix(match_prob)
      rownames(match_prob) <- names(match_freq)
      match_best <- apply(
        match_prob,
        1,
        function(x) names(x)[order(x, decreasing = TRUE)][1]
      )
    }
    srt_query[[paste0("predicted_", ref_group)]] <- match_best[colnames(
      srt_query
    )]
  }

  return(srt_query)
}
