run_integration_methods <- function(
  args,
  integration_methods,
  runner = RunIntegration
) {
  result <- NULL

  for (i in seq_along(integration_methods)) {
    method_args <- args
    if (i > 1L) {
      method_args[["srt_merge"]] <- result
      method_args[["srt_list"]] <- NULL
    }
    method_args[["integration_method"]] <- NULL
    method_args[["integration_methods"]] <- integration_methods[[i]]
    result <- do.call(runner, method_args)
  }

  result
}

collect_integration_metrics <- function(
  srt,
  reduction,
  batch_col,
  celltype_col = NULL,
  cluster_col = NULL,
  lisi_tool_name = NULL,
  lisi_prefix = NULL,
  k_graph = 15
) {
  emb <- Seurat::Embeddings(srt, reduction = reduction)
  summary_list <- list()
  if (!is.null(batch_col) && batch_col %in% colnames(srt@meta.data)) {
    summary_list[["batch_ASW_mixing"]] <- metric_silhouette(
      embeddings = emb,
      labels = srt[[batch_col, drop = TRUE]],
      maximize = FALSE
    )
  }
  if (!is.null(celltype_col) && celltype_col %in% colnames(srt@meta.data)) {
    summary_list[["celltype_ASW"]] <- metric_silhouette(
      embeddings = emb,
      labels = srt[[celltype_col, drop = TRUE]],
      maximize = TRUE
    )
    summary_list[["celltype_graph_connectivity"]] <- tryCatch(
      metric_graph_connectivity(
        embeddings = emb,
        labels = srt[[celltype_col, drop = TRUE]],
        k = k_graph
      ),
      error = function(e) NA_real_
    )
  }
  if (
    !is.null(cluster_col) &&
      cluster_col %in% colnames(srt@meta.data) &&
      !is.null(celltype_col) &&
      celltype_col %in% colnames(srt@meta.data)
  ) {
    metrics <- classification_metrics_compute(
      predicted = srt[[cluster_col, drop = TRUE]],
      truth = srt[[celltype_col, drop = TRUE]]
    )
    summary_list[["celltype_NMI"]] <- metrics[["nmi"]]
    summary_list[["celltype_ARI"]] <- metrics[["ari"]]
    summary_list[["celltype_purity"]] <- metrics[["purity"]]
  }
  if (!is.null(lisi_tool_name) && lisi_tool_name %in% names(srt@tools)) {
    lisi_res <- srt@tools[[lisi_tool_name]]
    if (!is.null(lisi_res$label_colnames) && !is.null(lisi_res$scores)) {
      for (label in lisi_res$label_colnames) {
        cols <- grep(
          paste0("_", label, "_LISI$"),
          colnames(lisi_res$scores),
          value = TRUE
        )
        if (!is.null(lisi_prefix)) {
          cols <- cols[grepl(
            paste0("(^|\\.)", make.names(lisi_prefix), "_"),
            cols
          )]
        }
        if (length(cols) > 0) {
          summary_list[[paste0(label, "_LISI_mean")]] <- mean(
            as.numeric(as.matrix(lisi_res$scores[, cols, drop = FALSE])),
            na.rm = TRUE
          )
        }
      }
    }
  }
  if (length(summary_list) == 0L) {
    raw_df <- data.frame(
      metric = character(),
      value = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    raw_df <- data.frame(
      metric = names(summary_list),
      value = unlist(summary_list, use.names = FALSE),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }
  format_integration_metrics(
    raw_df = raw_df,
    srt = srt,
    batch_col = batch_col,
    celltype_col = celltype_col
  )
}

resolve_cluster_algorithm_index <- function(cluster_algorithm) {
  cluster_algorithms <- c("louvain", "slm", "leiden")
  if (!cluster_algorithm %in% cluster_algorithms) {
    log_message(
      "{.arg cluster_algorithm} must be one of {.val {cluster_algorithms}}",
      message_type = "error"
    )
  }

  switch(
    EXPR = tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )
}

validate_nonlinear_reductions <- function(
  nonlinear_reduction,
  allowed = c(
    "umap",
    "umap-naive",
    "tsne",
    "dm",
    "phate",
    "pacmap",
    "trimap",
    "largevis",
    "fr"
  )
) {
  if (any(!nonlinear_reduction %in% allowed)) {
    log_message(
      "{.arg nonlinear_reduction} must be one of {.val {allowed}}",
      message_type = "error"
    )
  }

  invisible(nonlinear_reduction)
}

validate_integration_input_cells <- function(srt_list, srt_merge) {
  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "{.arg srt_list} and {.arg srt_merge} were all empty",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    list_cells <- sort(unique(unlist(lapply(srt_list, colnames))))
    merge_cells <- sort(unique(colnames(srt_merge)))
    if (!identical(list_cells, merge_cells)) {
      log_message(
        "{.arg srt_list} and {.arg srt_merge} have different cells",
        message_type = "error"
      )
    }
  }
  invisible(TRUE)
}

find_neighbors_and_clusters <- function(
  srt,
  reduction,
  dims_use,
  graph_prefix,
  graph_snn,
  cluster_colname,
  HVF,
  neighbor_metric,
  neighbor_k,
  cluster_algorithm,
  cluster_algorithm_index,
  cluster_resolution,
  run_find_neighbors = TRUE,
  verbose
) {
  neighbor_k_use <- min(as.integer(neighbor_k), max(1L, ncol(srt) - 1L))
  if (!identical(neighbor_k_use, as.integer(neighbor_k))) {
    log_message(
      "Adjust neighbor k from {.val {neighbor_k}} to {.val {neighbor_k_use}} for small-sample clustering",
      verbose = verbose
    )
  }
  srt <- tryCatch(
    {
      if (isTRUE(run_find_neighbors)) {
        srt <- FindNeighbors(
          object = srt,
          reduction = reduction,
          dims = dims_use,
          annoy.metric = neighbor_metric,
          k.param = neighbor_k_use,
          graph.name = paste0(graph_prefix, c("KNN", "SNN")),
          verbose = FALSE
        )
      }

      log_message(
        "Perform {.fn Seurat::FindClusters} with {.val {cluster_algorithm}}",
        verbose = verbose
      )
      srt <- FindClusters(
        object = srt,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        leiden_method = "igraph",
        graph.name = graph_snn,
        verbose = FALSE
      )
      log_message("Reorder clusters...")
      srt <- srt_reorder(
        srt,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srt[["seurat_clusters"]] <- NULL
      srt[[cluster_colname]] <- SeuratObject::Idents(srt)
      srt
    },
    error = function(error) {
      err_msg <- conditionMessage(error)
      err_msg <- gsub("{", "{{", err_msg, fixed = TRUE)
      err_msg <- gsub("}", "}}", err_msg, fixed = TRUE)
      log_message(err_msg, message_type = "warning", verbose = verbose)
      log_message(
        "Error when performing {.fn Seurat::FindClusters}. Skip this step",
        message_type = "warning",
        verbose = verbose
      )
      srt
    }
  )

  return(srt)
}

resolve_linear_dims_use <- function(
  srt,
  reduction,
  linear_reduction_dims_use = NULL,
  normalization_method = "LogNormalize",
  reduction_method = NULL,
  verbose = FALSE
) {
  if (!is.null(linear_reduction_dims_use)) {
    return(linear_reduction_dims_use)
  }
  RunDimsEstimate(
    object = srt,
    reduction = reduction,
    reduction_method = reduction_method,
    skip_first = normalization_method == "TFIDF",
    use_stored = TRUE,
    verbose = verbose
  )
}

run_nonlinear_reduction <- function(
  srt,
  prefix,
  reduction_use = NULL,
  reduction_dims = NULL,
  graph_use = NULL,
  neighbor_use = NULL,
  nonlinear_reduction,
  nonlinear_reduction_dims,
  nonlinear_reduction_params,
  force_nonlinear_reduction,
  seed,
  verbose
) {
  if (
    !is.null(reduction_use) && reduction_use %in% SeuratObject::Reductions(srt)
  ) {
    available_dims <- seq_len(
      ncol(Seurat::Embeddings(srt, reduction = reduction_use))
    )
    reduction_dims <- intersect(reduction_dims, available_dims)
    if (length(reduction_dims) == 0) {
      log_message(
        "No valid dimensions remain for {.arg reduction_use = '{reduction_use}'}",
        message_type = "warning",
        verbose = verbose
      )
      return(srt)
    }
  }
  srt <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        params_use <- nonlinear_reduction_params
        if (nr %in% c("fr")) {
          params_use[["n.neighbors"]] <- NULL
        }
        for (n in nonlinear_reduction_dims) {
          srt <- RunDimsReduction(
            srt,
            prefix = prefix,
            reduction_use = reduction_use,
            reduction_dims = reduction_dims,
            graph_use = graph_use,
            neighbor_use = neighbor_use,
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = params_use,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = verbose,
            seed = seed
          )
        }
      }
      srt
    },
    error = function(error) {
      err_msg <- conditionMessage(error)
      err_msg <- gsub("{", "{{", err_msg, fixed = TRUE)
      err_msg <- gsub("}", "}}", err_msg, fixed = TRUE)
      log_message(err_msg, message_type = "warning", verbose = verbose)
      log_message(
        "Error when performing nonlinear dimension reduction. Skip this step",
        message_type = "warning",
        verbose = verbose
      )
      srt
    }
  )

  return(srt)
}
metric_graph_connectivity <- function(
  embeddings,
  labels,
  k = 15,
  backend = "r"
) {
  backend <- match.arg(backend, "r")

  embeddings <- as.matrix(embeddings)
  storage.mode(embeddings) <- "double"
  if (nrow(embeddings) != length(labels)) {
    log_message(
      "{.arg embeddings} rows must match {.arg labels} length",
      message_type = "error"
    )
  }

  labels <- as.factor(labels)
  keep <- !is.na(labels)
  embeddings <- embeddings[keep, , drop = FALSE]
  labels <- droplevels(labels[keep])
  if (nrow(embeddings) < 3 || nlevels(labels) < 1) {
    return(NA_real_)
  }

  k_use <- min(as.integer(k), nrow(embeddings) - 1L)
  if (k_use < 1) {
    return(NA_real_)
  }

  edges <- graph_conn_edges_r(
    embeddings = embeddings,
    k = k_use
  )

  graph_conn_score(
    edges = edges,
    labels = labels
  )
}

graph_conn_edges_from_index <- function(
  index,
  k,
  remove_self = TRUE
) {
  rows <- seq_len(nrow(index))
  nn <- lapply(rows, function(i) {
    idx <- as.integer(index[i, ])
    idx <- idx[!is.na(idx) & idx >= 1L]
    if (isTRUE(remove_self)) {
      idx <- idx[idx != i]
    }
    idx[seq_len(min(k, length(idx)))]
  })

  cbind(
    rep(rows, lengths(nn)),
    unlist(nn, use.names = FALSE)
  )
}

graph_conn_edges_r <- function(embeddings, k) {
  if (isTRUE(all(unlist(check_r("BiocNeighbors", install = FALSE, verbose = FALSE), use.names = FALSE)))) {
    find_knn <- get_namespace_fun("BiocNeighbors", "findKNN")
    kmknn_param <- get_namespace_fun("BiocNeighbors", "KmknnParam")
    knn <- find_knn(
      embeddings,
      k = k,
      BNPARAM = kmknn_param(distance = "Euclidean"),
      num.threads = 1L
    )
    return(graph_conn_edges_from_index(
      index = knn$index,
      k = k,
      remove_self = FALSE
    ))
  }
  d <- as.matrix(stats::dist(embeddings))
  diag(d) <- Inf
  index <- t(apply(d, 1L, function(row) {
    as.integer(order(row)[seq_len(k)])
  }))
  graph_conn_edges_from_index(
    index = index,
    k = k,
    remove_self = FALSE
  )
}


graph_conn_score <- function(edges, labels) {
  if (is.null(edges) || nrow(edges) == 0L) {
    return(NA_real_)
  }

  graph <- igraph::graph_from_edgelist(edges, directed = FALSE)
  graph <- igraph::simplify(graph)
  per_label <- tapply(seq_along(labels), labels, function(idx) {
    if (length(idx) <= 1) {
      return(1)
    }
    comps <- igraph::components(igraph::induced_subgraph(graph, vids = idx))
    max(comps$csize) / length(idx)
  })
  mean(unlist(per_label), na.rm = TRUE)
}
metric_silhouette <- function(
  embeddings,
  labels,
  maximize = TRUE,
  n_max = 800L,
  seed = 11
) {
  check_r("cluster", verbose = FALSE)
  labels <- as.factor(labels)
  keep <- !is.na(labels)
  embeddings <- embeddings[keep, , drop = FALSE]
  labels <- droplevels(labels[keep])
  if (nrow(embeddings) < 3 || nlevels(labels) < 2) {
    return(NA_real_)
  }
  n_max <- as.integer(n_max)
  if (is.finite(n_max) && n_max >= 10L && nrow(embeddings) > n_max) {
    set.seed(seed)
    grouped <- split(seq_len(nrow(embeddings)), labels)
    idx <- unlist(lapply(grouped, function(cells) {
      take <- max(2L, as.integer(round(n_max * length(cells) / nrow(embeddings))))
      take <- min(length(cells), take)
      sample(cells, take)
    }), use.names = FALSE)
    if (length(idx) > n_max) {
      idx <- sample(idx, n_max)
    }
    embeddings <- embeddings[idx, , drop = FALSE]
    labels <- droplevels(labels[idx])
    if (nrow(embeddings) < 3 || nlevels(labels) < 2) {
      return(NA_real_)
    }
  }
  sil <- cluster::silhouette(
    x = as.integer(labels),
    dist = stats::dist(embeddings)
  )
  score <- mean(sil[, "sil_width"], na.rm = TRUE)
  if (isTRUE(maximize)) {
    return(score)
  }
  1 - abs(score)
}
