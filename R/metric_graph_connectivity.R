metric_graph_connectivity <- function(
  embeddings,
  labels,
  k = 15,
  backend = c("auto", "cpp", "bioc_kmknn", "fnn")
) {
  backend <- match.arg(backend)

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

  backend <- metric_graph_connectivity_resolve_backend(
    backend = backend
  )
  edges <- switch(backend,
    cpp = metric_graph_connectivity_edges_cpp(
      embeddings = embeddings,
      k = k_use
    ),
    bioc_kmknn = metric_graph_connectivity_edges_bioc_kmknn(
      embeddings = embeddings,
      k = k_use
    ),
    fnn = metric_graph_connectivity_edges_fnn(
      embeddings = embeddings,
      k = k_use
    )
  )

  metric_graph_connectivity_from_edges(
    edges = edges,
    labels = labels
  )
}

metric_graph_connectivity_resolve_backend <- function(
  backend
) {
  if (!identical(backend, "auto")) {
    return(backend)
  }
  "cpp"
}

metric_graph_connectivity_edges_from_index <- function(
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

metric_graph_connectivity_edges_fnn <- function(embeddings, k) {
  metric_graph_connectivity_check_fnn()

  knn <- FNN::get.knn(embeddings, k = k)
  cbind(
    rep(seq_len(nrow(embeddings)), each = k),
    as.vector(t(knn$nn.index))
  )
}

metric_graph_connectivity_check_fnn <- function() {
  if (!requireNamespace("FNN", quietly = TRUE)) {
    log_message(
      "{.pkg FNN} is required to compute graph connectivity metrics.",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

metric_graph_connectivity_edges_bioc_kmknn <- function(embeddings, k) {
  check_r("BiocNeighbors", verbose = FALSE)

  knn <- BiocNeighbors::findKNN(
    embeddings,
    k = k,
    BNPARAM = BiocNeighbors::KmknnParam(distance = "Euclidean"),
    num.threads = 1L
  )
  metric_graph_connectivity_edges_from_index(
    index = knn$index,
    k = k,
    remove_self = FALSE
  )
}

metric_graph_connectivity_edges_cpp <- function(embeddings, k) {
  knn <- run_knn_topk(
    reference = embeddings,
    k = k,
    metric = "euclidean",
    backend = "cpp",
    exclude_self = TRUE
  )
  metric_graph_connectivity_edges_from_index(
    index = knn[["idx"]],
    k = k,
    remove_self = FALSE
  )
}

metric_graph_connectivity_from_edges <- function(edges, labels) {
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
