metric_graph_connectivity <- function(embeddings, labels, k = 15) {
  if (!requireNamespace("FNN", quietly = TRUE)) {
    log_message(
      "{.pkg FNN} is required to compute graph connectivity metrics.",
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
  knn <- FNN::get.knn(embeddings, k = k_use)
  edges <- cbind(
    rep(seq_len(nrow(embeddings)), each = k_use),
    as.vector(knn$nn.index)
  )
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
