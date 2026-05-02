metric_graph_connectivity <- function(
  embeddings,
  labels,
  k = 15,
  backend = c("fnn", "torch"),
  device = c("auto", "mps", "cpu"),
  torch_chunk_size = 2048
) {
  backend <- match.arg(backend)
  device <- match.arg(device)
  if (identical(backend, "fnn")) {
    metric_graph_connectivity_check_fnn()
  }

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

  edges <- switch(backend,
    fnn = metric_graph_connectivity_edges_fnn(
      embeddings = embeddings,
      k = k_use
    ),
    torch = metric_graph_connectivity_edges_torch(
      embeddings = embeddings,
      k = k_use,
      device = device,
      chunk_size = torch_chunk_size
    )
  )

  metric_graph_connectivity_from_edges(
    edges = edges,
    labels = labels
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

metric_graph_connectivity_edges_torch <- function(
  embeddings,
  k,
  device = c("auto", "mps", "cpu"),
  chunk_size = 2048
) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    log_message(
      "{.pkg torch} is required for {.arg backend = 'torch'} graph connectivity.",
      message_type = "error"
    )
  }
  if (!isTRUE(torch::torch_is_installed())) {
    log_message(
      "{.pkg torch} is installed, but LibTorch is not available. Run {.fn torch::install_torch} first.",
      message_type = "error"
    )
  }

  device <- match.arg(device)
  torch_device <- metric_graph_connectivity_torch_device(device)
  n_cells <- nrow(embeddings)
  chunk_size <- max(1L, min(as.integer(chunk_size), n_cells))

  x <- torch::torch_tensor(
    embeddings,
    dtype = torch::torch_float32(),
    device = torch_device
  )

  topk_n <- min(k + 1L, n_cells)
  starts <- seq.int(1L, n_cells, by = chunk_size)
  edge_from <- vector("list", length(starts))
  edge_to <- vector("list", length(starts))
  index_base <- NULL

  for (i in seq_along(starts)) {
    start_i <- starts[[i]]
    end_i <- min(start_i + chunk_size - 1L, n_cells)
    len_i <- end_i - start_i + 1L

    x_chunk <- torch::torch_narrow(
      x,
      dim = 1L,
      start = as.integer(start_i),
      length = as.integer(len_i)
    )
    dist <- torch::torch_cdist(
      x_chunk,
      x,
      p = 2
    )
    idx <- torch::torch_topk(
      dist,
      k = as.integer(topk_n),
      dim = 2L,
      largest = FALSE,
      sorted = TRUE
    )[[2]]
    idx <- as.matrix(torch::as_array(idx$to(device = "cpu")))
    if (is.null(index_base)) {
      index_base <- if (any(idx == 0L, na.rm = TRUE)) 0L else 1L
    }
    if (identical(index_base, 0L)) {
      idx <- idx + 1L
    }

    rows <- seq.int(start_i, end_i)
    nn <- mapply(
      function(row_i, idx_i) {
        idx_i <- as.integer(idx_i)
        idx_i <- idx_i[!is.na(idx_i) & idx_i >= 1L & idx_i <= n_cells]
        idx_i <- idx_i[idx_i != row_i]
        idx_i[seq_len(min(k, length(idx_i)))]
      },
      row_i = rows,
      idx_i = split(idx, row(idx)),
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE
    )

    edge_from[[i]] <- rep(rows, lengths(nn))
    edge_to[[i]] <- unlist(nn, use.names = FALSE)
  }

  cbind(
    unlist(edge_from, use.names = FALSE),
    unlist(edge_to, use.names = FALSE)
  )
}

metric_graph_connectivity_torch_device <- function(device = c("auto", "mps", "cpu")) {
  device <- match.arg(device)
  mps_available <- "backends_mps_is_available" %in% getNamespaceExports("torch") &&
    isTRUE(torch::backends_mps_is_available())

  if (identical(device, "mps") && !mps_available) {
    log_message(
      "{.pkg torch} MPS backend is not available on this machine.",
      message_type = "error"
    )
  }
  if (identical(device, "auto") && mps_available) {
    return(torch::torch_device("mps"))
  }
  if (identical(device, "mps")) {
    return(torch::torch_device("mps"))
  }

  torch::torch_device("cpu")
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
