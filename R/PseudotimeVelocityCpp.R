# C++-accelerated pseudotime velocity computation
# Internal helper functions

run_pseudotime_velocity_knn <- function(x_emb, pseudotime, neighbors, normalize = TRUE) {
  if (!is.matrix(x_emb)) {
    x_emb <- as.matrix(x_emb)
  }
  storage.mode(x_emb) <- "double"
  if (!is.matrix(neighbors)) {
    neighbors <- as.matrix(neighbors)
  }
  storage.mode(neighbors) <- "integer"
  pseudotime_velocity_knn(
    x_emb = x_emb,
    pseudotime = as.numeric(pseudotime),
    neighbors = neighbors,
    normalize = isTRUE(normalize)
  )
}

run_pseudotime_velocity_gradient <- function(x_emb, pseudotime, neighbors, smooth = 0.5, normalize = TRUE) {
  if (!is.matrix(x_emb)) {
    x_emb <- as.matrix(x_emb)
  }
  storage.mode(x_emb) <- "double"
  if (!is.matrix(neighbors)) {
    neighbors <- as.matrix(neighbors)
  }
  storage.mode(neighbors) <- "integer"
  pseudotime_velocity_gradient(
    x_emb = x_emb,
    pseudotime = as.numeric(pseudotime),
    neighbors = neighbors,
    smooth = as.numeric(smooth),
    normalize = isTRUE(normalize)
  )
}

neighbors_list_to_matrix <- function(neighbors_list, k) {
  # Convert a list of neighbor index vectors to an n x k integer matrix
  n <- length(neighbors_list)
  k_use <- min(k, max(vapply(neighbors_list, length, integer(1)), 0))
  mat <- matrix(NA_integer_, nrow = n, ncol = k_use)
  for (i in seq_len(n)) {
    nb <- neighbors_list[[i]]
    if (length(nb) > 0) {
      n_take <- min(length(nb), k_use)
      mat[i, seq_len(n_take)] <- as.integer(nb[seq_len(n_take)])
    }
  }
  mat
}
