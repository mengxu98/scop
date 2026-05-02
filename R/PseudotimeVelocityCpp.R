# C++-accelerated pseudotime velocity computation
# Internal helper functions

run_pseudotime_velocity_knn_cpp <- function(x_emb, pseudotime, neighbors, normalize = TRUE) {
  if (!run_pseudotime_velocity_knn_cpp_available()) {
    log_message(
      "{.arg backend = 'cpp'} requires the compiled {.pkg scop} shared library. Reinstall the package to build native code.",
      message_type = "error"
    )
  }
  if (!is.matrix(x_emb)) {
    x_emb <- as.matrix(x_emb)
  }
  storage.mode(x_emb) <- "double"
  if (!is.matrix(neighbors)) {
    neighbors <- as.matrix(neighbors)
  }
  storage.mode(neighbors) <- "integer"
  pseudotime_velocity_knn_cpp(
    x_emb = x_emb,
    pseudotime = as.numeric(pseudotime),
    neighbors = neighbors,
    normalize = isTRUE(normalize)
  )
}

run_pseudotime_velocity_knn_cpp_available <- function() {
  exists("pseudotime_velocity_knn_cpp", mode = "function") &&
    isTRUE(is.loaded("_scop_pseudotime_velocity_knn_cpp"))
}

run_pseudotime_velocity_gradient_cpp <- function(x_emb, pseudotime, neighbors, smooth = 0.5, normalize = TRUE) {
  if (!run_pseudotime_velocity_gradient_cpp_available()) {
    log_message(
      "{.arg backend = 'cpp'} requires the compiled {.pkg scop} shared library. Reinstall the package to build native code.",
      message_type = "error"
    )
  }
  if (!is.matrix(x_emb)) {
    x_emb <- as.matrix(x_emb)
  }
  storage.mode(x_emb) <- "double"
  if (!is.matrix(neighbors)) {
    neighbors <- as.matrix(neighbors)
  }
  storage.mode(neighbors) <- "integer"
  pseudotime_velocity_gradient_cpp(
    x_emb = x_emb,
    pseudotime = as.numeric(pseudotime),
    neighbors = neighbors,
    smooth = as.numeric(smooth),
    normalize = isTRUE(normalize)
  )
}

run_pseudotime_velocity_gradient_cpp_available <- function() {
  exists("pseudotime_velocity_gradient_cpp", mode = "function") &&
    isTRUE(is.loaded("_scop_pseudotime_velocity_gradient_cpp"))
}

.neighbors_list_to_matrix <- function(neighbors_list, k) {
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
