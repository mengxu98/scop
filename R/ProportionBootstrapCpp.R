# C++-accelerated bootstrap for proportion testing
# Internal helper functions

run_proportion_bootstrap_log2fd_cpp <- function(v1, v2, n_bootstrap = 1000, pseudocount = 1e-5) {
  if (!run_proportion_bootstrap_log2fd_cpp_available()) {
    log_message(
      "{.arg backend = 'cpp'} requires the compiled {.pkg scop} shared library. Reinstall the package to build native code.",
      message_type = "error"
    )
  }
  proportion_bootstrap_log2fd_cpp(
    v1 = as.numeric(v1),
    v2 = as.numeric(v2),
    n_bootstrap = as.integer(n_bootstrap),
    pseudocount = as.numeric(pseudocount)
  )
}

run_proportion_bootstrap_log2fd_cpp_available <- function() {
  exists("proportion_bootstrap_log2fd_cpp", mode = "function") &&
    isTRUE(is.loaded("_scop_proportion_bootstrap_log2fd_cpp"))
}

run_proportion_bootstrap_stats_cpp <- function(v1, v2, n_bootstrap = 1000, pseudocount = 1e-5) {
  if (!run_proportion_bootstrap_stats_cpp_available()) {
    log_message(
      "{.arg backend = 'cpp'} requires the compiled {.pkg scop} shared library. Reinstall the package to build native code.",
      message_type = "error"
    )
  }
  proportion_bootstrap_stats_cpp(
    v1 = as.numeric(v1),
    v2 = as.numeric(v2),
    n_bootstrap = as.integer(n_bootstrap),
    pseudocount = as.numeric(pseudocount)
  )
}

run_proportion_bootstrap_stats_cpp_available <- function() {
  exists("proportion_bootstrap_stats_cpp", mode = "function") &&
    isTRUE(is.loaded("_scop_proportion_bootstrap_stats_cpp"))
}
