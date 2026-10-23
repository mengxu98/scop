check_python_dependencies <- function(python, packages, envname) {
  output <- tempfile("pip-check-")
  on.exit(unlink(output), add = TRUE)
  status <- get_namespace_fun("reticulate", "system2t")(
    python, c("-m", "pip", "check"), stdout = output, stderr = output
  )
  if (status == 0L) {
    return(invisible(TRUE))
  }
  diagnostics <- trimws(readLines(output, warn = FALSE))
  diagnostics <- diagnostics[nzchar(diagnostics)]
  closure <- python_dependency_closure(python, packages)
  pattern <- "^([[:alnum:]][[:alnum:]_.-]*) [^[:space:]]+ (has requirement |requires |is not supported on this platform)"
  recognized <- grepl(pattern, diagnostics)
  owners <- canonical_python_distribution_name(sub(" .*", "", diagnostics))
  relevant <- !recognized | owners %in% closure
  detail <- paste(diagnostics, collapse = "\n")
  if (status != 1L || !length(diagnostics) || any(relevant)) {
    log_message(
      "Python dependency check failed in environment {.val {envname}}; conflicts affecting the requested dependency closure or unclassified failures remain:\n{detail}",
      message_type = "error"
    )
  }
  log_message(
    "Unrelated legacy Python dependency conflicts in environment {.val {envname}}:\n{detail}\nCreate a separate environment for the affected module, or rebuild it with PrepareEnv(force = TRUE) after preserving needed data.",
    message_type = "warning",
    verbose = TRUE
  )
  invisible(TRUE)
}

python_dependency_closure <- function(python, packages) {
  roots <- vapply(seq_along(packages), function(i) {
    pkg <- packages[[i]]
    if (grepl("^git\\+", pkg)) {
      name <- parse_python_requirement(pkg, fallback_name = names(packages)[i])$name
      pkg <- paste(name, "@", pkg)
    }
    pkg
  }, "")
  output <- tempfile("dependency-closure-")
  errors <- tempfile("dependency-errors-")
  on.exit(unlink(c(output, errors)), add = TRUE)
  script <- runner_script_path("dependency_closure.py", "dependency metadata")
  status <- get_namespace_fun("reticulate", "system2t")(
    python, shQuote(c(script, roots)), stdout = output, stderr = errors
  )
  if (status != 0L) {
    detail <- paste(readLines(errors, warn = FALSE), collapse = "\n")
    log_message(
      "Unable to determine the requested Python dependency closure:\n{detail}",
      message_type = "error"
    )
  }
  closure <- readLines(output, warn = FALSE)
  if (!length(closure) || any(!grepl("^[a-z0-9][a-z0-9-]*$", closure))) {
    log_message("Invalid Python dependency closure metadata", message_type = "error")
  }
  closure
}

scvi_jax_requirements <- function() {
  c(
    jax = "jax[cpu]==0.4.38",
    jaxlib = "jaxlib==0.4.38",
    numpyro = "numpyro==0.17.0",
    optax = "optax==0.2.5",
    flax = "flax==0.10.4",
    chex = "chex==0.1.89",
    "orbax-checkpoint" = "orbax-checkpoint==0.8.0",
    tensorstore = "tensorstore==0.1.65"
  )
}

assert_jax_module_compatibility <- function(modules) {
  if (all(c("scvi", "moscot") %in% modules)) {
    log_message(
      "The moscot JAX profile is incompatible with scvi-tools 1.3.3 pinned to JAX 0.4.38. Prepare these modules in separate environments.",
      message_type = "error"
    )
  }
  invisible(NULL)
}
