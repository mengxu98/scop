#' @title Benchmark spatial domain clustering methods
#'
#' @description
#' Run supported spatial domain clustering methods from the same immutable
#' input, compare their aligned labels with a gold standard, and record
#' classification quality, elapsed time, and sampled peak process-tree memory.
#' Each method runs in an isolated R process so one failed backend does not
#' corrupt the input or prevent the remaining methods from being assessed.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt A spatial `Seurat` object.
#' @param gold_standard Either one metadata column in `srt` or a named vector
#' whose names match the spot names in `srt`.
#' @param methods Spatial domain methods to benchmark. `NULL` uses every
#' benchmarked producer (`BayesSpace`, `BANKSY`, and `SmoothClust`). Method
#' names may be written with or without the `Run` prefix.
#' @param method_params Named list of per-method argument lists. Arguments are
#' passed to the corresponding SCOP producer, never directly to its backend.
#' @param n_clusters Optional common number of domains. When `NULL`, methods
#' that require a cluster count use the number of non-missing gold-standard
#' classes. A method-specific `q` or `n_clusters` in `method_params` wins.
#' @param metrics Quality metrics selected by default when plotting. Supported
#' values are `"ARI"`, `"NMI"`, and `"purity"`. All three are retained in
#' the result summary.
#' @param keep_objects Whether to keep each successful method's full producer
#' result. The default keeps only aligned labels and compact run metadata.
#' @param install_missing Whether missing optional backends may enter their
#' producer's normal `check_r()` installation path. The default records them as
#' unavailable without changing the R library.
#' @param seed Seed set inside each isolated method process and forwarded to
#' producers with a public `seed` argument unless overridden in
#' `method_params`.
#' @param timeout Maximum wall time in seconds for each isolated run, including
#' process start and result serialization. `Inf` disables the timeout.
#' @param poll_interval Seconds between process-tree memory samples.
#'
#' @return A `benchmark_result` object. Use `result$summary` for the quality
#' table, [BenchmarkPlot()] for visualization, and `$predictions` for the
#' aligned labels.
#' @export
#'
#' @examples
#' \dontrun{
#' data(visium_human_pancreas_sub)
#' visium_human_pancreas_sub$gold_domain <- factor(
#'   paste0("domain_", (seq_len(ncol(visium_human_pancreas_sub)) - 1) %% 3 + 1)
#' )
#' bench <- RunBenchmark(
#'   visium_human_pancreas_sub,
#'   gold_standard = "gold_domain",
#'   method_params = list(
#'     BayesSpace = list(n.PCs = 5, n.HVGs = 200),
#'     BANKSY = list(layer = "counts"),
#'     SmoothClust = list(layer = "counts", min_spots = 1)
#'   )
#' )
#' bench
#' BenchmarkPlot(data = bench)
#' }
RunBenchmark <- function(
  srt,
  gold_standard,
  methods = NULL,
  method_params = list(),
  n_clusters = NULL,
  metrics = c("ARI", "NMI"),
  keep_objects = FALSE,
  install_missing = FALSE,
  seed = 1,
  timeout = Inf,
  poll_interval = 0.1,
  verbose = TRUE
) {
  benchmark_validate_input(srt)
  benchmark_assert_flag(keep_objects, "keep_objects")
  benchmark_assert_flag(install_missing, "install_missing")
  benchmark_assert_number(
    seed, "seed",
    lower = 0, upper = .Machine$integer.max, integer = TRUE
  )
  benchmark_assert_number(timeout, "timeout", lower = 0, allow_infinite = TRUE)
  benchmark_assert_number(poll_interval, "poll_interval", lower = 0, strict = TRUE)

  truth <- benchmark_resolve_truth(srt, gold_standard)
  n_truth_clusters <- length(unique(truth$labels[!is.na(truth$labels)]))
  if (is.null(n_clusters)) {
    n_clusters <- n_truth_clusters
  } else {
    benchmark_assert_number(
      n_clusters, "n_clusters",
      lower = 2,
      upper = .Machine$integer.max, integer = TRUE
    )
    n_clusters <- as.integer(n_clusters)
  }
  metrics <- benchmark_resolve_metrics(metrics)
  methods <- benchmark_resolve_methods(methods)
  method_params <- benchmark_normalize_method_params(method_params, methods)
  prepared_params <- lapply(methods, function(method) {
    benchmark_complete_method_params(
      method = method,
      params = method_params[[method]] %||% list(),
      n_clusters = n_clusters,
      seed = seed
    )
  })
  names(prepared_params) <- methods
  availability_results <- lapply(methods, benchmark_method_availability_safe)
  names(availability_results) <- methods
  will_execute <- isTRUE(install_missing) | vapply(
    availability_results,
    function(availability) identical(availability$status, "available"),
    logical(1)
  )

  work_dir <- NULL
  input_path <- NULL
  package_context <- NULL
  if (any(will_execute)) {
    benchmark_require_runtime()
    work_dir <- tempfile("scop-benchmark-")
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)
    input_path <- file.path(work_dir, "input.rds")
    saveRDS(srt, input_path, compress = FALSE)
    package_context <- benchmark_package_context()
  }

  run_results <- vector("list", length(methods))
  names(run_results) <- methods

  for (method in methods) {
    params <- prepared_params[[method]]
    availability <- availability_results[[method]]
    if (!will_execute[[method]]) {
      run_results[[method]] <- benchmark_unavailable_result(
        method = method,
        params = params,
        availability = availability,
        poll_interval = poll_interval
      )
      log_message(
        "Skip {.pkg {method}} because its backend is {.val {availability$status}}",
        message_type = "warning",
        verbose = verbose
      )
      next
    }

    log_message(
      "Benchmark {.pkg {method}} in an isolated R process",
      message_type = "running",
      verbose = verbose
    )
    run_dir <- file.path(work_dir, make.names(method))
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    run_results[[method]] <- benchmark_run_isolated(
      input_path = input_path,
      run_dir = run_dir,
      method = method,
      params = params,
      seed = seed,
      keep_object = keep_objects,
      timeout = timeout,
      poll_interval = poll_interval,
      package_context = package_context
    )
    if (identical(run_results[[method]]$status, "success")) {
      log_message(
        "Benchmark {.pkg {method}} completed",
        message_type = "success",
        verbose = verbose
      )
    } else {
      log_message(
        "Benchmark {.pkg {method}} ended with {.val {run_results[[method]]$status}}: {.val {run_results[[method]]$error}}",
        message_type = "warning",
        verbose = verbose
      )
    }
  }

  benchmark_build_result(
    srt = srt,
    truth = truth,
    methods = methods,
    run_results = run_results,
    metrics = metrics,
    n_clusters = n_clusters,
    keep_objects = keep_objects,
    install_missing = install_missing,
    seed = seed,
    timeout = timeout,
    poll_interval = poll_interval
  )
}

benchmark_method_aliases <- c(
  "bayesspace" = "BayesSpace",
  "runbayesspace" = "BayesSpace",
  "banksy" = "BANKSY",
  "runbanksy" = "BANKSY",
  "smoothclust" = "SmoothClust",
  "runsmoothclust" = "SmoothClust"
)

benchmark_method_runner <- function(method) {
  switch(method,
    BayesSpace = RunBayesSpace,
    BANKSY = RunBANKSY,
    SmoothClust = RunSmoothClust,
    log_message("Unknown benchmark method {.val {method}}", message_type = "error")
  )
}

benchmark_method_package <- function(method) {
  switch(method,
    BayesSpace = "BayesSpace",
    BANKSY = "Banksy",
    SmoothClust = "smoothclust",
    log_message("Unknown benchmark method {.val {method}}", message_type = "error")
  )
}

benchmark_resolve_methods <- function(methods = NULL) {
  defaults <- c("BayesSpace", "BANKSY", "SmoothClust")
  if (is.null(methods)) {
    return(defaults)
  }
  if (!is.character(methods) || length(methods) == 0L || anyNA(methods)) {
    log_message("{.arg methods} must contain one or more method names", message_type = "error")
  }
  resolved <- unname(benchmark_method_aliases[tolower(methods)])
  if (anyNA(resolved)) {
    unknown <- methods[is.na(resolved)]
    log_message(
      "Unknown benchmark methods: {.val {unknown}}. Available methods are {.val {defaults}}",
      message_type = "error"
    )
  }
  unique(resolved)
}

benchmark_normalize_method_params <- function(method_params, methods) {
  if (!is.list(method_params)) {
    log_message("{.arg method_params} must be a named list", message_type = "error")
  }
  out <- stats::setNames(rep(list(list()), length(methods)), methods)
  if (length(method_params) == 0L) {
    return(out)
  }
  nms <- names(method_params)
  if (is.null(nms) || any(!nzchar(nms))) {
    log_message("{.arg method_params} must be named by method", message_type = "error")
  }
  resolved <- unname(benchmark_method_aliases[tolower(nms)])
  if (anyNA(resolved)) {
    log_message(
      "Unknown methods in {.arg method_params}: {.val {nms[is.na(resolved)]}}",
      message_type = "error"
    )
  }
  if (anyDuplicated(resolved)) {
    log_message("{.arg method_params} contains duplicated method aliases", message_type = "error")
  }
  outside <- setdiff(resolved, methods)
  if (length(outside) > 0L) {
    log_message(
      "{.arg method_params} includes methods not selected by {.arg methods}: {.val {outside}}",
      message_type = "error"
    )
  }
  for (i in seq_along(method_params)) {
    value <- method_params[[i]]
    if (
      !is.list(value) ||
        (
          length(value) > 0L &&
            (
              is.null(names(value)) || any(!nzchar(names(value))) ||
                anyDuplicated(names(value))
            )
        )
    ) {
      log_message(
        "Parameters for {.val {nms[[i]]}} must be a uniquely named list",
        message_type = "error"
      )
    }
    if ("srt" %in% names(value)) {
      log_message("{.arg method_params} cannot replace the benchmark input {.arg srt}", message_type = "error")
    }
    out[[resolved[[i]]]] <- value
  }
  out
}

benchmark_complete_method_params <- function(method, params, n_clusters, seed) {
  cluster_arg <- switch(method,
    BayesSpace = "q",
    BANKSY = NA_character_,
    SmoothClust = "n_clusters",
    log_message("Unknown benchmark method {.val {method}}", message_type = "error")
  )
  if (!is.na(cluster_arg) && !cluster_arg %in% names(params)) {
    params[[cluster_arg]] <- as.integer(n_clusters)
  }
  seed_arg <- if (identical(method, "BayesSpace")) NA_character_ else "seed"
  if (!is.na(seed_arg) && !seed_arg %in% names(params)) {
    params[[seed_arg]] <- as.integer(seed)
  }
  if (!"verbose" %in% names(params)) params$verbose <- FALSE
  params
}

benchmark_validate_input <- function(srt) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  cells <- colnames(srt)
  if (length(cells) < 3L || anyNA(cells) || any(!nzchar(cells)) || anyDuplicated(cells)) {
    log_message(
      "{.arg srt} must contain at least three uniquely named spots",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

benchmark_resolve_truth <- function(srt, gold_standard) {
  cells <- colnames(srt)
  source <- NULL
  if (
    is.character(gold_standard) && length(gold_standard) == 1L &&
      !is.na(gold_standard) && gold_standard %in% colnames(srt@meta.data)
  ) {
    labels <- srt@meta.data[cells, gold_standard, drop = TRUE]
    source <- list(type = "metadata", column = gold_standard)
  } else {
    if (is.list(gold_standard) || length(gold_standard) == 0L) {
      log_message(
        "{.arg gold_standard} must be a metadata column or named label vector",
        message_type = "error"
      )
    }
    label_names <- names(gold_standard)
    if (is.null(label_names) || anyNA(label_names) || any(!nzchar(label_names)) || anyDuplicated(label_names)) {
      log_message(
        "A vector {.arg gold_standard} must have unique non-empty spot names",
        message_type = "error"
      )
    }
    missing <- setdiff(cells, label_names)
    extra <- setdiff(label_names, cells)
    if (length(missing) > 0L || length(extra) > 0L) {
      log_message(
        "{.arg gold_standard} names must exactly match {.arg srt} spots; missing {.val {missing}}, extra {.val {extra}}",
        message_type = "error"
      )
    }
    labels <- gold_standard[cells]
    source <- list(type = "vector", column = NA_character_)
  }
  labels <- as.character(labels)
  labels[is.na(labels) | !nzchar(labels)] <- NA_character_
  comparable <- labels[!is.na(labels)]
  if (length(comparable) < 2L || length(unique(comparable)) < 2L) {
    log_message(
      "{.arg gold_standard} must contain at least two non-missing classes",
      message_type = "error"
    )
  }
  list(labels = stats::setNames(labels, cells), source = source)
}

benchmark_resolve_metrics <- function(metrics) {
  allowed <- c(ari = "ARI", nmi = "NMI", purity = "purity")
  if (!is.character(metrics) || length(metrics) == 0L || anyNA(metrics)) {
    log_message("{.arg metrics} must select ARI, NMI, or purity", message_type = "error")
  }
  keys <- tolower(metrics)
  if (any(!keys %in% names(allowed))) {
    log_message(
      "Unsupported {.arg metrics}: {.val {metrics[!keys %in% names(allowed)]}}",
      message_type = "error"
    )
  }
  unique(unname(allowed[keys]))
}

benchmark_assert_flag <- function(x, arg) {
  validate_scalar_flag(x, arg)
}

benchmark_assert_number <- function(
  x,
  arg,
  lower = -Inf,
  upper = Inf,
  integer = FALSE,
  strict = FALSE,
  allow_infinite = FALSE
) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x)
  if (!allow_infinite) valid <- valid && is.finite(x)
  if (strict) valid <- valid && x > lower else valid <- valid && x >= lower
  valid <- valid && x <= upper
  if (integer) valid <- valid && isTRUE(all.equal(x, round(x)))
  if (!valid) {
    log_message("{.arg {arg}} has an invalid value", message_type = "error")
  }
  invisible(TRUE)
}

benchmark_require_runtime <- function() {
  missing <- c("callr", "ps")[!vapply(
    c("callr", "ps"),
    requireNamespace,
    logical(1),
    quietly = TRUE
  )]
  if (length(missing) > 0L) {
    log_message(
      "{.pkg {missing}} {?is/are} required for isolated benchmark execution and memory measurement",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

benchmark_method_availability <- function(method) {
  package <- benchmark_method_package(method)
  installed <- isTRUE(requireNamespace(package, quietly = TRUE))
  list(
    status = if (installed) "available" else "missing",
    detail = paste0(package, "=", if (installed) "available" else "missing"),
    versions = benchmark_backend_versions(package)
  )
}

benchmark_method_availability_safe <- function(method) {
  tryCatch(
    {
      availability <- benchmark_method_availability(method)
      if (
        !is.list(availability) ||
          !is.character(availability$status) || length(availability$status) != 1L ||
          is.na(availability$status) || !nzchar(availability$status) ||
          !is.character(availability$detail) || length(availability$detail) != 1L ||
          is.na(availability$detail)
      ) {
        log_message("backend availability diagnostic returned an invalid result", message_type = "error")
      }
      if (is.null(availability$versions)) availability$versions <- character()
      availability
    },
    error = function(error) {
      list(
        status = "diagnostic_error",
        detail = paste("Backend availability check failed:", conditionMessage(error)),
        versions = character()
      )
    }
  )
}

benchmark_backend_versions <- function(package) {
  version <- tryCatch(
    as.character(utils::packageVersion(package)),
    error = function(e) NA_character_
  )
  if (is.na(version)) {
    return(character())
  }
  stats::setNames(version, package)
}

benchmark_param_string <- function(params, name, default) {
  value <- params[[name]]
  if (
    is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
  ) {
    value
  } else {
    default
  }
}

benchmark_unavailable_result <- function(method, params, availability, poll_interval) {
  now <- Sys.time()
  list(
    status = "unavailable",
    error = availability$detail,
    prediction = NULL,
    object = NULL,
    runtime_s = NA_real_,
    baseline_memory_mb = NA_real_,
    peak_memory_mb = NA_real_,
    memory_delta_mb = NA_real_,
    started_at = now,
    finished_at = now,
    parameters = params,
    cluster_colname = benchmark_param_string(
      params, "cluster_colname", paste0(method, "_cluster")
    ),
    backend_versions = availability$versions,
    memory_method = "sampled process-tree RSS",
    poll_interval = poll_interval
  )
}

benchmark_package_context <- function() {
  namespace <- asNamespace("scop")
  path <- tryCatch(getNamespaceInfo(namespace, "path"), error = function(e) "")
  is_source <- nzchar(path) &&
    file.exists(file.path(path, "DESCRIPTION")) &&
    file.exists(file.path(path, ".Rbuildignore")) &&
    dir.exists(file.path(path, "R"))
  list(
    source_path = if (is_source) normalizePath(path, winslash = "/", mustWork = TRUE) else NULL,
    libpath = .libPaths()
  )
}

benchmark_run_isolated <- function(
  input_path,
  run_dir,
  method,
  params,
  seed,
  keep_object,
  timeout,
  poll_interval,
  package_context,
  child_entry = benchmark_child_entry
) {
  ready_path <- file.path(run_dir, "ready")
  go_path <- file.path(run_dir, "go")
  done_path <- file.path(run_dir, "done")
  result_path <- file.path(run_dir, "result.rds")
  stdout_path <- file.path(run_dir, "stdout.log")
  stderr_path <- file.path(run_dir, "stderr.log")
  started_at <- Sys.time()

  process <- tryCatch(
    callr::r_bg(
      func = child_entry,
      args = list(
        input_path = input_path,
        result_path = result_path,
        ready_path = ready_path,
        go_path = go_path,
        done_path = done_path,
        method = method,
        params = params,
        seed = seed,
        keep_object = keep_object,
        source_path = package_context$source_path
      ),
      libpath = package_context$libpath,
      stdout = stdout_path,
      stderr = stderr_path,
      supervise = TRUE,
      user_profile = FALSE,
      system_profile = FALSE
    ),
    error = function(error) error
  )
  if (inherits(process, "error")) {
    return(benchmark_failed_isolated_result(
      status = "failed", error = conditionMessage(process),
      method = method, params = params, started_at = started_at,
      poll_interval = poll_interval
    ))
  }

  timed_out <- FALSE
  start_clock <- proc.time()[["elapsed"]]
  deadline_reached <- function() {
    is.finite(timeout) && (proc.time()[["elapsed"]] - start_clock) > timeout
  }
  while (!file.exists(ready_path) && process$is_alive()) {
    if (deadline_reached()) {
      timed_out <- TRUE
      break
    }
    Sys.sleep(min(0.05, poll_interval))
  }
  if (timed_out) {
    benchmark_kill_process_tree(process)
    return(benchmark_failed_isolated_result(
      status = "timeout", error = "Timed out before backend execution started",
      method = method, params = params, started_at = started_at,
      poll_interval = poll_interval
    ))
  }
  if (!file.exists(ready_path)) {
    process$wait(timeout = 5000)
    return(benchmark_failed_isolated_result(
      status = "failed", error = benchmark_process_error(stderr_path, stdout_path),
      method = method, params = params, started_at = started_at,
      poll_interval = poll_interval
    ))
  }

  baseline_bytes <- benchmark_process_tree_rss(process$get_pid())
  peak_bytes <- baseline_bytes
  file.create(go_path)
  while (!file.exists(done_path) && process$is_alive()) {
    current <- benchmark_process_tree_rss(process$get_pid())
    if (is.finite(current)) peak_bytes <- max(peak_bytes, current, na.rm = TRUE)
    if (deadline_reached()) {
      timed_out <- TRUE
      break
    }
    Sys.sleep(poll_interval)
  }
  current <- benchmark_process_tree_rss(process$get_pid())
  if (is.finite(current)) peak_bytes <- max(peak_bytes, current, na.rm = TRUE)
  if (timed_out) {
    benchmark_kill_process_tree(process)
    return(benchmark_failed_isolated_result(
      status = "timeout", error = "Benchmark method exceeded timeout",
      method = method, params = params, started_at = started_at,
      baseline_bytes = baseline_bytes, peak_bytes = peak_bytes,
      poll_interval = poll_interval
    ))
  }

  while (process$is_alive()) {
    if (deadline_reached()) {
      timed_out <- TRUE
      break
    }
    Sys.sleep(min(0.05, poll_interval))
  }
  if (timed_out) {
    benchmark_kill_process_tree(process)
    return(benchmark_failed_isolated_result(
      status = "timeout", error = "Benchmark result serialization exceeded timeout",
      method = method, params = params, started_at = started_at,
      baseline_bytes = baseline_bytes, peak_bytes = peak_bytes,
      poll_interval = poll_interval
    ))
  }
  process$wait(timeout = 5000)
  if (!file.exists(result_path)) {
    return(benchmark_failed_isolated_result(
      status = "failed", error = benchmark_process_error(stderr_path, stdout_path),
      method = method, params = params, started_at = started_at,
      baseline_bytes = baseline_bytes, peak_bytes = peak_bytes,
      poll_interval = poll_interval
    ))
  }
  child <- tryCatch(readRDS(result_path), error = function(error) error)
  if (inherits(child, "error")) {
    return(benchmark_failed_isolated_result(
      status = "failed",
      error = paste("Failed to read isolated benchmark result:", conditionMessage(child)),
      method = method, params = params, started_at = started_at,
      baseline_bytes = baseline_bytes, peak_bytes = peak_bytes,
      poll_interval = poll_interval
    ))
  }
  child_error <- benchmark_child_result_error(child)
  if (!is.null(child_error)) {
    return(benchmark_failed_isolated_result(
      status = "failed", error = child_error,
      method = method, params = params, started_at = started_at,
      baseline_bytes = baseline_bytes, peak_bytes = peak_bytes,
      poll_interval = poll_interval
    ))
  }
  mb <- 1024^2
  child$status <- if (isTRUE(child$ok)) "success" else "failed"
  child$ok <- NULL
  child$baseline_memory_mb <- baseline_bytes / mb
  child$peak_memory_mb <- peak_bytes / mb
  child$memory_delta_mb <- max(0, peak_bytes - baseline_bytes) / mb
  child$started_at <- started_at
  child$finished_at <- Sys.time()
  child$parameters <- params
  child$backend_versions <- benchmark_backend_versions(benchmark_method_package(method))
  child$memory_method <- "sampled process-tree RSS"
  child$poll_interval <- poll_interval
  child
}

benchmark_child_result_error <- function(child) {
  if (!is.list(child)) {
    return("Isolated benchmark returned a non-list result")
  }
  if (!is.logical(child$ok) || length(child$ok) != 1L || is.na(child$ok)) {
    return("Isolated benchmark result has an invalid success flag")
  }
  if (!is.character(child$error) || length(child$error) != 1L || is.na(child$error)) {
    return("Isolated benchmark result has an invalid error field")
  }
  if (
    !is.numeric(child$runtime_s) || length(child$runtime_s) != 1L ||
      is.na(child$runtime_s) || !is.finite(child$runtime_s) || child$runtime_s < 0
  ) {
    return("Isolated benchmark result has an invalid runtime")
  }
  if (
    !is.character(child$cluster_colname) || length(child$cluster_colname) != 1L ||
      is.na(child$cluster_colname) || !nzchar(child$cluster_colname)
  ) {
    return("Isolated benchmark result has an invalid cluster column")
  }
  if (isTRUE(child$ok) && is.null(child$prediction)) {
    return("Isolated benchmark reported success without predictions")
  }
  NULL
}

benchmark_child_entry <- function(
  input_path,
  result_path,
  ready_path,
  go_path,
  done_path,
  method,
  params,
  seed,
  keep_object,
  source_path
) {
  if (!is.null(source_path)) {
    pkgload::load_all(
      source_path,
      quiet = TRUE,
      compile = FALSE,
      helpers = FALSE,
      export_all = FALSE,
      attach_testthat = FALSE
    )
  } else {
    loadNamespace("scop")
  }
  srt <- readRDS(input_path)
  file.create(ready_path)
  while (!file.exists(go_path)) Sys.sleep(0.01)
  set.seed(seed)
  start <- proc.time()[["elapsed"]]
  param_string <- get(
    "benchmark_param_string",
    envir = asNamespace("scop"), inherits = FALSE
  )
  result <- tryCatch(
    {
      execute <- get("benchmark_execute_method", envir = asNamespace("scop"), inherits = FALSE)
      value <- execute(srt = srt, method = method, params = params)
      list(
        ok = TRUE,
        error = "",
        prediction = value$prediction,
        object = if (isTRUE(keep_object)) value$object else NULL,
        cluster_colname = value$cluster_colname
      )
    },
    error = function(error) {
      list(
        ok = FALSE,
        error = conditionMessage(error),
        prediction = NULL,
        object = NULL,
        cluster_colname = param_string(
          params, "cluster_colname", paste0(method, "_cluster")
        )
      )
    }
  )
  result$runtime_s <- unname(proc.time()[["elapsed"]] - start)
  file.create(done_path)
  saveRDS(result, result_path, compress = FALSE)
  invisible(TRUE)
}

benchmark_execute_method <- function(srt, method, params) {
  function_name <- paste0("Run", method)
  fun <- benchmark_method_runner(method)
  args <- c(list(srt = srt), params)
  output <- do.call(fun, args)
  if (!inherits(output, "Seurat")) {
    log_message(function_name, " did not return a Seurat object",
      message_type = "error"
    )
  }
  cluster_colname <- benchmark_param_string(
    params, "cluster_colname", paste0(method, "_cluster")
  )
  if (!cluster_colname %in% colnames(output@meta.data)) {
    log_message(
      function_name, " did not create cluster column ", cluster_colname,
      message_type = "error"
    )
  }
  cells <- colnames(srt)
  output_cells <- rownames(output@meta.data)
  if (
    length(output_cells) != length(cells) || anyDuplicated(output_cells) ||
      !setequal(cells, output_cells)
  ) {
    log_message(
      function_name, " returned a non-identical set of spot identifiers",
      message_type = "error"
    )
  }
  prediction <- output@meta.data[cells, cluster_colname, drop = TRUE]
  prediction <- stats::setNames(as.character(prediction), cells)
  prediction[is.na(prediction) | !nzchar(prediction)] <- NA_character_
  if (all(is.na(prediction))) {
    log_message(function_name, " returned no usable cluster assignments",
      message_type = "error"
    )
  }
  list(object = output, prediction = prediction, cluster_colname = cluster_colname)
}

benchmark_process_tree_rss <- function(pid) {
  tryCatch(
    {
      root <- ps::ps_handle(pid)
      handles <- c(list(root), ps::ps_children(root, recursive = TRUE))
      values <- vapply(handles, function(handle) {
        tryCatch(as.numeric(ps::ps_memory_info(handle)[["rss"]]), error = function(e) NA_real_)
      }, numeric(1))
      if (all(is.na(values))) NA_real_ else sum(values, na.rm = TRUE)
    },
    error = function(e) NA_real_
  )
}

benchmark_kill_process_tree <- function(process) {
  try(process$kill_tree(), silent = TRUE)
  try(process$wait(timeout = 5000), silent = TRUE)
  invisible(TRUE)
}

benchmark_process_error <- function(stderr_path, stdout_path) {
  lines <- c(
    if (file.exists(stderr_path)) readLines(stderr_path, warn = FALSE) else character(),
    if (file.exists(stdout_path)) readLines(stdout_path, warn = FALSE) else character()
  )
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0L) "Isolated R process exited without a result" else utils::tail(lines, 1L)
}

benchmark_failed_isolated_result <- function(
  status,
  error,
  method,
  params,
  started_at,
  baseline_bytes = NA_real_,
  peak_bytes = NA_real_,
  poll_interval
) {
  mb <- 1024^2
  list(
    status = status,
    error = error,
    prediction = NULL,
    object = NULL,
    runtime_s = NA_real_,
    baseline_memory_mb = baseline_bytes / mb,
    peak_memory_mb = peak_bytes / mb,
    memory_delta_mb = if (is.finite(baseline_bytes) && is.finite(peak_bytes)) {
      max(0, peak_bytes - baseline_bytes) / mb
    } else {
      NA_real_
    },
    started_at = started_at,
    finished_at = Sys.time(),
    parameters = params,
    cluster_colname = benchmark_param_string(
      params, "cluster_colname", paste0(method, "_cluster")
    ),
    backend_versions = benchmark_backend_versions(benchmark_method_package(method)),
    memory_method = "sampled process-tree RSS",
    poll_interval = poll_interval
  )
}

benchmark_build_result <- function(
  srt,
  truth,
  methods,
  run_results,
  metrics,
  n_clusters,
  keep_objects,
  install_missing,
  seed,
  timeout,
  poll_interval
) {
  summary_rows <- vector("list", length(run_results))
  prediction_rows <- list()
  run_rows <- vector("list", length(run_results))
  objects <- list()
  names(summary_rows) <- names(run_results)
  names(run_rows) <- names(run_results)

  for (method in names(run_results)) {
    run <- run_results[[method]]
    values <- c(ARI = NA_real_, NMI = NA_real_, purity = NA_real_)
    n_evaluated <- 0L
    n_predicted_clusters <- NA_integer_
    if (identical(run$status, "success")) {
      prediction <- tryCatch(
        benchmark_align_prediction(run$prediction, colnames(srt)),
        error = function(error) error
      )
      if (inherits(prediction, "error")) {
        run$status <- "failed"
        run$error <- paste(
          "Prediction validation failed:", conditionMessage(prediction)
        )
      } else {
        keep <- !is.na(truth$labels) & !is.na(prediction)
        n_evaluated <- sum(keep)
        if (n_evaluated < 2L || length(unique(prediction[keep])) < 1L) {
          run$status <- "failed"
          run$error <- "Too few aligned non-missing predictions for benchmark metrics"
        } else {
          computed <- tryCatch(
            benchmark_compute_metrics(
              prediction = prediction[keep],
              truth = truth$labels[keep]
            ),
            error = function(error) error
          )
          if (inherits(computed, "error")) {
            run$status <- "failed"
            run$error <- paste("Metric computation failed:", conditionMessage(computed))
          } else {
            values <- computed
            n_predicted_clusters <- length(unique(prediction[keep]))
            prediction_rows[[method]] <- data.frame(
              spot_id = colnames(srt),
              gold_standard = unname(truth$labels[colnames(srt)]),
              method = method,
              prediction = unname(prediction),
              stringsAsFactors = FALSE
            )
            if (isTRUE(keep_objects)) objects[[method]] <- run$object
          }
        }
      }
    }
    summary_rows[[method]] <- data.frame(
      method = method,
      workflow = "SpatialDomain",
      ARI = values[["ARI"]],
      NMI = values[["NMI"]],
      purity = values[["purity"]],
      runtime_s = run$runtime_s,
      baseline_memory_mb = run$baseline_memory_mb,
      peak_memory_mb = run$peak_memory_mb,
      memory_delta_mb = run$memory_delta_mb,
      n_evaluated = as.integer(n_evaluated),
      n_clusters = as.integer(n_predicted_clusters),
      status = run$status,
      error = run$error %||% "",
      stringsAsFactors = FALSE
    )
    run_rows[[method]] <- data.frame(
      method = method,
      function_name = paste0("Run", method),
      backend_id = tolower(method),
      workflow = "SpatialDomain",
      tool_key = benchmark_param_string(run$parameters, "tool_name", method),
      cluster_colname = run$cluster_colname,
      status = run$status,
      error = run$error %||% "",
      runtime_s = run$runtime_s,
      baseline_memory_mb = run$baseline_memory_mb,
      peak_memory_mb = run$peak_memory_mb,
      memory_delta_mb = run$memory_delta_mb,
      started_at = as.character(run$started_at),
      finished_at = as.character(run$finished_at),
      backend_versions = paste(
        paste(names(run$backend_versions), run$backend_versions, sep = "="),
        collapse = ";"
      ),
      memory_method = run$memory_method,
      poll_interval = run$poll_interval,
      stringsAsFactors = FALSE
    )
    run_rows[[method]]$parameters <- I(list(run$parameters))
  }

  summary <- do.call(rbind, summary_rows)
  rownames(summary) <- NULL
  predictions <- if (length(prediction_rows) == 0L) {
    data.frame(
      spot_id = character(), gold_standard = character(),
      method = character(), prediction = character(),
      stringsAsFactors = FALSE
    )
  } else {
    out <- do.call(rbind, prediction_rows)
    rownames(out) <- NULL
    out
  }
  runs <- do.call(rbind, run_rows)
  rownames(runs) <- NULL
  metric_specs <- data.frame(
    metric = c("ARI", "NMI", "purity", "runtime_s", "peak_memory_mb"),
    direction = c("higher", "higher", "higher", "lower", "lower"),
    stringsAsFactors = FALSE
  )
  metric_rows <- lapply(seq_len(nrow(summary)), function(i) {
    data.frame(
      method = summary$method[[i]],
      metric = metric_specs$metric,
      value = as.numeric(unlist(summary[i, metric_specs$metric, drop = FALSE], use.names = FALSE)),
      workflow = "SpatialDomain",
      direction = metric_specs$direction,
      status = summary$status[[i]],
      stringsAsFactors = FALSE
    )
  })
  metric_table <- do.call(rbind, metric_rows)
  rownames(metric_table) <- NULL

  structure(
    list(
      method = "Benchmark",
      result_type = "benchmark",
      summary = summary,
      metrics = metric_table,
      predictions = predictions,
      runs = runs,
      objects = objects,
      source = list(
        gold_standard = truth$source,
        cells = colnames(srt),
        n_cells = ncol(srt)
      ),
      provenance = list(
        producer = "RunBenchmark",
        backend_id = paste(tolower(methods), collapse = ";"),
        scop_version = as.character(utils::packageVersion("scop"))
      ),
      parameters = list(
        methods = methods,
        metrics = metrics,
        n_clusters = n_clusters,
        keep_objects = keep_objects,
        install_missing = install_missing,
        seed = seed,
        timeout = timeout,
        poll_interval = poll_interval
      )
    ),
    class = c("benchmark_result", "list")
  )
}

benchmark_align_prediction <- function(prediction, cells) {
  if (
    is.null(prediction) || is.list(prediction) || is.matrix(prediction) ||
      is.data.frame(prediction) || length(prediction) == 0L
  ) {
    log_message("predictions must be a named atomic vector", message_type = "error")
  }
  ids <- names(prediction)
  if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    log_message("prediction spot identifiers must be unique and non-empty", message_type = "error")
  }
  missing <- setdiff(cells, ids)
  extra <- setdiff(ids, cells)
  if (length(missing) > 0L || length(extra) > 0L) {
    log_message("prediction spot identifiers must exactly match the benchmark input", message_type = "error")
  }
  aligned <- as.character(prediction[cells])
  aligned[is.na(aligned) | !nzchar(aligned)] <- NA_character_
  stats::setNames(aligned, cells)
}

benchmark_compute_metrics <- function(prediction, truth) {
  computed <- classification_metrics_compute(
    predicted = prediction,
    truth = truth
  )
  required <- c("ari", "nmi", "purity")
  if (!is.list(computed) && !is.atomic(computed)) {
    log_message("classification metric output must be a named vector or list", message_type = "error")
  }
  if (is.null(names(computed)) || !all(required %in% names(computed))) {
    log_message("classification metric output is missing ARI, NMI, or purity", message_type = "error")
  }
  values <- vapply(required, function(metric) {
    value <- computed[[metric]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      log_message("classification metric output contains a non-finite scalar", message_type = "error")
    }
    as.numeric(value)
  }, numeric(1))
  if (values[["ari"]] < -1 - 1e-8 || values[["ari"]] > 1 + 1e-8) {
    log_message("ARI is outside [-1, 1]", message_type = "error")
  }
  if (any(values[c("nmi", "purity")] < -1e-8 | values[c("nmi", "purity")] > 1 + 1e-8)) {
    log_message("NMI or purity is outside [0, 1]", message_type = "error")
  }
  c(ARI = values[["ari"]], NMI = values[["nmi"]], purity = values[["purity"]])
}
