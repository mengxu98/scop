#' @title Benchmark spatial domain clustering methods
#'
#' @description
#' Run registered spatial domain clustering methods from the same immutable
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
#' @param methods Spatial domain methods to benchmark. `NULL` uses every stable
#' registered domain producer. Use `"all"` to additionally include supported
#' legacy baselines such as GiottoCluster, or select them explicitly. Method
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
#' result, including standalone legacy results when explicitly selected. The
#' default keeps only aligned labels and compact run metadata.
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
#' @return A `scop_benchmark` object. Use `as.data.frame()` for its summary,
#' `plot()` or [BenchmarkPlot()] for visualization, and `$predictions` for the
#' aligned labels.
#' @concept spatial-producer
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
    seed, "seed", lower = 0, upper = .Machine$integer.max, integer = TRUE
  )
  benchmark_assert_number(timeout, "timeout", lower = 0, allow_infinite = TRUE)
  benchmark_assert_number(poll_interval, "poll_interval", lower = 0, strict = TRUE)

  truth <- benchmark_resolve_truth(srt, gold_standard)
  n_truth_clusters <- length(unique(truth$labels[!is.na(truth$labels)]))
  if (is.null(n_clusters)) {
    n_clusters <- n_truth_clusters
  } else {
    benchmark_assert_number(
      n_clusters, "n_clusters", lower = 2,
      upper = .Machine$integer.max, integer = TRUE
    )
    n_clusters <- as.integer(n_clusters)
  }
  metrics <- benchmark_resolve_metrics(metrics)
  adapters <- benchmark_resolve_adapters(methods)
  method_params <- benchmark_resolve_method_params(method_params, adapters)
  prepared_params <- lapply(names(adapters), function(method) {
    benchmark_complete_method_params(
      adapter = adapters[[method]],
      params = method_params[[method]] %||% list(),
      n_clusters = n_clusters,
      seed = seed
    )
  })
  names(prepared_params) <- names(adapters)
  availability_results <- lapply(adapters, benchmark_method_availability_safe)
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

  run_results <- vector("list", length(adapters))
  names(run_results) <- names(adapters)

  for (method in names(adapters)) {
    adapter <- adapters[[method]]
    params <- prepared_params[[method]]
    availability <- availability_results[[method]]
    if (!will_execute[[method]]) {
      run_results[[method]] <- benchmark_unavailable_result(
        adapter = adapter,
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
      adapter = adapter,
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
    adapters = adapters,
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

#' @export
print.scop_benchmark <- function(x, ...) {
  cat("<scop_benchmark>\n")
  cat("  workflow: spatial domain clustering\n")
  cat("  methods: ", nrow(x$summary), "\n", sep = "")
  status_text <- paste(
    names(table(x$summary$status)),
    as.integer(table(x$summary$status)),
    sep = "=",
    collapse = ", "
  )
  cat("  status: ", status_text, "\n", sep = "")
  print(x$summary, row.names = FALSE)
  invisible(x)
}

#' @export
as.data.frame.scop_benchmark <- function(x, row.names = NULL, optional = FALSE, ...) {
  as.data.frame(x$summary, row.names = row.names, optional = optional, ...)
}

#' @export
plot.scop_benchmark <- function(x, ...) {
  BenchmarkPlot(data = x, plot_type = "overview", ...)
}

benchmark_adapter_registry <- function() {
  list(
    BayesSpace = list(
      method = "BayesSpace",
      function_name = "RunBayesSpace",
      aliases = c("BayesSpace", "RunBayesSpace"),
      backend_id = "bayesspace",
      tool_key = "BayesSpace",
      cluster_colname = "BayesSpace_cluster",
      cluster_arg = "q",
      seed_arg = NA_character_
    ),
    BANKSY = list(
      method = "BANKSY",
      function_name = "RunBANKSY",
      aliases = c("BANKSY", "Banksy", "RunBANKSY", "RunBanksy"),
      backend_id = "banksy",
      tool_key = "BANKSY",
      cluster_colname = "BANKSY_cluster",
      cluster_arg = NA_character_,
      seed_arg = "seed"
    ),
    SmoothClust = list(
      method = "SmoothClust",
      function_name = "RunSmoothClust",
      aliases = c("SmoothClust", "smoothclust", "RunSmoothClust"),
      backend_id = "smoothclust",
      tool_key = "SmoothClust",
      cluster_colname = "SmoothClust_cluster",
      cluster_arg = "n_clusters",
      seed_arg = "seed",
      workflow = "SpatialDomain",
      tier = "stable",
      output_type = "seurat"
    ),
    GiottoCluster = list(
      method = "GiottoCluster",
      function_name = "RunGiottoCluster",
      aliases = c("GiottoCluster", "Giotto", "RunGiottoCluster"),
      backend_id = "giotto",
      tool_key = "GiottoCluster",
      cluster_colname = "Giotto_cluster",
      cluster_arg = NA_character_,
      seed_arg = "seed",
      workflow = "SpatialDomainLegacy",
      tier = "legacy",
      output_type = "giotto_cluster"
    )
  )
}

benchmark_resolve_adapters <- function(methods = NULL) {
  registry <- benchmark_adapter_registry()
  stable_domain <- spatial_method_registry()
  stable_domain <- stable_domain[
    stable_domain$task == "domain" &
      stable_domain$kind == "analysis" &
      stable_domain$status == "stable",
    "method",
    drop = TRUE
  ]
  available_adapters <- vapply(registry, `[[`, character(1), "function_name")
  missing_adapters <- setdiff(stable_domain, available_adapters)
  if (length(missing_adapters) > 0L) {
    log_message(
      "Stable spatial domain producers lack benchmark adapters: {.val {missing_adapters}}",
      message_type = "error"
    )
  }
  default_registry <- registry[available_adapters %in% stable_domain]
  if (is.null(methods)) return(default_registry)
  if (!is.character(methods) || length(methods) == 0L || anyNA(methods)) {
    log_message("{.arg methods} must contain one or more method names", message_type = "error")
  }
  if (any(tolower(methods) == "all")) {
    if (length(methods) != 1L) {
      log_message("{.val all} cannot be combined with other {.arg methods}", message_type = "error")
    }
    return(registry)
  }
  alias_map <- do.call(c, unname(lapply(registry, function(x) {
    stats::setNames(rep(x$method, length(x$aliases)), tolower(x$aliases))
  })))
  resolved <- unname(alias_map[tolower(methods)])
  if (anyNA(resolved)) {
    unknown <- methods[is.na(resolved)]
    log_message(
      "Unknown benchmark methods: {.val {unknown}}. Available methods are {.val {names(registry)}}",
      message_type = "error"
    )
  }
  resolved <- unique(resolved)
  registry[resolved]
}

benchmark_resolve_method_params <- function(method_params, adapters) {
  if (!is.list(method_params)) {
    log_message("{.arg method_params} must be a named list", message_type = "error")
  }
  if (length(method_params) == 0L) {
    return(stats::setNames(rep(list(list()), length(adapters)), names(adapters)))
  }
  if (is.null(names(method_params)) || any(!nzchar(names(method_params)))) {
    log_message("{.arg method_params} must be named by method", message_type = "error")
  }
  all_adapters <- benchmark_adapter_registry()
  alias_map <- do.call(c, unname(lapply(all_adapters, function(x) {
    stats::setNames(rep(x$method, length(x$aliases)), tolower(x$aliases))
  })))
  resolved_names <- unname(alias_map[tolower(names(method_params))])
  if (anyNA(resolved_names)) {
    log_message(
      "Unknown methods in {.arg method_params}: {.val {names(method_params)[is.na(resolved_names)]}}",
      message_type = "error"
    )
  }
  if (anyDuplicated(resolved_names)) {
    log_message("{.arg method_params} contains duplicated method aliases", message_type = "error")
  }
  outside <- setdiff(resolved_names, names(adapters))
  if (length(outside) > 0L) {
    log_message(
      "{.arg method_params} includes methods not selected by {.arg methods}: {.val {outside}}",
      message_type = "error"
    )
  }
  out <- stats::setNames(rep(list(list()), length(adapters)), names(adapters))
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
        "Parameters for {.val {names(method_params)[[i]]}} must be a uniquely named list",
        message_type = "error"
      )
    }
    if ("srt" %in% names(value)) {
      log_message("{.arg method_params} cannot replace the benchmark input {.arg srt}", message_type = "error")
    }
    out[[resolved_names[[i]]]] <- value
  }
  out
}

benchmark_complete_method_params <- function(adapter, params, n_clusters, seed) {
  cluster_arg <- adapter$cluster_arg
  if (!is.na(cluster_arg) && !cluster_arg %in% names(params)) {
    params[[cluster_arg]] <- as.integer(n_clusters)
  }
  seed_arg <- adapter$seed_arg
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
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    log_message("{.arg {arg}} must be TRUE or FALSE", message_type = "error")
  }
  invisible(TRUE)
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

benchmark_method_availability <- function(adapter) {
  status <- SpatialBackendStatus(
    method = adapter$function_name,
    api_check = TRUE,
    refresh = TRUE
  )
  if (nrow(status) == 0L) {
    return(list(status = "missing", detail = "backend is not registered", versions = character()))
  }
  available <- all(status$availability == "available")
  detail <- paste(
    paste0(status$backend_id, "=", status$availability),
    collapse = "; "
  )
  versions <- benchmark_backend_versions(adapter$backend_id)
  list(
    status = if (available) "available" else paste(unique(status$availability), collapse = ","),
    detail = detail,
    versions = versions
  )
}

benchmark_method_availability_safe <- function(adapter) {
  tryCatch(
    {
      availability <- benchmark_method_availability(adapter)
      if (
        !is.list(availability) ||
          !is.character(availability$status) || length(availability$status) != 1L ||
          is.na(availability$status) || !nzchar(availability$status) ||
          !is.character(availability$detail) || length(availability$detail) != 1L ||
          is.na(availability$detail)
      ) {
        stop("backend availability diagnostic returned an invalid result", call. = FALSE)
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

benchmark_backend_versions <- function(backend_id) {
  versions <- tryCatch(
    spatial_result_backend_versions(backend_id),
    error = function(error) character()
  )
  if (is.null(versions) || !is.atomic(versions)) return(character())
  out <- as.character(versions)
  names(out) <- names(versions)
  out
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

benchmark_unavailable_result <- function(adapter, params, availability, poll_interval) {
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
      params, "cluster_colname", adapter$cluster_colname
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
  adapter,
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
        adapter = adapter,
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
      adapter = adapter, params = params, started_at = started_at,
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
      adapter = adapter, params = params, started_at = started_at,
      poll_interval = poll_interval
    ))
  }
  if (!file.exists(ready_path)) {
    process$wait(timeout = 5000)
    return(benchmark_failed_isolated_result(
      status = "failed", error = benchmark_process_error(stderr_path, stdout_path),
      adapter = adapter, params = params, started_at = started_at,
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
      adapter = adapter, params = params, started_at = started_at,
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
      adapter = adapter, params = params, started_at = started_at,
      baseline_bytes = baseline_bytes, peak_bytes = peak_bytes,
      poll_interval = poll_interval
    ))
  }
  process$wait(timeout = 5000)
  if (!file.exists(result_path)) {
    return(benchmark_failed_isolated_result(
      status = "failed", error = benchmark_process_error(stderr_path, stdout_path),
      adapter = adapter, params = params, started_at = started_at,
      baseline_bytes = baseline_bytes, peak_bytes = peak_bytes,
      poll_interval = poll_interval
    ))
  }
  child <- tryCatch(readRDS(result_path), error = function(error) error)
  if (inherits(child, "error")) {
    return(benchmark_failed_isolated_result(
      status = "failed",
      error = paste("Failed to read isolated benchmark result:", conditionMessage(child)),
      adapter = adapter, params = params, started_at = started_at,
      baseline_bytes = baseline_bytes, peak_bytes = peak_bytes,
      poll_interval = poll_interval
    ))
  }
  child_error <- benchmark_child_result_error(child)
  if (!is.null(child_error)) {
    return(benchmark_failed_isolated_result(
      status = "failed", error = child_error,
      adapter = adapter, params = params, started_at = started_at,
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
  child$backend_versions <- benchmark_backend_versions(adapter$backend_id)
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
  adapter,
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
    "benchmark_param_string", envir = asNamespace("scop"), inherits = FALSE
  )
  result <- tryCatch(
    {
      execute <- get("benchmark_execute_adapter", envir = asNamespace("scop"), inherits = FALSE)
      value <- execute(srt = srt, adapter = adapter, params = params)
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
          params, "cluster_colname", adapter$cluster_colname
        )
      )
    }
  )
  result$runtime_s <- unname(proc.time()[["elapsed"]] - start)
  file.create(done_path)
  saveRDS(result, result_path, compress = FALSE)
  invisible(TRUE)
}

benchmark_execute_adapter <- function(srt, adapter, params) {
  fun <- get(adapter$function_name, envir = asNamespace("scop"), inherits = FALSE)
  args <- c(list(srt = srt), params)
  output <- do.call(fun, args)
  if (identical(adapter$output_type, "giotto_cluster")) {
    if (
      !inherits(output, "giotto2_cluster") ||
        !is.data.frame(output$clusters) ||
        !"cluster" %in% colnames(output$clusters)
    ) {
      stop(adapter$function_name, " did not return a valid giotto2_cluster result", call. = FALSE)
    }
    cells <- colnames(srt)
    output_cells <- rownames(output$clusters)
    if (
      is.null(output_cells) || length(output_cells) != length(cells) ||
        anyDuplicated(output_cells) || !setequal(cells, output_cells)
    ) {
      stop(
        adapter$function_name, " returned a non-identical set of spot identifiers",
        call. = FALSE
      )
    }
    prediction <- output$clusters[cells, "cluster", drop = TRUE]
    prediction <- stats::setNames(as.character(prediction), cells)
    prediction[is.na(prediction) | !nzchar(prediction)] <- NA_character_
    if (all(is.na(prediction))) {
      stop(adapter$function_name, " returned no usable cluster assignments", call. = FALSE)
    }
    return(list(
      object = output,
      prediction = prediction,
      cluster_colname = benchmark_param_string(
        params, "cluster_colname", adapter$cluster_colname
      )
    ))
  }
  if (!inherits(output, "Seurat")) {
    stop(adapter$function_name, " did not return a Seurat object", call. = FALSE)
  }
  cluster_colname <- benchmark_param_string(
    params, "cluster_colname", adapter$cluster_colname
  )
  if (!cluster_colname %in% colnames(output@meta.data)) {
    stop(
      adapter$function_name, " did not create cluster column ", cluster_colname,
      call. = FALSE
    )
  }
  cells <- colnames(srt)
  output_cells <- rownames(output@meta.data)
  if (
    length(output_cells) != length(cells) || anyDuplicated(output_cells) ||
      !setequal(cells, output_cells)
  ) {
    stop(
      adapter$function_name, " returned a non-identical set of spot identifiers",
      call. = FALSE
    )
  }
  prediction <- output@meta.data[cells, cluster_colname, drop = TRUE]
  prediction <- stats::setNames(as.character(prediction), cells)
  prediction[is.na(prediction) | !nzchar(prediction)] <- NA_character_
  if (all(is.na(prediction))) {
    stop(adapter$function_name, " returned no usable cluster assignments", call. = FALSE)
  }
  list(object = output, prediction = prediction, cluster_colname = cluster_colname)
}

benchmark_process_tree_rss <- function(pid) {
  tryCatch({
    root <- ps::ps_handle(pid)
    handles <- c(list(root), ps::ps_children(root, recursive = TRUE))
    values <- vapply(handles, function(handle) {
      tryCatch(as.numeric(ps::ps_memory_info(handle)[["rss"]]), error = function(e) NA_real_)
    }, numeric(1))
    if (all(is.na(values))) NA_real_ else sum(values, na.rm = TRUE)
  }, error = function(e) NA_real_)
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
  adapter,
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
      params, "cluster_colname", adapter$cluster_colname
    ),
    backend_versions = benchmark_backend_versions(adapter$backend_id),
    memory_method = "sampled process-tree RSS",
    poll_interval = poll_interval
  )
}

benchmark_build_result <- function(
  srt,
  truth,
  adapters,
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
      workflow = adapters[[method]]$workflow %||% "SpatialDomain",
      tier = adapters[[method]]$tier %||% "stable",
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
    adapter <- adapters[[method]]
    run_rows[[method]] <- data.frame(
      method = method,
      function_name = adapter$function_name,
      backend_id = adapter$backend_id,
      workflow = adapter$workflow %||% "SpatialDomain",
      tier = adapter$tier %||% "stable",
      tool_key = benchmark_param_string(run$parameters, "tool_name", adapter$tool_key),
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
      workflow = adapters[[summary$method[[i]]]]$workflow %||% "SpatialDomain",
      direction = metric_specs$direction,
      status = summary$status[[i]],
      stringsAsFactors = FALSE
    )
  })
  metric_table <- do.call(rbind, metric_rows)
  rownames(metric_table) <- NULL

  structure(
    list(
      schema_version = "1.0.0",
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
        backend_id = paste(vapply(adapters, `[[`, character(1), "backend_id"), collapse = ";"),
        scop_version = as.character(utils::packageVersion("scop"))
      ),
      parameters = list(
        methods = names(adapters),
        metrics = metrics,
        n_clusters = n_clusters,
        keep_objects = keep_objects,
        install_missing = install_missing,
        seed = seed,
        timeout = timeout,
        poll_interval = poll_interval
      )
    ),
    class = c("scop_benchmark", "list")
  )
}

benchmark_align_prediction <- function(prediction, cells) {
  if (
    is.null(prediction) || is.list(prediction) || is.matrix(prediction) ||
      is.data.frame(prediction) || length(prediction) == 0L
  ) {
    stop("predictions must be a named atomic vector", call. = FALSE)
  }
  ids <- names(prediction)
  if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("prediction spot identifiers must be unique and non-empty", call. = FALSE)
  }
  missing <- setdiff(cells, ids)
  extra <- setdiff(ids, cells)
  if (length(missing) > 0L || length(extra) > 0L) {
    stop("prediction spot identifiers must exactly match the benchmark input", call. = FALSE)
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
    stop("classification metric output must be a named vector or list", call. = FALSE)
  }
  if (is.null(names(computed)) || !all(required %in% names(computed))) {
    stop("classification metric output is missing ARI, NMI, or purity", call. = FALSE)
  }
  values <- vapply(required, function(metric) {
    value <- computed[[metric]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      stop("classification metric output contains a non-finite scalar", call. = FALSE)
    }
    as.numeric(value)
  }, numeric(1))
  if (values[["ari"]] < -1 - 1e-8 || values[["ari"]] > 1 + 1e-8) {
    stop("ARI is outside [-1, 1]", call. = FALSE)
  }
  if (any(values[c("nmi", "purity")] < -1e-8 | values[c("nmi", "purity")] > 1 + 1e-8)) {
    stop("NMI or purity is outside [0, 1]", call. = FALSE)
  }
  c(ARI = values[["ari"]], NMI = values[["nmi"]], purity = values[["purity"]])
}
