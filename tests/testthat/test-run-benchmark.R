make_benchmark_seurat <- function() {
  counts <- Matrix::Matrix(
    matrix(
      seq_len(24),
      nrow = 4,
      dimnames = list(paste0("gene", 1:4), paste0("spot", 1:6))
    ),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object$gold <- c("A", "A", "B", "B", "C", "C")
  object
}

mock_benchmark_runs <- function(code) {
  testthat::local_mocked_bindings(
    benchmark_method_availability = function(adapter) {
      list(status = "available", detail = "available", versions = c(mock = "1.0"))
    },
    benchmark_package_context = function() list(source_path = NULL, libpath = .libPaths()),
    benchmark_run_isolated = function(
      input_path, run_dir, adapter, params, seed, keep_object,
      timeout, poll_interval, package_context
    ) {
      object <- readRDS(input_path)
      truth <- object$gold
      prediction <- switch(adapter$method,
        BayesSpace = c("1", "1", "2", "2", "3", "3"),
        BANKSY = c("x", "x", "y", "y", "z", NA_character_),
        SmoothClust = c("1", "2", "1", "2", "1", "2")
      )
      prediction <- stats::setNames(prediction, colnames(object))
      list(
        status = "success", error = "", prediction = prediction,
        object = if (keep_object) object else NULL,
        runtime_s = switch(adapter$method, BayesSpace = 2, BANKSY = 1, SmoothClust = 4),
        baseline_memory_mb = 100,
        peak_memory_mb = switch(adapter$method, BayesSpace = 250, BANKSY = 150, SmoothClust = 500),
        memory_delta_mb = switch(adapter$method, BayesSpace = 150, BANKSY = 50, SmoothClust = 400),
        started_at = Sys.time(), finished_at = Sys.time(),
        parameters = params,
        cluster_colname = params$cluster_colname %||% adapter$cluster_colname,
        backend_versions = c(mock = "1.0"),
        memory_method = "sampled process-tree RSS",
        poll_interval = poll_interval
      )
    }
  )
  force(code)
}

test_that("benchmark adapters cover every stable spatial domain producer", {
  stable <- spatial_method_registry()
  stable <- stable$method[
    stable$task == "domain" & stable$kind == "analysis" & stable$status == "stable"
  ]
  adapters <- benchmark_adapter_registry()
  defaults <- benchmark_resolve_adapters()
  expect_setequal(vapply(defaults, `[[`, character(1), "function_name"), stable)
  expect_true(all(stable %in% vapply(adapters, `[[`, character(1), "function_name")))
  expect_identical(benchmark_resolve_adapters("all")$GiottoCluster$tier, "legacy")
})

test_that("method aliases and cluster-count parameters resolve deterministically", {
  adapters <- benchmark_resolve_adapters(c("runbayesspace", "Banksy", "smoothclust"))
  expect_named(adapters, c("BayesSpace", "BANKSY", "SmoothClust"))

  params <- benchmark_resolve_method_params(
    list(runbayesspace = list(q = 4), banksy = list(resolution = 0.4)),
    adapters
  )
  expect_identical(params$BayesSpace$q, 4)
  expect_identical(params$BANKSY$resolution, 0.4)
  expect_identical(benchmark_complete_method_params(adapters$BayesSpace, list(), 3, 17)$q, 3L)
  expect_identical(benchmark_complete_method_params(adapters$SmoothClust, list(), 3, 17)$n_clusters, 3L)
  expect_identical(benchmark_complete_method_params(adapters$SmoothClust, list(), 3, 17)$seed, 17L)
  expect_identical(benchmark_complete_method_params(adapters$BANKSY, list(), 3, 17)$seed, 17L)
  expect_identical(benchmark_complete_method_params(adapters$BANKSY, list(seed = 9), 3, 17)$seed, 9)
  expect_false("n_clusters" %in% names(benchmark_complete_method_params(adapters$BANKSY, list(), 3, 17)))

  expect_error(benchmark_resolve_adapters("unknown"), "Unknown benchmark methods")
  expect_named(benchmark_resolve_adapters("Giotto"), "GiottoCluster")
  expect_error(benchmark_resolve_adapters(c("all", "BANKSY")), "cannot be combined")
  expect_error(
    benchmark_resolve_method_params(list(BayesSpace = list(), RunBayesSpace = list()), adapters),
    "duplicated method aliases"
  )
  duplicated_params <- list(1, 2)
  names(duplicated_params) <- c("q", "q")
  expect_error(
    benchmark_resolve_method_params(list(BayesSpace = duplicated_params), adapters),
    "uniquely named list"
  )
})

test_that("Giotto legacy adapter validates and aligns standalone cluster output", {
  object <- make_benchmark_seurat()
  adapter <- benchmark_adapter_registry()$GiottoCluster
  giotto <- structure(
    list(clusters = data.frame(
      cluster = c("C", "A", "B", "A", "C", "B"),
      row.names = rev(colnames(object)),
      stringsAsFactors = FALSE
    )),
    class = c("giotto2_cluster", "giotto2_result", "list")
  )
  testthat::local_mocked_bindings(
    RunGiottoCluster = function(srt, ...) giotto
  )
  out <- benchmark_execute_adapter(object, adapter, list())
  expect_identical(names(out$prediction), colnames(object))
  expect_identical(unname(out$prediction), c("B", "C", "A", "B", "A", "C"))
  expect_identical(out$cluster_colname, "Giotto_cluster")
})

test_that("all-unavailable dry runs do not require subprocess dependencies", {
  object <- make_benchmark_seurat()
  testthat::local_mocked_bindings(
    benchmark_method_availability = function(adapter) {
      list(status = "missing", detail = paste0(adapter$backend_id, "=missing"), versions = character())
    },
    benchmark_require_runtime = function() stop("runtime must not be checked"),
    benchmark_package_context = function() stop("package context must not be resolved")
  )
  result <- RunBenchmark(object, "gold", verbose = FALSE)
  expect_identical(result$summary$status, rep("unavailable", 3))
  expect_true(all(is.na(result$summary$ARI)))
  expect_true(all(is.na(result$summary$runtime_s)))
  expect_equal(nrow(result$predictions), 0L)
  expect_equal(length(result$objects), 0L)
  expect_s3_class(BenchmarkPlot(data = result), "patchwork")
  expect_no_error(patchwork::patchworkGrob(BenchmarkPlot(data = result)))
  expect_silent(as.data.frame(result))
  expect_output(print(result), "unavailable=3")
})

test_that("availability diagnostic errors remain method-local and truthful", {
  object <- make_benchmark_seurat()
  testthat::local_mocked_bindings(
    benchmark_method_availability = function(adapter) stop("registry diagnostic broke"),
    benchmark_require_runtime = function() stop("runtime must not be checked")
  )
  result <- RunBenchmark(object, "gold", methods = c("BayesSpace", "BANKSY"), verbose = FALSE)
  expect_identical(result$summary$status, c("unavailable", "unavailable"))
  expect_true(all(grepl("availability check failed", result$summary$error)))
})

test_that("malformed availability diagnostics are treated as unavailable", {
  adapter <- benchmark_adapter_registry()$BayesSpace
  testthat::local_mocked_bindings(
    benchmark_method_availability = function(adapter) list(status = NULL)
  )
  availability <- benchmark_method_availability_safe(adapter)
  expect_identical(availability$status, "diagnostic_error")
  expect_match(availability$detail, "invalid result")
})

test_that("install_missing enters the producer path instead of pretending unavailable", {
  object <- make_benchmark_seurat()
  runtime_checks <- 0L
  testthat::local_mocked_bindings(
    benchmark_method_availability = function(adapter) {
      list(status = "missing", detail = "backend missing", versions = character())
    },
    benchmark_require_runtime = function() {
      runtime_checks <<- runtime_checks + 1L
      invisible(TRUE)
    },
    benchmark_run_isolated = function(
      input_path, run_dir, adapter, params, seed, keep_object,
      timeout, poll_interval, package_context
    ) {
      benchmark_failed_isolated_result(
        status = "failed", error = "producer installation failed",
        adapter = adapter, params = params, started_at = Sys.time(),
        poll_interval = poll_interval
      )
    }
  )
  result <- RunBenchmark(
    object, "gold", methods = c("BayesSpace", "BANKSY"),
    install_missing = TRUE, verbose = FALSE
  )
  expect_identical(runtime_checks, 1L)
  expect_identical(result$summary$status, c("failed", "failed"))
  expect_true(all(grepl("installation failed", result$summary$error)))
})

test_that("malformed predictions fail one method without aborting the batch", {
  object <- make_benchmark_seurat()
  testthat::local_mocked_bindings(
    benchmark_method_availability = function(adapter) {
      list(status = "available", detail = "available", versions = character())
    },
    benchmark_package_context = function() list(source_path = NULL, libpath = .libPaths()),
    benchmark_run_isolated = function(
      input_path, run_dir, adapter, params, seed, keep_object,
      timeout, poll_interval, package_context
    ) {
      object <- readRDS(input_path)
      prediction <- if (adapter$method == "BayesSpace") {
        unname(as.character(object$gold))
      } else {
        stats::setNames(as.character(object$gold), colnames(object))
      }
      list(
        status = "success", error = "", prediction = prediction, object = NULL,
        runtime_s = 0.1, baseline_memory_mb = 50, peak_memory_mb = 60,
        memory_delta_mb = 10, started_at = Sys.time(), finished_at = Sys.time(),
        parameters = params, cluster_colname = adapter$cluster_colname,
        backend_versions = character(), memory_method = "sampled process-tree RSS",
        poll_interval = poll_interval
      )
    }
  )
  result <- RunBenchmark(
    object, "gold", methods = c("BayesSpace", "BANKSY"), verbose = FALSE
  )
  expect_identical(result$summary$status, c("failed", "success"))
  expect_match(result$summary$error[[1]], "Prediction validation failed")
  expect_true(all(is.na(result$summary[1, c("ARI", "NMI", "purity")])))
  expect_equal(result$summary$ARI[[2]], 1)
  expect_identical(unique(result$predictions$method), "BANKSY")
})

test_that("prediction alignment rejects missing extra and duplicated spot IDs", {
  object <- make_benchmark_seurat()
  cells <- colnames(object)
  prediction <- stats::setNames(as.character(object$gold), cells)
  expect_identical(benchmark_align_prediction(prediction[rev(cells)], cells), prediction)
  expect_error(benchmark_align_prediction(unname(prediction), cells), "unique and non-empty")
  expect_error(benchmark_align_prediction(prediction[-1], cells), "exactly match")
  extra <- c(prediction, extra = "A")
  expect_error(benchmark_align_prediction(extra, cells), "exactly match")
  duplicated <- prediction
  names(duplicated)[[2]] <- names(duplicated)[[1]]
  expect_error(benchmark_align_prediction(duplicated, cells), "unique and non-empty")
})

test_that("metric computation handles permutation partial labels and unequal class counts", {
  perfect <- benchmark_compute_metrics(
    prediction = c("x", "x", "y", "y"),
    truth = c("A", "A", "B", "B")
  )
  expect_equal(perfect, c(ARI = 1, NMI = 1, purity = 1))
  unequal <- benchmark_compute_metrics(
    prediction = c("x", "x", "x", "y", "y", "y"),
    truth = c("A", "A", "B", "B", "C", "C")
  )
  expect_true(all(is.finite(unequal)))
  expect_true(unequal[["ARI"]] < 1)
  partial <- benchmark_align_prediction(
    stats::setNames(c("x", "x", NA, "y", "z", "z"), paste0("spot", 1:6)),
    paste0("spot", 1:6)
  )
  keep <- !is.na(partial)
  partial_metrics <- benchmark_compute_metrics(partial[keep], c("A", "A", "B", "C", "C"))
  expect_true(all(is.finite(partial_metrics)))
})

test_that("numeric controls reject overflow and non-finite values", {
  object <- make_benchmark_seurat()
  expect_error(RunBenchmark(object, "gold", seed = .Machine$integer.max + 1), "invalid value")
  expect_error(RunBenchmark(object, "gold", n_clusters = Inf), "invalid value")
  expect_error(RunBenchmark(object, "gold", timeout = -1), "invalid value")
  expect_error(RunBenchmark(object, "gold", poll_interval = 0), "invalid value")
})

test_that("RunBenchmark aligns spots and retains metrics and resource measurements", {
  object <- make_benchmark_seurat()
  mock_benchmark_runs({
    result <- RunBenchmark(
      object,
      gold_standard = "gold",
      methods = c("BayesSpace", "BANKSY"),
      keep_objects = TRUE,
      verbose = FALSE
    )
  })

  expect_s3_class(result, "scop_benchmark")
  expect_identical(result$summary$method, c("BayesSpace", "BANKSY"))
  expect_equal(result$summary$ARI[[1]], 1)
  expect_equal(result$summary$NMI[[1]], 1)
  expect_equal(result$summary$purity[[1]], 1)
  expect_equal(result$summary$n_evaluated, c(6L, 5L))
  expect_equal(result$summary$runtime_s, c(2, 1))
  expect_equal(result$summary$peak_memory_mb, c(250, 150))
  expect_identical(result$predictions$spot_id[1:6], colnames(object))
  expect_named(result$objects, c("BayesSpace", "BANKSY"))
  expect_s3_class(as.data.frame(result), "data.frame")
  expect_setequal(unique(result$metrics$direction), c("higher", "lower"))
  expect_true(all(c("parameters", "backend_versions", "cluster_colname") %in% names(result$runs)))
  expect_identical(result$runs$parameters[[2]]$seed, 1L)
})

test_that("benchmark execution never mutates the caller's Seurat object", {
  object <- make_benchmark_seurat()
  original <- object
  mock_benchmark_runs({
    invisible(RunBenchmark(
      object, "gold", methods = c("BayesSpace", "BANKSY"), verbose = FALSE
    ))
  })
  expect_identical(object, original)
})

test_that("custom producer result keys are recorded truthfully", {
  object <- make_benchmark_seurat()
  mock_benchmark_runs({
    result <- RunBenchmark(
      object, "gold", methods = "BANKSY",
      method_params = list(BANKSY = list(
        tool_name = "BANKSY_custom",
        cluster_colname = "BANKSY_custom_cluster"
      )),
      verbose = FALSE
    )
  })
  expect_identical(result$runs$tool_key, "BANKSY_custom")
  expect_identical(result$runs$cluster_colname, "BANKSY_custom_cluster")
  expect_identical(result$runs$parameters[[1]]$tool_name, "BANKSY_custom")
})

test_that("invalid backend-specific names cannot corrupt benchmark bookkeeping", {
  adapter <- benchmark_adapter_registry()$BANKSY
  failed <- benchmark_failed_isolated_result(
    status = "failed", error = "producer rejected arguments",
    adapter = adapter,
    params = list(cluster_colname = c("a", "b"), tool_name = c("x", "y")),
    started_at = Sys.time(), poll_interval = 0.1
  )
  expect_identical(failed$cluster_colname, adapter$cluster_colname)
  expect_identical(
    benchmark_param_string(failed$parameters, "tool_name", adapter$tool_key),
    adapter$tool_key
  )
})

test_that("gold-standard vectors require exact unique spot identifiers", {
  object <- make_benchmark_seurat()
  truth <- stats::setNames(object$gold, colnames(object))
  expect_identical(benchmark_resolve_truth(object, truth)$labels, truth)
  expect_error(benchmark_resolve_truth(object, unname(truth)), "unique non-empty spot names")
  names(truth)[[1]] <- names(truth)[[2]]
  expect_error(benchmark_resolve_truth(object, truth), "unique non-empty spot names")
  expect_error(benchmark_resolve_truth(object, stats::setNames(rep("A", 6), colnames(object))), "two non-missing classes")
})

test_that("failed and unavailable methods do not receive pseudo metrics", {
  object <- make_benchmark_seurat()
  testthat::local_mocked_bindings(
    benchmark_method_availability = function(adapter) {
      if (adapter$method == "BANKSY") {
        list(status = "missing", detail = "banksy=missing", versions = character())
      } else {
        list(status = "available", detail = "available", versions = character())
      }
    },
    benchmark_package_context = function() list(source_path = NULL, libpath = .libPaths()),
    benchmark_run_isolated = function(
      input_path, run_dir, adapter, params, seed, keep_object,
      timeout, poll_interval, package_context
    ) {
      list(
        status = "failed", error = "backend error", prediction = NULL, object = NULL,
        runtime_s = NA_real_, baseline_memory_mb = 80, peak_memory_mb = 90,
        memory_delta_mb = 10, started_at = Sys.time(), finished_at = Sys.time(),
        parameters = params, cluster_colname = adapter$cluster_colname,
        backend_versions = character(), memory_method = "sampled process-tree RSS",
        poll_interval = poll_interval
      )
    }
  )
  result <- RunBenchmark(
    object, "gold", methods = c("BayesSpace", "BANKSY"), verbose = FALSE
  )
  expect_identical(result$summary$status, c("failed", "unavailable"))
  expect_true(all(is.na(result$summary$ARI)))
  expect_true(all(is.na(result$summary$NMI)))
  expect_equal(nrow(result$predictions), 0L)
})

test_that("timeout status remains method-local in a mixed batch", {
  object <- make_benchmark_seurat()
  testthat::local_mocked_bindings(
    benchmark_method_availability = function(adapter) {
      list(status = "available", detail = "available", versions = character())
    },
    benchmark_package_context = function() list(source_path = NULL, libpath = .libPaths()),
    benchmark_run_isolated = function(
      input_path, run_dir, adapter, params, seed, keep_object,
      timeout, poll_interval, package_context
    ) {
      if (adapter$method == "BayesSpace") {
        return(benchmark_failed_isolated_result(
          status = "timeout", error = "method exceeded timeout",
          adapter = adapter, params = params, started_at = Sys.time(),
          poll_interval = poll_interval
        ))
      }
      object <- readRDS(input_path)
      list(
        status = "success", error = "",
        prediction = stats::setNames(as.character(object$gold), colnames(object)),
        object = NULL, runtime_s = 0.1, baseline_memory_mb = 50,
        peak_memory_mb = 60, memory_delta_mb = 10,
        started_at = Sys.time(), finished_at = Sys.time(), parameters = params,
        cluster_colname = adapter$cluster_colname, backend_versions = character(),
        memory_method = "sampled process-tree RSS", poll_interval = poll_interval
      )
    }
  )
  result <- RunBenchmark(
    object, "gold", methods = c("BayesSpace", "BANKSY"), verbose = FALSE
  )
  expect_identical(result$summary$status, c("timeout", "success"))
  expect_true(is.na(result$summary$ARI[[1]]))
  expect_equal(result$summary$ARI[[2]], 1)
  expect_identical(unique(result$predictions$method), "BANKSY")
})

test_that("adapter execution requires exact spot identity and expected output", {
  object <- make_benchmark_seurat()
  adapter <- benchmark_adapter_registry()$SmoothClust
  valid <- object[, rev(colnames(object))]
  valid$SmoothClust_cluster <- rev(c("1", "1", "2", "2", "3", "3"))
  testthat::local_mocked_bindings(
    RunSmoothClust = function(srt, ...) valid
  )
  aligned <- benchmark_execute_adapter(object, adapter, list())
  expect_identical(names(aligned$prediction), colnames(object))

  testthat::local_mocked_bindings(
    RunSmoothClust = function(srt, ...) "not Seurat"
  )
  expect_error(benchmark_execute_adapter(object, adapter, list()), "did not return a Seurat")

  missing_column <- object
  testthat::local_mocked_bindings(
    RunSmoothClust = function(srt, ...) missing_column
  )
  expect_error(benchmark_execute_adapter(object, adapter, list()), "did not create cluster column")

  missing_spot <- object[, -1]
  missing_spot$SmoothClust_cluster <- "1"
  testthat::local_mocked_bindings(
    RunSmoothClust = function(srt, ...) missing_spot
  )
  expect_error(benchmark_execute_adapter(object, adapter, list()), "non-identical set")
})

test_that("process-tree RSS is measured and supervised children are terminated", {
  expect_gt(benchmark_process_tree_rss(Sys.getpid()), 0)
  process <- callr::r_bg(function() Sys.sleep(30), supervise = TRUE)
  pid <- process$get_pid()
  expect_true(ps::ps_is_running(ps::ps_handle(pid)))
  benchmark_kill_process_tree(process)
  expect_false(process$is_alive())
})

benchmark_test_child <- function(
  input_path, result_path, ready_path, go_path, done_path,
  adapter, params, seed, keep_object, source_path
) {
  object <- readRDS(input_path)
  file.create(ready_path)
  while (!file.exists(go_path)) Sys.sleep(0.01)
  payload <- numeric(params$allocate)
  payload[[1]] <- 1
  Sys.sleep(params$sleep)
  result <- list(
    ok = TRUE, error = "",
    prediction = stats::setNames(as.character(object$gold), colnames(object)),
    object = NULL, runtime_s = params$sleep,
    cluster_colname = adapter$cluster_colname
  )
  file.create(done_path)
  saveRDS(result, result_path, compress = FALSE)
}

benchmark_invalid_child <- function(
  input_path, result_path, ready_path, go_path, done_path,
  adapter, params, seed, keep_object, source_path
) {
  file.create(ready_path)
  while (!file.exists(go_path)) Sys.sleep(0.01)
  saveRDS("not a benchmark result", result_path)
  file.create(done_path)
}

benchmark_corrupt_child <- function(
  input_path, result_path, ready_path, go_path, done_path,
  adapter, params, seed, keep_object, source_path
) {
  file.create(ready_path)
  while (!file.exists(go_path)) Sys.sleep(0.01)
  writeBin(charToRaw("not-an-rds"), result_path)
  file.create(done_path)
}

benchmark_crash_child <- function(
  input_path, result_path, ready_path, go_path, done_path,
  adapter, params, seed, keep_object, source_path
) {
  stop("child crashed before ready")
}

test_that("isolated monitor records runtime memory and timeout status", {
  object <- make_benchmark_seurat()
  input_path <- tempfile(fileext = ".rds")
  saveRDS(object, input_path)
  adapter <- benchmark_adapter_registry()$SmoothClust
  context <- list(source_path = NULL, libpath = .libPaths())
  success_dir <- tempfile("benchmark-success-")
  dir.create(success_dir)
  success <- benchmark_run_isolated(
    input_path, success_dir, adapter,
    params = list(allocate = 3e6, sleep = 0.25), seed = 1,
    keep_object = FALSE, timeout = 5, poll_interval = 0.02,
    package_context = context, child_entry = benchmark_test_child
  )
  expect_identical(success$status, "success")
  expect_gte(success$peak_memory_mb, success$baseline_memory_mb)
  expect_gte(success$memory_delta_mb, 0)
  expect_equal(success$runtime_s, 0.25)

  timeout_dir <- tempfile("benchmark-timeout-")
  dir.create(timeout_dir)
  timed <- benchmark_run_isolated(
    input_path, timeout_dir, adapter,
    params = list(allocate = 1, sleep = 30), seed = 1,
    keep_object = FALSE, timeout = 0.2, poll_interval = 0.02,
    package_context = context, child_entry = benchmark_test_child
  )
  expect_identical(timed$status, "timeout")
})

test_that("invalid corrupt and crashed child results become failed runs", {
  object <- make_benchmark_seurat()
  input_path <- tempfile(fileext = ".rds")
  saveRDS(object, input_path)
  adapter <- benchmark_adapter_registry()$SmoothClust
  context <- list(source_path = NULL, libpath = .libPaths())
  children <- list(
    invalid = benchmark_invalid_child,
    corrupt = benchmark_corrupt_child,
    crashed = benchmark_crash_child
  )
  results <- lapply(names(children), function(name) {
    run_dir <- tempfile(paste0("benchmark-", name, "-"))
    dir.create(run_dir)
    benchmark_run_isolated(
      input_path, run_dir, adapter,
      params = list(), seed = 1, keep_object = FALSE,
      timeout = 5, poll_interval = 0.02, package_context = context,
      child_entry = children[[name]]
    )
  })
  expect_identical(vapply(results, `[[`, character(1), "status"), rep("failed", 3))
  expect_match(results[[1]]$error, "non-list result")
  expect_match(results[[2]]$error, "Failed to read isolated benchmark result")
  expect_match(results[[3]]$error, "child crashed before ready")
})

test_that("real child error handling resolves internal helpers from the namespace", {
  object <- make_benchmark_seurat()
  input_path <- tempfile(fileext = ".rds")
  saveRDS(object, input_path)
  adapter <- list(
    method = "MissingProducer",
    function_name = "RunProducerThatDoesNotExist",
    backend_id = "core",
    cluster_colname = "missing_cluster"
  )
  run_dir <- tempfile("benchmark-real-child-error-")
  dir.create(run_dir)
  result <- benchmark_run_isolated(
    input_path, run_dir, adapter,
    params = list(), seed = 1, keep_object = FALSE,
    timeout = 30, poll_interval = 0.05,
    package_context = benchmark_package_context()
  )
  expect_identical(result$status, "failed")
  expect_match(result$error, "RunProducerThatDoesNotExist")
  expect_false(grepl("benchmark_param_string", result$error, fixed = TRUE))
  expect_identical(result$cluster_colname, "missing_cluster")
})

test_that("repeated timeouts leave no new R subprocesses", {
  object <- make_benchmark_seurat()
  input_path <- tempfile(fileext = ".rds")
  saveRDS(object, input_path)
  adapter <- benchmark_adapter_registry()$SmoothClust
  context <- list(source_path = NULL, libpath = .libPaths())
  r_processes <- function() {
    processes <- ps::ps()
    unique(processes$pid[grepl("^(r|rterm|rscript)(\\.exe)?$", processes$name, ignore.case = TRUE)])
  }
  before <- r_processes()
  results <- lapply(seq_len(3), function(i) {
    run_dir <- tempfile("benchmark-repeat-timeout-")
    dir.create(run_dir)
    benchmark_run_isolated(
      input_path, run_dir, adapter,
      params = list(allocate = 1, sleep = 30), seed = i,
      keep_object = FALSE, timeout = 0.15, poll_interval = 0.02,
      package_context = context, child_entry = benchmark_test_child
    )
  })
  Sys.sleep(0.2)
  after <- r_processes()
  expect_identical(vapply(results, `[[`, character(1), "status"), rep("timeout", 3))
  expect_length(setdiff(after, before), 0L)
})
