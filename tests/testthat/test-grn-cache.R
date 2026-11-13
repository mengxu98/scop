test_that("completed GRN cache reuses validated results and changes identity with inputs", {
  directory <- withr::local_tempdir()
  input <- list(expression = matrix(1:12, nrow = 3L), cells = c("c1", "c2", "c3", "c4"), groups = list(all = c("c1", "c2", "c3", "c4")))
  config <- list(method = "test", work_dir = directory, runtime = list(version = "test"), method_args = list(alpha = 1), resources = list(), seed = 3L)
  runs <- 0L
  compute <- function() {
    runs <<- runs + 1L
    structure(list(schema_version = "1.0", execution = list(status = "completed", artifacts = list()), groups = list(all = list(native = data.frame()))), class = c("GRNResult", "list"))
  }
  first <- grn_cached_run(input, config, compute)
  expect_false(first$execution$cache$hit)
  second <- grn_cached_run(input, config, compute)
  expect_true(second$execution$cache$hit)
  expect_identical(runs, 1L)
  changed <- input
  changed$expression[[1L]] <- 100L
  third <- grn_cached_run(changed, config, compute)
  expect_false(third$execution$cache$hit)
  expect_false(identical(third$execution$cache$key, first$execution$cache$key))
  config$method_args$alpha <- 2
  fourth <- grn_cached_run(input, config, compute)
  expect_false(identical(fourth$execution$cache$key, first$execution$cache$key))
  expect_identical(runs, 3L)
})

test_that("GRN cache rejects partial, corrupted and concurrently owned results", {
  directory <- withr::local_tempdir()
  input <- list(expression = matrix(1:4, 2), cells = c("a", "b"))
  config <- list(method = "test", work_dir = directory, resources = list())
  runs <- 0L
  compute <- function() {
    runs <<- runs + 1L
    structure(list(schema_version = "1.0", execution = list(status = "completed", artifacts = list()), groups = list()), class = "GRNResult")
  }
  first <- grn_cached_run(input, config, compute)
  cache_dir <- first$execution$cache$directory
  manifest <- readRDS(file.path(cache_dir, "completed.rds"))
  writeLines("partial", file.path(cache_dir, manifest$result_file))
  second <- grn_cached_run(input, config, compute)
  expect_false(second$execution$cache$hit)
  expect_identical(runs, 2L)
  dir.create(file.path(cache_dir, "running.lock"))
  config$force <- TRUE
  expect_error(grn_cached_run(input, config, compute), "Another GRN run owns")
  unlink(file.path(cache_dir, "running.lock"), recursive = TRUE)
  expect_error(grn_cached_run(input, config, function() stop("algorithm-sentinel")), "algorithm-sentinel")
  statuses <- lapply(list.files(cache_dir, pattern = "status[.]rds$", recursive = TRUE, full.names = TRUE), readRDS)
  expect_true(any(vapply(statuses, function(x) identical(x$status, "failed"), logical(1))))
  config$force <- FALSE
  expect_true(grn_cached_run(input, config, compute)$execution$cache$hit)
})

test_that("GRN cache validates directory artifacts before reuse", {
  directory <- withr::local_tempdir()
  output <- file.path(directory, "scientific-output")
  dir.create(output)
  input <- list(expression = matrix(1:4, 2), cells = c("a", "b"))
  config <- list(method = "test", work_dir = directory, resources = list())
  runs <- 0L
  compute <- function() {
    runs <<- runs + 1L
    writeLines(c("source\ttarget\tscore", "TF\tgene\t1"), file.path(output, "edges.tsv"))
    structure(list(schema_version = "1.0", execution = list(status = "completed", artifacts = list(list(path = output))), groups = list()), class = "GRNResult")
  }
  first <- grn_cached_run(input, config, compute)
  expect_true(grn_cached_run(input, config, compute)$execution$cache$hit)
  writeLines("corrupted", file.path(output, "edges.tsv"))
  expect_false(grn_cached_run(input, config, compute)$execution$cache$hit)
  unlink(file.path(output, "edges.tsv"))
  expect_false(grn_cached_run(input, config, compute)$execution$cache$hit)
  expect_identical(runs, 3L)
})
