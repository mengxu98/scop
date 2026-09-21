test_that("Python logger imports and runners support both thisutils layouts", {
  root <- tempfile("thisutils-layout-")
  dir.create(root)
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  for (layout in c("scripts", "python")) {
    dir.create(file.path(root, layout))
  }
  system_file <- base::system.file
  scope <- new.env(parent = asNamespace("scop"))
  scope$system.file <- function(..., package, mustWork = FALSE) {
    if (identical(package, "thisutils")) return(file.path(root, ...))
    system_file(..., package = package, mustWork = mustWork)
  }
  scope$system2 <- function(...) Sys.getenv("PYTHONPATH")
  import <- scop_python_import
  run <- runner_system2
  environment(import) <- environment(run) <- scope
  imported <- NULL
  testthat::local_mocked_bindings(
    py_available = function(initialize = FALSE) TRUE,
    py_eval = function(code, convert = TRUE) FALSE,
    import_from_path = function(module, path, convert = TRUE) {
      if (identical(module, "log_message")) imported <<- path
      invisible(NULL)
    },
    .package = "reticulate"
  )
  for (layout in c("python", "scripts")) {
    file.create(file.path(root, layout, "log_message.py"))
    import("functions")
    expect_identical(imported, file.path(root, layout))
    for (previous in c(NA_character_, "existing-path")) {
      withr::with_envvar(c(PYTHONPATH = previous), {
        observed <- run("unused", character(), character(), TRUE, TRUE)
        expected <- paste(c(file.path(root, layout), previous[!is.na(previous)]),
          collapse = .Platform$path.sep)
        expect_identical(observed, expected)
        expect_identical(Sys.getenv("PYTHONPATH", unset = NA_character_), previous)
      })
    }
  }
  unlink(file.path(root, c("scripts", "python"), "log_message.py"))
  expect_error(import("functions"), "log_message.py")
})

test_that("runner_system2 supports an empty environment override", {
  output <- tempfile()
  error_output <- tempfile()
  status <- getFromNamespace("runner_system2", "scop")(
    command = file.path(R.home("bin"), "Rscript"),
    args = c("-e", shQuote("cat('runner-ok')")),
    env = character(),
    stdout = output,
    stderr = error_output
  )

  expect_identical(status, 0L)
  expect_identical(readLines(output, warn = FALSE), "runner-ok")
})

test_that("runner locks are exclusive and only their owner releases them", {
  lock_path <- tempfile("runner_lock_")
  acquire <- getFromNamespace("runner_acquire_lock", "scop")
  release <- getFromNamespace("runner_release_lock", "scop")
  lock <- acquire(lock_path, backend = "test backend")

  expect_error(
    acquire(lock_path, backend = "test backend"),
    "Another.*run"
  )
  expect_false(release(list(path = lock_path, token = "not-the-owner")))
  expect_true(file.exists(lock_path))
  expect_true(release(lock))
  expect_false(file.exists(lock_path))
})

test_that("runner JSON writes leave only a complete target", {
  path <- tempfile("runner_json_", fileext = ".json")
  write_json <- getFromNamespace("runner_write_json", "scop")
  read_json <- getFromNamespace("runner_read_json", "scop")

  write_json(list(value = "first"), path)
  write_json(list(value = "second"), path)

  expect_identical(read_json(path)$value, "second")
  leftovers <- list.files(
    dirname(path),
    pattern = paste0("^\\.", basename(path), "\\."),
    full.names = TRUE
  )
  expect_length(leftovers, 0L)
})

test_that("runner JSON writes report failed atomic replacements", {
  target <- tempfile("runner-json-directory-")
  dir.create(target)

  expect_error(
    getFromNamespace("runner_write_json", "scop")(
      list(value = "new"),
      target
    ),
    "Unable to atomically write"
  )
  expect_true(dir.exists(target))
})

test_that("runner log tails stay bounded and retain the final lines", {
  output <- tempfile()
  writeLines(
    c(paste("line", seq_len(5000L)), "", "final sentinel"),
    output
  )

  tail_lines <- getFromNamespace("runner_tail_lines", "scop")(
    output,
    max_lines = 20L,
    chunk_size = 37L
  )

  expect_length(tail_lines, 20L)
  expect_identical(tail_lines[[20L]], "final sentinel")
  expect_false(any(!nzchar(tail_lines)))
})

test_that("Python distribution matching follows canonical package names", {
  canonical <- getFromNamespace("canonical_python_distribution_name", "scop")
  expect_identical(
    canonical(c("scCODA", "tf_keras", "scvi.tools", "python-igraph")),
    c("sccoda", "tf-keras", "scvi-tools", "python-igraph")
  )
})

test_that("Python requirements preserve names while removing duplicates", {
  unique_requirements <- getFromNamespace("unique_python_requirements", "scop")
  packages <- c(
    jax = "jax[cpu]==0.4.38",
    scanpy = "scanpy==1.11.3",
    jax_duplicate = "jax[cpu]==0.4.38"
  )

  expect_identical(
    unique_requirements(packages),
    c(jax = "jax[cpu]==0.4.38", scanpy = "scanpy==1.11.3")
  )
})

test_that("Python requirement parsing treats extras as distribution metadata", {
  parse <- getFromNamespace("parse_python_requirement", "scop")

  expect_identical(
    parse("jax[cpu]==0.4.38"),
    list(name = "jax", operator = "==", version = "0.4.38")
  )
  expect_identical(
    parse("requests[socks]"),
    list(name = "requests", operator = NA_character_, version = NA_character_)
  )
})
