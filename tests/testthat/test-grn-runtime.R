test_that("Python runner processes preserve argument boundaries and interpreter state", {
  directory <- withr::local_tempdir()
  path <- file.path(directory, "a file with spaces.txt")
  code <- "x<-commandArgs(TRUE);writeLines(x[[2]],x[[1]]);cat('stdout-marker\\n');message('stderr-marker')"
  before <- Sys.getenv(c("RETICULATE_PYTHON", "R_LIBS_USER"))
  result <- runner_run_process(file.path(R.home("bin"), "Rscript"), c("--vanilla", "-e", code, path, "literal $x; `command`"), work_dir = directory, timeout = 10)
  expect_identical(result$status, "completed")
  expect_identical(result$exit_status, 0L)
  expect_identical(readLines(path), "literal $x; `command`")
  expect_match(paste(readLines(result$stdout), collapse = ""), "stdout-marker")
  expect_match(paste(readLines(result$stderr), collapse = ""), "stderr-marker")
  expect_identical(Sys.getenv(c("RETICULATE_PYTHON", "R_LIBS_USER")), before)
  expect_identical(readRDS(result$receipt)$status, "completed")
})

test_that("Python runner failures and timeouts keep diagnostic receipts", {
  directory <- withr::local_tempdir()
  expect_error(runner_run_process(file.path(R.home("bin"), "Rscript"), c("--vanilla", "-e", "message('failure-sentinel');quit(status=7L)"), work_dir = directory, timeout = 5), "exit status 7.*Logs", fixed = FALSE)
  receipts <- lapply(list.files(directory, pattern = "[.]rds$", full.names = TRUE), readRDS)
  expect_true(any(vapply(receipts, function(x) identical(x$status, "failed") && identical(x$exit_status, 7L), logical(1))))
  expect_error(runner_run_process(file.path(R.home("bin"), "Rscript"), c("--vanilla", "-e", "Sys.sleep(30)"), work_dir = directory, timeout = 0.3), "timed out")
  receipts <- lapply(list.files(directory, pattern = "[.]rds$", full.names = TRUE), readRDS)
  expect_true(any(vapply(receipts, function(x) identical(x$status, "timed_out"), logical(1))))
})

test_that("Python runner cleanup removes descendants after parent exit or timeout", {
  directory <- withr::local_tempdir()
  parent_script <- file.path(directory, "parent.R")
  child_script <- file.path(directory, "child.R")
  writeLines(c("args <- commandArgs(TRUE)", "writeLines(as.character(Sys.getpid()), args[[1L]])", "Sys.sleep(60)"), child_script)
  writeLines(c(
    "args <- commandArgs(TRUE)",
    "child <- processx::process$new(file.path(R.home('bin'), 'Rscript'), c('--vanilla', args[[1L]], args[[2L]]), cleanup = FALSE, cleanup_tree = FALSE, stdout = args[[3L]], stderr = args[[4L]])",
    "while (!file.exists(args[[2L]])) { if(!child$is_alive()) stop('Child failed to start'); Sys.sleep(0.01) }",
    "if(args[[5L]] == 'timeout') Sys.sleep(60)",
    "quit(status = if(args[[5L]] == 'failure') 7L else 0L)"
  ), parent_script)
  for (kind in c("success", "failure", "timeout")) {
    pid_file <- file.path(directory, paste0(kind, ".pid"))
    arguments <- c(
      "--vanilla", parent_script, child_script, pid_file,
      file.path(directory, paste0(kind, "-child.stdout")), file.path(directory, paste0(kind, "-child.stderr")), kind
    )
    work_dir <- file.path(directory, kind)
    timeout <- if (kind == "timeout") 2 else 10
    if (kind == "success") {
      expect_identical(runner_run_process(file.path(R.home("bin"), "Rscript"), arguments, work_dir, timeout)$status, "completed")
    } else {
      expect_error(runner_run_process(file.path(R.home("bin"), "Rscript"), arguments, work_dir, timeout), if (kind == "failure") "exit status 7" else "timed out")
    }
    pid <- as.integer(readLines(pid_file))
    handle <- tryCatch(ps::ps_handle(pid), error = function(e) NULL)
    active <- function() {
      tryCatch(
        !is.null(handle) && ps::ps_is_running(handle) && !ps::ps_status(handle) %in% c("zombie", "dead"),
        error = function(e) FALSE
      )
    }
    deadline <- Sys.time() + 2
    while (active() && Sys.time() < deadline) Sys.sleep(0.02)
    alive <- active()
    if (alive) ps::ps_kill(handle)
    expect_false(alive)
  }
})

test_that("GRN dependency profiles remain isolated from the default environment", {
  modules <- c("grnboost2", "gniplr", "scsgl", "linger", "regdiffusion")
  expect_false(any(modules %in% default_env_modules()))
  expect_identical(env_module_environment("grnboost2"), "grnboost2_env")
  expect_identical(env_module_environment("linger"), "linger_env")
  expect_identical(env_module_environment("regdiffusion"), "regdiffusion_env")
  expect_identical(env_requirements(version = "3.9-1", modules = "grnboost2")$packages[["numpy"]], "numpy==1.24.3")
  expect_false("scanpy" %in% names(env_requirements(modules = "grnboost2")$packages))
  expect_error(env_requirements(modules = c("scenic", "regdiffusion")), "separate dependency profiles")
  expect_true("pybiomart" %in% names(env_requirements(modules = "scenicplus")$packages))
  expect_identical(env_r_packages("ppcor"), "ppcor")
})
