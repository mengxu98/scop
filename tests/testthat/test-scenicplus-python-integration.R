test_that("real SCENIC+ objects retain Python tables and AUC through R import", {
  python <- Sys.getenv("SCOP_TEST_PYTHON", unset = Sys.getenv("RETICULATE_PYTHON"))
  skip_if(!nzchar(python) || !file.exists(python), "Set SCOP_TEST_PYTHON to a SCENIC+ Python runtime")
  result_dir <- withr::local_tempdir()
  functions_file <- testthat::test_path("..", "..", "inst", "python", "functions.py")
  if (!file.exists(functions_file)) {
    functions_file <- system.file("python", "functions.py", package = "scop", mustWork = TRUE)
  }
  script <- testthat::test_path("..", "python", "test_scenicplus_result.py")
  log_file <- system.file("scripts", "log_message.py", package = "thisutils")
  if (!nzchar(log_file)) {
    log_file <- system.file("python", "log_message.py", package = "thisutils", mustWork = TRUE)
  }
  withr::local_envvar(PYTHONPATH = paste(dirname(log_file), Sys.getenv("PYTHONPATH"), sep = .Platform$path.sep))
  output <- suppressWarnings(system2(
    python,
    shQuote(c(script, functions_file, result_dir)),
    stdout = TRUE,
    stderr = TRUE
  ))
  expect_equal(attr(output, "status") %||% 0L, 0L, info = paste(output, collapse = "\n"))
  if (!is.null(attr(output, "status"))) {
    return(invisible(NULL))
  }
  cells <- c("c1", "c3", "c2")
  counts <- matrix(1:12, nrow = 4, dimnames = list(paste0("G", 1:4), cells))
  srt <- suppressWarnings(Seurat::CreateSeuratObject(counts))
  out <- suppressWarnings(RunSCENICPlus(
    srt,
    backend = "python",
    python_result_dir = file.path(result_dir, "result"),
    auto_prepare = FALSE,
    verbose = FALSE
  ))
  reference <- file.path(result_dir, "reference")
  result <- out@tools[["SCENICPlus"]]
  for (table in c("tf_gene", "region_gene", "tf_region", "triplets", "eregulons")) {
    expected <- utils::read.delim(file.path(reference, paste0(table, ".tsv")), check.names = FALSE)
    expect_equal(result[[table]], expected)
  }
  expected_auc <- utils::read.delim(file.path(reference, "auc.tsv"), row.names = 1, check.names = FALSE)
  expected_auc <- expected_auc[cells, , drop = FALSE]
  expect_equal(result$auc, expected_auc)
  expect_equal(unname(as.matrix(SeuratObject::GetAssayData(out[["scenicplus"]], layer = "data"))),
               unname(t(as.matrix(expected_auc))))
  expect_identical(out[["scenicplus"]][[]][["termnames"]], colnames(expected_auc))
  expect_identical(colnames(out[["scenicplus"]]), cells)
  expect_identical(unname(lengths(result$regulon_list)), c(2L, 2L))
  conda <- Sys.getenv("SCOP_TEST_CONDA")
  if (nzchar(conda)) {
    withr::local_envvar(KMP_DUPLICATE_LIB_OK = "FALSE")
    python_initialized <- reticulate::py_available(initialize = FALSE)
    native_input <- matrix(c(2, 1, 3, 1, 3, 2), nrow = 3)
    native_before <- scenicplus_rank_cpp(native_input, cores = 2L)
    extracted <- suppressWarnings(RunSCENICPlus(
      srt,
      backend = "python",
      scplus_object = file.path(result_dir, "scplus.h5mu"),
      work_dir = file.path(result_dir, "r_extracted"),
      envname = Sys.getenv("SCOP_TEST_SCENICPLUS_ENV", unset = "scenicplus_io_env"),
      conda = conda,
      verbose = FALSE
    ))
    fields <- setdiff(names(result), "parameters")
    expect_equal(extracted@tools[["SCENICPlus"]][fields], result[fields])
    expect_identical(reticulate::py_available(initialize = FALSE), python_initialized)
    expect_identical(Sys.getenv("KMP_DUPLICATE_LIB_OK"), "FALSE")
    expect_identical(scenicplus_rank_cpp(native_input, cores = 2L), native_before)
    other_python <- Sys.getenv("SCOP_TEST_OTHER_PYTHON")
    if (nzchar(other_python)) {
      child <- callr::r(function(root, python, envname, conda, object, srt) {
        if (dir.exists(file.path(root, "R"))) {
          suppressPackageStartupMessages(pkgload::load_all(root, compile = FALSE, quiet = TRUE))
        } else {
          suppressPackageStartupMessages(library(scop))
        }
        reticulate::use_python(python, required = TRUE)
        reticulate::py_run_string("import sys")
        configured <- reticulate::py_config()$python
        variables <- c("PATH", "PYTHONPATH", "RETICULATE_PYTHON", "LD_LIBRARY_PATH",
          "DYLD_LIBRARY_PATH", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "KMP_DUPLICATE_LIB_OK")
        before <- Sys.getenv(variables, unset = NA_character_)
        results <- lapply(seq_len(2L), function(i) {
          suppressWarnings(RunSCENICPlus(srt, backend = "python", scplus_object = object,
            work_dir = tempfile("scenicplus_initialized_"), envname = envname, conda = conda, verbose = FALSE))@tools$SCENICPlus$auc
        })
        list(auc = results, interpreter_unchanged = identical(reticulate::py_config()$python, configured),
          environment_unchanged = identical(Sys.getenv(variables, unset = NA_character_), before))
      }, args = list(root = getNamespaceInfo(asNamespace("scop"), "path"), python = other_python,
        envname = Sys.getenv("SCOP_TEST_SCENICPLUS_ENV", unset = "scenicplus_io_env"),
        conda = conda, object = file.path(result_dir, "scplus.h5mu"), srt = srt))
      expect_true(child$interpreter_unchanged)
      expect_true(child$environment_unchanged)
      expect_equal(child$auc, list(expected_auc, expected_auc))
    }
  }
})
