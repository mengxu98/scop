pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)

make_scmf_seurat <- function() {
  counts <- matrix(
    c(
      4, 0, 1,
      0, 3, 2,
      1, 2, 5
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:3), paste0("Cell", 1:3))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$x <- c(1, 2, 3)
  srt$y <- c(4, 5, 6)
  srt
}

with_mock_scmalignantfinder <- function(funs, code) {
  testthat::local_mocked_bindings(
    PrepareEnv = function(modules, ...) {
      expect_identical(modules, "scanpy")
      invisible(TRUE)
    },
    check_python = function(packages, ...) {
      expect_true(packages %in% c("scMalignantFinder", "xgboost"))
      invisible(TRUE)
    },
    scmf_python_classifier_available = function() {
      TRUE
    },
    srt_to_adata = function(srt, layer_x, ...) {
      expect_identical(layer_x, "counts")
      list(cells = colnames(srt))
    },
    import_scmalignantfinder = function(convert = TRUE) {
      funs
    }
  )
  force(code)
}

test_that("RunscMalignantFinder appends predictions by cell barcode", {
  srt <- make_scmf_seurat()
  obs <- data.frame(
    scMalignantFinder_prediction = c("Malignant", "Normal"),
    malignancy_probability = c(0.91, 0.12),
    row.names = c("Cell1", "Cell2")
  )
  funs <- list(
    run_scmalignantfinder = function(test_input, pretrain_dir, norm_type, ...) {
      expect_identical(test_input$cells, c("Cell1", "Cell2"))
      expect_identical(pretrain_dir, "model_dir")
      expect_true(norm_type)
      obs
    }
  )

  with_mock_scmalignantfinder(funs, {
    out <- RunscMalignantFinder(
      srt,
      cells = c("Cell1", "Cell2"),
      pretrain_dir = "model_dir",
      verbose = FALSE
    )
  })

  expect_equal(
    unname(out$scMalignantFinder_prediction),
    c("Malignant", "Normal", NA)
  )
  expect_equal(unname(out$malignancy_probability[1:2]), c(0.91, 0.12))
  expect_true(is.na(out$malignancy_probability[[3]]))
  expect_true("scMalignantFinder" %in% names(out@tools))
})

test_that("RunscMalignantFinder respects explicit norm_type and expands h5ad paths", {
  obs <- data.frame(
    scMalignantFinder_prediction = "Normal",
    malignancy_probability = 0.1,
    row.names = "Cell1"
  )
  funs <- list(
    run_scmalignantfinder = function(test_input, norm_type, ...) {
      expect_identical(test_input, path.expand("~/input.h5ad"))
      expect_false(norm_type)
      obs
    }
  )

  with_mock_scmalignantfinder(funs, {
    out <- RunscMalignantFinder(
      h5ad = "~/input.h5ad",
      pretrain_dir = "model_dir",
      norm_type = FALSE,
      return_seurat = FALSE,
      verbose = FALSE
    )
  })

  expect_equal(out$malignancy_probability, 0.1)
})

test_that("RunscMalignantFinder checks xgboost for scratch XGBoost training", {
  checks <- character()
  testthat::local_mocked_bindings(
    PrepareEnv = function(...) invisible(TRUE),
    scmf_python_classifier_available = function() {
      TRUE
    },
    check_python = function(packages, ...) {
      checks <<- c(checks, packages)
      invisible(TRUE)
    },
    import_scmalignantfinder = function(convert = TRUE) {
      list(
        run_scmalignantfinder = function(...) {
          data.frame(
            scMalignantFinder_prediction = "Normal",
            malignancy_probability = 0.1,
            row.names = "Cell1"
          )
        }
      )
    }
  )

  RunscMalignantFinder(
    h5ad = "input.h5ad",
    train_h5ad_path = "train.h5ad",
    feature_path = "features.tsv",
    model_method = "XGBoost",
    return_seurat = FALSE,
    verbose = FALSE
  )

  expect_true("xgboost" %in% checks)
})

test_that("RunscMalignantFinder respects an already usable explicit Python", {
  calls <- character()
  obs <- data.frame(
    scMalignantFinder_prediction = "Normal",
    malignancy_probability = 0.1,
    row.names = "Cell1"
  )
  withr::local_envvar(RETICULATE_PYTHON = "python")
  testthat::local_mocked_bindings(
    scmf_python_classifier_available = function() {
      TRUE
    },
    PrepareEnv = function(...) {
      calls <<- c(calls, "PrepareEnv")
      invisible(TRUE)
    },
    check_python = function(...) {
      calls <<- c(calls, "check_python")
      invisible(TRUE)
    },
    import_scmalignantfinder = function(convert = TRUE) {
      list(
        run_scmalignantfinder = function(...) obs
      )
    }
  )

  out <- RunscMalignantFinder(
    h5ad = "input.h5ad",
    pretrain_dir = "model_dir",
    return_seurat = FALSE,
    verbose = FALSE
  )

  expect_equal(out$malignancy_probability, 0.1)
  expect_length(calls, 0)
})

test_that("RunscMalignantFinder falls back to GitHub install when package check is insufficient", {
  checks <- character()
  installed <- FALSE
  obs <- data.frame(
    scMalignantFinder_prediction = "Normal",
    malignancy_probability = 0.1,
    row.names = "Cell1"
  )
  testthat::local_mocked_bindings(
    PrepareEnv = function(...) invisible(TRUE),
    check_python = function(packages, ...) {
      checks <<- c(checks, packages)
      invisible(FALSE)
    },
    scmf_python_classifier_available = function() {
      installed
    },
    scmf_install_python_github = function(...) {
      installed <<- TRUE
      invisible(TRUE)
    },
    import_scmalignantfinder = function(convert = TRUE) {
      list(run_scmalignantfinder = function(...) obs)
    }
  )

  out <- RunscMalignantFinder(
    h5ad = "input.h5ad",
    pretrain_dir = "model_dir",
    return_seurat = FALSE,
    verbose = FALSE
  )

  expect_true(installed)
  expect_true("scMalignantFinder" %in% checks)
  expect_equal(out$malignancy_probability, 0.1)
})

test_that("RunscMalignantRegion appends spatial region outputs", {
  srt <- make_scmf_seurat()
  gmt <- tempfile(fileext = ".gmt")
  writeLines(paste("Malignant_up", "na", "Gene1", "Gene2", sep = intToUtf8(9)), gmt)
  obs <- data.frame(
    cluster = c(1, 2, 1),
    region_prediction = c("Malignant", "Normal", "Malignant"),
    Malignant_up = c(0.8, 0.2, 0.7),
    row.names = paste0("Cell", 1:3)
  )
  funs <- list(
    run_scmalignant_region = function(
      test_input,
      signature_gmt,
      spatial_coordinates,
      spatial_key,
      ...
    ) {
      expect_identical(signature_gmt, gmt)
      expect_equal(dim(spatial_coordinates), c(3, 2))
      expect_identical(spatial_key, "custom_spatial")
      obs
    }
  )

  with_mock_scmalignantfinder(funs, {
    out <- RunscMalignantRegion(
      srt,
      signature_gmt = gmt,
      spatial.cols = c("x", "y"),
      spatial_key = "custom_spatial",
      verbose = FALSE
    )
  })

  expect_equal(
    unname(out$scMalignantFinder_region_prediction),
    c("Malignant", "Normal", "Malignant")
  )
  expect_equal(unname(out$scMalignantFinder_Malignant_up), c(0.8, 0.2, 0.7))
  expect_true("scMalignantFinder_region" %in% names(out@tools))
})

test_that("RunscMalignantRegion rejects invalid spatial metadata", {
  srt <- make_scmf_seurat()
  gmt <- tempfile(fileext = ".gmt")
  writeLines(paste("Malignant_up", "na", "Gene1", "Gene2", sep = intToUtf8(9)), gmt)
  srt$bad <- c(1, Inf, 3)

  expect_error(
    RunscMalignantRegion(
      srt,
      signature_gmt = gmt,
      spatial.cols = c("x", "bad"),
      verbose = FALSE
    ),
    "finite numeric"
  )
  expect_error(
    RunscMalignantRegion(
      srt,
      signature_gmt = gmt,
      spatial.cols = "x",
      verbose = FALSE
    ),
    "exactly two"
  )
})

test_that("RunscMalignantStates appends cancer state scores", {
  srt <- make_scmf_seurat()
  gmt <- tempfile(fileext = ".gmt")
  writeLines(
    paste("MP1 Cell Cycle - G2/M", "na", "Gene1", "Gene2", sep = intToUtf8(9)),
    gmt
  )
  obs <- data.frame(
    "MP1-A" = c(0.1, 0.2, 0.3),
    "MP1 A" = c(0.4, 0.5, 0.6),
    row.names = paste0("Cell", 1:3),
    check.names = FALSE
  )
  funs <- list(
    run_scmalignant_states = function(test_input, gene_sets, ...) {
      expect_identical(gene_sets, gmt)
      obs
    }
  )

  with_mock_scmalignantfinder(funs, {
    out <- RunscMalignantStates(
      srt,
      gene_sets = gmt,
      verbose = FALSE
    )
  })

  state_cols <- grep("^scMalignantState_MP", colnames(out[[]]), value = TRUE)
  expect_length(state_cols, 2)
  expect_identical(anyDuplicated(state_cols), 0L)
  expect_equal(unname(out[[state_cols[[1]]]][, 1]), c(0.1, 0.2, 0.3))
  expect_true("scMalignantFinder_states" %in% names(out@tools))
})

test_that("scMalignantFinder Python bridge is installed as an isolated module", {
  py_file <- c(
    system.file("python", "scmalignantfinder.py", package = "scop", mustWork = FALSE),
    file.path("inst", "python", "scmalignantfinder.py"),
    file.path(testthat::test_path("..", ".."), "inst", "python", "scmalignantfinder.py")
  )
  py_file <- py_file[nzchar(py_file) & file.exists(py_file)][1]
  skip_if(is.na(py_file), "scMalignantFinder Python bridge file is not available")
  source <- paste(readLines(py_file, warn = FALSE), collapse = intToUtf8(10))
  expect_match(source, "def run_scmalignantfinder", fixed = TRUE)
  expect_match(source, "def run_scmalignant_region", fixed = TRUE)
  expect_match(source, "def run_scmalignant_states", fixed = TRUE)

  python <- Sys.which("python3")
  skip_if(!nzchar(python), "python3 is not available")
  out <- system2(
    python,
    c("-m", "py_compile", py_file),
    stdout = TRUE,
    stderr = TRUE
  )
  expect_equal(attr(out, "status") %||% 0L, 0L, info = paste(out, collapse = "\n"))
})
