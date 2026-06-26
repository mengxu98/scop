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
      expect_identical(packages, "scMalignantFinder")
      invisible(TRUE)
    },
    srt_to_adata = function(srt, ...) {
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
    run_scmalignantfinder = function(test_input, pretrain_dir, ...) {
      expect_identical(test_input$cells, c("Cell1", "Cell2"))
      expect_identical(pretrain_dir, "model_dir")
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
      ...
    ) {
      expect_identical(signature_gmt, gmt)
      expect_equal(dim(spatial_coordinates), c(3, 2))
      obs
    }
  )

  with_mock_scmalignantfinder(funs, {
    out <- RunscMalignantRegion(
      srt,
      signature_gmt = gmt,
      spatial.cols = c("x", "y"),
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

test_that("RunscMalignantStates appends cancer state scores", {
  srt <- make_scmf_seurat()
  gmt <- tempfile(fileext = ".gmt")
  writeLines(
    paste("MP1 Cell Cycle - G2/M", "na", "Gene1", "Gene2", sep = intToUtf8(9)),
    gmt
  )
  obs <- data.frame(
    "MP1 Cell Cycle - G2/M" = c(0.1, 0.2, 0.3),
    "MP2 EMT" = c(0.4, 0.5, 0.6),
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
  expect_equal(unname(out[[state_cols[[1]]]][, 1]), c(0.1, 0.2, 0.3))
  expect_true("scMalignantFinder_states" %in% names(out@tools))
})

test_that("scMalignantFinder Python bridge is installed as an isolated module", {
  py_file <- system.file(
    "python",
    "scmalignantfinder.py",
    package = "scop",
    mustWork = TRUE
  )
  source <- paste(readLines(py_file, warn = FALSE), collapse = intToUtf8(10))
  expect_match(source, "def run_scmalignantfinder", fixed = TRUE)
  expect_match(source, "def run_scmalignant_region", fixed = TRUE)
  expect_match(source, "def run_scmalignant_states", fixed = TRUE)
})
