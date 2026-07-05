test_that("db_scds cxds path does not run hybrid backend", {
  skip_if_not_installed("scds")
  skip_if_not_installed("SingleCellExperiment")

  counts <- Matrix::Matrix(
    c(
      1, 0, 2, 0,
      0, 3, 0, 4,
      5, 0, 0, 6
    ),
    nrow = 3,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)

  testthat::local_mocked_bindings(
    cxds = function(sce, ...) {
      sce[["cxds_score"]] <- seq_len(ncol(sce))
      sce
    },
    cxds_bcds_hybrid = function(...) {
      stop("hybrid backend should not be called")
    },
    .package = "scds"
  )

  out <- suppressWarnings(
    db_scds(
      srt,
      method = "cxds",
      db_rate = 0.5,
      data_type = "raw_counts",
      verbose = FALSE
    )
  )

  expect_true("db.scds_cxds_score" %in% colnames(out[[]]))
  expect_true("db.scds_cxds_class" %in% colnames(out[[]]))
  expect_false("db.scds_hybrid_score" %in% colnames(out[[]]))
})

test_that("db_Scrublet passes integer PCA component arguments", {
  counts <- Matrix::Matrix(
    c(
      1, 0, 2, 0,
      0, 3, 0, 4,
      5, 0, 0, 6
    ),
    nrow = 3,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)

  observed <- new.env(parent = emptyenv())
  scrub <- new.env(parent = emptyenv())
  scrub$threshold_ <- 0.5
  scrub$scrub_doublets <- function(...) {
    observed$args <- list(...)
    list(rep(0.1, ncol(srt)), rep(FALSE, ncol(srt)))
  }
  scrublet <- list(
    Scrublet = function(...) {
      scrub
    }
  )

  testthat::local_mocked_bindings(
    PrepareEnv = function(...) invisible(NULL),
    check_python = function(...) invisible(TRUE),
    CheckDataType = function(...) "raw_counts",
    .package = "scop"
  )
  testthat::local_mocked_bindings(
    import = function(name, ...) {
      expect_identical(name, "scrublet")
      scrublet
    },
    py_has_attr = function(x, name) TRUE,
    py_to_r = function(x) x,
    .package = "reticulate"
  )

  out <- db_Scrublet(
    srt,
    db_rate = 0.01,
    n_prin_comps = 10,
    min_counts = 1,
    min_cells = 1,
    data_type = "raw_counts",
    verbose = FALSE
  )

  expect_type(observed$args$n_prin_comps, "integer")
  expect_type(observed$args$min_counts, "integer")
  expect_type(observed$args$min_cells, "integer")
  expect_true("db.Scrublet_score" %in% colnames(out[[]]))
})
