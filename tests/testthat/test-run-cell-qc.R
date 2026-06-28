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
