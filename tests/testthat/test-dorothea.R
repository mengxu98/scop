test_that("RunDorothea writes TF activity metadata when storing a new assay", {
  counts <- matrix(
    c(
      5, 0, 4,
      0, 3, 2
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Gene1", "Gene2"), paste0("Cell", 1:3))
  )
  srt <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  )
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  regulons <- data.frame(
    tf = c("TF1", "TF2"),
    target = c("Gene1", "Gene2"),
    mor = c(1, -1),
    confidence = c("A", "B")
  )

  testthat::local_mocked_bindings(
    check_r = function(...) TRUE,
    dorothea_get_run_fun = function(method) {
      function(...) {
        data.frame(
          source = rep(c("TF1", "TF2"), each = 3),
          condition = rep(colnames(srt), times = 2),
          score = c(1, 2, 3, -1, -2, -3)
        )
      }
    }
  )

  out <- RunDorothea(
    srt,
    regulons = regulons,
    confidence = c("A", "B"),
    method = "ulm",
    new_assay = TRUE,
    verbose = FALSE
  )

  expect_true("dorothea" %in% SeuratObject::Assays(out))
  expect_true(all(c("dorothea_TF1", "dorothea_TF2") %in% colnames(out@meta.data)))
  expect_equal(unname(out$dorothea_TF1), c(1, 2, 3))
  expect_equal(unname(out$dorothea_TF2), c(-1, -2, -3))
})

test_that("RunDorothea can skip TF activity metadata", {
  counts <- matrix(
    c(
      5, 0, 4,
      0, 3, 2
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Gene1", "Gene2"), paste0("Cell", 1:3))
  )
  srt <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  )
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  regulons <- data.frame(
    tf = c("TF1", "TF2"),
    target = c("Gene1", "Gene2"),
    mor = c(1, -1),
    confidence = c("A", "B")
  )

  testthat::local_mocked_bindings(
    check_r = function(...) TRUE,
    dorothea_get_run_fun = function(method) {
      function(...) {
        data.frame(
          source = rep(c("TF1", "TF2"), each = 3),
          condition = rep(colnames(srt), times = 2),
          score = c(1, 2, 3, -1, -2, -3)
        )
      }
    }
  )

  out <- RunDorothea(
    srt,
    regulons = regulons,
    confidence = c("A", "B"),
    method = "ulm",
    new_assay = TRUE,
    add_meta = FALSE,
    verbose = FALSE
  )

  expect_true("dorothea" %in% SeuratObject::Assays(out))
  expect_false(any(grepl("^dorothea_", colnames(out@meta.data))))
  expect_true("scores" %in% names(out@tools$Dorothea))
  expect_false(out@tools$Dorothea$parameters$add_meta)
})
