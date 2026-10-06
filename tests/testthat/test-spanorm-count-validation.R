test_that("SpaNorm rejects invalid counts without silently changing them", {
  srt <- SeuratObject::CreateSeuratObject(counts = matrix(1:8, 2,
    dimnames = list(c("g1", "g2"), paste0("s", 1:4))))
  srt$x <- 1:4
  srt$y <- c(0, 0, 1, 1)
  for (bad in c(-1, NA_real_, NaN, Inf, -Inf, 1.5)) {
    mat <- Matrix::Matrix(matrix(c(bad, 2:8), 2,
      dimnames = list(c("g1", "g2"), paste0("s", 1:4))), sparse = TRUE)
    testthat::local_mocked_bindings(GetAssayData5 = function(...) mat)
    expect_error(spanorm_prepare_input(srt, "RNA", "counts",
      coord.cols = c("x", "y"), coordinate_space = "raw"),
      "counts must contain only finite, non-negative integer values")
  }
})

test_that("SpaNorm retains valid sparse counts and aligns spot identifiers", {
  cells <- paste0("s", 1:4)
  mat <- Matrix::Matrix(matrix(c(0, 2:8), 2,
    dimnames = list(c("g2", "g1"), cells)), sparse = TRUE)
  srt <- SeuratObject::CreateSeuratObject(counts = mat)
  srt$x <- 1:4
  srt$y <- c(0, 0, 1, 1)
  testthat::local_mocked_bindings(GetAssayData5 = function(...) mat[, rev(cells)])
  out <- spanorm_prepare_input(srt, "RNA", "counts", coord.cols = c("x", "y"),
    coordinate_space = "raw")
  expect_equal(as.matrix(out$counts), as.matrix(mat))
  expect_identical(colnames(out$counts), rownames(out$coords))
})
