test_that("pancreas_sub nCount_RNA matches counts colSums used by the fast path", {
  data("pancreas_sub", package = "scop")
  expect_true("nCount_RNA" %in% colnames(pancreas_sub@meta.data))
  expect_equal(
    pancreas_sub@meta.data$nCount_RNA,
    unname(Matrix::colSums(GetAssayData5(pancreas_sub, layer = "counts")))
  )
})

test_that("GroupHeatmap libsize uses meta.data fast path and colSums fallback", {
  data("pancreas_sub", package = "scop")
  feats <- head(rownames(pancreas_sub), 5)

  p_fast <- GroupHeatmap(
    pancreas_sub,
    features = feats,
    group.by = "CellType",
    lib_normalize = TRUE
  )
  expect_false(is.null(p_fast))

  pancreas_sub@meta.data$nCount_RNA <- NULL
  p_fallback <- GroupHeatmap(
    pancreas_sub,
    features = feats,
    group.by = "CellType",
    lib_normalize = TRUE
  )
  expect_false(is.null(p_fallback))
})
