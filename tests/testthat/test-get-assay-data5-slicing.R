test_that("GetAssayData5 slices features and cells on Assay5 objects", {
  data("pancreas_sub", package = "scop")
  feats <- head(rownames(pancreas_sub), 3)
  cells <- colnames(pancreas_sub)[1:20]

  full <- GetAssayData5(pancreas_sub, layer = "counts")
  sliced <- GetAssayData5(
    pancreas_sub,
    layer = "counts",
    features = c(feats[2], feats[1], "NOT_A_GENE"),
    cells = cells
  )

  expect_identical(rownames(sliced), feats[c(2, 1)])
  expect_identical(colnames(sliced), cells)
  expect_equal(
    as.matrix(sliced),
    as.matrix(full[feats[c(2, 1)], cells, drop = FALSE])
  )
})

test_that("GetAssayData5 slices features and cells on Assay objects", {
  data("pancreas_sub", package = "scop")
  counts <- GetAssayData5(pancreas_sub, layer = "counts")
  assay_v3 <- SeuratObject::CreateAssayObject(counts = counts)
  feats <- head(rownames(assay_v3), 3)
  cells <- colnames(assay_v3)[11:30]

  sliced <- GetAssayData5(
    assay_v3,
    layer = "counts",
    features = c(feats[3], "NOT_A_GENE", feats[1]),
    cells = cells
  )

  expect_identical(rownames(sliced), feats[c(3, 1)])
  expect_identical(colnames(sliced), cells)
  expect_equal(
    as.matrix(sliced),
    as.matrix(counts[feats[c(3, 1)], cells, drop = FALSE])
  )
})

test_that("GetAssayData5.Seurat passes slicing through to the assay", {
  data("pancreas_sub", package = "scop")
  feats <- head(rownames(pancreas_sub), 2)
  cells <- colnames(pancreas_sub)[1:10]

  sliced <- GetAssayData5(
    pancreas_sub,
    layer = "counts",
    assay = "RNA",
    features = feats,
    cells = cells
  )

  expect_identical(dimnames(sliced), list(feats, cells))
})
