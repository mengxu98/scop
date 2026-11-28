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

test_that("GetAssayData5 complete sparse joins preserve split-layer axes and values", {
  counts_a <- Matrix::Matrix(
    matrix(c(0, 2, 3, 4), nrow = 2L,
           dimnames = list(c("g2", "g1"), c("c2", "c1"))),
    sparse = TRUE
  )
  counts_b <- Matrix::Matrix(
    matrix(c(5, 6, 7, 0), nrow = 2L,
           dimnames = list(c("g1", "g3"), c("c4", "c3"))),
    sparse = TRUE
  )
  counts <- list(A = counts_a, B = counts_b)
  assay <- SeuratObject::CreateAssay5Object(counts = counts, data = lapply(counts, log1p))
  SeuratObject::VariableFeatures(assay) <- c("g3", "g1")
  selected_features <- SeuratObject::Features(assay, layer = c("counts.A", "counts.B"))
  selected_cells <- SeuratObject::Cells(assay, layer = c("counts.A", "counts.B"))
  expect_identical(selected_features, rownames(assay))
  expect_identical(selected_cells, colnames(assay))
  initial_metadata <- GetFeaturesData(assay)

  selected <- subset(assay, features = selected_features, cells = selected_cells)
  expected_counts <- NULL
  for (layer in c("counts", "data")) {
    expected <- SeuratObject::JoinLayers(selected, layers = layer, new = layer)
    expected <- SeuratObject::GetAssayData(expected, layer = layer)
    observed <- GetAssayData5(assay, layer = layer)
    expect_identical(observed, expected)
    if (identical(layer, "counts")) expected_counts <- expected
  }

  feature_order <- c("g3", "g2", "g3")
  cell_order <- c("c3", "c2", "c1", "c4")
  expected <- expected_counts[feature_order, cell_order, drop = FALSE]
  expect_identical(
    suppressWarnings(GetAssayData5(assay, features = feature_order, cells = cell_order)),
    expected
  )
  expect_identical(SeuratObject::VariableFeatures(assay), c("g3", "g1"))
  expect_identical(GetFeaturesData(assay), initial_metadata)
  expect_true(methods::validObject(assay, test = TRUE))
})
