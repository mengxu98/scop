test_that("disk-backed preprocessing matches official Seurat without replacing counts", {
  skip_if_not_installed("BPCells")
  set.seed(221)
  counts <- Matrix::Matrix(matrix(rpois(120 * 60, 4), 120,
    dimnames = list(paste0("g", 1:120), paste0("c", 1:60))
  ), sparse = TRUE)
  directory <- tempfile("spatial-counts-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  write_matrix <- get_namespace_fun("BPCells", "write_matrix_dir")
  open_matrix <- get_namespace_fun("BPCells", "open_matrix_dir")
  write_matrix(counts, dir = directory)
  object <- SeuratObject::CreateSeuratObject(open_matrix(directory))
  original <- object
  expect_identical(CheckDataType(GetAssayData5(object), verbose = FALSE), "unknown")
  expected <- Seurat::NormalizeData(object, verbose = FALSE)
  actual <- NormalizeData(object, verbose = FALSE)
  expect_s4_class(GetAssayData5(actual, layer = "counts"), "IterableMatrix")
  expect_s4_class(GetAssayData5(actual, layer = "data"), "IterableMatrix")
  expect_equal(as.matrix(GetAssayData5(actual, layer = "data")),
    as.matrix(GetAssayData5(expected, layer = "data")),
    tolerance = 1e-10
  )
  actual <- FindVariableFeatures(actual, nfeatures = 25, verbose = FALSE)
  expected <- Seurat::FindVariableFeatures(expected, nfeatures = 25, verbose = FALSE)
  expect_identical(SeuratObject::VariableFeatures(actual), SeuratObject::VariableFeatures(expected))
  actual <- ScaleData(actual, verbose = FALSE)
  expected <- Seurat::ScaleData(expected, verbose = FALSE)
  expect_equal(as.matrix(GetAssayData5(actual, layer = "scale.data")),
    as.matrix(GetAssayData5(expected, layer = "scale.data")),
    tolerance = 1e-8
  )
  expect_identical(GetAssayData5(actual), GetAssayData5(original))
  kept <- CheckDataList(list(actual),
    batch = "", assay = "RNA", do_normalization = FALSE,
    do_HVF_finding = FALSE, HVF = SeuratObject::VariableFeatures(actual), nHVF = 25, verbose = FALSE
  )
  expect_identical(GetAssayData5(kept$srt_list[[1L]], layer = "data"), GetAssayData5(actual, layer = "data"))
  sliced <- GetAssayData5(actual, features = c("g8", "g2"), cells = c("c5", "c1"))
  expect_identical(dimnames(sliced), list(c("g8", "g2"), c("c5", "c1")))
  expect_s4_class(sliced, "IterableMatrix")
  expect_equal(as.matrix(sliced), as.matrix(counts[c("g8", "g2"), c("c5", "c1")]))
  actual$x <- rep(1:10, 6)
  actual$y <- rep(1:6, each = 10)
  result <- RunSpatialVariableFeatures(actual,
    features = rownames(actual)[1:5],
    k = 3, min_spots = 2, nperm = 0, backend = "cpp", verbose = FALSE
  )
  expect_s4_class(GetAssayData5(result), "IterableMatrix")
  expect_setequal(result@tools$SpatialVariableFeatures$result$feature, rownames(actual)[1:5])
})

test_that("split-layer extraction selects before joining and preserves requested order", {
  counts <- Matrix::Matrix(matrix(1:24, 4, dimnames = list(paste0("g", 1:4), paste0("c", 1:6))), sparse = TRUE)
  object <- SeuratObject::CreateAssay5Object(counts = list(a = counts[, 1:3], b = counts[, 4:6]))
  before <- object
  out <- GetAssayData5(object, features = c("g4", "g1"), cells = c("c6", "c2"))
  expect_equal(as.matrix(out), as.matrix(counts[c("g4", "g1"), c("c6", "c2")]))
  expect_identical(object, before)
  expect_equal(dim(GetAssayData5(object, features = "absent")), c(0L, 6L))
})
