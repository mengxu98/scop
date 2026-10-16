make_recovery_object <- function(values, storage = c("sparse", "dense", "legacy")) {
  storage <- match.arg(storage)
  sparse <- Matrix::Matrix(values, sparse = TRUE)
  assay <- if (storage == "legacy") {
    SeuratObject::CreateAssayObject(counts = sparse)
  } else {
    SeuratObject::CreateAssay5Object(counts = sparse)
  }
  object <- SeuratObject::CreateSeuratObject(assay)
  if (storage == "dense") {
    object <- SeuratObject::SetAssayData(object, layer = "counts", new.data = values)
  }
  object
}

test_that("RecoverCounts exp recovers sparse and dense assay inputs", {
  counts <- matrix(c(1, 2, 3, 1, 3, 5, 1, 4, 2), nrow = 3,
    dimnames = list(paste0("g", 1:3), paste0("c", 1:3)))
  normalized <- sweep(counts, 2, colSums(counts), "/") * 10000
  for (storage in c("sparse", "dense", "legacy")) {
    for (sf in list(10000, NULL)) {
      object <- make_recovery_object(log(normalized), storage)
      if (storage == "dense") {
        expect_true(is.matrix(GetAssayData5(object, layer = "counts")))
      }
      expect_identical(CheckDataType(GetAssayData5(object, layer = "counts"), verbose = FALSE), "log_normalized_counts")
      actual <- RecoverCounts(object, trans = "exp", sf = sf, min_count = 1, verbose = FALSE)
      recovered <- GetAssayData5(actual, layer = "counts")
      expect_s4_class(recovered, "dgCMatrix")
      expect_equal(as.matrix(recovered), counts)
      expect_identical(names(actual$nCount_RNA), colnames(counts))
      expect_equal(unname(actual$nCount_RNA), unname(colSums(counts)))
    }
  }
})

test_that("RecoverCounts exp expands implicit zeros to positive counts", {
  counts <- matrix(c(1, 2, 7, 1, 3, 6, 1, 4, 5), nrow = 3,
    dimnames = list(paste0("g", 1:3), paste0("c", 1:3)))
  object <- make_recovery_object(log(counts))
  input <- GetAssayData5(object, layer = "counts")
  expect_lt(length(input@x), length(counts))
  actual <- RecoverCounts(object, trans = "exp", sf = 10, min_count = 1, verbose = FALSE)
  expect_equal(as.matrix(GetAssayData5(actual, layer = "counts")), counts)
  expect_equal(unname(actual$nCount_RNA), rep(10, 3))
  expect_equal(unname(actual$nFeature_RNA), rep(3, 3))
})

test_that("RecoverCounts preserves single-cell dimensions and other assays", {
  counts <- matrix(c(1, 2, 3), ncol = 1, dimnames = list(paste0("g", 1:3), "only_cell"))
  normalized <- counts / sum(counts) * 10000
  object <- make_recovery_object(counts)
  object[["logged"]] <- SeuratObject::CreateAssay5Object(counts = Matrix::Matrix(log(normalized), sparse = TRUE))
  original_rna <- object[["RNA"]]
  original_idents <- SeuratObject::Idents(object)
  actual <- RecoverCounts(object, assay = "logged", trans = "exp", min_count = 1, verbose = FALSE)
  expect_equal(as.matrix(GetAssayData5(actual, assay = "logged", layer = "counts")), counts)
  expect_identical(actual[["RNA"]], original_rna)
  expect_identical(SeuratObject::DefaultAssay(actual), "RNA")
  expect_identical(SeuratObject::Idents(actual), original_idents)
  expect_equal(unname(actual$nCount_logged), 6)
})

test_that("RecoverCounts preserves expm1, none and raw-count behavior", {
  counts <- matrix(c(1, 0, 3, 1, 3, 0, 1, 4, 2), nrow = 3,
    dimnames = list(paste0("g", 1:3), paste0("c", 1:3)))
  normalized <- sweep(counts, 2, colSums(counts), "/") * 10000
  for (trans in c("expm1", "none")) {
    values <- if (trans == "expm1") log1p(normalized) else normalized
    object <- make_recovery_object(values)
    actual <- RecoverCounts(object, trans = trans, sf = 10000, min_count = 1, verbose = FALSE)
    expect_equal(as.matrix(GetAssayData5(actual, layer = "counts")), counts)
    expect_equal(unname(actual$nCount_RNA), unname(colSums(counts)))
    expect_identical(actual$nFeature_RNA, object$nFeature_RNA)
  }
  raw <- make_recovery_object(counts)
  expect_identical(RecoverCounts(raw, trans = "exp", verbose = FALSE), raw)
})

test_that("RecoverCounts keeps validation and unsuccessful recovery non-mutating", {
  counts <- matrix(c(1, 2, 3, 1, 3, 5, 1, 4, 2), nrow = 3,
    dimnames = list(paste0("g", 1:3), paste0("c", 1:3)))
  normalized <- sweep(counts, 2, colSums(counts), "/") * 10000
  object <- make_recovery_object(log(normalized))
  expect_error(RecoverCounts(object, trans = "invalid", verbose = FALSE), "arg")
  expect_warning(actual <- RecoverCounts(object, trans = "exp", sf = c(10000, 20000), verbose = TRUE), "Scale factor")
  expect_identical(actual, object)
  expect_warning(actual <- RecoverCounts(object, trans = "exp", sf = 10000, min_count = 0.2, tolerance = 0.01, verbose = TRUE), "presumed nCount")
  expect_identical(actual, object)
})
