test_that("prep_ref aligns an Assay5 subset before metadata filtering", {
  counts <- matrix(seq_len(24), 4, dimnames = list(paste0("g", 1:4), paste0("c", 1:6)))
  srt <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE))
  srt$label <- c("MG", "OL", NA, "MG", "OL", "MG")
  srt$sample <- c("a", "a", "b", "b", "c", "c")
  srt$state <- paste0("s", 1:6)
  srt[["subset"]] <- SeuratObject::CreateAssay5Object(counts = Matrix::Matrix(counts[, c(6, 3, 2), drop = FALSE], sparse = TRUE))
  out <- prep_ref(srt, group.by = "label", sample.by = "sample", cellstate.by = "state", assay = "subset", dense = FALSE)
  expect_identical(colnames(out$counts), c("c2", "c6"))
  expect_identical(rownames(out$meta), colnames(out$counts))
  expect_identical(out$meta$celltype_scop, c("OL", "MG"))
  expect_identical(out$meta$sample_scop, c("a", "c"))
  expect_identical(out$meta$cellstate_scop, c("s2", "s6"))
  expect_equal(as.matrix(out$counts), counts[, c(2, 6), drop = FALSE])
})

test_that("prep_ref retains full-reference behavior", {
  counts <- matrix(1:12, 3, dimnames = list(paste0("g",1:3), paste0("c",1:4)))
  srt <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE))
  srt$label <- c("MG", NA, "OL", "MG"); srt$sample <- "a"
  out <- prep_ref(srt, group.by = "label", sample.by = "sample")
  expect_equal(out$counts, counts[, c(1,3,4), drop = FALSE])
  expect_identical(rownames(out$meta), colnames(out$counts))
})


test_that("prep_ref filters shuffled counts using barcode metadata", {
  counts <- matrix(1:12, 3, dimnames = list(paste0("g",1:3), paste0("c",1:4)))
  srt <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE))
  srt$label <- c("MG", NA, "OL", "MG"); srt$sample <- c("a", "b", "c", "d")
  testthat::local_mocked_bindings(GetAssayData5 = function(...) counts[, c(4,2,1), drop = FALSE])
  out <- prep_ref(srt, group.by = "label", sample.by = "sample")
  expect_identical(colnames(out$counts), c("c4", "c1"))
  expect_identical(out$meta$sample_scop, c("d", "a"))
  expect_identical(rownames(out$meta), colnames(out$counts))
})
