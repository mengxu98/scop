make_atac_hvf_rank_object <- function(offset = 0L, prefix = "cell") {
  skip_if_not_installed("Signac")
  counts <- Matrix::Matrix(matrix(rep(12:1, 8), nrow = 12), sparse = TRUE)
  dimnames(counts) <- list(paste0("peak", 1:12), paste0(prefix, 1:8))
  ranges <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(seq(1, by = 1000, length.out = 12) + offset, width = 50)
  )
  assay <- Signac::CreateChromatinAssay(
    counts = counts, ranges = ranges, min.cells = 0, min.features = 0,
    verbose = FALSE
  )
  SeuratObject::CreateSeuratObject(assay, assay = "peaks")
}

check_atac_hvf_rank_objects <- function(objects, n = 3, minimum = 1, find = TRUE) {
  CheckDataList(
    objects,
    batch = paste0("batch", seq_along(objects)), assay = "peaks",
    normalization_method = "TFIDF", do_normalization = FALSE,
    do_HVF_finding = find, HVF_source = "separate", nHVF = n,
    HVF_min_intersection = minimum, verbose = FALSE
  )
}

test_that("single chromatin HVF truncation preserves native peak rank", {
  object <- make_atac_hvf_rank_object()
  before <- SeuratObject::LayerData(object[["peaks"]], layer = "counts")
  result <- check_atac_hvf_rank_objects(list(object))
  native <- Seurat::SelectIntegrationFeatures(
    result$srt_list,
    nfeatures = 3, assay = "peaks", verbose = FALSE
  )
  expect_identical(result$HVF, native)
  expect_identical(match(result$HVF, rownames(object)), 1:3)
  expect_identical(SeuratObject::LayerData(object[["peaks"]], layer = "counts"), before)
  expect_identical(colnames(result$srt_list[[1]]), colnames(object))
})

test_that("uniform toy coordinate translation does not change selected peak rank", {
  original <- make_atac_hvf_rank_object()
  translated <- make_atac_hvf_rank_object(offset = 10000000L)
  left <- check_atac_hvf_rank_objects(list(original))
  right <- check_atac_hvf_rank_objects(list(translated))
  expect_identical(match(left$HVF, rownames(original)), match(right$HVF, rownames(translated)))
  expect_identical(match(left$HVF, rownames(original)), 1:3)
  expect_identical(
    unname(SeuratObject::LayerData(original[["peaks"]], layer = "counts")),
    unname(SeuratObject::LayerData(translated[["peaks"]], layer = "counts"))
  )
})

test_that("multiple batches use native occurrence and median-rank ties", {
  first <- make_atac_hvf_rank_object(prefix = "first")
  second <- make_atac_hvf_rank_object(prefix = "second")
  ids <- rownames(first)
  SeuratObject::VariableFeatures(first) <- ids[c(1, 2, 3, 4)]
  SeuratObject::VariableFeatures(second) <- ids[c(2, 3, 5, 1)]
  result <- check_atac_hvf_rank_objects(list(first, second), find = FALSE)
  native <- Seurat::SelectIntegrationFeatures(
    result$srt_list,
    nfeatures = 3, assay = c("peaks", "peaks"), verbose = FALSE
  )
  expect_identical(result$HVF, native)
  expect_identical(result$HVF[[1]], ids[[2]])
  expect_setequal(result$HVF, ids[1:3])
  filtered <- check_atac_hvf_rank_objects(list(first, second), n = 12, minimum = 2, find = FALSE)
  expect_identical(filtered$HVF, native)
  capped <- check_atac_hvf_rank_objects(list(first, second), n = 12, find = FALSE)
  expect_length(capped$HVF, 4)
  expect_identical(capped$HVF, Seurat::SelectIntegrationFeatures(
    capped$srt_list,
    nfeatures = 4, assay = c("peaks", "peaks"), verbose = FALSE
  ))
  expect_error(
    check_atac_hvf_rank_objects(list(first, second), minimum = 3, find = FALSE),
    "[Nn]o HVF available"
  )
})

test_that("zero chromatin feature cap retains explicit empty-selection error", {
  object <- make_atac_hvf_rank_object()
  expect_error(check_atac_hvf_rank_objects(list(object), n = 0), "[Nn]o HVF available")
})
