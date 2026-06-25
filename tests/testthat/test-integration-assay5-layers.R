test_that("legacy integration scaling handles split Assay5 data layers", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  set.seed(1)
  counts <- Matrix::rsparsematrix(80, 60, density = 0.2)
  counts@x <- abs(round(counts@x * 10)) + 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))

  srt <- Seurat::CreateSeuratObject(counts)
  srt$batch <- rep(c("a", "b"), each = 30)
  srt_list <- Seurat::SplitObject(srt, split.by = "batch")
  srt_list <- lapply(srt_list, function(x) {
    NormalizeData(x, verbose = FALSE)
  })
  merged <- Reduce(merge, srt_list)
  assay <- SeuratObject::DefaultAssay(merged)

  expect_false("data" %in% SeuratObject::Layers(merged[[assay]], search = NA))
  expect_gt(length(SeuratObject::Layers(merged[[assay]], search = "data")), 1)

  expect_no_error(
    suppressWarnings(
      integration_scop(
        srt_merge = srt,
        batch = "batch",
        integration_method = "Uncorrected",
        nHVF = 20,
        linear_reduction_dims = 5,
        linear_reduction_dims_use = 1:3,
        nonlinear_reduction = NULL,
        neighbor_k = 5,
        verbose = FALSE
      )
    )
  )
})
