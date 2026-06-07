test_that("RunDimsReduction tolerates SCT scale.data missing some HVF", {
  set.seed(1)
  counts1 <- matrix(sample(0:5, 40, replace = TRUE), nrow = 5)
  counts2 <- matrix(sample(0:5, 40, replace = TRUE), nrow = 5)
  rownames(counts1) <- rownames(counts2) <- paste0("Gene", 1:5)
  colnames(counts1) <- paste0("A", 1:8)
  colnames(counts2) <- paste0("B", 1:8)

  s1 <- Seurat::CreateSeuratObject(counts = Matrix::Matrix(counts1, sparse = TRUE))
  s2 <- Seurat::CreateSeuratObject(counts = Matrix::Matrix(counts2, sparse = TRUE))
  s1$batch <- "A"
  s2$batch <- "B"

  srt <- merge(s1, s2)
  srt[["RNA"]] <- SeuratObject::JoinLayers(srt[["RNA"]])
  srt[["RNA"]] <- split(
    x = srt[["RNA"]],
    f = srt[["batch", drop = TRUE]]
  )
  srt <- Seurat::SCTransform(
    object = srt,
    assay = "RNA",
    variable.features.n = 2000,
    verbose = FALSE
  )

  hvf <- SeuratObject::VariableFeatures(srt, assay = "SCT", nfeatures = 2000)
  scale_features <- rownames(scop::GetAssayData5(srt, assay = "SCT", layer = "scale.data"))
  expect_true(length(setdiff(hvf, scale_features)) > 0)

  expect_silent(
    suppressWarnings(
      scop::RunDimsReduction(
        srt = srt,
        features = hvf,
        assay = "SCT",
        layer = "scale.data",
        linear_reduction = "pca",
        linear_reduction_dims = 2,
        force_linear_reduction = TRUE,
        verbose = FALSE,
        seed = 11
      )
    )
  )
})
