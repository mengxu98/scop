test_that("RunDEtest uses scop FindMarkers-compatible sparse Wilcoxon path", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  set.seed(1)
  counts <- Matrix::Matrix(rpois(80 * 30, lambda = 1), nrow = 80, sparse = TRUE)
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  SeuratObject::Idents(srt) <- rep(c("A", "B"), each = 15)

  cells1 <- colnames(srt)[1:15]
  cells2 <- colnames(srt)[16:30]
  features <- rownames(srt)[1:50]

  seurat_fm <- Seurat::FindMarkers(
    srt,
    ident.1 = "A",
    ident.2 = "B",
    features = features,
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = FALSE,
    verbose = FALSE
  )
  scop_fm <- scop::FindMarkers(
    srt,
    ident.1 = "A",
    ident.2 = "B",
    features = features,
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = FALSE,
    verbose = FALSE
  )
  fast_fm <- scop::FindMarkers(
    srt,
    .scop_backend = "sparse_wilcox",
    cells.1 = cells1,
    cells.2 = cells2,
    features = features,
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = FALSE
  )
  out <- RunDEtest(
    srt,
    cells1 = cells1,
    cells2 = cells2,
    features = features,
    fc.threshold = 1,
    min.pct = 0,
    only.pos = FALSE,
    verbose = FALSE
  )
  markers <- out@tools$DEtest_custom$AllMarkers_wilcox

  expect_equal(scop_fm, seurat_fm)
  expect_error(
    scop::FindMarkers(
      srt,
      cells.1 = cells1,
      cells.2 = cells2,
      features = features,
      logfc.threshold = 0,
      min.pct = 0,
      only.pos = FALSE,
      verbose = FALSE
    ),
    "At least 1 ident"
  )
  expect_identical(rownames(fast_fm), markers$gene)
  expect_equal(fast_fm$p_val, markers$p_val, tolerance = 1e-12)
  expect_equal(fast_fm$avg_log2FC, markers$avg_log2FC, tolerance = 0)
  expect_equal(fast_fm$pct.1, markers$pct.1, tolerance = 0)
  expect_equal(fast_fm$pct.2, markers$pct.2, tolerance = 0)
})
