data("pbmc_small", package = "SeuratObject")
srt_base <- suppressMessages(NormalizeData(pbmc_small, verbose = FALSE))
srt_vf <- suppressMessages(FindVariableFeatures(srt_base, verbose = FALSE))
Idents(srt_vf) <- srt_vf$groups
flat <- function(m) {
  as.numeric(if (inherits(m, "sparseMatrix")) methods::as(m, "matrix") else m)
}

test_that("NormalizeData covers CLR and RC like native Seurat", {
  skip_if_not_installed("Seurat")
  for (spec in list(c("LogNormalize", 1L), c("CLR", 2L), c("CLR", 1L), c("RC", 1L))) {
    method <- spec[1]; margin <- as.integer(spec[2])
    o1 <- suppressMessages(NormalizeData(
      srt_base,
      normalization.method = method,
      margin = margin,
      verbose = FALSE
    ))
    o2 <- suppressWarnings(Seurat::NormalizeData(
      srt_base,
      normalization.method = method,
      margin = margin,
      verbose = FALSE
    ))
    d1 <- flat(SeuratObject::GetAssayData(o1[["RNA"]], layer = "data"))
    d2 <- flat(SeuratObject::GetAssayData(o2[["RNA"]], layer = "data"))
    expect_lt(max(abs(d1 - d2)), 1e-8)
  }
})

test_that("ScaleData supports regression, split.by, and use.umi", {
  args_list <- list(
    list(vars.to.regress = "nCount_RNA"),
    list(split.by = "groups"),
    list(use.umi = TRUE),
    list(vars.to.regress = "nCount_RNA", use.umi = TRUE),
    list(do.scale = FALSE),
    list(scale.max = 3)
  )
  for (args in args_list) {
    o1 <- do.call(ScaleData.Seurat, c(list(object = srt_vf, verbose = FALSE), args))
    o2 <- do.call(Seurat::ScaleData, c(list(object = srt_vf, verbose = FALSE), args))
    A <- as.matrix(SeuratObject::GetAssayData(o1[["RNA"]], layer = "scale.data"))
    B <- as.matrix(SeuratObject::GetAssayData(o2[["RNA"]], layer = "scale.data"))
    common <- intersect(rownames(A), rownames(B))
    cells <- intersect(colnames(A), colnames(B))
    expect_lt(max(abs(A[common, cells] - B[common, cells])), 1e-6)
  }
})

test_that("FindVariableFeatures supports disp and mean.var.plot", {
  for (method in c("disp", "mean.var.plot")) {
    a_obj <- FindVariableFeatures(srt_base, selection.method = method, verbose = FALSE)
    b_obj <- suppressWarnings(Seurat::FindVariableFeatures(
      srt_base,
      selection.method = method,
      verbose = FALSE
    ))
    va <- SeuratObject::VariableFeatures(a_obj)
    vb <- SeuratObject::VariableFeatures(b_obj)
    jaccard <- length(intersect(va, vb)) /
      max(length(union(va, vb)), 1L)
    expect_gt(jaccard, 0.99)
  }
})

test_that("FoldChange dispatches to native behavior on Seurat objects", {
  o1 <- FoldChange(srt_vf, ident.1 = "g1")
  o2 <- Seurat::FoldChange(srt_vf, ident.1 = "g1")
  if (is.list(o1) && is.list(o2)) {
    expect_equal(o1[[1]], o2[[1]])
  } else {
    expect_equal(o1[[1]], o2[[1]])
  }
})
