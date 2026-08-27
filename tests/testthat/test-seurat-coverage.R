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

test_that("ScaleData accepts explicit latent.data covariates", {
  latent <- data.frame(
    custom_covariate = seq_len(ncol(srt_vf)),
    row.names = colnames(srt_vf)
  )

  actual <- ScaleData.Seurat(
    srt_vf,
    vars.to.regress = "custom_covariate",
    latent.data = latent,
    verbose = FALSE
  )
  expected_assay <- Seurat:::ScaleData.Assay(
    srt_vf[["RNA"]],
    vars.to.regress = "custom_covariate",
    latent.data = latent,
    verbose = FALSE
  )

  expect_equal(
    as.matrix(SeuratObject::GetAssayData(actual[["RNA"]], layer = "scale.data")),
    as.matrix(SeuratObject::GetAssayData(expected_assay, layer = "scale.data")),
    tolerance = 1e-8
  )
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

test_that("layered preprocessing aligns partially overlapping features", {
  set.seed(42)
  counts_a <- Matrix::Matrix(
    matrix(stats::rpois(80 * 30, lambda = 2), nrow = 80),
    sparse = TRUE,
    dimnames = list(paste0("g", 1:80), paste0("a", 1:30))
  )
  counts_b <- Matrix::Matrix(
    matrix(stats::rpois(80 * 30, lambda = 2), nrow = 80),
    sparse = TRUE,
    dimnames = list(paste0("g", 41:120), paste0("b", 1:30))
  )
  object <- suppressWarnings(suppressMessages(
    Seurat::CreateSeuratObject(counts = list(A = counts_a, B = counts_b))
  ))
  object <- suppressMessages(NormalizeData(object, verbose = FALSE))

  for (method in c("disp", "mean.var.plot")) {
    actual <- suppressWarnings(FindVariableFeatures(
      object,
      selection.method = method,
      nfeatures = 20,
      verbose = FALSE
    ))
    expected <- suppressWarnings(Seurat::FindVariableFeatures(
      object,
      selection.method = method,
      nfeatures = 20,
      verbose = FALSE
    ))
    expect_identical(
      SeuratObject::VariableFeatures(actual),
      SeuratObject::VariableFeatures(expected)
    )
  }

  variable_object <- suppressWarnings(Seurat::FindVariableFeatures(
    object,
    selection.method = "vst",
    nfeatures = 20,
    verbose = FALSE
  ))
  actual_scaled <- ScaleData(
    variable_object,
    features = SeuratObject::VariableFeatures(variable_object),
    verbose = FALSE
  )
  expected_scaled <- Seurat::ScaleData(
    variable_object,
    features = SeuratObject::VariableFeatures(variable_object),
    verbose = FALSE
  )
  expect_equal(
    as.matrix(SeuratObject::LayerData(actual_scaled, layer = "scale.data")),
    as.matrix(SeuratObject::LayerData(expected_scaled, layer = "scale.data")),
    tolerance = 1e-8
  )
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

test_that("RunCCA.default reaches the registered matrix helpers", {
  set.seed(9)
  object1 <- matrix(
    stats::rnorm(12 * 8),
    nrow = 12,
    dimnames = list(paste0("g", 1:12), paste0("a", 1:8))
  )
  object2 <- matrix(
    stats::rnorm(12 * 7),
    nrow = 12,
    dimnames = list(paste0("g", 1:12), paste0("b", 1:7))
  )

  result <- RunCCA.default(
    object1,
    object2,
    num.cc = 3,
    seed.use = 9,
    verbose = FALSE
  )

  expect_identical(dim(result$ccv), c(15L, 3L))
  expect_true(all(is.finite(result$ccv)))
  expect_true(all(is.finite(result$d)))
})
