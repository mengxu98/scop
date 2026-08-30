data("pbmc_small", package = "SeuratObject")
srt_base <- suppressMessages(NormalizeData(pbmc_small, verbose = FALSE))
srt_vf <- suppressMessages(FindVariableFeatures(srt_base, verbose = FALSE))
Idents(srt_vf) <- srt_vf$groups
flat <- function(m) {
  as.numeric(if (inherits(m, "sparseMatrix")) methods::as(m, "matrix") else m)
}

test_that("NormalizeData covers CLR and RC like native Seurat", {
  skip_if_not_installed("Seurat")
  for (spec in list(
    c("LogNormalize", 1L),
    c("CLR", 2L),
    c("CLR", 1L),
    c("RC", 1L),
    c("RC", 2L)
  )) {
    method <- spec[1]
    margin <- as.integer(spec[2])
    o1 <- suppressMessages(NormalizeData(
      srt_base,
      normalization.method = method,
      margin = margin,
      verbose = FALSE
    ))
    o2 <- suppressWarnings(seurat_reference_method(
      "NormalizeData",
      "Seurat",
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

test_that("matrix default preprocessing methods retain Seurat parity", {
  set.seed(20260903)
  matrix_input <- Matrix::rsparsematrix(100, 60, density = 0.12)
  matrix_input@x <- stats::rpois(length(matrix_input@x), lambda = 3) + 1
  rownames(matrix_input) <- paste0("g", seq_len(nrow(matrix_input)))
  colnames(matrix_input) <- paste0("c", seq_len(ncol(matrix_input)))

  for (method in c("LogNormalize", "CLR", "RC")) {
    for (margin in if (identical(method, "CLR")) c(1L, 2L) else 1L) {
      actual <- NormalizeData(matrix_input, normalization.method = method, margin = margin)
      expected <- seurat_reference_method(
        "NormalizeData",
        "default",
        matrix_input,
        normalization.method = method,
        margin = margin
      )
      expect_equal(as.matrix(actual), as.matrix(expected), tolerance = 1e-12)
    }
  }

  actual_vf <- suppressWarnings(FindVariableFeatures(
    matrix_input,
    nfeatures = 30,
    verbose = FALSE
  ))
  expected_vf <- suppressWarnings(seurat_reference_method(
    "FindVariableFeatures",
    "default",
    matrix_input,
    nfeatures = 30,
    verbose = FALSE
  ))
  expect_equal(actual_vf, expected_vf, tolerance = 1e-10)

  normalized <- seurat_reference_method(
    "NormalizeData",
    "default",
    matrix_input,
    verbose = FALSE
  )
  for (args in list(
    list(),
    list(do.center = FALSE),
    list(do.scale = FALSE),
    list(do.center = 1),
    list(do.scale = 1),
    list(use.umi = 1),
    list(split.by = rep(c("a", "b"), each = 30))
  )) {
    actual <- do.call(ScaleData, c(list(object = normalized, verbose = FALSE), args))
    expected <- do.call(
      function(object, ...) seurat_reference_method("ScaleData", "default", object, ...),
      c(list(object = normalized, verbose = FALSE), args)
    )
    expect_equal(actual, expected, tolerance = 1e-8)
  }

  pca_input <- as.matrix(ScaleData(normalized, verbose = FALSE))
  for (args in list(
    list(weight.by.var = 1),
    list(rev.pca = 1),
    list(approx = 1),
    list(npcs = 5.5)
  )) {
    actual <- suppressWarnings(do.call(
      RunPCA,
      utils::modifyList(list(object = pca_input, npcs = 5, verbose = FALSE), args)
    ))
    expected <- suppressWarnings(do.call(
      function(object, ...) seurat_reference_method("RunPCA", "default", object, ...),
      utils::modifyList(list(object = pca_input, npcs = 5, verbose = FALSE), args)
    ))
    expect_identical(actual@cell.embeddings, expected@cell.embeddings)
    expect_identical(actual@stdev, expected@stdev)
  }
})

test_that("default methods do not recurse through scop S3 registration", {
  set.seed(20260830)
  x <- Matrix::rsparsematrix(30, 20, density = 0.2)
  x@x <- stats::rpois(length(x@x), lambda = 2) + 1
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  colnames(x) <- paste0("c", seq_len(ncol(x)))
  expect_no_error(FindVariableFeatures(x, nfeatures = 10, verbose = FALSE))
  expect_no_error(FindMarkers(
    x,
    cells.1 = colnames(x)[1:10],
    cells.2 = colnames(x)[11:20],
    logfc.threshold = 0,
    min.pct = 0,
    verbose = FALSE
  ))
})

test_that("FindMarkers delegates non-logical only.pos exactly", {
  args <- list(
    object = srt_vf,
    ident.1 = levels(SeuratObject::Idents(srt_vf))[[1]],
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = 1,
    verbose = FALSE
  )
  actual <- suppressWarnings(do.call(FindMarkers, args))
  expected <- suppressWarnings(do.call(
    function(object, ...) seurat_reference_method("FindMarkers", "Seurat", object, ...),
    args
  ))
  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("NormalizeData delegates Assay5 layer and save overrides", {
  set.seed(20260830)
  counts <- Matrix::rsparsematrix(30, 20, density = 0.2)
  counts@x <- abs(round(counts@x * 5)) + 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  object <- Seurat::CreateSeuratObject(counts)
  object$batch <- rep(c("a", "b"), each = 10)
  object[["RNA"]] <- split(object[["RNA"]], f = object$batch)

  actual <- NormalizeData(
    object,
    layer = "counts.a",
    save = "custom",
    verbose = FALSE
  )
  expected <- seurat_reference_method(
    "NormalizeData",
    "Seurat",
    object,
    layer = "counts.a",
    save = "custom",
    verbose = FALSE
  )
  expect_setequal(
    SeuratObject::Layers(actual[["RNA"]], search = NA),
    SeuratObject::Layers(expected[["RNA"]], search = NA)
  )
  expect_equal(
    as.matrix(SeuratObject::LayerData(actual, layer = "custom")),
    as.matrix(SeuratObject::LayerData(expected, layer = "custom")),
    tolerance = 1e-12
  )
})

test_that("ScaleData supports regression, split.by, and use.umi", {
  args_list <- list(
    list(vars.to.regress = "nCount_RNA"),
    list(split.by = "groups"),
    list(use.umi = TRUE),
    list(vars.to.regress = "nCount_RNA", use.umi = TRUE),
    list(do.scale = FALSE),
    list(do.center = FALSE),
    list(do.scale = FALSE, do.center = FALSE),
    list(block.size = 17L, min.cells.to.block = 1L),
    list(scale.max = 3)
  )
  for (args in args_list) {
    o1 <- do.call(
      get("ScaleData.Seurat", asNamespace("scop")),
      c(list(object = srt_vf, verbose = FALSE), args)
    )
    o2 <- do.call(
      function(object, ...) seurat_reference_method("ScaleData", "Seurat", object, ...),
      c(list(object = srt_vf, verbose = FALSE), args)
    )
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

  actual <- get("ScaleData.Seurat", asNamespace("scop"))(
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

test_that("ScaleData delegates non-linear regression models", {
  features <- head(SeuratObject::VariableFeatures(srt_vf), 20L)
  actual <- suppressWarnings(ScaleData(
    srt_vf,
    features = features,
    vars.to.regress = "nCount_RNA",
    model.use = "poisson",
    verbose = FALSE
  ))
  expected <- suppressWarnings(seurat_reference_method(
    "ScaleData",
    "Seurat",
    srt_vf,
    features = features,
    vars.to.regress = "nCount_RNA",
    model.use = "poisson",
    verbose = FALSE
  ))
  expect_equal(
    as.matrix(SeuratObject::LayerData(actual, layer = "scale.data")),
    as.matrix(SeuratObject::LayerData(expected, layer = "scale.data")),
    tolerance = 1e-12
  )
})

test_that("FindVariableFeatures supports disp and mean.var.plot", {
  for (method in c("disp", "mean.var.plot")) {
    a_obj <- FindVariableFeatures(srt_base, selection.method = method, verbose = FALSE)
    b_obj <- suppressWarnings(seurat_reference_method(
      "FindVariableFeatures",
      "Seurat",
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

test_that("FindVariableFeatures covers VST clipping and equal-frequency MVP", {
  vst_actual <- FindVariableFeatures(
    srt_base,
    selection.method = "vst",
    clip.max = 2,
    nfeatures = 30L,
    verbose = FALSE
  )
  vst_expected <- seurat_reference_method(
    "FindVariableFeatures",
    "Seurat",
    srt_base,
    selection.method = "vst",
    clip.max = 2,
    nfeatures = 30L,
    verbose = FALSE
  )
  expect_identical(
    SeuratObject::VariableFeatures(vst_actual),
    SeuratObject::VariableFeatures(vst_expected)
  )

  mvp_actual <- FindVariableFeatures(
    srt_base,
    selection.method = "mean.var.plot",
    binning.method = "equal_frequency",
    num.bin = 10L,
    verbose = FALSE
  )
  mvp_expected <- suppressWarnings(seurat_reference_method(
    "FindVariableFeatures",
    "Seurat",
    srt_base,
    selection.method = "mean.var.plot",
    binning.method = "equal_frequency",
    num.bin = 10L,
    verbose = FALSE
  ))
  expect_identical(
    SeuratObject::VariableFeatures(mvp_actual),
    SeuratObject::VariableFeatures(mvp_expected)
  )
})

test_that("FindVariableFeatures accelerates explicit VST layer overrides", {
  set.seed(20260830)
  counts <- Matrix::rsparsematrix(60, 30, density = 0.2)
  counts@x <- abs(round(counts@x * 5)) + 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  object <- Seurat::CreateSeuratObject(counts)
  object$batch <- rep(c("a", "b"), each = 15)
  object[["RNA"]] <- split(object[["RNA"]], f = object$batch)

  actual <- suppressWarnings(FindVariableFeatures(
    object,
    selection.method = "vst",
    layer = "counts.a",
    nfeatures = 15L,
    verbose = FALSE
  ))
  expected <- suppressWarnings(seurat_reference_method(
    "FindVariableFeatures",
    "Seurat",
    object,
    selection.method = "vst",
    layer = "counts.a",
    nfeatures = 15L,
    verbose = FALSE
  ))
  expect_identical(
    SeuratObject::VariableFeatures(actual),
    SeuratObject::VariableFeatures(expected)
  )
  expect_equal(
    SeuratObject::HVFInfo(actual[["RNA"]], method = "vst", layer = "counts.a"),
    SeuratObject::HVFInfo(expected[["RNA"]], method = "vst", layer = "counts.a"),
    tolerance = 1e-10
  )
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
    expected <- suppressWarnings(seurat_reference_method(
      "FindVariableFeatures",
      "Seurat",
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

  variable_object <- suppressWarnings(seurat_reference_method(
    "FindVariableFeatures",
    "Seurat",
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
  expected_scaled <- seurat_reference_method(
    "ScaleData",
    "Seurat",
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

  result <- get("RunCCA.default", asNamespace("scop"))(
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
