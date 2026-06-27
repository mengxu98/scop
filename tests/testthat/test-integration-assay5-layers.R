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

test_that("CheckDataMerge reuses log-normalized RNA merge without changing HVF output", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  set.seed(2)
  counts <- Matrix::rsparsematrix(60, 40, density = 0.25)
  counts@x <- abs(round(counts@x * 10)) + 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))

  srt <- Seurat::CreateSeuratObject(counts)
  srt$batch <- rep(c("a", "b"), each = 20)
  srt <- NormalizeData(srt, verbose = FALSE)

  checked <- CheckDataMerge(
    srt_merge = srt,
    batch = "batch",
    nHVF = 15,
    verbose = FALSE
  )

  expect_s4_class(checked$srt_merge, "Seurat")
  expect_equal(dim(checked$srt_merge), dim(srt))
  expect_equal(colnames(checked$srt_merge), colnames(srt))
  expect_length(checked$HVF, 15)
  expect_equal(SeuratObject::VariableFeatures(checked$srt_merge), checked$HVF)
})

test_that("NormalizeData handles split Assay5 counts layers", {
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
  assay <- SeuratObject::DefaultAssay(srt)
  srt[[assay]] <- split(srt[[assay]], f = srt$batch)

  expect_setequal(
    SeuratObject::Layers(srt[[assay]], search = "counts"),
    c("counts.a", "counts.b")
  )
  expect_no_error(srt <- NormalizeData(srt, verbose = FALSE))
  expect_setequal(
    SeuratObject::Layers(srt[[assay]], search = "data"),
    c("data.a", "data.b")
  )
  expect_true(all(Matrix::colSums(SeuratObject::LayerData(srt, layer = "data.a")) > 0))
})

test_that("FindVariableFeatures handles split Assay5 counts layers", {
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
  assay <- SeuratObject::DefaultAssay(srt)
  srt[[assay]] <- split(srt[[assay]], f = srt$batch)

  expect_no_error(
    srt <- FindVariableFeatures(srt, nfeatures = 20, verbose = FALSE)
  )
  expect_length(SeuratObject::VariableFeatures(srt), 20)
})

test_that("ScaleData handles split Assay5 data layers", {
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
  assay <- SeuratObject::DefaultAssay(srt)
  srt_single <- srt
  srt_single <- NormalizeData(srt_single, verbose = FALSE)
  srt_single <- FindVariableFeatures(srt_single, nfeatures = 20, verbose = FALSE)
  srt_single <- ScaleData(
    srt_single,
    features = SeuratObject::VariableFeatures(srt_single),
    verbose = FALSE
  )

  srt[[assay]] <- split(srt[[assay]], f = srt$batch)
  srt <- NormalizeData(srt, verbose = FALSE)
  srt <- FindVariableFeatures(srt, nfeatures = 20, verbose = FALSE)

  expect_no_error(
    srt <- ScaleData(
      srt,
      features = SeuratObject::VariableFeatures(srt),
      verbose = FALSE
    )
  )
  expect_true("scale.data" %in% SeuratObject::Layers(srt[[assay]], search = NA))
  expect_equal(
    nrow(SeuratObject::LayerData(srt, layer = "scale.data")),
    20L
  )
  expect_equal(
    SeuratObject::VariableFeatures(srt),
    SeuratObject::VariableFeatures(srt_single)
  )
  expect_equal(
    as.matrix(SeuratObject::LayerData(srt, layer = "scale.data")),
    as.matrix(SeuratObject::LayerData(srt_single, layer = "scale.data")),
    tolerance = 1e-12
  )
})

test_that("srt_to_adata stacks split Assay5 layers into sparse X", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("reticulate")
  skip_if_not(reticulate::py_module_available("anndata"))
  skip_if_not(reticulate::py_module_available("scipy"))

  old_skip_python_prepare <- getOption("scop_skip_python_prepare", FALSE)
  options(scop_skip_python_prepare = TRUE)
  on.exit(options(scop_skip_python_prepare = old_skip_python_prepare), add = TRUE)

  set.seed(1)
  counts <- Matrix::rsparsematrix(20, 12, density = 0.2)
  counts@x <- abs(round(counts@x * 10)) + 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))

  srt <- Seurat::CreateSeuratObject(counts)
  srt$batch <- rep(c("a", "b"), each = 6)
  assay <- SeuratObject::DefaultAssay(srt)
  srt[[assay]] <- split(srt[[assay]], f = srt$batch)

  layers <- SeuratObject::Layers(srt[[assay]], search = "counts")
  expect_setequal(layers, c("counts.a", "counts.b"))

  adata <- srt_to_adata(
    srt,
    assay_x = assay,
    layer_x = "counts",
    verbose = FALSE
  )

  expect_identical(reticulate::py_to_r(adata$X$format), "csr")
  expect_equal(
    unlist(reticulate::py_to_r(adata$X$shape)),
    c(ncol(srt), nrow(srt[[assay]]))
  )
  expect_identical(
    reticulate::py_to_r(adata$obs_names$to_list()),
    colnames(srt)
  )
  expect_equal(nrow(reticulate::py_to_r(adata$obs)), ncol(srt))
})
