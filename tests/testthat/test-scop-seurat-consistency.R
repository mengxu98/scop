make_consistency_object <- function(seed = 21, n_genes = 60, n_cells = 50) {
  set.seed(seed)
  mat <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.15)
  mat@x <- abs(mat@x * 10)
  rownames(mat) <- paste0("g", seq_len(nrow(mat)))
  colnames(mat) <- paste0("c", seq_len(ncol(mat)))
  obj <- Seurat::CreateSeuratObject(counts = mat)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 40, verbose = FALSE)
  obj <- ScaleData(obj, features = SeuratObject::VariableFeatures(obj), verbose = FALSE)
  obj
}

test_that("RunPCA default output is identical to Seurat", {
  obj <- make_consistency_object()
  hvf <- SeuratObject::VariableFeatures(obj)

  scop_pca <- RunPCA(
    obj,
    features = hvf,
    npcs = 8,
    reduction.name = "pca_scop",
    seed.use = 42,
    verbose = FALSE
  )
  seurat_pca <- seurat_reference_method(
    "RunPCA",
    "Seurat",
    obj,
    features = hvf,
    npcs = 8,
    reduction.name = "pca_seurat",
    seed.use = 42,
    verbose = FALSE
  )

  expect_identical(
    SeuratObject::Embeddings(scop_pca[["pca_scop"]]),
    SeuratObject::Embeddings(seurat_pca[["pca_seurat"]])
  )
})

test_that("RunPCA cpp backend remains an explicit override", {
  obj <- make_consistency_object()
  hvf <- SeuratObject::VariableFeatures(obj)
  expect_no_error(
    RunPCA(obj, features = hvf, npcs = 8, backend = "cpp", verbose = FALSE)
  )
})

test_that("RunPCA native auto gate covers shape and argument branches", {
  gate <- get("pca_native_auto_supported", asNamespace("scop"))
  x <- matrix(0, nrow = 100L, ncol = 800L)
  expect_true(gate(x, npcs = 20L, extra_args = list()))
  expect_false(gate(x[, seq_len(799L), drop = FALSE], 20L, list()))
  expect_false(gate(x, 20L, list(tol = 1e-5)))
  expect_false(gate(Matrix::Matrix(x, sparse = TRUE), 20L, list()))
  expect_false(gate(matrix(0, nrow = 5001L, ncol = 1L), 20L, list()))
  expect_false(gate(x, 0L, list()))
})

test_that("RunPCA native eigengap gate rejects unstable trailing PCs", {
  acceptable <- get(
    "pca_native_candidate_acceptable",
    asNamespace("scop")
  )
  expect_true(acceptable(c(10, 8, 6, 5), npcs = 3L))
  expect_false(acceptable(c(10, 8, 6, 5.99), npcs = 3L))
  expect_false(acceptable(c(10, 8, NA_real_, 5), npcs = 3L))
})

test_that("RunPCA matrix method supplies Seurat's default assay contract", {
  set.seed(20260830)
  x <- matrix(stats::rnorm(30L * 20L), nrow = 30L)
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  colnames(x) <- paste0("c", seq_len(ncol(x)))
  result <- RunPCA(x, npcs = 5L, verbose = FALSE)
  expect_s4_class(result, "DimReduc")
  expect_identical(result@assay.used, "RNA")
})

test_that("RunPCA matrix fallback retains exact and reverse PCA branches", {
  set.seed(20260905)
  x <- matrix(stats::rnorm(40L * 30L), nrow = 40L)
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  colnames(x) <- paste0("c", seq_len(ncol(x)))
  for (args in list(
    list(approx = FALSE),
    list(rev.pca = TRUE, approx = FALSE),
    list(weight.by.var = FALSE, approx = FALSE)
  )) {
    actual <- suppressWarnings(do.call(
      RunPCA,
      c(list(object = x, npcs = 5L, verbose = FALSE), args)
    ))
    expected <- suppressWarnings(do.call(
      function(object, ...) seurat_reference_method("RunPCA", "default", object, ...),
      c(list(object = x, npcs = 5L, verbose = FALSE), args)
    ))
    expect_equal(
      SeuratObject::Embeddings(actual),
      SeuratObject::Embeddings(expected),
      tolerance = 1e-12
    )
    expect_equal(
      SeuratObject::Loadings(actual),
      SeuratObject::Loadings(expected),
      tolerance = 1e-12
    )
  }

  sparse <- methods::as(Matrix::Matrix(x, sparse = TRUE), "dgCMatrix")
  expect_error(
    RunPCA(sparse, npcs = 5L, verbose = FALSE),
    "unused argument"
  )
  expect_no_error(RunPCA(
    sparse,
    npcs = 5L,
    backend = "cpp",
    verbose = FALSE
  ))
})

test_that("RunPCA native auto path preserves a separated PCA subspace", {
  set.seed(20260830)
  features <- 80L
  cells <- 800L
  latent <- 12L
  x <- matrix(stats::rnorm(features * latent), nrow = features) %*%
    matrix(stats::rnorm(latent * cells), nrow = latent)
  x <- x + matrix(stats::rnorm(length(x), sd = 0.05), nrow = features)
  x <- x - rowMeans(x)
  rownames(x) <- paste0("g", seq_len(features))
  colnames(x) <- paste0("c", seq_len(cells))
  reference <- RunPCA(
    x,
    assay = "RNA",
    npcs = 10L,
    backend = "irlba",
    seed.use = 42,
    verbose = FALSE
  )
  candidate <- RunPCA(
    x,
    assay = "RNA",
    npcs = 10L,
    seed.use = 42,
    verbose = FALSE
  )
  reference_embedding <- SeuratObject::Embeddings(reference)
  candidate_embedding <- SeuratObject::Embeddings(candidate)
  signs <- sign(colSums(reference_embedding * candidate_embedding))
  signs[signs == 0] <- 1
  candidate_embedding <- sweep(candidate_embedding, 2L, signs, `*`)
  expect_lt(
    sqrt(mean((candidate_embedding - reference_embedding)^2)),
    1e-8
  )
  expect_equal(candidate@stdev, reference@stdev, tolerance = 1e-10)
})

test_that("FindNeighbors nn graph is identical to Seurat (self edges preserved)", {
  obj <- make_consistency_object()
  hvf <- SeuratObject::VariableFeatures(obj)
  obj <- RunPCA(obj, features = hvf, npcs = 8, verbose = FALSE)

  scop_nn <- FindNeighbors(
    obj,
    reduction = "pca",
    dims = 1:8,
    k.param = 10,
    graph.name = c("nn_scop", "snn_scop"),
    verbose = FALSE
  )
  seurat_nn <- seurat_reference_method(
    "FindNeighbors",
    "Seurat",
    obj,
    reduction = "pca",
    dims = 1:8,
    k.param = 10,
    graph.name = c("nn_seurat", "snn_seurat"),
    verbose = FALSE
  )

  expect_identical(
    as.matrix(scop_nn@graphs[["nn_scop"]]),
    as.matrix(seurat_nn@graphs[["nn_seurat"]])
  )
  expect_equal(
    sum(Matrix::diag(scop_nn@graphs[["nn_scop"]])),
    ncol(obj),
    label = "nn graph keeps self edges like Seurat"
  )
})

test_that("FindNeighbors SNN graph is identical to Seurat", {
  obj <- make_consistency_object()
  hvf <- SeuratObject::VariableFeatures(obj)
  obj <- RunPCA(obj, features = hvf, npcs = 8, verbose = FALSE)

  scop_nn <- FindNeighbors(
    obj,
    reduction = "pca",
    dims = 1:8,
    k.param = 10,
    graph.name = c("nn_scop", "snn_scop"),
    verbose = FALSE
  )
  seurat_nn <- seurat_reference_method(
    "FindNeighbors",
    "Seurat",
    obj,
    reduction = "pca",
    dims = 1:8,
    k.param = 10,
    graph.name = c("nn_seurat", "snn_seurat"),
    verbose = FALSE
  )

  expect_identical(
    as.matrix(scop_nn@graphs[["snn_scop"]]),
    as.matrix(seurat_nn@graphs[["snn_seurat"]])
  )
})

test_that("VST on split Assay5 layers matches Seurat's layered reference", {
  obj <- make_consistency_object()
  obj[["RNA"]] <- split(
    obj[["RNA"]],
    f = rep(c("A", "B"), length.out = ncol(obj))
  )
  joined <- SeuratObject::JoinLayers(obj, assay = "RNA")

  layered_hvf <- SeuratObject::VariableFeatures(
    FindVariableFeatures(obj, selection.method = "vst", nfeatures = 40, verbose = FALSE)
  )
  seurat_layered_hvf <- SeuratObject::VariableFeatures(
    seurat_reference_method(
      "FindVariableFeatures",
      "Seurat",
      obj,
      selection.method = "vst",
      nfeatures = 40,
      verbose = FALSE
    )
  )
  joined_hvf <- SeuratObject::VariableFeatures(
    FindVariableFeatures(joined, selection.method = "vst", nfeatures = 40, verbose = FALSE)
  )

  # scop follows Seurat's per-layer VST semantics: the layered HVF must be
  # byte-identical to the reference implementation on the same layered input.
  expect_identical(layered_hvf, seurat_layered_hvf)
  # Seurat's layered behaviour itself differs from the joined reference
  # (per-layer vst + consensus merge vs one-shot joined vst); scop mirrors
  # that behaviour, so the layered result may differ from the joined one.
  expect_identical(
    identical(layered_hvf, joined_hvf),
    identical(seurat_layered_hvf, joined_hvf)
  )
})

test_that("VST aligns Assay5 layers with partially overlapping features", {
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
  obj <- suppressWarnings(suppressMessages(
    Seurat::CreateSeuratObject(counts = list(A = counts_a, B = counts_b))
  ))

  scop_result <- suppressWarnings(FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 20,
    verbose = FALSE
  ))
  seurat_result <- suppressWarnings(seurat_reference_method(
    "FindVariableFeatures",
    "Seurat",
    obj,
    selection.method = "vst",
    nfeatures = 20,
    verbose = FALSE
  ))

  expect_identical(
    SeuratObject::VariableFeatures(scop_result),
    SeuratObject::VariableFeatures(seurat_result)
  )
  expect_identical(
    rownames(scop_result[[SeuratObject::DefaultAssay(scop_result)]]),
    rownames(seurat_result[[SeuratObject::DefaultAssay(seurat_result)]])
  )
  scop_assay <- scop_result[[SeuratObject::DefaultAssay(scop_result)]]
  seurat_assay <- seurat_result[[SeuratObject::DefaultAssay(seurat_result)]]
  scop_meta <- methods::slot(scop_assay, "meta.data")
  seurat_meta <- methods::slot(seurat_assay, "meta.data")
  hvf_columns <- grep("^vf_vst_", colnames(scop_meta), value = TRUE)
  expect_identical(rownames(scop_meta), rownames(seurat_meta))
  expect_identical(
    unname(scop_meta[, hvf_columns, drop = FALSE]),
    unname(seurat_meta[, hvf_columns, drop = FALSE])
  )
  expect_true(isTRUE(methods::validObject(scop_assay, test = TRUE)))

  reordered_obj <- suppressWarnings(suppressMessages(
    Seurat::CreateSeuratObject(
      counts = list(A = counts_a, B = counts_b[rev(seq_len(nrow(counts_b))), ])
    )
  ))
  reordered_scop <- suppressWarnings(FindVariableFeatures(
    reordered_obj,
    selection.method = "vst",
    nfeatures = 20,
    verbose = FALSE
  ))
  reordered_seurat <- suppressWarnings(seurat_reference_method(
    "FindVariableFeatures",
    "Seurat",
    reordered_obj,
    selection.method = "vst",
    nfeatures = 20,
    verbose = FALSE
  ))
  expect_identical(
    SeuratObject::VariableFeatures(reordered_scop),
    SeuratObject::VariableFeatures(reordered_seurat)
  )
})
