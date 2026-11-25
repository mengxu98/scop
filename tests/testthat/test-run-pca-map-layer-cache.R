make_pcamap_cache_objects <- function() {
  counts <- Matrix::Matrix(
    matrix(
      c(4, 1, 2, 3, 1, 5, 2, 4, 3, 6, 5, 2),
      nrow = 4,
      dimnames = list(paste0("Gene", 1:4), paste0("Cell", 1:3))
    ),
    sparse = TRUE
  )
  query <- Seurat::CreateSeuratObject(counts)
  reference <- Seurat::CreateSeuratObject(counts)
  loadings <- matrix(
    c(0.4, -0.2, 0.1, 0.3, 0.2, 0.5, -0.4, 0.1),
    nrow = 4,
    dimnames = list(rownames(reference), c("PC_1", "PC_2"))
  )
  reference[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      seq_len(6),
      ncol = 2,
      dimnames = list(colnames(reference), c("PC_1", "PC_2"))
    ),
    loadings = loadings,
    assay = "RNA",
    key = "PC_"
  )
  reference[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      seq_len(6),
      ncol = 2,
      dimnames = list(colnames(reference), c("UMAP_1", "UMAP_2"))
    ),
    assay = "RNA",
    key = "UMAP_"
  )
  list(query = query, reference = reference, loadings = loadings)
}

test_that("RunPCAMap reuses checked data layers for PCA projection", {
  objects <- make_pcamap_cache_objects()
  query_data <- matrix(
    c(1, 2, 3, 2, 4, 5, 3, 5, 7, 4, 6, 8),
    nrow = 4,
    dimnames = list(paste0("Gene", 1:4), colnames(objects$query))
  )
  ref_data <- matrix(
    c(2, 3, 5, 4, 6, 7, 8, 5, 4, 9, 3, 6),
    nrow = 4,
    dimnames = list(paste0("Gene", 1:4), colnames(objects$reference))
  )
  calls <- character()

  testthat::local_mocked_bindings(
    GetAssayData5 = function(object, layer, assay, ...) {
      calls <<- c(calls, paste(if (identical(object, objects$query)) "query" else "reference", layer, sep = ":"))
      if (identical(object, objects$query)) query_data else ref_data
    },
    CheckDataType = function(object, ...) {
      force(object)
      "lognormalized"
    },
    RunKNNMap = function(srt_query, ...) srt_query,
    .package = "scop"
  )

  out <- RunPCAMap(
    objects$query,
    objects$reference,
    ref_pca = "pca",
    ref_umap = "umap",
    ref_dims = 1:2,
    projection_method = "knn",
    verbose = FALSE
  )
  expected <- scale(
    t(query_data),
    center = rowMeans(ref_data),
    scale = apply(ref_data, 1, stats::sd)
  ) %*% objects$loadings

  expect_identical(calls, c("query:data", "reference:data"))
  expect_equal(SeuratObject::Embeddings(out[["ref.pca"]]), expected, tolerance = 1e-12)
})

make_pcamap_single_feature_objects <- function() {
  ref_counts <- matrix(
    c(2, 8, 4, 9, 1, 6, 3, 14, 12, 2, 7, 6, 5, 7, 3, 4),
    nrow = 4L, byrow = TRUE,
    dimnames = list(paste0("g", 1:4), c("r3", "r1", "r4", "r2"))
  )
  query_counts <- matrix(
    c(4, 8, 2, 2, 8, 16, 3, 2, 5, 1, 3, 7),
    nrow = 4L, byrow = TRUE,
    dimnames = list(c("qx", "g2", "g1", "qy"), c("q3", "q1", "q2"))
  )
  ref_data <- log1p(ref_counts)
  query_data <- log1p(query_counts[c("qx", "g2", "qy"), , drop = FALSE])
  fit <- stats::prcomp(t(ref_data), center = TRUE, scale. = TRUE, rank. = 2L)
  colnames(fit$x) <- colnames(fit$rotation) <- c("PC_1", "PC_2")

  make_object <- function(counts, data) {
    SeuratObject::CreateSeuratObject(SeuratObject::CreateAssay5Object(
      counts = Matrix::Matrix(counts, sparse = TRUE),
      data = Matrix::Matrix(data, sparse = TRUE)
    ))
  }
  reference <- make_object(ref_counts, ref_data)
  query <- make_object(query_counts, query_data)
  reference[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = fit$x, loadings = fit$rotation, assay = "RNA", key = "PC_"
  )
  umap <- matrix(
    c(1, 3, 5, 2, 4, 7, 6, 8), ncol = 2L,
    dimnames = list(colnames(reference), c("UMAP_1", "UMAP_2"))
  )
  reference[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = umap, assay = "RNA", key = "UMAP_"
  )
  list(query = query, reference = reference, fit = fit, query_data = query_data,
       reference_data = ref_data, umap = umap)
}

test_that("RunPCAMap projects one eligible shared feature without dropping axes", {
  objects <- make_pcamap_single_feature_objects()
  result <- RunPCAMap(
    objects$query, objects$reference, ref_pca = "pca", ref_umap = "umap",
    ref_dims = 1L, projection_method = "knn", nn_method = "raw", k = 1L,
    distance_metric = "euclidean", verbose = FALSE
  )

  q <- objects$query_data["g2", ]
  r <- objects$reference_data["g2", ]
  expected_pc <- outer((q - mean(r)) / stats::sd(r), objects$fit$rotation["g2", ])
  expect_identical(colnames(result), names(q))
  expect_equal(SeuratObject::Embeddings(result[["ref.pca"]]), expected_pc, tolerance = 1e-12)
  nearest <- apply(abs(outer(objects$fit$x[, 1L], expected_pc[, 1L], "-")), 2L, which.min)
  expected_umap <- objects$umap[nearest, , drop = FALSE]
  observed_umap <- SeuratObject::Embeddings(result[["ref.embeddings"]])[names(q), , drop = FALSE]
  expect_equal(unname(observed_umap), unname(expected_umap), tolerance = 0)
  expect_true(methods::validObject(result, test = TRUE))
})

test_that("RunPCAMap explains missing overlap or unusable sole reference variance", {
  objects <- make_pcamap_single_feature_objects()
  project <- function(query, reference) {
    RunPCAMap(
      query, reference, ref_pca = "pca", ref_umap = "umap", ref_dims = 1L,
      projection_method = "knn", nn_method = "raw", k = 1L,
      distance_metric = "euclidean", verbose = FALSE
    )
  }
  query <- subset(objects$query, features = c("qx", "qy"))
  expect_error(project(query, objects$reference), "No shared PCA features")

  reference <- objects$reference
  layer <- SeuratObject::LayerData(reference, assay = "RNA", layer = "data")
  layer["g2", ] <- layer["g2", 1L]
  SeuratObject::LayerData(reference, assay = "RNA", layer = "data") <- layer
  expect_error(project(objects$query, reference), "sole shared PCA feature")
})
