test_that("RunKNNMap reuses checked data layers for feature projection", {
  counts <- Matrix::Matrix(
    matrix(
      c(1, 3, 2, 4, 2, 1, 5, 3, 4, 2, 1, 6),
      nrow = 4,
      dimnames = list(paste0("Gene", 1:4), paste0("Cell", 1:3))
    ),
    sparse = TRUE
  )
  query <- Seurat::CreateSeuratObject(counts)
  reference <- Seurat::CreateSeuratObject(counts)
  SeuratObject::VariableFeatures(query) <- rownames(query)
  SeuratObject::VariableFeatures(reference) <- rownames(reference)
  umap <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      seq_len(6), ncol = 2,
      dimnames = list(colnames(reference), c("UMAP_1", "UMAP_2"))
    ),
    assay = "RNA",
    key = "UMAP_"
  )
  umap@misc$model <- list(
    layout = matrix(0, nrow = 3, ncol = 2),
    config = list(n_neighbors = 2)
  )
  reference[["umap"]] <- umap

  query_data <- matrix(
    seq_len(12), nrow = 4,
    dimnames = list(rownames(query), colnames(query))
  )
  ref_data <- matrix(
    seq_len(12) + 20, nrow = 4,
    dimnames = list(rownames(reference), colnames(reference))
  )
  calls <- character()
  projection_input <- NULL

  testthat::local_mocked_bindings(
    GetAssayData5 = function(object, layer, assay, ...) {
      calls <<- c(calls, if (identical(object, query)) "query" else "reference")
      if (identical(object, query)) query_data else ref_data
    },
    CheckDataType = function(...) "lognormalized",
    RunUMAP2 = function(object, ...) {
      projection_input <<- object
      SeuratObject::CreateDimReducObject(
        embeddings = matrix(
          0,
          nrow(object),
          2,
          dimnames = list(rownames(object), c("REF_1", "REF_2"))
        ),
        assay = "RNA",
        key = "REF_"
      )
    },
    .package = "scop"
  )

  out <- RunKNNMap(
    srt_query = query,
    srt_ref = reference,
    ref_umap = "umap",
    features = rownames(query),
    projection_method = "model",
    verbose = FALSE
  )

  expect_identical(calls, c("query", "reference"))
  expect_equal(projection_input, t(query_data), tolerance = 1e-12)
  expect_true("ref.embeddings" %in% SeuratObject::Reductions(out))
})

test_that("RunKNNMap fills missing neighbor distances like the legacy loop", {
  distances <- matrix(
    c(NA_real_, NA_real_, 0.2, NA_real_, Inf, NA_real_, 0.5, 0.1, 0.3),
    nrow = 3,
    byrow = TRUE
  )
  legacy <- distances
  global_max <- max(legacy, na.rm = TRUE)
  for (i in seq_len(nrow(legacy))) {
    missing <- is.na(legacy[i, ])
    if (any(missing)) {
      legacy[i, missing] <- if (all(missing)) {
        global_max
      } else {
        max(legacy[i, !missing], na.rm = TRUE)
      }
    }
  }
  testthat::local_mocked_bindings(check_r = function(...) TRUE, .package = "scop")

  expect_identical(knn_fill_missing_distances(distances), legacy)
})
