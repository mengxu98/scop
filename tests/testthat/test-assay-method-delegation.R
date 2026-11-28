test_that("RunGLMPCA Assay5 delegates with an unchanged public signature", {
  expect_identical(
    formals(RunGLMPCA.Assay5),
    formals(RunGLMPCA.Assay)
  )

  received <- NULL
  testthat::local_mocked_bindings(
    RunGLMPCA.Assay = function(...) {
      received <<- list(...)
      "glmpca-result"
    }
  )
  object <- structure(list(), class = "Assay5")
  result <- RunGLMPCA.Assay5(
    object,
    assay = "RNA",
    layer = "counts",
    features = c("g1", "g2"),
    L = 3,
    fam = "nb",
    rev.gmlpca = TRUE,
    ndims.print = 1:2,
    nfeatures.print = 4,
    reduction.key = "TEST_",
    verbose = FALSE,
    seed.use = 17,
    custom = "value"
  )

  expect_identical(result, "glmpca-result")
  expect_identical(received$object, object)
  expect_identical(received$features, c("g1", "g2"))
  expect_identical(received$seed.use, 17)
  expect_identical(received$custom, "value")
})

test_that("RunMDS Assay5 delegates with an unchanged public signature", {
  expect_identical(
    formals(RunMDS.Assay5),
    formals(RunMDS.Assay)
  )

  received <- NULL
  testthat::local_mocked_bindings(
    RunMDS.Assay = function(...) {
      received <<- list(...)
      "mds-result"
    }
  )
  object <- structure(list(), class = "Assay5")
  result <- RunMDS.Assay5(
    object,
    assay = "RNA",
    layer = "data",
    features = c("g1", "g2"),
    nmds = 3,
    dist.method = "manhattan",
    mds.method = "isoMDS",
    rev.mds = TRUE,
    reduction.key = "TEST_",
    verbose = FALSE,
    seed.use = 19,
    custom = "value"
  )

  expect_identical(result, "mds-result")
  expect_identical(received$object, object)
  expect_identical(received$features, c("g1", "g2"))
  expect_identical(received$seed.use, 19)
  expect_identical(received$custom, "value")
})

test_that("Assay5 delegation matches real GLMPCA and MDS results", {
  skip_if_not_installed("glmpca")

  counts <- Matrix::Matrix(
    matrix(
      c(
        1, 0, 3, 1, 4, 2, 0, 5,
        0, 2, 1, 3, 0, 4, 2, 1,
        5, 3, 0, 2, 1, 0, 4, 3,
        2, 1, 4, 0, 3, 2, 5, 0,
        0, 4, 2, 5, 1, 3, 0, 2,
        3, 0, 5, 2, 4, 1, 3, 1
      ),
      nrow = 6,
      byrow = TRUE,
      dimnames = list(
        paste0("gene", 1:6),
        paste0("cell", 1:8)
      )
    ),
    sparse = TRUE
  )
  assay5 <- SeuratObject::CreateAssay5Object(counts = counts)
  features <- rownames(counts)

  set.seed(23)
  glmpca_direct <- RunGLMPCA.Assay(
    assay5,
    assay = "RNA",
    features = features,
    L = 2,
    verbose = FALSE
  )
  set.seed(23)
  glmpca_delegated <- RunGLMPCA.Assay5(
    assay5,
    assay = "RNA",
    features = features,
    L = 2,
    verbose = FALSE
  )
  expect_equal(
    SeuratObject::Embeddings(glmpca_delegated),
    SeuratObject::Embeddings(glmpca_direct),
    tolerance = 1e-10
  )
  expect_equal(
    SeuratObject::Loadings(glmpca_delegated),
    SeuratObject::Loadings(glmpca_direct),
    tolerance = 1e-10
  )

  mds_direct <- RunMDS.Assay(
    assay5,
    assay = "RNA",
    layer = "counts",
    features = features,
    nmds = 2,
    seed.use = 29
  )
  mds_delegated <- RunMDS.Assay5(
    assay5,
    assay = "RNA",
    layer = "counts",
    features = features,
    nmds = 2,
    seed.use = 29
  )
  expect_equal(
    SeuratObject::Embeddings(mds_delegated),
    SeuratObject::Embeddings(mds_direct),
    tolerance = 1e-10
  )
})

test_that("RenameFeatures supports full unnamed vectors in a standard v3 assay", {
  counts <- Matrix::Matrix(
    matrix(0:8, nrow = 3L, dimnames = list(c("g1", "g2", "g3"), c("c1", "c2", "c3"))),
    sparse = TRUE
  )
  assay <- SeuratObject::CreateAssayObject(counts = counts)
  assay <- SeuratObject::AddMetaData(assay, data.frame(score = c(3, 1, 2), row.names = rownames(counts)))
  SeuratObject::VariableFeatures(assay) <- c("g3", "g1")
  assay@scale.data <- matrix(
    1:6, nrow = 2L, dimnames = list(c("g2", "g3"), colnames(counts))
  )
  object <- SeuratObject::CreateSeuratObject(assay)
  renamed <- RenameFeatures(object, c("z", "y", "x"), verbose = FALSE)
  result <- renamed[["RNA"]]

  expect_identical(rownames(renamed), c("z", "y", "x"))
  for (layer in c("counts", "data", "scale.data")) {
    old <- SeuratObject::LayerData(object[["RNA"]], layer = layer)
    new <- SeuratObject::LayerData(result, layer = layer)
    rownames(old) <- if (layer == "scale.data") c("y", "x") else c("z", "y", "x")
    expect_identical(new, old)
  }
  expect_identical(result[[]][["score"]], object[["RNA"]][[]][["score"]])
  expect_identical(SeuratObject::VariableFeatures(result), c("x", "z"))
  expect_identical(SeuratObject::Key(result), SeuratObject::Key(object[["RNA"]]))
  expect_identical(renamed[[]], object[[]])
  expect_true(methods::validObject(renamed, test = TRUE))
})

test_that("RenameFeatures preserves partially overlapping v5 layers and assay identity", {
  sparse <- function(values, features, cells) {
    Matrix::Matrix(matrix(values, nrow = length(features), dimnames = list(features, cells)),
                   sparse = TRUE)
  }
  layer_a <- sparse(c(0, 4, 9, 3), c("g3", "g1"), c("c1", "c2"))
  layer_b <- sparse(c(1, 12, 5, 1), c("g1", "g2"), c("c3", "c4"))
  assay <- SeuratObject::CreateAssay5Object(counts = list(A = layer_a, B = layer_b))
  SeuratObject::VariableFeatures(assay) <- c("g2", "g1")
  assay[["type"]] <- setNames(factor(c("III", "I", "II")), rownames(assay))
  object <- SeuratObject::CreateSeuratObject(assay)
  object[["other"]] <- SeuratObject::CreateAssay5Object(
    counts = sparse(seq_len(8L), c("different", "g2"), colnames(object))
  )
  result <- RenameFeatures(object, c(g1 = "alpha", g2 = "beta"), assays = "RNA", verbose = FALSE)
  renamed <- result[["RNA"]]

  expect_identical(rownames(renamed), c("g3", "alpha", "beta"))
  expect_identical(SeuratObject::Layers(renamed), c("counts.A", "counts.B"))
  rownames(layer_a) <- c("g3", "alpha")
  rownames(layer_b) <- c("alpha", "beta")
  expect_identical(SeuratObject::LayerData(renamed, layer = "counts.A"), layer_a)
  expect_identical(SeuratObject::LayerData(renamed, layer = "counts.B"), layer_b)
  expect_identical(SeuratObject::VariableFeatures(renamed), c("beta", "alpha"))
  expect_identical(SeuratObject::Key(renamed), SeuratObject::Key(object[["RNA"]]))
  expect_identical(renamed[[]][["var.features.rank"]], assay[[]][["var.features.rank"]])
  expect_identical(renamed[[]][["type"]], assay[[]][["type"]])
  expect_identical(rownames(renamed[[]]), rownames(renamed))
  expect_identical(result[["other"]], object[["other"]])
  expect_identical(result[[]], object[[]])
  expect_true(methods::validObject(result, test = TRUE))

  renamed_other <- RenameFeatures(object, c("first", "second"), assays = "other", verbose = FALSE)
  expect_identical(renamed_other[["RNA"]], object[["RNA"]])
  expect_identical(rownames(renamed_other[["other"]]), c("first", "second"))
  expect_error(RenameFeatures(object, letters[1:3], verbose = FALSE), "different feature IDs")
  expect_error(RenameFeatures(object, c(g1 = "g2"), verbose = FALSE), "must be unique")
  expect_error(RenameFeatures(object, c(g1 = "a", g1 = "b"), verbose = FALSE), "must be non-empty and unique")
})

test_that("RenameFeatures refuses stale feature-dependent reductions", {
  counts <- Matrix::Matrix(
    matrix(1:6, nrow = 2L, dimnames = list(c("g1", "g2"), c("c1", "c2", "c3"))),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(SeuratObject::CreateAssay5Object(counts = counts))
  embeddings <- matrix(1:3, ncol = 1L, dimnames = list(colnames(object), "PC_1"))
  loadings <- matrix(c(0.2, 0.8), ncol = 1L, dimnames = list(rownames(object), "PC_1"))
  object[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings, loadings = loadings, key = "PC_", assay = "RNA"
  )
  expect_error(RenameFeatures(object, c(g1 = "new"), verbose = FALSE), "contains feature loadings")
  expect_identical(RenameFeatures(object, c(absent = "new"), verbose = FALSE), object)
})
