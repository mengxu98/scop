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
