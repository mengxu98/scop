make_standard_spatial_atac_object <- function() {
  skip_if_not_installed("Signac")
  counts <- matrix(
    c(
      2, 0, 1, 0,
      0, 3, 0, 1,
      1, 1, 0, 2,
      0, 1, 2, 1,
      3, 0, 1, 1,
      1, 2, 0, 1
    ),
    nrow = 6,
    byrow = TRUE,
    dimnames = list(
      paste0("chr1-", seq(100, 600, by = 100), "-", seq(149, 649, by = 100)),
      paste0("spot", 1:4)
    )
  )
  chromatin <- Signac::CreateChromatinAssay(
    counts = Matrix::Matrix(counts, sparse = TRUE),
    sep = c("-", "-")
  )
  Seurat::CreateSeuratObject(counts = chromatin, assay = "peaks")
}

add_standard_spatial_atac_reduction <- function(srt, reduction) {
  embeddings <- matrix(
    seq_len(ncol(srt) * 3L),
    nrow = ncol(srt),
    dimnames = list(colnames(srt), paste0("DR_", 1:3))
  )
  srt[[reduction]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings,
    assay = "peaks",
    key = "DR_"
  )
  srt
}

test_that("ATAC defaults are applied before spatial output planning", {
  srt <- make_standard_spatial_atac_object()
  captured <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) {
      captured$nested <- list(...)
      add_standard_spatial_atac_reduction(srt, "ATACsvd")
    },
    RunBayesSpace = function(
      srt,
      cluster_colname,
      init_colname,
      use_reduction = NULL,
      ...
    ) {
      captured$bayes <- list(
        cluster_colname = cluster_colname,
        init_colname = init_colname,
        use_reduction = use_reduction
      )
      srt[[cluster_colname]] <- rep("domain1", ncol(srt))
      srt@tools[["BayesSpace"]] <- list(result = "fresh")
      srt
    },
    .package = "scop"
  )
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")

  out <- original(
    srt,
    assay = "peaks",
    do_spot_qc = FALSE,
    do_spatial_variable_features = FALSE,
    do_spatial_cluster = TRUE,
    spatial_q = 2,
    bayesspace_params = list(
      cluster_colname = "Standardclusters",
      init_colname = NULL
    ),
    do_deconvolution = FALSE,
    verbose = FALSE
  )
  stages <- out@tools$run_standard_spatial_workflow$stages
  domain <- stages[stages$stage == "spatial_clustering", , drop = FALSE]

  expect_identical(captured$nested$prefix, "ATAC")
  expect_true(captured$nested$do_normalization)
  expect_identical(captured$nested$normalization_method, "TFIDF")
  expect_true(captured$nested$do_HVF_finding)
  expect_identical(captured$nested$nHVF, 20000)
  expect_false(captured$nested$do_scaling)
  expect_identical(captured$nested$linear_reduction, "svd")
  expect_identical(captured$nested$linear_reduction_dims_use, 2:30)
  expect_identical(captured$nested$neighbor_metric, "cosine")
  expect_identical(captured$bayes$cluster_colname, "Standardclusters")
  expect_null(captured$bayes$init_colname)
  expect_identical(captured$bayes$use_reduction, "ATACsvd")
  expect_identical(domain$status, "completed")
  expect_identical(domain$result_metadata_key, "Standardclusters")
  expect_identical(
    out@tools$run_standard_spatial_workflow$parameters$prefix,
    "ATAC"
  )
})

test_that("ATAC preprocessing cluster targets reject BayesSpace collisions", {
  srt <- make_standard_spatial_atac_object()
  caller_before <- srt
  producer_calls <- character()
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) {
      producer_calls <<- c(producer_calls, "RunStandardWorkflow")
      srt
    },
    RunBayesSpace = function(srt, ...) {
      producer_calls <<- c(producer_calls, "RunBayesSpace")
      srt
    },
    .package = "scop"
  )
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")

  for (target in c("ATACclusters", "ATACsvdclusters", "ATAClsiclusters")) {
    for (bayes_arg in c("cluster_colname", "init_colname")) {
      bayes_params <- list(
        cluster_colname = "CustomCluster",
        init_colname = "CustomInit"
      )
      bayes_params[[bayes_arg]] <- target
      expect_error(
        original(
          srt,
          assay = "peaks",
          do_spot_qc = FALSE,
          do_spatial_variable_features = FALSE,
          do_spatial_cluster = TRUE,
          spatial_q = 2,
          bayesspace_params = bayes_params,
          do_deconvolution = FALSE,
          verbose = FALSE
        ),
        "metadata outputs collide"
      )
      expect_length(producer_calls, 0L)
      expect_identical(srt, caller_before)
    }
  }
})

test_that("custom ATAC workflow values retain atac_defaults semantics", {
  srt <- make_standard_spatial_atac_object()
  captured <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) {
      captured$nested <- list(...)
      add_standard_spatial_atac_reduction(srt, "Customsvd")
    },
    RunBayesSpace = function(
      srt,
      cluster_colname,
      init_colname,
      use_reduction = NULL,
      ...
    ) {
      captured$bayes_use_reduction <- use_reduction
      srt[[cluster_colname]] <- rep("domain1", ncol(srt))
      srt@tools[["BayesSpace"]] <- list(result = "fresh")
      srt
    },
    .package = "scop"
  )
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")

  out <- original(
    srt,
    prefix = "Custom",
    assay = "peaks",
    do_spot_qc = FALSE,
    do_spatial_variable_features = FALSE,
    do_spatial_cluster = TRUE,
    spatial_q = 2,
    bayesspace_params = list(
      cluster_colname = "Customicaclusters",
      init_colname = NULL
    ),
    do_deconvolution = FALSE,
    do_normalization = FALSE,
    normalization_method = "TFIDF",
    do_HVF_finding = FALSE,
    nHVF = 321,
    do_scaling = FALSE,
    linear_reduction = "ica",
    linear_reduction_dims_use = 1:2,
    neighbor_metric = "manhattan",
    verbose = FALSE
  )
  stages <- out@tools$run_standard_spatial_workflow$stages
  domain <- stages[stages$stage == "spatial_clustering", , drop = FALSE]

  expect_identical(captured$nested$prefix, "Custom")
  expect_false(captured$nested$do_normalization)
  expect_identical(captured$nested$normalization_method, "TFIDF")
  expect_false(captured$nested$do_HVF_finding)
  expect_identical(captured$nested$nHVF, 321)
  expect_false(captured$nested$do_scaling)
  expect_identical(captured$nested$linear_reduction, "ica")
  expect_identical(captured$nested$linear_reduction_dims_use, 1:2)
  expect_identical(captured$nested$neighbor_metric, "manhattan")
  expect_identical(captured$bayes_use_reduction, "Customsvd")
  expect_identical(domain$status, "completed")
  expect_identical(domain$result_metadata_key, "Customicaclusters")
})
