make_spatialqm_seurat <- function() {
  counts <- matrix(
    c(
      5, 0, 2, 1,
      0, 3, 0, 4,
      1, 1, 1, 1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:3), paste0("Cell", 1:4))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts, assay = "Spatial")
  srt$sample_id <- "sampleA"
  srt
}

with_mock_spatialqm <- function(code, fail_entropy = FALSE) {
  fake_metric <- function(name) {
    switch(name,
      getNcells = function(seu_obj, ...) {
        expect_true("RNA" %in% SeuratObject::Assays(seu_obj))
        ncol(seu_obj)
      },
      getTxPerCell = function(seu_obj, features = NULL, sample.p = NULL, ...) {
        expect_equal(features, c("Gene1", "Gene2"))
        expect_equal(sample.p, 0.25)
        mat <- GetAssayData5(seu_obj, assay = "RNA", layer = "counts")
        mean(Matrix::colSums(mat[features, , drop = FALSE]))
      },
      getSparsity = function(seu_obj, features = NULL, ...) 0.5,
      getEntropy = function(seu_obj, features = NULL, ...) {
        if (isTRUE(fail_entropy)) {
          stop("entropy failed", call. = FALSE)
        }
        1.25
      },
      stop("unexpected SpatialQM metric: ", name)
    )
  }
  testthat::local_mocked_bindings(
    check_r = function(packages, ...) {
      expect_identical(packages, "SpatialQM")
      invisible(TRUE)
    },
    get_namespace_fun = function(package, name) {
      expect_identical(package, "SpatialQM")
      fake_metric(name)
    }
  )
  force(code)
}

test_that("RunSpatialQM stores metric results and summary", {
  srt <- make_spatialqm_seurat()
  with_mock_spatialqm({
    out <- RunSpatialQM(
      srt,
      assay = "Spatial",
      layer = "counts",
      metrics = c("n_cells", "tx_per_cell", "sparsity", "entropy"),
      features = c("Gene1", "Gene2"),
      platform = "Xenium",
      sample.p = 0.25,
      verbose = FALSE
    )
  })

  expect_s4_class(out, "Seurat")
  expect_false("RNA" %in% SeuratObject::Assays(out))
  expect_true("SpatialQM" %in% names(out@tools))
  bundle <- out@tools$SpatialQM
  expect_named(bundle$results, c("n_cells", "tx_per_cell", "sparsity", "entropy"))
  expect_equal(bundle$results$n_cells, 4)
  expect_equal(bundle$results$sparsity, 0.5)
  expect_named(bundle$summary, c("metric", "status", "class", "rows", "columns", "value", "message"))
  expect_equal(bundle$summary$status, rep("ok", 4))
  expect_equal(bundle$parameters$assay, "Spatial")
  expect_equal(bundle$parameters$features, c("Gene1", "Gene2"))
  expect_equal(bundle$parameters$platform, "Xenium")
  expect_equal(bundle$parameters$backend_args$sample.p, 0.25)
})

test_that("RunSpatialQM can record failed metrics without stopping", {
  srt <- make_spatialqm_seurat()
  with_mock_spatialqm(fail_entropy = TRUE, {
    out <- RunSpatialQM(
      srt,
      assay = "Spatial",
      metrics = c("n_cells", "entropy"),
      on_error = "warning",
      verbose = FALSE
    )
  })

  expect_true("SpatialQM" %in% names(out@tools))
  expect_equal(names(out@tools$SpatialQM$results), "n_cells")
  expect_equal(out@tools$SpatialQM$summary$status, c("ok", "error"))
  expect_equal(out@tools$SpatialQM$errors$entropy, "entropy failed")
})

test_that("RunSpatialQM validates inputs clearly", {
  srt <- make_spatialqm_seurat()
  expect_error(
    RunSpatialQM(srt, assay = "missing", verbose = FALSE),
    "not present"
  )
  expect_error(
    RunSpatialQM(srt, metrics = "not_a_metric", verbose = FALSE),
    "Unknown SpatialQM metric"
  )
  expect_error(
    RunSpatialQM(srt, features = "missing_gene", verbose = FALSE),
    "features"
  )
})

test_that("RunSpatialQM can skip tool storage", {
  srt <- make_spatialqm_seurat()
  with_mock_spatialqm({
    out <- RunSpatialQM(
      srt,
      assay = "Spatial",
      metrics = "n_cells",
      store_results = FALSE,
      verbose = FALSE
    )
  })

  expect_false("SpatialQM" %in% names(out@tools))
})

test_that("RunSpatialQM accepts SpatialQM function-name metric aliases", {
  srt <- make_spatialqm_seurat()
  with_mock_spatialqm({
    out <- RunSpatialQM(
      srt,
      assay = "Spatial",
      metrics = c("getNcells", "n_cells"),
      verbose = FALSE
    )
  })

  expect_equal(names(out@tools$SpatialQM$results), "n_cells")
  expect_equal(out@tools$SpatialQM$results$n_cells, 4)
})
