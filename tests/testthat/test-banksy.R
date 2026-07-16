make_banksy_seurat <- function() {
  counts <- matrix(
    c(
      10, 8, 1, 0,
      0, 2, 9, 8,
      6, 0, 1, 0,
      1, 7, 2, 5
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:4), paste0("Spot", 1:4))
  )
  srt <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  )
  srt$col <- c(1, 2, 1, 2)
  srt$row <- c(1, 1, 2, 2)
  srt$sample <- c("S1", "S1", "S2", "S2")
  srt
}

with_mock_banksy <- function(code) {
  compute_fun <- function(se, assay_name, coord_names, compute_agf, M, k_geom, ...) {
    expect_s4_class(se, "SpatialExperiment")
    expect_equal(assay_name, "scop_input")
    expect_equal(coord_names, c("x", "y"))
    expect_false(compute_agf)
    expect_equal(M, 1)
    expect_equal(k_geom, 15)
    se
  }
  pca_fun <- function(se, assay_name, M, lambda, npcs, use_agf, group, seed, ...) {
    expect_equal(lambda, 0.2)
    expect_equal(npcs, 20)
    expect_false(use_agf)
    expect_equal(group, "sample")
    expect_equal(seed, 1)
    se
  }
  cluster_fun <- function(se, assay_name, M, lambda, use_agf, npcs, algo, k_neighbors, resolution, group, seed, ...) {
    expect_equal(npcs, 20)
    expect_equal(algo, "leiden")
    expect_equal(k_neighbors, 50)
    expect_equal(resolution, 0.6)
    cdata <- as.data.frame(SummarizedExperiment::colData(se))
    cdata$BANKSY_leiden <- c("1", "1", "2", "2")
    SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(cdata)
    se
  }
  cluster_names <- function(se) "BANKSY_leiden"
  testthat::local_mocked_bindings(
    check_r = function(packages, ...) {
      expect_true("Banksy" %in% packages)
      invisible(TRUE)
    },
    get_namespace_fun = function(package, name) {
      if (identical(package, "SpatialExperiment")) {
        return(getExportedValue(package, name))
      }
      expect_identical(package, "Banksy")
      switch(name,
        computeBanksy = compute_fun,
        runBanksyPCA = pca_fun,
        clusterBanksy = cluster_fun,
        clusterNames = cluster_names,
        stop("unexpected function")
      )
    }
  )
  force(code)
}

test_that("RunBANKSY writes cluster metadata and tool results", {
  testthat::skip_if_not_installed("SpatialExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("S4Vectors")
  srt <- make_banksy_seurat()
  with_mock_banksy({
    out <- RunBANKSY(srt, layer = "counts", group = "sample", verbose = FALSE)
  })

  expect_equal(unname(out$BANKSY_cluster), c("1", "1", "2", "2"))
  expect_true("BANKSY" %in% names(out@tools))
  expect_equal(out@tools$BANKSY$cluster_source, "BANKSY_leiden")
  expect_equal(out@tools$BANKSY$parameters$group, "sample")
  expect_identical(out@tools$BANKSY$source$coordinate_space, "raw")
  expect_named(out@tools$BANKSY$summary, c("n_spots", "domains"))
})

test_that("RunBANKSY validates inputs before backend work", {
  srt <- make_banksy_seurat()
  expect_error(
    RunBANKSY(matrix(1, nrow = 2), verbose = FALSE),
    "Seurat"
  )
  with_mock_banksy({
    expect_error(
      RunBANKSY(srt, layer = "counts", features = "AbsentGene", verbose = FALSE),
      "No features"
    )
    expect_error(
      RunBANKSY(srt, layer = "counts", group = "missing", verbose = FALSE),
      "group"
    )
    expect_error(
      RunBANKSY(srt, layer = "counts", run_pca_params = list(1), verbose = FALSE),
      "named arguments"
    )
  })
})

test_that("BANKSY clusters reuse SCOP SpatialSpotPlot", {
  srt <- make_banksy_seurat()
  srt$BANKSY_cluster <- c("1", "1", "2", "2")
  p <- SpatialSpotPlot(
    srt,
    group.by = "BANKSY_cluster",
    overlay_image = FALSE
  )
  expect_s3_class(p, "ggplot")
})
