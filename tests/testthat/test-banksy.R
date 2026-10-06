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

with_mock_banksy <- function(code, fail_sample = NULL, expected_group = "sample") {
  compute_fun <- function(se, assay_name, coord_names, compute_agf, M, k_geom, ...) {
    expect_s4_class(se, "SpatialExperiment")
    expect_equal(assay_name, "scop_input")
    expect_equal(coord_names, c("x", "y"))
    expect_false(compute_agf)
    expect_equal(M, 1)
    expect_equal(k_geom, 15)
    se
  }
  pca_fun <- function(se, assay_name, M, lambda, npcs, use_agf, group = NULL, seed, ...) {
    expect_equal(lambda, 0.2)
    expect_equal(npcs, 20)
    expect_false(use_agf)
    expect_equal(group, expected_group)
    expect_equal(seed, 1)
    se
  }
  cluster_fun <- function(se, assay_name, M, lambda, use_agf, npcs, algo, k_neighbors, resolution, group, seed, ...) {
    expect_equal(npcs, 20)
    expect_equal(algo, "leiden")
    expect_equal(k_neighbors, 50)
    expect_equal(resolution, 0.6)
    cdata <- as.data.frame(SummarizedExperiment::colData(se))
    if (!is.null(fail_sample) && "sample" %in% colnames(cdata) &&
      any(as.character(cdata$sample) == fail_sample)) {
      stop("mock backend failure", call. = FALSE)
    }
    cdata$BANKSY_leiden <- rep(c("1", "1", "2", "2"), length.out = nrow(cdata))
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
  expect_identical(out@tools$BANKSY$parameters$coordinate_space, "raw")
  expect_false("per_sample" %in% names(out@tools$BANKSY))
  expect_named(out@tools$BANKSY$summary, c("n_spots", "domains"))

  result <- GetSpatialResult(out, "BANKSY")
  expect_named(result, c("clusters", "parameters", "summary"))
  expect_equal(result$clusters, out@tools$BANKSY$clusters)
  expect_equal(result$parameters, out@tools$BANKSY$parameters)
  expect_equal(result$summary, out@tools$BANKSY$summary)
})

test_that("BANKSY chooses a spatial assay before RNA and respects explicit assays", {
  srt <- make_banksy_seurat()
  suppressWarnings(srt[["Spatial"]] <- srt[["RNA"]])
  SeuratObject::DefaultAssay(srt) <- "RNA"

  expect_identical(banksy_resolve_assay(srt)$assay, "Spatial")
  expect_identical(banksy_resolve_assay(srt, assay = "RNA")$assay, "RNA")
  expect_identical(banksy_resolve_assay(srt)$source, "Spatial assay")
})

test_that("BANKSY uses an image-associated assay and reports a missing layer", {
  data(visium_human_pancreas_sub)
  srt <- visium_human_pancreas_sub

  expect_identical(
    banksy_resolve_assay(srt, image = "slice1")$assay,
    "Spatial"
  )
  expect_error(
    RunBANKSY(make_banksy_seurat(), verbose = FALSE),
    "Layer .* is not present"
  )
})

test_that("RunBANKSY auto-selects its single image and associated assay", {
  data(visium_human_pancreas_sub)
  srt <- suppressWarnings(visium_human_pancreas_sub[, seq_len(4)])
  expected_coords <- suppressWarnings(attr(
    resolve_spatial_spot_coords(srt, colnames(srt), image = "slice1"),
    "spatial_source",
    exact = TRUE
  )$coord.cols)
  suppressWarnings(with_mock_banksy({
    out <- RunBANKSY(
      srt,
      layer = "counts",
      features = rownames(srt)[seq_len(5)],
      verbose = FALSE
    )
  }, expected_group = NULL))

  expect_identical(out@tools$BANKSY$parameters$assay, "Spatial")
  expect_identical(out@tools$BANKSY$parameters$image, "slice1")
  expect_identical(out@tools$BANKSY$parameters$coord.cols, expected_coords)
})

test_that("BANKSY keeps metadata clusters accessible without detailed storage", {
  srt <- make_banksy_seurat()
  suppressWarnings(with_mock_banksy({
    out <- RunBANKSY(
      srt,
      layer = "counts",
      group = "sample",
      store_results = FALSE,
      verbose = FALSE
    )
  }))

  expect_false("BANKSY" %in% names(out@tools))
  result <- GetSpatialResult(out, "BANKSY")
  expect_null(result$parameters)
  expect_equal(unname(result$clusters$BANKSY_cluster), unname(out$BANKSY_cluster))
  expect_equal(result$summary$n_spots, 4L)
  expect_equal(result$summary$domains$domain, c("1", "2"))
})

test_that("BANKSY receipt shows resolved inputs, stored results, and a plot call", {
  srt <- make_banksy_seurat()
  messages <- testthat::capture_messages(with_mock_banksy({
    RunBANKSY(srt, layer = "counts", group = "sample", verbose = TRUE)
  }))
  plain <- cli::ansi_strip(paste(messages, collapse = "\n"))

  expect_match(plain, "Using assay")
  expect_match(plain, "Using coordinates")
  expect_match(plain, "BANKSY completed")
  expect_match(plain, "GetSpatialResult")
  expect_match(plain, "SpatialSpotPlot")
})

make_banksy_multi_image_seurat <- function() {
  srt <- make_banksy_seurat()
  suppressWarnings(srt[["slice1"]] <- SeuratObject::CreateFOV(
    data.frame(x = c(1, 2), y = c(1, 1), row.names = c("Spot1", "Spot2")),
    type = "centroids",
    assay = "RNA",
    key = "banksy1_"
  ))
  suppressWarnings(srt[["slice2"]] <- SeuratObject::CreateFOV(
    data.frame(x = c(1, 2), y = c(2, 2), row.names = c("Spot3", "Spot4")),
    type = "centroids",
    assay = "RNA",
    key = "banksy2_"
  ))
  srt
}

test_that("BANKSY fits samples independently and returns prefixed labels", {
  srt <- make_banksy_multi_image_seurat()
  suppressWarnings(with_mock_banksy({
    out <- RunBANKSY(
      srt,
      layer = "counts",
      group = "sample",
      sample.by = "sample",
      image = c(S1 = "slice1", S2 = "slice2"),
      verbose = FALSE
    )
  }))

  expect_equal(
    unname(out$BANKSY_cluster),
    c("S1_1", "S1_1", "S2_1", "S2_1")
  )
  expect_named(out@tools$BANKSY$per_sample, c("S1", "S2"))
  expect_identical(out@tools$BANKSY$parameters$sample.by, "sample")
  expect_identical(out@tools$BANKSY$parameters$image, c(S1 = "slice1", S2 = "slice2"))
  expect_equal(
    unname(GetSpatialResult(out, "BANKSY", sample = "S1")$clusters$BANKSY_cluster),
    c("1", "1")
  )

  suppressWarnings(with_mock_banksy({
    out_without_details <- RunBANKSY(
      srt,
      layer = "counts",
      group = "sample",
      sample.by = "sample",
      image = c(S1 = "slice1", S2 = "slice2"),
      store_results = FALSE,
      verbose = FALSE
    )
  }))
  expect_false("BANKSY" %in% names(out_without_details@tools))
  sample_result <- GetSpatialResult(out_without_details, "BANKSY", sample = "S1")
  expect_equal(unname(sample_result$clusters$BANKSY_cluster), c("1", "1"))
  expect_equal(sample_result$summary$n_spots, 2L)
})

test_that("BANKSY sample failures are reported without returning partial results", {
  srt <- make_banksy_seurat()
  expect_error(
    with_mock_banksy({
      RunBANKSY(
        srt,
        layer = "counts",
        group = "sample",
        sample.by = "sample",
        verbose = FALSE
      )
    }, fail_sample = "S2"),
    "BANKSY failed for sample.*S2.*mock backend failure"
  )
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
