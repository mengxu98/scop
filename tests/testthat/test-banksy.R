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

skip_if_missing_banksy_test_dependencies <- function() {
  testthat::skip_if_not_installed("SpatialExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("S4Vectors")
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
    banksy_resolve_assay(srt, image = "missing_slice"),
    "is not present"
  )
  expect_error(
    RunBANKSY(make_banksy_seurat(), verbose = FALSE),
    "Layer .* is not present"
  )
})

test_that("RunBANKSY auto-selects its single image and associated assay", {
  skip_if_missing_banksy_test_dependencies()
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
  skip_if_missing_banksy_test_dependencies()
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

  expect_named(out@tools$BANKSY, "result_index")
  result <- GetSpatialResult(out, "BANKSY")
  expect_null(result$parameters)
  expect_equal(unname(result$clusters$BANKSY_cluster), unname(out$BANKSY_cluster))
  expect_equal(result$summary$n_spots, 4L)
  expect_equal(result$summary$domains$domain, c("1", "2"))
})

test_that("BANKSY receipt shows resolved inputs, stored results, and a plot call", {
  skip_if_missing_banksy_test_dependencies()
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
  skip_if_missing_banksy_test_dependencies()
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
  expect_named(out_without_details@tools$BANKSY, "result_index")
  sample_result <- GetSpatialResult(out_without_details, "BANKSY", sample = "S1")
  expect_equal(unname(sample_result$clusters$BANKSY_cluster), c("1", "1"))
  expect_equal(sample_result$summary$n_spots, 2L)
})

test_that("BANKSY sample failures are reported without returning partial results", {
  skip_if_missing_banksy_test_dependencies()
  srt <- make_banksy_seurat()
  srt$BANKSY_cluster <- "old"
  original <- srt
  for (store in c(TRUE, FALSE)) {
    expect_error(
      with_mock_banksy({
        RunBANKSY(srt, layer = "counts", group = "sample", sample.by = "sample",
          store_results = store, verbose = FALSE)
      }, fail_sample = "S2"),
      "BANKSY failed for sample.*S2.*mock backend failure"
    )
    expect_identical(srt, original)
  }
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

test_that("BANKSY lightweight results retain exact samples and custom columns", {
  skip_if_missing_banksy_test_dependencies()
  srt <- make_banksy_seurat()
  srt$sample <- c("S1", "S1", "S1_A", "S1_A")
  srt$BANKSY_cluster <- rep("stale", 4)
  # An earlier default result must not override the newest custom column.
  srt@tools$BANKSY <- list(clusters = data.frame(
    BANKSY_cluster = rep("stale", 4), row.names = colnames(srt)
  ))
  for (tool in c("BANKSY", "custom_fit")) {
    with_mock_banksy({
      detailed <- RunBANKSY(srt, layer = "counts", group = "sample",
        sample.by = "sample", cluster_colname = "domains", tool_name = tool,
        store_results = TRUE, verbose = FALSE)
      compact <- RunBANKSY(detailed, layer = "counts", group = "sample",
        sample.by = "sample", cluster_colname = "domains", tool_name = tool,
        store_results = FALSE, verbose = FALSE)
    })
    expect_named(compact@tools[[tool]], "result_index")
    expect_identical(compact@tools[[tool]]$result_index$samples,
      stats::setNames(srt$sample, colnames(srt)))
    expect_equal(GetSpatialResult(compact, tool)$clusters,
      GetSpatialResult(detailed, tool)$clusters)
    for (sample_name in c("S1", "S1_A")) {
      result <- GetSpatialResult(compact, tool, sample = sample_name)
      expected <- GetSpatialResult(detailed, tool, sample = sample_name)
      expect_equal(result$clusters, expected$clusters)
      expect_equal(result$summary, expected$summary)
      expect_equal(rownames(result$clusters), colnames(srt)[srt$sample == sample_name])
      expect_equal(unname(result$clusters$domains), c("1", "1"))
      expect_null(result$parameters)
    }
    expect_error(GetSpatialResult(compact, tool, sample = "S"), "No .* results")
    expect_error(GetSpatialResult(compact, tool, sample = "missing"), "No .* results")
    # Retrieval uses the captured run membership even after metadata changes.
    compact$sample <- "changed"
    expect_equal(nrow(GetSpatialResult(compact, tool, sample = "S1")$clusters), 2L)
  }
})

test_that("BANKSY compact single fits resolve custom columns without stale labels", {
  skip_if_missing_banksy_test_dependencies()
  for (previous in c(FALSE, TRUE)) {
    srt <- make_banksy_seurat()
    if (previous) srt$BANKSY_cluster <- "stale"
    with_mock_banksy({
      detailed <- RunBANKSY(srt, layer = "counts", group = "sample",
        cluster_colname = "domains", tool_name = "custom_fit", verbose = FALSE)
      compact <- RunBANKSY(detailed, layer = "counts", group = "sample",
        cluster_colname = "domains", tool_name = "custom_fit",
        store_results = FALSE, verbose = FALSE)
    })
    result <- GetSpatialResult(compact, "custom_fit")
    expect_equal(result$clusters, GetSpatialResult(detailed, "custom_fit")$clusters)
    expect_equal(result$summary, GetSpatialResult(detailed, "custom_fit")$summary)
    expect_named(compact@tools$custom_fit, "result_index")
    expect_null(result$parameters)
    expect_error(GetSpatialResult(compact, "custom_fit", sample = "S1"),
      "No exact sample membership")
  }
})

test_that("BANKSY reruns clear stale labels from unanalyzed zero-count cells", {
  skip_if_missing_banksy_test_dependencies()
  srt <- make_banksy_seurat()
  counts <- Seurat::GetAssayData(srt, layer = "counts")
  counts[, "Spot2"] <- 0
  srt <- Seurat::CreateSeuratObject(counts = counts, meta.data = srt[[]])
  srt$sample <- c("S1", "S1", "S1_A", "S1_A")
  srt$BANKSY_cluster <- "stale"
  with_mock_banksy({
    detailed <- RunBANKSY(srt, layer = "counts", group = "sample",
      sample.by = "sample", verbose = FALSE)
    compact <- RunBANKSY(srt, layer = "counts", group = "sample",
      sample.by = "sample", store_results = FALSE, verbose = FALSE)
  })
  expect_equal(GetSpatialResult(compact, "BANKSY", sample = "S1")$clusters,
    GetSpatialResult(detailed, "BANKSY", sample = "S1")$clusters)
  expect_identical(rownames(GetSpatialResult(compact, "BANKSY", sample = "S1")$clusters), "Spot1")
  expect_false("Spot2" %in% names(compact@tools$BANKSY$result_index$samples))
  expect_true(is.na(compact$BANKSY_cluster[["Spot2"]]))
  expect_equal(GetSpatialResult(compact, "BANKSY")$clusters,
    GetSpatialResult(detailed, "BANKSY")$clusters)
  expect_equal(GetSpatialResult(compact, "BANKSY")$summary$n_spots, 3L)
  with_mock_banksy({
    single <- RunBANKSY(srt, layer = "counts", group = "sample",
      store_results = FALSE, verbose = FALSE)
  })
  expect_true(is.na(single$BANKSY_cluster[["Spot2"]]))
  expect_identical(rownames(GetSpatialResult(single, "BANKSY")$clusters),
    c("Spot1", "Spot3", "Spot4"))

})

test_that("GetSpatialResult does not guess sample membership in legacy metadata", {
  srt <- make_banksy_seurat()
  srt$BANKSY_cluster <- c("S1_Domain_1", "S1_Domain_2", "S1_A_Domain_1", "S1_A_Domain_2")
  expect_equal(nrow(GetSpatialResult(srt, "BANKSY")$clusters), 4L)
  expect_error(GetSpatialResult(srt, "BANKSY", sample = "S1"),
    "No exact sample membership")
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

test_that("BANKSY getters follow current cells and order in both storage modes", {
  skip_if_missing_banksy_test_dependencies()
  for (by_sample in c(FALSE, TRUE)) {
    srt <- make_banksy_seurat()
    srt$sample <- c("S1", "S1", "S1_A", "S1_A")
    args <- list(layer = "counts", group = "sample", cluster_colname = "domains",
      tool_name = "custom_fit", verbose = FALSE)
    if (by_sample) args$sample.by <- "sample"
    with_mock_banksy({
      detailed <- do.call(RunBANKSY, c(list(object = srt), args))
      compact <- do.call(RunBANKSY, c(list(object = srt, store_results = FALSE), args))
    })
    for (rename in c(FALSE, TRUE)) {
      # Reversing the existing IDs is a permutation, not detectable from sets.
      outputs <- list(detailed, compact)
      if (rename) outputs <- lapply(outputs, SeuratObject::RenameCells,
        new.names = rev(colnames(srt)))
      outputs <- lapply(outputs, function(x) x[, colnames(x)[c(3L, 1L)]])
      # Seurat subset may retain its original order. Reconstruct a valid
      # reordered object to ensure retrieval follows current order as well.
      outputs <- lapply(outputs, function(x) {
        cells <- rev(colnames(x))
        current <- Seurat::CreateSeuratObject(
          counts = Seurat::GetAssayData(x, layer = "counts")[, cells],
          meta.data = x[[]][cells, , drop = FALSE])
        current@tools <- x@tools
        current
      })
      results <- lapply(outputs, GetSpatialResult, method = "custom_fit")
      expect_equal(results[[1L]]$clusters, results[[2L]]$clusters)
      expect_equal(results[[1L]]$summary, results[[2L]]$summary)
      for (i in seq_along(outputs)) {
        before <- outputs[[i]]
        result <- results[[i]]
        expect_identical(rownames(result$clusters), colnames(before))
        expect_identical(result$clusters$domains, as.character(before$domains))
        expect_equal(result$summary$n_spots, 2L)
        expect_equal(sum(result$summary$domains$count), 2L)
        expect_identical(GetSpatialResult(before, "custom_fit"), result)
        expect_identical(outputs[[i]], before)
        if (by_sample) {
          sample_result <- GetSpatialResult(before, "custom_fit", sample = "S1_A")
          expect_identical(rownames(sample_result$clusters),
            colnames(before)[before$sample == "S1_A"])
          expect_identical(sample_result$clusters$domains, "1")
          expect_equal(sample_result$summary$n_spots, 1L)
        }
      }
    }
  }
})

test_that("BANKSY getters return typed empty results for absent fitted cells", {
  skip_if_missing_banksy_test_dependencies()
  srt <- make_banksy_seurat()
  counts <- Seurat::GetAssayData(srt, layer = "counts")
  counts[, "Spot2"] <- 0
  srt <- Seurat::CreateSeuratObject(counts = counts, meta.data = srt[[]])
  for (store in c(TRUE, FALSE)) {
    with_mock_banksy({
      single <- RunBANKSY(srt, layer = "counts", group = "sample",
        store_results = store, verbose = FALSE)
      multi <- RunBANKSY(srt, layer = "counts", group = "sample",
        sample.by = "sample", store_results = store, verbose = FALSE)
    })
    empty_single <- GetSpatialResult(single[, "Spot2"], "BANKSY")
    empty_sample <- GetSpatialResult(multi[, c("Spot3", "Spot4")], "BANKSY", sample = "S1")
    for (result in list(empty_single, empty_sample)) {
      expect_equal(nrow(result$clusters), 0L)
      expect_identical(colnames(result$clusters), "BANKSY_cluster")
      expect_equal(result$summary$n_spots, 0L)
      expect_identical(result$summary$domains,
        data.frame(domain = character(), count = integer()))
    }
    # With detailed storage, neither assignments nor subsets require the label
    # metadata column. The identity marker is retained independently.
    if (store) {
      single$BANKSY_cluster <- NULL
      expect_equal(nrow(GetSpatialResult(single[, "Spot2"], "BANKSY")$clusters), 0L)
      expect_equal(nrow(GetSpatialResult(single[, "Spot1"], "BANKSY")$clusters), 1L)
    } else {
      single$BANKSY_cluster <- NULL
      expect_error(GetSpatialResult(single, "BANKSY"), "No cluster results")
    }
  }
})

test_that("BANKSY identity markers are collision safe and reused on reruns", {
  skip_if_missing_banksy_test_dependencies()
  srt <- make_banksy_seurat()
  srt$.scop_BANKSY_cell_id <- "user data"
  with_mock_banksy({
    out <- RunBANKSY(srt, layer = "counts", group = "sample", verbose = FALSE)
    column <- out@tools$BANKSY$result_index$cell_id_colname
    expect_false(identical(column, ".scop_BANKSY_cell_id"))
    expect_identical(out$.scop_BANKSY_cell_id, srt$.scop_BANKSY_cell_id)
    rerun <- RunBANKSY(out, layer = "counts", group = "sample",
      store_results = FALSE, verbose = FALSE)
    expect_identical(rerun@tools$BANKSY$result_index$cell_id_colname, column)
    expect_identical(colnames(rerun[[]]), colnames(out[[]]))
    another <- RunBANKSY(rerun, layer = "counts", group = "sample",
      cluster_colname = "other_domains", tool_name = "other_fit", verbose = FALSE)
    expect_false(identical(another@tools$other_fit$result_index$cell_id_colname, column))
    before_collision <- another
    expect_error(RunBANKSY(another, layer = "counts", group = "sample",
      cluster_colname = column, tool_name = "other_fit", verbose = FALSE),
      "identity metadata column")
    expect_identical(another, before_collision)
    renamed <- SeuratObject::RenameCells(another, add.cell.id = "new")
    expect_identical(GetSpatialResult(renamed, "BANKSY")$clusters$BANKSY_cluster,
      GetSpatialResult(renamed, "other_fit")$clusters$other_domains)
    renamed_rerun <- RunBANKSY(renamed, layer = "counts", group = "sample",
      cluster_colname = "other_domains", tool_name = "other_fit",
      store_results = FALSE, verbose = FALSE)
    expect_identical(colnames(renamed_rerun[[]]), colnames(renamed[[]]))
    expect_identical(GetSpatialResult(renamed_rerun, "BANKSY")$clusters,
      GetSpatialResult(renamed, "BANKSY")$clusters)
    expect_identical(rownames(GetSpatialResult(renamed_rerun, "other_fit")$clusters),
      colnames(renamed_rerun))
    # A requested output column that collides with the old identity locator
    # must retain its cluster labels and receive a new independent locator.
    collision <- RunBANKSY(out, layer = "counts", group = "sample",
      cluster_colname = column, verbose = FALSE)
    expect_false(identical(collision@tools$BANKSY$result_index$cell_id_colname, column))
    expect_identical(GetSpatialResult(collision, "BANKSY")$clusters[[column]],
      c("1", "1", "2", "2"))
  })
})

test_that("BANKSY rejects missing or invalid tracked cell identities", {
  skip_if_missing_banksy_test_dependencies()
  for (store in c(TRUE, FALSE)) {
    with_mock_banksy({
      out <- RunBANKSY(make_banksy_seurat(), layer = "counts", group = "sample",
        store_results = store, verbose = FALSE)
    })
    column <- out@tools$BANKSY$result_index$cell_id_colname
    for (invalid in list(NULL, rep("Spot1", 4L), c(NA, "Spot2", "Spot3", "Spot4"),
      paste0("unrelated", 1:4))) {
      broken <- out
      broken@meta.data[[column]] <- invalid
      expect_error(GetSpatialResult(broken, "BANKSY"), "identities.*stale or unverifiable")
    }
    broken <- out
    broken@tools$BANKSY$result_index$cell_id_colname <- NULL
    expect_error(GetSpatialResult(broken, "BANKSY"), "identities.*stale or unverifiable")
    broken <- out
    broken@tools$BANKSY$result_index$object_cells <- rep("Spot1", 4L)
    expect_error(GetSpatialResult(broken, "BANKSY"), "identities.*stale or unverifiable")
  }
})

test_that("legacy spatial results retain exact-ID subsets and diagnose stale names", {
  skip_if_missing_banksy_test_dependencies()
  for (store in c(TRUE, FALSE)) {
    with_mock_banksy({
      out <- RunBANKSY(make_banksy_seurat(), layer = "counts", group = "sample",
        store_results = store, verbose = FALSE)
    })
    column <- out@tools$BANKSY$result_index$cell_id_colname
    out@tools$BANKSY$result_index$cell_id_colname <- NULL
    out@tools$BANKSY$result_index$object_cells <- NULL
    out@meta.data[[column]] <- NULL
    current <- out[, c("Spot4", "Spot1")]
    result <- GetSpatialResult(current, "BANKSY")
    expect_identical(rownames(result$clusters), colnames(current))
    expect_identical(result$clusters$BANKSY_cluster,
      as.character(out$BANKSY_cluster[colnames(current)]))
    expect_equal(result$summary$n_spots, 2L)
    renamed <- SeuratObject::RenameCells(out, add.cell.id = "new")
    expect_error(GetSpatialResult(renamed, "BANKSY"), "identities.*stale or unverifiable")
    if (store) {
      out$BANKSY_cluster <- NULL
      expect_equal(nrow(GetSpatialResult(out[, "Spot1"], "BANKSY")$clusters), 1L)
      renamed$BANKSY_cluster <- NULL
      expect_error(GetSpatialResult(renamed, "BANKSY"), "identities.*stale or unverifiable")
    }
  }
})

test_that("other spatial summaries are not misrepresented after subsetting", {
  srt <- make_banksy_seurat()
  srt@tools$other <- list(clusters = data.frame(cluster = c("A", "A", "B", "B"),
    row.names = colnames(srt)), parameters = list(seed = 1),
    summary = list(n_spots = 4L, fit_statistic = 10))
  expect_identical(GetSpatialResult(srt, "other")$summary, srt@tools$other$summary)
  srt@tools$other$per_sample <- list(S1 = list(
    clusters = srt@tools$other$clusters[c("Spot1", "Spot2"), , drop = FALSE],
    summary = list(n_spots = 2L, fit_statistic = 5)))
  expect_identical(GetSpatialResult(srt, "other", sample = "S1")$summary,
    srt@tools$other$per_sample$S1$summary)
  expect_null(GetSpatialResult(srt[, "Spot1"], "other", sample = "S1")$summary)
  # A valid object with a reordered count matrix/metadata exercises current
  # cell order independently of Seurat versions that preserve subset order.
  cells <- c("Spot3", "Spot1")
  current <- Seurat::CreateSeuratObject(
    counts = Seurat::GetAssayData(srt, layer = "counts")[, cells],
    meta.data = srt[[]][cells, , drop = FALSE])
  current@tools <- srt@tools
  result <- GetSpatialResult(current, "other")
  expect_identical(rownames(result$clusters), cells)
  expect_identical(result$clusters$cluster, c("B", "A"))
  expect_identical(result$parameters, list(seed = 1))
  expect_null(result$summary)
})

test_that("a BANKSY cluster column named cell is not mistaken for an identity field", {
  skip_if_missing_banksy_test_dependencies()
  with_mock_banksy({
    out <- RunBANKSY(make_banksy_seurat(), layer = "counts", group = "sample",
      cluster_colname = "cell", verbose = FALSE)
  })
  expect_spatial_result_lifecycle(out, "BANKSY", "cell")
})
