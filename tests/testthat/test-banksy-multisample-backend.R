test_that("real BANKSY independent sample fits match standalone fits in both storage modes", {
  skip_if_not_installed("Banksy")
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("BiocNeighbors")

  set.seed(731)
  n_per_sample <- 48L
  counts <- matrix(rpois(30L * n_per_sample * 2L, lambda = 4), nrow = 30L)
  counts[1:10, rep(seq_len(n_per_sample) <= 24L, 2L)] <-
    counts[1:10, rep(seq_len(n_per_sample) <= 24L, 2L)] + 6L
  dimnames(counts) <- list(paste0("Gene", seq_len(nrow(counts))),
    paste0("Spot", seq_len(ncol(counts))))
  srt <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  )
  srt$sample <- rep(c("S1", "S1_A"), each = n_per_sample)
  # Samples deliberately share coordinates. They must be fitted separately.
  srt$col <- rep(rep(seq_len(8L), 6L), 2L)
  srt$row <- rep(rep(seq_len(6L), each = 8L), 2L)
  srt$BANKSY_cluster <- "old"
  args <- list(layer = "counts", k_geom = 4L, npcs = 5L, k_neighbors = 5L,
    algo = "leiden", resolution = 0.6, seed = 1L, cluster_colname = "domains",
    tool_name = "custom_fit", verbose = FALSE)
  detailed <- do.call(RunBANKSY, c(list(object = srt, sample.by = "sample"), args))
  compact <- do.call(RunBANKSY, c(list(object = detailed, sample.by = "sample",
    store_results = FALSE), args))

  expect_named(compact@tools$custom_fit, "result_index")
  expect_equal(GetSpatialResult(compact, "custom_fit")$clusters,
    GetSpatialResult(detailed, "custom_fit")$clusters)
  for (sample_name in c("S1", "S1_A")) {
    cells <- colnames(srt)[srt$sample == sample_name]
    standalone <- do.call(RunBANKSY, c(list(object = srt[, cells]), args))
    observed <- GetSpatialResult(compact, "custom_fit", sample = sample_name)
    stored <- GetSpatialResult(detailed, "custom_fit", sample = sample_name)
    expect_identical(rownames(observed$clusters), cells)
    expect_equal(observed$clusters, stored$clusters)
    expect_equal(observed$clusters, GetSpatialResult(standalone, "custom_fit")$clusters)
    expect_equal(observed$summary, stored$summary)
    expect_false(anyNA(observed$clusters$domains))
    expect_gt(length(unique(observed$clusters$domains)), 0L)
    expect_identical(as.character(compact$domains[cells]),
      paste0(sample_name, "_", observed$clusters$domains))
  }
  # Sample selection still uses exact fitted membership after a rename and
  # arbitrary current-object ordering, including an absent sample.
  renamed <- lapply(list(detailed, compact), SeuratObject::RenameCells,
    new.names = rev(colnames(srt)))
  positions <- c(70L, 2L, 53L, 10L)
  renamed <- lapply(renamed, function(x) x[, colnames(x)[positions]])
  for (sample_name in c("S1", "S1_A")) {
    results <- lapply(renamed, GetSpatialResult, method = "custom_fit", sample = sample_name)
    expect_equal(results[[1L]]$clusters, results[[2L]]$clusters)
    expect_equal(results[[1L]]$summary, results[[2L]]$summary)
    expected_cells <- colnames(renamed[[1L]])[renamed[[1L]]$sample == sample_name]
    expect_identical(rownames(results[[1L]]$clusters), expected_cells)
    expect_equal(results[[1L]]$summary$n_spots, 2L)
  }
  for (out in renamed) {
    only_s1 <- out[, colnames(out)[out$sample == "S1"]]
    empty <- GetSpatialResult(only_s1, "custom_fit", sample = "S1_A")
    expect_equal(nrow(empty$clusters), 0L)
    expect_equal(empty$summary$n_spots, 0L)
    expect_equal(nrow(empty$summary$domains), 0L)
  }
  # Prefixes distinguish independent fits; they do not imply aligned domains.
  expect_true(all(grepl("^S1_|^S1_A_", compact$domains)))
})

test_that("real BANKSY results survive subsets and exact renamed-cell mapping", {
  skip_if_not_installed("Banksy")
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("BiocNeighbors")

  set.seed(7)
  counts <- matrix(rpois(12L * 36L, lambda = 5), nrow = 12L,
    dimnames = list(paste0("Gene", 1:12), paste0("Spot", 1:36)))
  srt <- Seurat::CreateSeuratObject(counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix"))
  srt$x <- rep(1:6, 6)
  srt$y <- rep(1:6, each = 6)
  args <- list(layer = "counts", k_geom = 4L, npcs = 4L, k_neighbors = 4L,
    seed = 1L, verbose = FALSE)
  detailed <- do.call(RunBANKSY, c(list(object = srt), args))
  compact <- do.call(RunBANKSY, c(list(object = srt, store_results = FALSE), args))
  baseline <- GetSpatialResult(detailed, "BANKSY")
  expect_equal(GetSpatialResult(compact, "BANKSY")$clusters, baseline$clusters)
  expect_equal(baseline$summary$n_spots, 36L)
  expect_identical(baseline$summary, detailed@tools$BANKSY$summary)
  positions <- c(36L, 2L, 14L, 7L, 25L, 11L, 19L, 3L, 30L, 16L, 22L, 5L)
  for (rename in c("none", "new_names", "permutation")) {
    outputs <- list(detailed, compact)
    if (rename != "none") {
      new_names <- if (rename == "permutation") rev(colnames(srt)) else paste0("renamed", 1:36)
      outputs <- lapply(outputs, SeuratObject::RenameCells, new.names = new_names)
    }
    for (take_subset in c(FALSE, TRUE)) {
      current <- if (take_subset) lapply(outputs, function(x) x[, colnames(x)[positions]]) else outputs
      results <- lapply(current, GetSpatialResult, method = "BANKSY")
      expect_equal(results[[1L]]$clusters, results[[2L]]$clusters)
      expect_equal(results[[1L]]$summary, results[[2L]]$summary)
      for (i in seq_along(current)) {
        identity_column <- current[[i]]@tools$BANKSY$result_index$cell_id_colname
        fitted_ids <- current[[i]]@meta.data[[identity_column]]
        expected_labels <- baseline$clusters[fitted_ids, "BANKSY_cluster"]
        result <- results[[i]]
        expect_identical(rownames(result$clusters), colnames(current[[i]]))
        expect_identical(result$clusters$BANKSY_cluster, expected_labels)
        expect_equal(result$summary$n_spots, ncol(current[[i]]))
        expect_equal(sum(result$summary$domains$count), ncol(current[[i]]))
        expected_counts <- table(expected_labels)
        expect_equal(result$summary$domains$count,
          unname(as.integer(expected_counts[result$summary$domains$domain])))
        expect_identical(GetSpatialResult(current[[i]], "BANKSY"), result)
      }
    }
  }
})
