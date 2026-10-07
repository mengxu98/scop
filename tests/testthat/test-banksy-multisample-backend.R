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
  # Prefixes distinguish independent fits; they do not imply aligned domains.
  expect_true(all(grepl("^S1_|^S1_A_", compact$domains)))
})
