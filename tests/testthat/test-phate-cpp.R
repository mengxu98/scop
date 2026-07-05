test_that("phate_potential_distance_cpp returns a square matrix with landmarks", {
  set.seed(1)
  transition <- matrix(stats::runif(25), nrow = 5)
  transition <- transition / rowSums(transition)
  log_transition <- log(transition + 1e-10)

  dist <- phate_potential_distance_cpp(log_transition, n_landmarks = 2L)

  expect_equal(dim(dist), c(5L, 5L))
  expect_equal(diag(dist), rep(0, 5), tolerance = 1e-12)
  expect_equal(dist, t(dist), tolerance = 1e-12)
  expect_true(all(is.finite(dist)))
})

test_that("phate_metric_mds_cpp handles one-dimensional requests", {
  dist <- matrix(
    c(
      0, 1, 2,
      1, 0, 1,
      2, 1, 0
    ),
    nrow = 3,
    byrow = TRUE
  )

  embedding <- phate_metric_mds_cpp(dist, n_components = 1L)

  expect_equal(dim(embedding), c(3L, 1L))
  expect_true(all(is.finite(embedding)))
})

test_that("phate_metric_mds_cpp handles single-cell input", {
  embedding <- phate_metric_mds_cpp(matrix(0, nrow = 1, ncol = 1), n_components = 1L)

  expect_equal(dim(embedding), c(1L, 1L))
  expect_equal(embedding[1, 1], 0)
})

test_that("PHATE C++ helpers reject invalid matrix contracts", {
  expect_error(
    phate_potential_distance_cpp(matrix(0, nrow = 3, ncol = 2), n_landmarks = -1L),
    "square matrix"
  )
  expect_error(
    phate_metric_mds_cpp(matrix(0, nrow = 3, ncol = 2), n_components = 2L),
    "square distance matrix"
  )
  expect_error(
    phate_metric_mds_cpp(matrix(0, nrow = 2, ncol = 2), n_components = 0L),
    "at least 1"
  )
})

test_that("RunPHATE cpp backend stores a finite reduction", {
  data("pancreas_sub", package = "scop")
  srt <- suppressWarnings(standard_scop(pancreas_sub[, seq_len(60)], verbose = FALSE))

  out <- RunPHATE(
    srt,
    reduction = "Standardpca",
    dims = 1:5,
    backend = "cpp",
    n_components = 2,
    knn = 5,
    n_landmark = 60,
    t = 3,
    t_max = 10,
    mds = "classic",
    reduction.name = "phate_cpp",
    reduction.key = "PHATECPP_",
    verbose = FALSE
  )

  emb <- Seurat::Embeddings(out[["phate_cpp"]])
  expect_equal(dim(emb), c(60L, 2L))
  expect_true(all(is.finite(emb)))
  expect_equal(SeuratObject::Misc(out[["phate_cpp"]], slot = "backend"), "cpp")
  expect_equal(SeuratObject::Misc(out[["phate_cpp"]], slot = "t"), 3L)
})
