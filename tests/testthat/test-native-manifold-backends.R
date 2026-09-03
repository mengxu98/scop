test_that("native PaCMAP covers initialization and distance branches", {
  set.seed(20260730)
  data <- rbind(
    matrix(rnorm(30 * 6, mean = -2), ncol = 6),
    matrix(rnorm(30 * 6, mean = 2), ncol = 6)
  )
  rownames(data) <- paste0("Cell", seq_len(nrow(data)))

  for (metric in c("euclidean", "manhattan", "cosine", "angular", "hamming")) {
    for (init in c("pca", "random")) {
      first <- RunPaCMAP(
        data,
        backend = "cpp",
        distance_method = metric,
        init = init,
        n.neighbors = 6,
        num_iters = 20,
        seed.use = 17,
        verbose = FALSE
      )
      second <- RunPaCMAP(
        data,
        backend = "cpp",
        distance_method = metric,
        init = init,
        n.neighbors = 6,
        num_iters = 20,
        seed.use = 17,
        verbose = FALSE
      )
      expect_s4_class(first, "DimReduc")
      expect_equal(Seurat::Embeddings(first), Seurat::Embeddings(second))
      expect_true(all(is.finite(Seurat::Embeddings(first))))
    }
  }
})

test_that("native TriMap covers all optimizers and distance branches", {
  set.seed(20260731)
  data <- matrix(rnorm(50 * 8), nrow = 50)
  rownames(data) <- paste0("Cell", seq_len(nrow(data)))

  for (metric in c("euclidean", "manhattan", "cosine", "angular", "hamming")) {
    for (optimizer in c("dbd", "sd", "momentum")) {
      result <- RunTriMap(
        data,
        backend = "cpp",
        distance_method = metric,
        opt_method = optimizer,
        n_inliers = 5,
        n_outliers = 2,
        n_random = 1,
        n_iters = 20,
        seed.use = 19,
        verbose = FALSE
      )
      expect_s4_class(result, "DimReduc")
      expect_equal(dim(Seurat::Embeddings(result)), c(50L, 2L))
      expect_true(all(is.finite(Seurat::Embeddings(result))))
    }
  }
})

test_that("native BBKNN graph contains neighbors from every batch", {
  set.seed(20260801)
  embedding <- rbind(
    matrix(rnorm(20 * 5, mean = -1), ncol = 5),
    matrix(rnorm(20 * 5, mean = 1), ncol = 5)
  )
  rownames(embedding) <- paste0("Cell", seq_len(nrow(embedding)))
  batches <- rep(c("A", "B"), each = 20)
  native <- getFromNamespace("bbknn_native_matrix", "scop")(
    embedding,
    batches,
    params = list(
      neighbors_within_batch = 3,
      computation = "exact",
      metric = "euclidean",
      cores = 1
    )
  )
  distance <- native[[1]]

  for (cell in seq_len(nrow(distance))) {
    neighbors <- which(distance[cell, ] > 0)
    expect_true(all(c("A", "B") %in% unique(batches[neighbors])))
  }
  expect_true(Matrix::isSymmetric(native[[2]]))
  expect_identical(native[[3]]$backend, "cpp_exact")
})

test_that("BBKNN Annoy path uses Seurat and reports parameters", {
  set.seed(20260802)
  embedding <- matrix(rnorm(80 * 8), nrow = 80)
  rownames(embedding) <- paste0("Cell", seq_len(nrow(embedding)))
  batches <- rep(c("A", "B", "C", "D"), each = 20)
  params <- list(
    neighbors_within_batch = 3,
    computation = "annoy",
    annoy_n_trees = 10,
    metric = "euclidean",
    trim = 0,
    cores = 2
  )
  first <- getFromNamespace("bbknn_native_matrix", "scop")(
    embedding, batches, params
  )
  second <- getFromNamespace("bbknn_native_matrix", "scop")(
    embedding, batches, params
  )

  expect_equal(first[[1]], second[[1]], tolerance = 0)
  expect_equal(first[[2]], second[[2]], tolerance = 0)
  expect_identical(first[[3]]$backend, "seurat_annoy")
  expect_identical(first[[3]]$computation, "annoy")
  expect_identical(first[[3]]$annoy_n_trees, 10L)
})
