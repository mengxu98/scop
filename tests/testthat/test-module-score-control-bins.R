test_that("module score control bins preserve draw order, repeats and RNG", {
  expr <- matrix(
    c(1, 2, 3, 4, 5, 6, 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5),
    nrow = 6L,
    dimnames = list(paste0("g", 1:6), paste0("s", 1:3))
  )
  expr <- Matrix::Matrix(expr, sparse = TRUE)
  features <- list(first = c("g5", "g2", "g3", "g6", "g5"), second = c("g3", "g2"),
                   missing_pool = "g6")
  pool <- c("g5", "g4", "g2", "g3", "g1")
  seed <- 37L
  ctrl <- 1L
  set.seed(seed)
  avg <- Matrix::rowMeans(expr[pool, , drop = FALSE])
  avg <- avg[order(avg)]
  bins <- ggplot2::cut_number(avg + stats::rnorm(length(avg)) / 1e+30,
                              n = 2L, labels = FALSE, right = FALSE)
  names(bins) <- names(avg)
  features_unique <- lapply(features, unique)
  ref_controls <- lapply(features_unique, function(gs) {
    blocks <- lapply(gs, function(g) bins[which(bins == bins[g])])
    vector <- unlist(unname(blocks))
    match(names(sample(vector, size = min(ctrl * length(gs), length(vector)))), rownames(expr))
  })
  reference_seed <- .Random.seed
  result <- module_score_sparse(expr, lapply(features_unique, match, table = rownames(expr)),
                                ref_controls, n_threads = 1L)
  dimnames(result) <- list(colnames(expr), names(features))
  observed <- run_seurat_module_scores(expr, features, pool = pool, nbin = 2L,
                                       ctrl = ctrl, seed = seed, cores = 1L)
  expect_identical(observed, result, num.eq = FALSE)
  expect_identical(.Random.seed, reference_seed)
})
