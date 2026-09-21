
expect_close_to_reference <- function(actual, reference, atol = 1e-12, rtol = 1e-10) {
  reference <- as.numeric(reference)
  actual <- as.numeric(actual)
  expect_equal(length(actual), length(reference))
  ok <- abs(actual - reference) <= atol + rtol * abs(reference)
  bad <- which(!is.na(ok) & !ok)
  if (length(bad) > 0L) {
    fail(sprintf(
      "%d of %d values outside atol + rtol * |reference|; worst = %.6g vs %.6g (index %d)",
      length(bad), length(ok),
      actual[bad[which.max(abs(actual[bad] - reference[bad]))]],
      reference[bad[which.max(abs(actual[bad] - reference[bad]))]],
      bad[which.max(abs(actual[bad] - reference[bad]))]
    ))
  }
  expect_true(all(is.na(ok) == is.na(reference)))
}

test_that("Palantir Markov chain uses the reference bandwidth and start state", {
  set.seed(9)
  n <- 40L
  n_dims <- 4L
  wp_data <- matrix(stats::rnorm(n * n_dims), nrow = n)
  pseudotime <- as.numeric(scale(rowSums(wp_data)))
  knn <- 12L

  actual <- palantir_markov_chain_cpp(wp_data, knn, pseudotime)

  dists <- as.matrix(stats::dist(wp_data))
  idx <- t(apply(dists, 1L, order))[, seq_len(knn), drop = FALSE]
  d <- matrix(dists[cbind(rep(seq_len(n), each = knn), as.vector(t(idx)))],
    nrow = n, byrow = TRUE
  )
  adaptive_k <- max(1L, min(as.integer(floor(knn / 3)) - 1L, knn - 1L))
  adaptive_std <- apply(d, 1L, function(row) sort(row)[adaptive_k + 1L])

  reference_i <- integer()
  reference_j <- integer()
  reference_x <- numeric()
  for (i in seq_len(n)) {
    cutoff <- pseudotime[i] - adaptive_std[i]
    keep <- which(pseudotime[idx[i, ]] >= cutoff)
    if (length(keep) == 0L) {
      reference_i <- c(reference_i, i)
      reference_j <- c(reference_j, i)
      reference_x <- c(reference_x, 1)
      next
    }
    nb <- idx[i, keep]
    w <- exp(-0.5 * d[i, keep]^2 *
      (1 / adaptive_std[i]^2 + 1 / adaptive_std[nb]^2))
    reference_i <- c(reference_i, rep(i, length(nb)))
    reference_j <- c(reference_j, nb)
    reference_x <- c(reference_x, w / sum(w))
  }

  order_actual <- order(actual$T_i, actual$T_j)
  order_reference <- order(reference_i, reference_j)
  expect_identical(actual$T_i[order_actual], as.integer(reference_i[order_reference]))
  expect_identical(actual$T_j[order_actual], as.integer(reference_j[order_reference]))
  expect_close_to_reference(
    actual$T_x[order_actual], reference_x[order_reference],
    atol = 1e-12, rtol = 1e-10
  )
})

test_that("CellRank lineage drivers use the Fisher z p-value", {
  set.seed(31)
  n_genes <- 30L
  n_cells <- 120L
  n_lineages <- 2L
  expression <- matrix(stats::rnorm(n_genes * n_cells), nrow = n_genes)
  abs_probs <- matrix(stats::runif(n_cells * n_lineages), nrow = n_cells)

  actual <- cellrank_lineage_drivers_cpp(expression, abs_probs)

  corr <- as.numeric(actual$correlation)
  pval <- as.numeric(actual$pval)
  z <- atanh(corr) * sqrt(n_cells - 3)
  reference <- 2 * stats::pnorm(-abs(z))

  expect_close_to_reference(pval, reference, atol = 1e-12, rtol = 1e-10)
  expect_true(all(pval <= 1 & pval >= 0))
  expect_true(any(corr < 0))
  expect_true(all(pval[corr < 0] < 1))
})

test_that("scVelo stochastic gamma honours the reference percentile mask", {
  set.seed(12)
  n_cells <- 150L
  n_genes <- 6L

  reference_gamma <- function(expr_s, expr_u, expr_ss, expr_us, masked = TRUE) {
    vapply(seq_len(n_genes), function(g) {
      s <- expr_s[g, ]
      u <- expr_u[g, ]
      var_ss <- 2 * expr_ss[g, ] - s
      cov_us <- 2 * expr_us[g, ] + u
      std0 <- function(x) sqrt(mean((x - mean(x))^2))

      s_norm <- s / max(s, 1e-3)
      u_norm <- u / max(u, 1e-3)
      normalized <- s_norm + u_norm
      det_keep <- if (masked) {
        normalized >= stats::quantile(normalized, 0.95)
      } else {
        rep(TRUE, length(normalized))
      }
      gamma_det <- sum(s[det_keep] * u[det_keep]) / sum(s[det_keep]^2)
      if (!is.finite(gamma_det) || gamma_det < 0) gamma_det <- 0

      res_std <- std0(u - gamma_det * s)
      if (!is.finite(res_std) || res_std < 1e-8) {
        return(gamma_det)
      }

      gamma2 <- sum(var_ss * cov_us) / max(sum(var_ss^2), 1e-12)
      res2_std <- std0(cov_us - gamma2 * var_ss)
      if (!is.finite(res2_std) || res2_std < 1e-8) res2_std <- 1
      keep <- if (masked) {
        normalized >= stats::quantile(normalized, 0.95) |
          s >= stats::quantile(s, 0.95)
      } else {
        rep(TRUE, length(s))
      }
      x1 <- ifelse(keep, s / res_std, 0)
      y1 <- ifelse(keep, u / res_std, 0)
      x2 <- var_ss / res2_std
      y2 <- cov_us / res2_std
      den <- sum(x1^2) + sum(x2^2)
      if (den <= 1e-12) 0 else (sum(x1 * y1) + sum(x2 * y2)) / den
    }, numeric(1))
  }

  Ms <- matrix(0, n_genes, n_cells)
  Mu <- matrix(0, n_genes, n_cells)
  Mss <- matrix(0, n_genes, n_cells)
  Mus <- matrix(0, n_genes, n_cells)
  for (g in seq_len(n_genes)) {
    s <- stats::rgamma(n_cells, shape = 2, rate = 1)
    u <- 0.7 * s + stats::rnorm(n_cells, 0, 0.4)
    Ms[g, ] <- s
    Mu[g, ] <- u
    Mss[g, ] <- s^2 * 1.1
    Mus[g, ] <- s * u * 0.9
  }

  actual <- scanpy_stochastic_cpp(
    Ms = Ms, Mu = Mu, Mss = Mss, Mus = Mus,
    knn_idx = matrix(seq_len(n_cells), ncol = 1L),
    embedding = cbind(Ms[1L, ])
  )

  reference <- reference_gamma(Ms, Mu, Mss, Mus, masked = TRUE)
  unmasked <- reference_gamma(Ms, Mu, Mss, Mus, masked = FALSE)
  expect_gt(max(abs(reference - unmasked)), 1e-6)
  expect_close_to_reference(actual$gamma, reference, atol = 1e-10, rtol = 1e-8)
})

test_that("ssGSEA truncates tied average ranks like GSVA::ssgsea", {
  set.seed(5)
  expr <- matrix(0, nrow = 30L, ncol = 12L)
  expr[1:10, ] <- rep(c(0, 1, 1, 2, 2, 2), each = 10)[seq_len(10)]
  expr[11:30, ] <- matrix(stats::rpois(20L * 12L, 3), nrow = 20L)
  expr <- methods::as(Matrix::Matrix(expr, sparse = TRUE), "dgCMatrix")
  rownames(expr) <- paste0("g", seq_len(nrow(expr)))
  sets <- list(c(2L, 3L, 5L, 7L, 11L, 13L, 17L, 19L, 23L, 29L))

  actual <- ssgsea_rank_dense(expr, sets, alpha = 1, normalize = FALSE)

  n_genes <- nrow(expr)
  x <- as.numeric(expr[, 1L])
  set <- sets[[1L]]
  ord <- order(x, seq_len(n_genes))
  rank_by_gene <- numeric(n_genes)
  weight_by_gene <- numeric(n_genes)
  i <- 1L
  while (i <= n_genes) {
    j <- i
    while (j < n_genes && x[ord[j + 1L]] == x[ord[i]]) j <- j + 1L
    rv <- as.integer((i + j) / 2)
    rank_by_gene[ord[i:j]] <- rv
    weight_by_gene[ord[i:j]] <- abs(rv)^1
    i <- j + 1L
  }
  ranking <- order(-rank_by_gene, seq_len(n_genes))
  position <- integer(n_genes)
  position[ranking] <- seq_len(n_genes)
  inv <- n_genes - position[set] + 1L
  reference <- sum(weight_by_gene[set] * inv) / sum(weight_by_gene[set]) -
    (n_genes * (n_genes + 1) / 2 - sum(inv)) / (n_genes - length(set))

  expect_close_to_reference(actual[1L, 1L], reference, atol = 1e-12, rtol = 1e-10)
})

test_that("thread counts are named cores and forwarded under the backend name", {
  for (fn in list(RunSmoothClust, RunSecAct, RunscMalignantFinder, RunMERINGUE, RunSpatialEcoTyper, RunCHOIR)) {
    expect_true("cores" %in% names(formals(fn)))
  }
  for (fn in list(RunLargeVis.Seurat, RunLargeVis.default, RunUMAP2.Seurat, RunUMAP2.default)) {
    expect_true("cores" %in% names(formals(fn)))
  }
  expect_false(any(c("n_threads", "ncores", "n_thread") %in% names(formals(RunLargeVis.default))))
  expect_false(any(c("n_threads", "ncores", "n_thread") %in% names(formals(RunUMAP2.default))))
  expect_false("n_cores" %in% names(formals(RunCHOIR)))
  expect_false("n_cores" %in% names(formals(RunSpatialEcoTyper)))
})

test_that("RunPCA accepts cores and keeps the embedding stable", {
  skip_if_not_installed("Seurat")
  set.seed(17)
  X <- matrix(stats::rnorm(300L * 400L), nrow = 300L, ncol = 400L)
  object <- Seurat::CreateSeuratObject(
    Matrix::Matrix(abs(X), sparse = TRUE),
    assay = "RNA"
  )
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object <- Seurat::FindVariableFeatures(object, nfeatures = 40L, verbose = FALSE)
  object <- Seurat::ScaleData(object, verbose = FALSE)

  single <- RunPCA(object, npcs = 8L, cores = 1L, backend = "cpp", verbose = FALSE)
  multi <- RunPCA(object, npcs = 8L, cores = 2L, backend = "cpp", verbose = FALSE)

  scaled <- as.matrix(SeuratObject::LayerData(object[["RNA"]], layer = "scale.data"))
  expect_gte(ncol(scaled), 8L * nrow(scaled))
  native_one <- pca_backend_run(scaled, 8L, TRUE, n_threads = 1L)
  native_many <- pca_backend_run(scaled, 8L, TRUE, n_threads = 2L)
  expect_type(native_one, "list")
  expect_close_to_reference(
    as.numeric(native_many$embeddings),
    as.numeric(native_one$embeddings),
    atol = 1e-8, rtol = 1e-6
  )
  expect_close_to_reference(
    as.numeric(native_one$sdev),
    sqrt(pmax(as.numeric(native_one$eigvals), 0)) / sqrt(ncol(scaled) - 1),
    atol = 1e-10, rtol = 1e-10
  )

  expect_s4_class(single[["pca"]], "DimReduc")
  expect_close_to_reference(
    as.numeric(SeuratObject::Embeddings(multi[["pca"]])),
    as.numeric(SeuratObject::Embeddings(single[["pca"]])),
    atol = 1e-8, rtol = 1e-6
  )
})

test_that("RunSCVELO rejects transition matrices beyond max_dense_gib", {
  expect_error(
    assert_cpp_dense_budget(
      n_rows = 50000L, n_cols = 50000L, copies = 2,
      max_dense_gib = 8, context = "RunSCVELO(backend = \"cpp\")"
    ),
    "dense"
  )
  expect_silent(
    assert_cpp_dense_budget(
      n_rows = 2000L, n_cols = 2000L, copies = 2,
      max_dense_gib = 8, context = "RunSCVELO(backend = \"cpp\")"
    )
  )
})

test_that("RunCIBERSORT takes its thread count as cores only", {
  expect_true("cores" %in% names(formals(RunCIBERSORT)))
  expect_false("n_threads" %in% names(formals(RunCIBERSORT)))
  expect_false("n_threads" %in% names(formals(run_aucell_official_scores)))
})

test_that("ssGSEA truncates tied average ranks like GSVA::ssgsea", {
  set.seed(5)
  expr <- matrix(0, nrow = 30L, ncol = 12L)
  expr[1:10, ] <- rep(c(0, 1, 1, 2, 2, 2), each = 10)[seq_len(10)]
  expr[11:30, ] <- matrix(stats::rpois(20L * 12L, 3), nrow = 20L)
  expr <- methods::as(Matrix::Matrix(expr, sparse = TRUE), "dgCMatrix")
  rownames(expr) <- paste0("g", seq_len(nrow(expr)))
  sets <- list(c(2L, 3L, 5L, 7L, 11L, 13L, 17L, 19L, 23L, 29L))

  actual <- ssgsea_rank_dense(expr, sets, alpha = 1, normalize = FALSE)

  n_genes <- nrow(expr)
  x <- as.numeric(expr[, 1L])
  set <- sets[[1L]]
  ord <- order(x, seq_len(n_genes))
  rank_by_gene <- numeric(n_genes)
  weight_by_gene <- numeric(n_genes)
  i <- 1L
  while (i <= n_genes) {
    j <- i
    while (j < n_genes && x[ord[j + 1L]] == x[ord[i]]) j <- j + 1L
    rv <- as.integer((i + j) / 2)
    rank_by_gene[ord[i:j]] <- rv
    weight_by_gene[ord[i:j]] <- abs(rv)^1
    i <- j + 1L
  }
  ranking <- order(-rank_by_gene, seq_len(n_genes))
  position <- integer(n_genes)
  position[ranking] <- seq_len(n_genes)
  inv <- n_genes - position[set] + 1L
  reference <- sum(weight_by_gene[set] * inv) / sum(weight_by_gene[set]) -
    (n_genes * (n_genes + 1) / 2 - sum(inv)) / (n_genes - length(set))

  expect_close_to_reference(actual[1L, 1L], reference, atol = 1e-12, rtol = 1e-10)
})

test_that("thread counts are named cores and forwarded under the backend name", {
  for (fn in list(RunSmoothClust, RunSecAct, RunscMalignantFinder, RunMERINGUE, RunSpatialEcoTyper, RunCHOIR)) {
    expect_true("cores" %in% names(formals(fn)))
  }
  for (fn in list(RunLargeVis.Seurat, RunLargeVis.default, RunUMAP2.Seurat, RunUMAP2.default)) {
    expect_true("cores" %in% names(formals(fn)))
  }
  expect_false(any(c("n_threads", "ncores", "n_thread") %in% names(formals(RunLargeVis.default))))
  expect_false(any(c("n_threads", "ncores", "n_thread") %in% names(formals(RunUMAP2.default))))
  expect_false("n_cores" %in% names(formals(RunCHOIR)))
  expect_false("n_cores" %in% names(formals(RunSpatialEcoTyper)))
})

test_that("RunPCA accepts cores and keeps the embedding stable", {
  skip_if_not_installed("Seurat")
  set.seed(17)
  X <- matrix(stats::rnorm(300L * 400L), nrow = 300L, ncol = 400L)
  object <- Seurat::CreateSeuratObject(
    Matrix::Matrix(abs(X), sparse = TRUE),
    assay = "RNA"
  )
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object <- Seurat::FindVariableFeatures(object, nfeatures = 40L, verbose = FALSE)
  object <- Seurat::ScaleData(object, verbose = FALSE)

  single <- RunPCA(object, npcs = 8L, cores = 1L, backend = "cpp", verbose = FALSE)
  multi <- RunPCA(object, npcs = 8L, cores = 2L, backend = "cpp", verbose = FALSE)

  scaled <- as.matrix(SeuratObject::LayerData(object[["RNA"]], layer = "scale.data"))
  expect_gte(ncol(scaled), 8L * nrow(scaled))
  native_one <- pca_backend_run(scaled, 8L, TRUE, n_threads = 1L)
  native_many <- pca_backend_run(scaled, 8L, TRUE, n_threads = 2L)
  expect_type(native_one, "list")
  expect_close_to_reference(
    as.numeric(native_many$embeddings),
    as.numeric(native_one$embeddings),
    atol = 1e-8, rtol = 1e-6
  )
  expect_close_to_reference(
    as.numeric(native_one$sdev),
    sqrt(pmax(as.numeric(native_one$eigvals), 0)) / sqrt(ncol(scaled) - 1),
    atol = 1e-10, rtol = 1e-10
  )

  expect_s4_class(single[["pca"]], "DimReduc")
  expect_close_to_reference(
    as.numeric(SeuratObject::Embeddings(multi[["pca"]])),
    as.numeric(SeuratObject::Embeddings(single[["pca"]])),
    atol = 1e-8, rtol = 1e-6
  )
})

test_that("RunSCVELO rejects transition matrices beyond max_dense_gib", {
  expect_error(
    assert_cpp_dense_budget(
      n_rows = 50000L, n_cols = 50000L, copies = 2,
      max_dense_gib = 8, context = "RunSCVELO(backend = \"cpp\")"
    ),
    "dense"
  )
  expect_silent(
    assert_cpp_dense_budget(
      n_rows = 2000L, n_cols = 2000L, copies = 2,
      max_dense_gib = 8, context = "RunSCVELO(backend = \"cpp\")"
    )
  )
})

test_that("the ssGSEA walk keeps the GSVA integer rank convention", {
  set.seed(9)
  expr <- matrix(stats::rpois(40L * 8L, 2), nrow = 40L)
  expr[1:20, ] <- rep(c(0, 1, 1, 2), each = 10)[seq_len(20)]
  expr <- methods::as(Matrix::Matrix(expr, sparse = TRUE), "dgCMatrix")
  rownames(expr) <- paste0("g", seq_len(nrow(expr)))
  sets <- list(c(1L, 3L, 5L, 9L, 15L, 21L, 33L))
  actual <- ssgsea_rank_dense(expr, sets, alpha = 1, normalize = FALSE)

  n_genes <- nrow(expr)
  x <- as.numeric(expr[, 1L])
  ord <- order(x, seq_len(n_genes))
  rank_by_gene <- integer(n_genes)
  i <- 1L
  while (i <= n_genes) {
    j <- i
    while (j < n_genes && x[ord[j + 1L]] == x[ord[i]]) j <- j + 1L
    rank_by_gene[ord[i:j]] <- as.integer((i + j) / 2)
    i <- j + 1L
  }
  ranking <- order(-rank_by_gene, seq_len(n_genes))
  position <- integer(n_genes)
  position[ranking] <- seq_len(n_genes)
  inv <- n_genes - position[sets[[1L]]] + 1L
  w <- abs(rank_by_gene[sets[[1L]]])
  reference <- sum(w * inv) / sum(w) -
    (n_genes * (n_genes + 1) / 2 - sum(inv)) / (n_genes - length(sets[[1L]]))

  expect_close_to_reference(actual[1L, 1L], reference, atol = 1e-12, rtol = 1e-10)
})

test_that("scVelo terminal states and pseudotime match an independent dense solve", {
  set.seed(23)
  n_cells <- 70L
  n_dims <- 4L
  knn_k <- 8L
  embedding <- matrix(stats::rnorm(n_cells * n_dims), n_cells, n_dims)
  velocity_embedding <- matrix(stats::rnorm(n_cells * n_dims), n_cells, n_dims)
  knn_idx <- t(vapply(seq_len(n_cells), function(i) {
    d <- sqrt(rowSums((embedding - matrix(embedding[i, ], n_cells, n_dims, byrow = TRUE))^2))
    order(d)[seq_len(knn_k)]
  }, integer(knn_k)))

  graph <- scanpy_velocity_graph_cpp(
    Ms = matrix(stats::rnorm(2L * n_cells), nrow = 2L),
    Mu = matrix(stats::rnorm(2L * n_cells), nrow = 2L),
    residual = matrix(stats::rnorm(2L * n_cells), nrow = 2L),
    knn_idx = knn_idx, n_neighbors_velo = knn_k,
    sqrt_transform = FALSE, n_recurse_neighbors = 1L, n_threads = 1L
  )
  args <- list(
    graph_rows = graph$velocity_graph_rows, graph_cols = graph$velocity_graph_cols,
    graph_vals = graph$velocity_graph_vals,
    graph_neg_rows = graph$velocity_graph_neg_rows,
    graph_neg_cols = graph$velocity_graph_neg_cols,
    graph_neg_vals = graph$velocity_graph_neg_vals,
    knn_idx = knn_idx
  )

  eigen_oriented <- function(backward, scale = 10) {
    A <- matrix(0, n_cells, n_cells)
    for (k in seq_along(args$graph_rows)) {
      i <- args$graph_rows[k] + 1L; j <- args$graph_cols[k] + 1L
      A[i, j] <- A[i, j] + expm1(args$graph_vals[k] * scale)
    }
    for (k in seq_along(args$graph_neg_rows)) {
      i <- args$graph_neg_rows[k] + 1L; j <- args$graph_neg_cols[k] + 1L
      A[i, j] <- A[i, j] + exp(args$graph_neg_vals[k] * scale)
    }
    oriented <- if (backward) t(A) else A
    rs <- rowSums(oriented)
    for (i in seq_len(n_cells)) {
      if (rs[i] > 1e-12 && is.finite(rs[i])) oriented[i, ] <- oriented[i, ] / rs[i]
    }
    t(oriented)
  }

  leading_components <- function(E, eps = 1e-3, k = 10L) {
    eig <- eigen(E, symmetric = FALSE)
    values <- Re(eig$values)
    keep <- which(is.finite(values))
    keep <- keep[order(values[keep], decreasing = TRUE)][seq_len(min(k, length(keep)))]
    keep <- keep[values[keep] >= 1 - eps]
    if (length(keep) == 0L) return(matrix(0, nrow(E), 0L))
    out <- matrix(0, nrow(E), length(keep))
    for (c in seq_along(keep)) {
      values_abs <- abs(eig$vectors[, keep[c]])
      lower <- stats::quantile(values_abs, 0.02, names = FALSE)
      upper <- stats::quantile(values_abs, 0.98, names = FALSE)
      v <- ifelse(values_abs < lower, 0, pmin(values_abs, upper))
      if (max(v) > 0) v <- v / max(v)
      out[, c] <- v
    }
    out
  }

  smooth <- function(score) {
    vapply(seq_len(n_cells), function(i) {
      nb <- sort(unique(c(i, knn_idx[i, ])))
      mean(score[nb])
    }, numeric(1))
  }
  clip_scale <- function(x) {
    sorted <- sort(x[is.finite(x)])
    upper <- sorted[min(length(sorted), floor(0.98 * (length(sorted) - 1L)) + 1L)]
    out <- pmax(0, pmin(x, upper))
    rng <- max(out) - min(out)
    if (rng > 1e-12) out <- (out - min(out)) / rng
    out
  }

  ts <- do.call(scanpy_terminal_states_graph_cpp, args)
  expect_length(ts$root_cells, n_cells)
  expect_true(all(is.finite(ts$root_cells)))
  expect_true(all(ts$root_cells >= 0 & ts$root_cells <= 1))
  expect_lte(ts$n_root_regions, 10L)
  expect_lte(ts$n_end_regions, 10L)

  ref_roots <- clip_scale(smooth(rowSums(leading_components(eigen_oriented(TRUE)))))
  ref_ends <- clip_scale(smooth(rowSums(leading_components(eigen_oriented(FALSE)))))
  expect_close_to_reference(ts$root_cells, ref_roots, atol = 1e-6, rtol = 1e-6)
  expect_close_to_reference(ts$end_points, ref_ends, atol = 1e-6, rtol = 1e-6)

  pt <- do.call(scanpy_pseudotime_graph_cpp,
    c(args, list(root_cells = ts$root_cells, end_points = ts$end_points, n_dcs = 10L)))
  expect_length(pt$pseudotime, n_cells)
  expect_true(all(is.finite(pt$pseudotime)))
  expect_true(all(pt$pseudotime >= 0 & pt$pseudotime <= 1))
  expect_true(pt$root_cell >= 1L && pt$root_cell <= n_cells)
  expect_true(pt$end_cell >= 1L && pt$end_cell <= n_cells)
  expect_identical(dim(pt$diffusion_components), c(n_cells, 10L))

  direct <- scanpy_pseudotime_cpp(
    velocity_embedding = velocity_embedding, embedding = embedding,
    knn_idx = knn_idx, root_cells = ts$root_cells, end_points = ts$end_points,
    n_neighbors_velo = knn_k
  )
  expect_length(direct$pseudotime, n_cells)
  expect_true(all(is.finite(direct$pseudotime)))
  expect_true(all(direct$pseudotime >= 0 & direct$pseudotime <= 1))
})

test_that("the deprecated single_data alias warns before it is used", {
  expect_warning(
    tryCatch(
      RunscPagwas(single_data = 42, gwas_data = NULL, verbose = FALSE),
      error = function(e) NULL
    ),
    "deprecated"
  )
  expect_error(
    suppressWarnings(
      RunscPagwas(object = 1, single_data = 2, gwas_data = NULL, verbose = FALSE)
    ),
    "only one of"
  )
})

