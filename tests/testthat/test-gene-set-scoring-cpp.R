reference_ssgsea_scores <- function(expr, gene_sets, alpha = 0.25, normalize = FALSE) {
  n_genes <- nrow(expr)
  n_cells <- ncol(expr)
  scores <- matrix(NA_real_, n_cells, length(gene_sets))
  min_score <- Inf
  max_score <- -Inf

  for (cell in seq_len(n_cells)) {
    value_by_gene <- rep(0, n_genes)
    col <- expr[, cell]
    nonzero <- which(col != 0 & is.finite(col))
    value_by_gene[nonzero] <- col[nonzero]

    rank_by_gene <- integer(n_genes)
    rank_weight_by_gene <- numeric(n_genes)
    order_by_value <- order(value_by_gene, seq_len(n_genes), decreasing = FALSE)
    tie_start <- 1L
    while (tie_start <= n_genes) {
      tie_end <- tie_start
      tie_value <- value_by_gene[order_by_value[tie_start]]
      while (
        tie_end < n_genes &&
          value_by_gene[order_by_value[tie_end + 1L]] == tie_value
      ) {
        tie_end <- tie_end + 1L
      }
      rank_value <- as.integer((tie_start + tie_end) / 2)
      tied_genes <- order_by_value[tie_start:tie_end]
      rank_by_gene[tied_genes] <- rank_value
      rank_weight_by_gene[tied_genes] <- abs(rank_value)^alpha
      tie_start <- tie_end + 1L
    }

    ranking <- order(-rank_by_gene, seq_len(n_genes), decreasing = FALSE)
    position_by_gene <- integer(n_genes)
    position_by_gene[ranking] <- seq_len(n_genes)

    for (set_i in seq_along(gene_sets)) {
      set <- sort(unique(gene_sets[[set_i]]))
      set_size <- length(set)
      if (set_size <= 0 || set_size >= n_genes) {
        next
      }

      inverse_position <- n_genes - position_by_gene[set] + 1
      gene_weight <- rank_weight_by_gene[set]
      in_weight_sum <- sum(gene_weight)
      if (in_weight_sum <= 0) {
        next
      }

      total_position_sum <- n_genes * (n_genes + 1) / 2
      score <- sum(gene_weight * inverse_position) / in_weight_sum -
        (total_position_sum - sum(inverse_position)) / (n_genes - set_size)
      scores[cell, set_i] <- score
      if (is.finite(score)) {
        min_score <- min(min_score, score)
        max_score <- max(max_score, score)
      }
    }
  }

  if (isTRUE(normalize)) {
    scores <- scores / (max_score - min_score)
  }
  scores
}

reference_plage_scores <- function(expr, gene_sets, min_size = 1L, max_size = .Machine$integer.max) {
  n_genes <- nrow(expr)
  n_cells <- ncol(expr)
  scores <- matrix(NA_real_, n_cells, length(gene_sets))
  row_means <- rep(NA_real_, n_genes)
  row_sds <- rep(NA_real_, n_genes)
  for (gene in seq_len(n_genes)) {
    nonzero <- expr[gene, ] != 0 & is.finite(expr[gene, ])
    if (sum(nonzero) > 1L) {
      row_means[gene] <- mean(expr[gene, nonzero])
      row_sds[gene] <- stats::sd(expr[gene, nonzero])
    }
  }
  row_valid <- is.finite(row_sds) & row_sds > 0

  for (set_i in seq_along(gene_sets)) {
    set <- sort(unique(gene_sets[[set_i]]))
    set <- set[set >= 1L & set <= n_genes]
    valid <- set[row_valid[set]]
    effective_size <- length(valid)
    if (effective_size < min_size || effective_size > max_size || n_cells <= 1L) {
      next
    }

    z <- matrix(0, nrow = length(valid), ncol = n_cells)
    for (row in seq_along(valid)) {
      gene <- valid[[row]]
      nonzero <- expr[gene, ] != 0 & is.finite(expr[gene, ])
      z[row, nonzero] <- (expr[gene, nonzero] - row_means[gene]) / row_sds[gene]
    }
    first_v <- svd(z, nu = 0, nv = 1)$v[, 1]
    mean_z_by_cell <- colMeans(z)
    direction <- if (sum(first_v * mean_z_by_cell) < 0) -1 else 1
    scores[, set_i] <- direction * first_v
  }

  scores
}

test_that("ssGSEA sparse ranking matches full ranking with negative values", {
  expr <- matrix(
    c(
      -2, 0, 0, 2,
      0, 3, 0, 0,
      4, 0, 0, 0,
      0, -1, 0, 2,
      2, 0, 4, 0,
      0, 0, 3, 0
    ),
    nrow = 6,
    byrow = TRUE
  )
  gene_sets <- list(c(1L, 3L, 5L), c(2L, 4L, 6L), c(1L, 2L, 6L))

  cpp <- ssgsea_rank_dense(
    methods::as(Matrix::Matrix(expr, sparse = TRUE), "dgCMatrix"),
    gene_sets,
    alpha = 0.25,
    normalize = FALSE
  )
  ref <- reference_ssgsea_scores(expr, gene_sets, alpha = 0.25, normalize = FALSE)

  expect_equal(cpp, ref, tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("PLAGE gene-gene covariance path matches right singular vector scores", {
  set.seed(42)
  expr <- matrix(stats::rpois(40 * 320, lambda = 0.8), nrow = 40)
  expr[expr < 2] <- 0
  gene_sets <- list(1:8, 5:20, 1:35)

  cpp <- plage_dense(
    methods::as(Matrix::Matrix(expr, sparse = TRUE), "dgCMatrix"),
    gene_sets,
    min_size = 1L,
    max_size = 500L
  )
  ref <- reference_plage_scores(expr, gene_sets, min_size = 1L, max_size = 500L)

  expect_equal(cpp, ref, tolerance = 1e-12, ignore_attr = TRUE)
})
