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
  z <- matrix(0, nrow = n_genes, ncol = n_cells)
  row_valid <- logical(n_genes)
  for (gene in seq_len(n_genes)) {
    values <- expr[gene, ]
    if (sum(is.finite(values)) <= 1L) {
      next
    }
    nonzero <- values[is.finite(values) & values != 0]
    if (length(nonzero) > 0L && min(nonzero) == max(nonzero)) {
      # GSVA >= 2.6 drops genes whose stored (non-zero) values are constant,
      # even when the complete row still varies because of structural zeros.
      next
    }
    scaled <- scale(values)
    if (all(is.finite(scaled))) {
      z[gene, ] <- scaled
      row_valid[gene] <- TRUE
    }
  }
  orient_mean <- rowMeans(expr)
  orient_sd <- apply(expr, 1L, stats::sd)
  orient_valid <- is.finite(orient_sd) & orient_sd > 0

  for (set_i in seq_along(gene_sets)) {
    set <- sort(unique(gene_sets[[set_i]]))
    set <- set[set >= 1L & set <= n_genes]
    valid <- set[row_valid[set]]
    effective_size <- length(valid)
    if (effective_size < min_size || effective_size > max_size || n_cells <= 1L) {
      next
    }

    first_v <- svd(z[valid, , drop = FALSE], nu = 0, nv = 1)$v[, 1]
    orient <- set[orient_valid[set]]
    ref <- colMeans((expr[orient, , drop = FALSE] - orient_mean[orient]) / orient_sd[orient])
    direction <- if (sum(first_v * ref) < 0) -1 else 1
    scores[, set_i] <- direction * first_v
  }

  scores
}

test_that("dense gene-set variable-row filtering matches the R reference", {
  expr <- rbind(
    variable = c(0, 2, NA_real_, Inf),
    constant = c(3, 3, 3, 3),
    one_finite = c(NA_real_, 1, NA_real_, Inf),
    no_finite = c(NA_real_, NaN, Inf, -Inf)
  )
  legacy_keep <- apply(expr, 1L, function(x) {
    x <- x[is.finite(x)]
    length(x) > 1L && diff(range(x)) > 0
  })

  expect_identical(
    scop:::gene_set_scoring_keep_variable_rows(expr),
    expr[legacy_keep, , drop = FALSE]
  )
})

test_that("sparse gene-set variable-row filtering matches stored-value semantics", {
  expr <- methods::as(Matrix::Matrix(
    rbind(
      variable = c(0, 2, NA_real_, Inf),
      constant = c(3, 3, 3, 3),
      one_finite = c(NA_real_, 1, NA_real_, Inf),
      no_finite = c(NA_real_, NaN, Inf, -Inf),
      structural_zero = c(0, 4, 0, 0)
    ),
    sparse = TRUE
  ), "dgCMatrix")
  legacy_keep <- vapply(split(expr@x, expr@i + 1L), function(x) {
    x <- x[is.finite(x)]
    length(x) > 1L && diff(range(x)) > 0
  }, logical(1))
  expected_keep <- rep(FALSE, nrow(expr))
  expected_keep[as.integer(names(legacy_keep))] <- legacy_keep

  expect_identical(
    scop:::gene_set_scoring_keep_variable_rows(expr),
    expr[expected_keep, , drop = FALSE]
  )
})

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

  cpp <- scop:::ssgsea_rank_dense(
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

  cpp <- scop:::plage_dense(
    methods::as(Matrix::Matrix(expr, sparse = TRUE), "dgCMatrix"),
    gene_sets,
    min_size = 1L,
    max_size = 500L,
    dense_standardize = TRUE
  )
  ref <- reference_plage_scores(expr, gene_sets, min_size = 1L, max_size = 500L)

  expect_equal(cpp, ref, tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("PLAGE sparse standardization matches GSVA", {
  skip_if_not_installed("GSVA")
  set.seed(20260711)
  expr <- matrix(stats::rpois(18 * 41, lambda = 0.9), nrow = 18)
  expr[expr < 2] <- 0
  rownames(expr) <- paste0("g", seq_len(nrow(expr)))
  colnames(expr) <- paste0("c", seq_len(ncol(expr)))
  expr <- methods::as(Matrix::Matrix(expr, sparse = TRUE), "dgCMatrix")
  gene_sets <- list(a = rownames(expr)[1:7], b = rownames(expr)[5:15])

  cpp <- t(scop:::run_plage_scores(
    expr, gene_sets, min_gs_size = 1L, max_gs_size = 50L,
    dense_standardize = scop:::gene_set_scoring_plage_dense_standardize()
  ))
  reference <- GSVA::gsva(
    GSVA::plageParam(
      exprData = expr, geneSets = gene_sets, minSize = 1L, maxSize = 50L
    ),
    verbose = FALSE
  )
  correlations <- vapply(rownames(reference), function(term) {
    stats::cor(cpp[term, ], reference[term, ])
  }, numeric(1L))

  expect_equal(unname(abs(correlations)), rep(1, length(correlations)), tolerance = 1e-12)
})

test_that("PLAGE dense standardization matches the CellScoring GSVA contract", {
  skip_if_not_installed("GSVA")
  set.seed(20260714)
  expr <- matrix(stats::rpois(20 * 47, lambda = 0.8), nrow = 20)
  expr[expr < 2] <- 0
  rownames(expr) <- paste0("g", seq_len(nrow(expr)))
  colnames(expr) <- paste0("c", seq_len(ncol(expr)))
  gene_sets <- list(a = rownames(expr)[1:8], b = rownames(expr)[6:16])

  cpp <- t(scop:::run_plage_scores(
    methods::as(Matrix::Matrix(expr, sparse = TRUE), "dgCMatrix"),
    gene_sets,
    min_gs_size = 1L,
    max_gs_size = 50L,
    dense_standardize = TRUE
  ))
  cpp <- scop:::orient_plage_scores(cpp, expr, gene_sets)
  reference <- GSVA::gsva(
    GSVA::plageParam(
      exprData = expr, geneSets = gene_sets, minSize = 1L, maxSize = 50L
    ),
    verbose = FALSE
  )
  reference <- scop:::orient_plage_scores(reference, expr, gene_sets)

  expect_identical(dim(cpp), dim(reference))
  expect_equal(as.numeric(cpp), as.numeric(reference), tolerance = 1e-10)
})

test_that("z-score sparse standardization matches the RunGSVA GSVA contract", {
  skip_if_not_installed("GSVA")
  set.seed(20260714)
  expr <- matrix(stats::rpois(20 * 47, lambda = 0.8), nrow = 20)
  expr[expr < 2] <- 0
  rownames(expr) <- paste0("g", seq_len(nrow(expr)))
  colnames(expr) <- paste0("c", seq_len(ncol(expr)))
  expr_sparse <- methods::as(Matrix::Matrix(expr, sparse = TRUE), "dgCMatrix")
  gene_sets <- list(a = rownames(expr)[1:8], b = rownames(expr)[6:16])

  cpp <- t(scop:::run_zscore_scores(
    expr_sparse,
    gene_sets,
    min_gs_size = 1L,
    max_gs_size = 50L,
    sparse_standardize = !scop:::gene_set_scoring_zscore_sparse_standardize_full(),
    sparse_standardize_full = scop:::gene_set_scoring_zscore_sparse_standardize_full()
  ))
  reference <- GSVA::gsva(
    GSVA::zscoreParam(
      exprData = expr_sparse, geneSets = gene_sets, minSize = 1L, maxSize = 50L
    ),
    verbose = FALSE
  )

  expect_identical(dim(cpp), dim(reference))
  expect_equal(as.numeric(cpp), as.numeric(reference), tolerance = 1e-10)
})

test_that("GSVA sparse delegated kernel preserves the GSVA default contract", {
  skip_if_not_installed("GSVA")
  set.seed(20260714)
  expr <- matrix(stats::rpois(20 * 31, lambda = 0.8), nrow = 20)
  expr[expr < 2] <- 0
  rownames(expr) <- paste0("g", seq_len(nrow(expr)))
  colnames(expr) <- paste0("c", seq_len(ncol(expr)))
  expr_sparse <- methods::as(Matrix::Matrix(expr, sparse = TRUE), "dgCMatrix")
  gene_sets <- list(a = rownames(expr)[1:8], b = rownames(expr)[6:16])

  delegated <- t(scop:::run_gsva_scores(
    expr_sparse,
    gene_sets,
    kcdf = "Gaussian",
    min_gs_size = 1L,
    max_gs_size = 50L,
    kernel = "delegated"
  ))
  reference <- GSVA::gsva(
    GSVA::gsvaParam(
      exprData = expr_sparse,
      geneSets = gene_sets,
      minSize = 1L,
      maxSize = 50L,
      kcdf = "Gaussian"
    ),
    verbose = FALSE
  )

  expect_identical(dim(delegated), dim(reference))
  expect_equal(as.numeric(delegated), as.numeric(reference), tolerance = 1e-10)
})

test_that("native Gaussian GSVA is an opt-in result-compatible kernel", {
  skip_if_not_installed("GSVA")
  set.seed(20260714)
  expr <- matrix(stats::rlnorm(48 * 72, meanlog = -1, sdlog = 0.7), nrow = 48)
  expr[expr < 0.5] <- 0
  rownames(expr) <- paste0("g", seq_len(nrow(expr)))
  colnames(expr) <- paste0("c", seq_len(ncol(expr)))
  expr_sparse <- methods::as(Matrix::Matrix(expr, sparse = TRUE), "dgCMatrix")
  gene_sets <- list(a = rownames(expr)[1:12], b = rownames(expr)[14:32])

  native <- scop:::run_gsva_scores(
    expr_sparse,
    gene_sets,
    kcdf = "Gaussian",
    min_gs_size = 1L,
    max_gs_size = 50L,
    sparse = FALSE,
    kernel = "native"
  )
  delegated <- scop:::run_gsva_scores(
    expr_sparse,
    gene_sets,
    kcdf = "Gaussian",
    min_gs_size = 1L,
    max_gs_size = 50L,
    sparse = FALSE,
    kernel = "delegated"
  )

  expect_identical(dim(native), dim(delegated))
  correlations <- vapply(seq_len(ncol(native)), function(i) {
    stats::cor(native[, i], delegated[, i], method = "spearman")
  }, numeric(1L))
  expect_gte(min(correlations), 0.95)
  expect_error(
    scop:::run_gsva_scores(
      expr_sparse,
      gene_sets,
      kcdf = "Poisson",
      min_gs_size = 1L,
      max_gs_size = 50L,
      kernel = "native"
    ),
    "Gaussian"
  )
})

test_that("RunGSVA selects the native Gaussian kernel through backend", {
  set.seed(20260714)
  counts <- matrix(
    stats::rpois(50 * 30, lambda = 1),
    nrow = 50,
    dimnames = list(paste0("g", seq_len(50)), paste0("c", seq_len(30)))
  )
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)

  out <- RunGSVA(
    srt,
    features = list(a = rownames(srt)[1:15], b = rownames(srt)[16:30]),
    method = "gsva",
    backend = "cpp",
    kcdf = "Gaussian",
    new_assay = FALSE,
    store_metadata = TRUE,
    verbose = FALSE
  )

  expect_true(all(c("GSVA_A", "GSVA_B") %in% colnames(out[[]])))
})
