make_estimate_mock <- function(n_samples = 6, seed = 1) {
  set.seed(seed)
  sig <- scop:::estimate_get_signatures()
  genes <- unique(c(
    sig$stromal_signature[seq_len(20)],
    sig$immune_signature[seq_len(20)],
    sig$common_genes[seq_len(80)]
  ))
  mat <- matrix(
    stats::runif(length(genes) * n_samples, 0, 20),
    nrow = length(genes),
    ncol = n_samples,
    dimnames = list(genes, paste0("Sample", seq_len(n_samples)))
  )
  mat[sig$stromal_signature[seq_len(20)], seq_len(n_samples) <= n_samples / 2] <-
    mat[sig$stromal_signature[seq_len(20)], seq_len(n_samples) <= n_samples / 2] + 20
  mat[sig$immune_signature[seq_len(20)], seq_len(n_samples) > n_samples / 2] <-
    mat[sig$immune_signature[seq_len(20)], seq_len(n_samples) > n_samples / 2] + 20
  mat
}

estimate_ssgsea_scores_reference <- function(expr, gene_sets) {
  ranked <- apply(expr, 2L, function(x) {
    rank(as.numeric(x), ties.method = "average", na.last = "keep")
  })
  ranked <- as.matrix(ranked)
  ranked[!is.finite(ranked)] <- 0
  ranked <- 10000 * ranked / nrow(ranked)
  rownames(ranked) <- rownames(expr)
  colnames(ranked) <- colnames(expr)
  scores <- matrix(NA_real_, nrow = ncol(ranked), ncol = length(gene_sets), dimnames = list(colnames(ranked), names(gene_sets)))
  for (set_i in seq_along(gene_sets)) {
    common_genes <- intersect(gene_sets[[set_i]], rownames(ranked))
    if (length(common_genes) == 0L) next
    for (sample_i in seq_len(ncol(ranked))) {
      ord <- order(ranked[, sample_i], decreasing = TRUE)
      ordered <- ranked[ord, sample_i]
      names(ordered) <- rownames(ranked)[ord]
      hit_ind <- names(ordered) %in% common_genes
      no_hit_ind <- !hit_ind
      if (!any(hit_ind) || !any(no_hit_ind)) next
      ordered_weight <- ordered^0.25
      hit_exp <- ordered_weight[hit_ind]
      if (sum(hit_exp) <= 0) next
      no_hit_penalty <- cumsum(as.numeric(no_hit_ind) / sum(no_hit_ind))
      hit_reward <- cumsum((as.numeric(hit_ind) * ordered_weight) / sum(hit_exp))
      scores[sample_i, set_i] <- sum(hit_reward - no_hit_penalty)
    }
  }
  scores
}

test_that("ESTIMATE ssGSEA precomputation preserves reference scores", {
  expr <- rbind(
    c(0, 3, 3, 1),
    c(2, 0, 3, 1),
    c(2, 4, 0, 1),
    c(0, 4, 1, 1),
    c(1, 1, 1, 1)
  )
  rownames(expr) <- paste0("Gene", seq_len(nrow(expr)))
  colnames(expr) <- paste0("Sample", seq_len(ncol(expr)))
  gene_sets <- list(first = c("Gene1", "Gene3"), second = c("Gene2", "Gene4", "missing"))

  expect_equal(
    estimate_ssgsea_scores(expr, gene_sets),
    estimate_ssgsea_scores_reference(expr, gene_sets),
    tolerance = 1e-12
  )
})

test_that("RunESTIMATE returns fixed score columns for matrix input", {
  mat <- make_estimate_mock()
  out <- RunESTIMATE(
    count_matrix = mat,
    min_sig_genes = 5,
    verbose = FALSE
  )

  expect_equal(out$status, "success")
  expect_equal(colnames(out$scores), c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity"))
  expect_equal(rownames(out$scores), colnames(mat))
  expect_equal(
    out$scores$ESTIMATEScore,
    out$scores$StromalScore + out$scores$ImmuneScore,
    tolerance = 1e-10
  )
  expect_equal(
    out$scores$TumorPurity,
    cos(0.6049872018 + 0.0001467884 * out$scores$ESTIMATEScore),
    tolerance = 1e-12
  )
  expect_true(all(is.finite(out$scores$StromalScore)))
  expect_true(all(is.finite(out$scores$ImmuneScore)))
})

test_that("RunESTIMATE fails clearly with too few signature genes", {
  mat <- matrix(
    stats::runif(20),
    nrow = 5,
    dimnames = list(paste0("Gene", 1:5), paste0("Sample", 1:4))
  )
  expect_error(
    RunESTIMATE(count_matrix = mat, min_sig_genes = 3, verbose = FALSE),
    "No genes overlap|Too few ESTIMATE signature genes"
  )
})

test_that("RunESTIMATE stores SummarizedExperiment metadata", {
  skip_if_not_installed("SummarizedExperiment")
  mat <- make_estimate_mock()
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = mat),
    colData = S4Vectors::DataFrame(
      condition = rep(c("A", "B"), each = 3),
      row.names = colnames(mat)
    )
  )
  out <- RunESTIMATE(
    object = se,
    min_sig_genes = 5,
    verbose = FALSE
  )
  store <- S4Vectors::metadata(out)[["ESTIMATE"]]
  expect_equal(store$status, "success")
  expect_equal(rownames(store$scores), colnames(mat))
})

test_that("RunESTIMATE stores Seurat pseudo-bulk results", {
  skip_if_not_installed("Seurat")
  mat <- make_estimate_mock(n_samples = 4)
  counts <- matrix(
    rep(mat, each = 2),
    nrow = nrow(mat),
    dimnames = list(rownames(mat), paste0("Cell", seq_len(ncol(mat) * 2)))
  )
  meta <- data.frame(
    sample = rep(colnames(mat), each = 2),
    condition = rep(c("A", "B"), each = 4),
    row.names = colnames(counts)
  )
  srt <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix"),
    meta.data = meta
  )
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  out <- RunESTIMATE(
    object = srt,
    sample.by = "sample",
    group.by = "condition",
    min_sig_genes = 5,
    verbose = FALSE
  )

  expect_equal(out@tools$ESTIMATE$status, "success")
  expect_equal(rownames(out@tools$ESTIMATE$scores), colnames(mat))
  expect_equal(out@tools$ESTIMATE$details$sample_metadata$group, c("A", "A", "B", "B"))
})

test_that("ESTIMATE plot helpers return ggplot objects", {
  mat <- make_estimate_mock()
  out <- RunESTIMATE(
    count_matrix = mat,
    min_sig_genes = 5,
    verbose = FALSE
  )
  groups <- setNames(rep(c("A", "B"), each = 3), rownames(out$scores))

  p1 <- EstimateScorePlot(
    object = out,
    group.data = groups,
    plot_type = "violin",
    add_stat = FALSE
  )
  p2 <- EstimateGenePlot(
    object = out,
    gene.data = mat,
    features = rownames(mat)[1],
    plot_type = "scatter"
  )
  p3 <- EstimateGenePlot(
    object = out,
    gene.data = mat,
    features = rownames(mat)[1],
    plot_type = "violin",
    add_stat = FALSE
  )

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
})
