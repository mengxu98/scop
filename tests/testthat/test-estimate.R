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
