test_that("AUCell consistency strategy matches official AUCell scores", {
  skip_if_not_installed("AUCell")

  set.seed(42)
  expr <- Matrix::rsparsematrix(120, 24, density = 0.08)
  expr@x <- abs(round(expr@x * 5))
  rownames(expr) <- paste0("gene", seq_len(nrow(expr)))
  colnames(expr) <- paste0("cell", seq_len(ncol(expr)))
  gene_sets <- list(
    set_a = paste0("gene", 1:30),
    set_b = paste0("gene", 25:75),
    set_c = paste0("gene", 70:115)
  )

  set.seed(11)
  expected <- suppressWarnings(
    getFromNamespace("run_aucell_official_scores", "scop")(
      expr_counts = expr,
      gene_sets = gene_sets
    )
  )
  set.seed(11)
  observed <- suppressWarnings(
    getFromNamespace("run_aucell_official_scores", "scop")(
      expr_counts = expr,
      gene_sets = gene_sets
    )
  )

  expect_equal(observed, expected)
})

test_that("CellScoring AUCell backend switch controls R and C++ paths", {
  skip_if_not_installed("AUCell")
  skip_if_not_installed("Seurat")

  set.seed(43)
  counts <- Matrix::rsparsematrix(90, 18, density = 0.1)
  counts@x <- abs(round(counts@x * 4))
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  features <- list(
    set_a = paste0("gene", 1:25),
    set_b = paste0("gene", 20:55)
  )

  r_out <- suppressWarnings(CellScoring(
    srt,
    features = features,
    method = "AUCell",
    backend = "r",
    classification = FALSE,
    name = "auc_r",
    verbose = FALSE
  ))
  cpp_out <- suppressWarnings(CellScoring(
    srt,
    features = features,
    method = "AUCell",
    backend = "cpp",
    classification = FALSE,
    name = "auc_cpp",
    verbose = FALSE
  ))

  r_scores <- as.matrix(r_out@meta.data[, c("auc_r_set_a", "auc_r_set_b")])
  cpp_scores <- as.matrix(cpp_out@meta.data[, c("auc_cpp_set_a", "auc_cpp_set_b")])
  colnames(r_scores) <- colnames(cpp_scores)

  expect_equal(dim(cpp_scores), dim(r_scores))
  expect_true(all(is.finite(cpp_scores)))
})

test_that("CellScoring AUCell cpp backend keeps high consistency with R backend", {
  skip_if_not_installed("AUCell")
  skip_if_not_installed("Seurat")

  set.seed(44)
  counts <- Matrix::rsparsematrix(120, 24, density = 0.08)
  counts@x <- abs(round(counts@x * 4))
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  features <- list(
    set_a = paste0("gene", 1:35),
    set_b = paste0("gene", 25:85)
  )

  r_out <- suppressWarnings(CellScoring(
    srt,
    features = features,
    method = "AUCell",
    backend = "r",
    classification = FALSE,
    name = "auc_r",
    verbose = FALSE
  ))
  cpp_out <- suppressWarnings(CellScoring(
    srt,
    features = features,
    method = "AUCell",
    backend = "cpp",
    classification = FALSE,
    name = "auc_cpp",
    verbose = FALSE
  ))

  r_scores <- as.matrix(r_out@meta.data[, c("auc_r_set_a", "auc_r_set_b")])
  cpp_scores <- as.matrix(cpp_out@meta.data[, c("auc_cpp_set_a", "auc_cpp_set_b")])
  expect_gte(
    suppressWarnings(stats::cor(
      as.numeric(r_scores),
      as.numeric(cpp_scores),
      method = "spearman",
      use = "complete.obs"
    )),
    0.95
  )
})

test_that("RunMetabolism exposes backend without AUCell strategy parameter", {
  expect_true("backend" %in% names(formals(RunMetabolism)))
  expect_false(any(grepl("strategy", names(formals(RunMetabolism)), fixed = TRUE)))
})
