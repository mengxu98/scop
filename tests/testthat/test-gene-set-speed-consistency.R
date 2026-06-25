test_that("GSVA cpp default keeps R backend consistency on small input", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("Seurat")

  set.seed(20260619)
  counts <- Matrix::rsparsematrix(90, 24, density = 0.18)
  counts@x <- abs(round(counts@x * 5))
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  features <- list(
    set_a = paste0("gene", 1:30),
    set_b = paste0("gene", 25:60),
    set_c = paste0("gene", 50:85)
  )

  r_out <- suppressWarnings(CellScoring(
    srt,
    features = features,
    method = "GSVA",
    backend = "r",
    classification = FALSE,
    name = "gsva_r",
    verbose = FALSE
  ))
  cpp_out <- suppressWarnings(CellScoring(
    srt,
    features = features,
    method = "GSVA",
    backend = "cpp",
    classification = FALSE,
    name = "gsva_cpp",
    verbose = FALSE
  ))

  r_scores <- as.matrix(r_out@meta.data[, paste0("gsva_r_", names(features))])
  cpp_scores <- as.matrix(cpp_out@meta.data[, paste0("gsva_cpp_", names(features))])
  rho <- suppressWarnings(stats::cor(
    as.numeric(r_scores),
    as.numeric(cpp_scores),
    method = "spearman"
  ))
  expect_gte(rho, 0.95)
})

test_that("CytoTRACE2 defaults to cpp backend", {
  expect_identical(eval(formals(RunCytoTRACE.Seurat)$backend)[[1]], "cpp")
  expect_identical(eval(formals(RunCytoTRACE.default)$backend)[[1]], "cpp")
})
