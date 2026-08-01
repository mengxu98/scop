test_that("native UCell scores match official tie and missing-gene semantics", {
  skip_if_not_installed("UCell")
  set.seed(20260730)
  expr <- matrix(
    sample(0:5, 20 * 12, replace = TRUE),
    nrow = 20,
    dimnames = list(paste0("Gene", 1:20), paste0("Cell", 1:12))
  )
  expr <- methods::as(Matrix::Matrix(expr, sparse = TRUE), "dgCMatrix")
  features <- list(
    Positive = c("Gene1", "Gene2", "Gene3", "Missing"),
    Signed = c("Gene4+", "Gene5", "Gene10-", "Gene11-")
  )
  run_ucell <- getFromNamespace("run_ucell_scores", "scop")

  for (missing_genes in c("impute", "skip")) {
    for (ties_method in c("average", "min", "max", "dense", "first", "last")) {
      expected <- UCell::ScoreSignatures_UCell(
        matrix = expr,
        features = features,
        maxRank = 15,
        missing_genes = missing_genes,
        ties.method = ties_method,
        ncores = 1,
        BPPARAM = BiocParallel::SerialParam(),
        name = ""
      )
      observed <- run_ucell(
        expr,
        features,
        max_rank = 15,
        missing_genes = missing_genes,
        ties_method = ties_method
      )
      expect_equal(
        unname(observed),
        unname(expected),
        tolerance = 1e-12,
        info = paste(missing_genes, ties_method)
      )
    }
  }
})

test_that("CellScoring UCell backend agrees end to end", {
  skip_if_not_installed("UCell")
  set.seed(20260731)
  counts <- matrix(
    sample(0:8, 30 * 20, replace = TRUE),
    nrow = 30,
    dimnames = list(paste0("Gene", 1:30), paste0("Cell", 1:20))
  )
  srt <- Seurat::CreateSeuratObject(
    methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  )
  features <- list(
    State1 = paste0("Gene", 1:8),
    State2 = c(paste0("Gene", 9:14), "Gene15-", "Gene16-")
  )
  reference <- CellScoring(
    srt,
    features = features,
    layer = "counts",
    method = "UCell",
    backend = "r",
    classification = FALSE,
    name = "UC",
    maxRank = 25,
    BPPARAM = BiocParallel::SerialParam(),
    verbose = FALSE
  )
  candidate <- CellScoring(
    srt,
    features = features,
    layer = "counts",
    method = "UCell",
    backend = "cpp",
    classification = FALSE,
    name = "UC",
    maxRank = 25,
    verbose = FALSE
  )
  score_columns <- c("UC_State1", "UC_State2")
  expect_equal(
    as.matrix(candidate[[]][, score_columns]),
    as.matrix(reference[[]][, score_columns]),
    tolerance = 1e-12
  )
})
