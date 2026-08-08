# End-to-end consistency between the scop CCC wrappers and the original
# backend methods on real data (pancreas_sub). No mocks.

real_ccc_input <- function(n_cells = 300, seed = 42) {
  data(pancreas_sub)
  set.seed(seed)
  srt <- Seurat::NormalizeData(pancreas_sub, verbose = FALSE)
  srt <- srt[, sample(colnames(srt), n_cells)]
  srt$CellType <- factor(srt$CellType)
  srt
}

test_that("RunLIANA results match the original liana pipeline", {
  skip_on_cran()
  skip_if_not_installed("liana")

  srt <- real_ccc_input()

  wrapped <- RunLIANA(
    srt,
    group.by = "CellType",
    method = "logfc",
    species = "mouse",
    backend = "cpp",
    verbose = FALSE
  )
  long_table <- wrapped@tools$LIANA$long_table
  expect_true(nrow(long_table) > 0)
  expect_true(all(long_table$method == "LIANA"))

  # original pipeline with identical inputs
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = Seurat::GetAssayData(srt, assay = "RNA", layer = "counts"),
      logcounts = Seurat::GetAssayData(srt, assay = "RNA", layer = "data")
    ),
    colData = srt[[]]
  )
  original <- liana::liana_wrap(
    sce = sce,
    method = "logfc",
    resource = "MouseConsensus",
    idents_col = "CellType",
    assay = NULL,
    min_cells = 5,
    return_all = FALSE,
    verbose = FALSE,
    base = exp(1)
  )
  original_df <- as.data.frame(original)
  expect_true(nrow(original_df) > 0)

  long_key <- paste(long_table$sender, long_table$receiver,
    long_table$ligand.complex, long_table$receptor.complex)
  original_key <- paste(original_df$source, original_df$target,
    original_df$ligand.complex, original_df$receptor.complex)
  expect_gt(length(intersect(long_key, original_key)) / length(unique(original_key)), 0.9)

  # scores agree for shared interactions (logfc_comb is the logfc method score)
  merged <- merge(
    data.frame(key = original_key, score_original = original_df$logfc_comb),
    data.frame(key = long_key, score_wrapped = long_table$score),
    by = "key"
  )
  expect_gt(nrow(merged), 0)
  expect_equal(merged$score_wrapped, merged$score_original, tolerance = 1e-10)
  expect_identical(wrapped@tools$LIANA$parameters$backend, "cpp")
})

test_that("RunLIANA cpp and r backends produce identical tables", {
  skip_on_cran()
  skip_if_not_installed("liana")

  srt <- real_ccc_input(n_cells = 200)

  cpp_out <- RunLIANA(srt, group.by = "CellType", method = "logfc",
    species = "mouse", backend = "cpp", verbose = FALSE)
  r_out <- RunLIANA(srt, group.by = "CellType", method = "logfc",
    species = "mouse", backend = "r", verbose = FALSE)

  lt_cpp <- cpp_out@tools$LIANA$long_table
  lt_r <- r_out@tools$LIANA$long_table
  expect_identical(dim(lt_cpp), dim(lt_r))
  expect_equal(lt_cpp$score, lt_r$score, tolerance = 1e-12)
  expect_equal(lt_cpp$pvalue, lt_r$pvalue, tolerance = 1e-12)

  pt_cpp <- cpp_out@tools$LIANA$pair_table
  pt_r <- r_out@tools$LIANA$pair_table
  pt_cpp <- pt_cpp[order(pt_cpp$sender, pt_cpp$receiver), , drop = FALSE]
  pt_r <- pt_r[order(pt_r$sender, pt_r$receiver), , drop = FALSE]
  rownames(pt_cpp) <- NULL
  rownames(pt_r) <- NULL
  expect_equal(pt_cpp, pt_r, tolerance = 1e-12)
})
