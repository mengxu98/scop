test_that("C-SIDE summaries count records rather than unique genes", {
  srt <- Seurat::CreateSeuratObject(matrix(1:6, 2,
    dimnames = list(c("g1", "g2"), c("s1", "s2", "s3"))))
  tab <- data.frame(feature = c("g1", "g1", "g2"), significant = c(TRUE, TRUE, FALSE))
  meta <- cside_summary_metadata(srt, "custom", "single", tab)
  expect_identical(meta$custom_n_sig, rep(2L, 3))
  expect_identical(meta$custom_mode, rep("single", 3))
})
