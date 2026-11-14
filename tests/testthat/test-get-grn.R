test_that("GetGRN preserves negative and unknown scores without changing native", {
  native <- data.frame(regulator = c("01", "1"), target = c("G", "G"), weight = c(-2, 3), beta = c(-0.4, 0.8), source_type = c("original", "original"))
  entry <- list(
    cells = c("A", "B"), native = list(edges = native),
    layer_specs = list(tf_gene = list(status = "computed", table_path = "edges", source = "regulator", target = "target", score = "weight", source_type = "gene", target_type = "gene", score_type = "correlation", directed = FALSE, sign_column = "beta", sign_basis = "correlation", rank_by = "absolute_score")),
    execution = list(status = "completed")
  )
  result <- structure(list(schema_version = "1.0", groups = list(all = entry)), class = c("GRNResult", "list"))
  before <- serialize(result, NULL)
  edges <- GetGRN(result)
  expect_identical(edges$source, c("01", "1"))
  expect_identical(edges$raw_score, c(-2, 3))
  expect_identical(edges$sign, c(-1L, 1L))
  expect_false(any(edges$directed))
  expect_identical(edges$native_source_type, native$source_type)
  expect_identical(GetGRN(result, what = "native")$edges, native)
  expect_identical(GetGRN(result, what = "gene_sets"), list())
  expect_identical(before, serialize(result, NULL))
  result$groups$all$layer_specs$tf_gene$directed <- TRUE
  result$groups$all$layer_specs$tf_gene$source_type <- "TF"
  expect_named(GetGRN(result, what = "gene_sets"), c("01__negative", "1__positive"))
  result$groups$all$layer_specs$tf_gene$rows <- 2L
  expect_identical(GetGRN(result)$source, "1")
  expect_identical(GetGRN(result, what = "native")$edges, native)
  result$groups$all$layer_specs$tf_gene$rows <- 3L
  expect_error(GetGRN(result), "row selection")
  result$groups$all$layer_specs$tf_gene$rows <- NULL
  result$groups$all$layer_specs$tf_gene$sign_column <- NULL
  expect_true(all(is.na(GetGRN(result)$sign)))
  result$groups$all$layer_specs$tf_gene$status <- "not_computed"
  expect_error(GetGRN(result), "not_computed")
})
