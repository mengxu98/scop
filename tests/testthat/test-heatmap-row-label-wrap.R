test_that("heatmap row labels remain unchanged when wrapping is disabled", {
  labels <- stats::setNames(
    c("GO_LONG_PATHWAY_NAME", "Gene_symbol"),
    c("pathway", "gene")
  )

  expect_identical(heatmap_wrap_row_labels(labels), labels)
})

test_that("heatmap row labels wrap display text without changing identifiers", {
  labels <- stats::setNames(
    c(
      "GO_POSITIVE_REGULATION_OF_LEUKOCYTE_ACTIVATION",
      "REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX"
    ),
    c("pathway_1", "pathway_2")
  )
  original <- labels

  wrapped <- heatmap_wrap_row_labels(labels, width = 24)

  expect_identical(labels, original)
  expect_identical(names(wrapped), names(labels))
  expect_true(all(grepl("\n", wrapped, fixed = TRUE)))
  expect_false(any(grepl("_", wrapped, fixed = TRUE)))
})

test_that("heatmap row label wrapping validates its width", {
  expect_error(heatmap_wrap_row_labels("label", width = 0), "positive number")
  expect_error(heatmap_wrap_row_labels("label", width = c(10, 20)), "positive number")
  expect_error(heatmap_wrap_row_labels("label", width = "20"), "positive number")
})

test_that("wrapped row labels receive enough layout width", {
  labels <- heatmap_wrap_row_labels(
    "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
    width = 24
  )

  width <- heatmap_row_labels_max_width(labels)

  expect_s3_class(width, "unit")
  expect_gt(grid::convertWidth(width, "mm", valueOnly = TRUE), 0)
})

test_that("public heatmap functions expose optional row-name wrapping", {
  expect_null(formals(GroupHeatmap)[["row_names_wrap"]])
  expect_null(formals(FeatureHeatmap)[["row_names_wrap"]])
  expect_null(formals(NMFHeatmap)[["row_names_wrap"]])
})
