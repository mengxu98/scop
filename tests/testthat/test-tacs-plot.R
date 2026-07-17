test_that("GetSimilarFeatures min/max aggregation matches legacy apply", {
  correlations <- rbind(
    GeneA = c(0.2, -0.5, 0.4),
    GeneB = c(NA_real_, NA_real_, NA_real_),
    GeneC = c(0.8, 0.1, -0.2)
  )

  expect_identical(
    scop:::get_similar_features_aggregate(correlations, "min"),
    apply(correlations, 1, min)
  )
  expect_identical(
    scop:::get_similar_features_aggregate(correlations, "max"),
    apply(correlations, 1, max)
  )
})
