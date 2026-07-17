test_that("EnrichmentPlot overlap edges match pairwise shared-gene counts", {
  actual <- enrichment_overlap_edges(
    c("term_a", "term_b", "term_c"),
    list(c("A", "B", "B"), c("B", "C"), "D")
  )
  expected <- data.frame(
    from = "term_a",
    to = "term_b",
    weight = 1,
    stringsAsFactors = FALSE
  )
  expect_equal(actual, expected)
})

test_that("EnrichmentPlot overlap edges omit disjoint and singleton term sets", {
  actual <- enrichment_overlap_edges(c("term_a", "term_b"), list("A", "B"))
  expect_equal(actual, data.frame(
    from = character(), to = character(), weight = numeric()
  ))
})
