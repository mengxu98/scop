test_that("GetSimilarFeatures min/max aggregation matches legacy apply", {
  correlations <- rbind(
    GeneA = c(0.2, -0.5, 0.4),
    GeneB = c(NA_real_, NA_real_, NA_real_),
    GeneC = c(0.8, 0.1, -0.2)
  )

  expect_identical(
    get_similar_features_aggregate(correlations, "min"),
    apply(correlations, 1, min)
  )
  expect_identical(
    get_similar_features_aggregate(correlations, "max"),
    apply(correlations, 1, max)
  )
})

test_that("GetSimilarFeatures preserves feature names", {
  expression <- rbind(
    Query = c(1, 2, 3, 4, 5),
    Similar = c(1, 2, 3, 4, 5),
    Different = c(5, 1, 4, 2, 3),
    Opposite = c(5, 4, 3, 2, 1)
  )
  colnames(expression) <- paste0("Cell", seq_len(ncol(expression)))
  srt <- SeuratObject::CreateSeuratObject(counts = expression)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)

  similar <- GetSimilarFeatures(
    srt,
    features = "Query",
    n = 2,
    assay = "RNA",
    layer = "data",
    verbose = FALSE
  )

  expect_type(similar, "character")
  expect_length(similar, 2)
  expect_false(anyNA(similar))
  expect_false("Query" %in% similar)
  expect_true(all(similar %in% rownames(srt)))

  similar_multiple <- GetSimilarFeatures(
    srt,
    features = c("Query", "Different"),
    n = 1,
    assay = "RNA",
    layer = "data",
    verbose = FALSE
  )

  expect_type(similar_multiple, "character")
  expect_length(similar_multiple, 1)
  expect_false(anyNA(similar_multiple))
  expect_false(similar_multiple %in% c("Query", "Different"))
  expect_true(similar_multiple %in% rownames(srt))
})
