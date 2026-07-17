test_that("RunDynamicEnrichment avoids unused raw matrices and layer retrieval", {
  counts <- Matrix::Matrix(
    matrix(
      c(2, 0, 1, 3, 1, 4),
      nrow = 2,
      dimnames = list(c("Gene1", "Gene2"), paste0("Cell", 1:3))
    ),
    sparse = TRUE
  )
  srt <- Seurat::CreateSeuratObject(counts)
  srt@tools$DynamicFeatures_L1 <- list(
    DynamicFeatures = data.frame(
      features = "Gene1",
      exp_ncells = 10,
      r.sq = 1,
      dev.expl = 1,
      padjust = 0.001,
      row.names = "Gene1"
    )
  )
  observed_features <- NULL

  testthat::local_mocked_bindings(
    CellScoring = function(srt, name, ...) {
      srt[[name]] <- SeuratObject::CreateAssayObject(counts = counts)
      srt
    },
    RunDynamicFeatures = function(srt, features, ...) {
      observed_features <<- features
      srt
    },
    GetAssayData5 = function(...) stop("RunDynamicEnrichment should not retrieve the score layer only for row names"),
    .package = "scop"
  )

  expect_no_error(
    out <- RunDynamicEnrichment(
      srt,
      lineages = "L1",
      score_method = "AUCell",
      features = list(Set1 = c("Gene1", "Gene2")),
      minGSSize = 1,
      maxGSSize = 10,
      backend = "r",
      verbose = FALSE
    )
  )
  expect_s4_class(out, "Seurat")
  expect_identical(observed_features, c("Gene1", "Gene2"))
})
