legacy_feature_cor_geometric_mean <- function(x, log_normalized = FALSE) {
  if (isTRUE(log_normalized)) {
    return(apply(expm1(x), 2, function(values) log1p(exp(mean(log(values))))))
  }
  apply(x, 2, function(values) exp(mean(log(values))))
}

test_that("FeatureCorPlot sparse geometric mean matches legacy preprocessing", {
  counts <- Matrix::Matrix(
    matrix(
      c(1, 2, 0, 4, 5, 6, 7, 0, 9, 10, 11, 12),
      nrow = 3,
      dimnames = list(paste0("gene", 1:3), paste0("cell", 1:4))
    ),
    sparse = TRUE
  )

  expect_equal(
    feature_cor_geometric_mean(counts),
    legacy_feature_cor_geometric_mean(counts)
  )
  log_counts <- log1p(counts)
  expect_equal(
    feature_cor_geometric_mean(log_counts, log_normalized = TRUE),
    legacy_feature_cor_geometric_mean(log_counts, log_normalized = TRUE)
  )
})

test_that("FeatureCorPlot geometric mean retains legacy fallback for invalid values", {
  values <- matrix(c(1, -1, NA_real_, 3, 2, 4), nrow = 2)
  expect_equal(
    suppressWarnings(feature_cor_geometric_mean(values)),
    suppressWarnings(legacy_feature_cor_geometric_mean(values))
  )
})

test_that("FeatureCorPlot calculates co-expression on a public plot call", {
  counts <- Matrix::Matrix(
    matrix(
      c(1, 2, 3, 4, 4, 3, 2, 1, 1, 3, 2, 4),
      nrow = 3, byrow = TRUE,
      dimnames = list(paste0("gene", 1:3), paste0("cell", 1:4))
    ),
    sparse = TRUE
  )
  srt <- Seurat::CreateSeuratObject(counts)
  srt <- NormalizeData(srt, verbose = FALSE)
  srt$group <- factor(c("a", "a", "b", "b"))

  expect_no_error(
    FeatureCorPlot(
      srt,
      features = rownames(srt)[1:2],
      group.by = "group",
      calculate_coexp = TRUE,
      add_smooth = FALSE,
      add_r2 = FALSE,
      add_pvalue = FALSE,
      verbose = FALSE
    )
  )
})
