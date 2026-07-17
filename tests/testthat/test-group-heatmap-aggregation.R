group_heatmap_legacy_aggregate <- function(mat, groups, fun) {
  out <- Matrix::t(stats::aggregate(Matrix::t(mat), by = list(groups), FUN = fun))
  colnames(out) <- out[1, , drop = FALSE]
  out <- out[-1, , drop = FALSE]
  class(out) <- "numeric"
  out
}

test_that("GroupHeatmap default aggregation and expression fractions match aggregate", {
  mat <- matrix(
    c(1, 0, 2, 4, 3, 1, 5, 2, 4, 6, 0, 8),
    nrow = 3,
    dimnames = list(paste0("Gene", 1:3), paste0("Cell", 1:4))
  )
  groups <- factor(c("B", "A", "B", "A"), levels = c("A", "B", "unused"))
  names(groups) <- colnames(mat)
  legacy_mean <- group_heatmap_legacy_aggregate(mat, groups, base::mean)
  legacy_fraction <- group_heatmap_legacy_aggregate(
    mat,
    groups,
    function(x) sum(x > 2) / length(x)
  )

  expect_identical(group_heatmap_aggregate_matrix(mat, groups), legacy_mean)
  expect_identical(group_heatmap_expression_fraction(mat, groups, exp_cutoff = 2), legacy_fraction)
})

test_that("GroupHeatmap custom aggregation retains the aggregate fallback", {
  mat <- matrix(
    c(1, 0, 2, 4, 3, 1, 5, 2, 4, 6, 0, 8),
    nrow = 3,
    dimnames = list(paste0("Gene", 1:3), paste0("Cell", 1:4))
  )
  groups <- factor(c("B", "A", "B", "A"), levels = c("A", "B"))
  names(groups) <- colnames(mat)
  expect_identical(
    group_heatmap_aggregate_matrix(mat, groups, stats::median),
    group_heatmap_legacy_aggregate(mat, groups, stats::median)
  )
})
