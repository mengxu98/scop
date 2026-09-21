# Backend parity: the cpp aggregation and the reference R aggregation of
# observed label pairs must return identical pair tables, including row order
# (deterministic ordering added in spatial_neighborhood_observed_pairs).
# The internals are exercised directly; no spatial backend package is needed.

spatial_parity_cells <- function(n = 60L, seed = 11L) {
  set.seed(seed)
  data.frame(
    cell = paste0("c", seq_len(n)),
    x = stats::runif(n, 0, 10),
    y = stats::runif(n, 0, 10),
    group = rep(c("A", "B", "C"), length.out = n),
    sample = "s1",
    condition = "all",
    subject = "s1",
    stringsAsFactors = FALSE
  )
}

test_that("spatial neighborhood cpp and r backends return identical pair tables", {
  cells <- spatial_parity_cells()
  observed <- getFromNamespace("spatial_neighborhood_observed_pairs", "scop")
  cpp <- observed(cells = cells, k = 6L, backend = "cpp")
  ref <- observed(cells = cells, k = 6L, backend = "r")

  expect_identical(cpp$edge_table, ref$edge_table)
  expect_identical(cpp$pair_table, ref$pair_table)
  expect_identical(rownames(cpp$pair_table), rownames(ref$pair_table))
})

test_that("spatial neighborhood pair table rows are deterministically ordered", {
  cells <- spatial_parity_cells()
  observed <- getFromNamespace("spatial_neighborhood_observed_pairs", "scop")
  ref <- observed(cells = cells, k = 6L, backend = "r")

  pt <- ref$pair_table
  expect_identical(
    pt,
    pt[order(pt$sample, pt$condition, pt$subject, pt$from, pt$to), , drop = FALSE]
  )
  expect_identical(rownames(pt), as.character(seq_len(nrow(pt))))
})
