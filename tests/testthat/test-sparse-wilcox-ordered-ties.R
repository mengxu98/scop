test_that("sparse Wilcoxon keeps zero, signed values and ties in rank order", {
  native <- getFromNamespace("wilcox_rank_sum_sparse", "scop")
  values <- rbind(
    mixed = c(-4, -4, -1, 0, 0, -0, 0, 1, 1, 1),
    negative = c(-3, -3, -1, -1, 0, -2, -2, -2, 0, -2),
    invalid = c(0, 0, NA, NaN, -Inf, Inf, 0, 1, 2, 2),
    unsigned = c(0, 1, 0, 2, 0, 3, 0, 2, 0, 1),
    boundary = c(0.2, 0.2 + .Machine$double.eps, 0.2, 0.2, 0,
                 0.2, 0.2 - .Machine$double.eps, 0.2, 0.2, 1),
    all_zero = rep(0, 10)
  )
  positions <- which(is.na(values) | values != 0, arr.ind = TRUE)
  # Deliberately include two stored zero entries adjacent to implicit zeros.
  positions <- rbind(positions, c(1L, 4L), c(1L, 6L))
  x <- Matrix::sparseMatrix(
    i = positions[, 1L], j = positions[, 2L], x = values[positions],
    dims = dim(values), dimnames = dimnames(values)
  )
  expect_true(any(x@x == 0, na.rm = TRUE))

  for (cut in c(-Inf, -2, 0, 0.2, Inf)) {
    effective <- values
    effective[!is.finite(effective) | effective <= cut] <- 0
    for (n1 in c(1L, 5L, 7L)) {
      reference <- apply(effective, 1L, function(row) {
        if (length(unique(row)) == 1L) return(1)
        left <- seq_len(n1)
        stats::wilcox.test(row[left], row[-left], exact = FALSE)$p.value
      })
      expect_equal(native(x, n1, cut, 1L), unname(reference), tolerance = 1e-15)
    }
  }
})
