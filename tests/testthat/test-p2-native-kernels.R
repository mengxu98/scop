test_that("native KNN voting preserves probabilities and first-level ties", {
  vote <- getFromNamespace("knn_vote_labels", "scop")
  labels <- matrix(
    c("B", "A", "A", "B", "B", "C", "A", "C", "C", "B", "A", "C"),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("q", 1:4), NULL)
  )
  levels <- c("B", "A", "C")
  out <- vote(labels, levels = levels)
  expected <- t(vapply(seq_len(nrow(labels)), function(i) {
    tabulate(match(labels[i, ], levels), nbins = length(levels)) / ncol(labels)
  }, numeric(length(levels))))
  dimnames(expected) <- list(rownames(labels), levels)

  expect_equal(out$probability, expected)
  expect_identical(
    unname(out$best),
    levels[max.col(expected, ties.method = "first")]
  )
})

test_that("native VECTOR arrows agree with the R reference", {
  vector_field <- getFromNamespace("vector_field", "scop")
  set.seed(17)
  emb <- matrix(rnorm(240), ncol = 2)
  pca <- matrix(rnorm(120 * 12), nrow = 120)
  rownames(emb) <- rownames(pca) <- paste0("cell", seq_len(nrow(emb)))

  native <- vector_field(
    emb,
    pca,
    grid.n = 9,
    arrow.p = 0.88,
    arrow.ol = 1.3,
    backend = "cpp"
  )
  reference <- vector_field(
    emb,
    pca,
    grid.n = 9,
    arrow.p = 0.88,
    arrow.ol = 1.3,
    backend = "r"
  )

  expect_equal(native$score, reference$score, tolerance = 1e-14)
  expect_equal(native$grid, reference$grid, tolerance = 1e-14)
  expect_equal(native$arrows, reference$arrows, tolerance = 1e-12)
  expect_identical(native$cell_grid, reference$cell_grid)
})

test_that("native CytoSPACE fraction estimation agrees with the R reference", {
  estimate <- getFromNamespace("cytospace_estimate_fractions", "scop")
  set.seed(23)
  st <- matrix(rpois(80 * 16, lambda = 5), nrow = 80)
  ref <- matrix(rpois(80 * 24, lambda = 4), nrow = 80)
  labels <- factor(rep(c("A", "B", "C"), each = 8), levels = c("A", "B", "C"))
  weights <- sample(1:6, ncol(st), replace = TRUE)

  native <- estimate(
    st, ref, labels, levels(labels), weights, backend = "cpp"
  )
  reference <- estimate(
    st, ref, labels, levels(labels), weights, backend = "r"
  )
  expect_equal(native, reference, tolerance = 1e-12)
})

test_that("native spatial pair aggregation agrees with aggregate", {
  count_native <- getFromNamespace("spatial_pair_count_cpp", "scop")
  set.seed(29)
  edges <- data.frame(
    sample = sample(c("s1", "s2"), 300, TRUE),
    condition = sample(c("control", "case"), 300, TRUE),
    subject = sample(c("p1", "p2", "p3"), 300, TRUE),
    from = sample(LETTERS[1:4], 300, TRUE),
    to = sample(LETTERS[1:4], 300, TRUE),
    distance = runif(300),
    stringsAsFactors = FALSE
  )
  native <- count_native(
    edges$sample, edges$condition, edges$subject, edges$from, edges$to
  )
  reference <- stats::aggregate(
    distance ~ sample + condition + subject + from + to,
    data = edges,
    FUN = length
  )
  colnames(reference)[colnames(reference) == "distance"] <- "count"
  keys <- c("sample", "condition", "subject", "from", "to")
  order_rows <- function(x) {
    out <- x[do.call(order, x[keys]), c(keys, "count"), drop = FALSE]
    rownames(out) <- NULL
    out
  }
  expect_equal(order_rows(native), order_rows(reference))
})
