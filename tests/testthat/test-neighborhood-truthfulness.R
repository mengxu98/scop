test_that("neighborhood counts preserve saved scope, zero and unknown values", {
  srt <- Seurat::CreateSeuratObject(matrix(1:12, 2,
    dimnames = list(c("g1", "g2"), paste0("s", 1:6))))
  srt$x <- c(0, 1, 20, 30, 31, 90)
  srt$y <- rep(0, 6)
  srt$label <- c("A", "B", "A", "A", "B", NA)
  srt$condition <- c("one", "one", "one", "two", "two", "one")
  srt$sample <- c("x", "x", "x", "y", "y", "x")
  run <- function(radius = 2, ...) RunSpatialNeighborhood(srt, group.by = "label",
    method = "observed", sample.by = "sample", split.by = "condition",
    coord.cols = c("x", "y"), radius = radius, verbose = FALSE, ...)
  out <- run(backend = "r")
  local_mocked_bindings(SpatialSpotPlot = function(srt, values, ...) values, .package = "scop")
  plot <- function(object = out, ...) SpatialNeighborhoodPlot(object, plot_type = "spatial", ...)
  expected <- stats::setNames(c(1, NA, 0, 1, NA, NA), colnames(srt))
  expect_identical(plot(pair = c("A", "B")), expected)
  expect_identical(plot(pair = "A|B", value = "count"), expected)
  expected["s4"] <- NA_real_
  expect_identical(plot(pair = "A|B", condition = "one"), expected)
  expect_error(plot(pair = "A|B", condition = "missing"), "condition")
  expect_error(plot(pair = "typo|B"), "known labels")
  expect_error(plot(value = "estimate"), "only value")
  expect_error(plot(value = "fraction"), "only value")
  expect_error(plot(comparison = "one"), "comparison")
  sub <- out[, c("s3", "s1")]
  expect_identical(plot(sub, pair = "A|B"), c(s3 = 0, s1 = 1)[colnames(sub)])
  out$condition <- "changed"
  expect_identical(plot(pair = "A|B", condition = "one"), expected)
  empty <- run(radius = 0.01, backend = "r")
  expect_equal(empty@tools$SpatialNeighborhood$summary$n_edges, 0)
  expect_identical(plot(empty, pair = "A|B"), stats::setNames(c(0, NA, 0, 0, NA, NA), colnames(srt)))
  expect_error(plot(empty), "specify pair")
  expect_s3_class(SpatialNeighborhoodPlot(empty, plot_type = "stat"), "ggplot")
  limited <- run(from = "A", to = "A", backend = "r")
  expect_error(plot(limited, pair = "A|B"), "scope")
  expect_equal(limited@tools$SpatialNeighborhood$summary$n_edges, 0)
  expect_error(run(from = "typo", backend = "r"), "from/to")
  empty@tools$SpatialNeighborhood$methods$observed$input <- NULL
  expect_error(plot(empty, pair = "A|B"), "rerun")
})

test_that("empty native neighborhoods render and aggregation backends agree", {
  srt <- Seurat::CreateSeuratObject(Matrix::Matrix(matrix(1:8, 2,
    dimnames = list(c("g1", "g2"), paste0("s", 1:4))), sparse = TRUE))
  srt$x <- c(0, 1, 0, 1)
  srt$y <- rep(0, 4)
  srt$label <- c("A", "B", "A", "B")
  srt$sample <- c("x", "x", "y", "y")
  run <- function(backend, radius) RunSpatialNeighborhood(srt, group.by = "label",
    sample.by = "sample", coord.cols = c("x", "y"), radius = radius,
    backend = backend, verbose = FALSE)
  for (radius in c(0.01, 2)) {
    r <- run("r", radius)
    cpp <- run("cpp", radius)
    rb <- r@tools$SpatialNeighborhood$methods$observed
    cb <- cpp@tools$SpatialNeighborhood$methods$observed
    canonical <- function(x) {
      x <- x[order(x$sample, x$condition, x$from, x$to), , drop = FALSE]
      rownames(x) <- NULL
      x
    }
    expect_equal(canonical(rb$pair_table), canonical(cb$pair_table))
    expect_equal(rb$edge_table, cb$edge_table)
    expect_equal(nrow(cb$edge_table), if (radius < 1) 0 else 4)
    expect_s3_class(SpatialNeighborhoodPlot(cpp, plot_type = "spatial", pair = "A|B",
      coord.cols = c("x", "y"), overlay_image = FALSE, theme_use = NULL), "ggplot")
  }
})
