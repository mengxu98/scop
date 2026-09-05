profile_fixture <- function(n = 6L) {
  counts <- matrix(1L, 2, n, dimnames = list(c("g1", "g2"), paste0("c", seq_len(n))))
  input <- if (n == 1L) SeuratObject::CreateAssayObject(counts) else counts
  object <- suppressWarnings(SeuratObject::CreateSeuratObject(input))
  object$col <- seq_len(n) - 1
  object$row <- rep(0, n)
  object$label <- rep(c("T", "M"), length.out = n)
  object$sample <- rep("s1", n)
  object
}

profile_quiet <- function(...) SpatialNeighborhoodProfile(..., verbose = FALSE)

profile_oracle <- function(object, result, cumulative = FALSE) {
  meta <- object[[]]
  vapply(seq_len(nrow(result)), function(i) {
    q <- match(result$cell_id[i], rownames(meta))
    d <- sqrt((meta$col - meta$col[q])^2 + (meta$row - meta$row[q])^2)
    sum(seq_len(nrow(meta)) != q & meta$sample == result$sample[i] &
          meta$label == result$group[i] & d <= result$radius[i] &
          (cumulative | result$lower[i] == 0 | d > result$lower[i]))
  }, integer(1))
}

test_that("profiles match brute force and preserve context, identity and ordering", {
  object <- profile_fixture()
  object$col <- c(0, 0, 3, -3, 100, 0)
  object$row <- c(0, 0, 4, -4, 0, 0)
  object$sample <- c(rep("s1", 5), "s2")
  before <- serialize(object, NULL)
  targets <- c("c4", "c1", "c6")
  out <- profile_quiet(object, "label", c(1, 5), cells = targets, sample.by = "sample")
  expect_identical(out$count, profile_oracle(object, out))
  expect_identical(unique(out$cell_id), targets)
  expect_identical(serialize(object, NULL), before)
  expect_identical(attr(out, "source")$coordinate_space, "raw")
  expect_equal(attr(out, "parameters")$context_n, 6)
  expect_true(all(out$count[out$cell_id == "c6"] == 0))
  expect_true(all(is.na(out$fraction[out$total == 0])))
  expect_equal(out$count[out$cell_id == "c1" & out$group == "M" & out$radius == 1], 1L)
  full <- profile_quiet(object, "label", c(1, 5), sample.by = "sample")
  key <- function(x) paste(x$cell_id, x$group, x$radius)
  expect_identical(out$count, full$count[match(key(out), key(full))])
  shuffled <- object[, rev(colnames(object))]
  expect_identical(profile_quiet(shuffled, "label", c(1, 5), targets, "sample")$count, out$count)
  cumulative <- profile_quiet(object, "label", c(1, 5), targets, "sample", cumulative = TRUE)
  expect_identical(cumulative$count, profile_oracle(object, cumulative, TRUE))
  valid <- out[out$total > 0, ]
  expect_true(all(abs(aggregate(fraction ~ cell_id + radius, valid, sum)$fraction - 1) < 1e-12))
})

test_that("random coordinates and multiple samples agree with a separate oracle", {
  for (seed in 1:10) {
    set.seed(seed)
    object <- profile_fixture(100)
    object$col <- runif(100, -10, 10)
    object$row <- runif(100, -10, 10)
    object$sample <- sample(c("s1", "s2"), 100, TRUE)
    object$label <- sample(LETTERS[1:4], 100, TRUE)
    targets <- sample(colnames(object), 12)
    radii <- sort(runif(4, .1, 5))
    out <- profile_quiet(object, "label", radii, targets, "sample")
    expect_identical(out$count, profile_oracle(object, out))
    cumulative <- profile_quiet(object, "label", radii, targets, "sample", cumulative = TRUE)
    expect_identical(cumulative$count, profile_oracle(object, cumulative, TRUE))
    for (r in radii) {
      single <- profile_quiet(object, "label", r, targets, "sample")
      expect_identical(cumulative$count[cumulative$radius == r], single$count)
    }
  }
})

test_that("empty targets and single-cell contexts have stable table semantics", {
  object <- profile_fixture(1)
  one <- profile_quiet(object, "label", 1)
  expect_equal(nrow(one), 1L)
  expect_identical(one$count, 0L)
  expect_equal(one$total, 0)
  expect_true(is.na(one$fraction))
  empty <- profile_quiet(object, "label", c(1, 2), character())
  expect_s3_class(empty, "data.frame")
  expect_equal(nrow(empty), 0L)
  expect_identical(names(empty), names(one))
  expect_type(empty$count, "integer")
  expect_type(empty$fraction, "double")
  expect_equal(attr(empty, "parameters")$context_n, 1L)
})

test_that("invalid input errors rather than returning a partial profile", {
  object <- profile_fixture()
  before <- serialize(object, NULL)
  expect_error(profile_quiet(NULL, "label", 1), "Seurat")
  expect_error(profile_quiet(object, "missing", 1), "metadata")
  expect_error(profile_quiet(object, "label", 1, sample.by = "missing"), "metadata")
  expect_error(profile_quiet(object, "label", c(1, 1)), "increasing")
  expect_error(profile_quiet(object, "label", c(2, 1)), "increasing")
  expect_error(profile_quiet(object, "label", 0), "positive")
  expect_error(profile_quiet(object, "label", Inf), "finite")
  expect_error(profile_quiet(object, "label", 1e-200), "magnitude")
  expect_error(profile_quiet(object, "label", 1, cumulative = NA), "logical")
  expect_error(profile_quiet(object, "label", 1, c("c1", "c1")), "unique")
  expect_error(profile_quiet(object, "label", 1, "missing"), "resolved context")
  expect_identical(serialize(object, NULL), before)
  object$label[2] <- NA_character_
  expect_error(profile_quiet(object, "label", 1, "c1"), "nonmissing")
  object$label[2] <- "M"; object$sample[2] <- ""
  expect_error(profile_quiet(object, "label", 1, "c1", "sample"), "nonempty")
  object$col[2] <- NA_real_
  expect_error(profile_quiet(object, "label", 1), "non-finite")
})

test_that("multiple images require selection and restrict both targets and context", {
  object <- profile_fixture(4)
  assay <- SeuratObject::DefaultAssay(object)
  object[["slice1"]] <- SeuratObject::CreateFOV(
    data.frame(x = c(0, 1), y = c(0, 0), row.names = c("c1", "c2")),
    type = "centroids", assay = assay, key = "s1_"
  )
  object[["slice2"]] <- SeuratObject::CreateFOV(
    data.frame(x = c(0, 10), y = c(0, 0), row.names = c("c3", "c4")),
    type = "centroids", assay = assay, key = "s2_"
  )
  expect_error(profile_quiet(object, "label", 5), "Multiple spatial images")
  out <- profile_quiet(object, "label", c(5, 10), image = "slice2")
  expect_identical(unique(out$cell_id), c("c3", "c4"))
  expect_identical(attr(out, "source")$image, "slice2")
  expect_equal(attr(out, "parameters")$context_n, 2)
  expect_true(all(out$count[out$radius == 5] == 0))
  expect_equal(sum(out$count[out$radius == 10]), 2L)
  expect_error(profile_quiet(object, "label", 10, "c1", image = "slice2"), "resolved context")
})

test_that("the thin receipt describes returned, not saved, results", {
  recorded <- NULL
  local_mocked_bindings(spatial_run_receipt = function(...) recorded <<- list(...))
  out <- SpatialNeighborhoodProfile(profile_fixture(), "label", c(1, 3), cells = "c1")
  expect_s3_class(out, "data.frame")
  expect_match(recorded$scope, "1 target cells; 6 context cells")
  expect_match(recorded$inspect, "Returned data.frame")
  expect_null(recorded$saved)
})
