enrichment_serial_fixture <- function(n_groups) {
  genes <- paste0("G", seq_len(n_groups * 2L))
  list(
    geneID = genes,
    geneID_groups = rep(paste0("Cluster", seq_len(n_groups)), each = 2L),
    features = list(
      set1 = genes[1:3],
      set2 = genes[4:6]
    )
  )
}

test_that("RunEnrichment runs the group/database fan-out without spawning workers", {
  fixture <- enrichment_serial_fixture(6L)
  args <- list(
    geneID = fixture$geneID,
    geneID_groups = fixture$geneID_groups,
    features = fixture$features,
    minGSSize = 2,
    verbose = FALSE
  )
  reference <- do.call(RunEnrichment, args)
  expect_true(nrow(reference$enrichment) > 0)

  testthat::local_mocked_bindings(
    parallelize_fun = function(...) {
      stop("RunEnrichment must not spawn parallel workers")
    },
    .package = "scop"
  )
  serial <- do.call(RunEnrichment, args)

  expect_identical(serial$enrichment, reference$enrichment)
  expect_identical(serial$results, reference$results)
})

test_that("RunGSEA runs the group/database fan-out without spawning workers", {
  testthat::skip_if_not_installed("clusterProfiler")
  genes <- paste0("G", seq_len(6L))
  args <- list(
    geneID = rep(genes, 3L),
    geneScore = rep(c(6, 5, 4, 3, 2, 1), 3L),
    geneID_groups = rep(paste0("Cluster", seq_len(3L)), each = 6L),
    TERM2GENE = data.frame(
      Term = rep(c("Term1", "Term2"), each = 3L),
      symbol = genes
    ),
    minGSSize = 2,
    scoreType = "pos",
    verbose = FALSE
  )
  set.seed(11)
  reference <- do.call(RunGSEA, args)
  expect_true(nrow(reference$enrichment) > 0)

  testthat::local_mocked_bindings(
    parallelize_fun = function(...) {
      stop("RunGSEA must not spawn parallel workers")
    },
    .package = "scop"
  )
  set.seed(11)
  serial <- do.call(RunGSEA, args)

  # The returned `enrichResult` objects carry environments, which never compare
  # identical across two calls, so compare the tables they hold instead.
  extract_tables <- function(res) {
    lapply(res, function(x) x@result)
  }
  expect_identical(names(serial$results), names(reference$results))
  expect_identical(extract_tables(serial$results), extract_tables(reference$results))
  expect_identical(serial$enrichment, reference$enrichment)
})
