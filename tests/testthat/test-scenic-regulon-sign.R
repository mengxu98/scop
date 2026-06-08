test_that("C++ SCENIC module builder labels positive regulons by default", {
  expr <- matrix(
    c(
      1, 2, 3, 4, 5, 6,
      1, 2, 3, 4, 5, 6,
      6, 5, 4, 3, 2, 1,
      1, 1, 2, 2, 3, 3
    ),
    nrow = 6,
    ncol = 4,
    dimnames = list(paste0("Cell", 1:6), c("TF1", "GenePos", "GeneNeg", "GeneWeak"))
  )
  adjacency <- data.frame(
    TF = "TF1",
    target = c("GenePos", "GeneNeg", "GeneWeak"),
    importance = c(0.9, 0.8, 0.7),
    stringsAsFactors = FALSE
  )

  modules <- getFromNamespace("scenic_modules_from_adjacencies", "scop")(
    adjacency = adjacency,
    expr_mtx = expr,
    thresholds = 0.5,
    top_n_targets = 3,
    top_n_regulators = integer(0),
    min_genes = 2,
    keep_only_activating = TRUE,
    verbose = FALSE
  )

  expect_true(length(modules) > 0)
  expect_true(all(vapply(modules, `[[`, integer(1), "regulation") == 1L))
  expect_true(all(vapply(modules, `[[`, character(1), "suffix") == "(+)"))
  expect_true(all(vapply(modules, function(module) "GeneNeg" %in% module$genes, logical(1)) == FALSE))
})

test_that("C++ SCENIC module builder can include negative regulons", {
  expr <- matrix(
    c(
      1, 2, 3, 4, 5, 6,
      1, 2, 3, 4, 5, 6,
      6, 5, 4, 3, 2, 1
    ),
    nrow = 6,
    ncol = 3,
    dimnames = list(paste0("Cell", 1:6), c("TF1", "GenePos", "GeneNeg"))
  )
  adjacency <- data.frame(
    TF = "TF1",
    target = c("GenePos", "GeneNeg"),
    importance = c(0.9, 0.8),
    stringsAsFactors = FALSE
  )

  modules <- getFromNamespace("scenic_modules_from_adjacencies", "scop")(
    adjacency = adjacency,
    expr_mtx = expr,
    thresholds = 0.5,
    top_n_targets = 2,
    top_n_regulators = integer(0),
    min_genes = 2,
    keep_only_activating = FALSE,
    verbose = FALSE
  )

  suffixes <- unique(vapply(modules, `[[`, character(1), "suffix"))
  expect_setequal(suffixes, c("(-)", "(+)"))
  expect_true(any(vapply(modules, function(module) {
    identical(module$suffix, "(-)") && "GeneNeg" %in% module$genes
  }, logical(1))))
  expect_true(any(vapply(modules, function(module) {
    identical(module$suffix, "(+)") && "GenePos" %in% module$genes
  }, logical(1))))
})

test_that("C++ SCENIC regulon fallback uses pySCENIC-style positive names", {
  adjacency <- data.frame(
    TF = "TF1",
    target = c("Gene1", "Gene2"),
    importance = c(0.9, 0.8),
    stringsAsFactors = FALSE
  )

  regulons <- getFromNamespace("build_regulons", "scop")(
    adjacency = adjacency,
    max_targets = 2,
    min_regulon_size = 2
  )

  expect_equal(names(regulons), "TF1(+)")
})

test_that("C++ SCENIC AUCell scores are invariant to positive regulon naming", {
  expr <- matrix(
    rpois(30 * 24, lambda = 1),
    nrow = 30,
    ncol = 24,
    dimnames = list(paste0("Gene", seq_len(30)), paste0("Cell", seq_len(24)))
  )
  expr[1:4, 1:6] <- expr[1:4, 1:6] + 5
  genes <- c("Gene1", "Gene2", "Gene3", "Gene4")

  old_scores <- getFromNamespace("scenic_compute_aucell_score", "scop")(
    counts = expr,
    regulon_list = list(TF1 = genes),
    min_regulon_size = 1,
    backend = "cpp",
    cpp_strategy = "sparse",
    seed = 42,
    verbose = FALSE
  )
  new_scores <- getFromNamespace("scenic_compute_aucell_score", "scop")(
    counts = expr,
    regulon_list = list("TF1(+)" = genes),
    min_regulon_size = 1,
    backend = "cpp",
    cpp_strategy = "sparse",
    seed = 42,
    verbose = FALSE
  )

  expect_equal(unname(old_scores[[1]]), unname(new_scores[[1]]))
})

test_that("C++ SCENIC AUCell scores are reproducible with seed", {
  expr <- matrix(
    rpois(30 * 24, lambda = 1),
    nrow = 30,
    ncol = 24,
    dimnames = list(paste0("Gene", seq_len(30)), paste0("Cell", seq_len(24)))
  )
  genes <- c("Gene1", "Gene2", "Gene3", "Gene4")

  first <- getFromNamespace("scenic_compute_aucell_score", "scop")(
    counts = expr,
    regulon_list = list("TF1(+)" = genes),
    min_regulon_size = 1,
    backend = "cpp",
    cpp_strategy = "sparse",
    seed = 42,
    verbose = FALSE
  )
  second <- getFromNamespace("scenic_compute_aucell_score", "scop")(
    counts = expr,
    regulon_list = list("TF1(+)" = genes),
    min_regulon_size = 1,
    backend = "cpp",
    cpp_strategy = "sparse",
    seed = 42,
    verbose = FALSE
  )

  expect_equal(first, second)
})
