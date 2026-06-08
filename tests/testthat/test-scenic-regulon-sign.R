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

test_that("C++ SCENIC negative modules do not change positive module targets", {
  expr <- matrix(
    c(
      1, 2, 3, 4, 5, 6,
      1, 2, 3, 4, 5, 6,
      2, 3, 4, 5, 6, 7,
      3, 4, 5, 6, 7, 8,
      6, 5, 4, 3, 2, 1,
      7, 6, 5, 4, 3, 2
    ),
    nrow = 6,
    ncol = 6,
    dimnames = list(
      paste0("Cell", 1:6),
      c("TF1", "GenePos1", "GenePos2", "GenePos3", "GeneNeg1", "GeneNeg2")
    )
  )
  adjacency <- data.frame(
    TF = "TF1",
    target = c("GeneNeg1", "GeneNeg2", "GenePos1", "GenePos2", "GenePos3"),
    importance = c(0.99, 0.98, 0.50, 0.49, 0.48),
    stringsAsFactors = FALSE
  )

  positive_only <- getFromNamespace("scenic_modules_from_adjacencies", "scop")(
    adjacency = adjacency,
    expr_mtx = expr,
    thresholds = 0.5,
    top_n_targets = 3,
    top_n_regulators = integer(0),
    min_genes = 2,
    keep_only_activating = TRUE,
    verbose = FALSE
  )
  with_negative <- getFromNamespace("scenic_modules_from_adjacencies", "scop")(
    adjacency = adjacency,
    expr_mtx = expr,
    thresholds = 0.5,
    top_n_targets = 3,
    top_n_regulators = integer(0),
    min_genes = 2,
    keep_only_activating = FALSE,
    verbose = FALSE
  )
  with_negative_positive <- with_negative[
    vapply(with_negative, `[[`, character(1), "suffix") == "(+)"
  ]

  expect_equal(
    lapply(positive_only, `[[`, "genes"),
    lapply(with_negative_positive, `[[`, "genes")
  )
})

test_that("C++ SCENIC threshold modules use global pySCENIC weight quantiles", {
  expr <- matrix(
    c(
      1, 2, 3, 4, 5, 6,
      1, 2, 3, 4, 5, 6,
      2, 3, 4, 5, 6, 7,
      3, 4, 5, 6, 7, 8,
      6, 5, 4, 3, 2, 1,
      7, 6, 5, 4, 3, 2,
      8, 7, 6, 5, 4, 3
    ),
    nrow = 6,
    ncol = 7,
    dimnames = list(
      paste0("Cell", 1:6),
      c("TF1", "GenePos1", "GenePos2", "GenePos3", "GeneNeg1", "GeneNeg2", "GeneNeg3")
    )
  )
  adjacency <- data.frame(
    TF = "TF1",
    target = c("GenePos1", "GenePos2", "GenePos3", "GeneNeg1", "GeneNeg2", "GeneNeg3"),
    importance = c(0.99, 0.98, 0.97, 0.30, 0.29, 0.28),
    stringsAsFactors = FALSE
  )

  modules <- getFromNamespace("scenic_modules_from_adjacencies", "scop")(
    adjacency = adjacency,
    expr_mtx = expr,
    thresholds = 0.5,
    top_n_targets = 3,
    top_n_regulators = integer(0),
    min_genes = 2,
    keep_only_activating = FALSE,
    verbose = FALSE
  )

  negative_threshold_modules <- modules[
    vapply(modules, function(module) {
      identical(module$suffix, "(-)") && identical(module$context, "weight>50%")
    }, logical(1))
  ]
  expect_length(negative_threshold_modules, 0)
})

test_that("C++ SCENIC module builder matches pySCENIC module construction", {
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not(
    reticulate::py_module_available("pyscenic") &&
      reticulate::py_module_available("pandas"),
    "pySCENIC and pandas are required for SCENIC module-builder parity"
  )

  expr <- matrix(
    c(
      1, 2, 3, 4, 5, 6,
      1, 2, 3, 4, 5, 6,
      2, 3, 4, 5, 6, 7,
      3, 4, 5, 6, 7, 8,
      6, 5, 4, 3, 2, 1,
      7, 6, 5, 4, 3, 2,
      8, 7, 6, 5, 4, 3,
      6, 5, 4, 3, 2, 1,
      7, 6, 5, 4, 3, 2,
      8, 7, 6, 5, 4, 3
    ),
    nrow = 6,
    ncol = 10,
    dimnames = list(
      paste0("Cell", 1:6),
      c(
        "TF1", "GenePos1", "GenePos2", "GenePos3",
        "GeneNeg1", "GeneNeg2", "GeneNeg3",
        "GeneNeg4", "GeneNeg5", "GeneNeg6"
      )
    )
  )
  adjacency <- data.frame(
    TF = "TF1",
    target = c(
      "GenePos1", "GenePos2", "GenePos3",
      "GeneNeg1", "GeneNeg2", "GeneNeg3",
      "GeneNeg4", "GeneNeg5", "GeneNeg6"
    ),
    importance = c(0.99, 0.98, 0.97, 0.40, 0.39, 0.38, 0.37, 0.36, 0.35),
    stringsAsFactors = FALSE
  )

  cpp_modules <- getFromNamespace("scenic_modules_from_adjacencies", "scop")(
    adjacency = adjacency,
    expr_mtx = expr,
    thresholds = 0.5,
    top_n_targets = 3,
    top_n_regulators = integer(0),
    min_genes = 2,
    keep_only_activating = FALSE,
    verbose = FALSE
  )
  cpp_modules <- lapply(cpp_modules, function(module) {
    list(
      tf = module$tf,
      sign = if (identical(module$suffix, "(-)")) "repressing" else "activating",
      context = sub("\\.0%", "%", module$context),
      genes = sort(setdiff(module$genes, module$tf))
    )
  })

  pd <- reticulate::import("pandas", convert = FALSE)
  py_utils <- reticulate::import("pyscenic.utils", convert = FALSE)
  py_builtins <- reticulate::import_builtins(convert = FALSE)
  py_adj <- pd$DataFrame(reticulate::r_to_py(adjacency))
  expr_df <- as.data.frame(expr, check.names = FALSE)
  py_expr <- pd$DataFrame(reticulate::r_to_py(expr_df))
  py_expr$index <- reticulate::r_to_py(rownames(expr))
  py_modules <- py_utils$modules_from_adjacencies(
    py_adj,
    py_expr,
    thresholds = reticulate::tuple(0.5),
    top_n_targets = reticulate::tuple(3L),
    top_n_regulators = reticulate::tuple(),
    min_genes = 2L,
    rho_dichotomize = TRUE,
    keep_only_activating = FALSE,
    rho_threshold = 0.03
  )
  py_modules <- reticulate::py_to_r(py_modules)
  py_modules <- lapply(py_modules, function(module) {
    context <- as.character(reticulate::py_to_r(py_builtins$list(module$context)))
    sign <- if ("repressing" %in% context) "repressing" else "activating"
    context <- sub("\\.0%", "%", setdiff(context, c("activating", "repressing")))
    genes <- setdiff(reticulate::py_to_r(py_builtins$list(module$genes)), module$transcription_factor)
    list(
      tf = module$transcription_factor,
      sign = sign,
      context = context,
      genes = sort(as.character(genes))
    )
  })

  normalize_modules <- function(modules) {
    modules <- lapply(modules, function(module) {
      module$context <- sort(as.character(module$context))
      module$genes <- sort(as.character(module$genes))
      module
    })
    modules[order(
      vapply(modules, `[[`, character(1), "tf"),
      vapply(modules, `[[`, character(1), "sign"),
      vapply(modules, function(module) paste(module$context, collapse = "|"), character(1)),
      vapply(modules, function(module) paste(module$genes, collapse = "|"), character(1))
    )]
  }

  expect_equal(normalize_modules(cpp_modules), normalize_modules(py_modules))
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
