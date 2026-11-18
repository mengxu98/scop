test_that("Pando views retain interaction coefficients and unknown membership signs", {
  native <- list(
    coefficients = data.frame(tf = c("01", "01"), region = c("r1", "r2"), target = c("g1", "g1"), estimate = c(-2, 0.5), pval = c(0.01, 0.02), padj = c(0.02, 0.03)),
    goodness_of_fit = data.frame(target = "g1", rsq = 0.7, nvariables = 2),
    modules = data.frame(tf = "01", target = "g1", estimate = -1.5, pval = 0.01, padj = 0.02),
    tf_region = data.frame(tf = c("01", "01"), region = c("r1", "r2"), membership = c(1, 1)),
    region_gene = data.frame(region = c("r1", "r2"), target = c("g1", "g1"), membership = c(1, 1)),
    gene_sets = list("01(-)" = "g1"),
    parameters = list(p_thresh = 0.05),
    package_versions = c(Pando = "1.1.1")
  )
  native$triplets <- native$coefficients
  result <- grn_pando_group(native, list(cells = c("cell2", "cell1")), list())
  expect_identical(result$native, native)
  expect_identical(result$cells, c("cell2", "cell1"))
  expect_identical(result$native$triplets$estimate, c(-2, 0.5))
  expect_identical(result$layer_specs$tf_gene$score, "estimate")
  expect_identical(result$layer_specs$tf_gene$sign_column, "estimate")
  for (layer in c("tf_region", "region_gene")) {
    expect_identical(result$layer_specs[[layer]]$score_type, "membership")
    expect_null(result$layer_specs[[layer]]$sign_column)
    expect_identical(result$layer_specs[[layer]]$sign_basis, "unknown")
  }
  expect_identical(result$chains$table_path, "triplets")
  expect_identical(result$gene_sets$table_path, "gene_sets")
  expect_identical(result$execution$status, "completed")
})

test_that("Pando retains successful empty native layers", {
  native <- list(
    coefficients = data.frame(tf = character(), region = character(), target = character(), estimate = numeric(), pval = numeric(), padj = numeric()),
    goodness_of_fit = data.frame(target = character(), rsq = numeric(), nvariables = integer()),
    modules = data.frame(tf = character(), target = character(), estimate = numeric()),
    tf_region = data.frame(tf = character(), region = character(), membership = numeric()),
    region_gene = data.frame(region = character(), target = character(), membership = numeric()),
    gene_sets = list(), parameters = list(), package_versions = c(Pando = "1.1.1")
  )
  native$triplets <- native$coefficients
  result <- grn_pando_group(native, list(cells = "cell1"), list())
  expect_identical(result$execution$status, "completed")
  expect_true(all(vapply(result$layer_specs, function(spec) identical(spec$status, "computed"), logical(1))))
  expect_equal(nrow(result$native$modules), 0L)
  expect_equal(nrow(result$native$triplets), 0L)
})

test_that("Pando rejects incomplete native outputs and missing multiome input", {
  expect_error(grn_pando_group(list(), list(cells = "cell1"), list()), "incomplete native result")
  expect_error(RunGRN(matrix(1, 2, 2), methods = "pando"), "requires a Seurat object")
  counts <- Matrix::sparseMatrix(i = c(1L, 1L, 2L), j = c(1L, 2L, 2L), x = c(1, 2, 3), dims = c(2L, 2L), dimnames = list(c("a", "b"), c("c1", "c2")))
  object <- SeuratObject::CreateSeuratObject(counts)
  expect_error(infer_pando(list(object = object), list(resources = list())), "resources\\$genome")
  expect_error(RunGRN(object, methods = "pando", method_args = list(typo = 1)), "Unknown argument")
  expect_error(infer_pando(list(object = object), list(resources = list(genome = TRUE), method_args = list(interaction_term = "*"))), "interaction_term")
  expect_error(infer_pando(list(object = object), list(resources = list(genome = TRUE), method_args = list(extend = 5000))), "custom extend")
  expect_identical(grn_pando_defaults()$p_thresh, 0.05)
})

test_that("Pando module inference preserves multiCSN dispatch in both load orders", {
  skip_if_not_installed("callr")
  skip_if(!length(find.package("Pando", quiet = TRUE)), "Pando is not installed")
  skip_if_not_installed("multiCSN", minimum_version = "0.1.0")
  for (pando_first in c(FALSE, TRUE)) {
    result <- callr::r(function(pando_first) {
      if (pando_first) loadNamespace("Pando")
      loadNamespace("multiCSN")
      before <- search()
      network <- methods::new(methods::getClassDef("Network", where = asNamespace("multiCSN")))
      coefficients <- stats::coef(network)
      method <- getS3method("coef", "Network")
      functions <- getFromNamespace("grn_pando_functions", "scop")(1L)
      native <- methods::new(methods::getClassDef("Network", where = asNamespace("Pando")),
        fit = data.frame(target = "gene", rsq = 1, nvariables = 20L),
        coefs = data.frame(target = "gene", tf = "tf", region = "peak", estimate = 2, pval = 0, padj = 0),
        params = list(method = "glm")
      )
      fitted <- environment(functions$find_modules)$find_modules(native, verbose = FALSE)
      list(
        estimate = getExportedValue("Pando", "NetworkModules")(fitted)@meta$estimate,
        multiCSN = identical(stats::coef(network), coefficients) &&
          identical(getS3method("coef", "Network"), method),
        attached = !identical(search(), before)
      )
    }, args = list(pando_first), libpath = .libPaths())
    expect_identical(result, list(estimate = 2, multiCSN = TRUE, attached = FALSE))
  }
})
