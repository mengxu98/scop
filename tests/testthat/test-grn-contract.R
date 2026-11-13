grn_contract_expression <- function() {
  set.seed(1201)
  expression <- matrix(stats::rpois(12L * 48L, 3), 12L,
    dimnames = list(paste0("TF", seq_len(12L)), paste0("C", seq_len(48L))))
  expression[5L, ] <- 2 * expression[1L, ] + expression[2L, ]
  expression
}

test_that("default GRN calls preserve legacy numerical results and position", {
  expression <- grn_contract_expression()
  regulators <- rownames(expression)[1:4]
  targets <- rownames(expression)[5:10]
  actual <- RunGRN(expression, regulators, targets, "rows", "grnboost2", "cpp",
    n_rounds = 20, verbose = FALSE)
  direct <- RunGRNBoost2(expression, regulators = regulators, targets = targets,
    genes_in = "rows", n_rounds = 20, verbose = FALSE)
  expect_identical(actual, direct)
  expect_identical(names(formals(RunGRN.default))[1:5], c("object", "regulators", "targets", "genes_in", "methods"))
  old <- names(formals(grn_legacy_inference))
  old[old == "grn_method"] <- "methods"
  expect_identical(names(formals(RunGRN.default))[seq_along(old)], old)
  expect_identical(RunGRN(expression, regulators = regulators, targets = targets,
    methods = "GRNBoost2", n_rounds = 20, verbose = FALSE), direct)
})

test_that("result grouping and object storage preserve native scores", {
  expression <- grn_contract_expression()
  grouping <- stats::setNames(rep(c("A", "B"), each = 24L), colnames(expression))
  work <- tempfile("grn-contract-")
  result <- RunGRN(expression, regulators = rownames(expression)[1:4],
    targets = rownames(expression)[5:10], methods = "GRNBoost2", group.by = grouping,
    method_args = list(n_rounds = 20), return_data = TRUE, work_dir = work, verbose = FALSE)
  expect_s3_class(result, "GRNResult")
  expect_identical(names(result$groups), c("A", "B"))
  expect_identical(result$groups$A$cells, colnames(expression)[1:24])
  expect_error(GetGRN(result), "Select a group")
  expect_named(GetGRN(result, all_groups = TRUE), c("A", "B"))
  expect_identical(GetGRN(result, group = "A")$raw_score, result$groups$A$native$importance)
  expect_true(all(is.na(GetGRN(result, group = "A")$sign)))
  expect_error(RunGRN(expression, group.by = grouping, n_rounds = 20, verbose = FALSE), "multiple groups")
  object <- SeuratObject::CreateSeuratObject(methods::as(expression, "dgCMatrix"))
  object$group <- unname(grouping)
  stored <- RunGRN(object, methods = "grnboost2", group.by = "group", n_rounds = 20,
    return_seurat = TRUE, tool_name = "net", work_dir = tempfile(), verbose = FALSE)
  expect_s3_class(GetGRN(stored, tool_name = "net", what = "result"), "GRNResult")
  expect_error(RunGRN(stored, methods = "grnboost2", return_seurat = TRUE, tool_name = "net", verbose = FALSE), "already exists")
  path <- tempfile(fileext = ".rds")
  saveRDS(result, path)
  expect_identical(GetGRN(readRDS(path), group = "A"), GetGRN(result, group = "A"))
})

test_that("method boundaries reject unknown arguments and missing time", {
  expression <- grn_contract_expression()
  expect_error(RunGRN(expression, methods = c("pidc", "ppcor")), "exactly one")
  expect_error(RunGRN(expression, methods = "missing"), "Unknown GRN method")
  object <- SeuratObject::CreateSeuratObject(methods::as(expression, "dgCMatrix"))
  expect_error(RunGRN(expression, grn_method = "gniplr"), "Use 'methods'")
  expect_error(RunGRN(object, grn_method = "gniplr"), "Use 'methods'")
  expect_error(RunGRN(expression, return_type = "result"), "Use 'return_data'")
  expect_error(RunGRN(object, return_data = TRUE, return_seurat = TRUE), "Choose only one")
  expect_error(RunGRN(expression, return_data = NA), "TRUE or FALSE")
  expect_error(RunGRN(expression, methods = "ppcor", method_args = list(typo = 1)), "Unknown argument")
  expect_error(RunGRN(expression, methods = "grnboost2", n_rounds = 2,
    method_args = list(n_rounds = 3)), "repeat")
  expect_error(RunGRN(expression, methods = "leap", return_data = TRUE), "requires named pseudotime")
  expect_error(RunGRN(expression, methods = "multicsn", input_mode = "multiome"), "Seurat")
  expect_error(RunGRN(expression, methods = "ppcor", return_seurat = TRUE), "Seurat")
})

test_that("multiome resume rejects ambiguous groups and incomplete LINGER checkpoints", {
  expression <- grn_contract_expression()
  object <- SeuratObject::CreateSeuratObject(methods::as(expression, "dgCMatrix"))
  atac <- expression[seq_len(2L), , drop = FALSE]
  rownames(atac) <- c("chr1-1-50", "chr1-60-100")
  object[["peaks"]] <- SeuratObject::CreateAssayObject(methods::as(atac, "dgCMatrix"))
  object$sample <- rep(c("A", "B"), each = 24L)
  expect_error(RunGRN(object, methods = "scenicplus", return_data = TRUE,
    atac_assay = "peaks", group.by = "sample", method_args = list(resume_dir = tempfile())),
    "resume_dir requires exactly one selected group")
  expect_error(RunGRN(object, methods = "linger", return_data = TRUE,
    atac_assay = "peaks", group.by = "sample", method_args = list(resume_dir = tempfile())),
    "resume_dir requires exactly one selected group")
  directory <- withr::local_tempdir()
  expect_error(RunGRN(object, methods = "linger", return_data = TRUE,
    atac_assay = "peaks", group.by = "sample", groups = "A",
    method_args = list(resume_dir = directory, runtime = list(python = tempfile()))),
    "resume_dir must contain a completed prepare checkpoint")
})

test_that("new GNIPLR views restrict candidates and record effective defaults", {
  expression <- grn_contract_expression()
  regulators <- rownames(expression)[1:4]
  targets <- rownames(expression)[5:10]
  result <- RunGRN(expression, regulators = regulators, targets = targets,
    methods = "gniplr", method_args = list(correlation_threshold = 0, lasso_degree = 1, max_lag = 1),
    return_data = TRUE, work_dir = tempfile(), verbose = FALSE)
  expect_true(all(GetGRN(result)$source %in% regulators))
  expect_true(all(GetGRN(result)$target %in% targets))
  expect_identical(GetGRN(result, what = "parameters")$lasso_alpha, 0.1)
  expect_identical(GetGRN(result, what = "parameters")$max_edges_per_target, Inf)
  expect_identical(GetGRN(result, what = "parameters")$candidate_filter, "postfit_view")
})

test_that("GENIE3 accepts gene-by-cell input and preserves orientation", {
  skip_if_not_installed("GENIE3")
  expression <- grn_contract_expression()
  regulators <- rownames(expression)[1:4]
  targets <- rownames(expression)[5:10]
  set.seed(7)
  actual <- RunGRN(expression, regulators = regulators, targets = targets,
    methods = "genie3", nTrees = 20, verbose = FALSE)
  set.seed(7)
  weights <- GENIE3::GENIE3(exprMatrix = expression, regulators = regulators,
    targets = targets, nTrees = 20, nCores = 1L)
  direct <- GENIE3::getLinkList(weights)
  names(direct)[1:3] <- c("TF", "target", "importance")
  direct <- cap_grn_edges_per_target(direct)
  expect_identical(actual, direct)
  set.seed(7)
  columns <- RunGRN(t(expression), genes_in = "columns", regulators = regulators,
    targets = targets, methods = "genie3", nTrees = 20, verbose = FALSE)
  expect_identical(actual, columns)
})
