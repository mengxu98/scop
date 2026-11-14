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
  expect_identical(RunGRN(expression, regulators = regulators, targets = targets,
    methods = "GRNBoost2", n_rounds = 20, verbose = FALSE), direct)
})

test_that("result grouping and object storage preserve native scores", {
  expression <- grn_contract_expression()
  grouping <- stats::setNames(rep(c("A", "B"), each = 24L), colnames(expression))
  work <- tempfile("grn-contract-")
  result <- RunGRN(expression, regulators = rownames(expression)[1:4],
    targets = rownames(expression)[5:10], methods = "GRNBoost2", group.by = grouping,
    method_args = list(n_rounds = 20), return_type = "result", work_dir = work, verbose = FALSE)
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
    tool_name = "net", work_dir = tempfile(), verbose = FALSE)
  expect_s4_class(stored, "Seurat")
  stored@tools$unrelated <- list(network = FALSE)
  expect_identical(GetGRN(stored, what = "result"), GetGRN(stored, tool_name = "net", what = "result"))
  expect_error(GetGRN(stored), "Select a group")
  expect_s3_class(GRNPlot(stored, group = "A"), "ggplot")
  expect_error(GetGRN(object), "No GRNResult")
  expect_error(GetGRN(stored, tool_name = "unrelated"), "identify a stored GRNResult")
  expect_error(RunGRN(stored, methods = "grnboost2", tool_name = "net", verbose = FALSE), "already exists")
  edge_file <- tempfile()
  edge_table <- RunGRN(object, return_type = "edges", output_file = edge_file, n_rounds = 20, verbose = FALSE)
  expect_named(edge_table, c("TF", "target", "importance"))
  rownames(edge_table) <- NULL
  expect_equal(utils::read.delim(edge_file), edge_table)
  second <- RunGRN(stored, methods = "grnboost2", tool_name = "second", n_rounds = 20, verbose = FALSE)
  expect_error(GetGRN(second), "Specify tool_name: net, second")
  expect_error(GRNPlot(second), "Specify tool_name")
  expect_s3_class(GRNPlot(second, tool_name = "second"), "ggplot")
  expect_s3_class(GRNPlot(second, tool_name = "second", compare = stored,
    compare_group = "A", type = "overlap"), "ggplot")
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
  expect_error(RunGRN(expression, grn_meth = "gniplr"), "Use 'methods'")
  expect_error(RunGRN(expression, method = "grnboost2"), "Use 'methods'")
  expect_error(RunGRN(expression, methods = "grnboost2", meth = "gniplr"), "Use 'methods'")
  expect_error(RunGRN(expression, return_dat = TRUE), "Use 'return_type'")
  expect_error(RunGRN(object, return_seurat = TRUE), "Use 'return_type'")
  expect_error(RunGRN(expression, return_t = "result"), "Use 'return_type'")
  expect_error(RunGRN(object, return_type = NA), "must be one")
  expect_error(RunGRN(object, return_type = "r"), "must be one")
  expect_error(RunGRN(object, output_file = tempfile()), "output_file requires")
  expect_error(RunGRN(object, tool_name = ""), "one non-empty name")
  expect_error(RunGRN(object, tool_name = NA_character_), "one non-empty name")
  expect_error(RunGRN(expression, methods = "ppcor", method_args = list(typo = 1)), "Unknown argument")
  expect_error(RunGRN(expression, methods = "grnboost2", n_rounds = 2,
    method_args = list(n_rounds = 3)), "repeat")
  expect_error(RunGRN(expression, methods = "leap", return_type = "result"), "requires named pseudotime")
  expect_error(RunGRN(expression, methods = "multicsn", input_mode = "multiome"), "Seurat")
  expect_error(RunGRN(expression, methods = "ppcor", return_type = "seurat"), "Seurat")
})

test_that("multiome resume rejects ambiguous groups and incomplete LINGER checkpoints", {
  expression <- grn_contract_expression()
  object <- SeuratObject::CreateSeuratObject(methods::as(expression, "dgCMatrix"))
  atac <- expression[seq_len(2L), , drop = FALSE]
  rownames(atac) <- c("chr1-1-50", "chr1-60-100")
  object[["peaks"]] <- SeuratObject::CreateAssayObject(methods::as(atac, "dgCMatrix"))
  object$sample <- rep(c("A", "B"), each = 24L)
  expect_error(RunGRN(object, methods = "scenicplus", return_type = "result",
    atac_assay = "peaks", group.by = "sample", method_args = list(resume_dir = tempfile())),
    "resume_dir requires exactly one selected group")
  expect_error(RunGRN(object, methods = "linger", return_type = "result",
    atac_assay = "peaks", group.by = "sample", method_args = list(resume_dir = tempfile())),
    "resume_dir requires exactly one selected group")
  directory <- withr::local_tempdir()
  expect_error(RunGRN(object, methods = "linger", return_type = "result",
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
    return_type = "result", work_dir = tempfile(), verbose = FALSE)
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
    methods = "genie3", nTrees = 20, seed = 7, verbose = FALSE)
  set.seed(7)
  weights <- GENIE3::GENIE3(exprMatrix = expression, regulators = regulators,
    targets = targets, nTrees = 20, nCores = 1L)
  direct <- GENIE3::getLinkList(weights)
  names(direct)[1:3] <- c("TF", "target", "importance")
  direct <- cap_grn_edges_per_target(direct)
  expect_identical(actual, direct)
  set.seed(7)
  columns <- RunGRN(t(expression), genes_in = "columns", regulators = regulators,
    targets = targets, methods = "genie3", nTrees = 20, seed = 7, verbose = FALSE)
  expect_identical(actual, columns)
  result <- RunGRN(expression, regulators = regulators, targets = targets,
    methods = "genie3", nTrees = 20, seed = 7, return_type = "result", verbose = FALSE)
  expect_identical(result$implementation, "GENIE3::GENIE3")
  expect_identical(result$package_versions$all, c(GENIE3 = as.character(utils::packageVersion("GENIE3"))))
})

test_that("return format does not change inference or the GNIPLR candidate limit", {
  set.seed(1)
  expression <- matrix(rnorm(120), 4L, dimnames = list(LETTERS[1:4], paste0("c", 1:30)))
  args <- list(object = expression, methods = "gniplr", regulators = "D", targets = "B",
    correlation_threshold = 0, lasso_degree = 3L, max_lag = 1L,
    max_edges_per_target = 1L, verbose = FALSE)
  edges <- do.call(RunGRN, args)
  result <- GetGRN(do.call(RunGRN, c(args, list(return_type = "result"))))
  expect_identical(edges$TF, "D")
  expect_identical(edges$target, "B")
  expect_identical(result$source, edges$TF)
  expect_identical(result$raw_score, edges$importance)
})
