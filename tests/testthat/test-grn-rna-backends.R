test_that("RNA adapters reject invalid identifiers and method arguments before fitting", {
  expression <- matrix(seq_len(24), nrow = 3L, dimnames = list(c("01", "1", "target"), paste0("cell", seq_len(8L))))
  expect_error(
    RunGRN(expression, methods = "ppcor", method_args = list(typo = TRUE), verbose = FALSE),
    "Unknown argument.*ppcor.*typo"
  )
  rownames(expression)[[2L]] <- "01"
  expect_error(RunGRN(expression, methods = "ppcor", verbose = FALSE), "unique.*gene and cell")
  rownames(expression)[[2L]] <- "1"
  expression[1L, 1L] <- Inf
  expect_error(RunGRN(expression, methods = "ppcor", verbose = FALSE), "finite numbers")
})

test_that("PPCOR retains negative associations and numeric-like gene identifiers", {
  skip_if_not_installed("ppcor")
  set.seed(401)
  expression <- matrix(stats::rnorm(200),
    nrow = 4L,
    dimnames = list(c("01", "1", "TF-A", "target"), paste0("cell", seq_len(50L)))
  )
  expression[4L, ] <- -expression[1L, ] + stats::rnorm(50L, sd = 0.1)
  input <- list(expression = expression, regulators = c("01", "1"), targets = "target", cells = colnames(expression))
  config <- list(method = "ppcor", method_args = list(), work_dir = tempdir(), cores = 1L, seed = 1L, verbose = FALSE)
  result <- grn_run_rna(input, config)
  pearson <- RunGRN(expression,
    methods = "ppcor", method = "pearson", regulators = input$regulators,
    targets = input$targets, return_type = "result", verbose = FALSE
  )
  expect_identical(GetGRN(pearson, what = "parameters")$method, "pearson")
  object <- SeuratObject::CreateSeuratObject(methods::as(abs(expression), "dgCMatrix"))
  stored <- RunGRN(object, methods = "ppcor", method = "pearson", verbose = FALSE)
  expect_identical(GetGRN(stored, what = "parameters")$method, "pearson")
  expect_error(RunGRN(expression,
    methods = "ppcor", method = "pearson",
    method_args = list(method = "kendall")
  ), "named lists")
  expect_setequal(result$native$edges$regulator, c("01", "1"))
  expect_identical(result$native$edges$target, rep("target", 2L))
  expect_true(any(result$native$edges$weight < 0))
  expect_false(result$layer_specs$tf_gene$directed)
  expect_identical(result$layer_specs$tf_gene$sign_basis, "association")
  constant <- RunPPCOR(rbind(expression, constant = 0),
    method = "pearson",
    return_type = "result", verbose = FALSE
  )
  regular <- RunPPCOR(expression, method = "pearson", return_type = "result", verbose = FALSE)
  expect_identical(GetGRN(constant), GetGRN(regular))
  expected <- ppcor::pcor(t(expression[c("01", "1", "target"), ]), method = "spearman")
  indices <- cbind(match(result$native$edges$regulator, rownames(expected$estimate)), match(result$native$edges$target, colnames(expected$estimate)))
  expect_identical(result$native$edges$weight, as.numeric(expected$estimate[indices]))
  expect_identical(result$native$edges$p_value, as.numeric(expected$p.value[indices]))
})

test_that("PIDC maps rectangular native regulator rows to target columns", {
  skip_if_not_installed("PIDC")
  set.seed(402)
  expression <- matrix(stats::runif(300),
    nrow = 5L,
    dimnames = list(c("01", "1", "TF-A", "target-a", "target-b"), paste0("cell", seq_len(60L)))
  )
  input <- list(expression = expression, regulators = c("01", "1", "TF-A"), targets = c("target-a", "target-b"), cells = colnames(expression))
  result <- RunGRN(expression,
    methods = "pidc", regulators = input$regulators, targets = input$targets,
    return_type = "result", verbose = FALSE
  )$groups$all
  expect_setequal(result$native$edges$regulator, input$regulators)
  expect_setequal(result$native$edges$target, input$targets)
  expect_equal(nrow(result$native$edges), 6L)
  expect_false(result$layer_specs$tf_gene$directed)
  expect_null(result$layer_specs$tf_gene$sign_column)
  expected <- getExportedValue("PIDC", "PIDC")(
    expression, regulators = input$regulators, targets = input$targets, diag = "zero", verbose = FALSE
  )
  indices <- cbind(match(result$native$edges$regulator, rownames(expected)), match(result$native$edges$target, colnames(expected)))
  expect_identical(result$native$edges$weight, as.numeric(expected[indices]))
})
