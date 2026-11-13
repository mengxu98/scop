test_that("RNA adapters reject invalid identifiers and method arguments before fitting", {
  expression <- matrix(seq_len(24), nrow = 3L, dimnames = list(c("01", "1", "target"), paste0("cell", seq_len(8L))))
  input <- list(expression = expression, regulators = c("01", "1"), targets = "target", cells = colnames(expression))
  config <- list(method = "ppcor", method_args = list(typo = TRUE), work_dir = tempdir(), cores = 1L, seed = 1L, verbose = FALSE)
  expect_error(grn_run_rna(input, config), "Unknown ppcor argument")
  config$method_args <- list()
  rownames(input$expression)[[2L]] <- "01"
  expect_error(grn_run_rna(input, config), "unique gene and cell")
  input$expression <- expression
  input$expression[1L, 1L] <- Inf
  expect_error(grn_run_rna(input, config), "finite numeric")
})

test_that("PPCOR retains negative associations and numeric-like gene identifiers", {
  skip_if_not_installed("ppcor")
  set.seed(401)
  expression <- matrix(stats::rnorm(200), nrow = 4L,
    dimnames = list(c("01", "1", "TF-A", "target"), paste0("cell", seq_len(50L))))
  expression[4L, ] <- -expression[1L, ] + stats::rnorm(50L, sd = 0.1)
  input <- list(expression = expression, regulators = c("01", "1"), targets = "target", cells = colnames(expression))
  config <- list(method = "ppcor", method_args = list(), work_dir = tempdir(), cores = 1L, seed = 1L, verbose = FALSE)
  result <- grn_run_rna(input, config)
  expect_setequal(result$native$edges$regulator, c("01", "1"))
  expect_identical(result$native$edges$target, rep("target", 2L))
  expect_true(any(result$native$edges$weight < 0))
  expect_false(result$layer_specs$tf_gene$directed)
  expect_identical(result$layer_specs$tf_gene$sign_basis, "association")
  expected <- ppcor::pcor(t(expression[c("01", "1", "target"), ]), method = "spearman")
  indices <- cbind(match(result$native$edges$regulator, rownames(expected$estimate)), match(result$native$edges$target, colnames(expected$estimate)))
  expect_identical(result$native$edges$weight, as.numeric(expected$estimate[indices]))
  expect_identical(result$native$edges$p_value, as.numeric(expected$p.value[indices]))
})

test_that("PIDC maps rectangular native regulator rows to target columns", {
  skip_if_not_installed("pbapply")
  skip_if_not_installed("purrr")
  source_dir <- Sys.getenv("SCOP_GRN_SOURCE_DIR")
  skip_if(!nzchar(source_dir), "Set SCOP_GRN_SOURCE_DIR to a multiCSN manuscript checkout")
  set.seed(402)
  expression <- matrix(stats::runif(300), nrow = 5L,
    dimnames = list(c("01", "1", "TF-A", "target-a", "target-b"), paste0("cell", seq_len(60L))))
  input <- list(expression = expression, regulators = c("01", "1", "TF-A"), targets = c("target-a", "target-b"), cells = colnames(expression))
  config <- list(method = "pidc", method_args = list(), resources = list(source_dir = source_dir),
    work_dir = tempdir(), cores = 1L, seed = 1L, verbose = FALSE)
  result <- grn_run_rna(input, config)
  expect_setequal(result$native$edges$regulator, input$regulators)
  expect_setequal(result$native$edges$target, input$targets)
  expect_equal(nrow(result$native$edges), 6L)
  expect_false(result$layer_specs$tf_gene$directed)
  expect_null(result$layer_specs$tf_gene$sign_column)
  upstream <- new.env(parent = baseenv())
  upstream$is <- methods::is
  sys.source(file.path(source_dir, "code", "vendor", "PIDC", "PIDC.R"), envir = upstream)
  expected <- upstream$PIDC(expression, regulators = input$regulators, targets = input$targets, diag = "zero", verbose = FALSE)
  indices <- cbind(match(result$native$edges$regulator, rownames(expected)), match(result$native$edges$target, colnames(expected)))
  expect_identical(result$native$edges$weight, as.numeric(expected[indices]))
})
