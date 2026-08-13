test_that("scCODA preserves sample-condition pairs", {
  expect_true(isTRUE(formals(RunscCODA)$reuse_reverse_comparisons))
  build_inputs <- getFromNamespace("build_sccoda_sample_inputs", "scop")
  dat <- data.frame(
    cluster = c("T", "B", "T", "B", "T", "B"),
    condition = c("A", "A", "B", "B", "B", "B"),
    sample = c("donor1", "donor2", "donor1", "donor2", "donor1", "donor2"),
    stringsAsFactors = FALSE
  )

  result <- build_inputs(dat)

  expect_equal(nrow(result$metadata), 4L)
  expect_equal(nrow(result$counts), 4L)
  expect_identical(rownames(result$counts), result$metadata$sample)
  condition_counts <- table(result$metadata$condition)
  expect_identical(names(condition_counts), c("A", "B"))
  expect_identical(as.integer(condition_counts), c(2L, 2L))
  expect_true(all(rowSums(result$counts) > 0))
  expect_equal(
    nrow(result$mapping[result$mapping$sample == "donor1", , drop = FALSE]),
    2L
  )
})

test_that("scCODA Python helpers accept structured comparisons and quoted groups", {
  python <- unname(Sys.which(c("python3", "python")))
  python <- python[nzchar(python)][1]
  skip_if(is.na(python), "Python is not available")
  runner <- getFromNamespace("runner_script_path", "scop")(
    "functions.py",
    "scCODA"
  )
  script <- tempfile(fileext = ".py")
  writeLines(
    c(
      "import runpy",
      "import sys",
      "import types",
      "log_module = types.ModuleType('log_message')",
      "log_module.log_message = lambda *args, **kwargs: None",
      "sys.modules['log_message'] = log_module",
      "module = runpy.run_path(sys.argv[1], run_name='sccoda_test')",
      "parse = module['_sccoda_parse_comparison']",
      "formula = module['_sccoda_formula']",
      "reverse = module['_sccoda_reverse_result_rows']",
      "assert parse([\"10x 3' v3\", \"10x 5' v1\"]) == (\"10x 3' v3\", \"10x 5' v1\")",
      "assert parse('A_vs_B') == ('A', 'B')",
      "value = formula('condition', \"10x 3' v3\")",
      "compile(value, '<sccoda-formula>', 'eval')",
      "rows = [{'obs_log2FD': 2.0, 'boot_mean_log2FD': 1.5, 'boot_CI_2.5': -1.0, 'boot_CI_97.5': 3.0, 'pval': 0.2, 'credible': True}]",
      "rev = reverse(rows)[0]",
      "assert rev['obs_log2FD'] == -2.0 and rev['boot_mean_log2FD'] == -1.5",
      "assert rev['boot_CI_2.5'] == -3.0 and rev['boot_CI_97.5'] == 1.0",
      "assert rev['pval'] == 0.2 and rev['credible'] is True",
      "print('sccoda-comparison-ok')"
    ),
    script
  )
  output <- suppressWarnings(system2(
    python,
    c(script, runner),
    stdout = TRUE,
    stderr = TRUE
  ))
  expect_identical(attr(output, "status") %||% 0L, 0L)
  expect_true(any(output == "sccoda-comparison-ok"))
})
