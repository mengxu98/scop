test_that("Python SCENIC GRNBoost2 uses threaded local Dask cluster", {
  skip_if(
    identical(Sys.info()[["sysname"]], "Darwin"),
    "reticulate Python discovery is unstable on macOS CI"
  )
  skip_if_not_installed("reticulate")
  skip_if_not(reticulate::py_module_available("arboreto"))
  skip_if_not(reticulate::py_module_available("distributed"))

  functions <- reticulate::import_from_path(
    "functions",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )

  tmp_dir <- tempfile("scenic-python-grn-")
  dir.create(tmp_dir)
  expr_file <- file.path(tmp_dir, "expression.csv")
  regulators_file <- file.path(tmp_dir, "regulators.txt")
  adj_file <- file.path(tmp_dir, "adj.tsv")

  expr <- data.frame(
    TF1 = c(1, 2, 3, 4, 5, 6),
    TF2 = c(6, 5, 4, 3, 2, 1),
    Gene1 = c(1, 1, 2, 3, 5, 8),
    Gene2 = c(8, 5, 3, 2, 1, 1),
    Gene3 = c(0, 1, 0, 1, 0, 1),
    row.names = paste0("Cell", seq_len(6)),
    check.names = FALSE
  )
  utils::write.csv(expr, expr_file)
  writeLines(c("TF1", "TF2"), regulators_file)

  functions$RunSCENICGrn(
    expression_mtx = expr_file,
    regulators = regulators_file,
    adj_output = adj_file,
    n_rounds = 10L,
    early_stop_window_length = 5L,
    cores = 1L,
    seed = 1234L,
    force = TRUE,
    verbose = FALSE
  )

  expect_true(file.exists(adj_file))
  adj <- utils::read.delim(adj_file, stringsAsFactors = FALSE)
  expect_true(all(c("TF", "target", "importance") %in% colnames(adj)))
  expect_gt(nrow(adj), 0)
})
