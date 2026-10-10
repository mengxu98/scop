# Backend-boundary test: real Seurat storage/identity operations, mocked RareQ fit.
test_that("RareQ named assignments retain exact identities through Seurat operations", {
  counts <- matrix(seq_len(24), nrow = 4,
    dimnames = list(paste0("gene", 1:4), paste0("cell", 1:6)))
  object <- Seurat::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object@neighbors$RNA.nn <- methods::new("Neighbor",
    nn.idx = matrix(rep(seq_len(6), 2), ncol = 2),
    nn.dist = matrix(0, nrow = 6, ncol = 2), cell.names = colnames(object))
  testthat::local_mocked_bindings(.package = "scop",
    check_r = function(...) invisible(TRUE),
    get_namespace_fun = function(package, fun) {
      expect_identical(package, "RareQ")
      switch(fun,
        ComputeQ = function(...) seq_len(6) / 10,
        FindRare = function(...) c("A", "B", "A", "C", "C", "B"),
        stop("Unexpected backend function"))
    })
  out <- RunRareQ(object, k = 2, run_neighbors = FALSE,
    cluster_colname = "custom_rare", tool_name = "custom_fit", verbose = FALSE)
  expect_type(out@tools$custom_fit$clusters, "character")
  expect_spatial_result_lifecycle(out, "custom_fit", "custom_rare")
  index <- out@tools$custom_fit$result_index
  broken <- out
  broken@meta.data[[index$cell_id_colname]][2] <- broken@meta.data[[index$cell_id_colname]][1]
  expect_error(GetSpatialResult(broken, "custom_fit"), "identities.*stale")
  legacy <- out
  legacy@tools$custom_fit$result_index <- NULL
  expect_equal(GetSpatialResult(legacy[, c("cell1", "cell3")], "custom_fit")$clusters,
    c(cell1 = "A", cell3 = "A"))
  legacy <- SeuratObject::RenameCells(legacy, new.names = rev(colnames(legacy)))
  expect_error(GetSpatialResult(legacy, "custom_fit"), "identities.*stale")
  # Another producer must not overwrite this fit's private identity field.
  expect_error(RunRareQ(out, k = 2, run_neighbors = FALSE, tool_name = "other_fit",
    q_colname = index$cell_id_colname, verbose = FALSE), "identity metadata column")
  # Auxiliary outputs may deliberately reuse our own old marker's name.
  # Even Q values that happen to equal the old cell IDs must stay numeric.
  numbered <- SeuratObject::RenameCells(object, new.names = as.character(seq_len(6) / 10))
  numbered <- RunRareQ(numbered, k = 2, run_neighbors = FALSE, verbose = FALSE)
  old_column <- numbered@tools$RareQ$result_index$cell_id_colname
  numbered <- SeuratObject::RenameCells(numbered, new.names = paste0("new", seq_len(6)))
  rerun <- RunRareQ(numbered, k = 2, run_neighbors = FALSE,
    q_colname = old_column, verbose = FALSE)
  expect_identical(rerun@meta.data[[old_column]], seq_len(6) / 10)
  expect_false(identical(rerun@tools$RareQ$result_index$cell_id_colname, old_column))
  expect_identical(names(GetSpatialResult(rerun, "RareQ")$clusters), colnames(rerun))
})
