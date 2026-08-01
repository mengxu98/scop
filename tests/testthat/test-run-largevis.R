test_that("RunLargeVis maps cores to the uwot n_threads interface", {
  skip_if_not_installed("uwot")
  skip_if_not_installed("Seurat")

  set.seed(1)
  x <- matrix(
    stats::rnorm(40 * 5),
    nrow = 40,
    dimnames = list(paste0("cell", 1:40), paste0("feature", 1:5))
  )
  out <- RunLargeVis(
    x,
    assay = "RNA",
    perplexity = 3,
    n_neighbors = 10,
    n_epochs = 2,
    cores = 1,
    verbose = FALSE
  )

  expect_s4_class(out, "DimReduc")
  expect_identical(rownames(Seurat::Embeddings(out)), rownames(x))
  expect_equal(dim(Seurat::Embeddings(out)), c(40L, 2L))
})
