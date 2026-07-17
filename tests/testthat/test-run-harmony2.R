test_that("RunHarmony2 writes matrix-level component standard deviations", {
  counts <- Matrix::Matrix(
    matrix(
      c(1, 2, 3, 4, 2, 3, 4, 5),
      nrow = 2,
      dimnames = list(c("Gene1", "Gene2"), paste0("Cell", 1:4))
    ),
    sparse = TRUE
  )
  srt <- Seurat::CreateSeuratObject(counts)
  srt$batch <- c("A", "A", "B", "B")
  pca_embeddings <- matrix(
    seq_len(8), ncol = 2,
    dimnames = list(colnames(srt), c("PC_1", "PC_2"))
  )
  srt[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = pca_embeddings,
    assay = "RNA",
    key = "PC_"
  )
  corrected <- matrix(
    c(1, 3, 2, 5, 4, 2, 6, 8),
    nrow = 2,
    dimnames = list(c("H1", "H2"), colnames(srt))
  )

  testthat::local_mocked_bindings(
    check_r = function(...) TRUE,
    get_namespace_fun = function(package, name) {
      expect_identical(package, "harmony")
      expect_identical(name, "RunHarmony")
      function(...) list(
        Z_corr = corrected,
        R = matrix(1, nrow = 1, ncol = ncol(srt))
      )
    },
    .package = "scop"
  )

  out <- RunHarmony2(
    srt,
    group.by.vars = "batch",
    reduction = "pca",
    dims.use = 1:2,
    project.dim = FALSE,
    verbose = FALSE
  )
  expected_embed <- t(corrected)

  expect_equal(
    unname(out[["Harmony"]]@stdev),
    unname(apply(expected_embed, 2, stats::sd)),
    tolerance = 1e-12
  )
})
