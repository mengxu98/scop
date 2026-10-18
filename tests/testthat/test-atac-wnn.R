make_wnn_paired_object <- function() {
  skip_if_not_installed("Signac")
  withr::local_seed(42)
  cells <- paste0("cell", seq_len(80))
  rna <- Matrix::Matrix(matrix(rpois(300 * 80, 2), 300, 80,
    dimnames = list(paste0("gene", seq_len(300)), cells)
  ), sparse = TRUE)
  peaks <- Matrix::Matrix(matrix(rpois(200 * 80, 1), 200, 80,
    dimnames = list(paste0("chr1-", seq(1, by = 200, length.out = 200), "-", seq(100, by = 200, length.out = 200)), cells)
  ), sparse = TRUE)
  object <- SeuratObject::CreateSeuratObject(rna)
  object[["peaks"]] <- Signac::CreateChromatinAssay(peaks, verbose = FALSE)
  object
}

test_that("WNN selects RNA independently of derived default assays", {
  object <- make_wnn_paired_object()
  object[["chromvar"]] <- SeuratObject::CreateAssayObject(data = matrix(1, 3, 80,
    dimnames = list(paste0("motif", 1:3), colnames(object))
  ))
  SeuratObject::DefaultAssay(object) <- "chromvar"
  expect_identical(wnn_assays(object), list(rna = "RNA", atac = "peaks"))
  expect_error(wnn_assays(object, c("chromvar", "peaks")), "derived")
  object[["scores"]] <- object[["RNA"]]
  SeuratObject::Misc(object[["scores"]], "scop_gene_activity") <- list(source_assay = "peaks")
  expect_error(wnn_assays(object, c("scores", "peaks")), "derived")
  for (selection in list(character(), c("RNA", "missing"), c("RNA", "RNA"), c("RNA", NA_character_))) {
    expect_error(wnn_assays(object, selection), "[Aa]ssay for WNN")
  }
  object[["measured"]] <- object[["RNA"]]
  expect_identical(wnn_assays(object, c("measured", "peaks")), list(rna = "measured", atac = "peaks"))
  expect_error(wnn_assays(object, c("RNA", "measured")), "requires one RNA")
  object[["RNA"]] <- NULL
  expect_error(wnn_assays(object), "ambiguous")
})

test_that("WNN append and dispatch retain reusable modality results", {
  object <- make_wnn_paired_object()
  args <- list(
    srt_merge = object, assay = c("RNA", "peaks"),
    linear_reduction_dims = 10, linear_reduction_dims_use = 2:8,
    nHVF = 100, neighbor_k = 10, nonlinear_reduction_dims = 2, verbose = FALSE
  )
  withr::local_seed(109)
  before <- .Random.seed
  full <- suppressWarnings(do.call(WNN_integrate, c(args, list(append = FALSE))))
  expect_identical(.Random.seed, before)
  appended <- suppressWarnings(do.call(WNN_integrate, c(args, list(append = TRUE))))
  dispatched <- suppressWarnings(do.call(RunIntegration, c(
    args,
    list(append = TRUE, batch = NULL, integration_methods = "WNN")
  )))
  for (result in list(appended, dispatched)) {
    expect_true(all(result@misc$WNN_reduction_list %in% SeuratObject::Reductions(result)))
    expect_true(all(c("WNNKNN", "WNNSNN") %in% names(result@graphs)))
    expect_true("WNN" %in% names(result@neighbors))
    expect_identical(colnames(result), colnames(object))
    expect_equal(SeuratObject::Embeddings(result[["ATAClsi"]]), SeuratObject::Embeddings(full[["ATAClsi"]]))
    expect_equal(result$ATAC.weight, full$ATAC.weight)
    expect_equal(result$WNNclusters, full$WNNclusters)
    expect_equal(
      SeuratObject::LayerData(result[["peaks"]], layer = "data"),
      SeuratObject::LayerData(full[["peaks"]], layer = "data")
    )
  }
  # A supplied batch must be acknowledged before any analysis or input error.
  expect_warning(
    expect_error(WNN_integrate(batch = "sample", verbose = FALSE), "must be provided"),
    "does not perform batch correction"
  )
})
