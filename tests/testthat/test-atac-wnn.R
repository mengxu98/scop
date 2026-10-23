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
  object[["chromvar"]] <- SeuratObject::CreateAssayObject(data = matrix(1, 2, ncol(object),
    dimnames = list(c("m1", "m2"), colnames(object))
  ))
  SeuratObject::DefaultAssay(object) <- "chromvar"
  args$srt_merge <- object
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
    expect_identical(SeuratObject::DefaultAssay(result), "RNA")
    expect_equal(result[["RNA"]], full[["RNA"]])
    expect_equal(result[["peaks"]], full[["peaks"]])
    expect_equal(result[["chromvar"]], object[["chromvar"]])
    expect_true(all(result@misc$WNN_reduction_list %in% SeuratObject::Reductions(result)))
    expect_true(all(c("WNNKNN", "WNNSNN") %in% names(result@graphs)))
    expect_true("WNN" %in% names(result@neighbors))
    expect_identical(colnames(result), colnames(object))
    expect_equal(SeuratObject::Embeddings(result[["ATAClsi"]]), SeuratObject::Embeddings(full[["ATAClsi"]]))
    expect_equal(result$ATAC.weight, full$ATAC.weight)
    expect_equal(result$WNNclusters, full$WNNclusters)
    receipt <- SeuratObject::Misc(result[["ATAClsi"]], "scop_atac_representation")
    expect_identical(receipt$dimensions, 2:8)
    expect_equal(receipt, SeuratObject::Misc(full[["ATAClsi"]], "scop_atac_representation"))
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

test_that("WNN requires the same measured cells in both modalities", {
  object <- make_wnn_paired_object()
  object[["peaks"]] <- subset(object[["peaks"]], cells = colnames(object)[-1])
  expect_error(WNN_integrate(srt_merge = object, verbose = FALSE), "same.*cells|paired.*cells")
})

test_that("WNN honors graph-only requests and rejects unsupported embeddings", {
  object <- make_wnn_paired_object()
  expect_error(WNN_integrate(srt_merge = object, nonlinear_reduction = "tsne", verbose = FALSE), "supports.*umap")
  result <- WNN_integrate(
    srt_merge = object,
    linear_reduction_dims = 6, linear_reduction_dims_use = 2:5,
    nHVF = 80, neighbor_k = 10, nonlinear_reduction = character(), verbose = FALSE
  )
  expect_true(all(c("WNNKNN", "WNNSNN") %in% names(result@graphs)))
  expect_true("WNNclusters" %in% colnames(result[[]]))
  expect_false(any(grepl("^WNN(UMAP|FR)", SeuratObject::Reductions(result))))
})

test_that("WNN normalizes RNA and chromatin independently", {
  object <- make_wnn_paired_object()
  before <- SeuratObject::LayerData(object[["peaks"]], layer = "counts")
  result <- WNN_integrate(
    srt_merge = object, normalization_method = "SCT",
    linear_reduction_dims = 6, linear_reduction_dims_use = 2:5,
    nHVF = 80, neighbor_k = 10, nonlinear_reduction = character(), verbose = FALSE
  )
  expect_s4_class(result[["SCT"]], "SCTAssay")
  expect_identical(SeuratObject::DefaultAssay(result), "SCT")
  expect_identical(SeuratObject::DefaultAssay(result[["RNApca"]]), "SCT")
  native <- Signac::RunTFIDF(object[["peaks"]], verbose = FALSE)
  expect_equal(
    SeuratObject::LayerData(result[["peaks"]], layer = "data"),
    SeuratObject::LayerData(native, layer = "data")
  )
  expect_identical(SeuratObject::LayerData(result[["peaks"]], layer = "counts"), before)
  dispatched <- RunIntegration(srt_merge = object, batch = NULL,
    integration_methods = "WNN", assay = c("RNA", "peaks"), normalization_method = "SCT",
    linear_reduction_dims = 6, linear_reduction_dims_use = 2:5,
    nHVF = 80, neighbor_k = 10, nonlinear_reduction = character(), verbose = FALSE)
  expect_true(all(dispatched@misc$WNN_reduction_list %in% SeuratObject::Reductions(dispatched)))
  expect_equal(SeuratObject::Embeddings(dispatched[["RNApca"]]), SeuratObject::Embeddings(result[["RNApca"]]))
  expect_equal(dispatched$RNA.weight, result$RNA.weight)
  expect_equal(dispatched[["SCT"]], result[["SCT"]])
  selected <- WNN_integrate(
    srt_merge = object, HVF = rownames(object[["RNA"]])[1:80],
    linear_reduction_dims = 6, linear_reduction_dims_use = 2:5,
    nHVF = 80, neighbor_k = 10, nonlinear_reduction = character(), verbose = FALSE
  )
  expect_setequal(rownames(SeuratObject::Loadings(selected[["RNApca"]])), rownames(object[["RNA"]])[1:80])
  expect_true(all(rownames(SeuratObject::Loadings(selected[["ATAClsi"]])) %in% rownames(object[["peaks"]])))
})

test_that("WNN clustering failures cannot return stale cluster labels", {
  object <- make_wnn_paired_object()
  object$WNNclusters <- "previous"
  original <- FindClusters
  testthat::local_mocked_bindings(FindClusters = function(object, ..., graph.name) {
    if (identical(graph.name, "WNNSNN")) stop("WNN clustering failed")
    original(object, ..., graph.name = graph.name)
  }, .package = "scop")
  expect_error(WNN_integrate(
    srt_merge = object,
    linear_reduction_dims = 6, linear_reduction_dims_use = 2:5,
    nHVF = 80, neighbor_k = 10, nonlinear_reduction = character(), verbose = FALSE
  ), "WNN clustering failed")
  expect_true(all(object$WNNclusters == "previous"))
})
