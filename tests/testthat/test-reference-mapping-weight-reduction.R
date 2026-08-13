test_that("RunReferenceMapping transfers labels in the ATAC reduction", {
  skip_if_not_installed("Signac")

  genes <- paste0("Gene", seq_len(6))
  query_cells <- paste0("Query", seq_len(4))
  reference_cells <- paste0("Reference", seq_len(5))
  query_counts <- Matrix::Matrix(
    matrix(
      seq_len(length(genes) * length(query_cells)),
      nrow = length(genes),
      dimnames = list(genes, query_cells)
    ),
    sparse = TRUE
  )
  reference_counts <- Matrix::Matrix(
    matrix(
      seq_len(length(genes) * length(reference_cells)),
      nrow = length(genes),
      dimnames = list(genes, reference_cells)
    ),
    sparse = TRUE
  )
  query <- Seurat::CreateSeuratObject(query_counts)
  peaks <- paste0("chr1-", seq(1, 61, by = 10), "-", seq(10, 70, by = 10))
  peak_counts <- Matrix::Matrix(
    matrix(1, nrow = length(peaks), ncol = length(query_cells),
      dimnames = list(peaks, query_cells)),
    sparse = TRUE
  )
  query[["peaks"]] <- Signac::CreateChromatinAssay(counts = peak_counts)
  SeuratObject::DefaultAssay(query) <- "peaks"
  query[["lsi"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      seq_len(length(query_cells) * 4),
      nrow = length(query_cells),
      dimnames = list(query_cells, paste0("LSI_", seq_len(4)))
    ),
    assay = "peaks",
    key = "LSI_"
  )

  reference <- Seurat::CreateSeuratObject(reference_counts)
  SeuratObject::VariableFeatures(reference) <- genes
  reference$cell_type <- rep(c("A", "B"), length.out = ncol(reference))
  reference[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      seq_len(length(reference_cells) * 2),
      nrow = length(reference_cells),
      dimnames = list(reference_cells, paste0("PC_", seq_len(2)))
    ),
    assay = "RNA",
    key = "PC_"
  )
  reference[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      seq_len(length(reference_cells) * 2),
      nrow = length(reference_cells),
      dimnames = list(reference_cells, paste0("UMAP_", seq_len(2)))
    ),
    assay = "RNA",
    key = "UMAP_"
  )

  transfer_args <- NULL
  testthat::local_mocked_bindings(
    FindTransferAnchors = function(...) list(),
    IntegrateEmbeddings = function(query, ...) {
      query[["ref.pca"]] <- SeuratObject::CreateDimReducObject(
        embeddings = matrix(
          seq_len(ncol(query) * 2),
          nrow = ncol(query),
          dimnames = list(colnames(query), paste0("REFPC_", seq_len(2)))
        ),
        assay = "peaks",
        key = "REFPC_"
      )
      query
    },
    .package = "Seurat"
  )
  testthat::local_mocked_bindings(
    atac_k_weight = function(...) 1L,
    atac_transfer_labels = function(
      srt, linear_reduction_dims_use, weight_reduction, ...
    ) {
      transfer_args <<- list(
        dims = linear_reduction_dims_use,
        reduction = weight_reduction
      )
      srt$predicted_predicted.id <- rep("A", ncol(srt))
      srt
    },
    RunKNNMap = function(srt_query, ...) srt_query,
    .package = "scop"
  )

  out <- RunReferenceMapping(
    srt = query,
    reference = reference,
    assay = "peaks",
    reference_assay = "RNA",
    reference_reduction = "pca",
    ref_umap = "umap",
    reference_dims = 1:2,
    reference_label = "cell_type",
    gene_activity_assay = "RNA",
    dims = 2:4,
    features = genes,
    verbose = FALSE
  )

  expect_s4_class(out, "Seurat")
  expect_identical(transfer_args$reduction, "lsi")
  expect_identical(transfer_args$dims, 2:4)
})
