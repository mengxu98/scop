test_that("RunRareQ names rebuilt neighbors for the requested assay", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("RareQ")

  set.seed(1)
  counts <- Matrix::Matrix(
    rpois(4000, lambda = 5),
    nrow = 100,
    sparse = TRUE
  )
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))

  srt <- Seurat::CreateSeuratObject(counts = counts, assay = "RNA")
  astrocyte_assay <- SeuratObject::CreateAssay5Object(counts = counts)
  SeuratObject::Key(astrocyte_assay) <- "Astrocyte_"
  srt[["Astrocyte"]] <- astrocyte_assay
  srt <- Seurat::NormalizeData(srt, assay = "Astrocyte", verbose = FALSE)
  srt <- Seurat::FindVariableFeatures(
    srt,
    assay = "Astrocyte",
    verbose = FALSE
  )
  srt <- Seurat::ScaleData(srt, assay = "Astrocyte", verbose = FALSE)
  SeuratObject::DefaultAssay(srt) <- "Astrocyte"
  srt <- Seurat::RunPCA(
    srt,
    assay = "Astrocyte",
    npcs = 10,
    reduction.name = "Astrocytepca",
    verbose = FALSE
  )
  SeuratObject::DefaultAssay(srt) <- "RNA"

  srt <- RunRareQ(
    srt = srt,
    assay = "RNA",
    reduction = "Astrocytepca",
    dims = 1:5,
    k.param = 8,
    k = 3,
    run_neighbors = TRUE,
    force_recalc = TRUE,
    prefix = "AstrocyteRareQ",
    verbose = FALSE
  )

  expect_true("RNA.nn" %in% names(srt@neighbors))
  expect_false("Astrocyte.nn" %in% names(srt@neighbors))
  expect_true(all(
    c(
      "AstrocyteRareQ_cluster",
      "AstrocyteRareQ_Q",
      "AstrocyteRareQ_cluster_size",
      "AstrocyteRareQ_is_rare"
    ) %in% colnames(srt@meta.data)
  ))
  expect_identical(srt@tools$RareQ$neighbor_slot, "RNA.nn")
})
