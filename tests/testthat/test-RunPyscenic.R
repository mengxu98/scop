test_that("pySCENIC environment requirements use the isolated stack", {
  req <- env_requirements(version = "3.10-1", modules = "pyscenic")
  req_py311 <- env_requirements(version = "3.11-1", modules = "pyscenic")

  expect_equal(req_py311$python, "3.10-1")
  expect_equal(req$packages[["numpy"]], "numpy==1.23.5")
  expect_equal(req$packages[["pyscenic"]], "pyscenic==0.12.1")
  expect_equal(req$packages[["arboreto"]], "arboreto==0.1.6")
  expect_equal(req$packages[["ctxcore"]], "ctxcore==0.2.0")
  expect_equal(req$packages[["dask"]], "dask==2024.2.1")
  expect_equal(req$packages[["distributed"]], "distributed==2024.2.1")
  expect_true("pyarrow" %in% names(req$packages))
  expect_false("scanpy" %in% names(req$packages))
  expect_false("leidenalg" %in% names(req$packages))
  expect_false("pyscenic" %in% scop:::default_env_modules())
  expect_error(
    env_requirements(modules = c("pyscenic", "scanpy")),
    "standalone"
  )
})

test_that("GRN input matrix is cells by genes after filtering", {
  counts <- Matrix::Matrix(
    c(
      1, 0, 2,
      0, 0, 0,
      3, 1, 0,
      0, 4, 5
    ),
    nrow = 4,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(counts) <- paste0("g", 1:4)
  colnames(counts) <- paste0("cell", 1:3)

  grn <- scop:::pyscenic_prepare_grn_matrix(
    counts = counts,
    ranking_genes = c("g1", "g3", "g4"),
    min_expr_cells = 2,
    verbose = FALSE
  )

  expect_equal(rownames(grn), colnames(counts))
  expect_equal(colnames(grn), c("g1", "g3", "g4"))
  expect_equal(dim(grn), c(3, 3))
})

test_that("default metacell targets scale by cell count", {
  expect_equal(scop:::pyscenic_default_target_metacells(500), 80L)
  expect_equal(scop:::pyscenic_default_target_metacells(1000), 100L)
  expect_equal(scop:::pyscenic_default_target_metacells(10000), 500L)
  expect_equal(scop:::pyscenic_default_target_metacells(100000), 1500L)
  expect_equal(scop:::pyscenic_default_target_metacells(300000), 3000L)
  expect_equal(scop:::pyscenic_default_target_metacells(1000000), 10000L)
  expect_equal(scop:::pyscenic_default_target_metacells(2000000), 10000L)
})

test_that("metacell resolution selection favors the closest and larger tied result", {
  resolution_summary <- data.frame(
    resolution = c(5, 10, 20),
    cluster_col = c("RNA_snn_res.5", "RNA_snn_res.10", "RNA_snn_res.20"),
    n_metacells = c(80, 120, 70),
    stringsAsFactors = FALSE
  )

  selected <- scop:::pyscenic_select_metacell_resolution(
    resolution_summary = resolution_summary,
    target_metacells = 100
  )

  expect_equal(selected$selected_resolution, 10)
  expect_equal(selected$selected_cluster_col, "RNA_snn_res.10")
  expect_equal(selected$resolution_summary$n_metacells[[1]], 120)
  expect_true("target_distance" %in% colnames(selected$resolution_summary))
})

test_that("metacell resolution summary counts grouped metacells", {
  counts <- Matrix::Matrix(
    matrix(1, nrow = 5, ncol = 6),
    sparse = TRUE
  )
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(
    counts = counts
  )
  srt@meta.data[["sample"]] <- c("A", "A", "B", "B", "B", "B")
  srt@meta.data[["RNA_snn_res.1"]] <- c(0, 0, 0, 0, 1, 1)
  srt@meta.data[["RNA_snn_res.2"]] <- c(0, 1, 0, 1, 0, 1)

  resolution_summary <- scop:::pyscenic_resolution_summary(
    srt = srt,
    cluster_cols = c("RNA_snn_res.1", "RNA_snn_res.2"),
    resolution_candidates = c(1, 2),
    group_df = srt@meta.data[, "sample", drop = FALSE]
  )

  expect_equal(resolution_summary$n_metacells, c(3L, 4L))
  expect_equal(resolution_summary$cluster_col, c("RNA_snn_res.1", "RNA_snn_res.2"))
})

test_that("RunPyscenic exposes default metacell reduction", {
  expect_equal(formals(RunPyscenic)$metacell_reduction, "pca")
})

test_that("metacell builder records the default pca reduction", {
  set.seed(1)
  counts <- Matrix::Matrix(
    matrix(rpois(80 * 25, lambda = 3), nrow = 80, ncol = 25),
    sparse = TRUE
  )
  rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)

  suppressWarnings(
    metacell_result <- scop:::pyscenic_build_metacell_counts(
      srt = srt,
      counts = counts,
      assay = "RNA",
      metacell_resolution = 0.1,
      metacell_dims = 1:30,
      verbose = FALSE
    )
  )

  expect_equal(metacell_result$info$reduction, "pca")
  expect_equal(metacell_result$info$dims, 1:24)
  expect_equal(rownames(metacell_result$counts), rownames(counts))
})

test_that("metacell builder rejects missing grouping annotations", {
  counts <- Matrix::Matrix(
    matrix(rpois(60 * 25, lambda = 3), nrow = 60, ncol = 25),
    sparse = TRUE
  )
  rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt@meta.data[["sample"]] <- rep(c("A", "B"), length.out = ncol(srt))
  srt@meta.data[["sample"]][[1]] <- NA_character_

  expect_error(
    scop:::pyscenic_build_metacell_counts(
      srt = srt,
      counts = counts,
      assay = "RNA",
      metacell.by = "sample",
      metacell_resolution = 0.1,
      metacell_dims = 1:3,
      verbose = FALSE
    ),
    "missing values"
  )
})

test_that("metacell builder validates custom reduction names and dimensions", {
  counts <- Matrix::Matrix(
    matrix(rpois(60 * 25, lambda = 3), nrow = 60, ncol = 25),
    sparse = TRUE
  )
  rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)

  expect_error(
    scop:::pyscenic_build_metacell_counts(
      srt = srt,
      counts = counts,
      assay = "RNA",
      metacell_reduction = "Harmony",
      metacell_resolution = 0.1,
      metacell_dims = 1:3,
      verbose = FALSE
    ),
    "metacell_reduction.*Harmony.*not present"
  )

  harmony_embeddings <- matrix(rnorm(ncol(srt) * 3), nrow = ncol(srt), ncol = 3)
  rownames(harmony_embeddings) <- colnames(srt)
  colnames(harmony_embeddings) <- paste0("Harmony_", 1:3)
  srt[["Harmony"]] <- Seurat::CreateDimReducObject(
    embeddings = harmony_embeddings,
    assay = "RNA",
    key = "Harmony_"
  )

  expect_error(
    scop:::pyscenic_build_metacell_counts(
      srt = srt,
      counts = counts,
      assay = "RNA",
      metacell_reduction = "Harmony",
      metacell_resolution = 0.1,
      metacell_dims = 1:4,
      verbose = FALSE
    ),
    "metacell_dims.*4.*Harmony.*3"
  )

  metacell_result <- scop:::pyscenic_build_metacell_counts(
    srt = srt,
    counts = counts,
    assay = "RNA",
    metacell_reduction = "Harmony",
    metacell_resolution = 0.1,
    metacell_dims = 1:3,
    verbose = FALSE
  )

  expect_equal(metacell_result$info$reduction, "Harmony")
  expect_equal(metacell_result$info$dims, 1:3)
  expect_equal(rownames(metacell_result$counts), rownames(counts))
})

test_that("regulon text table is converted to AUCell gene sets", {
  txt_file <- tempfile(fileext = ".txt")
  writeLines(
    c(
      "BCLAF1(12g)\thttp://motif/logo.png\tGeneA,GeneB,GeneC",
      "IRF9(10g)\thttp://motif/logo2.png\tGeneD,GeneE"
    ),
    txt_file
  )

  regulon_tbl <- scop:::pyscenic_read_regulon_txt(txt_file)
  regulon_list <- scop:::pyscenic_regulon_list(regulon_tbl)

  expect_equal(names(regulon_list), c("BCLAF1(+)", "IRF9(+)"))
  expect_equal(regulon_list[["BCLAF1(+)"]], c("GeneA", "GeneB", "GeneC"))
})

test_that("AUCell scoring passes batch cores and preserves cells and regulons", {
  counts <- Matrix::Matrix(
    matrix(
      sample(0:5, 120, replace = TRUE),
      nrow = 30,
      ncol = 4
    ),
    sparse = TRUE
  )
  rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))

  regulons <- list(
    TF1 = paste0("Gene", 1:10),
    TF2 = paste0("Gene", 11:20)
  )

  testthat::local_mocked_bindings(
    check_r = function(...) invisible(TRUE),
    pyscenic_calc_auc_batch = function(counts, regulon_list) {
      out <- matrix(
        seq_len(ncol(counts) * length(regulon_list)),
        nrow = ncol(counts),
        ncol = length(regulon_list)
      )
      rownames(out) <- colnames(counts)
      colnames(out) <- names(regulon_list)
      out
    },
    .package = "scop"
  )

  scores <- scop:::pyscenic_compute_aucell_score(
    counts = counts,
    regulon_list = regulons,
    min_regulon_size = 5,
    batch_size = 2,
    cores = if (.Platform$OS.type == "unix") 2 else 1,
    verbose = FALSE
  )

  expect_equal(rownames(scores), colnames(counts))
  expect_equal(colnames(scores), names(regulons))
})
