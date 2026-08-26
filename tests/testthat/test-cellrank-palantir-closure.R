make_cellrank_closure_srt <- function() {
  counts <- Matrix::Diagonal(4, x = rep(1, 4))
  rownames(counts) <- paste0("g", 1:4)
  colnames(counts) <- paste0("c", 1:4)
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt@meta.data$palantir_pseudotime <- seq(0, 1, length.out = 4)
  srt@tools$CellRank <- list(
    fate_probabilities = matrix(
      c(.9, .1, .8, .2, .2, .8, .1, .9),
      ncol = 2,
      dimnames = list(colnames(srt), c("A", "B"))
    ),
    lineage_drivers = data.frame(
      A_corr = c(.9, .4, -.2, .1),
      row.names = rownames(srt)
    )
  )
  srt
}

test_that("trajectory dependency profiles follow the requested Python version", {
  legacy <- env_requirements(
    version = "3.10-1",
    modules = c("scanpy", "palantir", "cellrank")
  )
  middle <- env_requirements(
    version = "3.11-1",
    modules = c("scanpy", "palantir", "cellrank")
  )
  modern <- env_requirements(
    version = "3.12-1",
    modules = c("scanpy", "palantir", "cellrank")
  )
  expect_identical(legacy$packages[["scanpy"]], "scanpy==1.11.3")
  expect_identical(legacy$packages[["cellrank"]], "cellrank==2.0.7")
  expect_identical(legacy$packages[["numpy"]], "numpy==1.26.4")
  expect_identical(middle$packages[["scanpy"]], "scanpy==1.11.3")
  expect_identical(middle$packages[["cellrank"]], "cellrank==2.0.7")
  expect_identical(modern$packages[["scanpy"]], "scanpy==1.12.2")
  expect_identical(modern$packages[["cellrank"]], "cellrank==2.3.2")
  expect_identical(modern$packages[["numpy"]], "numpy==2.4.6")
})

test_that("package checks use the selected environment interpreter profile", {
  testthat::local_mocked_bindings(
    .package = "scop",
    conda_python = function(...) "C:/fake/python.exe",
    python_minor_version = function(...) "3.12",
    resolve_conda = function(x) "C:/fake/conda.exe",
    is_windows = function() TRUE
  )
  req <- unique_python_requirements(
    "cellrank",
    envname = "scop_env",
    conda = "C:/fake/conda.exe"
  )
  expect_true("cellrank==2.3.2" %in% req)
})

test_that("RunCellRankTrends stores normalized trend payload", {
  srt <- make_cellrank_closure_srt()
  testthat::local_mocked_bindings(
    .package = "scop",
    PrepareEnv = function(...) NULL,
    check_python = function(...) TRUE,
    srt_to_adata = function(...) "adata",
    scop_python_import = function(...) {
      list(CellRankTrends = function(...) {
        list(
          trend_matrix = matrix(1:6, nrow = 2),
          genes = c("g1", "g2"),
          time = c(0, 1, 2),
          clusters = c("0", "1"),
          parameters = list(n_points = 3L)
        )
      })
    }
  )
  out <- RunCellRankTrends(
    srt,
    lineage = "A",
    top_n = 2,
    n_points = 3,
    min_expressed_cells = 1L,
    verbose = FALSE
  )
  expect_true("A" %in% names(out@tools$CellRank$trends))
  expect_equal(out@tools$CellRank$trends$A$cluster_table$cluster, c("0", "1"))
  expect_equal(out@tools$CellRank$trends$A$heatmap_genes, c("g1", "g2"))
})

test_that("RunCellRankTrends rejects fate matrices without cell names", {
  srt <- make_cellrank_closure_srt()
  rownames(srt@tools$CellRank$fate_probabilities) <- NULL
  expect_error(
    RunCellRankTrends(srt, lineage = "A", verbose = FALSE),
    "unique cell row names"
  )
})

test_that("RunCellRankEnrichment records empty and non-empty module results", {
  srt <- make_cellrank_closure_srt()
  srt@tools$CellRank$trends <- list(A = list(
    cluster_table = data.frame(gene = c("g1", "g2"), cluster = c("0", "1"))
  ))
  testthat::local_mocked_bindings(
    .package = "scop",
    check_r = function(...) TRUE,
    get_namespace_fun = function(pkg, fun) {
      function(...) {
        args <- list(...)
        expect_identical(names(args$TERM2GENE), c("term", "gene"))
        expect_true(all(args$TERM2GENE$gene %in% c("g1", "g2")))
        data.frame(ID = "t1", Description = "term", p.adjust = 0.01, Count = 2)
      }
    },
    PrepareDB = function(...) list(Mus_musculus = list(
      TEST = list(
        TERM2GENE = data.frame(
          term = "t1",
          entrez_id = c("1", "2"),
          symbol = c("g1", "g2")
        ),
        TERM2NAME = data.frame(term = "t1", name = "term")
      )
    ))
  )
  out <- RunCellRankEnrichment(srt, lineage = "A", db = "TEST", verbose = FALSE)
  expect_true("A" %in% names(out@tools$CellRank$enrichment))
  expect_true("TEST" %in% names(out@tools$CellRank$enrichment$A$databases))
  expect_true(nrow(out@tools$CellRank$enrichment$A$databases$TEST) >= 1)
  expect_s3_class(
    CellRankPlot(out, plot_type = "enrichment", lineage = "A", database = "TEST"),
    "ggplot"
  )
})

test_that("CellRankPlot circular consumes stored fate probabilities", {
  srt <- make_cellrank_closure_srt()
  p <- CellRankPlot(srt, plot_type = "circular")
  expect_s3_class(p, "ggplot")
})

test_that("CellRankPlot drivers consumes matrix driver payloads", {
  srt <- make_cellrank_closure_srt()
  p <- CellRankPlot(srt, plot_type = "drivers", lineage = "A", top_n = 2)
  expect_s3_class(p, "ggplot")
})

test_that("CellRankPlot projection consumes the stored transition matrix", {
  srt <- make_cellrank_closure_srt()
  coords <- matrix(
    c(0, 0, 1, 0, 1, 1, 0, 1),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(colnames(srt), c("UMAP_1", "UMAP_2"))
  )
  srt[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = coords,
    key = "UMAP_",
    assay = "RNA"
  )
  transition <- Matrix::Diagonal(4, x = 1)
  rownames(transition) <- colnames(srt)
  colnames(transition) <- colnames(srt)
  srt@tools$CellRank$transition_matrix <- transition
  p <- CellRankPlot(srt, plot_type = "projection", reduction = "umap")
  expect_s3_class(p, "ggplot")
  walks <- CellRankPlot(
    srt,
    plot_type = "random_walks",
    reduction = "umap",
    n_sims = 2L,
    max_iter = 3L
  )
  expect_s3_class(walks, "ggplot")
})

test_that("RunCellRankEnrichment distinguishes backend failures", {
  srt <- make_cellrank_closure_srt()
  srt@tools$CellRank$trends <- list(A = list(
    cluster_table = data.frame(gene = c("g1", "g2"), cluster = c("0", "1"))
  ))
  testthat::local_mocked_bindings(
    .package = "scop",
    check_r = function(...) TRUE,
    get_namespace_fun = function(pkg, fun) {
      function(...) stop("malformed TERM2GENE")
    },
    PrepareDB = function(...) list(Mus_musculus = list(
      TEST = list(TERM2GENE = data.frame(term = "t1", gene = "g1"), TERM2NAME = data.frame(term = "t1", name = "term"))
    ))
  )
  expect_error(
    RunCellRankEnrichment(srt, lineage = "A", db = "TEST", verbose = FALSE),
    "malformed TERM2GENE"
  )
})

test_that("C++ backend rejects CellRank reference-only options", {
  expect_error(
    RunCellRank(
      backend = "cpp",
      allow_approximate = TRUE,
      terminal_states = "A",
      show_plot = FALSE,
      verbose = FALSE
    ),
    "terminal_states"
  )
})

test_that("RunCellRank stores named fate and driver payloads", {
  srt <- make_cellrank_closure_srt()
  srt[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(runif(8), nrow = 4, dimnames = list(colnames(srt), c("PC_1", "PC_2"))),
    assay = "RNA", key = "PC_"
  )
  srt$cluster <- factor(c("a", "a", "b", "b"))
  adata_out <- srt
  adata_out$macrostates_fwd <- factor(c("A", "A", "B", "B"))
  adata_out$term_states_fwd <- factor(c("A", "transient", "transient", "B"))
  adata_out$term_states_fwd_probs <- c(.9, .4, .3, .8)
  payload <- list(
    fate_probabilities = list(
      values = srt@tools$CellRank$fate_probabilities,
      index = colnames(srt),
      columns = c("A", "B")
    ),
    lineage_drivers = list(
      values = matrix(c(.8, .1, .2, .05), nrow = 2, dimnames = list(c("g1", "g2"), c("A_corr", "A_qval"))),
      index = c("g1", "g2"),
      columns = c("A_corr", "A_qval")
    ),
    transition_key = "cellrank_transition",
    versions = list(cellrank = "2.3.2")
  )
  captured <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    .package = "scop",
    PrepareEnv = function(...) NULL,
    check_python = function(...) TRUE,
    py_to_r2 = function(x) x,
    srt_to_adata = function(...) list(obs = data.frame(cluster = factor(c("a", "a", "b", "b")))),
    palette_colors = function(...) c(a = "#111111", b = "#222222"),
    scop_python_import = function(...) list(CellRank = function(...) {
      captured$recompute_neighbors <- list(...)$recompute_neighbors
      list("adata", "estimator", "kernel", payload)
    }),
    adata_to_srt = function(...) adata_out
  )
  out <- RunCellRank(
    srt = srt,
    group.by = "cluster",
    linear_reduction = "pca",
    nonlinear_reduction = "pca",
    backend = "python",
    return_seurat = TRUE,
    recompute_neighbors = FALSE,
    verbose = FALSE
  )
  expect_equal(colnames(out@tools$CellRank$fate_probabilities), c("A", "B"))
  expect_true(all(c("cellrank_fate_A", "cellrank_fate_B") %in% colnames(out@meta.data)))
  expect_equal(out@tools$CellRank$versions$cellrank, "2.3.2")
  expect_false(captured$recompute_neighbors)
})

test_that("CellRank trend fallback uses the CellRank normal distribution name", {
  expect_identical(formals(RunCellRankTrends)$fallback_distribution, "normal")
})
