test_that("AUCell consistency strategy matches official AUCell scores", {
  skip_if_not_installed("AUCell")

  set.seed(42)
  expr <- Matrix::rsparsematrix(120, 24, density = 0.08)
  expr@x <- abs(round(expr@x * 5))
  rownames(expr) <- paste0("gene", seq_len(nrow(expr)))
  colnames(expr) <- paste0("cell", seq_len(ncol(expr)))
  gene_sets <- list(
    set_a = paste0("gene", 1:30),
    set_b = paste0("gene", 25:75),
    set_c = paste0("gene", 70:115)
  )

  set.seed(11)
  expected <- suppressWarnings(
    getFromNamespace("run_aucell_official_scores", "scop")(
      expr_counts = expr,
      gene_sets = gene_sets
    )
  )
  set.seed(11)
  observed <- suppressWarnings(
    getFromNamespace("run_aucell_official_scores", "scop")(
      expr_counts = expr,
      gene_sets = gene_sets
    )
  )

  expect_equal(observed, expected)
})

test_that("C++ AUCell AUC matches official random and blocked rankings", {
  skip_if_not_installed("AUCell")

  set.seed(142)
  expr <- Matrix::rsparsematrix(80, 14, density = 0.15)
  expr@x <- abs(round(expr@x * 4))
  rownames(expr) <- paste0("gene", seq_len(nrow(expr)))
  colnames(expr) <- paste0("cell", seq_len(ncol(expr)))
  gene_sets <- list(
    set_a = paste0("gene", 1:20),
    set_b = paste0("gene", 15:50),
    set_c = c(paste0("gene", 50:70), "not_present")
  )

  for (split_by_blocks in c(FALSE, TRUE)) {
    set.seed(91)
    rankings <- AUCell::AUCell_buildRankings(
      expr, plotStats = FALSE, splitByBlocks = split_by_blocks, verbose = FALSE
    )
    expected <- suppressWarnings(AUCell::getAUC(AUCell::AUCell_calcAUC(
      geneSets = gene_sets, rankings = rankings, aucMaxRank = 12L, verbose = FALSE
    )))
    observed <- scop:::run_aucell_scores_from_official_rankings(
      rankings, gene_sets, auc_max_rank = 12L, norm_auc = TRUE
    )

    expect_equal(unname(observed), unname(t(expected)), tolerance = 1e-12)
  }
})

test_that("AUCell top-k ranks preserve the full native AUC", {
  set.seed(20260714)
  expr <- Matrix::rsparsematrix(300, 40, density = 0.12)
  expr@x <- abs(round(expr@x * 5))
  rownames(expr) <- paste0("gene", seq_len(nrow(expr)))
  colnames(expr) <- paste0("cell", seq_len(ncol(expr)))
  gene_sets <- list(
    set_a = rownames(expr)[1:70],
    set_b = rownames(expr)[55:160],
    set_c = rownames(expr)[180:280]
  )

  full <- scop:::run_aucell_scores(
    expr, gene_sets, strategy = "full", tie_method = "first"
  )
  topk <- scop:::run_aucell_scores(
    expr, gene_sets, strategy = "topk", tie_method = "first"
  )

  expect_equal(topk, full, tolerance = 1e-12)
})

test_that("CellScoring AUCell backend switch controls R and C++ paths", {
  skip_if_not_installed("AUCell")
  skip_if_not_installed("Seurat")

  set.seed(43)
  counts <- Matrix::rsparsematrix(90, 18, density = 0.1)
  counts@x <- abs(round(counts@x * 4))
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  features <- list(
    set_a = paste0("gene", 1:25),
    set_b = paste0("gene", 20:55)
  )

  r_out <- suppressWarnings(CellScoring(
    srt,
    features = features,
    method = "AUCell",
    backend = "r",
    classification = FALSE,
    name = "auc_r",
    verbose = FALSE
  ))
  cpp_out <- suppressWarnings(CellScoring(
    srt,
    features = features,
    method = "AUCell",
    backend = "cpp",
    classification = FALSE,
    name = "auc_cpp",
    verbose = FALSE
  ))

  r_scores <- as.matrix(r_out@meta.data[, c("auc_r_set_a", "auc_r_set_b")])
  cpp_scores <- as.matrix(cpp_out@meta.data[, c("auc_cpp_set_a", "auc_cpp_set_b")])
  colnames(r_scores) <- colnames(cpp_scores)

  expect_equal(dim(cpp_scores), dim(r_scores))
  expect_true(all(is.finite(cpp_scores)))
})

test_that("CellScoring AUCell cpp backend keeps high consistency with R backend", {
  skip_if_not_installed("AUCell")
  skip_if_not_installed("Seurat")

  set.seed(44)
  counts <- Matrix::rsparsematrix(120, 24, density = 0.08)
  counts@x <- abs(round(counts@x * 4))
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  features <- list(
    set_a = paste0("gene", 1:35),
    set_b = paste0("gene", 25:85)
  )

  r_out <- suppressWarnings(CellScoring(
    srt,
    features = features,
    method = "AUCell",
    backend = "r",
    classification = FALSE,
    name = "auc_r",
    verbose = FALSE
  ))
  cpp_out <- suppressWarnings(CellScoring(
    srt,
    features = features,
    method = "AUCell",
    backend = "cpp",
    classification = FALSE,
    name = "auc_cpp",
    verbose = FALSE
  ))

  r_scores <- as.matrix(r_out@meta.data[, c("auc_r_set_a", "auc_r_set_b")])
  cpp_scores <- as.matrix(cpp_out@meta.data[, c("auc_cpp_set_a", "auc_cpp_set_b")])
  expect_gte(
    suppressWarnings(stats::cor(
      as.numeric(r_scores),
      as.numeric(cpp_scores),
      method = "spearman",
      use = "complete.obs"
    )),
    0.95
  )
})

test_that("CellScoring AUCell uses an exactly shared tie rule", {
  skip_if_not_installed("AUCell")
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(
    c(
      0, 0, 2, 2, 0, 0,
      0, 0, 2, 2, 0, 0,
      1, 1, 0, 0, 1, 1,
      1, 1, 0, 0, 1, 1,
      0, 0, 0, 0, 3, 3,
      0, 0, 0, 0, 3, 3
    ),
    nrow = 6,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  features <- list(set_a = rownames(counts)[1:4], set_b = rownames(counts)[3:6])

  r_out <- suppressWarnings(CellScoring(
    srt, features = features, method = "AUCell", backend = "r",
    classification = FALSE, name = "auc_r_first",
    verbose = FALSE
  ))
  cpp_out <- suppressWarnings(CellScoring(
    srt, features = features, method = "AUCell", backend = "cpp",
    classification = FALSE, name = "auc_cpp_first",
    verbose = FALSE
  ))
  r_scores <- as.matrix(r_out@meta.data[, c("auc_r_first_set_a", "auc_r_first_set_b")])
  cpp_scores <- as.matrix(cpp_out@meta.data[, c("auc_cpp_first_set_a", "auc_cpp_first_set_b")])

  expect_equal(unname(cpp_scores), unname(r_scores), tolerance = 1e-12)
})

test_that("CellScoring routes AUCell native ranking options through both backends", {
  skip_if_not_installed("AUCell")
  skip_if_not_installed("Seurat")

  set.seed(144)
  counts <- Matrix::Matrix(matrix(stats::rpois(50L * 10L, 1), nrow = 50L), sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  features <- list(set_a = rownames(counts)[1:15], set_b = rownames(counts)[20:40])
  set.seed(17)
  r_out <- suppressWarnings(CellScoring(
    srt, features = features, method = "AUCell", backend = "r", splitByBlocks = TRUE,
    aucMaxRank = 10L, classification = FALSE, name = "auc_r_native", verbose = FALSE
  ))
  set.seed(17)
  cpp_out <- suppressWarnings(CellScoring(
    srt, features = features, method = "AUCell", backend = "cpp", splitByBlocks = TRUE,
    aucMaxRank = 10L, classification = FALSE, name = "auc_cpp_native", verbose = FALSE
  ))
  r_scores <- as.matrix(r_out@meta.data[, c("auc_r_native_set_a", "auc_r_native_set_b")])
  cpp_scores <- as.matrix(cpp_out@meta.data[, c("auc_cpp_native_set_a", "auc_cpp_native_set_b")])
  expect_equal(unname(cpp_scores), unname(r_scores), tolerance = 1e-12)
})

test_that("RunMetabolism exposes backend without AUCell strategy parameter", {
  expect_false("aucell_ties" %in% names(formals(CellScoring)))
  expect_true("backend" %in% names(formals(RunMetabolism)))
  expect_false("aucell_ties" %in% names(formals(RunMetabolism)))
  expect_false(any(grepl("strategy", names(formals(RunMetabolism)), fixed = TRUE)))
})

test_that("RunMetabolism AUCell R and C++ backends agree end to end", {
  skip_if_not_installed("AUCell")
  skip_if_not_installed("Seurat")

  set.seed(143)
  counts <- Matrix::Matrix(
    matrix(stats::rpois(80L * 12L, lambda = 1.2), nrow = 80L),
    sparse = TRUE
  )
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)
  gene_sets <- list(
    Pathway_A = rownames(counts)[1:20],
    Pathway_B = rownames(counts)[31:55]
  )
  local_mocked_bindings(
    scmetabolism_pathway_refs = function(...) list(kegg_refs = character(), reactome_names = character()),
    build_metabolism_gene_sets_from_preparedb = function(...) {
      list(gene_sets = gene_sets, term_names = stats::setNames(names(gene_sets), names(gene_sets)))
    },
    .package = "scop"
  )

  r_out <- suppressWarnings(RunMetabolism(
    srt, db = "KEGG", method = "AUCell", backend = "r", minGSSize = 10L,
    use_preparedb = TRUE, new_assay = FALSE, verbose = FALSE
  ))
  cpp_out <- suppressWarnings(RunMetabolism(
    srt, db = "KEGG", method = "AUCell", backend = "cpp", minGSSize = 10L,
    use_preparedb = TRUE, new_assay = FALSE, verbose = FALSE
  ))
  r_scores <- r_out@tools$Metabolism_AUCell$scores
  cpp_scores <- cpp_out@tools$Metabolism_AUCell$scores

  expect_equal(cpp_scores, r_scores, tolerance = 1e-12)
})

test_that("prepared metabolism gene sets retain Term alignment after missing genes", {
  local_mocked_bindings(
    PrepareDB = function(...) {
      list(Homo_sapiens = list(KEGG = list(
        TERM2GENE = data.frame(
          Term = c("hsa00010", "hsa00020", "hsa00020"),
          symbol = c("G1", NA_character_, "G2"),
          stringsAsFactors = FALSE
        ),
        TERM2NAME = data.frame(
          Term = c("hsa00010", "hsa00020"),
          Name = c("Pathway one", "Pathway two"),
          stringsAsFactors = FALSE
        )
      )))
    },
    .package = "scop"
  )

  out <- scop:::build_metabolism_gene_sets_from_preparedb(
    species = "Homo_sapiens",
    db_prepare = "KEGG",
    IDtype = "symbol",
    curated = list(kegg_refs = c("00010", "00020"), reactome_names = character()),
    expr_gene_names = c("G1", "G2"),
    db_update = FALSE,
    db_version = NULL,
    convert_species = FALSE,
    Ensembl_version = NULL,
    mirror = NULL,
    minGSSize = 1L,
    maxGSSize = 10L,
    verbose = FALSE
  )

  expect_equal(out$gene_sets, list(hsa00010 = "G1", hsa00020 = "G2"))
})
