make_abundance_contract_srt <- function(n = 80L) {
  counts <- Matrix::sparseMatrix(
    i = rep(1L, n), j = seq_len(n), x = 1,
    dims = c(2L, n),
    dimnames = list(c("gene1", "gene2"), paste0("cell", seq_len(n)))
  )
  srt <- SeuratObject::CreateSeuratObject(counts)
  srt$Condition <- rep(c("A", "B"), each = n / 2L)
  srt$CellType <- c(rep("T", 22L), rep("B", 18L),
    rep("T", 18L), rep("B", 22L))
  srt
}

test_that("permutation seed controls the public result", {
  srt <- make_abundance_contract_srt()
  run <- function() RunProportionTest(srt,
    group.by = "CellType", split.by = "Condition",
    comparison = list(c("A", "B")), proportion_method = "permutation",
    n_permutations = 40L, seed = 11L, verbose = FALSE
  )@tools$ProportionTest$results$A_vs_B
  set.seed(1)
  first <- run()
  set.seed(2)
  second <- run()
  expect_identical(first$pval, second$pval)
  expect_identical(first$boot_CI_2.5, second$boot_CI_2.5)
})

test_that("paired donor IDs retain both conditions in sample-level results", {
  make_cells <- function(donor, condition, t_count) {
    data.frame(CellType = c(rep("T", t_count), rep("B", 10L - t_count)),
      Condition = condition, Donor = donor)
  }
  dat <- do.call(rbind, list(
    make_cells("d1", "A", 8), make_cells("d1", "B", 3),
    make_cells("d2", "A", 6), make_cells("d2", "B", 2),
    make_cells("d3", "A", 7), make_cells("d3", "B", 4),
    make_cells("d4", "A", 9), make_cells("d4", "B", 5)
  ))
  helper <- getFromNamespace("sample_level_proportion_test", "scop")
  out <- helper(dat, "CellType", "Condition", "Donor", "A", "B",
    n_bootstrap = 40L, seed = 11L)
  t_row <- out[out$clusters == "T", , drop = FALSE]
  expect_true(is.finite(t_row$obs_log2FD))
  expect_true(is.finite(t_row$pval))
  expect_true(is.finite(t_row$boot_CI_2.5))
  expect_true(is.finite(t_row$boot_CI_97.5))
  expect_gt(t_row$obs_log2FD, 0)
  expect_equal(t_row$pval, stats::t.test(
    c(0.8, 0.6, 0.7, 0.9), c(0.3, 0.2, 0.4, 0.5), paired = TRUE
  )$p.value)
  partial <- dat[!(dat$Donor == "d4" & dat$Condition == "B"), , drop = FALSE]
  expect_error(helper(partial, "CellType", "Condition", "Donor",
    "A", "B", n_bootstrap = 0L), "fully paired")
})

test_that("Milo neighborhood uses SpatialFDR and keeps ordinary FDR", {
  raw <- data.frame(clusters = "nhood_1", neighborhood = "nhood_1",
    logFC = 2, PValue = 0.001, FDR = 0.8, SpatialFDR = 0.01)
  std <- getFromNamespace("standardize_proportion_result", "scop")(
    raw, "A", "B", "A_vs_B", "milo")
  expect_equal(std$FDR, 0.01)
  expect_equal(std$BH_FDR, 0.8)
  expect_identical(std$FDR_source, "SpatialFDR")
})

test_that("scCODA credibility drives effect direction without an FDR", {
  raw <- data.frame(clusters = "T", obs_log2FD = 2, pval = NA_real_,
    FDR = NA_real_, inclusion_prob = 0.99, credible = TRUE)
  std <- getFromNamespace("standardize_proportion_result", "scop")(
    raw, "A", "B", "A_vs_B", "sccoda")
  prepared <- getFromNamespace("prepare_proportion_plot_data", "scop")(
    std, FDR_threshold = 0.05, log2FD_threshold = log2(1.5),
    order_by = "value", nlabel = 1, features_label = NULL, label = FALSE)
  expect_identical(as.character(prepared$direction), "Increased")
  expect_identical(as.character(prepared$significance), "Credible effect")
  srt <- make_abundance_contract_srt()
  srt@tools$ProportionTest <- list(active_method = "sccoda",
    parameters = list(group.by = "CellType"),
    methods = list(sccoda = list(results = list(A_vs_B = raw),
      parameters = list(group.by = "CellType", credible_effect_threshold = 0.95))))
  plot <- ProportionTestPlot(srt, proportion_method = "sccoda", combine = FALSE)
  expect_identical(as.character(plot[[1]]$data$direction), "Increased")
  embedding <- cbind(UMAP_1 = seq_len(ncol(srt)), UMAP_2 = seq_len(ncol(srt)))
  rownames(embedding) <- colnames(srt)
  srt[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding, key = "UMAP_", assay = "RNA")
  umap <- ProportionTestPlot(srt, proportion_method = "sccoda",
    plot_type = "umap", combine = FALSE)
  direction <- stats::setNames(
    as.character(umap[[1]]$data$.proportion_da_direction_A_vs_B),
    rownames(umap[[1]]$data)
  )
  expect_identical(unname(direction[c("cell1", "cell30")]),
    c("Increased", "Uncovered"))
})

test_that("Milo neighborhood UMAP uses stored memberships and marks overlaps", {
  srt <- make_abundance_contract_srt()
  embedding <- cbind(UMAP_1 = seq_len(ncol(srt)), UMAP_2 = seq_len(ncol(srt)))
  rownames(embedding) <- colnames(srt)
  srt[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding, key = "UMAP_", assay = "RNA")
  nhood <- data.frame(neighborhood = c("nhood_1", "nhood_2"),
    clusters = c("nhood_1", "nhood_2"), obs_log2FD = c(2, -2),
    pval = c(0.001, 0.001), FDR = c(0.01, 0.01))
  members <- list(nhood_1 = c("cell1", "cell2"),
    nhood_2 = c("cell2", "cell3"))
  srt@tools$ProportionTest <- list(active_method = "milo",
    parameters = list(group.by = "CellType"),
    methods = list(milo = list(
      results = list(A_vs_B = nhood),
      neighborhood_results = list(A_vs_B = nhood),
      parameters = list(group.by = "CellType"),
      details = list(milo_graph_data = list(.metadata = list(members = members)))
    )))
  effect_plot <- ProportionTestPlot(srt, proportion_method = "milo",
    result_level = "neighborhood", plot_type = "effect", combine = FALSE)
  expect_identical(effect_plot[[1]]$labels$x, "Milo neighborhood")
  expect_identical(effect_plot[[1]]$labels$y, "Milo log2 FC")
  plot <- ProportionTestPlot(srt, proportion_method = "milo",
    result_level = "neighborhood", plot_type = "umap", combine = FALSE)
  direction <- stats::setNames(
    as.character(plot[[1]]$data$.proportion_da_direction_A_vs_B),
    rownames(plot[[1]]$data)
  )
  expect_identical(unname(direction[paste0("cell", seq_len(4L))]),
    c("Increased", "Mixed", "Decreased", "Uncovered"))
  continuous <- ProportionTestPlot(srt, proportion_method = "milo",
    result_level = "neighborhood", plot_type = "umap",
    umap_mode = "continuous", combine = FALSE)
  effect <- stats::setNames(
    continuous[[1]]$data$.proportion_da_log2fd_A_vs_B,
    rownames(continuous[[1]]$data)
  )
  expect_equal(unname(effect[paste0("cell", seq_len(3L))]), c(2, 0, -2))
  srt@tools$ProportionTest$methods$milo$details$milo_graph_data$.metadata$members <- NULL
  expect_error(ProportionTestPlot(srt, proportion_method = "milo",
    result_level = "neighborhood", plot_type = "umap", combine = FALSE),
    "membership")
})

test_that("virtual samples require opt-in and suppress inferential columns", {
  srt <- make_abundance_contract_srt()
  expect_error(RunProportionTest(srt,
    group.by = "CellType", split.by = "Condition",
    proportion_method = "propeller", n_bootstrap = 0L, verbose = FALSE),
    "sample.by")
  out <- suppressWarnings(RunProportionTest(srt,
    group.by = "CellType", split.by = "Condition",
    proportion_method = "propeller", allow_pseudo_samples = TRUE,
    n_bootstrap = 0L, verbose = FALSE))
  rows <- out@tools$ProportionTest$results$A_vs_B
  expect_true(all(is.na(rows$pval)))
  expect_true(all(is.na(rows$FDR)))
  expect_false(out@tools$ProportionTest$methods$propeller$inference_valid)
  plot <- ProportionTestPlot(out, proportion_method = "propeller", combine = FALSE)
  expect_true(all(as.character(plot[[1]]$data$direction) == "Not tested"))
  expect_error(ProportionTestPlot(out, proportion_method = "propeller",
    result_level = "neighborhood", combine = FALSE), "not available")
})

test_that("Milo C++ keeps paired sample-condition units", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")
  skip_if_not_installed("BiocNeighbors")
  set.seed(42)
  n <- 160L
  counts <- Matrix::sparseMatrix(i = rep(1L, n), j = seq_len(n), x = 1,
    dims = c(2L, n),
    dimnames = list(c("gene1", "gene2"), paste0("cell", seq_len(n))))
  srt <- SeuratObject::CreateSeuratObject(counts)
  srt$Sample <- rep(rep(paste0("donor", seq_len(4L)), each = 20L), times = 2L)
  srt$Condition <- rep(c("A", "B"), each = 80L)
  srt$CellType <- sample(c("T", "B"), n, replace = TRUE)
  embedding <- cbind(PC_1 = rnorm(n), PC_2 = rnorm(n))
  rownames(embedding) <- colnames(srt)
  srt[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding, key = "PC_", assay = "RNA")
  result <- RunProportionTest(srt, group.by = "CellType",
    split.by = "Condition", sample.by = "Sample",
    comparison = list(c("A", "B")), proportion_method = "milo",
    backend = "cpp", reduction = "pca", milo_k = 10L, milo_d = 2L,
    n_bootstrap = 0L, verbose = FALSE)
  bundle <- result@tools$ProportionTest$methods$milo
  expect_true(nrow(bundle$neighborhood_results$A_vs_B) > 0L)
  expect_equal(bundle$neighborhood_results$A_vs_B$FDR,
    bundle$neighborhood_results$A_vs_B$SpatialFDR)
  expect_true(all(is.finite(bundle$results$A_vs_B$obs_log2FD)))
})
