make_proportion_srt <- function() {
  n_a_t <- 40L
  n_a_b <- 10L
  n_b_t <- 10L
  n_b_b <- 40L
  n_cells <- n_a_t + n_a_b + n_b_t + n_b_b
  counts <- Matrix::sparseMatrix(
    i = rep(1L, n_cells),
    j = seq_len(n_cells),
    x = 1,
    dims = c(2L, n_cells),
    dimnames = list(c("gene1", "gene2"), paste0("cell", seq_len(n_cells)))
  )
  srt <- SeuratObject::CreateSeuratObject(counts = counts)
  srt$CellType <- c(
    rep("T", n_a_t),
    rep("B", n_a_b),
    rep("T", n_b_t),
    rep("B", n_b_b)
  )
  srt$Condition <- c(
    rep("A", n_a_t + n_a_b),
    rep("B", n_b_t + n_b_b)
  )
  srt$Sample <- c(
    rep("s1", 25L),
    rep("s2", 25L),
    rep("s3", 25L),
    rep("s4", 25L)
  )
  srt
}

test_that("permutation obs_log2FD is log2(group1 / group2)", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  srt <- make_proportion_srt()
  out <- RunProportionTest(
    srt,
    group.by = "CellType",
    split.by = "Condition",
    comparison = list(c("A", "B")),
    proportion_method = "permutation",
    n_permutations = 20,
    include_all_cells = FALSE,
    verbose = FALSE
  )
  res <- out@tools$ProportionTest$results[["A_vs_B"]]
  t_row <- res[res$clusters == "T", , drop = FALSE]
  b_row <- res[res$clusters == "B", , drop = FALSE]

  expect_equal(as.numeric(t_row$obs_log2FD), 2, tolerance = 1e-6)
  expect_equal(as.numeric(b_row$obs_log2FD), -2, tolerance = 1e-6)
  expect_gt(as.numeric(t_row$obs_log2FD), 0)
  expect_lt(as.numeric(b_row$obs_log2FD), 0)
})

test_that("sample-level obs_log2FD is log2(group1 / group2)", {
  sample_level <- getFromNamespace("sample_level_proportion_test", "scop")
  meta_data <- data.frame(
    CellType = c(
      "T", "T", "T", "B",
      "T", "T", "T", "B",
      "T", "B", "B", "B",
      "T", "B", "B", "B"
    ),
    Condition = c(rep("A", 8), rep("B", 8)),
    Sample = rep(c("s1", "s2", "s3", "s4"), each = 4),
    stringsAsFactors = FALSE
  )

  res <- sample_level(
    meta_data = meta_data,
    group.by = "CellType",
    split.by = "Condition",
    sample.by = "Sample",
    cluster_1 = "A",
    cluster_2 = "B",
    n_bootstrap = 20,
    seed = 1,
    verbose = FALSE
  )
  t_row <- res[res$clusters == "T", , drop = FALSE]
  b_row <- res[res$clusters == "B", , drop = FALSE]

  expect_gt(as.numeric(t_row$obs_log2FD), 0)
  expect_lt(as.numeric(b_row$obs_log2FD), 0)
  expect_equal(
    as.numeric(t_row$obs_log2FD),
    log2((0.75 + 1e-5) / (0.25 + 1e-5)),
    tolerance = 1e-8
  )
})

test_that("bootstrap log2FD uses group1 as numerator", {
  boot <- getFromNamespace("proportion_bootstrap_log2fd", "scop")
  set.seed(11)
  values <- boot(
    v1 = c(0.8, 0.8, 0.8),
    v2 = c(0.2, 0.2, 0.2),
    n_bootstrap = 30,
    pseudocount = 1e-8,
    verbose = FALSE
  )
  expect_true(all(values > 0))
  expect_equal(mean(values), log2(0.8 / 0.2), tolerance = 1e-6)
})

test_that("proportion plot title follows group1 vs group2", {
  title_fun <- getFromNamespace("proportion_title", "scop")
  df <- data.frame(group1 = "A", group2 = "B", stringsAsFactors = FALSE)
  expect_identical(title_fun(df), "A vs B")
})

test_that("ProportionTestPlot accepts Milo neighborhood result level", {
  result_levels <- eval(formals(ProportionTestPlot)$result_level)

  expect_identical(result_levels, c("group", "neighborhood"))
  expect_identical(match.arg("neighborhood", result_levels), "neighborhood")
})

make_neighborhood_plot_srt <- function(members = TRUE) {
  counts <- matrix(1, nrow = 2, ncol = 4, dimnames = list(c("g1", "g2"), paste0("cell", 1:4)))
  srt <- Seurat::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE))
  srt$CellType <- c("T", "T", "B", "NK")
  nhood <- data.frame(
    clusters = c("nhood_1", "nhood_2", "nhood_3"),
    neighborhood = c("nhood_1", "nhood_2", "nhood_3"),
    obs_log2FD = c(2, 0.1, -2),
    FDR = c(0.001, 0.8, 0.001),
    pval = c(0.001, 0.8, 0.001),
    stringsAsFactors = FALSE
  )
  group <- data.frame(
    clusters = c("T", "B", "NK"),
    obs_log2FD = c(2, 0.1, -2),
    FDR = c(0.001, 0.8, 0.001),
    pval = c(0.001, 0.8, 0.001),
    stringsAsFactors = FALSE
  )
  stored_members <- if (isTRUE(members)) {
    list(
      nhood_1 = c("cell1", "cell2"),
      nhood_2 = "cell3",
      nhood_3 = "cell4"
    )
  } else {
    NULL
  }
  srt@tools[["ProportionTest"]] <- list(
    active_method = "milo",
    parameters = list(group.by = "CellType"),
    methods = list(
      milo = list(
        results = list(Treat_vs_Ctrl = group),
        neighborhood_results = list(Treat_vs_Ctrl = nhood),
        parameters = list(group.by = "CellType"),
        details = list(
          milo_graph_data = list(.metadata = list(members = stored_members))
        )
      )
    )
  )
  srt
}

test_that("neighborhood DA maps onto member cells instead of cell-group labels", {
  project <- getFromNamespace("project_neighborhood_da_to_cells", "scop")
  srt <- make_neighborhood_plot_srt()
  nhood <- srt@tools$ProportionTest$methods$milo$neighborhood_results$Treat_vs_Ctrl
  members <- srt@tools$ProportionTest$methods$milo$details$milo_graph_data$.metadata$members
  members$nhood_1 <- c("cell1", "cell2", "cell4")

  projected <- project(
    nhood,
    members,
    colnames(srt),
    FDR_threshold = 0.05,
    log2FD_threshold = log2(1.5)
  )

  expect_identical(as.character(projected$direction), c("Increased", "Increased", "NS", "Increased"))
  expect_equal(projected$obs_log2FD, c(2, 2, 0.1, 2))
})

test_that("missing Milo neighborhood results are not replaced with group results", {
  get_results <- getFromNamespace("get_proportion_plot_results", "scop")
  srt <- make_neighborhood_plot_srt()
  srt@tools$ProportionTest$methods$milo$neighborhood_results <- NULL

  expect_error(
    get_results(srt, result_level = "neighborhood"),
    "Neighborhood-level results are not available"
  )
})

test_that("neighborhood UMAP requires stored neighborhood membership", {
  srt <- make_neighborhood_plot_srt(members = FALSE)

  expect_error(
    ProportionTestPlot(
      srt,
      result_level = "neighborhood",
      plot_type = "umap",
      verbose = FALSE
    ),
    "RunMilo"
  )
})
