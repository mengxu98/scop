make_dynamic_plot_test_object <- function() {
  cells <- paste0("cell", seq_len(12))
  counts <- Matrix::Matrix(
    matrix(
      seq_len(12),
      nrow = 1,
      dimnames = list("Gene1", cells)
    ),
    sparse = TRUE
  )
  srt <- Seurat::CreateSeuratObject(counts)
  srt$Lineage1 <- rep(seq(0, 1, length.out = 6), 2)
  srt$condition <- rep(c("control", "HUA"), each = 6)
  srt
}

test_that("DynamicPlot fits each fit.by group independently", {
  srt <- make_dynamic_plot_test_object()
  fitted_cells <- list()

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      cells <- rownames(srt@meta.data)[is.finite(srt@meta.data[[lineages]])]
      fitted_cells[[length(fitted_cells) + 1L]] <<- cells
      fit_value <- if (all(srt$condition[cells] == "control")) 1 else 2
      pseudotime <- srt@meta.data[cells, lineages, drop = TRUE]
      raw <- cbind(
        pseudotime = pseudotime,
        Gene1 = as.numeric(srt[["RNA"]]$counts["Gene1", cells])
      )
      rownames(raw) <- cells
      fitted <- cbind(
        pseudotime = pseudotime,
        Gene1 = rep(fit_value, length(cells))
      )
      rownames(fitted) <- cells
      upr <- lwr <- fitted
      upr[, "Gene1"] <- upr[, "Gene1"] + 0.1
      lwr[, "Gene1"] <- lwr[, "Gene1"] - 0.1
      srt@tools[[paste0("DynamicFeatures_", lineages)]] <- list(
        raw_matrix = raw,
        fitted_matrix = fitted,
        upr_matrix = upr,
        lwr_matrix = lwr
      )
      srt
    },
    .package = "scop"
  )

  plots <- DynamicPlot(
    srt,
    lineages = "Lineage1",
    features = "Gene1",
    group.by = "condition",
    fit.by = "condition",
    layer = "counts",
    exp_method = "raw",
    lib_normalize = FALSE,
    add_point = FALSE,
    add_rug = FALSE,
    add_interval = FALSE,
    combine = FALSE,
    verbose = FALSE
  )

  expect_length(fitted_cells, 2)
  expect_setequal(fitted_cells[[1]], paste0("cell", 1:6))
  expect_setequal(fitted_cells[[2]], paste0("cell", 7:12))

  built <- ggplot2::ggplot_build(plots[[1]])
  line_data <- built$data[[1]]
  expect_equal(length(unique(line_data$group)), 2)
  expect_equal(sort(unique(line_data$y)), c(1, 2))
})

test_that("DynamicPlot group_use limits independently fitted groups", {
  srt <- make_dynamic_plot_test_object()
  fitted_cells <- list()

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      cells <- rownames(srt@meta.data)[is.finite(srt@meta.data[[lineages]])]
      fitted_cells[[length(fitted_cells) + 1L]] <<- cells
      pseudotime <- srt@meta.data[cells, lineages, drop = TRUE]
      matrix_out <- cbind(
        pseudotime = pseudotime,
        Gene1 = rep(1, length(cells))
      )
      rownames(matrix_out) <- cells
      srt@tools[[paste0("DynamicFeatures_", lineages)]] <- list(
        raw_matrix = matrix_out,
        fitted_matrix = matrix_out,
        upr_matrix = matrix_out,
        lwr_matrix = matrix_out
      )
      srt
    },
    .package = "scop"
  )

  expect_no_error(DynamicPlot(
    srt,
    lineages = "Lineage1",
    features = "Gene1",
    group.by = "condition",
    group_use = "HUA",
    fit.by = "condition",
    exp_method = "raw",
    lib_normalize = FALSE,
    add_point = FALSE,
    add_rug = FALSE,
    add_interval = FALSE,
    combine = FALSE,
    verbose = FALSE
  ))

  expect_length(fitted_cells, 1)
  expect_setequal(fitted_cells[[1]], paste0("cell", 7:12))
})

test_that("DynamicPlot preserves the existing overall fit by default", {
  srt <- make_dynamic_plot_test_object()
  cells <- rownames(srt@meta.data)
  pseudotime <- srt$Lineage1
  raw <- cbind(
    pseudotime = pseudotime,
    Gene1 = as.numeric(srt[["RNA"]]$counts["Gene1", cells])
  )
  fitted <- cbind(pseudotime = pseudotime, Gene1 = rep(3, length(cells)))
  rownames(raw) <- rownames(fitted) <- cells
  srt@tools$DynamicFeatures_Lineage1 <- list(
    raw_matrix = raw,
    fitted_matrix = fitted,
    upr_matrix = fitted,
    lwr_matrix = fitted
  )

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(...) {
      stop("The cached overall fit should be reused")
    },
    .package = "scop"
  )

  plots <- DynamicPlot(
    srt,
    lineages = "Lineage1",
    features = "Gene1",
    group.by = "condition",
    exp_method = "raw",
    lib_normalize = FALSE,
    add_point = FALSE,
    add_rug = FALSE,
    add_interval = FALSE,
    combine = FALSE,
    verbose = FALSE
  )

  line_data <- ggplot2::ggplot_build(plots[[1]])$data[[1]]
  expect_equal(length(unique(line_data$group)), 1)
  expect_equal(unique(line_data$y), 3)
})

test_that("DynamicPlot validates fit.by", {
  srt <- make_dynamic_plot_test_object()

  expect_error(
    DynamicPlot(
      srt,
      lineages = "Lineage1",
      features = "Gene1",
      fit.by = "missing_column",
      verbose = FALSE
    ),
    "not in the meta.data"
  )
})
