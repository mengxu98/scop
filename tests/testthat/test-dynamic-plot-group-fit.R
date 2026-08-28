make_dynamic_plot_test_object <- function(
  n_per_group = 12,
  factor_condition = FALSE,
  gene1_values = NULL
) {
  n_cells <- 2 * n_per_group
  cells <- paste0("cell", seq_len(n_cells))
  gene1_values <- gene1_values %||% seq_len(n_cells)
  counts <- Matrix::Matrix(
    rbind(
      Gene1 = gene1_values,
      Gene2 = rev(seq_len(n_cells))
    ),
    sparse = TRUE
  )
  colnames(counts) <- cells
  srt <- Seurat::CreateSeuratObject(counts)
  srt$Lineage1 <- rep(seq(0, 1, length.out = n_per_group), 2)
  srt$Lineage2 <- rep(seq(1, 0, length.out = n_per_group), 2)
  condition <- rep(c("control", "HUA"), each = n_per_group)
  srt$condition <- if (isTRUE(factor_condition)) {
    factor(condition, levels = c("HUA", "control"))
  } else {
    condition
  }
  srt
}

dynamic_plot_built_layer <- function(plot, type = c("point", "line")) {
  type <- match.arg(type)
  layers <- ggplot2::ggplot_build(plot)$data
  matches <- vapply(layers, function(x) {
    if (type == "point") {
      "shape" %in% names(x)
    } else {
      "linewidth" %in% names(x) && !"ymin" %in% names(x)
    }
  }, logical(1))
  layers[[which(matches)[1]]]
}

mock_dynamic_fit_result <- function(
  srt,
  lineages,
  features,
  fitted_by_group = c(control = 1, HUA = 2)
) {
  cells <- rownames(srt@meta.data)[is.finite(srt@meta.data[[lineages]])]
  pseudotime <- srt@meta.data[cells, lineages, drop = TRUE]
  raw_values <- as.matrix(srt[["RNA"]]$counts[features, cells, drop = FALSE])
  raw <- cbind(pseudotime = pseudotime, t(raw_values))
  group <- as.character(srt$condition[cells][1])
  fit_base <- fitted_by_group[[group]]
  fitted_values <- vapply(
    seq_along(features),
    function(i) rep(fit_base + i - 1, length(cells)),
    numeric(length(cells))
  )
  fitted <- cbind(pseudotime = pseudotime, fitted_values)
  colnames(fitted) <- c("pseudotime", features)
  rownames(raw) <- rownames(fitted) <- cells
  upr <- lwr <- fitted
  upr[, features] <- upr[, features, drop = FALSE] + 0.1
  lwr[, features] <- lwr[, features, drop = FALSE] - 0.1
  srt@tools[[paste0("DynamicFeatures_", lineages)]] <- list(
    raw_matrix = raw,
    fitted_matrix = fitted,
    upr_matrix = upr,
    lwr_matrix = lwr
  )
  srt
}

test_that("DynamicPlot fits each fit.by group independently", {
  srt <- make_dynamic_plot_test_object()
  fitted_cells <- list()

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      cells <- rownames(srt@meta.data)[is.finite(srt@meta.data[[lineages]])]
      fitted_cells[[length(fitted_cells) + 1L]] <<- cells
      mock_dynamic_fit_result(srt, lineages, features)
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
  expect_setequal(fitted_cells[[1]], paste0("cell", 1:12))
  expect_setequal(fitted_cells[[2]], paste0("cell", 13:24))

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
      mock_dynamic_fit_result(srt, lineages, features)
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
  expect_setequal(fitted_cells[[1]], paste0("cell", 13:24))
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

test_that("DynamicPlot keeps raw points when a group lacks enough GAM support", {
  srt <- make_dynamic_plot_test_object()
  srt$Lineage1[paste0("cell", 10:12)] <- NA_real_
  fitted_groups <- character()

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      cells <- rownames(srt@meta.data)[is.finite(srt@meta.data[[lineages]])]
      fitted_groups <<- c(fitted_groups, as.character(srt$condition[cells][1]))
      mock_dynamic_fit_result(srt, lineages, features)
    },
    .package = "scop"
  )

  plots <- DynamicPlot(
    srt,
    lineages = "Lineage1",
    features = "Gene1",
    group.by = "condition",
    fit.by = "condition",
    exp_method = "raw",
    lib_normalize = FALSE,
    add_interval = FALSE,
    add_rug = FALSE,
    combine = FALSE,
    verbose = FALSE
  )

  expect_identical(fitted_groups, "HUA")
  point_data <- dynamic_plot_built_layer(plots[[1]], "point")
  expect_equal(nrow(point_data), 21)
  expect_equal(length(unique(point_data$group)), 2)
})

test_that("DynamicPlot subsets custom library sizes for grouped raw points", {
  srt <- make_dynamic_plot_test_object()
  custom_libsize <- seq_len(ncol(srt))

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(srt, lineages, features)
    },
    .package = "scop"
  )

  plots <- DynamicPlot(
    srt,
    lineages = "Lineage1",
    features = "Gene1",
    group.by = "condition",
    fit.by = "condition",
    libsize = custom_libsize,
    exp_method = "raw",
    add_line = FALSE,
    add_interval = FALSE,
    add_rug = FALSE,
    combine = FALSE,
    verbose = FALSE
  )

  point_data <- dynamic_plot_built_layer(plots[[1]], "point")
  expect_equal(unique(point_data$y), stats::median(custom_libsize))
})

test_that("DynamicPlot uses shared expression transforms across fit groups", {
  values <- c(rep(1, 12), rep(3, 12))
  srt <- make_dynamic_plot_test_object(gene1_values = values)
  expected <- list(
    fc = c(0.5, 1.5),
    log2fc = log2(c(0.5, 1.5)),
    zscore = sort(unique(as.numeric(scale(values))))
  )

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(
        srt,
        lineages,
        features,
        fitted_by_group = c(control = 1, HUA = 3)
      )
    },
    .package = "scop"
  )

  for (method in names(expected)) {
    plots <- DynamicPlot(
      srt,
      lineages = "Lineage1",
      features = "Gene1",
      fit.by = "condition",
      exp_method = method,
      lib_normalize = FALSE,
      add_point = FALSE,
      add_interval = FALSE,
      add_rug = FALSE,
      combine = FALSE,
      verbose = FALSE
    )
    line_data <- dynamic_plot_built_layer(plots[[1]], "line")
    expect_equal(
      sort(unique(line_data$y)),
      expected[[method]],
      tolerance = 1e-8,
      info = method
    )
  }
})

test_that("DynamicPlot retains series identity in grouped combined panels", {
  srt <- make_dynamic_plot_test_object()

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(srt, lineages, features)
    },
    .package = "scop"
  )

  plots <- DynamicPlot(
    srt,
    lineages = c("Lineage1", "Lineage2"),
    features = c("Gene1", "Gene2"),
    fit.by = "condition",
    compare_lineages = TRUE,
    compare_features = TRUE,
    exp_method = "raw",
    lib_normalize = FALSE,
    add_point = FALSE,
    add_interval = FALSE,
    add_rug = FALSE,
    combine = FALSE,
    verbose = FALSE
  )

  line_data <- dynamic_plot_built_layer(plots[[1]], "line")
  expect_equal(length(unique(line_data$linetype)), 4)
})

test_that("DynamicPlot keeps a grouped interval legend without fitted lines", {
  srt <- make_dynamic_plot_test_object()
  captured_plot <- NULL

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(srt, lineages, features)
    },
    get_legend = function(plot) {
      captured_plot <<- plot
      grid::nullGrob()
    },
    .package = "scop"
  )

  DynamicPlot(
    srt,
    lineages = "Lineage1",
    features = "Gene1",
    fit.by = "condition",
    exp_method = "raw",
    lib_normalize = FALSE,
    add_point = FALSE,
    add_line = FALSE,
    add_interval = TRUE,
    add_rug = FALSE,
    combine = FALSE,
    verbose = FALSE
  )

  fill_scales <- Filter(
    function(x) any(grepl("fill", x$aesthetics, fixed = TRUE)),
    captured_plot$scales$scales
  )
  expect_true(any(vapply(
    fill_scales,
    function(x) !identical(x$guide, "none"),
    logical(1)
  )))
})

test_that("DynamicPlot uses factor levels for matching point and line colors", {
  values <- c(rep(1, 12), rep(2, 12))
  srt <- make_dynamic_plot_test_object(
    factor_condition = TRUE,
    gene1_values = values
  )

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(srt, lineages, features)
    },
    .package = "scop"
  )

  plots <- DynamicPlot(
    srt,
    lineages = "Lineage1",
    features = "Gene1",
    group.by = "condition",
    fit.by = "condition",
    point_palcolor = c("#D73027", "#2C7FB8"),
    exp_method = "raw",
    lib_normalize = FALSE,
    add_interval = FALSE,
    add_rug = FALSE,
    combine = FALSE,
    verbose = FALSE
  )

  point_data <- dynamic_plot_built_layer(plots[[1]], "point")
  line_data <- dynamic_plot_built_layer(plots[[1]], "line")
  point_colors <- vapply(
    split(point_data$colour, point_data$y),
    function(x) unique(x)[1],
    character(1)
  )
  line_colors <- vapply(
    split(line_data$colour, line_data$y),
    function(x) unique(x)[1],
    character(1)
  )
  expect_identical(line_colors[names(point_colors)], point_colors)
})

test_that("DynamicPlot preserves metadata named FitGroup", {
  srt <- make_dynamic_plot_test_object()
  srt$FitGroup <- rep(c("meta-A", "meta-B"), times = 12)

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(srt, lineages, features)
    },
    .package = "scop"
  )

  plots <- DynamicPlot(
    srt,
    lineages = "Lineage1",
    features = "Gene1",
    group.by = "FitGroup",
    exp_method = "raw",
    lib_normalize = FALSE,
    add_line = FALSE,
    add_interval = FALSE,
    add_rug = FALSE,
    combine = FALSE,
    verbose = FALSE
  )

  point_data <- dynamic_plot_built_layer(plots[[1]], "point")
  expect_equal(length(unique(point_data$colour)), 2)
})

test_that("DynamicPlot applies custom transforms over the shared lineage", {
  srt <- make_dynamic_plot_test_object(
    gene1_values = c(rep(1, 12), rep(3, 12))
  )

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(
        srt,
        lineages,
        features,
        fitted_by_group = c(control = 1, HUA = 3)
      )
    },
    .package = "scop"
  )

  plots <- DynamicPlot(
    srt,
    lineages = "Lineage1",
    features = "Gene1",
    group.by = "condition",
    fit.by = "condition",
    exp_method = function(x) x - mean(x, na.rm = TRUE),
    lib_normalize = FALSE,
    add_line = FALSE,
    add_interval = FALSE,
    add_rug = FALSE,
    combine = FALSE,
    verbose = FALSE
  )

  point_data <- dynamic_plot_built_layer(plots[[1]], "point")
  group_means <- sort(as.numeric(tapply(
    point_data$y,
    point_data$colour,
    mean
  )))
  expect_equal(group_means, c(-1, 1))
})

test_that("DynamicPlot caps log2fc infinities over the shared lineage", {
  values <- c(0, 1:11, 0, rep(100, 11))
  srt <- make_dynamic_plot_test_object(gene1_values = values)

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(srt, lineages, features)
    },
    .package = "scop"
  )

  plots <- DynamicPlot(
    srt,
    lineages = "Lineage1",
    features = "Gene1",
    group.by = "condition",
    fit.by = "condition",
    exp_method = "log2fc",
    lib_normalize = FALSE,
    add_line = FALSE,
    add_interval = FALSE,
    add_rug = FALSE,
    combine = FALSE,
    verbose = FALSE
  )

  point_data <- dynamic_plot_built_layer(plots[[1]], "point")
  group_minima <- tapply(point_data$y, point_data$colour, min)
  expect_equal(length(unique(round(group_minima, 8))), 1)
})

test_that("DynamicPlot distinguishes interval-only combined series", {
  srt <- make_dynamic_plot_test_object()

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(srt, lineages, features)
    },
    .package = "scop"
  )

  plots <- DynamicPlot(
    srt,
    lineages = c("Lineage1", "Lineage2"),
    features = c("Gene1", "Gene2"),
    fit.by = "condition",
    compare_lineages = TRUE,
    compare_features = TRUE,
    exp_method = "raw",
    lib_normalize = FALSE,
    add_point = FALSE,
    add_line = FALSE,
    add_interval = TRUE,
    add_rug = FALSE,
    combine = FALSE,
    verbose = FALSE
  )

  built <- ggplot2::ggplot_build(plots[[1]])$data
  ribbon_data <- built[[which(vapply(
    built,
    function(x) "ymin" %in% names(x),
    logical(1)
  ))[1]]]
  expect_equal(length(unique(ribbon_data$linetype)), 4)
})

test_that("DynamicPlot rejects more than six grouped combined series", {
  srt <- make_dynamic_plot_test_object()
  for (i in 3:7) {
    srt[[paste0("Lineage", i)]] <- srt$Lineage1
  }

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(srt, lineages, features)
    },
    .package = "scop"
  )

  expect_error(
    DynamicPlot(
      srt,
      lineages = paste0("Lineage", 1:7),
      features = "Gene1",
      fit.by = "condition",
      compare_lineages = TRUE,
      exp_method = "raw",
      lib_normalize = FALSE,
      add_point = FALSE,
      add_interval = FALSE,
      add_rug = FALSE,
      combine = FALSE,
      verbose = FALSE
    ),
    "at most six"
  )
})

test_that("DynamicPlot builds grouped legends from every panel", {
  srt <- make_dynamic_plot_test_object()
  srt$Lineage1[paste0("cell", 13:24)] <- NA_real_
  srt$Lineage2[paste0("cell", 1:12)] <- NA_real_

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(srt, lineages, features)
    },
    .package = "scop"
  )

  plots <- DynamicPlot(
    srt,
    lineages = c("Lineage1", "Lineage2"),
    features = "Gene1",
    fit.by = "condition",
    compare_lineages = FALSE,
    exp_method = "raw",
    lib_normalize = FALSE,
    add_point = FALSE,
    add_interval = FALSE,
    add_rug = FALSE,
    combine = FALSE,
    verbose = FALSE
  )

  built <- ggplot2::ggplot_build(plots[[1]])
  fit_scales <- Filter(
    function(x) identical(x$name, "condition"),
    built$plot$scales$scales
  )
  expect_length(fit_scales, 1)
  expect_setequal(fit_scales[[1]]$get_breaks(), c("control", "HUA"))
})

test_that("DynamicPlot deduplicates points across compared lineages", {
  srt <- make_dynamic_plot_test_object()

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(srt, lineages, features)
    },
    .package = "scop"
  )

  plots <- DynamicPlot(
    srt,
    lineages = c("Lineage1", "Lineage2"),
    features = "Gene1",
    group.by = "condition",
    fit.by = "condition",
    compare_lineages = TRUE,
    exp_method = "raw",
    lib_normalize = FALSE,
    add_line = FALSE,
    add_interval = FALSE,
    add_rug = FALSE,
    combine = FALSE,
    verbose = FALSE
  )

  point_data <- dynamic_plot_built_layer(plots[[1]], "point")
  expect_equal(nrow(point_data), ncol(srt))
})
