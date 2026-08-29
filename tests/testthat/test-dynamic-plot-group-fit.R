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
  gene_features <- features[features %in% rownames(srt[["RNA"]])]
  meta_features <- features[features %in% colnames(srt@meta.data)]
  raw_values <- matrix(
    NA_real_,
    nrow = length(features),
    ncol = length(cells),
    dimnames = list(features, cells)
  )
  if (length(gene_features) > 0) {
    raw_values[gene_features, ] <- as.matrix(
      srt[["RNA"]]$counts[gene_features, cells, drop = FALSE]
    )
  }
  if (length(meta_features) > 0) {
    raw_values[meta_features, ] <- t(as.matrix(
      srt@meta.data[cells, meta_features, drop = FALSE]
    ))
  }
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
  custom_libsize[24] <- NA_real_

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(srt, lineages, features)
    },
    .package = "scop"
  )

  plots <- suppressWarnings(DynamicPlot(
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
  ))

  point_data <- dynamic_plot_built_layer(plots[[1]], "point")
  expect_equal(sum(is.finite(point_data$y)), 23)
  expect_equal(
    unique(point_data$y[is.finite(point_data$y)]),
    stats::median(custom_libsize, na.rm = TRUE)
  )
})

test_that("DynamicPlot masks zero library sizes before shared transforms", {
  srt <- make_dynamic_plot_test_object()
  custom_libsize <- rep(1, ncol(srt))
  custom_libsize[24] <- 0

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(srt, lineages, features)
    },
    .package = "scop"
  )

  plots <- suppressWarnings(DynamicPlot(
    srt,
    lineages = "Lineage1",
    features = "Gene1",
    group.by = "condition",
    fit.by = "condition",
    libsize = custom_libsize,
    exp_method = "zscore",
    add_line = FALSE,
    add_interval = FALSE,
    add_rug = FALSE,
    combine = FALSE,
    verbose = FALSE
  ))

  point_data <- dynamic_plot_built_layer(plots[[1]], "point")
  expect_equal(sum(is.finite(point_data$y)), 23)
  expect_false(any(is.infinite(point_data$y)))
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

test_that("DynamicPlot shared transforms ignore partial missing values", {
  srt <- make_dynamic_plot_test_object()
  srt$score <- seq_len(ncol(srt))
  srt$score[3] <- NA_real_

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(srt, lineages, features)
    },
    .package = "scop"
  )

  for (method in c("zscore", "fc", "log2fc")) {
    plots <- suppressWarnings(DynamicPlot(
      srt,
      lineages = "Lineage1",
      features = "score",
      fit.by = "condition",
      exp_method = method,
      lib_normalize = FALSE,
      add_line = FALSE,
      add_interval = FALSE,
      add_rug = FALSE,
      combine = FALSE,
      verbose = FALSE
    ))
    point_data <- dynamic_plot_built_layer(plots[[1]], "point")
    expect_equal(sum(is.finite(point_data$y)), 23)
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

test_that("DynamicPlot keeps a point guide for globally raw-only fit groups", {
  cells <- paste0("cell", seq_len(25))
  counts <- Matrix::Matrix(
    rbind(Gene1 = seq_len(25)),
    sparse = TRUE
  )
  colnames(counts) <- cells
  srt <- Seurat::CreateSeuratObject(counts)
  srt$condition <- c(rep("A", 10), rep("B", 10), rep("C", 5))
  srt$Lineage1 <- c(
    seq(0, 1, length.out = 10),
    seq(0, 1, length.out = 10),
    rep(NA_real_, 5)
  )
  srt$Lineage2 <- c(
    seq(0, 1, length.out = 10),
    rep(NA_real_, 10),
    seq(0, 1, length.out = 5)
  )
  legend_plot <- NULL

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(
        srt,
        lineages,
        features,
        fitted_by_group = c(A = 1, B = 2, C = 3)
      )
    },
    get_legend = function(plot) {
      legend_plot <<- plot
      grid::nullGrob()
    },
    .package = "scop"
  )

  DynamicPlot(
    srt,
    lineages = c("Lineage1", "Lineage2"),
    features = "Gene1",
    group.by = "condition",
    fit.by = "condition",
    compare_lineages = FALSE,
    exp_method = "raw",
    lib_normalize = FALSE,
    add_interval = FALSE,
    add_rug = FALSE,
    combine = FALSE,
    verbose = FALSE
  )

  scales <- ggplot2::ggplot_build(legend_plot)$plot$scales$scales
  has_raw_only_guide <- any(vapply(scales, function(x) {
    "C" %in% x$get_breaks() && !identical(x$guide, "none")
  }, logical(1)))
  expect_true(has_raw_only_guide)
})

test_that("DynamicPlot keeps a rug guide for globally raw-only fit groups", {
  cells <- paste0("cell", seq_len(25))
  counts <- Matrix::Matrix(rbind(Gene1 = seq_len(25)), sparse = TRUE)
  colnames(counts) <- cells
  srt <- Seurat::CreateSeuratObject(counts)
  srt$condition <- c(rep("A", 10), rep("B", 10), rep("C", 5))
  srt$Lineage1 <- c(
    seq(0, 1, length.out = 10),
    seq(0, 1, length.out = 10),
    rep(NA_real_, 5)
  )
  srt$Lineage2 <- c(
    seq(0, 1, length.out = 10),
    rep(NA_real_, 10),
    seq(0, 1, length.out = 5)
  )
  legend_plot <- NULL

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(
        srt,
        lineages,
        features,
        fitted_by_group = c(A = 1, B = 2, C = 3)
      )
    },
    get_legend = function(plot) {
      legend_plot <<- plot
      grid::nullGrob()
    },
    .package = "scop"
  )

  DynamicPlot(
    srt,
    lineages = c("Lineage1", "Lineage2"),
    features = "Gene1",
    group.by = "condition",
    fit.by = "condition",
    compare_lineages = FALSE,
    exp_method = "raw",
    lib_normalize = FALSE,
    add_point = FALSE,
    add_rug = TRUE,
    add_interval = FALSE,
    combine = FALSE,
    verbose = FALSE
  )

  built <- ggplot2::ggplot_build(legend_plot)$plot
  rug_layers <- Filter(
    function(x) inherits(x$geom, "GeomRug"),
    built$layers
  )
  expect_length(rug_layers, 1)
  expect_true(rug_layers[[1]]$show.legend)
  expect_true(any(vapply(built$scales$scales, function(x) {
    "C" %in% x$get_breaks() && !identical(x$guide, "none")
  }, logical(1))))
})

test_that("DynamicPlot disambiguates combined series legend labels", {
  srt <- make_dynamic_plot_test_object()
  srt@meta.data[["A - B"]] <- srt$Lineage1
  srt@meta.data[["A"]] <- srt$Lineage2
  srt@meta.data[["C"]] <- seq_len(ncol(srt))
  srt@meta.data[["B - C"]] <- rev(seq_len(ncol(srt)))

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(srt, lineages, features)
    },
    .package = "scop"
  )

  plots <- DynamicPlot(
    srt,
    lineages = c("A - B", "A"),
    features = c("C", "B - C"),
    fit.by = "condition",
    compare_lineages = TRUE,
    compare_features = TRUE,
    exp_method = "raw",
    lib_normalize = FALSE,
    add_point = FALSE,
    add_rug = FALSE,
    add_interval = FALSE,
    combine = FALSE,
    verbose = FALSE
  )

  scales <- ggplot2::ggplot_build(plots[[1]])$plot$scales$scales
  linetype_scale <- Filter(
    function(x) "linetype" %in% x$aesthetics,
    scales
  )[[1]]
  expect_length(unique(linetype_scale$get_labels()), 4)
  expect_setequal(
    linetype_scale$get_labels(),
    c(
      "Lineage: A - B; feature: C",
      "Lineage: A; feature: C",
      "Lineage: A - B; feature: B - C",
      "Lineage: A; feature: B - C"
    )
  )
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

test_that("DynamicPlot keeps cells with missing fit groups as raw points", {
  srt <- make_dynamic_plot_test_object()
  srt$condition[c("cell1", "cell13")] <- c(NA_character_, "")

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

test_that("DynamicPlot checks finite support for each fitted feature", {
  srt <- make_dynamic_plot_test_object()
  srt$score <- seq_len(ncol(srt))
  srt$score[paste0("cell", 4:12)] <- NA_real_
  fitted_features <- list()

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      group <- as.character(srt$condition[
        is.finite(srt@meta.data[[lineages]])
      ][1])
      fitted_features[[group]] <<- features
      mock_dynamic_fit_result(srt, lineages, features)
    },
    .package = "scop"
  )

  plots <- DynamicPlot(
    srt,
    lineages = "Lineage1",
    features = c("Gene1", "score"),
    fit.by = "condition",
    exp_method = "raw",
    lib_normalize = FALSE,
    add_point = FALSE,
    add_interval = FALSE,
    add_rug = FALSE,
    combine = FALSE,
    verbose = FALSE
  )

  expect_identical(fitted_features[["control"]], "Gene1")
  expect_setequal(fitted_features[["HUA"]], c("Gene1", "score"))
  score_line <- dynamic_plot_built_layer(plots[["Lineage1.score"]], "line")
  expect_equal(length(unique(score_line$colour)), 1)
})

test_that("DynamicPlot appends fit.by after the positional API", {
  srt <- make_dynamic_plot_test_object()

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(srt, lineages, features)
    },
    .package = "scop"
  )

  expect_identical(names(formals(DynamicPlot))[[6]], "cells")
  expect_identical(tail(names(formals(DynamicPlot)), 1), "fit.by")
  expect_no_error(do.call(
    DynamicPlot,
    list(
      srt,
      "Lineage1",
      "Gene1",
      NULL,
      NULL,
      colnames(srt),
      add_point = FALSE,
      add_interval = FALSE,
      add_rug = FALSE,
      combine = FALSE,
      verbose = FALSE
    )
  ))
})

test_that("DynamicPlot keeps grouped series keys collision safe", {
  srt <- make_dynamic_plot_test_object()
  srt[["A-B"]] <- seq_len(ncol(srt))
  srt[["A"]] <- rev(seq_len(ncol(srt)))
  srt$condition <- rep(c("C", "B-C"), each = 12)

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(
        srt,
        lineages,
        features,
        fitted_by_group = c(C = 1, `B-C` = 2)
      )
    },
    .package = "scop"
  )

  plots <- DynamicPlot(
    srt,
    lineages = "Lineage1",
    features = c("A-B", "A"),
    fit.by = "condition",
    compare_features = TRUE,
    exp_method = "raw",
    lib_normalize = FALSE,
    add_interval = FALSE,
    add_point = FALSE,
    add_rug = FALSE,
    combine = FALSE,
    verbose = FALSE
  )

  line_data <- dynamic_plot_built_layer(plots[[1]], "line")
  expect_equal(length(unique(line_data$group)), 4)
})

test_that("DynamicPlot keeps feature linetypes stable across lineage panels", {
  srt <- make_dynamic_plot_test_object()

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(
        srt,
        lineages,
        features,
        fitted_by_group = c(control = 0, HUA = 10)
      )
    },
    .package = "scop"
  )

  plots <- DynamicPlot(
    srt,
    lineages = c("Lineage1", "Lineage2"),
    features = c("Gene1", "Gene2"),
    fit.by = "condition",
    compare_lineages = FALSE,
    compare_features = TRUE,
    exp_method = "raw",
    lib_normalize = FALSE,
    add_interval = FALSE,
    add_point = FALSE,
    add_rug = FALSE,
    combine = FALSE,
    verbose = FALSE
  )

  linetypes <- lapply(plots, function(plot) {
    line_data <- dynamic_plot_built_layer(plot, "line")
    vapply(
      split(line_data$linetype, line_data$y),
      function(x) unique(x)[1],
      character(1)
    )
  })
  expect_identical(linetypes[[1]], linetypes[[2]])
})

test_that("DynamicPlot includes usable library sizes in feature support", {
  srt <- make_dynamic_plot_test_object()
  custom_libsize <- rep(1, ncol(srt))
  custom_libsize[10:12] <- 0
  fitted_groups <- character()

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      cells <- rownames(srt@meta.data)[is.finite(srt@meta.data[[lineages]])]
      fitted_groups <<- c(fitted_groups, as.character(srt$condition[cells][1]))
      mock_dynamic_fit_result(srt, lineages, features)
    },
    .package = "scop"
  )

  DynamicPlot(
    srt,
    lineages = "Lineage1",
    features = "Gene1",
    fit.by = "condition",
    libsize = custom_libsize,
    exp_method = "raw",
    lib_normalize = FALSE,
    add_interval = FALSE,
    add_point = FALSE,
    add_rug = FALSE,
    combine = FALSE,
    verbose = FALSE
  )

  expect_identical(fitted_groups, "HUA")
})

test_that("DynamicPlot preserves metadata named LineagesFeaturesFitGroups", {
  srt <- make_dynamic_plot_test_object()
  srt$LineagesFeaturesFitGroups <- rep(c("meta-A", "meta-B"), times = 12)

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
    group.by = "LineagesFeaturesFitGroups",
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

test_that("DynamicPlot keeps global linetypes when panel fits differ", {
  srt <- make_dynamic_plot_test_object()
  srt$condition <- "all"
  srt$Lineage1[13:24] <- NA_real_
  srt$Lineage2[1:12] <- NA_real_
  srt$score1 <- seq_len(ncol(srt))
  srt$score1[1:12] <- NA_real_
  srt$score2 <- rev(seq_len(ncol(srt)))

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, ...) {
      mock_dynamic_fit_result(
        srt,
        lineages,
        features,
        fitted_by_group = c(all = 1)
      )
    },
    .package = "scop"
  )

  for (add_line in c(TRUE, FALSE)) {
    plots <- DynamicPlot(
      srt,
      lineages = c("Lineage1", "Lineage2"),
      features = c("score1", "score2"),
      fit.by = "condition",
      compare_lineages = FALSE,
      compare_features = TRUE,
      exp_method = "raw",
      lib_normalize = FALSE,
      add_line = add_line,
      add_interval = !add_line,
      add_point = FALSE,
      add_rug = FALSE,
      combine = FALSE,
      verbose = FALSE
    )

    plot_layers <- lapply(plots, function(plot) {
      layers <- ggplot2::ggplot_build(plot)$data
      is_target <- vapply(layers, function(x) {
        if (add_line) {
          "linewidth" %in% names(x) && !"ymin" %in% names(x)
        } else {
          "ymin" %in% names(x)
        }
      }, logical(1))
      layers[[which(is_target)[1]]]
    })
    lineage1_score2 <- as.character(unique(plot_layers[[1]]$linetype))
    score2_rows <- if (add_line) {
      plot_layers[[2]]$y == 2
    } else {
      plot_layers[[2]]$group == max(plot_layers[[2]]$group)
    }
    lineage2_score2 <- as.character(unique(
      plot_layers[[2]]$linetype[score2_rows]
    ))
    expect_identical(lineage1_score2, lineage2_score2)
  }
})

test_that("DynamicPlot preserves unnamed family order before feature regrouping", {
  srt <- make_dynamic_plot_test_object()
  srt$score <- seq_len(ncol(srt))
  fitted_families <- list()

  testthat::local_mocked_bindings(
    RunDynamicFeatures = function(srt, lineages, features, family, ...) {
      cells <- rownames(srt@meta.data)[is.finite(srt@meta.data[[lineages]])]
      group <- as.character(srt$condition[cells][1])
      fitted_families[[group]] <<- family
      mock_dynamic_fit_result(srt, lineages, features)
    },
    .package = "scop"
  )

  DynamicPlot(
    srt,
    lineages = "Lineage1",
    features = c("score", "Gene1"),
    fit.by = "condition",
    family = c("gaussian", "nb"),
    exp_method = "raw",
    lib_normalize = FALSE,
    add_line = TRUE,
    add_interval = FALSE,
    add_point = FALSE,
    add_rug = FALSE,
    combine = FALSE,
    verbose = FALSE
  )

  expected <- c(Gene1 = "nb", score = "gaussian")
  expect_identical(fitted_families[["control"]], expected)
  expect_identical(fitted_families[["HUA"]], expected)
})
