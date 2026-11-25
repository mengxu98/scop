make_feature_stat_meta <- function(score, group, group_levels = unique(group)) {
  cells <- paste0("Cell", seq_along(score))
  data.frame(
    cells = cells,
    score = score,
    group = factor(group, levels = group_levels),
    row.names = cells
  )
}

test_that("ExpressionStatPlot checks its bar summary dependency before rendering", {
  requested <- character()
  original_check <- get("check_r", asNamespace("scop"))
  testthat::local_mocked_bindings(
    check_r = function(packages, ...) {
      requested <<- c(requested, packages)
      if ("Hmisc" %in% packages) stop("bar summary dependency requested")
      original_check(packages, ...)
    },
    .package = "scop"
  )
  meta <- make_feature_stat_meta(seq_len(12), rep(c("A", "B"), each = 6))
  expect_no_error(ExpressionStatPlot(
    meta.data = meta, stat.by = "score", group.by = "group",
    plot_type = "box", force = TRUE
  ))
  expect_false("Hmisc" %in% requested)
  expect_error(ExpressionStatPlot(
    meta.data = meta, stat.by = "score", group.by = "group",
    plot_type = "bar", force = TRUE
  ), "bar summary dependency requested")
})

test_that("ExpressionStatPlot retains fill colours without unused colour warnings", {
  meta <- make_feature_stat_meta(seq_len(12), rep(c("A", "B"), each = 6))
  colors <- c(A = "#123456", B = "#FEDCBA")
  for (type in c("violin", "box", "bar", "dot", "col")) {
    plot <- ExpressionStatPlot(
      meta.data = meta, stat.by = "score", group.by = "group",
      plot_type = type, palcolor = colors, force = TRUE
    )[[1L]]
    expect_no_warning(ggplot2::ggplotGrob(plot))
    built <- ggplot2::ggplot_build(plot)
    expect_identical(
      unname(built$plot$scales$get_scales("fill")$map(c("A", "B"))),
      unname(colors)
    )
  }
})

get_stat_compare_layer <- function(plot) {
  compare_layers <- Filter(
    function(layer) inherits(layer$stat, "StatSignif"),
    plot$layers
  )
  if (length(compare_layers) == 0) {
    return(NULL)
  }
  compare_layers[[1]]
}

get_boxplot_layer <- function(plot) {
  boxplot_layers <- Filter(
    function(layer) inherits(layer$geom, "GeomBoxplot"),
    plot$layers
  )
  if (length(boxplot_layers) == 0) {
    return(NULL)
  }
  boxplot_layers[[1]]
}

get_stat_label_layer <- function(plot) {
  built <- suppressWarnings(ggplot2::ggplot_build(plot))
  label_layers <- Filter(function(data) "label" %in% names(data), built$data)
  if (length(label_layers) == 0) {
    return(NULL)
  }
  label_layers[[1]]
}

test_that("ExpressionStatPlot honors box_width for standalone box plots", {
  meta <- make_feature_stat_meta(
    score = seq_len(6),
    group = rep(c("A", "B"), each = 3)
  )

  narrow_plot <- ExpressionStatPlot(
    meta.data = meta,
    stat.by = "score",
    group.by = "group",
    plot_type = "box",
    box_width = 0.2,
    force = TRUE
  )[[1]]
  wide_plot <- ExpressionStatPlot(
    meta.data = meta,
    stat.by = "score",
    group.by = "group",
    plot_type = "box",
    box_width = 0.9,
    force = TRUE
  )[[1]]

  expect_equal(get_boxplot_layer(narrow_plot)$aes_params$width, 0.2)
  expect_equal(get_boxplot_layer(wide_plot)$aes_params$width, 0.9)
})

test_that("ExpressionStatPlot supports sparse co-expression", {
  cells <- paste0("Cell", seq_len(8))
  expr <- methods::as(Matrix::Matrix(
    rbind(
      gene1 = c(1, 0, 2, 0, 3, 0, 1, 0),
      gene2 = c(2, 0, 1, 0, 4, 0, 2, 0)
    ),
    sparse = TRUE,
    dimnames = list(c("gene1", "gene2"), cells)
  ), "dgCMatrix")
  meta <- data.frame(
    cells = cells,
    group = rep(c("A", "B"), each = 4),
    row.names = cells
  )

  expect_no_error(ExpressionStatPlot(
    exp.data = expr,
    meta.data = meta,
    stat.by = c("gene1", "gene2"),
    group.by = "group",
    calculate_coexp = TRUE,
    plot_type = "box",
    force = TRUE
  ))
})

test_that("auto_comparison selects the group with the highest median", {
  skip_if_not_installed("ggpubr")

  meta <- make_feature_stat_meta(
    score = c(1, 2, 3, 8, 9, 10, 4, 5, 6),
    group = rep(c("A", "B", "C"), each = 3),
    group_levels = c("A", "B", "C")
  )

  plots <- ExpressionStatPlot(
    meta.data = meta,
    stat.by = "score",
    group.by = "group",
    plot_type = "box",
    auto_comparison = TRUE,
    sig_label = "p.format",
    force = TRUE
  )

  compare_layer <- get_stat_compare_layer(plots[[1]])
  expect_equal(
    compare_layer$stat_params$comparisons,
    list(c("B", "A"), c("B", "C"))
  )
})

test_that("multiple-group p.format labels build with current ggpubr", {
  skip_if_not_installed("ggpubr")

  meta <- make_feature_stat_meta(
    score = c(1, 2, 3, 8, 9, 10, 4, 5, 6),
    group = rep(c("A", "B", "C"), each = 3),
    group_levels = c("A", "B", "C")
  )
  plots <- ExpressionStatPlot(
    meta.data = meta,
    stat.by = "score",
    group.by = "group",
    plot_type = "box",
    multiplegroup_comparisons = TRUE,
    sig_label = "p.format",
    force = TRUE
  )

  label_data <- get_stat_label_layer(plots[[1L]])
  expect_false(is.null(label_data))
  expect_equal(as.character(label_data$method), "Kruskal-Wallis")
  expect_equal(
    label_data$p,
    stats::kruskal.test(score ~ group, data = meta)$p.value,
    tolerance = 1e-8
  )
  expect_identical(label_data$label, paste0("p = ", label_data$p.format))
  expect_identical(label_data$label, "p = 0.027")
})

test_that("pairwise p.format labels build with current ggpubr", {
  skip_if_not_installed("ggpubr")

  meta <- make_feature_stat_meta(
    score = c(1, 2, 3, 8, 9, 10, 4, 5, 6),
    group = rep(c("A", "B", "C"), each = 3),
    group_levels = c("A", "B", "C")
  )
  plots <- ExpressionStatPlot(
    meta.data = meta,
    stat.by = "score",
    group.by = "group",
    plot_type = "box",
    comparisons = list(c("A", "B"), c("A", "C")),
    sig_label = "p.format",
    force = TRUE
  )

  expect_no_error(ggplot2::ggplot_build(plots[[1L]]))
})

test_that("auto_comparison breaks median ties by mean", {
  skip_if_not_installed("ggpubr")

  meta <- make_feature_stat_meta(
    score = c(1, 4, 6, 9, 2, 3, 7, 16, 0, 0.5, 1, 1.5),
    group = rep(c("A", "B", "C"), each = 4),
    group_levels = c("A", "B", "C")
  )

  plots <- ExpressionStatPlot(
    meta.data = meta,
    stat.by = "score",
    group.by = "group",
    plot_type = "box",
    auto_comparison = TRUE,
    sig_label = "p.format",
    force = TRUE
  )

  compare_layer <- get_stat_compare_layer(plots[[1]])
  expect_equal(
    compare_layer$stat_params$comparisons,
    list(c("B", "A"), c("B", "C"))
  )
})

test_that("auto_comparison does not add comparisons with fewer than two valid groups", {
  skip_if_not_installed("ggpubr")

  meta <- make_feature_stat_meta(
    score = c(1, 2, NA, NA),
    group = rep(c("A", "B"), each = 2),
    group_levels = c("A", "B")
  )

  plots <- ExpressionStatPlot(
    meta.data = meta,
    stat.by = "score",
    group.by = "group",
    plot_type = "box",
    auto_comparison = TRUE,
    sig_label = "p.format",
    force = TRUE
  )

  expect_null(get_stat_compare_layer(plots[[1]]))
})

test_that("auto_comparison is disabled for split comparisons", {
  meta <- make_feature_stat_meta(
    score = c(1, 2, 3, 8, 9, 10, 4, 5, 6),
    group = rep(c("A", "B", "C"), each = 3),
    group_levels = c("A", "B", "C")
  )
  meta$split <- factor(rep(c("S1", "S2", "S1"), times = 3))

  expect_warning(
    plots <- ExpressionStatPlot(
      meta.data = meta,
      stat.by = "score",
      group.by = "group",
      split.by = "split",
      plot_type = "box",
      auto_comparison = TRUE,
      sig_label = "p.format",
      force = TRUE
    ),
    "only supported"
  )

  expect_null(get_stat_compare_layer(plots[[1]]))
})

test_that("explicit comparisons take precedence over auto_comparison", {
  skip_if_not_installed("ggpubr")

  meta <- make_feature_stat_meta(
    score = c(1, 2, 3, 8, 9, 10, 4, 5, 6),
    group = rep(c("A", "B", "C"), each = 3),
    group_levels = c("A", "B", "C")
  )

  expect_warning(
    plots <- ExpressionStatPlot(
      meta.data = meta,
      stat.by = "score",
      group.by = "group",
      plot_type = "box",
      comparisons = list(c("A", "C")),
      auto_comparison = TRUE,
      sig_label = "p.format",
      force = TRUE
    ),
    "explicit.*comparisons"
  )

  compare_layer <- get_stat_compare_layer(plots[[1]])
  expect_equal(compare_layer$stat_params$comparisons, list(c("A", "C")))
})

test_that("explicit ref_group is used by auto_comparison", {
  skip_if_not_installed("ggpubr")

  meta <- make_feature_stat_meta(
    score = c(1, 2, 3, 8, 9, 10, 4, 5, 6),
    group = rep(c("A", "B", "C"), each = 3),
    group_levels = c("A", "B", "C")
  )

  plots <- suppressMessages(ExpressionStatPlot(
    meta.data = meta,
    stat.by = "score",
    group.by = "group",
    plot_type = "box",
    ref_group = "A",
    auto_comparison = TRUE,
    sig_label = "p.format",
    force = TRUE
  ))

  compare_layer <- get_stat_compare_layer(plots[[1]])
  expect_equal(
    compare_layer$stat_params$comparisons,
    list(c("A", "B"), c("A", "C"))
  )
})

test_that("expression fill retains per-split medians and theme callbacks", {
  cells <- paste0("c", seq_len(48L))
  meta <- data.frame(
    cells = cells,
    group = factor(rep(c("B", "A"), each = 24L), levels = c("A", "B")),
    batch = factor(rep(rep(c("S2", "S1"), each = 12L), 2L)),
    a = rep(seq_len(12L), 4L),
    b = seq_along(cells) / 10,
    row.names = cells
  )
  seen <- numeric()
  theme_observer <- function(...) {
    seen <<- c(seen, stats::runif(1))
    ggplot2::theme_bw()
  }
  withr::local_seed(21L)
  rng_before <- .Random.seed
  plots <- ExpressionStatPlot(
    meta.data = meta, stat.by = c("a", "b"), group.by = "group",
    split.by = "batch", fill.by = "expression", palette = "Blues",
    individual = TRUE, plot_type = "box", theme_use = theme_observer
  )
  expect_length(plots, 8L)
  expect_length(seen, length(plots))
  expect_identical(.Random.seed, rng_before)
  for (panel in plots) {
    d <- panel$data
    f <- unique(d$features)
    expected <- stats::median(meta[[f]][meta$group == as.character(d$group.by[[1L]]) &
      meta$batch == as.character(d$split.by[[1L]])])
    expect_equal(d$fill.by, rep(expected, nrow(d)))
  }
})

test_that("grouped expression plots accept several grouping columns and omit absent individual panels", {
  meta <- make_feature_stat_meta(seq_len(12L), rep(c("A", "B"), each = 6L))
  meta$other <- rep(c("C", "D"), each = 3L, length.out = nrow(meta))
  meta$batch <- factor(rep(c("u", "v"), each = 6L), levels = c("v", "u"))
  grouped <- ExpressionStatPlot(
    meta.data = meta, stat.by = "score", group.by = c("group", "other"),
    split.by = "batch", plot_type = "box"
  )
  expect_named(grouped, c("score:group:A,B:v,u", "score:other:C,D:v,u"))
  expect_identical(lapply(grouped, function(p) sort(p$data$score)),
    list("score:group:A,B:v,u" = meta$score, "score:other:C,D:v,u" = meta$score)
  )

  individual <- ExpressionStatPlot(
    meta.data = meta, stat.by = "score", group.by = "group", split.by = "batch",
    plot_type = "box", individual = TRUE, fill.by = "expression", palette = "Blues"
  )
  expect_named(individual, c("score:group:B:v", "score:group:A:u"))
  expect_identical(lapply(individual, function(p) unique(p$data$fill.by)),
    list("score:group:B:v" = 9.5, "score:group:A:u" = 3.5)
  )
  expect_no_warning(invisible(lapply(individual, ggplot2::ggplot_build)))
})
