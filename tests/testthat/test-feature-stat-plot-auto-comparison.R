make_feature_stat_meta <- function(score, group, group_levels = unique(group)) {
  cells <- paste0("Cell", seq_along(score))
  data.frame(
    cells = cells,
    score = score,
    group = factor(group, levels = group_levels),
    row.names = cells
  )
}

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

  expect_no_error(ggplot2::ggplot_build(plots[[1L]]))
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

  plots <- ExpressionStatPlot(
    meta.data = meta,
    stat.by = "score",
    group.by = "group",
    split.by = "split",
    plot_type = "box",
    auto_comparison = TRUE,
    sig_label = "p.format",
    force = TRUE
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

  plots <- ExpressionStatPlot(
    meta.data = meta,
    stat.by = "score",
    group.by = "group",
    plot_type = "box",
    comparisons = list(c("A", "C")),
    auto_comparison = TRUE,
    sig_label = "p.format",
    force = TRUE
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
