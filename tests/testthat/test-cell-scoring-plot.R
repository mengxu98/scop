make_scoreplot_srt <- function(n_cells = 24L, n_genes = 6L) {
  counts <- matrix(rpois(n_genes * n_cells, lambda = 5), nrow = n_genes)
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  srt <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  srt$celltype <- factor(rep(c("A", "B", "C"), length.out = n_cells))
  srt$AUC_one <- seq(0.01, 0.49, length.out = n_cells)
  srt$AUC_two <- rev(srt$AUC_one)
  angle <- seq(0, 2 * pi, length.out = n_cells + 1L)[-1L]
  embedding <- cbind(UMAP_1 = cos(angle), UMAP_2 = sin(angle))
  rownames(embedding) <- colnames(srt)
  srt[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    key = "UMAP_",
    assay = SeuratObject::DefaultAssay(srt)
  )
  srt
}

test_that("CellScoringPlot returns group/feature component pairs", {
  srt <- make_scoreplot_srt()
  original_meta <- srt@meta.data

  plots <- CellScoringPlot(
    srt,
    features = c("Program one" = "AUC_one", "Program two" = "AUC_two"),
    group.by = "celltype",
    reduction = "umap",
    add_bar = FALSE,
    point.fraction = 1,
    combine = FALSE,
    verbose = FALSE
  )

  expect_length(plots, 6L)
  expect_named(
    plots,
    c(
      "group_umap", "group_legend",
      "umap_Program one", "umap_Program two",
      "stats_Program one", "stats_Program two"
    )
  )
  expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))
  expect_equal(nrow(ggplot2::layer_data(plots[["stats_Program one"]], 2)), ncol(srt))
  expect_identical(srt@meta.data, original_meta)

  combined <- CellScoringPlot(
    srt,
    features = c("AUC_one", "AUC_two"),
    group.by = "celltype",
    reduction = "umap",
    thresholds = 0.25,
    ncol = 3,
    nrow = 2,
    verbose = FALSE
  )
  expect_s3_class(combined, "patchwork")
  expect_no_error(patchwork::patchworkGrob(combined))
})

test_that("CellScoringPlot infers a compact default layout", {
  expect_null(formals(CellScoringPlot)$ncol)
  expect_identical(formals(CellScoringPlot)$add_bar, TRUE)
  expect_identical(
    eval(formals(CellScoringPlot)$highlight.mark.type),
    c("ellipse", "hull")
  )

  srt <- make_scoreplot_srt()
  scores <- cbind(
    Alpha = srt$AUC_one,
    Beta = srt$AUC_two,
    Ductal = (srt$AUC_one + srt$AUC_two) / 2
  )
  plot <- CellScoringPlot(
    srt,
    scores = scores,
    features = colnames(scores),
    group.by = "celltype",
    add_bar = FALSE,
    verbose = FALSE
  )

  expect_identical(plot$patches$layout$ncol, 4L)
  expect_equal(plot$patches$layout$nrow, 2L)
})

test_that("CellScoringPlot infers highlights without redundant display parameters", {
  args <- names(formals(CellScoringPlot))
  expect_true("method" %in% args)
  expect_false("highlight.groups" %in% args)
  expect_false("theme.base.size" %in% args)
  expect_false("feature.palette" %in% args)
  expect_false("feature.palcolor" %in% args)
  expect_false("legend.labels" %in% args)
  expect_false(any(c(
    "auc", "auc_palette", "auc_palcolor", "show.auc.legend"
  ) %in% args))
  expect_true(all(c(
    "scores", "score.palette", "score.palcolor", "show.score.legend"
  ) %in% args))
  expect_true(all(c("bar.text.size", "bar.text.color") %in% args))

  highlights <- scoreplot_highlights(
    labels = c("Display A", "B", "Unmatched"),
    sources = c("A", "score_b", "score_c"),
    groups = c("A", "B", "C")
  )
  expect_identical(
    highlights,
    list("Display A" = "A", "B" = "B", "Unmatched" = character())
  )

  groups <- factor(rep(c("A", "B", "C"), c(1000, 20, 3)))
  expect_identical(
    scoreplot_group_labels(
      groups,
      labels = c("Display A", "B"),
      sources = c("A", "score_b")
    ),
    c(A = "Display A (n = 1,000)", B = "B (n = 20)", C = "C (n = 3)")
  )
})

test_that("CellScoringPlot discovers the latest CellScoring AUCell result", {
  srt <- make_scoreplot_srt(n_genes = 120L)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  programs <- list(
    "Program A" = paste0("gene", 1:35),
    "Program B" = paste0("gene", 36:70)
  )
  scored <- suppressWarnings(CellScoring(
    srt,
    features = programs,
    method = "AUCell",
    backend = "cpp",
    classification = FALSE,
    name = "score",
    verbose = FALSE
  ))

  record <- tail(scored@misc$CellScoring$records, 1L)[[1L]]
  expect_identical(record$method, "AUCell")
  expect_identical(
    record$features,
    c("Program A" = "score_Program.A", "Program B" = "score_Program.B")
  )

  plots <- CellScoringPlot(
    scored,
    group.by = "celltype",
    reduction = "umap",
    add_box = FALSE,
    add_bar = FALSE,
    combine = FALSE,
    verbose = FALSE
  )
  expect_named(
    plots,
    c(
      "group_umap", "group_legend",
      "umap_Program A", "umap_Program B",
      "stats_Program A", "stats_Program B"
    )
  )
  expect_identical(plots[["umap_Program A"]]$labels$title, "Program A")

  assay_scored <- suppressWarnings(CellScoring(
    srt,
    features = programs,
    method = "AUCell",
    backend = "cpp",
    classification = FALSE,
    name = "program_auc",
    new_assay = TRUE,
    store_metadata = FALSE,
    verbose = FALSE
  ))
  expect_false(any(
    c("program_auc_Program.A", "program_auc_Program.B") %in%
      colnames(assay_scored@meta.data)
  ))
  assay_plots <- CellScoringPlot(
    assay_scored,
    group.by = "celltype",
    reduction = "umap",
    add_box = FALSE,
    add_bar = FALSE,
    combine = FALSE,
    verbose = FALSE
  )
  expect_true(all(c("umap_Program A", "umap_Program B") %in% names(
    assay_plots
  )))
})

test_that("CellScoringPlot selects and labels the requested scoring method", {
  srt <- make_scoreplot_srt(n_genes = 120L)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  programs <- list(
    "Program A" = paste0("gene", 1:35),
    "Program B" = paste0("gene", 36:70)
  )
  scored <- suppressWarnings(CellScoring(
    srt,
    features = programs,
    method = c("AUCell", "Seurat"),
    backend = "cpp",
    classification = FALSE,
    name = "score",
    verbose = FALSE
  ))

  plots <- CellScoringPlot(
    scored,
    method = "Seurat",
    group.by = "celltype",
    reduction = "umap",
    add_box = TRUE,
    add_bar = FALSE,
    combine = FALSE,
    verbose = FALSE
  )
  expect_true(all(c("umap_Program A", "umap_Program B") %in% names(plots)))
  expect_identical(plots[["umap_Program A"]]$labels$colour, "Seurat")
  expect_identical(plots[["stats_Program A"]]$labels$y, "Seurat")

  explicit_scores <- as.matrix(scored@meta.data[, unname(
    tail(scored@misc$CellScoring$records, 1L)[[1L]]$features
  ), drop = FALSE])
  explicit_plots <- CellScoringPlot(
    scored,
    method = "Seurat",
    scores = explicit_scores,
    group.by = "celltype",
    reduction = "umap",
    add_box = FALSE,
    add_bar = FALSE,
    combine = FALSE,
    verbose = FALSE
  )
  expect_identical(
    explicit_plots[["umap_score_Seurat_Program.A"]]$labels$colour,
    "Seurat"
  )

  expect_error(
    CellScoringPlot(
      scored,
      method = "Seurat",
      group.by = "celltype",
      reduction = "umap",
      add_bar = TRUE,
      combine = FALSE,
      verbose = FALSE
    ),
    "Automatic thresholds are only available for AUCell"
  )
})

test_that("CellScoringPlot preserves paired visual primitives", {
  srt <- make_scoreplot_srt()
  plots <- CellScoringPlot(
    srt,
    features = c("A" = "AUC_one"),
    group.by = "celltype",
    reduction = "umap",
    thresholds = 0.25,
    show.score.legend = TRUE,
    combine = FALSE,
    verbose = FALSE
  )

  score_umap <- plots[["umap_A"]]
  expect_equal(score_umap$labels$title, "A")
  expect_equal(score_umap$labels$colour, "AUC")
  expect_match(score_umap$theme$plot.title$colour, "^#[0-9A-Fa-f]{6}$")
  expect_equal(score_umap$theme$legend.position, "bottom")
  expect_gte(length(score_umap$layers), 2L)
  expect_false(any(vapply(
    score_umap$layers,
    function(x) grepl("Text|Label", class(x$geom)[[1]]),
    logical(1)
  )))
  score_layers <- ggplot2::ggplot_build(score_umap)$data
  score_scattermore <- vapply(
    score_umap$layers,
    function(layer) inherits(layer$geom, "GeomScattermore"),
    logical(1)
  )
  expect_false(any(score_scattermore))
  score_point_layers <- vapply(
    score_umap$layers,
    function(layer) inherits(layer$geom, "GeomPoint"),
    logical(1)
  )
  expect_true(any(score_point_layers))
  expect_equal(
    sum(vapply(score_layers[score_point_layers], nrow, integer(1))),
    ncol(srt)
  )
  expect_true(is.numeric(score_umap$data$group.by))

  score_point_sizes <- vapply(
    score_umap$layers[score_point_layers],
    function(layer) layer$aes_params$size,
    numeric(1)
  )
  expect_identical(score_umap$theme$legend.direction, "horizontal")
  bottom_guide <- score_umap$guides$guides$colour$params$theme
  expect_gt(
    as.numeric(bottom_guide$legend.key.width),
    as.numeric(bottom_guide$legend.key.height)
  )

  group_umap <- plots[["group_umap"]]
  expect_equal(group_umap$theme$legend.position, "none")
  group_scattermore <- vapply(
    group_umap$layers,
    function(layer) inherits(layer$geom, "GeomScattermore"),
    logical(1)
  )
  expect_false(any(group_scattermore))
  group_point_layers <- vapply(
    group_umap$layers,
    function(layer) inherits(layer$geom, "GeomPoint"),
    logical(1)
  )
  expect_true(any(group_point_layers))
  group_point_sizes <- vapply(
    group_umap$layers[group_point_layers],
    function(layer) layer$aes_params$size,
    numeric(1)
  )
  expect_identical(group_point_sizes, score_point_sizes)

  group_legend <- plots[["group_legend"]]
  expect_s3_class(group_legend, "ggplot")
  expect_lte(group_legend$layers[[1]]$aes_params$size, 2.5)
  expect_lte(group_legend$layers[[2]]$aes_params$size, 2.1)
})

test_that("CellScoringPlot applies separate main and AUC UMAP themes", {
  srt <- make_scoreplot_srt()
  main_theme <- function(base_size = 11) {
    ggplot2::theme_void(base_size) +
      ggplot2::theme(plot.background = ggplot2::element_rect(fill = "pink"))
  }
  auc_theme <- function(base_size = 11) {
    ggplot2::theme_minimal(base_size) +
      ggplot2::theme(plot.background = ggplot2::element_rect(fill = "lightblue"))
  }
  plots <- CellScoringPlot(
    srt,
    features = c("A" = "AUC_one"),
    group.by = "celltype",
    thresholds = 0.25,
    theme_use = main_theme,
    auc_theme_use = auc_theme,
    combine = FALSE,
    verbose = FALSE
  )

  expect_identical(plots[["group_umap"]]$theme$plot.background$fill, "pink")
  expect_identical(plots[["umap_A"]]$theme$plot.background$fill, "lightblue")
})

test_that("CellScoringPlot exposes bar text size and color", {
  srt <- make_scoreplot_srt()
  bar <- CellScoringPlot(
    srt,
    features = c("A" = "AUC_one"),
    group.by = "celltype",
    thresholds = 0.25,
    add_box = FALSE,
    add_bar = TRUE,
    bar.text.size = 4.1,
    bar.text.color = "navy",
    combine = FALSE,
    verbose = FALSE
  )[["stats_A"]]

  expect_equal(bar$layers[[2]]$aes_params$size, 4.1)
  expect_identical(bar$layers[[2]]$aes_params$colour, "navy")
})

test_that("CellScoringPlot orients collected legends from their position", {
  srt <- make_scoreplot_srt()
  for (position in c("bottom", "top", "left", "right")) {
    plot <- CellScoringPlot(
      srt,
      features = c("A" = "AUC_one"),
      group.by = "celltype",
      thresholds = 0.25,
      legend.position = position,
      combine = FALSE,
      verbose = FALSE
    )[["umap_A"]]

    horizontal <- position %in% c("bottom", "top")
    expect_identical(plot$theme$legend.position, position)
    expect_identical(
      plot$theme$legend.direction,
      if (horizontal) "horizontal" else "vertical"
    )
    guide_theme <- plot$guides$guides$colour$params$theme
    expect_identical(
      as.numeric(guide_theme$legend.key.width) >
        as.numeric(guide_theme$legend.key.height),
      horizontal
    )
    expect_identical(
      guide_theme$legend.title.position,
      if (horizontal) "left" else "top"
    )
    expect_identical(
      guide_theme$legend.title$vjust,
      if (horizontal) 1 else 0.5
    )
    expect_identical(
      as.numeric(guide_theme$legend.text$margin[[1]]),
      if (horizontal) 1 else 0
    )
  }
})

test_that("CellScoringPlot anchors combined titles to UMAP panels", {
  srt <- make_scoreplot_srt()
  plot <- CellScoringPlot(
    srt,
    features = c("A" = "AUC_one"),
    group.by = "celltype",
    thresholds = 0.25,
    combine = FALSE,
    verbose = FALSE
  )[["umap_A"]]
  anchored <- scoreplot_panel_title(plot)

  expect_null(anchored$labels$title)
  expect_identical(anchored$coordinates$clip, "off")
  expect_identical(class(anchored$coordinates), class(plot$coordinates))
  title_grob <- anchored$layers[[length(anchored$layers)]]$geom_params$grob
  expect_equal(title_grob$x, grid::unit(0, "npc"))
  expect_equal(title_grob$y, grid::unit(1, "npc") + grid::unit(2, "pt"))
  expect_identical(title_grob$just, c("left", "bottom"))
  expect_identical(anchored$theme$panel.border$linewidth, 0.5)
})

test_that("CellScoringPlot uses compact threshold legends", {
  srt <- make_scoreplot_srt()
  bar <- CellScoringPlot(
    srt,
    features = c("A" = "AUC_one"),
    group.by = "celltype",
    thresholds = 0.25,
    add_box = FALSE,
    add_bar = TRUE,
    combine = FALSE,
    verbose = FALSE
  )[["stats_A"]]

  expect_lte(as.numeric(bar$theme$legend.key.size), 7)
  expect_lte(bar$theme$legend.text$size, 6)
  expect_lte(bar$theme$legend.title$size, 7)
})

test_that("CellScoringPlot supports benchmark-style hull highlights", {
  srt <- make_scoreplot_srt()
  plots <- CellScoringPlot(
    srt,
    features = c("A" = "AUC_one"),
    group.by = "celltype",
    reduction = "umap",
    thresholds = 0.25,
    highlight.mark.type = "hull",
    combine = FALSE,
    verbose = FALSE
  )

  geoms <- vapply(
    plots[["umap_A"]]$layers,
    function(layer) class(layer$geom)[[1]],
    character(1)
  )
  expect_true("GeomMarkHull" %in% geoms)
})

test_that("CellScoringPlot controls UMAP and statistic row proportions", {
  srt <- make_scoreplot_srt()
  combined <- CellScoringPlot(
    srt,
    features = c("AUC_one", "AUC_two"),
    group.by = "celltype",
    reduction = "umap",
    thresholds = 0.25,
    ncol = 3,
    nrow = 2,
    row.heights = c(0.47, 0.53),
    verbose = FALSE
  )

  expect_equal(combined$patches$layout$heights, c(0.47, 0.53))
  expect_error(
    CellScoringPlot(
      srt,
      features = "AUC_one",
      group.by = "celltype",
      thresholds = 0.25,
      row.heights = c(1, 0),
      verbose = FALSE
    ),
    "two positive finite values"
  )
})

test_that("CellScoringPlot uses theme_scop without default in-situ labels", {
  srt <- make_scoreplot_srt()
  plots <- CellScoringPlot(
    srt,
    features = c("A" = "AUC_one"),
    group.by = "celltype",
    thresholds = 0.25,
    combine = FALSE,
    verbose = FALSE
  )
  group_umap <- plots[["group_umap"]]
  group_geoms <- vapply(group_umap$layers, function(x) class(x$geom)[[1]], character(1))
  expect_false(any(grepl("Text|Label", group_geoms)))
  expect_s3_class(group_umap$theme$axis.line, "element_blank")
})

test_that("CellScoringPlot has an independent continuous AUCell palette", {
  srt <- make_scoreplot_srt()
  default_plot <- CellScoringPlot(
    srt,
    features = c("A" = "AUC_one"),
    group.by = "celltype",
    thresholds = 0.25,
    combine = FALSE,
    verbose = FALSE
  )[["umap_A"]]
  default_scale <- default_plot$scales$get_scales("colour")
  expected <- unname(thisplot::palette_colors(
    n = 1000, palette = "Spectral", type = "continuous", matched = FALSE
  ))
  expect_equal(default_scale$palette(c(0, 1)), expected[c(1, length(expected))])

  custom_plot <- CellScoringPlot(
    srt,
    features = "AUC_one",
    group.by = "celltype",
    thresholds = 0.25,
    score.palette = "viridis",
    score.palcolor = c("white", "navy"),
    combine = FALSE,
    verbose = FALSE
  )[["umap_AUC_one"]]
  custom_scale <- custom_plot$scales$get_scales("colour")
  expect_equal(custom_scale$palette(c(0, 1)), c("#FFFFFF", "#000080"))
})

test_that("CellScoringPlot uses matched group colors for feature titles and outlines", {
  srt <- make_scoreplot_srt()
  plots <- CellScoringPlot(
    srt,
    features = c("A" = "AUC_one"),
    group.by = "celltype",
    reduction = "umap",
    thresholds = 0.25,
    group.palcolor = c(C = "#333333", A = "#111111", B = "#222222"),
    combine = FALSE,
    verbose = FALSE
  )

  expect_equal(
    plots[["umap_A"]]$theme$plot.title$colour,
    "#111111"
  )
  umap_data <- ggplot2::ggplot_build(plots[["umap_A"]])$data
  expect_true(any(vapply(
    umap_data,
    function(layer) any(layer$colour == "#111111"),
    logical(1)
  )))
})

test_that("CellScoringPlot accepts both score-matrix orientations", {
  srt <- make_scoreplot_srt()
  feature_by_cell <- rbind(program_a = srt$AUC_one, program_b = srt$AUC_two)
  colnames(feature_by_cell) <- colnames(srt)

  from_columns <- CellScoringPlot(
    srt,
    scores = feature_by_cell,
    features = c("A" = "program_a", "B" = "program_b"),
    group.by = "celltype",
    thresholds = 0.25,
    combine = FALSE,
    verbose = FALSE
  )
  from_rows <- CellScoringPlot(
    srt,
    scores = t(feature_by_cell),
    features = c("A" = "program_a", "B" = "program_b"),
    group.by = "celltype",
    thresholds = 0.25,
    combine = FALSE,
    verbose = FALSE
  )

  expect_equal(
    ggplot2::layer_data(from_columns[["stats_A"]], 1)$middle,
    ggplot2::layer_data(from_rows[["stats_A"]], 1)$middle
  )
})

test_that("CellScoringPlot adds explicit threshold panels and validates paired layout", {
  srt <- make_scoreplot_srt()
  plots <- CellScoringPlot(
    srt,
    features = c("A" = "AUC_one"),
    group.by = "celltype",
    thresholds = 0.25,
    point.fraction = 0,
    combine = FALSE,
    verbose = FALSE
  )

  expect_s3_class(plots[["stats_A"]], "patchwork")
  expect_no_error(patchwork::patchworkGrob(plots[["stats_A"]]))

  expect_error(
    CellScoringPlot(
      srt,
      features = "AUC_one",
      group.by = "celltype",
      thresholds = 0.25,
      nrow = 3,
      verbose = FALSE
    ),
    "must be even"
  )
  expect_error(
    CellScoringPlot(
      srt,
      features = c("AUC_one", "AUC_two"),
      group.by = "celltype",
      thresholds = 0.25,
      ncol = 1,
      nrow = 2,
      verbose = FALSE
    ),
    "paired slots"
  )
})

test_that("CellScoringPlot forwards CellDimPlot functionality without moving legends", {
  srt <- make_scoreplot_srt()
  plots <- CellScoringPlot(
    srt,
    features = c("A" = "AUC_one"),
    group.by = "celltype",
    thresholds = 0.25,
    label = TRUE,
    label.fg = "purple",
    cells.highlight = colnames(srt)[1:2],
    cols.highlight = "orange",
    add_mark = TRUE,
    mark_type = "ellipse",
    combine = FALSE,
    verbose = FALSE
  )

  expect_s3_class(plots[["group_umap"]], "ggplot")
  expect_gte(length(plots[["group_umap"]]$layers), 5L)
  expect_equal(plots[["group_umap"]]$theme$legend.position, "none")
  expect_equal(plots[["umap_A"]]$theme$legend.position, "bottom")
})

test_that("CellScoringPlot controls box and threshold bar independently", {
  srt <- make_scoreplot_srt()

  box_only <- CellScoringPlot(
    srt,
    features = c("A" = "AUC_one"),
    group.by = "celltype",
    thresholds = 0.25,
    add_box = TRUE,
    add_bar = FALSE,
    combine = FALSE,
    verbose = FALSE
  )[["stats_A"]]
  expect_s3_class(box_only, "ggplot")
  mean_layer <- box_only$layers[[length(box_only$layers)]]
  mean_grob <- mean_layer$geom_params$grob
  expect_s3_class(mean_layer$geom, "GeomCustomAnn")
  expect_equal(mean_grob$x, grid::unit(5, "pt"))
  expect_equal(mean_grob$y, grid::unit(1, "npc") - grid::unit(5, "pt"))
  expect_equal(mean_grob$just, c("left", "top"))
  expect_equal(mean_grob$gp$col, "grey50")

  bar_only <- CellScoringPlot(
    srt,
    features = c("A" = "AUC_one"),
    group.by = "celltype",
    thresholds = 0.25,
    add_box = FALSE,
    add_bar = TRUE,
    combine = FALSE,
    verbose = FALSE
  )[["stats_A"]]
  expect_s3_class(bar_only, "ggplot")
  expect_equal(bar_only$labels$y, "Proportion")

  neither <- CellScoringPlot(
    srt,
    features = c("A" = "AUC_one"),
    group.by = "celltype",
    add_box = FALSE,
    add_bar = FALSE,
    combine = FALSE,
    verbose = FALSE
  )[["stats_A"]]
  expect_s3_class(neither, "spacer")

  skip_if_not_installed("AUCell")
  auto_bar <- suppressWarnings(CellScoringPlot(
    srt,
    features = c("A" = "AUC_one"),
    group.by = "celltype",
    add_box = FALSE,
    combine = FALSE,
    verbose = FALSE
  )[["stats_A"]])
  explicit_auto_bar <- suppressWarnings(CellScoringPlot(
    srt,
    features = c("A" = "AUC_one"),
    group.by = "celltype",
    thresholds = "auto",
    add_box = FALSE,
    combine = FALSE,
    verbose = FALSE
  )[["stats_A"]])
  expect_s3_class(auto_bar, "ggplot")
  expect_equal(
    ggplot2::layer_data(auto_bar, 1)$y,
    ggplot2::layer_data(explicit_auto_bar, 1)$y
  )
})

test_that("CellScoringPlot highlight colors are independently configurable", {
  srt <- make_scoreplot_srt()
  defaults <- CellScoringPlot(
    srt,
    features = c("A" = "AUC_one"),
    group.by = "celltype",
    thresholds = 0.25,
    combine = FALSE,
    verbose = FALSE
  )
  default_umap_data <- ggplot2::ggplot_build(defaults[["umap_A"]])$data
  expected_mark <- unname(thisplot::palette_colors(
    levels(srt$celltype),
    palette = "Chinese", type = "discrete", matched = TRUE
  )[["A"]])
  expect_true(any(vapply(
    default_umap_data, function(x) any(x$colour == expected_mark), logical(1)
  )))
  expect_false(any(vapply(
    default_umap_data, function(x) any(x$colour == "grey40"), logical(1)
  )))

  plots <- CellScoringPlot(
    srt,
    features = c("A" = "AUC_one"),
    group.by = "celltype",
    highlight.label = TRUE,
    highlight.mark.color = "orange",
    highlight.label.color = "grey60",
    highlight.mean.color = "purple",
    add_bar = FALSE,
    combine = FALSE,
    verbose = FALSE
  )
  umap_data <- ggplot2::ggplot_build(plots[["umap_A"]])$data
  expect_true(any(vapply(umap_data, function(x) any(x$colour == "orange"), logical(1))))
  expect_true(any(vapply(umap_data, function(x) any(x$colour == "grey60"), logical(1))))
  stat_layers <- plots[["stats_A"]]$layers
  expect_equal(
    stat_layers[[3L]]$aes_params$colour,
    "purple"
  )
  expect_equal(
    stat_layers[[length(stat_layers)]]$geom_params$grob$gp$col,
    "grey50"
  )
})
