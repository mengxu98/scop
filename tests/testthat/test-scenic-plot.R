make_scenic_plot_mock <- function(seed = 1) {
  set.seed(seed)
  counts <- Matrix::Matrix(matrix(
    stats::rpois(8 * 12, lambda = 5),
    nrow = 8,
    ncol = 12,
    dimnames = list(
      paste0("Gene", seq_len(8)),
      paste0("Cell", seq_len(12))
    )
  ), sparse = TRUE)
  counts[1, ] <- seq_len(ncol(counts))
  counts[2, ] <- rev(seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$CellType <- rep(c("A", "B", "C"), each = 4)

  auc <- matrix(
    stats::runif(6 * 12),
    nrow = 6,
    ncol = 12,
    dimnames = list(
      paste0(c("Jun", "Atf3", "Fos", "Klf4", "Sox9", "Birc3"), "(+)"),
      colnames(srt)
    )
  )
  srt@tools$SCENIC <- list(scores_cells_by_regulon = t(auc))
  list(srt = srt, auc = auc)
}

test_that("SCENICPlot rss heatmap returns drawable plot object", {
  dat <- make_scenic_plot_mock()
  out <- SCENICPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "rss_heatmap",
    features = rownames(dat$auc)[1:4],
    verbose = FALSE
  )

  expect_s3_class(out$plot, "ggplot")
  expect_s3_class(out$plots[[1]], "ggplot")
  expect_false("matrix_list" %in% names(out$plot))
  expect_s3_class(out$heatmap$plot, "ggplot")
  expect_true("matrix_list" %in% names(out$heatmap))
})

test_that("SCENICPlot activity heatmap returns drawable plot object", {
  dat <- make_scenic_plot_mock(seed = 2)
  out <- SCENICPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "activity_heatmap",
    features = rownames(dat$auc)[1:4],
    verbose = FALSE
  )

  expect_s3_class(out$plot, "ggplot")
  expect_s3_class(out$plots[[1]], "ggplot")
  expect_false("matrix_list" %in% names(out$plot))
  expect_s3_class(out$heatmap$plot, "ggplot")
  expect_true("matrix_list" %in% names(out$heatmap))
})

test_that("SCENICPlot activity heatmap can order rows by RSS source group", {
  dat <- make_scenic_plot_mock(seed = 22)
  dat$auc[,] <- 0
  dat$auc["Jun(+)", dat$srt$CellType == "A"] <- 4
  dat$auc["Atf3(+)", dat$srt$CellType == "A"] <- 3
  dat$auc["Fos(+)", dat$srt$CellType == "B"] <- 4
  dat$auc["Klf4(+)", dat$srt$CellType == "B"] <- 3
  dat$auc["Sox9(+)", dat$srt$CellType == "C"] <- 4
  dat$auc["Birc3(+)", dat$srt$CellType == "C"] <- 3
  dat$srt@tools$SCENIC <- list(scores_cells_by_regulon = t(dat$auc))

  out <- SCENICPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "activity_heatmap",
    heatmap_order = "group",
    heatmap_cluster_rows = FALSE,
    top_n = 1,
    verbose = FALSE
  )

  regulons <- rev(levels(out$plot_data$regulon))
  expect_length(regulons, 3)
  expect_equal(
    as.character(unique(out$plot_data$rss_group[match(regulons, out$plot_data$regulon)])),
    c("A", "B", "C")
  )
  expect_true(all(c("rss_group", "rss_rank") %in% colnames(out$plot_data)))
})

test_that("SCENICPlot activity heatmap explicit features are not replaced by top_n", {
  dat <- make_scenic_plot_mock(seed = 23)
  out <- SCENICPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "activity_heatmap",
    features = rownames(dat$auc)[1:4],
    heatmap_order = "input",
    top_n = 1,
    verbose = FALSE
  )

  expect_equal(rev(levels(out$plot_data$regulon)), rownames(dat$auc)[1:4])
  expect_true(all(c("rss_group", "rss_rank") %in% colnames(out$plot_data)))
})

test_that("SCENICPlot rss rank keeps requested labels by default", {
  dat <- make_scenic_plot_mock(seed = 3)
  out <- SCENICPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "rss_rank",
    top_n = 4,
    verbose = FALSE
  )

  label_layers <- Filter(
    function(layer) inherits(layer$geom, "GeomTextRepel"),
    out$plots[[1]]$layers
  )
  expect_length(label_layers, 1)
  expect_equal(label_layers[[1]]$geom_params$max.overlaps, Inf)
  expect_equal(
    sum(out$top_table[["group"]] == unique(out$top_table[["group"]])[[1]]),
    4
  )
})

test_that("SCENICPlot activity correlation dumbbell returns plot and statistics", {
  dat <- make_scenic_plot_mock(seed = 5)
  dat$auc[1, ] <- seq_len(ncol(dat$auc))
  dat$auc[2, ] <- rev(seq_len(ncol(dat$auc)))
  dat$srt@tools$SCENIC <- list(scores_cells_by_regulon = t(dat$auc))
  dat$srt$AFP_score <- as.numeric(dat$auc[1, ])
  dat$srt$CYP_signature <- as.numeric(dat$auc[2, ])

  out <- SCENICPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "activity_cor_dumbbell",
    features = rownames(dat$auc)[1:2],
    cor.features = c("AFP_score", "CYP_signature"),
    cor.feature.labels = c("AFP expression", "CYP signature"),
    cor_label = FALSE,
    p_cutoff = 0.05,
    verbose = FALSE
  )

  expect_s3_class(out$plot, "ggplot")
  expect_s3_class(out$plots[[1]], "ggplot")
  expect_equal(nrow(out$plot_data), 4)
  expect_setequal(
    as.character(unique(out$plot_data[["target"]])),
    c("AFP expression", "CYP signature")
  )
  expect_true(all(c(
    "regulon", "TF", "target", "target_feature",
    "target_source", "cor", "p_val", "significant"
  ) %in% colnames(out$plot_data)))
  expect_equal(
    out$plot_data[["significant"]],
    out$plot_data[["p_val"]] < 0.05
  )
})

test_that("SCENICPlot activity correlation dumbbell resolves TF and gene targets", {
  dat <- make_scenic_plot_mock(seed = 6)
  rownames(dat$auc)[1:2] <- c("Jun(+)", "Jun(-)")
  dat$auc[1, ] <- seq_len(ncol(dat$auc))
  dat$auc[2, ] <- rev(seq_len(ncol(dat$auc)))
  dat$srt@tools$SCENIC <- list(scores_cells_by_regulon = t(dat$auc))
  dat$srt$AFP_score <- as.numeric(dat$auc[1, ])

  out <- SCENICPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "activity_cor_dumbbell",
    features = "Jun",
    cor.features = c("AFP_score", "Gene1"),
    layer = "counts",
    cor_label = FALSE,
    p_cutoff = 0.2,
    verbose = FALSE
  )

  expect_setequal(
    as.character(unique(out$plot_data[["regulon"]])),
    c("Jun(+)", "Jun(-)")
  )
  expect_equal(nrow(out$plot_data), 4)
  expect_true("metadata" %in% out$plot_data[["target_source"]])
  expect_true("assay:RNA" %in% out$plot_data[["target_source"]])
  expect_equal(
    out$plot_data[["significant"]],
    out$plot_data[["p_val"]] < 0.2
  )
})

test_that("SCENICPlot resolves positive and negative regulon suffixes", {
  dat <- make_scenic_plot_mock(seed = 4)
  rownames(dat$auc)[1:2] <- c("Jun(+)", "Jun(-)")
  dat$srt@tools$SCENIC <- list(scores_cells_by_regulon = t(dat$auc))

  out <- SCENICPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "rss_dotplot",
    features = "Jun",
    verbose = FALSE
  )

  expect_setequal(as.character(unique(out$plot_data[["regulon"]])), c("Jun(-)", "Jun(+)"))
  expect_equal(
    sort(unique(out$rank_table[out$rank_table[["regulon"]] %in% c("Jun(+)", "Jun(-)"), "TF"])),
    "Jun"
  )
})
