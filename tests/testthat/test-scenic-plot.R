make_scenic_plot_mock <- function(seed = 1, add_embedding = FALSE) {
  set.seed(seed)
  counts <- Matrix::Matrix(matrix(
    stats::rpois(8 * 12, lambda = 5),
    nrow = 8,
    ncol = 12,
    dimnames = list(
      c("Gene1", "Gene2", "Jun", "Atf3", "Fos", "Klf4", "Sox9", "Birc3"),
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
  srt@tools$SCENIC <- list(
    scores_cells_by_regulon = t(auc),
    regulon_list = list(
      `Jun(+)` = c("Gene1", "Gene2", "Atf3"),
      `Atf3(+)` = c("Gene1", "Fos"),
      `Fos(+)` = c("Gene2", "Klf4"),
      `Klf4(+)` = c("Gene1", "Sox9"),
      `Sox9(+)` = c("Klf4", "Birc3"),
      `Birc3(+)` = c("Gene2")
    ),
    adjacency = data.frame(
      TF = c("Jun", "Jun", "Atf3", "Fos", "Klf4", "Sox9", "Birc3"),
      target = c("Gene1", "Gene2", "Fos", "Klf4", "Sox9", "Birc3", "Gene2"),
      importance = c(5, 4, 3, 2.5, 2, 3.5, 1),
      stringsAsFactors = FALSE
    )
  )
  if (isTRUE(add_embedding)) {
    emb <- cbind(
      UMAP_1 = as.numeric(as.integer(factor(srt$CellType))),
      UMAP_2 = seq_len(ncol(srt)) / ncol(srt)
    )
    rownames(emb) <- colnames(srt)
    srt[["umap"]] <- Seurat::CreateDimReducObject(
      embeddings = emb,
      key = "UMAP_",
      assay = "RNA"
    )
  }
  list(srt = srt, auc = auc)
}

make_scenicplus_plot_mock <- function(seed = 11) {
  dat <- make_scenic_plot_mock(seed = seed, add_embedding = TRUE)
  regulons <- c("Jun(+)", "Atf3(+)", "Sox9(+)")
  auc <- dat$auc[regulons, , drop = FALSE]
  scores <- Matrix::Matrix(auc, sparse = TRUE)
  regulon_list <- list(
    `Jun(+)` = c("Gene1", "Gene2", "Atf3"),
    `Atf3(+)` = c("Gene1", "Fos"),
    `Sox9(+)` = c("Klf4", "Birc3")
  )
  tf_gene <- data.frame(
    TF = c("Jun", "Jun", "Atf3", "Sox9"),
    target = c("Gene1", "Gene2", "Gene1", "Klf4"),
    importance = c(4, 3, 2, 5),
    stringsAsFactors = FALSE
  )
  peak_names <- c(
    "chr1-1000000-1000499",
    "chr1-1002000-1002499",
    "chr1-1004000-1004499",
    "chr1-1008000-1008499",
    "chr1-1012000-1012499"
  )
  triplets <- data.frame(
    TF = c("Jun", "Jun", "Jun", "Sox9"),
    region = peak_names[c(1, 2, 3, 5)],
    gene = c("Gene1", "Gene1", "Gene2", "Klf4"),
    score = c(0.8, 0.6, 0.5, 0.9),
    stringsAsFactors = FALSE
  )
  region_gene <- triplets[, c("region", "gene", "score", "TF")]
  region_auc <- auc
  region_auc[] <- region_auc[] * 0.7
  dat$srt[["scenicplus"]] <- Seurat::CreateAssayObject(counts = as.matrix(auc))
  peak_counts <- matrix(
    stats::rpois(length(peak_names) * ncol(dat$srt), lambda = 4),
    nrow = length(peak_names),
    ncol = ncol(dat$srt),
    dimnames = list(peak_names, colnames(dat$srt))
  )
  peak_counts[1:3, dat$srt$CellType == "A"] <- peak_counts[1:3, dat$srt$CellType == "A"] + 8
  chromatin <- Signac::CreateChromatinAssay(
    counts = Matrix::Matrix(peak_counts, sparse = TRUE),
    sep = c("-", "-")
  )
  dat$srt[["peaks"]] <- chromatin
  dat$srt@tools$SCENICPlus <- list(
    scores = scores,
    auc = t(as.matrix(auc)),
    regulon_list = regulon_list,
    tf_gene = tf_gene,
    triplets = triplets,
    region_gene = region_gene,
    auc_regions = region_auc
  )
  dat$srt@tools$SCENIC <- NULL
  list(srt = dat$srt, auc = auc, regulons = regulons)
}

test_that("SCENIC row scaling matches the legacy apply standard deviations", {
  mat <- rbind(
    c(1, 2, 3, 4),
    c(4, NA, 6, 8),
    c(5, 5, 5, 5),
    c(NA, NA, NA, NA)
  )
  rownames(mat) <- paste0("Regulon", seq_len(nrow(mat)))
  old_scale <- function(x) {
    means <- rowMeans(x, na.rm = TRUE)
    sds <- apply(x, 1, stats::sd, na.rm = TRUE)
    variable <- is.finite(sds) & sds > sqrt(.Machine$double.eps)
    out <- x
    if (any(variable)) {
      out[variable, ] <- sweep(sweep(x[variable, , drop = FALSE], 1, means[variable], "-"), 1, sds[variable], "/")
    }
    if (any(!variable)) out[!variable, ] <- 0
    attr(out, "constant_rows") <- rownames(x)[!variable]
    out
  }

  expect_equal(scenic_scale_rows(mat), old_scale(mat), tolerance = 1e-12)
})

test_that("SCENIC group heatmap order matches legacy row-wise maxima", {
  mat <- rbind(
    RegulonA = c(A = 0.8, B = 0.8, C = 0.2),
    RegulonB = c(A = NA_real_, B = 0.4, C = 0.9),
    RegulonC = c(A = NA_real_, B = NA_real_, C = NA_real_)
  )
  features <- rownames(mat)
  legacy <- function(x, features) {
    x <- as.matrix(x[features, , drop = FALSE])
    x[!is.finite(x)] <- NA_real_
    max_group_idx <- apply(x, 1L, function(values) {
      if (all(is.na(values))) {
        return(NA_integer_)
      }
      which.max(replace(values, is.na(values), -Inf))
    })
    max_value <- apply(x, 1L, function(values) {
      out <- suppressWarnings(max(values, na.rm = TRUE))
      if (is.finite(out)) out else NA_real_
    })
    max_group <- colnames(x)[max_group_idx]
    features[order(
      factor(max_group, levels = colnames(x)),
      -max_value,
      features,
      na.last = TRUE
    )]
  }

  expect_identical(
    scenic_order_heatmap_features(mat, features, heatmap_order = "group"),
    legacy(mat, features)
  )
})

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
  dat$auc[, ] <- 0
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

test_that("SCENICPlot heatmap_dotplot colors dots by TF expression", {
  dat <- make_scenic_plot_mock(seed = 7)
  out <- SCENICPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "heatmap_dotplot",
    features = rownames(dat$auc)[1:4],
    verbose = FALSE
  )

  expect_s3_class(out$plot, "ggplot")
  expect_true(all(c("RSS", "TF_expr", "TF") %in% colnames(out$plot_data)))
  expect_false(all(is.na(out$plot_data[["TF_expr"]])))
  expect_setequal(
    as.character(unique(out$plot_data[["TF"]])),
    c("Jun", "Atf3", "Fos", "Klf4")
  )
})

test_that("SCENICPlot activity_dim can compare TF expression", {
  dat <- make_scenic_plot_mock(seed = 8, add_embedding = TRUE)
  out <- SCENICPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "activity_dim",
    features = "Jun(+)",
    compare_expression = TRUE,
    reduction = "umap",
    return_data = TRUE,
    verbose = FALSE
  )

  expect_true(inherits(out$plot, c("ggplot", "patchwork")))
  expect_equal(as.character(out$plot_data[["regulon"]]), "Jun(+)")
  expect_equal(as.character(out$plot_data[["TF"]]), "Jun")
})

test_that("SCENICPlusPlot defaults to heatmap_dotplot and SCENIC+ slots", {
  dat <- make_scenicplus_plot_mock()
  out <- SCENICPlusPlot(
    dat$srt,
    group.by = "CellType",
    features = "Jun",
    verbose = FALSE
  )

  expect_equal(out$plot_type, "heatmap_dotplot")
  expect_s3_class(out$plot, "ggplot")
  expect_true("Jun(+)" %in% as.character(out$plot_data[["regulon"]]))
  expect_equal(unique(as.character(out$plot_data[["TF"]])), "Jun")
})

test_that("SCENICPlusPlot eregulon_dim includes region AUC when stored", {
  dat <- make_scenicplus_plot_mock()
  out <- SCENICPlusPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "eregulon_dim",
    features = "Jun",
    reduction = "umap",
    verbose = FALSE
  )

  expect_true(inherits(out$plot, c("ggplot", "patchwork")))
  expect_equal(as.character(out$plot_data[["regulon"]]), "Jun(+)")
  expect_gt(length(out$plots), 1)
})

test_that("SCENIC network palette remaps RdYlBu to Set1", {
  expect_equal(scenic_network_palette("RdYlBu"), "Set1")
  expect_equal(scenic_network_palette("Chinese"), "Chinese")
  shared <- palette_colors(c("Jun", "Atf3", "Fos"), palette = "Set1")
  one <- palette_colors("Atf3", palette = "Set1", palcolor = shared)
  expect_equal(unname(one[["Atf3"]]), unname(shared[["Atf3"]]))
  expect_false(identical(unname(shared[["Jun"]]), unname(shared[["Atf3"]])))
})

test_that("SCENICPlusPlot network uses TF-region-gene hubs from triplets", {
  dat <- make_scenicplus_plot_mock()
  out <- SCENICPlusPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "network",
    features = "Jun",
    verbose = FALSE
  )

  expect_s3_class(out$plot, "ggplot")
  expect_true("Jun" %in% out$plot_data$nodes$name)
  expect_true("region" %in% as.character(out$plot_data$nodes$node_type))
  expect_true("gene" %in% as.character(out$plot_data$nodes$node_type))
  tf_node <- out$plot_data$nodes[out$plot_data$nodes$name == "Jun", ]
  expect_equal(tf_node$x, 0, tolerance = 1e-8)
  expect_equal(tf_node$y, 0, tolerance = 1e-8)
  expect_true(all(c("tf_region", "region_gene") %in% out$plot_data$edges$edge_type))
})

test_that("SCENICPlot network places the TF at the hub center", {
  dat <- make_scenic_plot_mock()
  out <- SCENICPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "network",
    features = "Jun",
    verbose = FALSE
  )

  tf_node <- out$plot_data$nodes[out$plot_data$nodes$name == "Jun", ]
  gene_nodes <- out$plot_data$nodes[out$plot_data$nodes$node_type == "gene", ]
  expect_equal(tf_node$x, 0, tolerance = 1e-8)
  expect_equal(tf_node$y, 0, tolerance = 1e-8)
  expect_true(all(gene_nodes$x^2 + gene_nodes$y^2 > 1))
  expect_true(all(out$plot_data$edges$from == "Jun"))
})

test_that("SCENICPlusPlot egrn draws TF, region, and gene nodes", {
  dat <- make_scenicplus_plot_mock()
  out <- SCENICPlusPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "egrn",
    features = "Jun",
    verbose = FALSE
  )

  nodes <- out$plot_data$nodes
  expect_true(all(c("TF", "region", "gene") %in% as.character(nodes$node_type)))
  expect_true("Jun" %in% nodes$name)
})

test_that("SCENIC+ region diamonds are larger than the previous pin-head default", {
  dat <- make_scenicplus_plot_mock()
  out <- SCENICPlusPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "egrn",
    features = "Jun",
    verbose = FALSE
  )
  styled <- scenic_network_style_data(
    node_data = out$plot_data$nodes,
    edge_plot = out$plot_data$edges,
    edge_data = out$plot_data$edges
  )
  region_size <- unique(styled$nodes$node_size[as.character(styled$nodes$node_type) == "region"])
  gene_size <- unique(styled$nodes$node_size[as.character(styled$nodes$node_type) == "gene"])
  expect_gte(region_size, 3)
  expect_lt(region_size, gene_size)
})

test_that("SCENICPlusPlot egrn can use a tripartite layout", {
  dat <- make_scenicplus_plot_mock()
  out <- SCENICPlusPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "egrn",
    features = "Jun",
    network_layout = "tripartite",
    verbose = FALSE
  )

  nodes <- out$plot_data$nodes
  tf_x <- unique(nodes$x[as.character(nodes$node_type) == "TF"])
  region_x <- unique(nodes$x[as.character(nodes$node_type) == "region"])
  gene_x <- unique(nodes$x[as.character(nodes$node_type) == "gene"])
  expect_lt(tf_x, region_x)
  expect_lt(region_x, gene_x)
})

test_that("SCENICPlot overlap returns pairwise Jaccard values", {
  dat <- make_scenic_plot_mock()
  out <- SCENICPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "overlap",
    features = c("Jun(+)", "Atf3(+)", "Fos(+)"),
    verbose = FALSE
  )

  expect_s3_class(out$plot, "ggplot")
  expect_true(all(c("regulon_1", "regulon_2", "jaccard") %in% colnames(out$plot_data)))
  diag <- out$plot_data$regulon_1 == out$plot_data$regulon_2
  expect_true(all(out$plot_data$jaccard[diag] == 1))
})

test_that("SCENIC network_graph hub layout keeps TFs off a single line", {
  dat <- make_scenic_plot_mock()
  out <- SCENICPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "network_graph",
    network_layout = "hub",
    max_targets = 5,
    verbose = FALSE
  )

  tf_nodes <- out$plot_data$nodes[as.character(out$plot_data$nodes$node_type) == "TF", ]
  expect_gt(nrow(tf_nodes), 2)
  expect_gt(sd(tf_nodes$x), 0.2)
  expect_gt(sd(tf_nodes$y), 0.2)
})

test_that("SCENICPlusPlot network can disable region nodes", {
  dat <- make_scenicplus_plot_mock()
  out <- SCENICPlusPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "network",
    features = "Jun",
    network_include_regions = FALSE,
    verbose = FALSE
  )

  expect_false("region" %in% as.character(out$plot_data$nodes$node_type))
  expect_true(all(out$plot_data$edges$from == "Jun"))
})

test_that("scenic_tf_from_regulon parses SCENIC and SCENIC+ names", {
  expect_equal(scenic_tf_from_regulon("Jun(+)"), "Jun")
  expect_equal(scenic_tf_from_regulon("Jun(-)"), "Jun")
  expect_equal(scenic_tf_from_regulon("Sox9_direct_+/+_(12g)"), "Sox9")
  expect_equal(scenic_tf_from_regulon("SOX9_+_+"), "SOX9")
})

test_that("scenic_match_regulon_features matches SCENIC+ eRegulon names", {
  available <- c("Jun_direct_+/+_(3g)", "Sox9_direct_+/+_(2g)")
  expect_equal(
    scenic_match_regulon_features("Jun", available),
    "Jun_direct_+/+_(3g)"
  )
})

test_that("scenic_parse_genomic_region accepts Signac and SCENIC+ coordinates", {
  parsed <- scenic_parse_genomic_region(c(
    "chr1-100-200",
    "chr1:300-400",
    "chr2_500_600",
    "chr1-1e+06-1000499"
  ))
  expect_equal(parsed$seqnames, c("chr1", "chr1", "chr2", "chr1"))
  expect_equal(parsed$start, c(100, 300, 500, 1e6))
  expect_equal(parsed$end, c(200, 400, 600, 1000499))
})

test_that("SCENICPlusPlot coverage returns region-gene tracks", {
  skip_if_not_installed("Signac")
  dat <- make_scenicplus_plot_mock()
  out <- SCENICPlusPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "coverage",
    features = "Gene1",
    extend.upstream = 1000,
    extend.downstream = 1000,
    verbose = FALSE
  )

  expect_true(inherits(out$plot, c("ggplot", "patchwork")))
  expect_true(all(out$plot_data$gene == "Gene1"))
  expect_true(all(c("region", "start", "end", "score") %in% colnames(out$plot_data)))
  expect_gt(nrow(out$plot_data), 0)
})

test_that("SCENICPlusPlot coverage resolves TF features to target genes", {
  skip_if_not_installed("Signac")
  dat <- make_scenicplus_plot_mock()
  out <- SCENICPlusPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "coverage",
    features = "Jun",
    extend.upstream = 1000,
    extend.downstream = 1000,
    verbose = FALSE
  )

  expect_true(all(out$plot_data$TF == "Jun"))
  expect_true("Gene1" %in% out$plot_data$gene)
})
