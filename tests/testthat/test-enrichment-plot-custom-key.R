test_that("EnrichmentPlot finds results stored without a grouping variable", {
  skip_if_not_installed("Seurat")

  data("pancreas_sub", package = "scop")
  cell_types <- as.character(pancreas_sub$CellType)
  cells1 <- head(colnames(pancreas_sub)[cell_types == "Ductal"], 80L)
  cells2 <- head(colnames(pancreas_sub)[cell_types == "Endocrine"], 80L)

  srt <- RunDEtest(
    pancreas_sub,
    cells1 = cells1,
    cells2 = cells2,
    only.pos = FALSE,
    verbose = FALSE
  )
  expect_true("DEtest_custom" %in% names(srt@tools))

  de_df <- srt@tools[["DEtest_custom"]][["AllMarkers_wilcox"]]
  threshold <- "p_val_adj < 0.05 & avg_log2FC > 0.5"
  query_genes <- unique(as.character(
    de_df[["gene"]][
      de_df[["p_val_adj"]] < 0.05 & de_df[["avg_log2FC"]] > 0.5
    ]
  ))
  expect_gt(length(query_genes), 100L)
  term2gene <- rbind(
    data.frame(
      Term = rep(paste0("custom", seq_len(10)), each = 10),
      symbol = query_genes[seq_len(100)]
    ),
    data.frame(
      Term = "background",
      symbol = setdiff(unique(as.character(de_df[["gene"]])), query_genes)
    )
  )

  srt <- RunEnrichment(
    srt,
    db = "custom",
    TERM2GENE = term2gene,
    DE_threshold = threshold,
    minGSSize = 10,
    verbose = FALSE
  )
  expect_true("Enrichment_custom_wilcox" %in% names(srt@tools))
  enrichment <- srt@tools[["Enrichment_custom_wilcox"]][["enrichment"]]
  expect_gt(nrow(enrichment), 0L)
  expect_true(any(enrichment[["p.adjust"]] < 0.05))

  plot <- EnrichmentPlot(
    srt,
    db = "custom",
    plot_type = "bar",
    topTerm = 10,
    verbose = FALSE
  )
  expect_s3_class(plot, "ggplot")
  expect_setequal(plot$data[["ID"]], paste0("custom", seq_len(10)))

  expect_error(
    EnrichmentPlot(
      srt,
      db = "custom",
      group.by = "CellType",
      plot_type = "bar",
      verbose = FALSE
    ),
    "Enrichment_custom_wilcox"
  )
})

test_that("GSEAPlot finds results stored without a grouping variable", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("clusterProfiler")

  data("pancreas_sub", package = "scop")
  cell_types <- as.character(pancreas_sub$CellType)
  cells1 <- head(colnames(pancreas_sub)[cell_types == "Ductal"], 80L)
  cells2 <- head(colnames(pancreas_sub)[cell_types == "Endocrine"], 80L)

  srt <- RunDEtest(
    pancreas_sub,
    cells1 = cells1,
    cells2 = cells2,
    only.pos = FALSE,
    verbose = FALSE
  )
  de_df <- srt@tools[["DEtest_custom"]][["AllMarkers_wilcox"]]
  genes <- unique(as.character(
    de_df[["gene"]][
      de_df[["p_val_adj"]] < 0.05 & de_df[["avg_log2FC"]] > 0.5
    ]
  ))
  scores <- stats::setNames(
    de_df[["avg_log2FC"]][match(genes, de_df[["gene"]])],
    genes
  )
  term2gene <- rbind(
    data.frame(
      Term = rep(paste0("custom", seq_len(10)), each = 10),
      symbol = genes[seq_len(100)]
    ),
    data.frame(
      Term = "background",
      symbol = setdiff(unique(as.character(de_df[["gene"]])), genes)
    )
  )

  srt <- RunGSEA(
    srt,
    geneID = genes,
    geneScore = scores,
    geneID_groups = rep("group1", length(genes)),
    scoreType = "pos",
    db = "custom",
    TERM2GENE = term2gene,
    minGSSize = 10,
    verbose = FALSE
  )
  expect_true("GSEA_custom_wilcox" %in% names(srt@tools))

  plot <- GSEAPlot(
    srt,
    db = "custom",
    plot_type = "bar",
    verbose = FALSE
  )
  expect_s3_class(plot, "ggplot")
  expect_true(all(
    as.character(plot$data[["ID"]]) %in%
      as.character(srt@tools[["GSEA_custom_wilcox"]][["enrichment"]][["ID"]])
  ))

  expect_error(
    GSEAPlot(
      srt,
      db = "custom",
      group.by = "CellType",
      plot_type = "bar",
      verbose = FALSE
    ),
    "GSEA_custom_wilcox"
  )
})

test_that("EnrichmentPlot and GSEAPlot plot results passed through `res`", {
  skip_if_not_installed("clusterProfiler")
  term2gene <- data.frame(
    Term = c(rep("Endocrine markers", 5), rep("Ductal markers", 5)),
    symbol = c(
      "INS", "GCG", "SST", "IAPP", "PCSK1",
      "KRT19", "SOX9", "MUC1", "CFTR", "KRT7"
    )
  )
  genes <- unique(term2gene$symbol)

  enrich_out <- RunEnrichment(
    geneID = c(
      "INS", "GCG", "SST", "IAPP", "PRSS1", "CPA1",
      "KRT19", "SOX9", "MUC1", "CFTR", "KRT7", "REG1A"
    ),
    geneID_groups = rep(c("Cluster1", "Cluster2"), each = 6),
    TERM2GENE = term2gene,
    minGSSize = 2,
    verbose = FALSE
  )
  expect_true(any(enrich_out[["enrichment"]][["p.adjust"]] < 0.05))

  plot <- EnrichmentPlot(
    res = enrich_out,
    db = "custom",
    plot_type = "comparison",
    verbose = FALSE
  )
  expect_s3_class(plot, "ggplot")
  expect_setequal(plot$data[["Groups"]], c("Cluster1", "Cluster2"))

  scores <- stats::setNames(seq_along(genes) - 6, genes)
  gsea_out <- RunGSEA(
    geneID = genes,
    geneScore = scores,
    geneID_groups = rep("Cluster1", length(genes)),
    TERM2GENE = term2gene,
    minGSSize = 2,
    verbose = FALSE
  )
  gsea_plot <- GSEAPlot(
    res = gsea_out,
    db = "custom",
    plot_type = "comparison",
    verbose = FALSE
  )
  expect_s3_class(gsea_plot, "ggplot")

  expect_error(EnrichmentPlot(db = "custom", verbose = FALSE), "Either")
  expect_error(GSEAPlot(db = "custom", verbose = FALSE), "Either")
})
