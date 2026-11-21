test_that("CPP ORA respects an explicit universe and reports its effective size", {
  term2gene <- data.frame(
    Term = rep(c("A", "B", "C"), each = 4L),
    symbol = paste0("G", seq_len(12L))
  )
  term2name <- data.frame(
    Term = c("A", "B", "C"),
    Name = c("Term A", "Term B", "Term C")
  )

  default <- run_ora_result(
    gene = c("G1", "G2", "G5"),
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    IDtype = "symbol",
    min_gs_size = 2L
  )
  restricted <- run_ora_result(
    gene = c("G1", "G2", "G5"),
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    IDtype = "symbol",
    min_gs_size = 2L,
    universe = c(paste0("G", seq_len(10L)), "unannotated")
  )

  expect_true(all(default$UniverseSize == 12L))
  expect_true(all(restricted$UniverseSize == 10L))
  expect_equal(restricted$GeneRatio[restricted$ID == "A"], "2/3")
  expect_equal(restricted$BgRatio[restricted$ID == "A"], "4/10")
})

test_that("CPP ORA supports comparison-specific universes", {
  term2gene <- data.frame(
    Term = rep(c("A", "B", "C"), each = 4L),
    symbol = paste0("G", seq_len(12L))
  )
  term2name <- data.frame(
    Term = c("A", "B", "C"),
    Name = c("Term A", "Term B", "Term C")
  )
  genes <- list(
    Group1 = c("G1", "G2"),
    Group2 = c("G3", "G5")
  )
  universe <- list(
    Group1 = paste0("G", seq_len(10L)),
    Group2 = paste0("G", 3:12)
  )

  result <- run_ora_result(
    gene = genes,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    IDtype = "symbol",
    min_gs_size = 2L,
    universe = universe
  )

  expect_named(result, names(genes))
  expect_true(all(result$Group1$UniverseSize == 10L))
  expect_true(all(result$Group2$UniverseSize == 10L))
  expect_equal(result$Group1$BgRatio[result$Group1$ID == "A"], "4/10")
  expect_equal(result$Group2$BgRatio[result$Group2$ID == "A"], "2/10")
})

test_that("RunEnrichment uses the explicit universe for custom gene sets", {
  term2gene <- data.frame(
    Term = rep(c("A", "B", "C"), each = 4L),
    symbol = paste0("G", seq_len(12L))
  )
  term2name <- data.frame(
    Term = c("A", "B", "C"),
    Name = c("Term A", "Term B", "Term C")
  )

  result <- RunEnrichment(
    geneID = c("G1", "G2", "G5"),
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    universe = paste0("G", seq_len(10L)),
    minGSSize = 2L,
    verbose = FALSE
  )

  expect_true(nrow(result$enrichment) > 0L)
  expect_true(all(result$enrichment$UniverseSize == 10L))
  expect_equal(result$enrichment$BgRatio[result$enrichment$ID == "A"], "4/10")
})

test_that("RunEnrichment derives the background before applying DE_threshold", {
  skip_if_not_installed("SeuratObject")

  counts <- matrix(
    1,
    nrow = 8L,
    ncol = 4L,
    dimnames = list(paste0("G", seq_len(8L)), paste0("Cell", seq_len(4L)))
  )
  srt <- SeuratObject::CreateSeuratObject(counts = counts)
  de <- data.frame(
    gene = paste0("G", seq_len(8L)),
    group1 = rep("Cluster1", 8L),
    group2 = rep("rest", 8L),
    avg_log2FC = c(2, 1.5, rep(0, 6L)),
    p_val_adj = c(0.001, 0.01, rep(0.5, 6L)),
    stringsAsFactors = FALSE
  )
  srt@tools[["DEtest_test"]] <- list(AllMarkers_wilcox = de)

  term2gene <- data.frame(
    Term = rep(c("A", "B"), each = 4L),
    symbol = paste0("G", seq_len(8L))
  )
  term2name <- data.frame(
    Term = c("A", "B"),
    Name = c("Term A", "Term B")
  )

  srt <- RunEnrichment(
    object = srt,
    group.by = "test",
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    DE_threshold = "avg_log2FC > 1 & p_val_adj < 0.05",
    minGSSize = 2L,
    verbose = FALSE
  )
  result <- srt@tools[["Enrichment_test_wilcox"]][["enrichment"]]

  expect_true(nrow(result) > 0L)
  expect_true(all(result$UniverseSize == 8L))
  expect_true(all(result$BgRatio == "4/8"))
})

test_that("R and CPP backends use the same explicit universe", {
  skip_if_not_installed("clusterProfiler")

  term2gene <- data.frame(
    Term = rep(c("A", "B", "C"), each = 4L),
    symbol = paste0("G", seq_len(12L))
  )
  term2name <- data.frame(
    Term = c("A", "B", "C"),
    Name = c("Term A", "Term B", "Term C")
  )
  args <- list(
    geneID = c("G1", "G2", "G5", "G6"),
    geneID_groups = c("Group1", "Group1", "Group2", "Group2"),
    universe = list(
      Group1 = paste0("G", seq_len(10L)),
      Group2 = paste0("G", 3:12)
    ),
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    minGSSize = 2L,
    verbose = FALSE
  )

  cpp <- do.call(RunEnrichment, c(args, list(backend = "cpp")))
  r <- do.call(RunEnrichment, c(args, list(backend = "r")))

  select_columns <- function(x) {
    out <- x$enrichment
    out <- out[order(out$Groups, out$ID), , drop = FALSE]
    out[, c(
      "Groups", "ID", "GeneRatio", "BgRatio", "Count",
      "UniverseSize", "pvalue", "p.adjust"
    ), drop = FALSE]
  }

  expect_equal(select_columns(cpp), select_columns(r), tolerance = 1e-10)
})
