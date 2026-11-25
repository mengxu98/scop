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

test_that("excluding a query gene leaves the ORA background unchanged", {
  term2gene <- data.frame(Term = rep(LETTERS[1:3], each = 4L), symbol = paste0("G", 1:12))
  args <- list(
    geneID = c("G1", "G2", "G5"),
    TERM2GENE = term2gene,
    universe = term2gene$symbol,
    minGSSize = 2L,
    verbose = FALSE
  )
  for (backend in c("cpp", if (requireNamespace("clusterProfiler", quietly = TRUE)) "r")) {
    full <- do.call(RunEnrichment, c(args, list(backend = backend)))$enrichment
    excluded <- do.call(RunEnrichment, c(args, list(backend = backend, geneID_exclude = "G2")))$enrichment
    full <- full[full$ID == "A", , drop = FALSE]
    excluded <- excluded[excluded$ID == "A", , drop = FALSE]
    expect_equal(full$GeneRatio, "2/3")
    expect_equal(excluded$GeneRatio, "1/2")
    expect_equal(excluded$BgRatio, full$BgRatio)
    expect_equal(excluded$UniverseSize, 12L)
  }
})

test_that("real DE marker storage is not the earlier ORA eligibility set", {
  skip_if_not_installed("Seurat")

  genes <- paste0("G", 1:12)
  counts <- matrix(8L, nrow = 12L, ncol = 48L, dimnames = list(genes, paste0("C", 1:48)))
  counts[1:3, 1:24] <- counts[4:6, 25:48] <- 18L
  counts[1:3, 25:48] <- counts[4:6, 1:24] <- 1L
  counts[12, ] <- 0L
  counts[12, c(1, 25)] <- 1L
  srt <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE))
  srt$condition <- rep(c("Left", "Right"), each = 24)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)

  markers <- RunDEtest(srt, group.by = "condition", verbose = FALSE)
  eligible <- RunDEtest(srt,
    group.by = "condition", only.pos = FALSE, fc.threshold = 1, verbose = FALSE
  )
  de_by_group <- function(x) {
    tab <- x@tools[["DEtest_condition"]][["AllMarkers_wilcox"]]
    split(tab$gene, tab$group1)
  }
  background <- de_by_group(eligible)
  expect_setequal(de_by_group(markers)$Left, genes[1:3])
  expect_setequal(de_by_group(markers)$Right, genes[4:6])
  expect_setequal(background$Left, genes[1:11])
  expect_setequal(background$Right, genes[1:11])

  args <- list(
    object = markers, group.by = "condition",
    TERM2GENE = data.frame(Term = rep(LETTERS[1:3], each = 4L), symbol = genes),
    minGSSize = 2L, verbose = FALSE
  )
  automatic <- do.call(RunEnrichment, args)@tools[["Enrichment_condition_wilcox"]]
  explicit <- do.call(RunEnrichment, c(args, list(universe = background)))@tools[["Enrichment_condition_wilcox"]]
  expect_identical(automatic$input, explicit$input)
  expect_setequal(automatic$enrichment$UniverseSize, 3L)
  expect_setequal(explicit$enrichment$UniverseSize, 11L)
})

test_that("Seurat ORA reruns refresh only the requested comparison and database", {
  skip_if_not_installed("SeuratObject")

  genes <- paste0("G", 1:12)
  counts <- matrix(1L, nrow = 12L, ncol = 4L, dimnames = list(genes, paste0("C", 1:4)))
  srt <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE))
  annotations <- list(
    TERM2GENE = data.frame(Term = rep(LETTERS[1:3], each = 4L), symbol = genes),
    TERM2NAME = data.frame(Term = LETTERS[1:3], Name = LETTERS[1:3]),
    version = "fixture"
  )
  testthat::local_mocked_bindings(
    PrepareDB = function(species, db, ...) {
      stats::setNames(list(stats::setNames(rep(list(annotations), length(db)), db)), species)
    },
    .package = "scop"
  )
  enrich <- function(object, groups, db, universe, minGSSize = 2L) {
    RunEnrichment(
      object, group.by = "test", geneID = rep(genes[c(1, 2, 5)], length(groups)),
      geneID_groups = rep(groups, each = 3L), db = db, universe = universe,
      minGSSize = minGSSize, verbose = FALSE
    )
  }
  key <- "Enrichment_test_wilcox"
  srt <- enrich(srt, c("Alpha", "Beta"), c("GO_BP", "OtherDB"), genes)

  # Add a previously cached simplified view to check invalidation without
  # changing the actual C++ calculations or requiring the optional R backend.
  result <- srt@tools[[key]]
  sims <- lapply(c("Alpha-GO_BP", "Beta-GO_BP"), function(nm) {
    tab <- result$results[[nm]]
    tab$Database <- "GO_BP_sim"
    tab
  })
  names(sims) <- c("Alpha-GO_BP_sim", "Beta-GO_BP_sim")
  result$results <- c(result$results, sims)
  result$enrichment <- do.call(rbind, unname(result$results))
  srt@tools[[key]] <- result

  term <- function(result) {
    tab <- result$enrichment
    tab[tab$Groups == "Alpha" & tab$Database == "GO_BP" & tab$ID == "A", , drop = FALSE]
  }
  expect_equal(term(result)$BgRatio, "4/12")
  expect_equal(term(result)$pvalue, 52 / 220)

  rerun <- enrich(srt, "Alpha", "GO_BP", genes[1:10])
  refreshed <- rerun@tools[[key]]
  expect_equal(term(refreshed)$BgRatio, "4/10")
  expect_equal(term(refreshed)$pvalue, 40 / 120)
  expect_equal(sum(names(refreshed$results) == "Alpha-GO_BP"), 1L)
  expect_false("Alpha-GO_BP_sim" %in% names(refreshed$results))
  expect_false(any(refreshed$enrichment$Groups == "Alpha" & refreshed$enrichment$Database == "GO_BP_sim"))
  untouched <- c("Beta-GO_BP", "Beta-GO_BP_sim", "Alpha-OtherDB", "Beta-OtherDB")
  expect_identical(refreshed$results[untouched], result$results[untouched])

  empty <- enrich(rerun, "Alpha", "GO_BP", genes[1:10], minGSSize = 5L)@tools[[key]]
  expect_equal(nrow(term(empty)), 0L)
  expect_false("Alpha-GO_BP" %in% names(empty$results))
  expect_identical(empty$results[untouched], result$results[untouched])
})
