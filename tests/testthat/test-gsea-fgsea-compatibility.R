test_that("GSEA retains fgsea statistics and complete gene sets across dispatcher versions", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("fgsea")
  if (!"by" %in% names(formals(clusterProfiler::GSEA))) {
    skip_if_not_installed("enrichit")
  }
  scores <- stats::setNames(seq(3, -3, length.out = 40), paste0("g", seq_len(40)))
  sets <- list(up = names(scores)[1:7], down = names(scores)[34:40], middle = names(scores)[17:24])
  terms <- data.frame(term = rep(names(sets), lengths(sets)), gene = unlist(sets, use.names = FALSE))
  labels <- data.frame(term = names(sets), name = paste("Pathway", names(sets)))
  for (type in c("std", "pos", "neg")) {
    native <- withr::with_seed(19, fgsea::fgsea(
      sets, scores,
      minSize = 2, maxSize = 30, nPermSimple = 1000,
      eps = 0, scoreType = type, gseaParam = 1, nproc = 0
    ))
    result <- withr::with_seed(19, gsea_fgsea_result(
      scores, 2, 30, 1000, 0, type, "BH", Inf, terms, labels, FALSE
    ))
    expect_s4_class(result, "gseaResult")
    index <- match(result@result$ID, native$pathway)
    expect_equal(result@result$enrichmentScore, native$ES[index])
    expect_equal(result@result$NES, native$NES[index])
    expect_equal(result@result$pvalue, native$pval[index])
    expect_equal(result@result$p.adjust, native$padj[index])
    expect_equal(result@geneSets[names(sets)], sets, ignore_attr = TRUE)
    expect_identical(result@geneList, scores)
    expect_identical(result@result$Description, labels$name[match(result@result$ID, labels$term)])
    curve <- gsInfo(result, result@result$ID[[1]])
    expect_equal(nrow(curve), length(scores))
    expect_true(all(is.finite(curve$runningScore)))
  }
  empty <- withr::with_seed(19, gsea_fgsea_result(
    scores, 30, 35, 1000, 0, "std", "BH", Inf, terms, labels, FALSE
  ))
  expect_true(is.null(empty) || nrow(empty@result) == 0L)
})

test_that("one-to-many output IDs preserve original GSEA ranks", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("fgsea")
  if (!"by" %in% names(formals(clusterProfiler::GSEA))) {
    skip_if_not_installed("enrichit")
  }
  genes <- paste0("g", seq_len(40))
  genes[39:40] <- c("prefix", "prefix2")
  mapping <- data.frame(from_geneID = genes, row.names = genes)
  mapping$ensembl_id <- I(lapply(seq_along(genes), function(i) {
    paste0(genes[[i]], "_", seq_len(if (i == 1L) 3L else 1L))
  }))
  mapping <- mapping[-39, , drop = FALSE]
  testthat::local_mocked_bindings(
    GeneConvert = function(...) list(geneID_collapse = mapping),
    .package = "scop"
  )
  args <- list(
    geneID = genes, geneScore = seq(3, -3, length.out = 40),
    TERM2GENE = data.frame(
      Term = rep(c("up", "down", "middle"), c(7, 7, 8)),
      symbol = c(genes[1:7], genes[34:40], genes[17:24])
    ),
    minGSSize = 2, maxGSSize = 30, verbose = FALSE
  )
  set.seed(19)
  expected <- do.call(RunGSEA, args)
  expected_rng <- .Random.seed
  set.seed(19)
  actual <- do.call(RunGSEA, c(args, list(result_IDtype = "ensembl_id")))
  expect_identical(.Random.seed, expected_rng)
  result <- actual$results[[1]]
  control <- expected$results[[1]]
  expect_identical(result@geneList, control@geneList)
  expect_identical(result@geneSets, control@geneSets)
  columns <- setdiff(names(result@result), "core_enrichment")
  expect_identical(result@result[, columns], control@result[, columns])
  expect_length(result@gene2Symbol, length(genes))
  expect_identical(actual$geneMap$ensembl_id, mapping$ensembl_id)
  expect_equal(nrow(actual$input), length(genes) + 2L)
  expect_identical(actual$input$ensembl_id[actual$input$geneID == "prefix"], "prefix")
  expect_identical(actual$input$ensembl_id[actual$input$geneID == genes[[1]]], mapping$ensembl_id[[1]])
})
