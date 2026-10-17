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
