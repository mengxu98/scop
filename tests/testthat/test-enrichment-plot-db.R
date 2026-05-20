test_that("plot database aliases resolve to available child databases", {
  enrichment <- data.frame(
    Database = c("GO_BP", "KEGG"),
    stringsAsFactors = FALSE
  )

  expect_equal(scop:::resolve_enrichment_plot_db("GO", enrichment), "GO_BP")
  expect_equal(
    scop:::resolve_enrichment_plot_db(c("KEGG", "GO"), enrichment),
    c("KEGG", "GO_BP")
  )
})

test_that("simplified database aliases resolve to available child databases", {
  enrichment <- data.frame(
    Database = c("GO_CC_sim", "Reactome"),
    stringsAsFactors = FALSE
  )

  expect_equal(scop:::resolve_enrichment_plot_db("GO_sim", enrichment), "GO_CC_sim")
})

test_that("plot database alias resolution is reusable for non-GO databases", {
  enrichment <- data.frame(
    Database = c("MSigDB_H", "MSigDB_C2"),
    stringsAsFactors = FALSE
  )
  aliases <- list(MSigDB = c("MSigDB_H", "MSigDB_C2", "MSigDB_C5"))

  expect_equal(
    scop:::resolve_enrichment_plot_db("MSigDB", enrichment, aliases = aliases),
    c("MSigDB_H", "MSigDB_C2")
  )
})

test_that("GSEA result names can be matched after GO alias resolution", {
  enrichment <- data.frame(
    Groups = "g1",
    Database = "GO_BP",
    stringsAsFactors = FALSE
  )
  results <- list("g1-GO_BP" = TRUE)

  db <- scop:::resolve_enrichment_plot_db("GO", enrichment)
  comb <- expand.grid(unique(enrichment[["Groups"]]), db)
  use <- names(results)[names(results) %in% paste(comb$Var1, comb$Var2, sep = "-")]

  expect_equal(use, "g1-GO_BP")
})

test_that("group filtering reports missing selected enrichment groups before db errors", {
  enrichment <- data.frame(
    Groups = "Acinar",
    Database = "GO_BP",
    stringsAsFactors = FALSE
  )

  expect_equal(scop:::resolve_enrichment_plot_db("GO_BP", enrichment), "GO_BP")
  expect_error(
    scop:::filter_enrichment_plot_groups(
      enrichment,
      group_use = c("Ductal", "Endocrine")
    ),
    "No enrichment result found for selected groups"
  )
})
