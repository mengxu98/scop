test_that("legacy MSigDB ID columns are normalized to symbol", {
  term2gene <- data.frame(
    Term = "HALLMARK_TP53_PATHWAY",
    check.names = FALSE
  )
  term2gene[["symbol.ensembl_id"]] <- "TP53"

  normalized <- scop:::preparedb_normalize_term2gene_id_columns(term2gene)

  expect_equal(colnames(normalized), c("Term", "symbol"))
  expect_equal(normalized[["symbol"]], "TP53")
})

test_that("source ID type resolution returns one existing TERM2GENE column", {
  term2gene <- data.frame(
    Term = "HALLMARK_TP53_PATHWAY",
    symbol = "TP53",
    stringsAsFactors = FALSE
  )
  default_id_types <- list(MSigDB = c("symbol", "ensembl_id"))

  expect_equal(
    scop:::preparedb_source_idtype("MSigDB", term2gene, default_id_types),
    "symbol"
  )
})

test_that("local OrgDb mapper rejects vector source ID types without switch errors", {
  expect_null(
    scop:::preparedb_local_orgdb_id_map(
      geneID = "TP53",
      geneID_from_IDtype = c("symbol", "ensembl_id"),
      geneID_to_IDtype = "symbol",
      org_sp = "org.Hs.eg.db",
      org_key = "ENTREZID",
      verbose = FALSE
    )
  )
})
