test_that("RunGSEA and RunEnrichment tolerate a database with an empty version", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::local_mocked_bindings(
    PrepareDB = function(species, db, db_update = FALSE, db_version = "latest",
                         db_IDtypes = "symbol", convert_species = TRUE,
                         Ensembl_version = NULL, mirror = NULL, verbose = TRUE, ...) {
      list(
        "Homo_sapiens" = list(
          KEGG = list(
            TERM2GENE = data.frame(
              Term = rep(c("Term1", "Term2"), each = 3),
              symbol = paste0("G", 1:6)
            ),
            TERM2NAME = data.frame(
              Term = c("Term1", "Term2"),
              Name = c("Term One", "Term Two")
            ),
            version = character(0)
          )
        )
      )
    },
    .package = "scop"
  )

  genes <- paste0("G", 1:6)
  groups <- rep("Cluster1", length(genes))

  gsea <- RunGSEA(
    geneID = genes,
    geneScore = 6:1,
    geneID_groups = groups,
    species = "Homo_sapiens",
    db = "KEGG",
    IDtype = "symbol",
    result_IDtype = "symbol",
    minGSSize = 2,
    verbose = FALSE
  )
  expect_true(nrow(gsea$enrichment) > 0)
  expect_identical(unique(gsea$enrichment$Version), "")

  enrich <- RunEnrichment(
    geneID = genes,
    geneID_groups = groups,
    species = "Homo_sapiens",
    db = "KEGG",
    IDtype = "symbol",
    result_IDtype = "symbol",
    backend = "r",
    minGSSize = 2,
    verbose = FALSE
  )
  expect_true(nrow(enrich$enrichment) > 0)
  expect_identical(unique(enrich$enrichment$Version), "")
})
