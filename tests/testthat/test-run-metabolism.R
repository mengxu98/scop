test_that("metabolism GMT sources do not depend on an installed scMetabolism package", {
  sources <- scmetabolism_gmt_sources()
  expect_named(sources, c("KEGG", "Reactome"))
  expect_true(all(startsWith(sources, "https://raw.githubusercontent.com/mengxu98/datasets/")))
})

test_that("RunMetabolism uses the current gene conversion on each call", {
  counts <- Matrix::Matrix(matrix(
    c(5, 3, 1, 2, 2, 1, 4, 3),
    nrow = 4,
    dimnames = list(c("Gapdh", "MouseA", "MouseB", "Other"), c("cell1", "cell2"))
  ), sparse = TRUE)
  object <- SeuratObject::CreateSeuratObject(counts)
  calls <- list()
  local_mocked_bindings(
    load_scmetabolism_gmt = function(...) {
      list(
        gene_sets = list(KEGG.TEST = c("GAPDH", "HUMAN_ONLY")),
        term_info = data.frame(Name = "Test pathway", row.names = "KEGG.TEST")
      )
    },
    GeneConvert = function(...) {
      calls[[length(calls) + 1L]] <<- list(...)
      list(geneID_expand = data.frame(
        from_geneID = "HUMAN_ONLY",
        symbol = c("MouseA", "MouseB")[[length(calls)]]
      ))
    },
    run_aucell_scores = function(expr_counts, gene_sets, ...) {
      matrix(
        0,
        nrow = ncol(expr_counts), ncol = length(gene_sets),
        dimnames = list(colnames(expr_counts), names(gene_sets))
      )
    },
    .package = "scop"
  )
  run <- function() {
    RunMetabolism(
      object,
      assay = "RNA", db = "KEGG", species = "Mus_musculus",
      use_preparedb = FALSE, minGSSize = 1, new_assay = FALSE,
      Ensembl_version = 116, biomart = "ensembl", mirror = "useast",
      max_tries = 2, verbose = FALSE
    )
  }
  first <- run()
  second <- run()

  expect_length(calls, 2L)
  expect_identical(calls[[1L]], list(
    geneID = c("GAPDH", "HUMAN_ONLY"),
    geneID_from_IDtype = "symbol", geneID_to_IDtype = "symbol",
    species_from = "Homo_sapiens", species_to = "Mus_musculus",
    Ensembl_version = 116, biomart = "ensembl", mirror = "useast",
    max_tries = 2, verbose = FALSE
  ))
  expect_setequal(first@tools$Metabolism_AUCell$gene_sets[[1L]], c("Gapdh", "MouseA"))
  expect_setequal(second@tools$Metabolism_AUCell$gene_sets[[1L]], c("Gapdh", "MouseB"))
})
