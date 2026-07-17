test_that("RunMetabolism caches repeated BioMart conversion inputs", {
  cache <- scop:::.scmetabolism_conversion_cache
  rm(list = ls(envir = cache, all.names = TRUE), envir = cache)
  calls <- 0L
  expected <- list(
    geneID_expand = data.frame(
      from_geneID = c("GAPDH", "ACTB"),
      symbol = c("Gapdh", "Actb")
    )
  )
  local_mocked_bindings(
    GeneConvert = function(...) {
      calls <<- calls + 1L
      expected
    },
    .package = "scop"
  )

  first <- scop:::scmetabolism_convert_genes(
    all_human_genes = c("GAPDH", "ACTB"),
    species = "Mus_musculus",
    Ensembl_version = 116,
    verbose = FALSE
  )
  second <- scop:::scmetabolism_convert_genes(
    all_human_genes = c("GAPDH", "ACTB"),
    species = "Mus_musculus",
    Ensembl_version = 116,
    verbose = FALSE
  )

  expect_equal(calls, 1L)
  expect_identical(first, expected)
  expect_identical(second, expected)
})
