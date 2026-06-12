test_that("species names accept common Latin-name separators", {
  expect_identical(
    normalize_species_name(c("Bos taurus", "Bos_taurus", "bos.taurus")),
    rep("Bos_taurus", 3)
  )
  expect_identical(normalize_species_name("  Mus musculus  "), "Mus_musculus")
  expect_identical(normalize_species_name("Homo-sapiens"), "Homo_sapiens")
})

test_that("normalized species names produce Bioconductor org package names", {
  species <- normalize_species_name("Bos taurus")
  sp <- unlist(strsplit(species, split = "_"))
  org_sp <- paste0(
    "org.",
    paste0(substring(sp, 1, 1), collapse = ""),
    ".eg.db"
  )
  expect_identical(org_sp, "org.Bt.eg.db")
})
