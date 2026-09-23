test_that("MSigDB nested collection names keep every prefix", {
  expect_identical(
    msigdb_collection_db_names("M2:CP:BIOCARTA"),
    c("MSigDB_M2", "MSigDB_M2_CP", "MSigDB_M2_CP_BIOCARTA")
  )
  expect_identical(
    msigdb_collection_db_names(c("H", "MH", "C2:CGP")),
    c("MSigDB_H", "MSigDB_MH", "MSigDB_C2", "MSigDB_C2_CGP")
  )
  expect_identical(msigdb_collection_db_names(c(NA, "")), character(0))
})

test_that("colon MSigDB names normalize to underscores", {
  expect_identical(
    normalize_msigdb_db_names(c(
      "MSigDB_M2:CGP",
      "MSigDB_M2:CP:BIOCARTA",
      "MSigDB_M2",
      "MSigDB",
      "GO_BP",
      "MSigDB_.*"
    )),
    c(
      "MSigDB_M2_CGP",
      "MSigDB_M2_CP_BIOCARTA",
      "MSigDB_M2",
      "MSigDB",
      "GO_BP",
      "MSigDB_.*"
    )
  )
  expect_null(normalize_msigdb_db_names(NULL))
})

test_that("MSigDB collection subsets keep parent unions and exclude siblings", {
  term2name <- data.frame(
    Term = c("CGP1", "BIO1", "REACT1", "H1"),
    Name = c("cgp set", "biocarta set", "reactome set", "hallmark"),
    Collection = c("M2:CGP", "M2:CP:BIOCARTA", "M2:CP:REACTOME", "MH"),
    stringsAsFactors = FALSE
  )
  term2gene <- data.frame(
    Term = c("CGP1", "CGP1", "BIO1", "REACT1", "H1"),
    symbol = c("A", "B", "C", "D", "E"),
    stringsAsFactors = FALSE
  )
  subsets <- preparedb_msigdb_collection_subsets(term2gene, term2name)

  expect_setequal(
    names(subsets),
    c(
      "MSigDB_M2",
      "MSigDB_M2_CGP",
      "MSigDB_M2_CP",
      "MSigDB_M2_CP_BIOCARTA",
      "MSigDB_M2_CP_REACTOME",
      "MSigDB_MH"
    )
  )
  expect_setequal(
    subsets[["MSigDB_M2"]][["TERM2GENE"]][["Term"]],
    c("CGP1", "BIO1", "REACT1")
  )
  expect_identical(subsets[["MSigDB_M2_CGP"]][["TERM2GENE"]][["symbol"]], c("A", "B"))
  expect_setequal(subsets[["MSigDB_M2_CP"]][["TERM2GENE"]][["Term"]], c("BIO1", "REACT1"))
  expect_false("CGP1" %in% subsets[["MSigDB_M2_CP"]][["TERM2GENE"]][["Term"]])
  expect_identical(subsets[["MSigDB_M2_CP_BIOCARTA"]][["TERM2GENE"]][["Term"]], "BIO1")
  expect_identical(subsets[["MSigDB_MH"]][["TERM2GENE"]][["Term"]], "H1")
})

test_that("unknown MSigDB names fail with the available collections", {
  expect_error(
    preparedb_require_msigdb_names(
      db_list = list(Mus_musculus = list(MSigDB_M2 = list())),
      species = "Mus_musculus",
      db = "MSigDB_M2_NOPE"
    ),
    "MSigDB_M2_NOPE"
  )
  expect_error(
    preparedb_require_msigdb_names(
      db_list = list(Mus_musculus = list(MSigDB_M2 = list())),
      species = "Mus_musculus",
      db = "MSigDB_M2_NOPE"
    ),
    "Available MSigDB databases"
  )
  expect_null(
    preparedb_require_msigdb_names(
      db_list = list(Mus_musculus = list()),
      species = "Mus_musculus",
      db = "GO_BP"
    )
  )
})

test_that("PrepareDB keeps requested MSigDB names and drops siblings", {
  db_list <- list(
    Mus_musculus = list(
      MSigDB = list(TERM2GENE = "all"),
      MSigDB_M2 = list(TERM2GENE = "m2"),
      MSigDB_M2_CGP = list(TERM2GENE = "cgp"),
      MSigDB_M2_CP = list(TERM2GENE = "cp"),
      GO_BP = list(TERM2GENE = "go"),
      "MSigDB_M2:CGP" = list(TERM2GENE = "cgp")
    ),
    Homo_sapiens = list(MSigDB = list(TERM2GENE = "human")),
    CytoTRACE2 = list(version = "1.1.0")
  )
  out <- preparedb_keep_requested_dbs(
    db_list = db_list,
    species = "Mus_musculus",
    db_names = c("MSigDB_M2_CGP", "MSigDB_M2:CGP", "GO_BP")
  )
  expect_identical(
    names(out[["Mus_musculus"]]),
    c("MSigDB_M2_CGP", "GO_BP", "MSigDB_M2:CGP")
  )
  expect_false("Homo_sapiens" %in% names(out))
  expect_identical(out[["CytoTRACE2"]][["version"]], "1.1.0")
})

test_that("MSigDB parent map is used only when it has the requested IDs", {
  parent <- data.frame(
    Term = c("CGP1", "CGP1"),
    symbol = c("A", "B"),
    stringsAsFactors = FALSE
  )
  expect_null(preparedb_msigdb_parent_term2gene(NULL, "entrez_id"))
  expect_null(preparedb_msigdb_parent_term2gene(parent, "entrez_id"))
  expect_identical(
    preparedb_msigdb_parent_term2gene(parent, "symbol"),
    parent
  )
})

test_that("colon MSigDB requests alias the underscore database", {
  entry <- list(TERM2GENE = data.frame(Term = "CGP1", symbol = "A"))
  db_list <- list(Mus_musculus = list(MSigDB_M2_CGP = entry))
  out <- preparedb_alias_requested_msigdb(
    db_list = db_list,
    species = "Mus_musculus",
    db_requested = "MSigDB_M2:CGP",
    db_normalized = "MSigDB_M2_CGP"
  )
  expect_identical(
    out[["Mus_musculus"]][["MSigDB_M2:CGP"]],
    out[["Mus_musculus"]][["MSigDB_M2_CGP"]]
  )
})
