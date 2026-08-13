# ListCCCDB and PrepareCCCDB unit tests.

test_that("ListCCCDB enumerates CCC databases with a unified schema", {
  liana_loaded <- "liana" %in% loadedNamespaces()
  tcltk_loaded <- "tcltk" %in% loadedNamespaces()
  resources <- ListCCCDB()
  expect_identical("liana" %in% loadedNamespaces(), liana_loaded)
  expect_identical("tcltk" %in% loadedNamespaces(), tcltk_loaded)
  expect_identical(colnames(resources), c("Database", "Resource", "Species", "Default", "Status", "Description"))
  expect_setequal(unique(resources$Database), c("LIANA", "CellChat", "Nichenetr", "CellphoneDB"))
  available <- resources[resources$Status == "available", , drop = FALSE]
  expect_true(all(available$Species %in% c("human", "mouse", "zebrafish")))
  cellchat_rows <- resources[resources$Database == "CellChat", , drop = FALSE]
  testthat::skip_if_not_installed("CellChat")
  expect_setequal(
    cellchat_rows$Resource,
    c("CellChatDB.human", "CellChatDB.mouse", "CellChatDB.zebrafish")
  )
  expect_true(all(cellchat_rows$Default))
  expect_true(all(resources$Status %in% c(
    "available",
    "backend_not_installed",
    "requires_python_cellphonedb"
  )))
})

test_that("ListCCCDB returns a plain data frame", {
  resources <- ListCCCDB()
  expect_false(inherits(resources, "ccc_db_table"))
  df <- as.data.frame(resources)
  expect_identical(colnames(df), colnames(resources))
})

test_that("ListCCCDB filters by db and species", {
  resources <- ListCCCDB(db = c("CellChat", "LIANA"))
  expect_setequal(unique(resources$Database), c("CellChat", "LIANA"))
  testthat::skip_if_not_installed("CellChat")
  mouse <- ListCCCDB(species = "mouse")
  expect_true(all(mouse$Species == "mouse"))
  expect_true(any(mouse$Resource == "CellChatDB.mouse"))
  expect_false(any(mouse$Resource == "CellChatDB.human"))
})

test_that("ListDB returns unified Database/Species/Version/Date columns", {
  testthat::skip_if_not_installed("R.cache")
  dbinfo <- ListDB(species = "Homo_sapiens", db = "CellChat")
  expect_identical(colnames(dbinfo), c("Database", "Species", "Version", "Date"))
  expect_false(inherits(dbinfo, "db_table"))
  expect_true(all(dbinfo$Database == "CellChat"))
  expect_true(all(nzchar(dbinfo$Version)))
  expect_true(all(nzchar(dbinfo$Date)))
})

test_that("ListCCCDB reports unavailable backends without failing", {
  resources <- ListCCCDB(db = "CellphoneDB")
  expect_true(nrow(resources) >= 1L)
  expect_true(all(resources$Database == "CellphoneDB"))
  if (!all(resources$Status == "available")) {
    expect_true(any(resources$Status %in% c("backend_not_installed", "requires_python_cellphonedb")))
  }
})

test_that("PrepareCCCDB validates db and species arguments", {
  expect_error(PrepareCCCDB(db = "KEGG"), "should be one of")
  expect_error(PrepareCCCDB(species = "Canis_lupus"), "should be one of")
})

test_that("PrepareDB delegates CellChat and CellTalk to PrepareCCCDB", {
  testthat::skip_if_not_installed("biomaRt")
  mock_env <- environment()
  captured <- list()
  testthat::local_mocked_bindings(
    PrepareCCCDB = function(species, db, convert_species = TRUE, data_dir = NULL,
                            db_version = "latest", db_update = FALSE, verbose = TRUE, ...) {
      local_captured <- get0("captured", envir = mock_env, inherits = FALSE)
      local_captured$species <- species
      local_captured$db <- db
      assign("captured", local_captured, envir = mock_env)
      out <- list()
      out[[species]] <- stats::setNames(
        lapply(db, function(term) {
          list(
            TERM2GENE = data.frame(Term = paste0("ligand_", term), symbol = "GENE"),
            TERM2NAME = data.frame(Term = paste0("ligand_", term), Name = term),
            version = "mock"
          )
        }),
        db
      )
      out
    },
    .package = "scop"
  )

  out <- PrepareDB(
    species = "Homo_sapiens",
    db = "CellChat",
    db_update = TRUE,
    verbose = FALSE
  )
  expect_equal(captured$species, "Homo_sapiens")
  expect_equal(captured$db, "CellChat")
  expect_true("CellChat" %in% names(out$Homo_sapiens))
  expect_true(all(c("TERM2GENE", "TERM2NAME", "version") %in% names(out$Homo_sapiens$CellChat)))
})
