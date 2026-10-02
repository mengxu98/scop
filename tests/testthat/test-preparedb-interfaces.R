preparedb_mock_cache <- function() {
  annotation <- function(db) {
    list(
      TERM2GENE = data.frame(Term = paste0(db, "_term"), symbol = "G1"),
      TERM2NAME = data.frame(Term = paste0(db, "_term"), Name = "Example"),
      version = "test"
    )
  }
  local_mocked_bindings(
    list_db_cache_entries = function(species, db) {
      expand.grid(Species = species, DB = db, stringsAsFactors = FALSE) |>
        transform(
          timestamp = as.POSIXct("2026-01-01", tz = "UTC"),
          file = paste(Species, DB, sep = "@")
        )
    },
    .package = "scop", .env = parent.frame()
  )
  local_mocked_bindings(
    readCacheHeader = function(pathname, ...) {
      fields <- strsplit(pathname, "@", fixed = TRUE)[[1]]
      list(
        comment = paste0("test nterm:1|", fields[[1]], "-", fields[[2]]),
        timestamp = as.POSIXct("2026-01-01", tz = "UTC")
      )
    },
    loadCache = function(pathname, ...) annotation(strsplit(pathname, "@", fixed = TRUE)[[1]][[2]]),
    .package = "R.cache", .env = parent.frame()
  )
  invisible(annotation)
}

test_that("source interfaces prepare cached annotations independently of PrepareDB", {
  skip_if_not_installed("R.cache")
  annotation <- preparedb_mock_cache()
  local_mocked_bindings(
    PrepareDB = function(...) log_message("Unexpected dispatcher recursion", message_type = "error"),
    .package = "scop"
  )
  selectors <- c(
    "GO", "KEGG", "WikiPathway", "Reactome", "CORUM", "MP", "DO", "HPO", "PFAM",
    "Chromosome", "GeneType", "Enzyme", "TF", "CSPA", "Surfaceome", "SPRomeDB", "VerSeDa",
    "TFLink", "hTFtarget", "TRRUST", "JASPAR", "ENCODE", "MSigDB", "CellTalk", "CellChat"
  )
  for (db in selectors) {
    name <- if (db == "hTFtarget") "PrepareHTFtarget" else paste0("Prepare", db)
    prepare <- get(name, envir = asNamespace("scop"))
    result <- prepare(species = "Mus_musculus", db = db, db_IDtypes = "symbol", verbose = FALSE)
    expect_identical(result, list(Mus_musculus = setNames(list(annotation(db)), db)), info = name)
    expect_error(prepare(db = "Unknown", verbose = FALSE), "selector", info = name)
  }
})

test_that("combined preparation preserves cached order and MSigDB aliases", {
  skip_if_not_installed("R.cache")
  annotation <- preparedb_mock_cache()
  requested <- c("TF", "GO_CC", "MSigDB_M2:CP:BIOCARTA", "KEGG")
  result <- PrepareDB(species = "Mus_musculus", db = requested, db_IDtypes = "symbol", verbose = FALSE)
  expect_identical(names(result[["Mus_musculus"]]), c("TF", "GO_CC", "MSigDB_M2_CP_BIOCARTA", "KEGG", "MSigDB_M2:CP:BIOCARTA"))
  expect_identical(result[["Mus_musculus"]][["MSigDB_M2:CP:BIOCARTA"]], annotation("MSigDB_M2_CP_BIOCARTA"))
  go <- PrepareGO(species = "Mus_musculus", db = c("GO_CC", "GO_BP"), db_IDtypes = "symbol", verbose = FALSE)
  expect_identical(names(go[["Mus_musculus"]]), c("GO_CC", "GO_BP"))
  expect_identical(go[["Mus_musculus"]][["GO_BP"]], annotation("GO_BP"))
})

test_that("related selectors share one source preparation call", {
  skip_if_not_installed("R.cache")
  annotation <- preparedb_mock_cache()
  original <- PrepareGO
  calls <- 0L
  local_mocked_bindings(PrepareGO = function(...) {
    calls <<- calls + 1L
    original(...)
  }, .package = "scop")
  requested <- c("GO_CC", "GO_BP", "GO_MF", "GO")
  result <- PrepareDB(species = "Mus_musculus", db = requested, db_IDtypes = "symbol", verbose = FALSE)
  expect_identical(calls, 1L)
  expect_identical(names(result$Mus_musculus), requested)
  expect_identical(result$Mus_musculus, setNames(lapply(requested, annotation), requested))
})

test_that("source parsing and cross-species conversion use the shared annotation contract", {
  skip_if_not_installed("R.cache")
  skip_if_not_installed("clusterProfiler")
  directory <- tempfile("corum-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  writeLines("Complex\tDescription\tG1\tG2", file.path(directory, "gene_set_library_crisp.gmt"))
  saved <- list()
  local_mocked_bindings(saveCache = function(object, key, comment, ...) {
    saved[[length(saved) + 1L]] <<- list(object = object, key = key, comment = comment)
  }, .package = "R.cache")
  local_mocked_bindings(
    preparedb_local_orgdb_id_map = function(...) NULL,
    GeneConvert = function(geneID, geneID_from_IDtype, geneID_to_IDtype, species_from,
                           species_to, Ensembl_version, mirror, biomart, max_tries) {
      data <- data.frame(ensembl_id = paste0("E", seq_along(geneID)), row.names = geneID)
      list(geneID_res = data, geneID_collapse = data)
    }, .package = "scop"
  )
  human <- PrepareCORUM(species = "Homo_sapiens", db_IDtypes = "symbol", data_dir = directory, db_update = TRUE, verbose = FALSE)
  expect_identical(human[["Homo_sapiens"]]$CORUM$TERM2GENE$symbol, c("G1", "G2"))
  expect_identical(saved[[1]]$key, list("Harmonizome 3.0", "Homo_sapiens", "CORUM"))
  expect_identical(saved[[1]]$comment, "Harmonizome 3.0 nterm:1|Homo_sapiens-CORUM")
  mouse <- PrepareCORUM(species = "Mus_musculus", db_IDtypes = "ensembl_id", data_dir = directory, db_update = TRUE, verbose = FALSE)
  expect_identical(names(mouse), "Mus_musculus")
  expect_identical(mouse[["Mus_musculus"]]$CORUM$TERM2GENE$ensembl_id, c("E1", "E2"))
  expect_identical(mouse[["Mus_musculus"]]$CORUM$version, "Harmonizome 3.0(converted from Homo_sapiens)")
  expect_identical(saved[[2]]$key, list("Harmonizome 3.0(converted from Homo_sapiens)", "Mus_musculus", "CORUM"))
  expect_error(PrepareCORUM(species = "Mus_musculus", convert_species = FALSE, db_update = TRUE, data_dir = directory, verbose = FALSE), "Stop")
})

test_that("version selection and fallback share one cache policy", {
  skip_if_not_installed("R.cache")
  local_mocked_bindings(list_db_cache_entries = function(species, db) {
    data.frame(
      Species = species, DB = db, db_version = c("v2", "v1"),
      timestamp = as.POSIXct(c("2026-02-01", "2026-01-01"), tz = "UTC"), file = c("v2", "v1")
    )
  }, .package = "scop")
  local_mocked_bindings(
    readCacheHeader = function(pathname, ...) list(comment = paste0(pathname, "|Mus_musculus-KEGG"), timestamp = as.POSIXct("2026-01-01", tz = "UTC")),
    loadCache = function(pathname, ...) list(TERM2GENE = data.frame(Term = "T", symbol = "G1"), TERM2NAME = data.frame(Term = "T", Name = "Example"), version = pathname),
    .package = "R.cache"
  )
  selected <- PrepareKEGG(species = "Mus_musculus", db_IDtypes = "symbol", db_version = "v1", verbose = FALSE)
  expect_identical(selected$Mus_musculus$KEGG$version, "v1")
  latest <- PrepareKEGG(species = "Mus_musculus", db_IDtypes = "symbol", verbose = FALSE)
  expect_identical(latest$Mus_musculus$KEGG$version, "v2")
  expect_warning(fallback <- PrepareKEGG(species = "Mus_musculus", db_IDtypes = "symbol", db_version = "v3", verbose = TRUE), "no.*v3")
  expect_identical(fallback$Mus_musculus$KEGG$version, "v2")
})

test_that("an unavailable species can skip preparation without losing earlier resources", {
  skip_if_not_installed("R.cache")
  skip_if_not_installed("httr")
  local_mocked_bindings(
    kegg_get = function(url) {
      if (grepl("list/genome", url)) {
        return(data.frame(ID = "T1", Name = "hsa; Homo sapiens (human)"))
      }
      if (grepl("link/hsa", url)) {
        return(data.frame(Pathway = "path:hsa00010", Gene = "hsa:1"))
      }
      if (grepl("conv/ncbi-geneid", url)) {
        return(data.frame(Gene = "hsa:1", Entrez = "ncbi-geneid:1"))
      }
      data.frame(Pathway = "path:hsa00010", Name = "Glycolysis - Homo sapiens (human)")
    }, .package = "scop"
  )
  local_mocked_bindings(GET = function(...) "Release 119.0", content = identity, .package = "httr")
  local_mocked_bindings(saveCache = function(...) NULL, .package = "R.cache")
  result <- PrepareKEGG(species = c("Homo_sapiens", "Mus_musculus"), db_IDtypes = "entrez_id", db_update = TRUE, verbose = FALSE)
  expect_identical(names(result), "Homo_sapiens")
  expect_identical(result$Homo_sapiens$KEGG$TERM2GENE$entrez_id, "1")
  expect_error(PrepareKEGG(species = "Mus_musculus", convert_species = FALSE, db_update = TRUE, verbose = FALSE), "Stop")
})

test_that("model resources retain the existing cache key and top-level contract", {
  skip_if_not_installed("R.cache")
  directory <- tempfile("model-cache-")
  dir.create(file.path(directory, "CytoTRACE2"), recursive = TRUE)
  directory <- normalizePath(directory, winslash = "/", mustWork = TRUE)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  files <- c("model_parameters.rds", "features_model_training_17.csv", "mt_dict_human_to_mouse.csv", "mt_human_alias.csv", "mt_mouse_alias.csv")
  for (file in files) writeLines("cached", file.path(directory, "CytoTRACE2", file))
  model <- list(data_dir = file.path(directory, "CytoTRACE2"), files = files, version = "1.1.0")
  local_mocked_bindings(R_user_dir = function(package, which) directory, .package = "tools")
  saved <- list()
  local_mocked_bindings(
    loadCache = function(key, ...) {
      expect_identical(key, list("1.1.0", "CytoTRACE2", "CytoTRACE2"))
      model
    },
    saveCache = function(object, key, comment, ...) saved <<- list(object = object, key = key, comment = comment),
    .package = "R.cache"
  )
  expect_identical(PrepareCytoTRACE2(verbose = FALSE), list(CytoTRACE2 = model))
  expect_identical(PrepareDB(db = "CytoTRACE2", verbose = FALSE), list(CytoTRACE2 = model))
  downloaded <- character()
  local_mocked_bindings(download.file = function(url, destfile, mode, quiet) {
    downloaded <<- c(downloaded, basename(url))
    writeLines(url, destfile)
    0L
  }, .package = "utils")
  before <- getOption("timeout")
  refreshed <- PrepareCytoTRACE2(db_update = TRUE, verbose = FALSE)
  expect_identical(refreshed, list(CytoTRACE2 = model))
  expect_identical(downloaded, files)
  expect_identical(saved$key, list("1.1.0", "CytoTRACE2", "CytoTRACE2"))
  expect_identical(saved$comment, "1.1.0 nterm:5|CytoTRACE2-CytoTRACE2")
  expect_identical(getOption("timeout"), before)
})

test_that("custom preparation has an independent public entry point", {
  skip_if_not_installed("R.cache")
  saved <- NULL
  local_mocked_bindings(saveCache = function(object, key, comment, ...) saved <<- list(key = key, comment = comment), .package = "R.cache")
  mappings <- data.frame(Term = c("Cycle", "Cycle"), symbol = c("G1", "G2"))
  result <- PrepareCustomDB(
    species = "Mus_musculus", db = "Cycle", db_IDtypes = "symbol",
    custom_TERM2GENE = mappings, custom_species = "Mus_musculus", custom_IDtype = "symbol", custom_version = "v1", verbose = FALSE
  )
  expect_identical(names(result), "Mus_musculus")
  expect_identical(result$Mus_musculus$Cycle$TERM2GENE, mappings)
  expect_identical(saved$key, list("v1", "Mus_musculus", "Cycle"))
  expect_error(PrepareCustomDB(db = c("Cycle", "Other"), custom_TERM2GENE = mappings, verbose = FALSE), "length")
})

test_that("custom annotations can be prepared alongside model resources", {
  skip_if_not_installed("R.cache")
  model <- list(data_dir = tempdir(), files = "model_parameters.rds", version = "1.1.0")
  local_mocked_bindings(
    PrepareCytoTRACE2 = function(db_update, verbose) {
      expect_true(db_update)
      expect_false(verbose)
      list(CytoTRACE2 = model)
    }, .package = "scop"
  )
  saved <- NULL
  local_mocked_bindings(saveCache = function(object, key, comment, ...) saved <<- list(key = key, comment = comment), .package = "R.cache")
  mappings <- data.frame(Term = c("Cycle", "Cycle"), symbol = c("G1", "G2"))
  custom <- PrepareCustomDB(
    species = "Mus_musculus", db = "Cycle", db_IDtypes = "symbol", db_update = TRUE,
    custom_TERM2GENE = mappings, custom_species = "Mus_musculus", custom_IDtype = "symbol", custom_version = "v1", verbose = FALSE
  )
  expect_identical(PrepareDB(
    species = "Mus_musculus", db = "Cycle", db_IDtypes = "symbol", db_update = TRUE,
    custom_TERM2GENE = mappings, custom_species = "Mus_musculus", custom_IDtype = "symbol", custom_version = "v1", verbose = FALSE
  ), custom)
  for (selectors in list(c("CytoTRACE2", "Cycle"), c("Cycle", "CytoTRACE2"))) {
    combined <- PrepareDB(
      species = "Mus_musculus", db = selectors, db_IDtypes = "symbol", db_update = TRUE,
      custom_TERM2GENE = mappings, custom_species = "Mus_musculus", custom_IDtype = "symbol", custom_version = "v1", verbose = FALSE
    )
    expect_identical(combined, c(list(CytoTRACE2 = model), custom))
    expect_identical(saved$key, list("v1", "Mus_musculus", "Cycle"))
  }
})
