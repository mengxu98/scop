species_compatibility_object <- function() {
  counts <- matrix(seq_len(24 * 6), nrow = 24,
    dimnames = list(paste0("G", 1:24), paste0("C", 1:6)))
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object$Sample_Origin <- rep(c("A", "B"), each = 3)
  object
}

species_compatibility_db <- function(species, db, ...) {
  entry <- list(
    TERM2GENE = data.frame(Term = "T00010", symbol = paste0("G", 1:12)),
    TERM2NAME = data.frame(Term = "T00010", Name = "Test pathway"),
    version = "fixture"
  )
  setNames(list(setNames(rep(list(entry), length(db)), db)), normalize_species_name(species))
}

test_that("RunGSVA accepts species aliases in grouped and single-cell database modes", {
  local_mocked_bindings(PrepareDB = species_compatibility_db, .package = "scop")
  object <- species_compatibility_object()
  for (group in list("Sample_Origin", NULL)) {
    run <- function(species, combine = FALSE) RunGSVA(object,
      group.by = group, db = "MSigDB_H", species = species,
      db_combine = combine, method = "zscore", layer = "counts",
      new_assay = FALSE, store_metadata = FALSE, verbose = FALSE)
    canonical <- run("Homo_sapiens")
    key <- if (is.null(group)) "GSVA_cell_zscore" else "GSVA_Sample_Origin_zscore"
    for (alias in c("homo_sapiens", "HOMO SAPIENS", " Homo.sapiens ", "Homo-sapiens")) {
      out <- run(alias)
      expect_equal(out@tools[[key]]$scores, canonical@tools[[key]]$scores)
      expect_identical(out@tools[[key]]$species, "Homo_sapiens")
    }
    expect_equal(run("homo_sapiens", TRUE)@tools[[key]]$scores,
      run("Homo_sapiens", TRUE)@tools[[key]]$scores)
  }
})

test_that("CellScoring species aliases retain database-derived scores", {
  local_mocked_bindings(PrepareDB = species_compatibility_db, .package = "scop")
  object <- species_compatibility_object()
  run <- function(species) CellScoring(object, species = species,
    db = "MSigDB_H", method = "zscore", layer = "counts",
    classification = FALSE, name = "score", verbose = FALSE)
  expect_equal(run("homo sapiens")@meta.data, run("Homo_sapiens")@meta.data)
})

test_that("AnnotateFeatures species aliases append the same database annotations", {
  local_mocked_bindings(PrepareDB = species_compatibility_db, .package = "scop")
  object <- species_compatibility_object()
  run <- function(species) AnnotateFeatures(object, species = species,
    db = "MSigDB_H", verbose = FALSE)
  out <- GetFeaturesData(run("homo.sapiens"))
  expect_equal(out, GetFeaturesData(run("Homo_sapiens")))
  expect_true("MSigDB_H" %in% colnames(out))
  expect_equal(out["G1", "MSigDB_H"], "Test pathway")
})

test_that("CellRank enrichment species aliases do not mark available databases unavailable", {
  local_mocked_bindings(
    PrepareDB = species_compatibility_db,
    check_r = function(...) TRUE,
    get_namespace_fun = function(...) function(...) data.frame(ID = "T00010", Description = "Test pathway"),
    .package = "scop"
  )
  object <- species_compatibility_object()
  object@tools$CellRank <- list(trends = list(A = list(
    cluster_table = data.frame(gene = paste0("G", 1:12), cluster = "0")
  )))
  run <- function(species) RunCellRankEnrichment(object, lineage = "A",
    db = "MSigDB_H", species = species, verbose = FALSE)
  out <- run("mus-musculus")@tools$CellRank$enrichment$A
  expect_equal(out, run("Mus_musculus")@tools$CellRank$enrichment$A)
  expect_identical(out$manifest$MSigDB_H$status, "ok")
})

test_that("dynamic enrichment species aliases retain pathway feature lists", {
  calls <- list()
  local_mocked_bindings(
    PrepareDB = species_compatibility_db,
    CellScoring = function(object, features, name, ...) {
      calls[[length(calls) + 1L]] <<- features
      object[[name]] <- SeuratObject::CreateAssayObject(counts = matrix(1,
        nrow = length(features), ncol = ncol(object),
        dimnames = list(names(features), colnames(object))))
      object
    },
    RunDynamicFeatures = function(object, ...) object,
    .package = "scop"
  )
  object <- species_compatibility_object()
  object@tools$DynamicFeatures_A <- list(DynamicFeatures = data.frame(
    features = paste0("G", 1:12), exp_ncells = 30, r.sq = 0.9,
    dev.expl = 0.9, padjust = 0.001
  ))
  for (species in c("Homo_sapiens", "homo sapiens")) {
    RunDynamicEnrichment(object, lineages = "A", species = species,
      db = "MSigDB_H", verbose = FALSE)
  }
  expect_length(calls, 2L)
  expect_identical(calls[[1]], calls[[2]])
  expect_identical(calls[[2]][["Test pathway"]], paste0("G", 1:12))
})

test_that("metabolism database loading accepts species aliases without dropping pathways", {
  local_mocked_bindings(PrepareDB = species_compatibility_db, .package = "scop")
  run <- function(species) build_metabolism_gene_sets_from_preparedb(
    species = species, db_prepare = "KEGG", IDtype = "symbol",
    curated = list(kegg_refs = "00010"), expr_gene_names = paste0("G", 1:24),
    db_update = FALSE, db_version = "latest", convert_species = TRUE,
    Ensembl_version = NULL, mirror = NULL, minGSSize = 10, maxGSSize = 500,
    verbose = FALSE
  )
  out <- run("HOMO SAPIENS")
  expect_equal(out, run("Homo_sapiens"))
  expect_identical(out$gene_sets$T00010, paste0("G", 1:12))
})

test_that("RunMetabolism stores canonical species with unchanged pathway scores", {
  local_mocked_bindings(
    PrepareDB = species_compatibility_db,
    scmetabolism_pathway_refs = function(...) list(kegg_refs = "00010"),
    .package = "scop"
  )
  object <- species_compatibility_object()
  run <- function(species) RunMetabolism(object, species = species,
    db = "KEGG", method = "ssGSEA", new_assay = FALSE, verbose = FALSE)
  out <- run("homo-sapiens")@tools$Metabolism_ssGSEA
  expect_equal(out, run("Homo_sapiens")@tools$Metabolism_ssGSEA)
  expect_identical(out$species, "Homo_sapiens")
  expect_length(out$gene_sets, 1L)
})
