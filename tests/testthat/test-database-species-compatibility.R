species_compatibility_object <- function() {
  counts <- matrix(seq_len(24 * 6),
    nrow = 24,
    dimnames = list(paste0("G", 1:24), paste0("C", 1:6))
  )
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
    run <- function(species, combine = FALSE) {
      RunGSVA(object,
        group.by = group, db = "MSigDB_H", species = species,
        db_combine = combine, method = "ssgsea", backend = "cpp", layer = "counts",
        new_assay = FALSE, store_metadata = FALSE, verbose = FALSE
      )
    }
    canonical <- run("Homo_sapiens")
    key <- if (is.null(group)) "GSVA_cell_ssgsea" else "GSVA_Sample_Origin_ssgsea"
    for (alias in c("homo_sapiens", "HOMO SAPIENS", " Homo.sapiens ", "Homo-sapiens")) {
      out <- run(alias)
      expect_equal(out@tools[[key]]$scores, canonical@tools[[key]]$scores)
      expect_identical(out@tools[[key]]$species, "Homo_sapiens")
    }
    expect_equal(
      run("homo_sapiens", TRUE)@tools[[key]]$scores,
      run("Homo_sapiens", TRUE)@tools[[key]]$scores
    )
  }
})

test_that("CellScoring species aliases retain database-derived scores", {
  local_mocked_bindings(PrepareDB = species_compatibility_db, .package = "scop")
  object <- species_compatibility_object()
  run <- function(species) {
    CellScoring(object,
      species = species,
      db = "MSigDB_H", method = "zscore", layer = "counts",
      classification = FALSE, name = "score", verbose = FALSE
    )
  }
  expect_equal(run("homo sapiens")@meta.data, run("Homo_sapiens")@meta.data)
})

test_that("AnnotateFeatures species aliases append the same database annotations", {
  local_mocked_bindings(PrepareDB = species_compatibility_db, .package = "scop")
  object <- species_compatibility_object()
  run <- function(species) {
    AnnotateFeatures(object,
      species = species,
      db = "MSigDB_H", verbose = FALSE
    )
  }
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
  run <- function(species) {
    RunCellRankEnrichment(object,
      lineage = "A",
      db = "MSigDB_H", species = species, verbose = FALSE
    )
  }
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
        dimnames = list(names(features), colnames(object))
      ))
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
    RunDynamicEnrichment(object,
      lineages = "A", species = species,
      db = "MSigDB_H", verbose = FALSE
    )
  }
  expect_length(calls, 2L)
  expect_identical(calls[[1]], calls[[2]])
  expect_identical(calls[[2]][["Test pathway"]], paste0("G", 1:12))
})

test_that("metabolism database loading accepts species aliases without dropping pathways", {
  local_mocked_bindings(PrepareDB = species_compatibility_db, .package = "scop")
  run <- function(species) {
    build_metabolism_gene_sets_from_preparedb(
      species = species, db_prepare = "KEGG", IDtype = "symbol",
      curated = list(kegg_refs = "00010"), expr_gene_names = paste0("G", 1:24),
      db_update = FALSE, db_version = "latest", convert_species = TRUE,
      Ensembl_version = NULL, mirror = NULL, minGSSize = 10, maxGSSize = 500,
      verbose = FALSE
    )
  }
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
  run <- function(species) {
    RunMetabolism(object,
      species = species,
      db = "KEGG", method = "ssGSEA", new_assay = FALSE, verbose = FALSE
    )
  }
  out <- run("homo-sapiens")@tools$Metabolism_ssGSEA
  expect_equal(out, run("Homo_sapiens")@tools$Metabolism_ssGSEA)
  expect_identical(out$species, "Homo_sapiens")
  expect_length(out$gene_sets, 1L)
})

test_that("AnnotateFeatures joins standard GTF gene records exactly by symbol", {
  object <- species_compatibility_object()
  gtf <- withr::local_tempfile(fileext = ".gtf")
  writeLines(c(
    "# example gene annotation",
    'chr1\tsource\ttranscript\t13\t18\t.\t+\t.\tgene_id "ID2.1"; gene_name "G2";',
    'chr2\tsource\tgene\t11\t21\t.\t+\t.\tgene_id "ID2.1"; gene_name "G2"; gene_type "protein coding; curated #1";',
    'chr2\tsource\tgene\t11\t21\t.\t+\t.\tgene_id "ID2.1"; gene_name "G2";',
    'chr3\tsource\tgene\t31\t41\t.\t-\t.\tgene_id "ID2.1"; gene_name "G2";',
    "# an interleaved comment",
    'chr4\tsource\tgene\t51\t61\t.\t+\t.\tgene_id "ID12.2"; gene_name "G12";',
    'chr5\tsource\tgene\t71\t81\t.\t+\t.\tgene_id "ID71";'
  ), gtf)
  annotated <- AnnotateFeatures(object, gtf = gtf, verbose = FALSE)
  result <- GetFeaturesData(annotated)
  expect_identical(rownames(result), rownames(object))
  expect_identical(
    unname(unlist(result["G2", ], use.names = FALSE)),
    c("chr2;chr3", "gene", "11;31", "21;41", "+;-", "ID2.1", "G2", "protein coding; curated #1")
  )
  expect_identical(result["G12", "gene_id"], "ID12.2")
  expect_identical(result["G12", "gene_type"], NA_character_)
  expect_true(all(is.na(result["G1", ])))
  expect_true(methods::validObject(annotated, test = TRUE))
  expect_identical(SeuratObject::LayerData(annotated, layer = "counts"),
                   SeuratObject::LayerData(object, layer = "counts"))
})

test_that("AnnotateFeatures accepts GTF column selection and an independent key", {
  counts <- Matrix::sparseMatrix(
    i = 1:3, j = 1:3, x = 1L, dims = c(3L, 3L),
    dimnames = list(c("ID2.1", "ID1", "ID1.6"), c("C1", "C2", "C3"))
  )
  object <- SeuratObject::CreateSeuratObject(counts)
  object <- AddFeaturesData(
    object,
    data.frame(gene_name = rep("existing", 3L), keep = 1:3, row.names = rownames(counts))
  )
  gtf <- withr::local_tempfile(fileext = ".gtf.gz")
  conn <- gzfile(gtf, "wt")
  tryCatch(writeLines(c(
    "# transcript-only fallback",
    'chr1\tensembl\ttranscript\t1\t3\t.\t+\t.\tgene_id "ID1.6"; gene_name "first gene"; gene_biotype "type A";',
    'chr1\tensembl\ttranscript\t4\t5\t.\t+\t.\tgene_biotype "type B"; gene_name "second gene"; gene_id "ID2.1";',
    'chr1\tensembl\texon\t2\t2\t.\t+\t.\tgene_id "ID1.6"; gene_name "not used";'
  ), conn), finally = close(conn))

  run <- function(overwrite) {
    GetFeaturesData(AnnotateFeatures(
      object, gtf = gtf, IDtype = "ensembl_id", merge_gtf_by = "gene_id",
      columns = c("source", "gene_name", "gene_biotype"), overwrite = overwrite,
      verbose = FALSE
    ))
  }
  preserved <- run(FALSE)
  expect_identical(preserved$gene_name, rep("existing", 3L))
  replaced <- run(TRUE)
  expect_identical(rownames(replaced), rownames(counts))
  expect_identical(names(replaced), c("gene_name", "keep", "source", "gene_biotype"))
  expect_identical(replaced$keep, 1:3)
  expect_identical(replaced$source, c("ensembl", NA, "ensembl"))
  expect_identical(replaced$gene_name, c("second gene", NA, "first gene"))
  expect_identical(replaced$gene_biotype, c("type B", NA, "type A"))
  expect_error(
    AnnotateFeatures(object, gtf = gtf, merge_gtf_by = "absent_key", verbose = FALSE),
    "^No GTF records contain values"
  )
})
