test_that("PrepareDB reuses matching SCENICPlus reference files", {
  directory <- withr::local_tempdir()
  annotation <- withr::local_tempdir()
  writeLines(c(
    "Chromosome\tStart\tEnd\tStrand\tGene\tTranscription_Start_Site",
    "chr1\t10\t20\t+\tA\t10"
  ), file.path(annotation, "genome_annotation.tsv"))
  writeLines(c("Chromosome\tStart\tEnd", "chr1\t0\t100"), file.path(annotation, "chromsizes.tsv"))
  prefix <- "hg38_screen_v10_clust.regions_vs_motifs"
  for (type in c("rankings", "scores")) {
    writeBin(charToRaw("ARROW1fixtureARROW1"), file.path(directory, paste0(prefix, ".", type, ".feather")))
  }
  motif <- file.path(directory, "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")
  writeLines(c("#motif_id\tgene_name", "m1\tTF1"), motif)
  source <- list(SCENICPlus = list(
    path = directory,
    genome_annotation = file.path(annotation, "genome_annotation.tsv"),
    chromsizes = file.path(annotation, "chromsizes.tsv")
  ))
  testthat::local_mocked_bindings(scenic_download_reference_file = function(...) {
    stop("unexpected download")
  }, .package = "scop")
  result <- PrepareDB(db = "SCENICPlus", data_dir = source, verbose = FALSE)$SCENICPlus
  expect_named(result, c("ranking_db", "scores_db", "motif_annotation", "genome_annotation", "chromsizes"))
  expect_identical(result$motif_annotation, normalizePath(motif))
  expect_identical(result$genome_annotation, normalizePath(source$SCENICPlus$genome_annotation))
  expect_error(PrepareDB(species = "Mus_musculus", db = "SCENICPlus"), "Homo_sapiens")
  writeLines("bad", source$SCENICPlus$chromsizes)
  expect_error(PrepareDB(db = "SCENICPlus", data_dir = source, verbose = FALSE), "chromsizes columns")
})

test_that("PrepareDB checks a SCENICPlus download before replacing a cached file", {
  directory <- withr::local_tempdir()
  writeLines(
    "Chromosome\tStart\tEnd\tStrand\tGene\tTranscription_Start_Site",
    file.path(directory, "genome_annotation.tsv")
  )
  writeLines("Chromosome\tStart\tEnd", file.path(directory, "chromsizes.tsv"))
  ranking <- file.path(directory, "hg38_screen_v10_clust.regions_vs_motifs.rankings.feather")
  writeBin(charToRaw("ARROW1oldARROW1"), ranking)
  download <- character()
  testthat::local_mocked_bindings(download.file = function(url, destfile, mode, quiet) {
    download <<- c(download, url)
    writeBin(charToRaw("partial download"), destfile)
    0L
  }, .package = "utils")
  expect_error(PrepareDB(
    db = "SCENICPlus", data_dir = directory, db_update = TRUE,
    verbose = FALSE
  ), "Invalid SCENICPlus Feather file")
  expect_length(download, 1L)
  expect_match(download[[1L]], "/region_based/.*rankings.feather$")
  expect_identical(readBin(ranking, "raw", n = 64L), charToRaw("ARROW1oldARROW1"))
  expect_length(list.files(directory), 3L)
})

test_that("PrepareDB extracts and reuses a local LINGER archive", {
  tar <- Sys.which("tar")
  skip_if(!nzchar(tar))
  directory <- withr::local_tempdir()
  withr::local_envvar(R_USER_DATA_DIR = file.path(directory, "cache"))
  bulk <- file.path(directory, "data_bulk")
  dir.create(bulk)
  files <- c(
    "TFName.txt", "Match2.txt", "bulk_gene_all.txt", "all_hg19.txt",
    "motifWeight.txt", "RE_gene_corr_hg19.bed", "hg38_hg19_pair.bed",
    "MotifTarget_Matrix_chr1.txt", "MotifTarget_matrix_chr1.bed", "hg38_Peaks_chr1.bed",
    "chr1_gene.txt", "chr1_index.txt", "chr1_index_all.txt", "fisher_chr1.pt", "all_models_chr1.pt"
  )
  for (name in files) writeLines(name, file.path(bulk, name))
  writeLines("unused", file.path(bulk, "unrelated.txt"))
  local <- PrepareDB(db = "LINGER", data_dir = directory, verbose = FALSE)$LINGER
  expect_identical(local$bulk_reference, normalizePath(bulk))
  archive <- file.path(directory, "reference.tar.gz")
  processx::run(tar, c("-czf", archive, "-C", directory, "./data_bulk"))
  extracted <- PrepareDB(db = "LINGER", data_dir = archive, verbose = FALSE)$LINGER
  expect_setequal(list.files(extracted$bulk_reference), files)
  expect_identical(PrepareDB(db = "LINGER", data_dir = archive, verbose = FALSE)$LINGER, extracted)
  broken <- file.path(directory, "broken.tgz")
  writeLines("<html>not an archive</html>", broken)
  expect_error(PrepareDB(db = "LINGER", data_dir = broken, verbose = FALSE), "Could not list")
})
