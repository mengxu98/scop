test_that("RunCisTarget exposes the native cisTarget backend", {
  work_dir <- tempfile("cistarget_cpp_")
  dir.create(work_dir)
  ranking <- file.path(work_dir, "ranking.feather")
  motif <- file.path(work_dir, "motifs.tsv")
  file.create(ranking, motif)
  adj <- data.frame(
    TF = rep(c("TF1", "TF2"), each = 3),
    target = paste0("Gene", 1:6),
    stringsAsFactors = FALSE
  )

  testthat::local_mocked_bindings(
    .package = "scop",
    scenic_reference = function(...) {
      list(ranking_dbs = ranking, motif_annotations = motif)
    },
    cistarget2 = function(adjacency,
                          expr_mtx,
                          ranking_dbs,
                          motif_annotations,
                          min_regulon_size,
                          cores,
                          parallel_backend,
                          verbose) {
      expect_true("importance" %in% colnames(adjacency))
      expect_null(expr_mtx)
      expect_identical(ranking_dbs, ranking)
      expect_identical(motif_annotations, motif)
      expect_identical(min_regulon_size, 2L)
      expect_identical(parallel_backend, "psock")
      list("TF1(+)" = c("Gene1", "Gene2", "Gene3"))
    }
  )

  result <- RunCisTarget(
    adj,
    backend = "cpp",
    ranking_dbs = ranking,
    motif_annotations = motif,
    min_regulon_size = 2,
    parallel_backend = "psock",
    work_dir = work_dir,
    verbose = FALSE
  )
  expect_identical(names(result$regulons), "TF1(+)")
  expect_true(all(file.exists(unlist(result[-1]))))
  ctx <- utils::read.csv(result$ctx_file)
  expect_identical(ctx$target, c("Gene1", "Gene2", "Gene3"))
})

test_that("RunCisTarget forwards cores to the RcisTarget backend", {
  work_dir <- tempfile("cistarget_r_")
  dir.create(work_dir)
  ranking <- file.path(work_dir, "ranking.feather")
  motif <- file.path(work_dir, "motifs.tsv")
  file.create(ranking)
  writeLines("motif\tTF", motif)
  adj <- data.frame(
    TF = rep("TF1", 3),
    target = paste0("Gene", 1:3),
    importance = 1,
    stringsAsFactors = FALSE
  )

  testthat::local_mocked_bindings(
    .package = "scop",
    check_r = function(...) invisible(TRUE),
    # scenic_reference() downloads the species TF list when regulators
    # are missing; mock it so the test never depends on network access.
    scenic_reference = function(...) {
      list(ranking_dbs = ranking, motif_annotations = motif)
    },
    get_namespace_fun = function(package, fun) {
      expect_identical(package, "RcisTarget")
      if (identical(fun, "importRankings")) {
        return(function(path) {
          expect_identical(path, ranking)
          structure(list(), class = "mock_rankings")
        })
      }
      expect_identical(fun, "cisTarget")
      function(gene_sets, rankings, motifAnnot, nCores, ...) {
        expect_identical(nCores, 1L)
        expect_s3_class(motifAnnot, "data.table")
        data.frame(
          geneSet = "TF1",
          NES = 4,
          enrichedGenes = "Gene1;Gene2;Gene3",
          stringsAsFactors = FALSE
        )
      }
    }
  )

  result <- RunCisTarget(
    adj,
    backend = "r",
    ranking_dbs = ranking,
    motif_annotations = motif,
    min_regulon_size = 2,
    cores = 1,
    work_dir = work_dir,
    verbose = FALSE
  )
  expect_identical(names(result$regulons), "TF1(+)")
  expect_identical(result$regulons[[1L]], paste0("Gene", 1:3))
})
