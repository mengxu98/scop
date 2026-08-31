make_scenic_parallel_fixture <- function() {
  fixture_dir <- tempfile("scenic-parallel-fixture-")
  dir.create(fixture_dir)
  genes <- c("TF1", paste0("Gene", seq_len(59L)))
  motifs <- paste0("M", seq_len(100L))
  set.seed(42)

  ranking_paths <- file.path(
    fixture_dir,
    c("ranking-a.feather", "ranking-b.feather")
  )
  for (ranking_path in ranking_paths) {
    ranking <- replicate(
      length(genes),
      sample.int(length(genes), length(motifs), replace = TRUE) - 1L
    )
    colnames(ranking) <- genes
    ranking[1L, ] <- 0L
    ranking <- as.data.frame(ranking, check.names = FALSE)
    ranking[["motifs"]] <- motifs
    arrow::write_feather(ranking, ranking_path)
  }

  annotations <- data.frame(
    motif_id = motifs,
    gene_name = c("TF1", rep("OTHER", length(motifs) - 1L)),
    motif_similarity_qvalue = 0,
    orthologous_identity = 1,
    stringsAsFactors = FALSE
  )
  annotation_path <- file.path(fixture_dir, "motifs.tsv")
  utils::write.table(
    annotations,
    annotation_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  adjacency <- data.frame(
    TF = "TF1",
    target = genes[-1L],
    importance = rev(seq_along(genes[-1L])),
    stringsAsFactors = FALSE
  )
  list(
    adjacency = adjacency,
    ranking_paths = ranking_paths,
    annotation_path = annotation_path
  )
}

run_scenic_parallel_fixture <- function(fixture, cores, parallel_backend) {
  suppressWarnings(getFromNamespace("cistarget2", "scop")(
    adjacency = fixture$adjacency,
    ranking_dbs = fixture$ranking_paths,
    motif_annotations = fixture$annotation_path,
    min_regulon_size = 5L,
    rank_threshold = 50L,
    nes_threshold = 2,
    auc_threshold = 0.5,
    cores = cores,
    parallel_backend = parallel_backend,
    diagnostics = TRUE,
    verbose = FALSE
  ))
}

scenic_regulon_payload <- function(x) {
  attr(x, "profile_sec") <- NULL
  attr(x, "diagnostics") <- NULL
  x
}

test_that("RunSCENIC exposes and forwards the native cisTarget parallel backend", {
  expect_identical(
    eval(formals(RunSCENIC)$parallel_backend)[[1L]],
    "auto"
  )

  testthat::local_mocked_bindings(
    .package = "scop",
    scenic_cpp = function(...) {
      args <- list(...)
      args[c("cores", "parallel_backend")]
    }
  )
  result <- RunSCENIC(
    structure(list(), class = "Seurat"),
    backend = "cpp",
    cores = 3L,
    parallel_backend = "psock"
  )
  expect_identical(result$cores, 3L)
  expect_identical(result$parallel_backend, "psock")
})

test_that("native cisTarget parallel backend validates Windows fork explicitly", {
  testthat::local_mocked_bindings(
    .package = "base",
    .Platform = list(OS.type = "windows")
  )
  expect_error(
    getFromNamespace("cistarget2", "scop")(
      adjacency = NULL,
      ranking_dbs = character(),
      motif_annotations = "",
      parallel_backend = "fork"
    ),
    "fork parallel backend is unavailable on Windows"
  )
})

test_that("cisTarget worker closure contains only minimal state", {
  skip_if_not_installed("arrow")
  fixture <- make_scenic_parallel_fixture()
  captured <- new.env(parent = emptyenv())
  original_parallelize_fun <- thisutils::parallelize_fun
  testthat::local_mocked_bindings(
    .package = "thisutils",
    parallelize_fun = function(x, fun, ..., backend) {
      captured$backend <- backend
      captured$environment_names <- ls(environment(fun), all.names = TRUE)
      captured$serialized_bytes <- length(serialize(fun, NULL))
      original_parallelize_fun(x, fun, ..., backend = backend)
    }
  )

  result <- run_scenic_parallel_fixture(
    fixture,
    cores = 1L,
    parallel_backend = "auto"
  )
  expect_identical(captured$backend, "auto")
  expect_identical(
    captured$environment_names,
    "state"
  )
  expect_gt(captured$serialized_bytes, 0L)
  expect_identical(names(result), "TF1(+)")
})

test_that("serial, PSOCK, and fork cisTarget results have identical payloads", {
  skip_if_not_installed("arrow")
  fixture <- make_scenic_parallel_fixture()
  serial <- run_scenic_parallel_fixture(fixture, 1L, "psock")
  psock <- run_scenic_parallel_fixture(fixture, 2L, "psock")
  results <- list(serial = serial, psock = psock)
  if (.Platform$OS.type != "windows") {
    results$fork <- run_scenic_parallel_fixture(fixture, 2L, "fork")
  }

  payloads <- lapply(results, scenic_regulon_payload)
  expect_true(all(vapply(
    payloads[-1L],
    identical,
    logical(1L),
    payloads[[1L]]
  )))
  diagnostics <- lapply(results, attr, "diagnostics")
  expect_true(all(vapply(
    diagnostics[-1L],
    identical,
    logical(1L),
    diagnostics[[1L]]
  )))
  profile_names <- lapply(results, function(x) names(attr(x, "profile_sec")))
  expect_true(all(vapply(
    profile_names[-1L],
    identical,
    logical(1L),
    profile_names[[1L]]
  )))
})
