test_that("Presto checks existing installs without installing and uses the repository for installation", {
  requests <- list()
  testthat::local_mocked_bindings(
    check_r = function(packages, install = TRUE, verbose = TRUE) {
      requests[[length(requests) + 1L]] <<- list(
        packages = packages,
        install = install,
        verbose = verbose
      )
      c(presto = TRUE)
    },
    .package = "scop"
  )

  expect_true(presto_check_r(install = FALSE))
  expect_true(presto_check_r(install = TRUE))
  expect_identical(requests[[1]]$packages, "presto")
  expect_false(requests[[1]]$install)
  expect_false(requests[[1]]$verbose)
  expect_identical(requests[[2]]$packages, "immunogenomics/presto")
  expect_true(requests[[2]]$install)
  expect_false(requests[[2]]$verbose)
})

test_that("presto_get_fun checks and resolves the function exactly once", {
  calls <- c(check = 0L, resolve = 0L)
  wilcoxauc <- function(...) list(ok = TRUE)
  testthat::local_mocked_bindings(
    presto_check_r = function(install = FALSE) {
      calls[["check"]] <<- calls[["check"]] + 1L
      expect_false(install)
      TRUE
    },
    get_namespace_fun = function(package, fun) {
      calls[["resolve"]] <<- calls[["resolve"]] + 1L
      expect_identical(package, "presto")
      expect_identical(fun, "wilcoxauc")
      wilcoxauc
    },
    .package = "scop"
  )

  resolved <- presto_get_fun(
    install = FALSE,
    error_on_missing = FALSE
  )
  expect_identical(resolved, wilcoxauc)
  expect_equal(calls, c(check = 1L, resolve = 1L))
})

test_that("presto_get_fun supports a non-installing optional fallback", {
  resolved <- FALSE
  testthat::local_mocked_bindings(
    presto_check_r = function(install = FALSE) {
      expect_false(install)
      FALSE
    },
    get_namespace_fun = function(...) {
      resolved <<- TRUE
      stop("resolver should not run")
    },
    .package = "scop"
  )

  expect_null(presto_get_fun(
    install = FALSE,
    error_on_missing = FALSE
  ))
  expect_false(resolved)
  expect_error(
    presto_get_fun(install = FALSE),
    "optional.*presto.*unavailable"
  )
})

test_that("presto_get_fun treats a missing backend symbol as unavailable", {
  testthat::local_mocked_bindings(
    presto_check_r = function(...) TRUE,
    get_namespace_fun = function(...) stop("missing wilcoxauc"),
    .package = "scop"
  )

  expect_null(presto_get_fun(error_on_missing = FALSE))
  expect_error(presto_get_fun(), "does not provide.*wilcoxauc")
})

test_that("CellChat requests Presto only for a supported fast path", {
  installs <- logical()
  testthat::local_mocked_bindings(
    presto_check_r = function(install = FALSE) {
      installs <<- c(installs, install)
      TRUE
    },
    .package = "scop"
  )

  expect_false(cellchat_check_presto(FALSE, c("object", "do.fast")))
  expect_false(cellchat_check_presto(TRUE, "object"))
  expect_true(cellchat_check_presto(TRUE, c("object", "do.fast")))
  expect_true(cellchat_check_presto(TRUE, c("object", "...")))
  expect_equal(installs, c(TRUE, TRUE))
})

test_that("CellChat fails before its fast path when Presto installation fails", {
  testthat::local_mocked_bindings(
    presto_check_r = function(install = FALSE) {
      expect_true(install)
      FALSE
    },
    .package = "scop"
  )

  expect_error(
    cellchat_check_presto(TRUE, c("object", "do.fast")),
    "optional.*presto.*do.fast.*unavailable"
  )
})

test_that("RunCellChat routes its supported do.fast path through the Presto helper", {
  presto_installs <- logical()
  testthat::local_mocked_bindings(
    check_r = function(...) invisible(TRUE),
    get_namespace_fun = function(package, fun) {
      expect_identical(package, "CellChat")
      expect_identical(fun, "identifyOverExpressedGenes")
      function(object, do.fast = FALSE) object
    },
    presto_check_r = function(install = FALSE) {
      presto_installs <<- c(presto_installs, install)
      FALSE
    },
    validate_cc_input = function(...) stop("input validation reached"),
    .package = "scop"
  )

  expect_error(
    RunCellChat(NULL, group.by = "group", do.fast = FALSE, verbose = FALSE),
    "input validation reached"
  )
  expect_length(presto_installs, 0L)
  expect_error(
    RunCellChat(NULL, group.by = "group", do.fast = TRUE, verbose = FALSE),
    "optional.*presto.*do.fast.*unavailable"
  )
  expect_identical(presto_installs, TRUE)
})

test_that("automatic marker fast paths fall back without installing Presto", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  set.seed(6)
  counts <- Matrix::rsparsematrix(
    30,
    16,
    density = 0.2,
    rand.x = function(n) stats::rpois(n, 2) + 1
  )
  dimnames(counts) <- list(paste0("g", 1:30), paste0("c", 1:16))
  srt <- Seurat::CreateSeuratObject(counts)
  srt <- get("NormalizeData.Seurat", asNamespace("Seurat"))(
    srt,
    verbose = FALSE
  )
  SeuratObject::Idents(srt) <- factor(rep(c("A", "B"), each = 8))

  pair_fallback <- data.frame(path = "pair")
  all_fallback <- data.frame(path = "all")
  requests <- list()
  testthat::local_mocked_bindings(
    presto_get_fun = function(
      fun = "wilcoxauc",
      install = FALSE,
      error_on_missing = TRUE
    ) {
      requests[[length(requests) + 1L]] <<- list(
        fun = fun,
        install = install,
        error_on_missing = error_on_missing
      )
      NULL
    },
    get_namespace_fun = function(package, fun) {
      expect_identical(package, "Seurat")
      switch(fun,
        FindMarkers.Seurat = function(...) pair_fallback,
        FindAllMarkers = function(...) all_fallback,
        stop("unexpected Seurat fallback")
      )
    },
    .package = "scop"
  )

  pair <- scop::FindMarkers(
    srt,
    cells.1 = colnames(srt)[1:8],
    cells.2 = colnames(srt)[9:16],
    logfc.threshold = 0,
    min.pct = 0,
    verbose = FALSE
  )
  expect_length(requests, 1L)
  all_markers <- scop::FindAllMarkers(
    srt,
    logfc.threshold = 0,
    min.pct = 0,
    return.thresh = 1,
    verbose = FALSE
  )
  expect_length(requests, 2L)

  expect_identical(pair, pair_fallback)
  expect_identical(all_markers, all_fallback)
  expect_true(all(vapply(requests, function(x) identical(x$fun, "wilcoxauc"), logical(1))))
  expect_true(all(!vapply(requests, `[[`, logical(1), "install")))
  expect_true(all(!vapply(requests, `[[`, logical(1), "error_on_missing")))
})

test_that("Presto remains declared only as a runtime remote", {
  description_path <- system.file("DESCRIPTION", package = "scop")
  expect_true(nzchar(description_path))
  description <- read.dcf(description_path)[1, ]
  dependency_fields <- intersect(
    c("Depends", "Imports", "Suggests"),
    names(description)
  )
  dependencies <- unlist(strsplit(
    paste(description[dependency_fields], collapse = ","),
    ",",
    fixed = TRUE
  ))
  dependencies <- trimws(sub("\\s*\\(.*$", "", dependencies))
  remotes <- trimws(unlist(strsplit(description[["Remotes"]], ",", fixed = TRUE)))

  expect_false("presto" %in% dependencies)
  expect_true("immunogenomics/presto" %in% remotes)
})
