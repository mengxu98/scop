test_that("FindAllMarkers probes Presto only after support checks and before marker context", {
  native_all_markers <- get("FindAllMarkers.Seurat", asNamespace("scop"))
  fallback <- data.frame(path = "Seurat", stringsAsFactors = FALSE)
  events <- character()

  testthat::local_mocked_bindings(
    marker_assay_is_chromatin = function(...) {
      events <<- c(events, "chromatin")
      FALSE
    },
    marker_all_supported = function(test.use, ...) {
      events <<- c(events, "supported")
      identical(test.use, "wilcox")
    },
    check_r = function(packages, install = TRUE, verbose = TRUE) {
      if (identical(packages, "presto")) {
        events <<- c(events, "presto")
        return(c(presto = FALSE))
      }
      TRUE
    },
    marker_context = function(...) {
      events <<- c(events, "context")
      stop("marker context should not be materialized without Presto")
    },
    get_namespace_fun = function(package, fun) {
      expect_identical(package, "Seurat")
      expect_identical(fun, "FindAllMarkers")
      function(...) fallback
    },
    .package = "scop"
  )

  supported <- native_all_markers(
    object = NULL,
    test.use = "wilcox",
    verbose = FALSE
  )
  expect_identical(supported, fallback)
  expect_identical(events, c("chromatin", "supported", "presto"))

  events <- character()
  unsupported <- native_all_markers(
    object = NULL,
    test.use = "roc",
    verbose = FALSE
  )
  expect_identical(unsupported, fallback)
  expect_identical(events, c("chromatin", "supported"))
})

test_that("RunCellChat requests Presto only for a supported fast path", {
  installs <- logical()
  testthat::local_mocked_bindings(
    check_r = function(packages, install = TRUE, verbose = TRUE) {
      if (identical(packages, "immunogenomics/presto")) {
        installs <<- c(installs, install)
        return(c(presto = TRUE))
      }
      TRUE
    },
    get_namespace_fun = function(package, fun) {
      expect_identical(package, "CellChat")
      expect_identical(fun, "identifyOverExpressedGenes")
      function(object, do.fast = FALSE) object
    },
    validate_cc_input = function(...) stop("input validation reached"),
    .package = "scop"
  )

  expect_error(
    RunCellChat(NULL, group.by = "group", do.fast = FALSE, verbose = FALSE),
    "input validation reached"
  )
  expect_length(installs, 0L)
  expect_error(
    RunCellChat(NULL, group.by = "group", do.fast = TRUE, verbose = FALSE),
    "input validation reached"
  )
  expect_identical(installs, TRUE)
})

test_that("RunCellChat fails before its fast path when Presto installation fails", {
  testthat::local_mocked_bindings(
    check_r = function(packages, install = TRUE, verbose = TRUE) {
      if (identical(packages, "immunogenomics/presto")) {
        expect_true(install)
        return(c(presto = FALSE))
      }
      TRUE
    },
    get_namespace_fun = function(package, fun) {
      expect_identical(package, "CellChat")
      function(object, do.fast = FALSE) object
    },
    validate_cc_input = function(...) stop("input validation reached"),
    .package = "scop"
  )

  expect_error(
    RunCellChat(NULL, group.by = "group", do.fast = TRUE, verbose = FALSE),
    "optional.*presto.*do.fast.*unavailable"
  )
})

test_that("RunCellChat skips the Presto request when the backend has no fast path", {
  installs <- logical()
  testthat::local_mocked_bindings(
    check_r = function(packages, install = TRUE, verbose = TRUE) {
      if (identical(packages, "immunogenomics/presto")) {
        installs <<- c(installs, install)
        return(c(presto = TRUE))
      }
      TRUE
    },
    get_namespace_fun = function(package, fun) {
      expect_identical(package, "CellChat")
      expect_identical(fun, "identifyOverExpressedGenes")
      function(object) object
    },
    validate_cc_input = function(...) stop("input validation reached"),
    .package = "scop"
  )

  expect_error(
    RunCellChat(NULL, group.by = "group", do.fast = TRUE, verbose = FALSE),
    "input validation reached"
  )
  expect_length(installs, 0L)
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
    check_r = function(packages, install = TRUE, verbose = TRUE) {
      if (identical(packages, "presto")) {
        requests[[length(requests) + 1L]] <<- list(
          packages = packages,
          install = install
        )
        return(c(presto = FALSE))
      }
      TRUE
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
  expect_true(all(vapply(requests, function(x) identical(x$packages, "presto"), logical(1))))
  expect_true(all(!vapply(requests, `[[`, logical(1), "install")))
})

test_that("the installed runtime Presto backend resolves and executes its current API", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("presto")

  presto_fun <- get_namespace_fun("presto", "wilcoxauc")
  skip_if(!is.function(presto_fun), "The runtime-optional Presto backend is unavailable")

  expression <- Matrix::Matrix(
    matrix(
      c(
        8, 7, 9, 1, 0, 2,
        0, 1, 0, 7, 8, 9,
        2, 2, 3, 2, 3, 2
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(paste0("g", 1:3), paste0("c", 1:6))
    ),
    sparse = TRUE
  )
  result <- presto_fun(
    X = expression,
    y = factor(rep(c("A", "B"), each = 3)),
    verbose = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_setequal(unique(as.character(result$group)), c("A", "B"))
  expect_true(all(c("feature", "group", "pval", "padj") %in% colnames(result)))
  expect_setequal(unique(result$feature), rownames(expression))
})

test_that("the installed runtime Presto backend drives SCOP marker entry points", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("presto")

  presto_fun <- get_namespace_fun("presto", "wilcoxauc")
  skip_if(!is.function(presto_fun), "The runtime-optional Presto backend is unavailable")

  set.seed(16)
  counts <- Matrix::rsparsematrix(
    40,
    20,
    density = 0.25,
    rand.x = function(n) stats::rpois(n, 2) + 1
  )
  dimnames(counts) <- list(paste0("g", 1:40), paste0("c", 1:20))
  object <- Seurat::CreateSeuratObject(counts)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  SeuratObject::Idents(object) <- factor(rep(c("A", "B"), each = 10))

  testthat::local_mocked_bindings(
    FindMarkers.Seurat = function(...) stop("Seurat pairwise fallback was reached"),
    FindAllMarkers = function(...) stop("Seurat all-markers fallback was reached"),
    .package = "Seurat"
  )

  pair <- scop::FindMarkers(
    object,
    cells.1 = colnames(object)[1:10],
    cells.2 = colnames(object)[11:20],
    logfc.threshold = 0,
    min.pct = 0,
    verbose = FALSE
  )
  all_markers <- scop::FindAllMarkers(
    object,
    logfc.threshold = 0,
    min.pct = 0,
    return.thresh = 1,
    verbose = FALSE
  )

  expect_s3_class(pair, "data.frame")
  expect_s3_class(all_markers, "data.frame")
  expect_true(all(c("p_val", "p_val_adj") %in% colnames(pair)))
  expect_true(all(c("gene", "cluster", "p_val", "p_val_adj") %in% colnames(all_markers)))
})

test_that("optional GitHub backends stay out of the package metadata", {
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
  expect_false(any(grepl("presto", remotes, fixed = TRUE)))
  expect_setequal(remotes, c("mengxu98/thisplot", "mengxu98/thisutils"))
})
