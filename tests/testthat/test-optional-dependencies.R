test_that("loading scop does not eagerly add optional backend namespaces", {
  skip_if_not_installed("callr")
  loaded <- callr::r(function() {
    suppressPackageStartupMessages(loadNamespace("thisplot"))
    suppressPackageStartupMessages(loadNamespace("thisutils"))
    before <- loadedNamespaces()
    suppressPackageStartupMessages(loadNamespace("scop"))
    setdiff(loadedNamespaces(), before)
  }, libpath = .libPaths())
  expect_false("Signac" %in% loaded)
  expect_false("ComplexHeatmap" %in% loaded)
  expect_false("ggforce" %in% loaded)
})

test_that("optional backends are checked before backend resolution", {
  requested <- character()
  local_mocked_bindings(
    check_r = function(packages, verbose) {
      expect_false(verbose)
      requested <<- c(requested, packages)
      stop("Dependency checkpoint")
    },
    get_namespace_fun = function(...) stop("Backend lookup before dependency check"),
    .package = "scop"
  )
  expect_error(ReadATACData(), "Dependency checkpoint")
  expect_error(getFromNamespace("cnv_plot_heatmap", "scop")(NULL), "Dependency checkpoint")
  expect_error(FerrisWheelPlot(), "Dependency checkpoint")
  expect_identical(requested, c("Signac", "ComplexHeatmap", "ggforce"))
})

test_that("RNA normalization does not require the optional ATAC backend", {
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(matrix(
    c(1L, 3L, 2L, 4L, 2L, 1L),
    nrow = 3L,
    dimnames = list(paste0("g", 1:3), c("c1", "c2"))
  ), sparse = TRUE))
  local_mocked_bindings(
    check_r = function(...) stop("Unexpected optional dependency check"),
    .package = "scop"
  )
  normalized <- suppressWarnings(getFromNamespace("RunDEtestNormalizeDataLayer", "scop")(
    object, "RNA", "Missing normalized layer.", FALSE
  ))
  expected <- log1p(sweep(as.matrix(SeuratObject::LayerData(object, layer = "counts")), 2, c(6, 7), "/") * 10000)
  expect_equal(as.matrix(SeuratObject::LayerData(normalized, layer = "data")), expected)
})

test_that("optional SVD dispatch preserves Signac embeddings", {
  skip_if_not_installed("Signac")
  set.seed(17)
  counts <- Matrix::Matrix(matrix(
    rpois(40L * 20L, 3),
    nrow = 40L,
    dimnames = list(
      paste0("chr1-", seq_len(40L) * 200L, "-", seq_len(40L) * 200L + 100L),
      paste0("c", seq_len(20L))
    )
  ), sparse = TRUE)
  create_assay <- get_namespace_fun("Signac", "CreateChromatinAssay")
  object <- SeuratObject::CreateSeuratObject(create_assay(counts), assay = "peaks")
  object <- get_namespace_fun("Signac", "RunTFIDF")(object, verbose = FALSE)
  features <- rownames(object)
  set.seed(11)
  expected <- get_namespace_fun("Signac", "RunSVD")(
    object, features = features, n = 3, reduction.name = "svd",
    reduction.key = "LSI_", verbose = FALSE, seed.use = 11
  )
  observed <- RunDimsReduction(
    object,
    features = features, linear_reduction = "svd",
    linear_reduction_dims = 3, verbose = FALSE, seed = 11
  )
  expect_equal(
    SeuratObject::Embeddings(observed, "svd"),
    SeuratObject::Embeddings(expected, "svd"),
    tolerance = 1e-10
  )
})

test_that("optional heatmap resolution preserves the plotted matrix", {
  skip_if_not_installed("ComplexHeatmap")
  df <- data.frame(
    cell_type = factor(c("B", "A", "B", "A"), levels = c("A", "B")),
    sample = factor(c("s2", "s2", "s1", "s1"), levels = c("s1", "s2")),
    proportion = c(0.8, 0.2, 0.6, 0.4)
  )
  ht <- getFromNamespace("deconvolution_heatmap", "scop")(NULL, df)
  expect_s4_class(ht, "Heatmap")
  expect_equal(
    ht@matrix,
    matrix(c(0.4, 0.6, 0.2, 0.8), 2L,
      dimnames = list(cell_type = c("A", "B"), sample = c("s1", "s2"))
    )
  )
})
