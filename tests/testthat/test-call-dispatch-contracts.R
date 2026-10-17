make_dispatch_object <- function(nfeatures = 40L, ncells = 32L, seed = 511L, prefix = "c") {
  set.seed(seed)
  counts <- matrix(rpois(nfeatures * ncells, 5) + 1, nrow = nfeatures,
    dimnames = list(paste0("g", seq_len(nfeatures)), paste0(prefix, seq_len(ncells))))
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object <- NormalizeData(object, verbose = FALSE)
  ScaleData(object, features = rownames(object), verbose = FALSE)
}

test_that("Seurat fallbacks preserve distinct canonical command history", {
  object <- make_dispatch_object()
  object@commands <- list()
  normalize <- get("NormalizeData.Seurat", asNamespace("Seurat"))
  scale <- get("ScaleData.Seurat", asNamespace("Seurat"))
  variable <- get("FindVariableFeatures.Seurat", asNamespace("Seurat"))

  actual <- NormalizeData(object, block.size = 50, verbose = FALSE)
  expected <- normalize(object, block.size = 50, verbose = FALSE)
  expect_identical(actual@assays, expected@assays)
  expect_true("NormalizeData.RNA" %in% names(actual@commands))

  actual <- suppressWarnings(ScaleData(actual, nothing = "ignored", verbose = FALSE))
  expected <- suppressWarnings(scale(expected, nothing = "ignored", verbose = FALSE))
  expect_identical(actual@assays, expected@assays)
  expect_true(all(c("NormalizeData.RNA", "ScaleData.RNA") %in% names(actual@commands)))

  actual <- suppressWarnings(FindVariableFeatures(actual, nfeatures = 10, nothing = "ignored", verbose = FALSE))
  expected <- suppressWarnings(variable(expected, nfeatures = 10, nothing = "ignored", verbose = FALSE))
  expect_identical(actual@assays, expected@assays)
  expect_true(all(c("NormalizeData.RNA", "ScaleData.RNA", "FindVariableFeatures.RNA") %in% names(actual@commands)))
  expect_false(any(c(".RNA", ".fn.RNA") %in% names(actual@commands)))
  expect_identical(SeuratObject::Command(actual, "NormalizeData.RNA", value = "normalization.method"), "LogNormalize")
  expect_true(all(vapply(actual@commands, function(x) sum(nchar(x@call.string)) < 2000, logical(1))))
})

test_that("nonlinear ScaleData fallback keeps canonical history and values", {
  object <- make_dispatch_object()
  object@commands <- list()
  scale <- get("ScaleData.Seurat", asNamespace("Seurat"))
  actual <- suppressWarnings(ScaleData(object, model.use = "poisson", verbose = FALSE))
  expected <- suppressWarnings(scale(object, model.use = "poisson", verbose = FALSE))
  expect_identical(actual@assays, expected@assays)
  expect_true("ScaleData.RNA" %in% names(actual@commands))
  expect_identical(SeuratObject::Command(actual, "ScaleData.RNA", value = "model.use"), "poisson")
})

test_that("dynamic PCA keeps numerical results and compact command text", {
  object <- make_dispatch_object()
  testthat::local_mocked_bindings(RunDimsEstimate = function(...) seq_len(5L), .package = "scop")
  actual <- RunDimsReduction(object, features = rownames(object), linear_reduction = "pca",
    linear_reduction_dims = 5, linear_reduction_params = list(npcs = 5), verbose = FALSE)
  expected <- RunPCA(object, features = rownames(object), npcs = 5, slot = "data", seed.use = 11, verbose = FALSE)
  expect_identical(SeuratObject::Embeddings(actual, "pca"), SeuratObject::Embeddings(expected, "pca"))
  expect_identical(SeuratObject::Loadings(actual, "pca"), SeuratObject::Loadings(expected, "pca"))
  command <- tail(actual@commands, 1L)[[1L]]
  expect_lt(sum(nchar(command@call.string)), 2000)
})

test_that("RareQ neighbor dispatch keeps compact commands in both verbosity branches", {
  object <- make_dispatch_object()
  object <- RunPCA(object, features = rownames(object), npcs = 5, verbose = FALSE)
  n <- ncol(object)
  original_get <- get_namespace_fun
  testthat::local_mocked_bindings(
    check_r = function(...) invisible(TRUE),
    get_namespace_fun = function(package, fun) {
      if (identical(package, "RareQ")) {
        return(switch(fun,
          ComputeQ = function(...) rep(0.1, n),
          FindRare = function(...) rep(c("A", "B"), length.out = n),
          stop("Unexpected RareQ backend function")))
      }
      original_get(package, fun)
    }, .package = "scop")
  for (verbose in c(FALSE, TRUE)) {
    actual <- suppressMessages(RunRareQ(object, dims = 1:5, k.param = 3, k = 2, verbose = verbose))
    expect_length(actual@tools$RareQ$clusters, n)
    expect_identical(names(actual@tools$RareQ$clusters), colnames(object))
    command <- tail(actual@commands, 1L)[[1L]]
    expect_lt(sum(nchar(command@call.string)), 2000)
  }
})

test_that("metacell grouping treats metadata names as data", {
  object <- make_dispatch_object(nfeatures = 12L, ncells = 4L)
  object[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(seq_len(8), nrow = 4L,
      dimnames = list(colnames(object), c("PC_1", "PC_2"))), key = "PC_", assay = "RNA")
  counts <- GetAssayData5(object, layer = "counts")
  original_get <- get_namespace_fun
  testthat::local_mocked_bindings(
    check_r = function(...) invisible(TRUE),
    PrepareEnv = function(...) invisible(TRUE),
    scop_python_import = function(...) list(),
    python_anndata_from_matrix = function(x, embeddings) list(obsm = embeddings),
    get_namespace_fun = function(package, fun) {
      if (identical(package, "SuperCell")) {
        return(switch(fun,
          SCimplify = function(ge, ...) list(membership = rep(1L, ncol(ge))),
          supercell_GE = function(ge, membership) metacell_aggregate_counts(ge, membership),
          stop("Unexpected SuperCell backend function")))
      }
      original_get(package, fun)
    }, .package = "scop")
  groups <- data.frame(condition = c("A", "A", "B", "B"), row.names = colnames(object))
  runners <- list(
    function(group_df) metacell_supercell(counts, gamma = 2, group_df = group_df, verbose = FALSE),
    function(group_df) metacell_seacells(object, counts, gamma = 2, group_df = group_df, verbose = FALSE),
    function(group_df) metacell_metacell(counts, gamma = 2, group_df = group_df, verbose = FALSE)
  )
  for (runner in runners) {
    expected <- runner(groups)
    for (name in c("drop", "sep", "lex.order")) {
      renamed <- groups
      names(renamed) <- name
      actual <- runner(renamed)
      expect_identical(actual$membership, expected$membership)
      expect_equal(actual$counts, expected$counts)
    }
  }
})

test_that("LIANA sample columns cannot consume paste options in either backend", {
  native <- ccc_aggregate_liana_table_cpp
  native_results <- list()
  testthat::local_mocked_bindings(ccc_aggregate_liana_table_cpp = function(...) {
    result <- native(...)
    native_results[[length(native_results) + 1L]] <<- result
    result
  }, .package = "scop")
  tab <- data.frame(source = "A", target = "B", ligand_complex = "L", receptor_complex = "R",
    score = rep(c(1, 2), each = 500L), pvalue = 0.1,
    sample_id = rep(c("s1", "s2"), each = 500L), stringsAsFactors = FALSE)
  for (backend in c("r", "cpp")) {
    expected <- ccc_aggregate_liana_table(tab, sample_col = "sample_id", backend = backend)
    for (name in c("sep", "collapse", "recycle0")) {
      renamed <- tab
      names(renamed)[names(renamed) == "sample_id"] <- name
      actual <- ccc_aggregate_liana_table(renamed, sample_col = name, backend = backend)
      names(actual)[names(actual) == name] <- "sample_id"
      expect_equal(actual, expected)
      expect_equal(nrow(actual), 2L)
    }
  }
  expect_length(native_results, 4L)
  expect_true(all(vapply(native_results, function(x) is.data.frame(x) && nrow(x) == 2L, logical(1))))
})

test_that("gene immune correlations preserve each statistical method", {
  set.seed(512)
  x <- matrix(rnorm(80), ncol = 2, dimnames = list(NULL, c("g1", "g2")))
  y <- matrix(rnorm(80), ncol = 2, dimnames = list(NULL, c("i1", "i2")))
  for (method in c("pearson", "spearman", "kendall")) {
    actual <- gene_immune_cor_edges(x, y, method = method)
    for (i in seq_len(nrow(actual))) {
      expected <- stats::cor.test(x[, actual$gene[i]], y[, actual$path[i]], method = method, exact = FALSE)
      expect_identical(actual$r[i], unname(expected$estimate))
      expect_identical(actual$p[i], expected$p.value)
    }
  }
})

test_that("Seurat integration retains canonical anchor and integration commands", {
  first <- make_dispatch_object(nfeatures = 60L, ncells = 55L, seed = 513L, prefix = "a")
  second <- make_dispatch_object(nfeatures = 60L, ncells = 55L, seed = 514L, prefix = "b")
  for (method in c("LogNormalize", "Log")) {
    actual <- suppressWarnings(suppressMessages(Seurat_integrate(
      srt_list = list(first, second), batch = c("first", "second"), do_normalization = FALSE, do_HVF_finding = FALSE,
      HVF = rownames(first), do_scaling = TRUE, linear_reduction_dims = 5,
      linear_reduction_dims_use = 1:5, nonlinear_reduction = NULL, neighbor_k = 5,
      FindIntegrationAnchors_params = list(dims = 1:5, k.anchor = 3, k.filter = NA, k.score = 5),
      IntegrateData_params = list(dims = 1:5, k.weight = 5, normalization.method = method), verbose = FALSE)))
    expect_identical(SeuratObject::Command(actual, "FindIntegrationAnchors")@name, "FindIntegrationAnchors")
    expect_identical(SeuratObject::Command(actual, "IntegrateData")@name, "IntegrateData")
    expect_match(paste(SeuratObject::Command(actual, "IntegrateData")@call.string, collapse = " "), "^IntegrateData\\(")
    expect_false(".fn" %in% names(actual@commands))
    expect_identical(SeuratObject::Command(actual, "IntegrateData", value = "normalization.method"), "LogNormalize")
  }
})

test_that("TFIDF integration retains its canonical anchor command", {
  skip_if_not_installed("Signac")
  peaks <- paste0("chr1-", seq(1001, by = 1000, length.out = 120), "-", seq(1500, by = 1000, length.out = 120))
  make_object <- function(seed, prefix) {
    set.seed(seed)
    counts <- matrix(rbinom(120 * 60, 2, 0.3), nrow = 120,
      dimnames = list(peaks, paste0(prefix, seq_len(60))))
    assay <- Signac::CreateChromatinAssay(Matrix::Matrix(counts, sparse = TRUE), min.cells = 0, min.features = 0)
    SeuratObject::CreateSeuratObject(assay, assay = "ATAC")
  }
  actual <- suppressWarnings(suppressMessages(Seurat_integrate(
    srt_list = list(make_object(519, "a"), make_object(520, "b")), batch = c("one", "two"), assay = "ATAC",
    normalization_method = "TFIDF", do_normalization = TRUE, do_HVF_finding = FALSE,
    HVF = peaks, do_scaling = FALSE, linear_reduction_dims = 5, linear_reduction_dims_use = 2:5,
    nonlinear_reduction = NULL, neighbor_k = 5,
    FindIntegrationAnchors_params = list(dims = 2:5, k.anchor = 3, k.filter = NA, k.score = 5),
    IntegrateEmbeddings_params = list(dims.to.integrate = 1:5, k.weight = 5), verbose = FALSE)))
  expect_identical(SeuratObject::Command(actual, "FindIntegrationAnchors")@name, "FindIntegrationAnchors")
  expect_identical(dim(SeuratObject::Embeddings(actual, "Seuratlsi")), c(120L, 5L))
  expect_false(".fn" %in% names(actual@commands))
})
