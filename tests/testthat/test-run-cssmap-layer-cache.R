make_cssmap_cache_objects <- function() {
  counts <- methods::as(Matrix::Matrix(
    matrix(
      c(5, 0, 2, 1, 3, 0, 4, 1, 2, 1, 3, 2),
      nrow = 4,
      dimnames = list(paste0("Gene", 1:4), paste0("Cell", 1:3))
    ),
    sparse = TRUE
  ), "dgCMatrix")
  query <- Seurat::NormalizeData(Seurat::CreateSeuratObject(counts), verbose = FALSE)
  reference <- Seurat::NormalizeData(Seurat::CreateSeuratObject(counts), verbose = FALSE)
  css <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(seq_len(6), ncol = 2, dimnames = list(colnames(reference), c("CSS_1", "CSS_2"))),
    assay = "RNA",
    key = "CSS_"
  )
  css@misc$model <- list(sim2profiles = list())
  reference[["css"]] <- css
  reference[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(seq_len(6), ncol = 2, dimnames = list(colnames(reference), c("UMAP_1", "UMAP_2"))),
    assay = "RNA",
    key = "UMAP_"
  )
  list(query = query, reference = reference)
}

test_that("RunCSSMap reuses the checked data layers for projection", {
  skip_if_not_installed("simspec")
  objects <- make_cssmap_cache_objects()
  query_data <- Matrix::Matrix(matrix(seq_len(12), nrow = 4), sparse = TRUE)
  ref_data <- Matrix::Matrix(matrix(seq_len(12) + 10, nrow = 4), sparse = TRUE)
  rownames(query_data) <- rownames(ref_data) <- paste0("Gene", 1:4)
  colnames(query_data) <- colnames(objects$query)
  colnames(ref_data) <- colnames(objects$reference)
  calls <- character()

  testthat::local_mocked_bindings(
    css_project = function(object, model) {
      expect_identical(object, query_data)
      expect_identical(model, objects$reference[["css"]]@misc$model)
      matrix(seq_len(6), ncol = 2)
    },
    .package = "simspec"
  )
  testthat::local_mocked_bindings(
    check_r = function(...) TRUE,
    GetAssayData5 = function(object, layer, assay, ...) {
      calls <<- c(calls, paste(if (identical(object, objects$query)) "query" else "reference", layer, sep = ":"))
      if (identical(object, objects$query)) query_data else ref_data
    },
    CheckDataType = function(...) "lognormalized",
    RunKNNMap = function(srt_query, ...) srt_query,
    .package = "scop"
  )

  out <- RunCSSMap(
    objects$query,
    objects$reference,
    ref_css = "css",
    ref_umap = "umap",
    projection_method = "knn",
    verbose = FALSE
  )

  expect_identical(calls, c("query:data", "reference:data"))
  expect_true("css_proj" %in% SeuratObject::Reductions(out))
})

test_that("RunCSSMap preserves the real simspec projection", {
  skip_if_not_installed("simspec")
  objects <- make_cssmap_cache_objects()
  query_data <- GetAssayData5(objects$query, layer = "data")
  profiles <- matrix(c(1, 4, 2, 3, 4, 1, 3, 2),
    nrow = 4,
    dimnames = list(rownames(query_data), c("p1", "p2"))
  )
  model <- list(model = list(
    spectrum_type = "corr_raw", profiles = list(batch = profiles),
    args = c(corr_method = "pearson")
  ))
  objects$reference[["css"]]@misc$model <- model
  expected <- invoke_fun(simspec::css_project, list(object = query_data, model = model))
  rownames(expected) <- colnames(objects$query)
  colnames(expected) <- paste0("CSSPROJ_", seq_len(ncol(expected)))
  testthat::local_mocked_bindings(
    check_r = function(...) invisible(TRUE),
    RunKNNMap = function(srt_query, ...) srt_query,
    .package = "scop"
  )
  actual <- RunCSSMap(objects$query, objects$reference,
    ref_css = "css", ref_umap = "umap",
    projection_method = "knn", verbose = FALSE
  )
  expect_true(all(is.finite(expected)))
  expect_identical(SeuratObject::Embeddings(actual, "css_proj"), expected)
})

test_that("CSS integration preserves the Seurat SNN model", {
  skip_if_not_installed("simspec")
  set.seed(11)
  counts <- matrix(stats::rpois(60 * 100, lambda = 1),
    nrow = 60,
    dimnames = list(paste0("Gene", seq_len(60)), paste0("Cell", seq_len(100)))
  )
  object <- Seurat::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object <- Seurat::ScaleData(object, features = rownames(object), verbose = FALSE)
  object$batch <- rep(c("A", "B"), each = 50)
  original_graph <- getFromNamespace("build_knn_graph", "simspec")
  testthat::local_mocked_bindings(
    check_r = function(...) invisible(TRUE),
    find_neighbors_and_clusters = function(srt, ...) srt,
    run_nonlinear_reduction = function(srt, ...) srt,
    .package = "scop"
  )
  actual <- CSS_integrate(
    srt_merge = object, batch = "batch", do_HVF_finding = FALSE,
    HVF = rownames(object), linear_reduction_dims = 5,
    linear_reduction_dims_use = seq_len(5),
    CSS_params = list(k = 5, min_cluster_num = 1), verbose = FALSE
  )
  expected <- actual
  graph <- environment(css_simspec_method())$build_knn_graph
  expect_identical(body(graph), body(original_graph))
  clusters <- stats::setNames(rep(NA_character_, ncol(actual)), colnames(actual))
  for (batch in c("A", "B")) {
    cells <- colnames(actual)[actual$batch == batch]
    dr <- SeuratObject::Embeddings(actual, "CSSpca")[cells, , drop = FALSE]
    nn <- RANN::nn2(dr, dr, k = 6)$nn.idx[, -1, drop = FALSE]
    compute_snn <- getFromNamespace("ComputeSNN", "Seurat")
    args <- list(nn_ranked = nn, prune = 1 / 15)
    if ("nthreads" %in% names(formals(compute_snn))) args$nthreads <- 1L
    snn <- do.call(compute_snn, args)
    dimnames(snn) <- list(cells, cells)
    expect_equal(graph(t(dr), k = 5, verbose = FALSE), snn, tolerance = 1e-15)
    clusters[cells] <- as.character(Seurat::FindClusters(
      SeuratObject::as.Graph(snn),
      resolution = 0.6, verbose = FALSE
    )[, 1])
  }
  expected$css_test_cluster <- clusters[colnames(expected)]
  expected <- simspec::cluster_sim_spectrum(
    expected,
    label_tag = "batch", cluster_col = "css_test_cluster",
    var_genes = rownames(object), use_dr = "CSSpca", dims_use = seq_len(5),
    k = 5, min_cluster_num = 1, reduction.name = "CSS", reduction.key = "CSS_",
    verbose = FALSE
  )
  expect_equal(SeuratObject::Embeddings(actual, "CSS"),
    SeuratObject::Embeddings(expected, "CSS"),
    tolerance = 1e-12
  )
  expect_equal(actual[["CSS"]]@misc$model, expected[["CSS"]]@misc$model, tolerance = 1e-12)
  expect_identical(colnames(actual), colnames(object))
  expect_identical(getFromNamespace("build_knn_graph", "simspec"), original_graph)
})
