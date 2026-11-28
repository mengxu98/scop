test_that("RunKNNPredict probability maxima retain row labels and missing values", {
  probabilities <- matrix(
    c(0.1, 0.7, NA_real_, 0.5, 0.4, NA_real_),
    nrow = 3,
    dimnames = list(c("Q1", "Q2", "Q3"), c("A", "B"))
  )
  testthat::local_mocked_bindings(check_r = function(...) TRUE, .package = "scop")

  expect_identical(
    knn_match_prob_max(probabilities),
    apply(probabilities, 1, max)
  )
})

test_that("matrix neighbor dispatch remains available to KNN reference mapping", {
  set.seed(17)
  reference <- matrix(
    stats::rnorm(80),
    nrow = 20,
    dimnames = list(paste0("ref", seq_len(20)), paste0("dim", seq_len(4)))
  )
  query <- matrix(
    stats::rnorm(40),
    nrow = 10,
    dimnames = list(paste0("query", seq_len(10)), paste0("dim", seq_len(4)))
  )

  observed <- Seurat::FindNeighbors(
    object = reference,
    query = query,
    k.param = 3,
    nn.method = "annoy",
    annoy.metric = "cosine",
    return.neighbor = TRUE,
    verbose = FALSE
  )

  expect_s4_class(observed, "Neighbor")
  expect_identical(dim(observed@nn.idx), c(10L, 3L))
  expect_identical(observed@cell.names, rownames(query))
  expect_true(all(observed@nn.idx >= 1L & observed@nn.idx <= nrow(reference)))
})

test_that("KNN reference mapping and prediction preserve a single selected dimension", {
  make_object <- function(cells, pca) {
    counts <- Matrix::Matrix(
      matrix(1, nrow = 3L, ncol = length(cells),
             dimnames = list(c("g1", "g2", "g3"), cells)),
      sparse = TRUE
    )
    object <- SeuratObject::CreateSeuratObject(counts)
    dimnames(pca) <- list(cells, c("PC_1", "PC_2"))
    object[["pca"]] <- SeuratObject::CreateDimReducObject(
      embeddings = pca, assay = "RNA", key = "PC_"
    )
    object
  }
  reference_cells <- c("r4", "r1", "r2", "r3")
  query_cells <- c("q3", "q1", "q2")
  reference <- make_object(
    reference_cells, cbind(c(-4, 0, 10, 30), c(900, 100, -800, 8))
  )
  reference$kind <- c("A", "B", "C", "D")
  reference_umap <- cbind(UMAP_1 = 1:4, UMAP_2 = c(8, 3, 9, 6))
  rownames(reference_umap) <- reference_cells
  reference[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = reference_umap, assay = "RNA", key = "UMAP_"
  )
  query <- make_object(query_cells, cbind(c(2, 9, 28), c(870, 5, -750)))
  args <- list(
    srt_query = query, srt_ref = reference, ref_group = "kind",
    query_reduction = "pca", ref_reduction = "pca",
    query_dims = 1L, ref_dims = 1L, distance_metric = "euclidean",
    k = 1L, nn_method = "raw", verbose = FALSE
  )

  mapped <- do.call(
    RunKNNMap, c(args, list(ref_umap = "umap", projection_method = "knn"))
  )
  predicted <- do.call(RunKNNPredict, args)

  expect_identical(colnames(mapped), query_cells)
  expect_identical(colnames(predicted), query_cells)
  expect_identical(as.character(mapped$predicted_ref_group), c("B", "C", "D"))
  expect_identical(as.character(predicted$KNNPredict_classification), c("B", "C", "D"))
  expect_equal(unname(predicted$KNNPredict_dist), c(2, 1, 2), tolerance = 0)
  expect_equal(
    unname(SeuratObject::Embeddings(mapped[["ref.embeddings"]])[query_cells, , drop = FALSE]),
    unname(reference_umap[c("r1", "r2", "r3"), , drop = FALSE]),
    tolerance = 0
  )
})
