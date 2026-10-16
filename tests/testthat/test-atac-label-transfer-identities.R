make_label_transfer_fixture <- function(reverse_query = FALSE) {
  testthat::skip_if_not_installed("Signac")
  testthat::skip_if_not_installed("withr")
  withr::with_seed(91, {
    genes <- paste0("GENE", seq_len(150))
    make_counts <- function(n, prefix) {
      mu <- matrix(8, 150, n)
      mu[1:30, seq_len(n / 2)] <- 24
      mu[31:60, (n / 2 + 1):n] <- 24
      Matrix::Matrix(matrix(rnbinom(length(mu), mu = mu, size = 10), 150,
        dimnames = list(genes, paste0(prefix, seq_len(n)))
      ), sparse = TRUE)
    }
    reference <- SeuratObject::CreateSeuratObject(make_counts(60, "ref"))
    reference$label <- rep(c("classA", "classB"), each = 30)
    reference <- Seurat::NormalizeData(reference, verbose = FALSE)
    reference <- Seurat::FindVariableFeatures(reference, nfeatures = 150, verbose = FALSE)
    reference <- Seurat::ScaleData(reference, features = genes, verbose = FALSE)
    reference <- Seurat::RunPCA(reference, features = genes, npcs = 10, verbose = FALSE)
    activity <- make_counts(40, "query")
    if (reverse_query) activity <- activity[, rev(colnames(activity)), drop = FALSE]
    peaks <- Matrix::Matrix(matrix(rpois(200 * 40, 4), 200,
      dimnames = list(paste0("chr1-", seq(1, by = 200, length.out = 200), "-", seq(100, by = 200, length.out = 200)), colnames(activity))
    ), sparse = TRUE)
    chromatin <- Signac::CreateChromatinAssay(counts = peaks, genome = NULL)
    query <- SeuratObject::CreateSeuratObject(chromatin, assay = "peaks")
    query[["ACTIVITY"]] <- SeuratObject::CreateAssayObject(counts = activity)
    query <- Seurat::NormalizeData(query, assay = "ACTIVITY", verbose = FALSE)
    query <- Signac::RunTFIDF(query, assay = "peaks", verbose = FALSE)
    query <- Signac::FindTopFeatures(query, assay = "peaks", min.cutoff = "q0", verbose = FALSE)
    query <- Signac::RunSVD(query, assay = "peaks", n = 10, verbose = FALSE)
    list(reference = reference, query = query, features = genes)
  })
}
run_label_wrapper <- function(x, ...) {
  RunLabelTransfer(x$query, x$reference,
    assay = "peaks", reference_assay = "RNA", reference_reduction = "pca", reference_dims = 1:5,
    reference_label = "label", add_gene_activity = FALSE, gene_activity_assay = "ACTIVITY",
    weight_reduction = "lsi", dims = 2:5, features = x$features, k.weight = 5, verbose = FALSE, ...
  )
}
run_label_native <- function(x) {
  anchors <- Seurat::FindTransferAnchors(
    reference = x$reference, query = x$query,
    reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "pcaproject",
    reference.reduction = "pca", dims = 1:5, features = x$features, k.anchor = 5, k.filter = 39, k.score = 30, verbose = FALSE
  )
  Seurat::TransferData(
    anchorset = anchors, refdata = x$reference$label,
    weight.reduction = x$query[["lsi"]], dims = 2:5, k.weight = 5, verbose = FALSE
  )
}

label_prediction_frame <- function(object, columns, prefix = "predicted_") {
  out <- object[[]][, paste0(prefix, columns), drop = FALSE]
  colnames(out) <- columns
  out
}

test_that("disjoint reference-query transfer matches native and preserves sources", {
  x <- make_label_transfer_fixture()
  expect_length(intersect(colnames(x$query), colnames(x$reference)), 0L)
  x$reference$label <- factor(x$reference$label, levels = c("classA", "classB", "unused"))
  expected <- withr::with_seed(29, run_label_native(x))
  before <- list(query = serialize(x$query, NULL), reference = serialize(x$reference, NULL))
  set.seed(101)
  seed <- .Random.seed
  actual <- run_label_wrapper(x)
  expect_identical(.Random.seed, seed)
  expect_equal(label_prediction_frame(actual, colnames(expected)), expected, tolerance = 1e-12)
  expect_identical(serialize(x$query, NULL), before$query)
  expect_identical(serialize(x$reference, NULL), before$reference)
  expect_equal(actual[["peaks"]], x$query[["peaks"]])
  expect_equal(actual[["ACTIVITY"]], x$query[["ACTIVITY"]])
  expect_identical(actual[["lsi"]], x$query[["lsi"]])
  expect_equal(actual[[]][, colnames(x$query[[]]), drop = FALSE], x$query[[]])
  expect_identical(colnames(actual), colnames(x$query))
  expect_identical(SeuratObject::DefaultAssay(actual), "peaks")
})

test_that("query token IDs and cross-object repeated IDs retain complete predictions", {
  x <- make_label_transfer_fixture()
  expected <- withr::with_seed(29, run_label_native(x))
  for (case in c("mixed", "collision", "shared")) {
    ids <- paste0("cell", 1:40)
    if (case == "mixed") ids[1] <- "cell_query1"
    if (case == "collision") ids[2] <- "cell_query1"
    if (case == "shared") ids <- colnames(x$reference)[1:40]
    input <- x
    input$query <- SeuratObject::RenameCells(input$query, new.names = ids)
    actual <- withr::with_seed(29, run_label_wrapper(input))
    want <- expected
    rownames(want) <- ids
    expect_identical(colnames(actual), ids)
    expect_false(anyNA(actual$predicted_predicted.id))
    expect_equal(label_prediction_frame(actual, colnames(want)), want, tolerance = 1e-12)
  }
})

test_that("actual feature receipts and literal prefixes describe only current predictions", {
  x <- make_label_transfer_fixture()
  x$features <- NULL
  out <- withr::with_seed(29, run_label_wrapper(x, prediction_prefix = "map+"))
  receipt <- out@tools$ATAC_Seurat_LabelTransfer
  expect_setequal(receipt$features, SeuratObject::VariableFeatures(x$reference))
  expect_identical(
    receipt$probability_cols,
    colnames(out[[]])[startsWith(colnames(out[[]]), "map+prediction.score.")]
  )
  expect_length(receipt$probability_cols, 3L)
  x$features <- c("not-a-gene", rownames(x$reference)[1:100])
  out <- withr::with_seed(29, run_label_wrapper(x))
  expect_identical(out@tools$ATAC_Seurat_LabelTransfer$features, rownames(x$reference)[1:100])
})

test_that("reruns replace only recorded result columns and protect unrelated metadata", {
  x <- make_label_transfer_fixture()
  out <- withr::with_seed(29, run_label_wrapper(x))
  out$predicted_user_note <- "retain"
  x$query <- out
  x$reference$label <- paste0("new_", x$reference$label)
  updated <- withr::with_seed(29, run_label_wrapper(x))
  expect_false(any(c("predicted_prediction.score.classA", "predicted_prediction.score.classB") %in% colnames(updated[[]])))
  expect_true(all(c("predicted_prediction.score.new_classA", "predicted_prediction.score.new_classB") %in% colnames(updated[[]])))
  expect_identical(updated$predicted_user_note, out$predicted_user_note)
  expect_length(updated@tools$ATAC_Seurat_LabelTransfer$probability_cols, 3L)
  fresh <- make_label_transfer_fixture()
  fresh$query$predicted_predicted.id <- "user-owned"
  seed <- .Random.seed
  expect_error(run_label_wrapper(fresh), "outside this tool's recorded results")
  expect_identical(.Random.seed, seed)
  expect_true(all(fresh$query$predicted_predicted.id == "user-owned"))
})

test_that("explicit existing activity and dimensions retain clear failure boundaries", {
  x <- make_label_transfer_fixture()
  bad <- x
  bad$features <- "missing"
  expect_error(run_label_wrapper(bad), "No shared features")
  expect_error(RunLabelTransfer(x$query, x$reference,
    assay = "peaks",
    add_gene_activity = FALSE, gene_activity_assay = "ACTIVITY", reference_label = "label",
    weight_reduction = "lsi", dims = 999, verbose = FALSE
  ), "No valid dimensions")
  expect_error(RunLabelTransfer(x$query, x$reference,
    assay = "peaks",
    add_gene_activity = FALSE, gene_activity_assay = "ACTIVITY", reference_label = "absent",
    weight_reduction = "lsi", dims = 2:5, verbose = FALSE
  ), "reference_label")
  withr::local_preserve_seed()
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  invisible(run_label_wrapper(x))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("shuffled valid query cells stay aligned with native predictions", {
  x <- make_label_transfer_fixture(reverse_query = TRUE)
  SeuratObject::Idents(x$query) <- rep(c("query-group-B", "query-group-A"), each = 20)
  requested <- paste0("query", 40:1)
  expect_identical(colnames(x$query), requested)
  expect_true(validObject(x$query))
  expected <- withr::with_seed(29, run_label_native(x))
  before <- serialize(x$query, NULL)
  actual <- withr::with_seed(29, run_label_wrapper(x))
  expect_identical(colnames(actual), colnames(x$query))
  expect_identical(SeuratObject::Idents(actual), SeuratObject::Idents(x$query))
  expect_equal(label_prediction_frame(actual, colnames(expected)), expected, tolerance = 1e-12)
  expect_identical(serialize(x$query, NULL), before)
})

test_that("shared transfer helper retains the object-returning mapping contract", {
  x <- make_label_transfer_fixture()
  expected <- withr::with_seed(29, run_label_native(x))
  actual <- withr::with_seed(29, scop:::atac_transfer_labels(
    srt = x$query, assay = "peaks", prefix = "ATAC",
    linear_reduction_dims_use = 2:5, reference = x$reference,
    reference_assay = "RNA", reference_reduction = "pca", reference_dims = 1:5,
    reference_label = "label", add_gene_activity = FALSE,
    gene_activity_assay = "ACTIVITY", weight_reduction = "lsi", features = x$features,
    prediction_prefix = "predicted_", k.weight = 5, verbose = FALSE
  ))
  expect_s4_class(actual, "Seurat")
  expect_equal(label_prediction_frame(actual, colnames(expected)), expected, tolerance = 1e-12)
  expect_identical(actual@tools, x$query@tools)
})
