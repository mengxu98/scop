make_pseudobulk_object <- function(seed = 20260830) {
  set.seed(seed)
  counts <- Matrix::rsparsematrix(200, 120, density = 0.08)
  counts@x <- stats::rpois(length(counts@x), lambda = 4) + 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  object <- Seurat::CreateSeuratObject(counts)
  object$group <- factor(rep(c("A", "B", "C"), length.out = ncol(object)))
  object$batch <- factor(rep(c("x", "y"), each = ncol(object) / 2L))
  SeuratObject::Idents(object) <- object$group
  NormalizeData(object, verbose = FALSE)
}

test_that("AggregateExpression delegation matches Seurat", {
  object <- make_pseudobulk_object()
  features <- c("g17", "g3", "g101", "g8")
  expected <- Seurat::AggregateExpression(
    object,
    assays = "RNA",
    features = features,
    group.by = "group",
    verbose = FALSE
  )
  actual <- AggregateExpression(
    object,
    assays = "RNA",
    features = features,
    group.by = "group",
    verbose = FALSE
  )
  expect_identical(names(actual), "RNA")
  expect_equal(actual$RNA, expected$RNA, tolerance = 0)
  expect_identical(rownames(actual$RNA), rownames(expected$RNA))
})

test_that("AverageExpression delegation matches Seurat", {
  object <- make_pseudobulk_object()
  expected <- suppressMessages(Seurat::AverageExpression(
    object,
    assays = "RNA",
    group.by = "group",
    verbose = FALSE
  ))
  counts_expected <- suppressMessages(Seurat::AverageExpression(
    object,
    assays = "RNA",
    group.by = "ident",
    layer = "counts",
    verbose = FALSE
  ))
  actual <- suppressMessages(AverageExpression(
    object,
    assays = "RNA",
    group.by = "group",
    verbose = FALSE
  ))
  expect_equal(actual$RNA, expected$RNA, tolerance = 1e-12)
  counts_actual <- suppressMessages(AverageExpression(
    object,
    assays = "RNA",
    group.by = "ident",
    layer = "counts",
    verbose = FALSE
  ))
  expect_equal(counts_actual$RNA, counts_expected$RNA, tolerance = 0)
})

test_that("Pseudobulk delegation covers multi-group and split layers", {
  object <- make_pseudobulk_object()
  aggregate_actual <- AggregateExpression(
    object,
    group.by = c("group", "batch"),
    verbose = FALSE
  )
  aggregate_expected <- Seurat::AggregateExpression(
    object,
    group.by = c("group", "batch"),
    verbose = FALSE
  )
  expect_equal(aggregate_actual, aggregate_expected, tolerance = 0)

  seurat_actual <- AggregateExpression(
    object,
    group.by = "group",
    return.seurat = TRUE,
    verbose = FALSE
  )
  seurat_expected <- Seurat::AggregateExpression(
    object,
    group.by = "group",
    return.seurat = TRUE,
    verbose = FALSE
  )
  expect_s4_class(seurat_actual, "Seurat")
  expect_equal(
    as.matrix(SeuratObject::LayerData(seurat_actual, layer = "counts")),
    as.matrix(SeuratObject::LayerData(seurat_expected, layer = "counts")),
    tolerance = 0
  )

  split_object <- object
  split_object[["RNA"]] <- split(split_object[["RNA"]], f = split_object$batch)
  average_actual <- suppressMessages(AverageExpression(
    split_object,
    group.by = "group",
    verbose = FALSE
  ))
  average_expected <- suppressMessages(Seurat::AverageExpression(
    split_object,
    group.by = "group",
    verbose = FALSE
  ))
  expect_equal(average_actual, average_expected, tolerance = 1e-12)
})

test_that("Pseudobulk delegation covers multiple assays and feature lists", {
  object <- make_pseudobulk_object(seed = 20260906)
  object[["ADT"]] <- SeuratObject::CreateAssay5Object(
    counts = SeuratObject::LayerData(object[["RNA"]], layer = "counts")[1:40, ]
  )
  args <- list(
    object = object,
    assays = c("RNA", "ADT"),
    features = list(paste0("g", 1:30), paste0("g", 1:20)),
    group.by = c("group", "batch"),
    verbose = FALSE
  )
  actual <- do.call(AggregateExpression, args)
  expected <- do.call(Seurat::AggregateExpression, args)
  expect_identical(names(actual), names(expected))
  expect_equal(actual, expected, tolerance = 0)
})
