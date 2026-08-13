test_that("RunLISI forwards and records resource controls", {
  set.seed(41)
  counts <- matrix(
    stats::rpois(120, lambda = 2),
    nrow = 6,
    dimnames = list(paste0("gene", 1:6), paste0("cell", 1:20))
  )
  srt <- SeuratObject::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE)
  )
  embedding <- matrix(
    stats::rnorm(80),
    nrow = 20,
    dimnames = list(colnames(srt), paste0("DEMO_", 1:4))
  )
  srt[["demo"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    key = "DEMO_",
    assay = SeuratObject::DefaultAssay(srt)
  )
  srt$batch <- rep(c("A", "B"), each = 10)

  result <- RunLISI(
    srt,
    reductions = "demo",
    label_colnames = "batch",
    perplexity = 5,
    knn_algorithm = "brute_force",
    cores = 1,
    max_dense_bytes = 1024^2,
    verbose = FALSE
  )

  expect_true("demo_batch_LISI" %in% colnames(result@meta.data))
  expect_identical(result@tools$demo_LISI$knn_algorithm, "brute_force")
  expect_identical(result@tools$demo_LISI$cores, 1L)
  expect_identical(result@tools$demo_LISI$max_dense_bytes, 1024^2)
  expect_lt(result@tools$demo_LISI$estimated_dense_bytes$demo, 1024^2)
})

test_that("RunLISI enforces the dense-memory limit before thisutils", {
  set.seed(42)
  counts <- Matrix::Matrix(matrix(stats::rpois(120, 2), nrow = 6), sparse = TRUE)
  srt <- SeuratObject::CreateSeuratObject(counts = counts)
  embedding <- matrix(
    stats::rnorm(80),
    nrow = 20,
    dimnames = list(colnames(srt), paste0("DEMO_", 1:4))
  )
  srt[["demo"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    key = "DEMO_",
    assay = SeuratObject::DefaultAssay(srt)
  )
  srt$batch <- rep(c("A", "B"), each = 10)

  expect_error(
    RunLISI(
      srt,
      reductions = "demo",
      label_colnames = "batch",
      max_dense_bytes = 100,
      verbose = FALSE
    ),
    "exceeds.*max_dense_bytes"
  )
})

test_that("RunLISI keeps automatic resource defaults", {
  expect_null(formals(RunLISI)$cores)
  expect_identical(eval(formals(RunLISI)$max_dense_bytes), Inf)
})
