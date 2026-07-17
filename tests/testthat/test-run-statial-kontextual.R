make_statial_kontextual_seurat <- function() {
  counts <- matrix(
    c(
      5, 0, 2, 1, 3, 0,
      0, 3, 0, 4, 1, 2,
      1, 1, 1, 1, 1, 1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:3), paste0("Cell", 1:6))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$cell_type <- c("A", "B", "A", "C", "B", "C")
  srt$sample <- c("S1", "S1", "S1", "S2", "S2", "S2")
  srt$x <- c(0, 1, 2, 0, 1, 2)
  srt$y <- c(0, 0, 0, 1, 1, 1)
  srt
}

with_mock_statial_kontextual <- function(code) {
  fake_kontextual <- function(
    cells,
    r,
    parentDf,
    from,
    to,
    parent,
    image,
    inhom,
    edgeCorrect,
    window,
    window.length,
    includeOriginal,
    spatialCoords,
    cellType,
    imageID,
    cores,
    ...
  ) {
    expect_true(all(c("imageID", "cellType", "x", "y") %in% colnames(cells)))
    expect_equal(spatialCoords, c("x", "y"))
    expect_identical(cellType, "cellType")
    expect_identical(imageID, "imageID")
    expect_equal(sort(unique(cells$imageID)), c("S1", "S2"))
    expect_equal(r, c(10, 20))
    expect_equal(from, "A")
    expect_equal(to, "B")
    expect_equal(parent, c("A", "B", "C"))
    expect_false(inhom)
    expect_true(edgeCorrect)
    expect_identical(window, "convex")
    expect_true(includeOriginal)
    expect_identical(cores, 1L)
    data.frame(
      imageID = c("S1", "S2"),
      test = "A__B",
      original = c(0.2, -0.1),
      kontextual = c(1.5, -0.7),
      r = c(10, 20),
      inhomL = FALSE,
      stringsAsFactors = FALSE
    )
  }
  testthat::local_mocked_bindings(
    check_r = function(packages, ...) {
      expect_identical(packages, "Statial")
      invisible(TRUE)
    },
    get_namespace_fun = function(package, name) {
      expect_identical(package, "Statial")
      expect_identical(name, "Kontextual")
      fake_kontextual
    }
  )
  force(code)
}

test_that("RunStatialKontextual stores Kontextual results and summary", {
  srt <- make_statial_kontextual_seurat()
  with_mock_statial_kontextual({
    out <- RunStatialKontextual(
      srt,
      group.by = "cell_type",
      r = c(10, 20),
      from = "A",
      to = "B",
      parent = c("A", "B", "C"),
      sample.by = "sample",
      coord.cols = c("x", "y"),
      store_input = TRUE,
      verbose = FALSE
    )
  })

  expect_s4_class(out, "Seurat")
  expect_true("StatialKontextual" %in% names(out@tools))
  bundle <- out@tools$StatialKontextual
  expect_identical(bundle$source$coordinate_space, "raw")
  expect_named(bundle, c(
    "table", "raw", "summary", "parameters", "input", "method",
    "schema_version", "result_type", "source", "provenance"
  ))
  expect_named(bundle$table[1:4], c("imageID", "test", "original", "kontextual"))
  expect_equal(bundle$table$kontextual, c(1.5, -0.7))
  expect_equal(bundle$summary$n_records, 2)
  expect_equal(bundle$summary$n_images, 2)
  expect_equal(bundle$summary$n_localized, 1)
  expect_equal(bundle$summary$n_dispersed, 1)
  expect_equal(bundle$parameters$group.by, "cell_type")
  expect_equal(bundle$parameters$r, c(10, 20))
  expect_equal(nrow(bundle$input), 6)
})

test_that("RunStatialKontextual validates relationship and input fields", {
  srt <- make_statial_kontextual_seurat()
  expect_error(
    RunStatialKontextual(
      srt,
      group.by = "missing",
      r = 10,
      from = "A",
      to = "B",
      parent = c("A", "B"),
      coord.cols = c("x", "y"),
      verbose = FALSE
    ),
    "group.by"
  )
  expect_error(
    RunStatialKontextual(
      srt,
      group.by = "cell_type",
      r = 10,
      coord.cols = c("x", "y"),
      verbose = FALSE
    ),
    "from"
  )
  expect_error(
    RunStatialKontextual(
      srt,
      group.by = "cell_type",
      r = -1,
      from = "A",
      to = "B",
      parent = c("A", "B"),
      coord.cols = c("x", "y"),
      verbose = FALSE
    ),
    "positive finite"
  )
})

test_that("RunStatialKontextual can skip result storage", {
  srt <- make_statial_kontextual_seurat()
  with_mock_statial_kontextual({
    out <- RunStatialKontextual(
      srt,
      group.by = "cell_type",
      r = c(10, 20),
      from = "A",
      to = "B",
      parent = c("A", "B", "C"),
      sample.by = "sample",
      coord.cols = c("x", "y"),
      store_results = FALSE,
      verbose = FALSE
    )
  })

  expect_false("StatialKontextual" %in% names(out@tools))
})
