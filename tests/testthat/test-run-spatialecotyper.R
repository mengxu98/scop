make_spatialecotyper_seurat <- function() {
  counts <- matrix(
    c(
      4, 0, 1, 2,
      0, 3, 0, 5,
      2, 1, 4, 0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:3), paste0("Cell", 1:4))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$CellType <- c("T", "B", "T", "Myeloid")
  srt$X <- c(1, 2, 3, 4)
  srt$Y <- c(5, 6, 7, 8)
  srt
}

with_mock_spatialecotyper <- function(fun, code) {
  testthat::local_mocked_bindings(
    check_r = function(packages, ...) {
      expect_identical(packages, "digitalcytometry/SpatialEcoTyper")
      invisible(TRUE)
    },
    get_namespace_fun = function(package, name) {
      expect_identical(package, "SpatialEcoTyper")
      expect_identical(name, "SpatialEcoTyper")
      fun
    }
  )
  force(code)
}

test_that("RunSpatialEcoTyper writes SE labels and tool results", {
  srt <- make_spatialecotyper_seurat()
  fake_fun <- function(normdata, metadata, ...) {
    expect_identical(rownames(metadata), colnames(normdata))
    expect_named(metadata, c("X", "Y", "CellType"))
    expect_equal(metadata$CellType, unname(srt$CellType))
    list(
      metadata = data.frame(
        SE = c("SE1", "SE1", "SE2", "SE2"),
        row.names = colnames(normdata)
      ),
      obj = list()
    )
  }

  with_mock_spatialecotyper(fake_fun, {
    out <- RunSpatialEcoTyper(
      srt,
      celltype.by = "CellType",
      prefix = "SET",
      tool_name = "SET_tool",
      verbose = FALSE
    )
  })

  expect_equal(unname(out$SET_SE), c("SE1", "SE1", "SE2", "SE2"))
  expect_true("SET_tool" %in% names(out@tools))
  expect_equal(out@tools$SET_tool$parameters$prefix, "SET")
  expect_equal(out@tools$SET_tool$parameters$tool_name, "SET_tool")
})

test_that("RunSpatialEcoTyper errors on partial returned metadata by default", {
  srt <- make_spatialecotyper_seurat()
  fake_fun <- function(normdata, metadata, ...) {
    list(
      metadata = data.frame(
        SE = c("SE1", "SE2", "SE2"),
        row.names = colnames(normdata)[1:3]
      )
    )
  }

  with_mock_spatialecotyper(fake_fun, {
    expect_error(
      RunSpatialEcoTyper(srt, celltype.by = "CellType", verbose = FALSE),
      "returned no SE label"
    )
  })
})

test_that("RunSpatialEcoTyper can keep partial labels when requested", {
  srt <- make_spatialecotyper_seurat()
  fake_fun <- function(normdata, metadata, ...) {
    list(
      metadata = data.frame(
        SE = c("SE1", "SE2", "SE2"),
        row.names = colnames(normdata)[1:3]
      )
    )
  }

  with_mock_spatialecotyper(fake_fun, {
    out <- RunSpatialEcoTyper(
      srt,
      celltype.by = "CellType",
      allow_partial = TRUE,
      verbose = FALSE
    )
  })
  expect_equal(unname(out$SpatialEcoTyper_SE[1:3]), c("SE1", "SE2", "SE2"))
  expect_true(is.na(out$SpatialEcoTyper_SE[[4]]))
})
