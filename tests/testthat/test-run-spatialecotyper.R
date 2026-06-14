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
  srt$sample <- c("S1", "S1", "S2", "S2")
  srt$Region <- c("R1", "R1", "R2", "R2")
  srt$X <- c(1, 2, 3, 4)
  srt$Y <- c(5, 6, 7, 8)
  srt
}

with_mock_spatialecotyper <- function(funs, code) {
  if (is.function(funs)) {
    funs <- list(SpatialEcoTyper = funs)
  }
  testthat::local_mocked_bindings(
    check_r = function(packages, ...) {
      expect_identical(packages, "digitalcytometry/SpatialEcoTyper")
      invisible(TRUE)
    },
    get_namespace_fun = function(package, name) {
      expect_identical(package, "SpatialEcoTyper")
      expect_true(name %in% names(funs))
      funs[[name]]
    }
  )
  force(code)
}

test_that("RunSpatialEcoTyper writes single-sample SE labels and tool results", {
  srt <- make_spatialecotyper_seurat()
  fake_fun <- function(normdata, metadata, ...) {
    expect_s4_class(normdata, "dgCMatrix")
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
  expect_equal(out@tools$SET_tool$parameters$mode, "single")
  expect_equal(out@tools$SET_tool$parameters$prefix, "SET")
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
      "returned no SE result"
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

test_that("RunSpatialEcoTyper supports multi-sample conserved SE discovery", {
  srt <- make_spatialecotyper_seurat()
  fake_multi <- function(
    data_list,
    metadata_list,
    outdir,
    radius,
    ncores,
    Region,
    npcs,
    min.cells,
    iterations,
    grid.size,
    k,
    k.sn,
    dropcell,
    ...
  ) {
    expect_named(data_list, c("S1", "S2"))
    expect_named(metadata_list, c("S1", "S2"))
    expect_identical(colnames(data_list$S1), c("Cell1", "Cell2"))
    expect_identical(colnames(data_list$S2), c("Cell3", "Cell4"))
    expect_equal(metadata_list$S1$Region, c("R1", "R1"))
    expect_equal(metadata_list$S2$Region, c("R2", "R2"))
    expect_true(dir.exists(outdir))
    expect_identical(radius, 50)
    expect_identical(ncores, 4)
    expect_identical(Region, "Region")
    expect_identical(npcs, 11)
    expect_identical(min.cells, 2)
    expect_identical(iterations, 3)
    expect_identical(grid.size, 25)
    expect_identical(k, 5)
    expect_identical(k.sn, 7)
    expect_false(dropcell)
    data.frame(
      CID = paste0("Cell", 1:4),
      Sample = c("S1", "S1", "S2", "S2"),
      InitSE = c("Init1", "Init1", "Init2", "Init2"),
      SE = c("CSE1", "CSE1", "CSE2", "CSE2")
    )
  }

  with_mock_spatialecotyper(list(MultiSpatialEcoTyper = fake_multi), {
    out <- RunSpatialEcoTyper(
      srt,
      mode = "multi",
      celltype.by = "CellType",
      sample.by = "sample",
      Region = "Region",
      npcs = 11,
      min.cells = 2,
      iterations = 3,
      grid.size = 25,
      k = 5,
      k.sn = 7,
      dropcell = FALSE,
      prefix = "SET",
      verbose = FALSE
    )
  })

  expect_equal(unname(out$SET_InitSE), c("Init1", "Init1", "Init2", "Init2"))
  expect_equal(unname(out$SET_SE), c("CSE1", "CSE1", "CSE2", "CSE2"))
  expect_equal(out@tools$SpatialEcoTyper$parameters$mode, "multi")
})

test_that("RunSpatialEcoTyper recovers pretrained SE labels", {
  srt <- make_spatialecotyper_seurat()
  fake_recover <- function(dat, celltypes, scale, ncell.per.run, min.score, ...) {
    expect_s4_class(dat, "dgCMatrix")
    expect_identical(celltypes, as.character(srt$CellType))
    expect_true(scale)
    expect_identical(ncell.per.run, 500)
    expect_identical(min.score, 0.6)
    data.frame(
      CID = paste0("Cell", 1:4),
      CellType = celltypes,
      InitSE = c("Init1", "Init1", "Init2", "Init2"),
      SE = c("SE1", "SE1", "SE2", "SE2"),
      PredScore = c(0.91, 0.83, 0.77, 0.69)
    )
  }

  with_mock_spatialecotyper(list(RecoverSE = fake_recover), {
    out <- RunSpatialEcoTyper(
      srt,
      mode = "recover",
      celltype.by = "CellType",
      prefix = "SET",
      verbose = FALSE
    )
  })

  expect_equal(unname(out$SET_SE), c("SE1", "SE1", "SE2", "SE2"))
  expect_equal(unname(out$SET_InitSE), c("Init1", "Init1", "Init2", "Init2"))
  expect_equal(unname(out$SET_PredScore), c(0.91, 0.83, 0.77, 0.69))
})

test_that("RunSpatialEcoTyper aligns named recovery celltypes before length checks", {
  srt <- make_spatialecotyper_seurat()
  dat <- GetAssayData5(srt, layer = "data")[, c("Cell1", "Cell3"), drop = FALSE]
  full_celltypes <- stats::setNames(as.character(srt$CellType), colnames(srt))
  fake_recover <- function(dat, celltypes, ...) {
    expect_identical(unname(celltypes), c("T", "T"))
    data.frame(
      CID = colnames(dat),
      CellType = celltypes,
      InitSE = c("Init1", "Init2"),
      SE = c("SE1", "SE2"),
      PredScore = c(0.91, 0.82)
    )
  }

  with_mock_spatialecotyper(list(RecoverSE = fake_recover), {
    out <- RunSpatialEcoTyper(
      srt,
      mode = "recover",
      dat = dat,
      celltypes = full_celltypes,
      prefix = "SET",
      allow_partial = TRUE,
      verbose = FALSE
    )
  })

  expect_equal(out$SET_SE[["Cell1"]], "SE1")
  expect_equal(out$SET_SE[["Cell3"]], "SE2")
  expect_true(is.na(out$SET_SE[["Cell2"]]))
})

test_that("RunSpatialEcoTyper writes deconvoluted SE abundance to metadata", {
  srt <- make_spatialecotyper_seurat()
  fake_deconv <- function(dat, scale, nsample.per.run, sum2one, ...) {
    expect_s4_class(dat, "dgCMatrix")
    expect_true(scale)
    expect_identical(nsample.per.run, 500)
    expect_true(sum2one)
    matrix(
      c(
        0.8, 0.7, 0.2, 0.1,
        0.2, 0.3, 0.8, 0.9
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("SE1", "SE2"), paste0("Cell", 1:4))
    )
  }

  with_mock_spatialecotyper(list(DeconvoluteSE = fake_deconv), {
    out <- RunSpatialEcoTyper(
      srt,
      mode = "deconvolute",
      prefix = "SET",
      verbose = FALSE
    )
  })

  expect_equal(unname(out$SET_Abundance_SE1), c(0.8, 0.7, 0.2, 0.1))
  expect_equal(unname(out$SET_Abundance_SE2), c(0.2, 0.3, 0.8, 0.9))
  expect_equal(unname(out$SET_DominantSE), c("SE1", "SE1", "SE2", "SE2"))
})

test_that("RunSpatialEcoTyper returns deconvolution matrix for matrix input", {
  expr <- matrix(
    1:6,
    nrow = 3,
    dimnames = list(paste0("Gene", 1:3), paste0("Bulk", 1:2))
  )
  fake_deconv <- function(dat, ...) {
    expect_identical(dat, expr)
    matrix(
      c(0.4, 0.6, 0.5, 0.5),
      nrow = 2,
      dimnames = list(c("SE1", "SE2"), paste0("Bulk", 1:2))
    )
  }

  with_mock_spatialecotyper(list(DeconvoluteSE = fake_deconv), {
    out <- RunSpatialEcoTyper(expr, mode = "deconvolute", verbose = FALSE)
  })

  expect_equal(dim(out), c(2, 2))
  expect_equal(rownames(out), c("SE1", "SE2"))
})

test_that("SpatialEcoTyper plotting helpers return ggplot objects", {
  srt <- make_spatialecotyper_seurat()
  srt$SpatialEcoTyper_SE <- c("SE1", "SE1", "SE2", "SE2")
  p1 <- SpatialEcoTyperSpatialPlot(
    srt,
    group.by = "SpatialEcoTyper_SE",
    overlay_image = FALSE
  )
  p2 <- SpatialEcoTyperCompositionPlot(
    srt,
    se.by = "SpatialEcoTyper_SE",
    verbose = FALSE
  )
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})
