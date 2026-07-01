make_secact_seurat <- function() {
  counts <- matrix(
    c(
      4, 0, 1,
      0, 3, 2,
      1, 2, 5
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:3), paste0("Cell", 1:3))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$celltype <- c("A", "A", "B")
  srt$condition <- c("case", "control", "case")
  srt
}

with_mock_secact <- function(funs, code) {
  testthat::local_mocked_bindings(
    .package = "scop",
    check_r = function(packages, ...) {
      expect_identical(packages, "data2intelligence/SecAct")
      invisible(TRUE)
    },
    secact_namespace_available = function() {
      TRUE
    },
    get_namespace_fun = function(pkg, fun) {
      expect_identical(pkg, "SecAct")
      funs[[fun]]
    }
  )
  force(code)
}

test_that("RunSecAct stores matrix activity as a Seurat assay", {
  srt <- make_secact_seurat()
  expr <- matrix(
    1,
    nrow = 3,
    ncol = 3,
    dimnames = list(paste0("Gene", 1:3), paste0("Cell", 1:3))
  )
  activity <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 2,
    dimnames = list(c("TGFB1", "LY86"), paste0("Cell", 1:3))
  )
  funs <- list(
    SecAct.activity.inference = function(inputProfile, is.differential, ...) {
      expect_identical(inputProfile, expr)
      expect_true(is.differential)
      list(beta = activity, se = activity, zscore = activity, pvalue = activity)
    }
  )

  with_mock_secact(funs, {
    out <- RunSecAct(
      srt = srt,
      inputProfile = expr,
      mode = "matrix",
      is.differential = TRUE,
      verbose = FALSE
    )
  })

  expect_true("SecAct" %in% names(out@assays))
  expect_true("SecAct" %in% names(out@tools))
  expect_equal(rownames(out[["SecAct"]]), c("TGFB1", "LY86"))
})

test_that("RunSecAct does not store non-matrix activity as an assay", {
  srt <- make_secact_seurat()
  expr <- matrix(
    1,
    nrow = 3,
    ncol = 3,
    dimnames = list(paste0("Gene", 1:3), paste0("Cell", 1:3))
  )
  funs <- list(
    SecAct.activity.inference = function(...) {
      list(beta = 1:3, se = 1:3, zscore = 1:3, pvalue = 1:3)
    }
  )

  with_mock_secact(funs, {
    out <- suppressWarnings(RunSecAct(
      srt = srt,
      inputProfile = expr,
      mode = "matrix",
      verbose = FALSE
    ))
  })

  expect_false("SecAct" %in% names(out@assays))
  expect_true("SecAct" %in% names(out@tools))
})

test_that("RunSecAct scRNAseq dispatches to SecAct and stores single-cell assay", {
  srt <- make_secact_seurat()
  activity <- matrix(
    c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
    nrow = 2,
    dimnames = list(c("TGFB1", "LY86"), paste0("Cell", 1:3))
  )
  funs <- list(
    SecAct.activity.inference.scRNAseq = function(
      inputProfile,
      cellType_meta,
      is.singleCellLevel,
      ...
    ) {
      expect_s4_class(inputProfile, "Seurat")
      expect_null(cellType_meta)
      expect_true(is.singleCellLevel)
      inputProfile@misc$SecAct_output$SecretedProteinActivity <- list(
        beta = activity,
        se = activity,
        zscore = activity,
        pvalue = activity
      )
      inputProfile
    }
  )

  with_mock_secact(funs, {
    out <- RunSecAct(
      srt = srt,
      mode = "scRNAseq",
      is.singleCellLevel = TRUE,
      verbose = FALSE
    )
  })

  expect_true("SecAct" %in% names(out@assays))
  expect_true("SecAct" %in% names(out@tools))
})

test_that("RunSecAct scRNAseq return_seurat = FALSE returns activity only", {
  srt <- make_secact_seurat()
  activity <- matrix(
    c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
    nrow = 2,
    dimnames = list(c("TGFB1", "LY86"), paste0("Cell", 1:3))
  )
  funs <- list(
    SecAct.activity.inference.scRNAseq = function(inputProfile, ...) {
      inputProfile@misc$SecAct_output$SecretedProteinActivity <- list(
        beta = activity,
        se = activity,
        zscore = activity,
        pvalue = activity
      )
      inputProfile
    }
  )

  with_mock_secact(funs, {
    out <- RunSecAct(
      srt = srt,
      mode = "scRNAseq",
      is.singleCellLevel = TRUE,
      return_seurat = FALSE,
      verbose = FALSE
    )
  })

  expect_type(out, "list")
  expect_equal(out$zscore, activity)
})

test_that("RunSecAct rejects non-RNA Seurat assays explicitly", {
  srt <- make_secact_seurat()
  srt[["ADT"]] <- SeuratObject::CreateAssayObject(counts = SeuratObject::GetAssayData(srt, assay = "RNA", layer = "counts"))
  SeuratObject::DefaultAssay(srt) <- "ADT"

  with_mock_secact(list(), {
    expect_error(
      RunSecAct(
        srt = srt,
        mode = "scRNAseq",
        is.singleCellLevel = TRUE,
        verbose = FALSE
      ),
      "assay = 'RNA'"
    )
  })
})

test_that("RunSecActCCC records scRNAseq CCC runs", {
  srt <- make_secact_seurat()
  funs <- list(
    SecAct.CCC.scRNAseq = function(
      Seurat_obj,
      cellType_meta,
      condition_meta,
      conditionCase,
      conditionControl,
      ...
    ) {
      expect_identical(cellType_meta, "celltype")
      expect_identical(condition_meta, "condition")
      expect_identical(conditionCase, "case")
      expect_identical(conditionControl, "control")
      Seurat_obj
    }
  )

  with_mock_secact(funs, {
    out <- RunSecActCCC(
      srt = srt,
      mode = "scRNAseq",
      cellType_meta = "celltype",
      condition_meta = "condition",
      conditionCase = "case",
      conditionControl = "control",
      verbose = FALSE
    )
  })

  expect_true("SecAct_CCC" %in% names(out@tools))
})

test_that("RunSecActCCC allows condition_meta NULL and ignores case labels", {
  srt <- make_secact_seurat()
  funs <- list(
    SecAct.CCC.scRNAseq = function(
      Seurat_obj,
      condition_meta,
      conditionCase,
      conditionControl,
      ...
    ) {
      expect_null(condition_meta)
      expect_null(conditionCase)
      expect_null(conditionControl)
      Seurat_obj
    }
  )

  with_mock_secact(funs, {
    expect_warning(
      out <- RunSecActCCC(
        srt = srt,
        mode = "scRNAseq",
        cellType_meta = "celltype",
        conditionCase = "case",
        conditionControl = "control",
        verbose = FALSE
      ),
      "ignored"
    )
  })

  expect_true("SecAct_CCC" %in% names(out@tools))
})

test_that("SecAct SpaCET downstream wrappers dispatch to SecAct functions", {
  spacet <- structure(list(id = "mock"), class = "SpaCET")
  funs <- list(
    SecAct.signaling.pattern = function(SpaCET_obj, k, ...) {
      expect_identical(SpaCET_obj, spacet)
      expect_equal(k, 3)
      structure(list(pattern = TRUE), class = "SpaCET")
    },
    SecAct.signaling.pattern.gene = function(SpaCET_obj, n) {
      expect_equal(n, 2)
      matrix(1, nrow = 1, dimnames = list("TGFB1", "pattern2"))
    },
    SecAct.signaling.velocity.spotST = function(SpaCET_obj, gene, signalMode, radius, ...) {
      expect_identical(gene, "TGFB1")
      expect_identical(signalMode, "receiving")
      expect_equal(radius, 200)
      ggplot2::ggplot()
    },
    SecAct.signaling.velocity.scST = function(SpaCET_obj, colors, radius, ...) {
      expect_equal(radius, 20)
      expect_false(is.null(colors))
      ggplot2::ggplot()
    }
  )

  with_mock_secact(funs, {
    pattern_obj <- RunSecActSignalingPattern(spacet, k = 3, verbose = FALSE)
    genes <- RunSecActPatternGenes(pattern_obj, n = 2, verbose = FALSE)
    p <- RunSecActVelocity(
      spacet,
      mode = "spotST",
      gene = "TGFB1",
      verbose = FALSE
    )
    p2 <- RunSecActVelocity(
      spacet,
      mode = "scST",
      sender = "A",
      secretedProtein = "TGFB1",
      receiver = "B",
      cellType_meta = "celltype",
      verbose = FALSE
    )
  })

  expect_s3_class(pattern_obj, "SpaCET")
  expect_true(is.matrix(genes))
  expect_s3_class(p, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("RunSecActVelocity rejects unsupported spotST signalMode", {
  spacet <- structure(list(id = "mock"), class = "SpaCET")
  with_mock_secact(list(), {
    expect_error(
      RunSecActVelocity(
        spacet,
        mode = "spotST",
        gene = "TGFB1",
        signalMode = "both",
        verbose = FALSE
      ),
      "should be one of"
    )
  })
})
