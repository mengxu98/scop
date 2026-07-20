test_that("RunCellQC initializes count metrics from one counts-layer read", {
  counts <- methods::as(Matrix::Matrix(
    matrix(c(1, 0, 2, 0, 3, 0, 0, 4, 1, 2, 0, 5), nrow = 3),
    sparse = TRUE
  ), "dgCMatrix")
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt_missing <- srt
  srt_missing@meta.data[, c("nCount_RNA", "nFeature_RNA")] <- NULL

  calls <- 0L
  testthat::local_mocked_bindings(
    GetAssayData5 = function(...) {
      calls <<- calls + 1L
      counts
    },
    .package = "scop"
  )

  out <- cellqc_initialize_count_metrics(srt_missing, assay = "RNA")

  expect_equal(calls, 1L)
  expect_equal(unname(out$nCount_RNA), as.numeric(Matrix::colSums(counts)))
  expect_equal(unname(out$nFeature_RNA), as.numeric(Matrix::colSums(counts > 0)))
})

cellqc_feature_test_object <- function() {
  counts <- Matrix::Matrix(
    matrix(
      c(
        10, 0, 5, 0,
        0, 0, 0, 0,
        90, 0, 0, 0,
        0, 5, 20, 0,
        0, 95, 75, 100
      ),
      nrow = 5,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(counts) <- c("HBA1", "HBB", "HBP1", "CUSTOM1", "OTHER")
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  Seurat::CreateSeuratObject(counts = counts)
}

test_that("RunCellQC computes HB percentage without changing default filtering", {
  srt <- cellqc_feature_test_object()

  measured <- suppressWarnings(RunCellQC(
    srt,
    qc_metrics = character(),
    verbose = FALSE
  ))

  expect_equal(unname(measured$percent.hb), c(10, 0, 5, 0))
  expect_false("hb_qc" %in% colnames(measured[[]]))
  expect_true(all(measured$CellQC == "Pass"))

  filtered <- suppressWarnings(RunCellQC(
    srt,
    qc_metrics = "hb",
    verbose = FALSE
  ))

  expect_identical(as.character(filtered$hb_qc), c("Fail", "Pass", "Pass", "Pass"))
  expect_identical(as.character(filtered$CellQC), c("Fail", "Pass", "Pass", "Pass"))

  no_hb <- srt[c("CUSTOM1", "OTHER"), ]
  no_hb_measured <- suppressWarnings(RunCellQC(
    no_hb,
    qc_metrics = character(),
    verbose = FALSE
  ))
  expect_equal(unname(no_hb_measured$percent.hb), rep(0, ncol(no_hb)))
  expect_error(
    suppressWarnings(RunCellQC(
      no_hb,
      qc_metrics = "hb",
      verbose = FALSE
    )),
    "did not match"
  )
})

test_that("RunCellQC supports HB feature overrides", {
  srt <- cellqc_feature_test_object()

  out <- expect_warning(
    suppressMessages(RunCellQC(
      srt,
      qc_metrics = character(),
      hb_gene = c("HBP1", "MISSING"),
      verbose = FALSE
    )),
    "hb_gene ignored missing features"
  )

  expect_equal(unname(out$percent.hb), c(90, 0, 0, 0))
})

test_that("RunCellQC computes and filters named qc_features rules", {
  srt <- cellqc_feature_test_object()
  rules <- list(
    ambient = list(features = "CUSTOM1", range = c(0, 10)),
    hb_like = list(pattern = "^HB[AB]", range = c(0, 5))
  )

  measured <- suppressWarnings(RunCellQC(
    srt,
    qc_metrics = character(),
    qc_features = rules,
    verbose = FALSE
  ))
  expect_equal(unname(measured$percent.ambient), c(0, 5, 20, 0))
  expect_equal(unname(measured$percent.hb_like), c(10, 0, 5, 0))
  expect_false(any(c("ambient_qc", "hb_like_qc") %in% colnames(measured[[]])))

  filtered <- suppressWarnings(RunCellQC(
    srt,
    qc_metrics = c("ambient", "hb_like"),
    qc_features = rules,
    verbose = FALSE
  ))
  expect_identical(as.character(filtered$ambient_qc), c("Pass", "Pass", "Fail", "Pass"))
  expect_identical(as.character(filtered$hb_like_qc), c("Fail", "Pass", "Pass", "Pass"))
  expect_identical(as.character(filtered$CellQC), c("Fail", "Pass", "Fail", "Pass"))
})

test_that("RunCellQC preserves feature QC across splits and filtered returns", {
  srt <- cellqc_feature_test_object()
  srt$batch <- rep(c("a", "b"), each = 2)
  rules <- list(
    ambient = list(features = "CUSTOM1", range = c(0, 10))
  )

  out <- suppressWarnings(RunCellQC(
    srt,
    split.by = "batch",
    return_filtered = TRUE,
    qc_metrics = c("hb", "ambient"),
    qc_features = rules,
    verbose = FALSE
  ))

  expect_identical(colnames(out), c("c2", "c4"))
  expect_true(all(c("percent.hb", "percent.ambient") %in% colnames(out[[]])))
  expect_false(any(c("hb_qc", "ambient_qc", "CellQC") %in% colnames(out[[]])))
})

test_that("RunCellQC applies HB patterns per species", {
  counts <- Matrix::Matrix(
    matrix(
      c(
        10, 0,
        0, 10,
        90, 90,
        0, 0
      ),
      nrow = 4,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(counts) <- c(
    "human-HBA1", "mouse-Hba-a1", "human-OTHER", "mouse-Other"
  )
  colnames(counts) <- c("c1", "c2")
  srt <- Seurat::CreateSeuratObject(counts = counts)

  out <- suppressWarnings(RunCellQC(
    srt,
    qc_metrics = character(),
    species = c("human", "mouse"),
    species_gene_prefix = c("human", "mouse"),
    verbose = FALSE
  ))

  expect_equal(unname(out$percent.hb.human), c(10, 0))
  expect_equal(unname(out$percent.hb.mouse), c(0, 10))
})

test_that("qc_features validates rule names, definitions, and matches", {
  srt <- cellqc_feature_test_object()

  expect_error(
    RunCellQC(
      srt,
      qc_metrics = character(),
      qc_features = list(hb = list(features = "HBA1", range = c(0, 5))),
      verbose = FALSE
    ),
    "conflict with built-in"
  )
  expect_error(
    RunCellQC(
      srt,
      qc_metrics = character(),
      qc_features = list(bad = list(
        features = "HBA1",
        pattern = "^HB",
        range = c(0, 5)
      )),
      verbose = FALSE
    ),
    "exactly one"
  )
  expect_error(
    RunCellQC(
      srt,
      qc_metrics = character(),
      qc_features = list(bad = list(pattern = "^NOT_PRESENT", range = c(0, 5))),
      verbose = FALSE
    ),
    "did not match"
  )
  expect_error(
    RunCellQC(
      srt,
      qc_metrics = character(),
      qc_features = list(bad = list(features = "HBA1", range = c(5, 0))),
      verbose = FALSE
    ),
    "lower <= upper"
  )
  expect_warning(
    prepared <- cellqc_prepare_qc_features(
      qc_features = list(ok = list(
        features = c("HBA1", "NOT_PRESENT"),
        range = c(0, 5)
      )),
      assay_features = rownames(srt)
    ),
    "ignored missing features"
  )
  expect_identical(prepared$ok$features, "HBA1")
})

test_that("db_scds cxds path does not run hybrid backend", {
  skip_if_not_installed("scds")
  skip_if_not_installed("SingleCellExperiment")

  counts <- Matrix::Matrix(
    c(
      1, 0, 2, 0,
      0, 3, 0, 4,
      5, 0, 0, 6
    ),
    nrow = 3,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)

  testthat::local_mocked_bindings(
    cxds = function(sce, ...) {
      sce[["cxds_score"]] <- seq_len(ncol(sce))
      sce
    },
    cxds_bcds_hybrid = function(...) {
      stop("hybrid backend should not be called")
    },
    .package = "scds"
  )

  out <- suppressWarnings(
    db_scds(
      srt,
      method = "cxds",
      db_rate = 0.5,
      data_type = "raw_counts",
      verbose = FALSE
    )
  )

  expect_true("db.scds_cxds_score" %in% colnames(out[[]]))
  expect_true("db.scds_cxds_class" %in% colnames(out[[]]))
  expect_false("db.scds_hybrid_score" %in% colnames(out[[]]))
})

test_that("db_Scrublet passes integer PCA component arguments", {
  counts <- Matrix::Matrix(
    c(
      1, 0, 2, 0,
      0, 3, 0, 4,
      5, 0, 0, 6
    ),
    nrow = 3,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)

  observed <- new.env(parent = emptyenv())
  scrub <- new.env(parent = emptyenv())
  scrub$threshold_ <- 0.5
  scrub$scrub_doublets <- function(...) {
    observed$args <- list(...)
    list(rep(0.1, ncol(srt)), rep(FALSE, ncol(srt)))
  }
  scrublet <- list(
    Scrublet = function(...) {
      scrub
    }
  )

  testthat::local_mocked_bindings(
    PrepareEnv = function(...) invisible(NULL),
    check_python = function(...) invisible(TRUE),
    CheckDataType = function(...) "raw_counts",
    .package = "scop"
  )
  testthat::local_mocked_bindings(
    import = function(name, ...) {
      expect_identical(name, "scrublet")
      scrublet
    },
    py_has_attr = function(x, name) TRUE,
    py_to_r = function(x) x,
    .package = "reticulate"
  )

  out <- db_Scrublet(
    srt,
    db_rate = 0.01,
    n_prin_comps = 10,
    min_counts = 1,
    min_cells = 1,
    data_type = "raw_counts",
    verbose = FALSE
  )

  expect_type(observed$args$n_prin_comps, "integer")
  expect_type(observed$args$min_counts, "integer")
  expect_type(observed$args$min_cells, "integer")
  expect_true("db.Scrublet_score" %in% colnames(out[[]]))
})
