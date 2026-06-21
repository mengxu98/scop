make_cside_seurat <- function() {
  counts <- matrix(
    c(
      5, 0, 3, 2,
      0, 4, 1, 6,
      2, 1, 0, 3
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:3), paste0("Spot", 1:4))
  )
  srt <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  )
  srt$condition <- c("A", "B", "A", "B")
  srt$region <- c("left", "middle", "right", "right")
  srt$col <- c(1, 2, 3, 4)
  srt$row <- c(1, 1, 2, 2)
  srt@tools$RCTD <- list(
    object = structure(list(id = "mock-rctd"), class = "RCTD"),
    backend_api = "SpatialRNA/Reference/create.RCTD/run.RCTD",
    parameters = list(rctd_mode = "full")
  )
  srt
}

mock_cside_result <- function() {
  list(
    de_results = list(
      all_gene_list = list(
        Alpha = data.frame(
          Z_score = c(3.0, 1.2),
          log_fc = c(0.8, 0.1),
          se = c(0.2, 0.3),
          paramindex_best = c(2, 2),
          p_val = c(0.001, 0.2),
          row.names = c("Gene1", "Gene2")
        )
      ),
      sig_gene_list = list(
        Alpha = data.frame(
          Z_score = 3.0,
          log_fc = 0.8,
          p_val = 0.001,
          row.names = "Gene1"
        )
      ),
      gene_fits = list(converged = TRUE)
    ),
    internal_vars_de = list(
      cell_types = "Alpha",
      params_to_test = 2
    )
  )
}

with_mock_cside <- function(funs, code) {
  testthat::local_mocked_bindings(
    check_r = function(packages, ...) {
      expect_identical(packages, "spacexr")
      invisible(TRUE)
    },
    get_namespace_fun = function(package, name) {
      expect_identical(package, "spacexr")
      if (!name %in% names(funs)) {
        stop("missing backend")
      }
      funs[[name]]
    }
  )
  force(code)
}

test_that("RunCSIDE validates Seurat and stored old RCTD inputs", {
  expect_error(
    RunCSIDE(matrix(1, nrow = 2), verbose = FALSE),
    "Seurat"
  )
  srt <- make_cside_seurat()
  srt@tools$RCTD <- NULL
  expect_error(
    RunCSIDE(srt, verbose = FALSE),
    "No RCTD object"
  )
  srt@tools$RCTD <- list(object = list(new_api = TRUE))
  expect_error(
    RunCSIDE(srt, verbose = FALSE),
    "old-api"
  )
})

test_that("RunCSIDE condition.by dispatches to single backend", {
  srt <- make_cside_seurat()
  fake_single <- function(myRCTD, explanatory.variable, cell_types, gene_threshold, ...) {
    expect_s3_class(myRCTD, "RCTD")
    expect_identical(cell_types, "Alpha")
    expect_equal(unname(explanatory.variable[colnames(srt)]), c(0, 1, 0, 1))
    expect_identical(gene_threshold, 0.01)
    mock_cside_result()
  }

  with_mock_cside(list("run.CSIDE.single" = fake_single), {
    out <- RunCSIDE(
      srt,
      condition.by = "condition",
      celltypes = "Alpha",
      gene_threshold = 0.01,
      verbose = FALSE
    )
  })

  expect_true("CSIDE" %in% names(out@tools))
  expect_named(
    out@tools$CSIDE$result_table,
    c(
      "feature", "celltype", "parameter", "logFC", "statistic",
      "p_value", "q_value", "significant", "method"
    )
  )
  expect_equal(out@tools$CSIDE$result_table$feature, c("Gene1", "Gene2"))
  expect_equal(out@tools$CSIDE$result_table$celltype, c("Alpha", "Alpha"))
  expect_equal(out@tools$CSIDE$result_table$significant, c(TRUE, FALSE))
  expect_equal(unique(out$CSIDE_n_sig), 1)
  expect_equal(unique(out$CSIDE_mode), "single")
})

test_that("RunCSIDE group.by builds region_list and dispatches to regions backend", {
  srt <- make_cside_seurat()
  fake_regions <- function(myRCTD, region_list, cell_types, ...) {
    expect_identical(list(...)$doublet_mode, FALSE)
    expect_named(region_list, c("left", "middle", "right"))
    expect_equal(region_list$left, "Spot1")
    expect_equal(region_list$middle, "Spot2")
    expect_equal(region_list$right, c("Spot3", "Spot4"))
    expect_null(cell_types)
    mock_cside_result()
  }

  with_mock_cside(list("run.CSIDE.regions" = fake_regions), {
    out <- RunCSIDE(
      srt,
      group.by = "region",
      prefix = "REG",
      tool_name = "REG_tool",
      verbose = FALSE
    )
  })

  expect_true("REG_tool" %in% names(out@tools))
  expect_equal(unique(out$REG_mode), "regions")
  expect_equal(out@tools$REG_tool$barcodes, paste0("Spot", 1:4))
})

test_that("RunCSIDE design dispatches to general backend and keeps explicit backend args", {
  srt <- make_cside_seurat()
  design <- data.frame(
    intercept = 1,
    condition = c(0, 1, 0, 1),
    row.names = colnames(srt)
  )
  fake_general <- function(
    myRCTD,
    X,
    barcodes,
    cell_types,
    cell_type_specific,
    params_to_test,
    ...
  ) {
    expect_equal(rownames(X), c("Spot2", "Spot4"))
    expect_equal(barcodes, c("Spot2", "Spot4"))
    expect_identical(cell_types, c("Alpha", "Beta"))
    expect_identical(cell_type_specific, c(FALSE, TRUE))
    expect_identical(params_to_test, 2L)
    mock_cside_result()
  }

  with_mock_cside(list("run.CSIDE" = fake_general), {
    out <- RunCSIDE(
      srt,
      design = design,
      barcodes = c("Spot2", "Spot4"),
      celltypes = c("Alpha", "Beta"),
      cell_type_specific = c(FALSE, TRUE),
      params_to_test = 2L,
      verbose = FALSE
    )
  })

  expect_equal(unique(out$CSIDE_mode), "general")
  expect_equal(out@tools$CSIDE$parameters$backend_args$params_to_test, 2L)
})

test_that("RunCSIDE intercept mode and backend availability errors are clear", {
  srt <- make_cside_seurat()
  fake_intercept <- function(myRCTD, barcodes, cell_types, ...) {
    expect_equal(barcodes, c("Spot1", "Spot3"))
    expect_null(cell_types)
    mock_cside_result()
  }

  with_mock_cside(list("run.CSIDE.intercept" = fake_intercept), {
    out <- RunCSIDE(
      srt,
      barcodes = c("Spot1", "Spot3"),
      verbose = FALSE
    )
  })
  expect_equal(unique(out$CSIDE_mode), "intercept")

  with_mock_cside(list(), {
    expect_error(
      RunCSIDE(srt, condition.by = "condition", verbose = FALSE),
      "dmcable/spacexr"
    )
  })
})
