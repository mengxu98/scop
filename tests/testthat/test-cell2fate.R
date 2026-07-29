make_cell2fate_srt <- function() {
  cells <- paste0("Cell", 1:8)
  genes <- paste0("Gene", 1:6)
  spliced <- matrix(
    c(
      5, 4, 3, 2, 1, 2, 3, 4,
      2, 3, 4, 5, 2, 3, 4, 5,
      1, 1, 2, 2, 3, 3, 4, 4,
      6, 5, 4, 3, 2, 1, 2, 3,
      2, 2, 2, 3, 3, 3, 4, 4,
      4, 3, 2, 1, 4, 3, 2, 1
    ),
    nrow = 6,
    byrow = TRUE,
    dimnames = list(genes, cells)
  )
  unspliced <- matrix(
    c(
      2, 1, 1, 1, 1, 1, 2, 2,
      1, 1, 2, 2, 1, 1, 2, 2,
      0, 1, 0, 1, 1, 2, 1, 2,
      3, 2, 2, 1, 1, 1, 1, 2,
      1, 1, 1, 1, 2, 2, 2, 2,
      2, 2, 1, 1, 2, 2, 1, 1
    ),
    nrow = 6,
    byrow = TRUE,
    dimnames = list(genes, cells)
  )
  srt <- Seurat::CreateSeuratObject(
    methods::as(Matrix::Matrix(spliced + unspliced, sparse = TRUE), "dgCMatrix")
  )
  srt[["spliced"]] <- SeuratObject::CreateAssayObject(
    counts = methods::as(Matrix::Matrix(spliced, sparse = TRUE), "dgCMatrix")
  )
  srt[["unspliced"]] <- SeuratObject::CreateAssayObject(
    counts = methods::as(Matrix::Matrix(unspliced, sparse = TRUE), "dgCMatrix")
  )
  srt$celltype <- rep(c("Alpha", "Beta"), each = 4)
  srt
}

test_that("cell2fate environment is isolated and pinned to the upstream stack", {
  modules <- getFromNamespace("normalize_env_modules", "scop")("cell2fate")
  expect_identical(modules, "cell2fate")

  req <- env_requirements(modules = "cell2fate")
  expect_identical(req$python, "3.9-1")
  expect_identical(unname(req$packages[["jaxlib"]]), "jaxlib==0.4.10")
  expect_match(
    unname(req$packages[["cell2fate"]]),
    "BayraktarLab/cell2fate.git@c03d1ca0"
  )
  expect_error(
    getFromNamespace("normalize_env_modules", "scop")(
      c("cell2fate", "scanpy")
    ),
    "standalone"
  )
})

test_that("cell2fate environment check uses the pinned GitHub backend", {
  checked <- list()
  prepared <- NULL
  python <- Sys.which("python3")
  skip_if(!nzchar(python))

  testthat::local_mocked_bindings(
    .package = "scop",
    PrepareEnv = function(envname, version, modules, verbose) {
      prepared <<- list(
        envname = envname,
        version = version,
        modules = modules
      )
      invisible(NULL)
    },
    check_python = function(packages, envname, pip, verbose) {
      checked[[length(checked) + 1L]] <<- list(
        packages = packages,
        pip = pip
      )
      TRUE
    },
    conda_python = function(...) python,
    resolve_conda = function(...) "mamba"
  )

  found <- getFromNamespace("cell2fate_check_python", "scop")(
    envname = "cell2fate_test",
    verbose = FALSE
  )

  expect_identical(prepared$version, "3.9-1")
  expect_identical(prepared$modules, "cell2fate")
  expect_identical(checked[[1]]$packages, "jaxlib==0.4.10")
  expect_false(checked[[1]]$pip)
  expect_match(checked[[2]]$packages[[1]], "^git\\+https://")
  expect_match(checked[[2]]$packages[[1]], "cell2fate.git@c03d1ca0")
  expect_true(checked[[2]]$pip)
  expect_identical(found, normalizePath(python, winslash = "/", mustWork = TRUE))
})

test_that("Cell2fate subprocess environment supports paths with spaces", {
  rscript <- file.path(R.home("bin"), "Rscript")
  value <- tempfile("cell2fate path with spaces ")
  dir.create(value)
  output <- tempfile()
  error_output <- tempfile()
  old_value <- Sys.getenv("CELL2FATE_SPACE_TEST", unset = NA_character_)
  status <- getFromNamespace("cell2fate_run_system2", "scop")(
    command = rscript,
    args = c(
      "-e",
      shQuote("cat(Sys.getenv('CELL2FATE_SPACE_TEST'))")
    ),
    env = c(CELL2FATE_SPACE_TEST = value),
    stdout = output,
    stderr = error_output
  )

  expect_identical(status, 0L)
  expect_identical(paste(readLines(output, warn = FALSE), collapse = "\n"), value)
  expect_identical(
    Sys.getenv("CELL2FATE_SPACE_TEST", unset = NA_character_),
    old_value
  )
})

test_that("cell2fate input preparation aligns raw spliced and unspliced counts", {
  srt <- make_cell2fate_srt()
  prepared <- getFromNamespace("cell2fate_prepare_input", "scop")(
    srt = srt,
    spliced_assay = "spliced",
    unspliced_assay = "unspliced",
    spliced_layer = "counts",
    unspliced_layer = "counts",
    cluster.by = "celltype",
    features = c("Gene4", "Gene2", "Gene1")
  )

  expect_identical(prepared$features, c("Gene4", "Gene2", "Gene1"))
  expect_identical(rownames(prepared$spliced), prepared$features)
  expect_identical(dim(prepared$spliced), dim(prepared$unspliced))
  expect_identical(colnames(prepared$spliced), colnames(srt))
  expect_identical(unname(prepared$clusters), as.character(srt$celltype))
  expect_identical(names(prepared$clusters), colnames(srt))
})

test_that("cell2fate input conversion restores its temporary option", {
  srt <- make_cell2fate_srt()
  prepared <- getFromNamespace("cell2fate_prepare_input", "scop")(
    srt = srt,
    spliced_assay = "spliced",
    unspliced_assay = "unspliced",
    spliced_layer = "counts",
    unspliced_layer = "counts",
    cluster.by = "celltype",
    features = NULL
  )
  old_options <- options(scop_skip_python_prepare = NULL)
  on.exit(options(old_options), add = TRUE)

  testthat::local_mocked_bindings(
    .package = "scop",
    srt_to_h5ad = function(...) {
      expect_true(getOption("scop_skip_python_prepare"))
      invisible(NULL)
    }
  )

  getFromNamespace("cell2fate_write_input", "scop")(
    prepared,
    path = tempfile(fileext = ".h5ad"),
    verbose = FALSE
  )
  expect_null(getOption("scop_skip_python_prepare"))
})

test_that("cell2fate rejects normalized or misaligned velocity inputs", {
  srt <- make_cell2fate_srt()
  srt[["spliced"]] <- SeuratObject::SetAssayData(
    srt[["spliced"]],
    layer = "data",
    new.data = SeuratObject::LayerData(srt, assay = "spliced", layer = "counts") / 3
  )
  expect_error(
    getFromNamespace("cell2fate_prepare_input", "scop")(
      srt = srt,
      spliced_assay = "spliced",
      unspliced_assay = "unspliced",
      spliced_layer = "data",
      unspliced_layer = "counts",
      cluster.by = "celltype",
      features = NULL
    ),
    "raw integer counts"
  )

  srt$celltype[1] <- NA_character_
  expect_error(
    getFromNamespace("cell2fate_prepare_input", "scop")(
      srt = srt,
      spliced_assay = "spliced",
      unspliced_assay = "unspliced",
      spliced_layer = "counts",
      unspliced_layer = "counts",
      cluster.by = "celltype",
      features = NULL
    ),
    "missing"
  )
})

test_that("RunCell2fate maps posterior summaries back to Seurat", {
  srt <- make_cell2fate_srt()
  result_dir <- tempfile("cell2fate result with spaces ")
  dir.create(result_dir)

  testthat::local_mocked_bindings(
    .package = "scop",
    cell2fate_check_python = function(...) "python",
    cell2fate_runner_path = function() "cell2fate_runner.py",
    cell2fate_write_input = function(prepared, path, ...) {
      file.create(path)
      invisible(path)
    },
    cell2fate_write_json = function(...) invisible(NULL),
    cell2fate_run_system2 = function(command, args, env, stdout, stderr) {
      files <- getFromNamespace("cell2fate_result_files", "scop")(result_dir)
      dir.create(dirname(files$cell_metadata), recursive = TRUE, showWarnings = FALSE)
      dir.create(dirname(files$posterior), recursive = TRUE, showWarnings = FALSE)
      dir.create(files$model, recursive = TRUE, showWarnings = FALSE)
      metadata <- data.frame(
        Cell2fate_time = seq_len(8),
        Cell2fate_time_uncertainty = seq(0.1, 0.8, 0.1),
        Cell2fate_module_0_activation = seq(0.2, 0.9, 0.1),
        Cell2fate_module_0_state = rep(c("ON", "OFF"), 4),
        row.names = paste0("Cell", 1:8),
        check.names = FALSE
      )
      velocity <- matrix(
        seq_len(24) / 10,
        nrow = 8,
        dimnames = list(paste0("Cell", 1:8), paste0("Gene", 1:3))
      )
      utils::write.csv(metadata, files$cell_metadata)
      utils::write.csv(velocity, files$velocity)
      file.create(files$posterior)
      writeLines(
        '{"status":"complete","n_modules":1,"versions":{"cell2fate":"0.1a0"}}',
        files$manifest
      )
      file.create(stdout)
      file.create(stderr)
      0L
    }
  )

  out <- RunCell2fate(
    srt,
    result_dir = result_dir,
    cluster.by = "celltype",
    n_modules = 1L,
    train_params = list(max_epochs = 2L),
    store_velocity = TRUE,
    verbose = FALSE
  )

  expect_equal(unname(out$Cell2fate_time), seq_len(8))
  expect_equal(unname(out$Cell2fate_module_0_activation), seq(0.2, 0.9, 0.1))
  expect_true(all(out$Cell2fate_selected))
  expect_true("Cell2fate" %in% names(out@tools))
  expect_identical(out@tools$Cell2fate$method, "Cell2fate")
  expect_identical(out@tools$Cell2fate$result_type, "rna_velocity")
  expect_identical(out@tools$Cell2fate$cells, colnames(srt))
  expect_identical(dim(out@tools$Cell2fate$velocity), c(8L, 3L))
  expect_identical(out@tools$Cell2fate$manifest$status, "complete")
  expect_identical(
    out@tools$Cell2fate$provenance$backend_commit,
    "c03d1ca0bb963f550001c6070d4986a61ec8456a"
  )
})

test_that("Cell2fate marks cells excluded from backend training", {
  metadata <- data.frame(
    Cell2fate_time = c(1, 2),
    row.names = c("Cell2", "Cell4")
  )
  expanded <- getFromNamespace("cell2fate_expand_metadata", "scop")(
    metadata,
    cells = paste0("Cell", 1:4),
    prefix = "Cell2fate"
  )

  expect_equal(expanded$Cell2fate_time, c(NA, 1, NA, 2))
  expect_identical(
    expanded$Cell2fate_selected,
    c(FALSE, TRUE, FALSE, TRUE)
  )
})

test_that("RunCell2fate leaves the input unchanged when Python fails", {
  srt <- make_cell2fate_srt()
  before <- srt
  result_dir <- tempfile("cell2fate failed ")

  testthat::local_mocked_bindings(
    .package = "scop",
    cell2fate_check_python = function(...) "python",
    cell2fate_runner_path = function() "cell2fate_runner.py",
    cell2fate_write_input = function(prepared, path, ...) {
      file.create(path)
      invisible(path)
    },
    cell2fate_write_json = function(...) invisible(NULL),
    cell2fate_run_system2 = function(command, args, env, stdout, stderr) {
      writeLines("synthetic Cell2fate backend failure", stderr)
      file.create(stdout)
      2L
    }
  )

  expect_error(
    RunCell2fate(
      srt,
      result_dir = result_dir,
      cluster.by = "celltype",
      verbose = FALSE
    ),
    "[Ss]ynthetic Cell2fate backend failure"
  )
  expect_identical(before@meta.data, srt@meta.data)
  expect_identical(before@tools, srt@tools)
})

test_that("RunCell2fate validates parameter lists and output paths", {
  srt <- make_cell2fate_srt()
  expect_error(
    RunCell2fate(
      srt,
      result_dir = tempfile(),
      cluster.by = "celltype",
      train_params = 10
    ),
    "train_params"
  )
  expect_error(
    RunCell2fate(
      srt,
      result_dir = tempfile(),
      cluster.by = "missing"
    ),
    "cluster.by"
  )
})
