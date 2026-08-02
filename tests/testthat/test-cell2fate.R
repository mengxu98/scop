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
  expect_identical(
    unname(req$packages[["setuptools"]]),
    "setuptools<81"
  )
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
  expect_error(
    env_requirements(version = "3.9-1"),
    "only supported"
  )
  expect_error(
    env_requirements(version = "3.9-1", modules = "scanpy"),
    "only supported"
  )
})

test_that("cell2fate environment check uses the pinned GitHub backend", {
  checked <- list()
  prepared <- NULL
  python <- Sys.which("python3")
  skip_if(!nzchar(python))
  srt <- make_cell2fate_srt()
  result_dir <- tempfile("cell2fate check ")

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
    resolve_conda = function(...) "mamba",
    runner_script_path = function(...) "cell2fate.py",
    runner_request_id = function() "test-request",
    runner_write_json = function(...) invisible(NULL),
    cell2fate_write_input = function(prepared, path, ...) {
      file.create(path)
      invisible(path)
    },
    runner_system2 = function(command, args, env, stdout, stderr) {
      files <- list(
        input = file.path(result_dir, "inputs", "input.h5ad"),
        cell_metadata = file.path(result_dir, "tables", "cell_metadata.csv"),
        posterior = file.path(result_dir, "posterior", "cell2fate_posterior.h5ad"),
        model = file.path(result_dir, "model"),
        manifest = file.path(result_dir, "manifest.json")
      )
      dir.create(dirname(files$input), recursive = TRUE, showWarnings = FALSE)
      dir.create(dirname(files$cell_metadata), recursive = TRUE, showWarnings = FALSE)
      dir.create(dirname(files$posterior), recursive = TRUE, showWarnings = FALSE)
      dir.create(files$model, recursive = TRUE, showWarnings = FALSE)
      writeLines("model", file.path(files$model, "model.pt"))
      writeLines("model data", file.path(files$model, "adata.h5ad"))
      writeLines("input", files$input)
      metadata <- data.frame(
        Cell2fate_time = seq_len(8),
        Cell2fate_time_uncertainty = seq(0.1, 0.8, 0.1),
        Cell2fate_module_0_activation = seq(0.2, 0.9, 0.1),
        Cell2fate_module_0_state = rep(c("ON", "OFF"), 4),
        row.names = paste0("Cell", 1:8),
        check.names = FALSE
      )
      utils::write.csv(metadata, files$cell_metadata)
      writeLines("posterior", files$posterior)
      artifact_paths <- c(
        "inputs/input.h5ad" = files$input,
        "model/model.pt" = file.path(files$model, "model.pt"),
        "model/adata.h5ad" = file.path(files$model, "adata.h5ad"),
        "posterior/cell2fate_posterior.h5ad" = files$posterior,
        "tables/cell_metadata.csv" = files$cell_metadata
      )
      sha256sum <- get0(
        "sha256sum",
        envir = asNamespace("tools"),
        mode = "function",
        inherits = FALSE
      )
      artifacts <- lapply(
        artifact_paths,
        function(path) {
          list(
            size = unname(file.info(path)$size),
            sha256 = if (is.null(sha256sum)) {
              paste(rep("0", 64L), collapse = "")
            } else {
              unname(sha256sum(path))
            }
          )
        }
      )
      manifest <- list(
        producer = "RunCell2fate",
        runner_schema_version = 1L,
        backend_commit = "c03d1ca0bb963f550001c6070d4986a61ec8456a",
        status = "complete",
        request_id = "test-request",
        features = paste0("Gene", 1:3),
        cells = paste0("Cell", 1:8),
        artifacts = artifacts,
        n_modules = 1L,
        versions = list(cell2fate = "0.1a0")
      )
      to_json <- getExportedValue(
        "thisutils",
        "get_namespace_fun"
      )("jsonlite", "toJSON")
      writeLines(
        as.character(to_json(
          manifest,
          auto_unbox = TRUE,
          null = "null",
          digits = NA
        )),
        files$manifest
      )
      writeLines(
        paste0(
          '{"producer":"RunCell2fate","runner_schema_version":1,',
          '"backend_commit":"c03d1ca0bb963f550001c6070d4986a61ec8456a"}'
        ),
        file.path(result_dir, ".cell2fate.json")
      )
      file.create(file.path(result_dir, ".complete"))
      file.create(stdout)
      file.create(stderr)
      0L
    }
  )

  out <- RunCell2fate(
    srt,
    result_dir = result_dir,
    cluster.by = "celltype",
    envname = "cell2fate_test",
    n_modules = 1L,
    verbose = FALSE
  )

  expect_identical(prepared$version, "3.9-1")
  expect_identical(prepared$modules, "cell2fate")
  expect_identical(prepared$envname, "cell2fate_test")
  expect_setequal(
    checked[[1]]$packages,
    c("setuptools<81", "jaxlib==0.4.10")
  )
  expect_false(checked[[1]]$pip)
  expect_match(checked[[2]]$packages[[1]], "^git\\+https://")
  expect_match(checked[[2]]$packages[[1]], "cell2fate.git@c03d1ca0")
  expect_true(checked[[2]]$pip)
  expect_equal(unname(out$Cell2fate_time), seq_len(8))
})

test_that("Cell2fate subprocess environment supports paths with spaces", {
  rscript <- file.path(R.home("bin"), "Rscript")
  value <- tempfile("cell2fate path with spaces ")
  dir.create(value)
  output <- tempfile()
  error_output <- tempfile()
  old_value <- Sys.getenv("CELL2FATE_SPACE_TEST", unset = NA_character_)
  status <- getFromNamespace("runner_system2", "scop")(
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
    PrepareEnv = function(...) invisible(NULL),
    check_python = function(...) TRUE,
    conda_python = function(...) Sys.which("python3"),
    resolve_conda = function(...) "mamba",
    runner_script_path = function(...) "cell2fate.py",
    cell2fate_write_input = function(prepared, path, ...) {
      file.create(path)
      invisible(path)
    },
    runner_request_id = function() "test-request",
    runner_write_json = function(...) invisible(NULL),
    runner_system2 = function(command, args, env, stdout, stderr) {
      files <- list(
        input = file.path(result_dir, "inputs", "input.h5ad"),
        cell_metadata = file.path(result_dir, "tables", "cell_metadata.csv"),
        velocity = file.path(result_dir, "tables", "velocity.csv"),
        posterior = file.path(result_dir, "posterior", "cell2fate_posterior.h5ad"),
        model = file.path(result_dir, "model"),
        manifest = file.path(result_dir, "manifest.json")
      )
      dir.create(dirname(files$input), recursive = TRUE, showWarnings = FALSE)
      dir.create(dirname(files$cell_metadata), recursive = TRUE, showWarnings = FALSE)
      dir.create(dirname(files$posterior), recursive = TRUE, showWarnings = FALSE)
      dir.create(files$model, recursive = TRUE, showWarnings = FALSE)
      writeLines("model", file.path(files$model, "model.pt"))
      writeLines("model data", file.path(files$model, "adata.h5ad"))
      writeLines("input", files$input)
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
      writeLines("posterior", files$posterior)
      artifact_paths <- c(
        "inputs/input.h5ad" = files$input,
        "model/model.pt" = file.path(files$model, "model.pt"),
        "model/adata.h5ad" = file.path(files$model, "adata.h5ad"),
        "posterior/cell2fate_posterior.h5ad" = files$posterior,
        "tables/cell_metadata.csv" = files$cell_metadata,
        "tables/velocity.csv" = files$velocity
      )
      sha256sum <- get0(
        "sha256sum",
        envir = asNamespace("tools"),
        mode = "function",
        inherits = FALSE
      )
      artifacts <- lapply(
        artifact_paths,
        function(path) {
          list(
            size = unname(file.info(path)$size),
            sha256 = if (is.null(sha256sum)) {
              paste(rep("0", 64L), collapse = "")
            } else {
              unname(sha256sum(path))
            }
          )
        }
      )
      manifest <- list(
        producer = "RunCell2fate",
        runner_schema_version = 1L,
        backend_commit = "c03d1ca0bb963f550001c6070d4986a61ec8456a",
        status = "complete",
        request_id = "test-request",
        features = paste0("Gene", 1:3),
        cells = paste0("Cell", 1:8),
        artifacts = artifacts,
        n_modules = 1L,
        versions = list(cell2fate = "0.1a0")
      )
      to_json <- getExportedValue(
        "thisutils",
        "get_namespace_fun"
      )("jsonlite", "toJSON")
      writeLines(
        as.character(to_json(
          manifest,
          auto_unbox = TRUE,
          null = "null",
          digits = NA
        )),
        files$manifest
      )
      writeLines(
        paste0(
          '{"producer":"RunCell2fate","runner_schema_version":1,',
          '"backend_commit":"c03d1ca0bb963f550001c6070d4986a61ec8456a"}'
        ),
        file.path(result_dir, ".cell2fate.json")
      )
      file.create(file.path(result_dir, ".complete"))
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
  expect_identical(out@tools$Cell2fate$features, paste0("Gene", 1:3))
  expect_identical(out@tools$Cell2fate$manifest$status, "complete")
  expect_identical(
    out@tools$Cell2fate$provenance$backend_commit,
    "c03d1ca0bb963f550001c6070d4986a61ec8456a"
  )
  expect_false(file.exists(file.path(result_dir, ".cell2fate.lock")))
})

test_that("Cell2fate CSV readers preserve character cell names", {
  metadata_path <- tempfile(fileext = ".csv")
  velocity_path <- tempfile(fileext = ".csv")
  writeLines(
    c(
      ",Cell2fate_time,Cell2fate_module_0_state",
      "NA,0.5,Induction",
      "001,1.5,ON",
      "010,2.5,OFF"
    ),
    metadata_path
  )
  writeLines(
    c(
      ",Gene1,Gene2",
      "NA,0.5,1.5",
      "001,1.5,2.5",
      "010,3.5,4.5"
    ),
    velocity_path
  )

  metadata <- getFromNamespace("runner_read_csv", "scop")(
    metadata_path,
    "posterior metadata",
    backend = "Cell2fate"
  )
  velocity <- getFromNamespace("runner_read_numeric_csv", "scop")(
    velocity_path,
    "velocity",
    backend = "Cell2fate"
  )

  expect_identical(rownames(metadata), c("NA", "001", "010"))
  expect_identical(metadata$Cell2fate_time, c(0.5, 1.5, 2.5))
  expect_identical(rownames(velocity), c("NA", "001", "010"))
  expect_identical(unname(velocity[, "Gene1"]), c(0.5, 1.5, 3.5))
})

test_that("Cell2fate manifest features are always a character vector", {
  manifest <- list(features = as.list(c("Gene1", "Gene2")))
  features <- getFromNamespace("cell2fate_manifest_features", "scop")(
    manifest,
    input_features = c("Gene1", "Gene2", "Gene3")
  )

  expect_identical(features, c("Gene1", "Gene2"))
  expect_error(
    getFromNamespace("cell2fate_manifest_features", "scop")(
      manifest,
      input_features = c("Gene1", "Gene3")
    ),
    "absent from the prepared input"
  )
})

test_that("Cell2fate validates artifact sizes and checksums", {
  result_dir <- tempfile("cell2fate artifact ")
  path <- file.path(result_dir, "inputs", "input.h5ad")
  dir.create(dirname(path), recursive = TRUE)
  writeBin(charToRaw("first"), path)
  sha256sum <- get0(
    "sha256sum",
    envir = asNamespace("tools"),
    mode = "function",
    inherits = FALSE
  )
  manifest <- list(
    artifacts = list(
      "inputs/input.h5ad" = list(
        size = unname(file.info(path)$size),
        sha256 = if (is.null(sha256sum)) {
          paste(rep("0", 64L), collapse = "")
        } else {
          unname(sha256sum(path))
        }
      )
    )
  )
  validate <- getFromNamespace("cell2fate_validate_artifacts", "scop")
  expect_invisible(validate(manifest, result_dir))

  writeBin(charToRaw("other"), path)
  if (!is.null(sha256sum)) {
    expect_error(
      validate(manifest, result_dir),
      "failed integrity validation"
    )
  }
  writeBin(charToRaw("longer"), path)
  expect_error(
    validate(manifest, result_dir),
    "failed integrity validation"
  )
})

test_that("Cell2fate rejects non-finite or inconsistent result tables", {
  validate <- getFromNamespace("cell2fate_validate_result_tables", "scop")
  metadata <- data.frame(
    Cell2fate_time = c(1, 2),
    Cell2fate_time_uncertainty = c(0.1, 0.2),
    Cell2fate_module_0_activation = c(0.2, 0.3),
    Cell2fate_module_0_state = c("ON", "OFF"),
    row.names = c("Cell1", "Cell2")
  )
  expect_invisible(validate(
    metadata,
    velocity = NULL,
    manifest_cells = c("Cell1", "Cell2"),
    manifest_features = c("Gene1", "Gene2"),
    n_modules = 1L,
    prefix = "Cell2fate"
  ))

  invalid <- metadata
  invalid$Cell2fate_time[[2L]] <- NaN
  expect_error(
    validate(
      invalid,
      velocity = NULL,
      manifest_cells = c("Cell1", "Cell2"),
      manifest_features = c("Gene1", "Gene2"),
      n_modules = 1L,
      prefix = "Cell2fate"
    ),
    "inconsistent posterior metadata"
  )
  character_time <- metadata
  character_time$Cell2fate_time <- c("early", "late")
  expect_error(
    validate(
      character_time,
      velocity = NULL,
      manifest_cells = c("Cell1", "Cell2"),
      manifest_features = c("Gene1", "Gene2"),
      n_modules = 1L,
      prefix = "Cell2fate"
    ),
    "inconsistent posterior metadata"
  )
  expect_error(
    validate(
      metadata,
      velocity = matrix(
        1,
        nrow = 2,
        ncol = 2,
        dimnames = list(c("Cell2", "Cell1"), c("Gene1", "Gene2"))
      ),
      manifest_cells = c("Cell1", "Cell2"),
      manifest_features = c("Gene1", "Gene2"),
      n_modules = 1L,
      prefix = "Cell2fate"
    ),
    "inconsistent with its manifest"
  )
  missing_state <- metadata
  missing_state$Cell2fate_module_0_state <- NULL
  expect_error(
    validate(
      missing_state,
      velocity = NULL,
      manifest_cells = c("Cell1", "Cell2"),
      manifest_features = c("Gene1", "Gene2"),
      n_modules = 1L,
      prefix = "Cell2fate"
    ),
    "inconsistent posterior metadata"
  )
})

test_that("Cell2fate recognizes an older backend owner for safe replacement", {
  result_dir <- tempfile("older cell2fate owner ")
  dir.create(result_dir)
  writeLines(
    paste0(
      '{"producer":"RunCell2fate","runner_schema_version":1,',
      '"backend_commit":"older-backend"}'
    ),
    file.path(result_dir, ".cell2fate.json")
  )

  expect_true(
    getFromNamespace("cell2fate_result_dir_is_owned", "scop")(result_dir)
  )
})

test_that("Cell2fate refuses a non-empty unowned result directory", {
  srt <- make_cell2fate_srt()
  result_dir <- tempfile("unowned cell2fate ")
  logs_dir <- file.path(result_dir, "logs")
  dir.create(logs_dir, recursive = TRUE)
  stdout_path <- file.path(logs_dir, "cell2fate_stdout.log")
  writeLines("user-owned log", stdout_path)
  python_checked <- FALSE

  testthat::local_mocked_bindings(
    .package = "scop",
    PrepareEnv = function(...) {
      python_checked <<- TRUE
      invisible(NULL)
    }
  )

  expect_error(
    RunCell2fate(
      srt,
      result_dir = result_dir,
      cluster.by = "celltype",
      verbose = FALSE
    ),
    "empty or owned"
  )
  expect_false(python_checked)
  expect_identical(readLines(stdout_path), "user-owned log")
})

test_that("Cell2fate rejects symbolic links in managed result paths", {
  # file.symlink() creates junctions on Windows, which Sys.readlink()
  # does not report, so the POSIX symlink scenario cannot be simulated.
  skip_on_os("windows")

  result_dir <- tempfile("cell2fate symlink ")
  outside <- tempfile("cell2fate outside ")
  dir.create(result_dir)
  dir.create(outside)
  linked <- suppressWarnings(
    file.symlink(outside, file.path(result_dir, "logs"))
  )
  skip_if(!isTRUE(linked), "Symbolic links are not available")

  expect_error(
    getFromNamespace("cell2fate_prepare_result_dir", "scop")(result_dir),
    "symbolic links"
  )
  expect_length(list.files(outside, all.files = TRUE, no.. = TRUE), 0L)
})

test_that("Cell2fate runner errors retain stderr when stdout is long", {
  stdout_path <- tempfile()
  stderr_path <- tempfile()
  writeLines(paste("stdout line", seq_len(40)), stdout_path)
  writeLines(
    c("Traceback (most recent call last):", "cell2fate traceback sentinel"),
    stderr_path
  )

  expect_error(
    getFromNamespace("runner_error", "scop")(
      2L,
      stdout_path,
      stderr_path,
      backend = "Cell2fate"
    ),
    "cell2fate traceback\\s+sentinel"
  )
})

test_that("Cell2fate runner treats a malformed resume manifest as a cache miss", {
  python <- unname(Sys.which(c("python3", "python")))
  python <- python[nzchar(python)][1]
  skip_if(is.na(python), "Python is not available")
  python_version <- suppressWarnings(
    system2(python, "--version", stdout = TRUE, stderr = TRUE)
  )
  skip_if(
    length(python_version) == 0L || is.na(python_version[[1L]]),
    "Python version could not be determined"
  )
  version_match <- regmatches(
    python_version[[1L]],
    regexec("^Python ([0-9]+)\\.([0-9]+)", python_version[[1L]])
  )[[1]]
  skip_if(length(version_match) < 3L, "Python version could not be determined")
  python_major <- as.integer(version_match[[2L]])
  python_minor <- as.integer(version_match[[3L]])
  skip_if(
    python_major < 3L || (python_major == 3L && python_minor < 9L),
    "Python 3.9 or newer is required"
  )

  runner <- getFromNamespace("runner_script_path", "scop")(
    "cell2fate.py",
    "Cell2fate"
  )
  script <- tempfile(fileext = ".py")
  output <- tempfile()
  writeLines(
    c(
      "import runpy",
      "import json",
      "import os",
      "import sys",
      "import tempfile",
      "import types",
      "from concurrent.futures import ThreadPoolExecutor",
      "from pathlib import Path",
      "for name in ('anndata', 'numpy', 'pandas'):",
      "    sys.modules[name] = types.ModuleType(name)",
      "class FakeArray:",
      "    def __init__(self, shape, finite=True):",
      "        self.shape = shape",
      "        self.finite = finite",
      "    def __getitem__(self, key):",
      "        return FakeArray((1, 1), self.finite)",
      "    def __mul__(self, other):",
      "        return FakeArray((1, 1), self.finite and other.finite)",
      "    def __sub__(self, other):",
      "        return FakeArray((1, 1), self.finite and other.finite)",
      "numpy = sys.modules['numpy']",
      "numpy.asarray = lambda value, dtype=None: value",
      "numpy.isfinite = lambda value: types.SimpleNamespace(all=lambda: value.finite)",
      "module = runpy.run_path(sys.argv[1], run_name='cell2fate_test')",
      "model = types.SimpleNamespace(samples={'post_sample_means': {",
      "    'mu_expression': FakeArray((1, 1, 2), finite=False),",
      "    'beta_g': FakeArray((1, 1)),",
      "    'gamma_g': FakeArray((1, 1)),",
      "}})",
      "adata = types.SimpleNamespace(n_obs=1, n_vars=1, layers={})",
      "try:",
      "    module['_add_velocity_from_posterior'](adata, model)",
      "except RuntimeError as error:",
      "    assert 'non-finite' in str(error)",
      "else:",
      "    raise AssertionError('non-finite velocity was accepted')",
      "with tempfile.TemporaryDirectory() as temporary:",
      "    paths = module['_output_paths'](Path(temporary))",
      "    with module['_result_lock'](Path(temporary)):",
      "        try:",
      "            with module['_result_lock'](Path(temporary)):",
      "                pass",
      "        except RuntimeError as error:",
      "            assert 'Another Cell2fate run' in str(error)",
      "        else:",
      "            raise AssertionError('concurrent result lock was accepted')",
      "    foreign_lock = Path(temporary) / '.cell2fate.lock'",
      "    with module['_result_lock'](Path(temporary)):",
      "        foreign_lock.unlink()",
      "        foreign_lock.write_text('foreign-owner\\n', encoding='utf-8')",
      "    assert foreign_lock.read_text(encoding='utf-8').strip() == 'foreign-owner'",
      "    foreign_lock.unlink()",
      "    external_lock = Path(temporary) / '.cell2fate.lock'",
      "    external_lock.write_text('owner-token\\n', encoding='utf-8')",
      "    module['_assert_external_lock'](Path(temporary), 'owner-token')",
      "    try:",
      "        module['_assert_external_lock'](Path(temporary), 'other-token')",
      "    except RuntimeError as error:",
      "        assert 'owner changed' in str(error)",
      "    else:",
      "        raise AssertionError('foreign external lock was accepted')",
      "    stale_stage = Path(temporary) / '.cell2fate-stage-stale'",
      "    stale_stage.mkdir()",
      "    (stale_stage / 'partial').write_text('partial', encoding='utf-8')",
      "    module['_cleanup_stale_stages'](",
      "        Path(temporary), 'owner-token'",
      "    )",
      "    assert not stale_stage.exists()",
      "    if hasattr(os, 'symlink') and os.name == 'posix':",
      "        outside_stage = Path(temporary) / 'outside-stage'",
      "        outside_stage.mkdir()",
      "        outside_marker = outside_stage / 'keep'",
      "        outside_marker.write_text('keep', encoding='utf-8')",
      "        stage_link = Path(temporary) / '.cell2fate-stage-link'",
      "        try:",
      "            stage_link.symlink_to(outside_stage, target_is_directory=True)",
      "        except OSError:",
      "            pass",
      "        else:",
      "            module['_cleanup_stale_stages'](",
      "                Path(temporary), 'owner-token'",
      "            )",
      "            assert not stage_link.exists()",
      "            assert outside_marker.read_text(encoding='utf-8') == 'keep'",
      "    paths['model'].mkdir(parents=True)",
      "    old_marker = paths['model'] / 'old-result'",
      "    old_marker.write_text('keep', encoding='utf-8')",
      "    with module['_staged_output_paths'](Path(temporary)) as staged:",
      "        for key in ('inputs', 'model', 'posterior', 'tables'):",
      "            staged[key].mkdir(parents=True)",
      "        staged['manifest'].write_text('{}', encoding='utf-8')",
      "        staged['complete'].write_text('complete', encoding='utf-8')",
      "        external_lock.write_text('foreign-owner\\n', encoding='utf-8')",
      "        try:",
      "            module['_publish_staged_outputs'](",
      "                staged, paths, Path(temporary), 'owner-token'",
      "            )",
      "        except RuntimeError as error:",
      "            assert 'owner changed' in str(error)",
      "        else:",
      "            raise AssertionError('lost lock published staged output')",
      "    assert old_marker.read_text(encoding='utf-8') == 'keep'",
      "    old_marker.unlink()",
      "    paths['model'].rmdir()",
      "    external_lock.unlink()",
      "    if os.name == 'posix':",
      "        # os.replace() of the same target from multiple threads is",
      "        # atomic on POSIX but can fail with PermissionError when the",
      "        # target is transiently open on Windows.",
      "        concurrent_json = Path(temporary) / 'concurrent.json'",
      "        with ThreadPoolExecutor(max_workers=8) as executor:",
      "            list(executor.map(",
      "                lambda value: module['_write_json']({'value': value}, concurrent_json),",
      "                range(40),",
      "            ))",
      "        assert isinstance(json.loads(concurrent_json.read_text()), dict)",
      "    paths['model'].mkdir(parents=True)",
      "    paths['inputs'].mkdir(parents=True)",
      "    paths['posterior'].mkdir(parents=True)",
      "    paths['tables'].mkdir(parents=True)",
      "    paths['owner'].write_text('[]', encoding='utf-8')",
      "    assert module['_is_owned_result'](paths) is False",
      "    module['_write_json']({",
      "        'producer': module['PRODUCER'],",
      "        'runner_schema_version': module['RUNNER_SCHEMA_VERSION'],",
      "        'backend_commit': 'older-backend',",
      "    }, paths['owner'])",
      "    assert module['_is_owned_result'](paths) is True",
      "    (paths['model'] / 'model.pt').write_bytes(b'model')",
      "    (paths['inputs'] / 'input.h5ad').write_bytes(b'input')",
      "    (paths['model'] / 'adata.h5ad').write_bytes(b'adata')",
      "    (paths['posterior'] / 'cell2fate_posterior.h5ad').write_bytes(b'posterior')",
      "    (paths['tables'] / 'cell_metadata.csv').write_bytes(b'metadata')",
      "    paths['complete'].write_text('complete')",
      "    paths['manifest'].write_text('{', encoding='utf-8')",
      "    assert module['_can_resume'](",
      "        paths, 'fingerprint', {}, require_velocity=False",
      "    ) is False",
      "    module['_write_json']({",
      "        'producer': module['PRODUCER'],",
      "        'runner_schema_version': module['RUNNER_SCHEMA_VERSION'],",
      "        'backend_commit': module['BACKEND_COMMIT'],",
      "        'input_sha256': 'fingerprint',",
      "        'parameters': {},",
      "        'status': 'complete',",
      "        'artifacts': module['_artifact_records'](",
      "            paths, require_velocity=False",
      "        ),",
      "    }, paths['manifest'])",
      "    assert module['_can_resume'](",
      "        paths, 'fingerprint', {}, require_velocity=False",
      "    ) is True",
      "    strict_parameters = {'flag': False}",
      "    module['_write_json']({",
      "        'producer': module['PRODUCER'],",
      "        'runner_schema_version': module['RUNNER_SCHEMA_VERSION'],",
      "        'backend_commit': module['BACKEND_COMMIT'],",
      "        'input_sha256': 'fingerprint',",
      "        'parameters': strict_parameters,",
      "        'parameters_sha256': module['_parameter_fingerprint'](strict_parameters),",
      "        'status': 'complete',",
      "        'artifacts': module['_artifact_records'](",
      "            paths, require_velocity=False",
      "        ),",
      "    }, paths['manifest'])",
      "    assert module['_can_resume'](",
      "        paths, 'fingerprint', {'flag': 0}, require_velocity=False",
      "    ) is False",
      "    if hasattr(os, 'symlink') and os.name == 'posix':",
      "        outside = Path(temporary) / 'outside'",
      "        outside.mkdir()",
      "        try:",
      "            paths['logs'].symlink_to(outside, target_is_directory=True)",
      "        except OSError:",
      "            pass",
      "        else:",
      "            try:",
      "                module['_assert_safe_result_paths'](Path(temporary))",
      "            except RuntimeError as error:",
      "                assert 'symbolic link' in str(error)",
      "            else:",
      "                raise AssertionError('managed symlink was accepted')",
      "            paths['logs'].unlink()",
      "    (paths['model'] / 'model.pt').write_bytes(b'corrupt')",
      "    assert module['_can_resume'](",
      "        paths, 'fingerprint', strict_parameters, require_velocity=False",
      "    ) is False"
    ),
    script
  )
  status <- system2(
    python,
    c(shQuote(script), shQuote(runner)),
    stdout = output,
    stderr = output
  )

  expect_identical(
    status,
    0L,
    info = paste(readLines(output, warn = FALSE), collapse = "\n")
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

test_that("Cell2fate only replaces its own prior Seurat output", {
  srt <- make_cell2fate_srt()
  metadata <- data.frame(
    Cell2fate_time = seq_len(ncol(srt)),
    Cell2fate_module_0_activation = seq_len(ncol(srt)) / 10,
    Cell2fate_module_0_state = rep(c("ON", "OFF"), 4),
    row.names = colnames(srt),
    check.names = FALSE
  )
  srt <- Seurat::AddMetaData(srt, metadata)
  srt$Cell2fate_selected <- TRUE
  srt@tools$Cell2fate <- list(
    method = "Cell2fate",
    cell_metadata = metadata,
    parameters = list(prefix = "Cell2fate"),
    provenance = list(producer = "RunCell2fate")
  )

  cleaned <- getFromNamespace(
    "cell2fate_prepare_seurat_output",
    "scop"
  )(
    srt,
    prefix = "Cell2fate",
    tool_name = "Cell2fate"
  )
  expect_length(
    grep("^Cell2fate_", colnames(cleaned@meta.data), value = TRUE),
    0L
  )
  expect_false("Cell2fate" %in% names(cleaned@tools))
  expect_true("celltype" %in% colnames(cleaned@meta.data))

  unrelated <- make_cell2fate_srt()
  unrelated@tools$Cell2fate <- list(user = "keep")
  expect_error(
    getFromNamespace("cell2fate_prepare_seurat_output", "scop")(
      unrelated,
      prefix = "Cell2fate",
      tool_name = "Cell2fate"
    ),
    "unrelated result"
  )

  unowned_metadata <- make_cell2fate_srt()
  unowned_metadata$Cell2fate_time <- seq_len(ncol(unowned_metadata))
  expect_error(
    getFromNamespace("cell2fate_prepare_seurat_output", "scop")(
      unowned_metadata,
      prefix = "Cell2fate",
      tool_name = "Cell2fate"
    ),
    "without matching provenance"
  )
})

test_that("RunCell2fate leaves the input unchanged when Python fails", {
  srt <- make_cell2fate_srt()
  before <- srt
  result_dir <- tempfile("cell2fate failed ")

  testthat::local_mocked_bindings(
    .package = "scop",
    PrepareEnv = function(...) invisible(NULL),
    check_python = function(...) TRUE,
    conda_python = function(...) Sys.which("python3"),
    resolve_conda = function(...) "mamba",
    runner_script_path = function(...) "cell2fate.py",
    cell2fate_write_input = function(prepared, path, ...) {
      file.create(path)
      invisible(path)
    },
    runner_write_json = function(...) invisible(NULL),
    runner_system2 = function(command, args, env, stdout, stderr) {
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
    "[Ss]ynthetic\\s+Cell2fate\\s+backend\\s+failure"
  )
  expect_identical(before@meta.data, srt@meta.data)
  expect_identical(before@tools, srt@tools)
  expect_false(file.exists(file.path(result_dir, ".cell2fate.lock")))
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
