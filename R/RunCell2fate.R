#' @title Run Cell2fate RNA velocity analysis
#'
#' @description
#' Run the official Python Cell2fate model on raw spliced and unspliced counts
#' from a `Seurat` object. Cell2fate uses an isolated Python 3.9 environment
#' because its upstream dependency stack is not compatible with the default
#' scvi-tools environment. Inputs, model files, posterior output, logs, and a
#' reproducible manifest are persisted under `result_dir`.
#'
#' A returned object can be passed back to the function for a matching resumed
#' run. Output previously recorded by the same `tool_name` and `prefix` is
#' replaced in the returned copy; unrelated metadata or tool entries are never
#' overwritten.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt A `Seurat` object containing raw spliced and unspliced counts.
#' @param result_dir Empty directory, or a directory owned by an earlier
#' `RunCell2fate()` run, used to persist inputs, model files, posterior output,
#' per-attempt logs, and the run manifest.
#' @param spliced_assay,unspliced_assay Assays containing raw spliced and
#' unspliced counts.
#' @param spliced_layer,unspliced_layer Raw-count layers in the corresponding
#' assays.
#' @param cluster.by Metadata column containing cell-state labels used for
#' Cell2fate training-data selection.
#' @param features Optional features to consider before Cell2fate filtering and
#' variable-gene selection.
#' @param remove_clusters Optional cluster labels to remove before training.
#' @param cells_per_cluster Maximum cells retained per cluster. Use `NULL` to
#' retain every cell. Cells excluded by this sampling receive `NA` posterior
#' values and `FALSE` in the generated `<prefix>_selected` metadata column.
#' @param min_shared_counts Minimum total shared spliced and unspliced counts
#' required for a gene.
#' @param n_var_genes Number of variable genes retained for model fitting.
#' @param n_modules Number of Cell2fate modules. If `NULL`, the upstream
#' `get_max_modules()` heuristic is used.
#' @param model_params Named arguments passed to
#' `Cell2fate_DynamicalModel()`.
#' @param train_params Named arguments passed to the model `train()` method.
#' @param posterior_params Named arguments passed to
#' `export_posterior()` through its `sample_kwargs` argument.
#' @param seed Random seed used by Python, NumPy, PyTorch, and scvi-tools.
#' @param envname Name of the isolated Cell2fate environment. If `NULL`,
#' `"cell2fate_env"` is used.
#' @param resume Reuse a completed run only when its input fingerprint,
#' parameters, and artifact hashes match.
#' @param overwrite Permit replacement of incompatible artifacts in an owned
#' `result_dir`.
#' @param prefix Prefix used for Cell2fate metadata columns. Existing columns
#'   without matching Cell2fate provenance are rejected.
#' @param tool_name Name of the `srt@tools` result entry. An existing unrelated
#'   entry is rejected.
#' @param store_velocity Whether to write the dense posterior velocity matrix
#' to CSV and read it into `srt@tools`. The posterior `.h5ad`, including its
#' velocity layer, is always retained on disk.
#'
#' @return A `Seurat` object with Cell2fate time, uncertainty, module activation,
#' module-state, and training-cell-selection metadata. Cells not selected for
#' training have `NA` posterior values. Detailed provenance and optional
#' velocity values are stored in `srt@tools[[tool_name]]`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- RunCell2fate(
#'   pancreas_sub,
#'   result_dir = "pancreas_cell2fate",
#'   cluster.by = "SubCellType",
#'   n_modules = 10
#' )
#' FeatureDimPlot(pancreas_sub, "Cell2fate_time")
#' }
RunCell2fate <- function(
  srt,
  result_dir,
  spliced_assay = "spliced",
  unspliced_assay = "unspliced",
  spliced_layer = "counts",
  unspliced_layer = "counts",
  cluster.by,
  features = NULL,
  remove_clusters = NULL,
  cells_per_cluster = 100L,
  min_shared_counts = 10L,
  n_var_genes = 2000L,
  n_modules = NULL,
  model_params = list(),
  train_params = list(
    max_epochs = 500L,
    batch_size = 1000L,
    train_size = 1,
    lr = 0.01,
    accelerator = "auto"
  ),
  posterior_params = list(
    num_samples = 30L,
    batch_size = NULL,
    use_gpu = FALSE,
    return_samples = FALSE
  ),
  seed = 1L,
  envname = NULL,
  resume = TRUE,
  overwrite = FALSE,
  prefix = "Cell2fate",
  tool_name = "Cell2fate",
  store_velocity = FALSE,
  verbose = TRUE
) {
  validate_scalar_flag(resume, "resume")
  validate_scalar_flag(overwrite, "overwrite")
  validate_scalar_flag(store_velocity, "store_velocity")
  validate_scalar_string(result_dir, "result_dir")
  validate_scalar_string(cluster.by, "cluster.by")
  validate_scalar_string(prefix, "prefix")
  validate_scalar_string(tool_name, "tool_name")
  validate_named_list(model_params, "model_params")
  validate_named_list(train_params, "train_params")
  validate_named_list(posterior_params, "posterior_params")
  min_shared_counts <- validate_scalar_integer(
    min_shared_counts,
    "min_shared_counts",
    minimum = 0L,
    message = "must be a non-negative integer"
  )
  n_var_genes <- validate_scalar_integer(n_var_genes, "n_var_genes")
  seed <- validate_scalar_integer(
    seed,
    "seed",
    minimum = 0L,
    message = "must be a non-negative integer"
  )
  if (!is.null(cells_per_cluster)) {
    cells_per_cluster <- validate_scalar_integer(
      cells_per_cluster,
      "cells_per_cluster"
    )
  }
  if (!is.null(n_modules)) {
    n_modules <- validate_scalar_integer(n_modules, "n_modules")
  }
  if (!is.null(remove_clusters)) {
    remove_clusters <- unique(as.character(remove_clusters))
    remove_clusters <- remove_clusters[!is.na(remove_clusters) & nzchar(remove_clusters)]
  }

  prepared <- cell2fate_prepare_input(
    srt = srt,
    spliced_assay = spliced_assay,
    unspliced_assay = unspliced_assay,
    spliced_layer = spliced_layer,
    unspliced_layer = unspliced_layer,
    cluster.by = cluster.by,
    features = features
  )
  srt <- cell2fate_prepare_seurat_output(
    srt = srt,
    prefix = prefix,
    tool_name = tool_name
  )

  expanded_result_dir <- path.expand(result_dir)
  if (cell2fate_path_is_symlink(expanded_result_dir)) {
    log_message(
      "{.arg result_dir} must not be a symbolic link",
      message_type = "error"
    )
  }
  result_dir <- normalizePath(
    expanded_result_dir,
    winslash = "/",
    mustWork = FALSE
  )
  cell2fate_prepare_result_dir(result_dir)
  cell2fate_claim_result_dir(result_dir)
  result_lock <- runner_acquire_lock(
    file.path(result_dir, ".cell2fate.lock"),
    backend = "Cell2fate"
  )
  on.exit(runner_release_lock(result_lock), add = TRUE)
  cell2fate_assert_safe_result_paths(result_dir)
  request_id <- runner_request_id()

  python <- cell2fate_check_python(envname = envname, verbose = verbose)
  workdir <- tempfile("cell2fate_inputs_")
  dir.create(workdir, recursive = TRUE)
  on.exit(unlink(workdir, recursive = TRUE, force = TRUE), add = TRUE)

  input_path <- file.path(workdir, "input.h5ad")
  cell2fate_write_input(
    prepared = prepared,
    path = input_path,
    verbose = verbose
  )

  config <- list(
    result_dir = result_dir,
    input_path = normalizePath(input_path, winslash = "/", mustWork = TRUE),
    cluster_column = "cell2fate_cluster",
    remove_clusters = remove_clusters %||% character(),
    cells_per_cluster = cells_per_cluster,
    min_shared_counts = min_shared_counts,
    n_var_genes = min(n_var_genes, length(prepared$features)),
    n_modules = n_modules,
    model_params = model_params,
    train_params = train_params,
    posterior_params = posterior_params,
    seed = seed,
    resume = resume,
    overwrite = overwrite,
    prefix = prefix,
    store_velocity = store_velocity,
    request_id = request_id,
    lock_token = result_lock$token
  )
  config_path <- file.path(workdir, "config.json")
  runner_write_json(config, config_path)

  runner <- runner_script_path("cell2fate_runner.py", "Cell2fate")
  logs_dir <- file.path(result_dir, "logs")
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  attempt_dir <- tempfile("attempt_", tmpdir = logs_dir)
  dir.create(attempt_dir)
  stdout_path <- file.path(attempt_dir, "stdout.log")
  stderr_path <- file.path(attempt_dir, "stderr.log")
  cache_dir <- file.path(workdir, ".cache")
  dir.create(file.path(cache_dir, "matplotlib"), recursive = TRUE, showWarnings = FALSE)

  log_message(
    "Run {.pkg Cell2fate} with {.val {length(prepared$features)}} genes and {.val {ncol(prepared$spliced)}} cells",
    verbose = verbose
  )
  status <- runner_system2(
    command = python,
    args = c(shQuote(runner), "--config", shQuote(config_path)),
    env = c(
      PYTHONNOUSERSITE = "1",
      MPLCONFIGDIR = file.path(cache_dir, "matplotlib")
    ),
    stdout = stdout_path,
    stderr = stderr_path
  )
  if (!identical(status, 0L)) {
    runner_error(
      status,
      stdout_path,
      stderr_path,
      backend = "Cell2fate"
    )
  }

  cell2fate_assert_safe_result_paths(result_dir)
  files <- cell2fate_result_files(result_dir)
  required_names <- c("input", "cell_metadata", "posterior", "manifest")
  if (isTRUE(store_velocity)) {
    required_names <- c(required_names, "velocity")
  }
  required_files <- c(
    unlist(files[required_names]),
    file.path(files$model, c("model.pt", "adata.h5ad")),
    file.path(result_dir, c(".cell2fate.json", ".complete"))
  )
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0L) {
    log_message(
      "{.pkg Cell2fate} did not produce all required output artifacts",
      message_type = "error"
    )
  }

  manifest <- cell2fate_validate_manifest(
    runner_read_json(files$manifest),
    request_id = request_id,
    store_velocity = store_velocity
  )
  cell2fate_validate_artifacts(
    manifest = manifest,
    result_dir = result_dir
  )
  manifest_cells <- cell2fate_manifest_cells(manifest)
  result_features <- cell2fate_manifest_features(
    manifest,
    input_features = prepared$features
  )
  metadata <- runner_read_csv(
    files$cell_metadata,
    "posterior metadata",
    backend = "Cell2fate"
  )
  metadata <- cell2fate_align_cells(
    metadata,
    cells = colnames(srt),
    label = "posterior metadata"
  )

  velocity <- NULL
  if (isTRUE(store_velocity)) {
    velocity <- runner_read_numeric_csv(
      files$velocity,
      "velocity",
      backend = "Cell2fate"
    )
    velocity <- cell2fate_align_cells(
      velocity,
      cells = colnames(srt),
      label = "velocity"
    )
  }
  cell2fate_validate_result_tables(
    metadata = metadata,
    velocity = velocity,
    manifest_cells = manifest_cells,
    manifest_features = result_features,
    n_modules = manifest$n_modules,
    prefix = prefix
  )
  full_metadata <- cell2fate_expand_metadata(
    metadata,
    cells = colnames(srt),
    prefix = prefix
  )
  collisions <- intersect(colnames(full_metadata), colnames(srt@meta.data))
  if (length(collisions) > 0L) {
    log_message(
      "Cell2fate metadata columns already exist: {.val {collisions}}. Use a different {.arg prefix}.",
      message_type = "error"
    )
  }
  srt_out <- Seurat::AddMetaData(srt, metadata = full_metadata)
  stored_files <- files[
    c(
      "cell_metadata",
      "input",
      "posterior",
      "model",
      "manifest",
      if (isTRUE(store_velocity)) "velocity"
    )
  ]
  srt_out@tools[[tool_name]] <- list(
    method = "Cell2fate",
    result_type = "rna_velocity",
    cells = rownames(metadata),
    features = result_features,
    cell_metadata = metadata,
    velocity = velocity,
    manifest = manifest,
    paths = lapply(
      stored_files,
      normalizePath,
      winslash = "/",
      mustWork = FALSE
    ),
    parameters = config[setdiff(
      names(config),
      c("input_path", "result_dir", "request_id", "lock_token")
    )],
    provenance = list(
      producer = .cell2fate_producer,
      backend_id = "cell2fate",
      backend_repository = .cell2fate_repository,
      backend_commit = manifest$backend_commit
    )
  )

  log_message(
    "{.pkg Cell2fate} analysis completed",
    message_type = "success",
    verbose = verbose
  )
  srt_out
}

.cell2fate_repository <- "BayraktarLab/cell2fate"
.cell2fate_commit <- "c03d1ca0bb963f550001c6070d4986a61ec8456a"
.cell2fate_producer <- "RunCell2fate"
.cell2fate_schema_version <- 1L

cell2fate_result_dir_is_owned <- function(result_dir) {
  owner_path <- file.path(result_dir, ".cell2fate.json")
  if (
    !file.exists(owner_path) ||
      cell2fate_path_is_symlink(owner_path)
  ) {
    return(FALSE)
  }
  owner <- tryCatch(
    runner_read_json(owner_path),
    error = function(...) NULL
  )
  if (!is.list(owner)) {
    return(FALSE)
  }
  schema <- owner[["runner_schema_version"]]
  identical(owner[["producer"]], .cell2fate_producer) &&
    is.numeric(schema) &&
    length(schema) == 1L &&
    !is.na(schema) &&
    schema == .cell2fate_schema_version
}

cell2fate_prepare_result_dir <- function(result_dir) {
  cell2fate_assert_safe_result_paths(result_dir)
  if (dir.exists(result_dir)) {
    entries <- list.files(
      result_dir,
      all.files = TRUE,
      no.. = TRUE
    )
    if (length(entries) > 0L && !cell2fate_result_dir_is_owned(result_dir)) {
      log_message(
        "{.arg result_dir} must be empty or owned by an earlier {.fn RunCell2fate} run",
        message_type = "error"
      )
    }
    return(invisible(result_dir))
  }
  if (!dir.create(result_dir, recursive = TRUE)) {
    log_message(
      "Unable to create {.arg result_dir}: {.file {result_dir}}",
      message_type = "error"
    )
  }
  invisible(result_dir)
}

cell2fate_claim_result_dir <- function(result_dir) {
  if (!cell2fate_result_dir_is_owned(result_dir)) {
    entries <- list.files(
      result_dir,
      all.files = TRUE,
      no.. = TRUE
    )
    if (length(entries) > 0L) {
      log_message(
        "{.arg result_dir} changed before it could be claimed by {.fn RunCell2fate}",
        message_type = "error"
      )
    }
    runner_write_json(
      list(
        producer = .cell2fate_producer,
        runner_schema_version = .cell2fate_schema_version,
        backend_commit = .cell2fate_commit
      ),
      file.path(result_dir, ".cell2fate.json")
    )
  }
  invisible(result_dir)
}

cell2fate_path_is_symlink <- function(path) {
  link <- Sys.readlink(path)
  length(link) == 1L && !is.na(link) && nzchar(link)
}

cell2fate_assert_safe_result_paths <- function(result_dir) {
  paths <- file.path(
    result_dir,
    c(
      ".cell2fate.json",
      ".cell2fate.lock",
      ".complete",
      "inputs",
      file.path("inputs", "input.h5ad"),
      "model",
      file.path("model", "model.pt"),
      file.path("model", "adata.h5ad"),
      "posterior",
      file.path("posterior", "cell2fate_posterior.h5ad"),
      "tables",
      file.path("tables", "cell_metadata.csv"),
      file.path("tables", "velocity.csv"),
      "logs",
      "manifest.json"
    )
  )
  unsafe <- paths[vapply(paths, cell2fate_path_is_symlink, logical(1))]
  if (length(unsafe) > 0L) {
    log_message(
      "{.arg result_dir} contains symbolic links in paths managed by {.fn RunCell2fate}: {.file {unsafe}}",
      message_type = "error"
    )
  }
  invisible(result_dir)
}

cell2fate_validate_manifest <- function(
  manifest,
  request_id,
  store_velocity
) {
  schema <- if (is.list(manifest)) {
    manifest[["runner_schema_version"]]
  } else {
    NULL
  }
  artifacts <- if (is.list(manifest)) {
    manifest[["artifacts"]]
  } else {
    NULL
  }
  n_modules <- if (is.list(manifest)) {
    manifest[["n_modules"]]
  } else {
    NULL
  }
  expected_artifacts <- c(
    "inputs/input.h5ad",
    "model/model.pt",
    "model/adata.h5ad",
    "posterior/cell2fate_posterior.h5ad",
    "tables/cell_metadata.csv",
    if (isTRUE(store_velocity)) "tables/velocity.csv"
  )
  valid <- is.list(manifest) &&
    identical(manifest[["producer"]], .cell2fate_producer) &&
    is.numeric(schema) &&
    length(schema) == 1L &&
    !is.na(schema) &&
    schema == .cell2fate_schema_version &&
    identical(manifest[["status"]], "complete") &&
    identical(manifest[["backend_commit"]], .cell2fate_commit) &&
    identical(manifest[["request_id"]], request_id) &&
    is.list(artifacts) &&
    !is.null(names(artifacts)) &&
    all(nzchar(names(artifacts))) &&
    !anyDuplicated(names(artifacts)) &&
    setequal(names(artifacts), expected_artifacts) &&
    is.numeric(n_modules) &&
    length(n_modules) == 1L &&
    is.finite(n_modules) &&
    n_modules >= 1 &&
    n_modules == floor(n_modules)
  if (!valid) {
    log_message(
      "{.pkg Cell2fate} returned an invalid run manifest",
      message_type = "error"
    )
  }
  manifest$n_modules <- as.integer(n_modules)
  manifest
}

cell2fate_validate_artifacts <- function(manifest, result_dir) {
  artifacts <- manifest[["artifacts"]]
  sha256sum <- get0(
    "sha256sum",
    envir = asNamespace("tools"),
    mode = "function",
    inherits = FALSE
  )
  valid <- vapply(
    names(artifacts),
    function(name) {
      record <- artifacts[[name]]
      path <- file.path(result_dir, name)
      size <- if (is.list(record)) record[["size"]] else NULL
      sha256 <- if (is.list(record)) record[["sha256"]] else NULL
      if (
        !is.list(record) ||
          !setequal(names(record), c("size", "sha256")) ||
          !is.numeric(size) ||
          length(size) != 1L ||
          !is.finite(size) ||
          size < 1 ||
          size != floor(size) ||
          !is.character(sha256) ||
          length(sha256) != 1L ||
          is.na(sha256) ||
          !grepl("^[[:xdigit:]]{64}$", sha256) ||
          !file.exists(path)
      ) {
        return(FALSE)
      }
      observed_size <- unname(file.info(path)$size)
      if (
        length(observed_size) != 1L ||
          is.na(observed_size) ||
          as.numeric(observed_size) != as.numeric(size)
      ) {
        return(FALSE)
      }
      if (!is.null(sha256sum)) {
        observed_sha256 <- unname(sha256sum(path))
        if (!identical(tolower(observed_sha256), tolower(sha256))) {
          return(FALSE)
        }
      }
      TRUE
    },
    logical(1)
  )
  if (!length(valid) || !all(valid)) {
    log_message(
      "{.pkg Cell2fate} output artifacts failed integrity validation",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

cell2fate_manifest_features <- function(manifest, input_features = NULL) {
  features <- as.character(
    unlist(manifest[["features"]] %||% character(), use.names = FALSE)
  )
  if (
    length(features) < 2L ||
      anyNA(features) ||
      any(!nzchar(features)) ||
      anyDuplicated(features)
  ) {
    log_message(
      "{.pkg Cell2fate} returned invalid manifest features",
      message_type = "error"
    )
  }
  if (
    !is.null(input_features) &&
      length(setdiff(features, input_features)) > 0L
  ) {
    log_message(
      "{.pkg Cell2fate} returned features absent from the prepared input",
      message_type = "error"
    )
  }
  features
}

cell2fate_manifest_cells <- function(manifest) {
  cells <- as.character(
    unlist(manifest[["cells"]] %||% character(), use.names = FALSE)
  )
  if (
    length(cells) < 2L ||
      anyNA(cells) ||
      any(!nzchar(cells)) ||
      anyDuplicated(cells)
  ) {
    log_message(
      "{.pkg Cell2fate} returned invalid manifest cells",
      message_type = "error"
    )
  }
  cells
}

cell2fate_validate_result_tables <- function(
  metadata,
  velocity,
  manifest_cells,
  manifest_features,
  n_modules,
  prefix
) {
  metadata_names <- colnames(metadata)
  module_ids <- seq.int(0L, n_modules - 1L)
  expected_names <- paste0(
    prefix,
    "_",
    c(
      "time",
      "time_uncertainty",
      paste0("module_", rep(module_ids, each = 2L), "_", c("activation", "state"))
    )
  )
  if (!setequal(metadata_names, expected_names)) {
    log_message(
      "{.pkg Cell2fate} returned inconsistent posterior metadata",
      message_type = "error"
    )
  }
  continuous_names <- paste0(
    prefix,
    "_",
    c(
      "time",
      "time_uncertainty",
      paste0("module_", module_ids, "_activation")
    )
  )
  invalid_continuous <-
    any(!vapply(metadata[continuous_names], is.numeric, logical(1))) ||
    any(vapply(
      metadata[continuous_names],
      function(x) any(!is.finite(x)),
      logical(1)
    ))
  numeric_columns <- vapply(metadata, is.numeric, logical(1))
  invalid_numeric <- any(vapply(
    metadata[numeric_columns],
    function(x) any(!is.finite(x)),
    logical(1)
  ))
  invalid_values <- any(vapply(
    metadata,
    function(x) {
      anyNA(x) || (is.character(x) && any(!nzchar(x)))
    },
    logical(1)
  ))
  if (
    !identical(rownames(metadata), manifest_cells) ||
      anyDuplicated(colnames(metadata)) ||
      invalid_continuous ||
      invalid_numeric ||
      invalid_values
  ) {
    log_message(
      "{.pkg Cell2fate} returned inconsistent posterior metadata",
      message_type = "error"
    )
  }
  if (
    !is.null(velocity) &&
      (
        !identical(rownames(velocity), manifest_cells) ||
          !identical(colnames(velocity), manifest_features)
      )
  ) {
    log_message(
      "{.pkg Cell2fate} returned velocity dimensions inconsistent with its manifest",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

cell2fate_managed_metadata <- function(metadata_names, prefix) {
  prefix_with_separator <- paste0(prefix, "_")
  prefixed <- startsWith(metadata_names, prefix_with_separator)
  suffixes <- substring(
    metadata_names[prefixed],
    nchar(prefix_with_separator) + 1L
  )
  metadata_names[prefixed][
    suffixes %in% c("time", "time_uncertainty", "selected") |
      grepl("^module_[0-9]+_(activation|state)$", suffixes)
  ]
}

cell2fate_prepare_seurat_output <- function(srt, prefix, tool_name) {
  existing_columns <- cell2fate_managed_metadata(
    colnames(srt@meta.data),
    prefix = prefix
  )
  existing_tool <- srt@tools[[tool_name]]
  if (is.null(existing_tool)) {
    if (length(existing_columns) > 0L) {
      log_message(
        "Cell2fate metadata columns already exist without matching provenance: {.val {existing_columns}}. Use a different {.arg prefix}.",
        message_type = "error"
      )
    }
    return(srt)
  }
  tool_is_owned <- is.list(existing_tool) &&
    identical(existing_tool[["method"]], "Cell2fate") &&
    is.list(existing_tool[["provenance"]]) &&
    identical(
      existing_tool[["provenance"]][["producer"]],
      .cell2fate_producer
    )
  if (!tool_is_owned) {
    log_message(
      "{.arg tool_name} {.val {tool_name}} is already used by an unrelated result. Use a different {.arg tool_name}.",
      message_type = "error"
    )
  }

  parameters <- existing_tool[["parameters"]]
  previous_prefix <- if (is.list(parameters)) {
    parameters[["prefix"]]
  } else {
    NULL
  }
  if (
    !is.character(previous_prefix) ||
      length(previous_prefix) != 1L ||
      is.na(previous_prefix) ||
      !identical(previous_prefix, prefix)
  ) {
    log_message(
      "Existing {.pkg Cell2fate} tool {.val {tool_name}} uses a different or unverifiable metadata prefix. Use a different {.arg tool_name}.",
      message_type = "error"
    )
  }
  recorded_columns <- colnames(existing_tool[["cell_metadata"]])
  if (is.null(recorded_columns)) {
    log_message(
      "Existing {.pkg Cell2fate} tool {.val {tool_name}} has invalid output provenance",
      message_type = "error"
    )
  }
  owned_columns <- unique(c(
    recorded_columns,
    paste0(prefix, "_selected")
  ))
  if (
    !setequal(
      recorded_columns,
      cell2fate_managed_metadata(recorded_columns, prefix = prefix)
    ) ||
      length(setdiff(existing_columns, owned_columns)) > 0L
  ) {
    log_message(
      "Existing Cell2fate metadata cannot be safely matched to {.arg tool_name} {.val {tool_name}}",
      message_type = "error"
    )
  }
  for (column in intersect(owned_columns, colnames(srt@meta.data))) {
    srt@meta.data[[column]] <- NULL
  }
  srt@tools[[tool_name]] <- NULL
  srt
}

cell2fate_prepare_input <- function(
  srt,
  spliced_assay,
  unspliced_assay,
  spliced_layer,
  unspliced_layer,
  cluster.by,
  features
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  if (anyDuplicated(colnames(srt))) {
    log_message("{.arg srt} cell names must be unique", message_type = "error")
  }
  if (!cluster.by %in% colnames(srt@meta.data)) {
    log_message(
      "{.arg cluster.by} {.val {cluster.by}} was not found in metadata",
      message_type = "error"
    )
  }
  clusters <- as.character(srt@meta.data[[cluster.by]])
  if (anyNA(clusters) || any(!nzchar(trimws(clusters)))) {
    log_message(
      "{.arg cluster.by} contains missing or empty labels",
      message_type = "error"
    )
  }

  spliced <- cell2fate_counts(
    srt,
    assay = spliced_assay,
    layer = spliced_layer,
    label = "spliced"
  )
  unspliced <- cell2fate_counts(
    srt,
    assay = unspliced_assay,
    layer = unspliced_layer,
    label = "unspliced"
  )
  shared <- intersect(rownames(spliced), rownames(unspliced))
  if (!is.null(features)) {
    features <- unique(as.character(features))
    missing_features <- setdiff(features, shared)
    if (length(missing_features) > 0L) {
      log_message(
        "{.arg features} contains genes absent from the two velocity assays: {.val {missing_features}}",
        message_type = "error"
      )
    }
    shared <- features
  }
  if (length(shared) < 2L) {
    log_message(
      "Cell2fate requires at least two shared genes",
      message_type = "error"
    )
  }
  cells <- colnames(srt)
  if (!setequal(colnames(spliced), cells) || !setequal(colnames(unspliced), cells)) {
    log_message(
      "Spliced and unspliced counts must contain exactly the Seurat cells",
      message_type = "error"
    )
  }
  spliced <- spliced[shared, cells, drop = FALSE]
  unspliced <- unspliced[shared, cells, drop = FALSE]
  if (ncol(spliced) < 2L) {
    log_message("Cell2fate requires at least two cells", message_type = "error")
  }

  list(
    spliced = spliced,
    unspliced = unspliced,
    features = shared,
    clusters = stats::setNames(clusters, cells)
  )
}

cell2fate_counts <- function(srt, assay, layer, label) {
  validate_scalar_string(assay, paste0(label, "_assay"))
  validate_scalar_string(layer, paste0(label, "_layer"))
  if (!assay %in% names(srt@assays)) {
    log_message(
      "{.arg {label}_assay} {.val {assay}} was not found",
      message_type = "error"
    )
  }
  layer_names <- SeuratObject::Layers(srt[[assay]], search = layer)
  if (length(layer_names) != 1L) {
    log_message(
      "{.arg {label}_layer} must identify one joined layer in assay {.val {assay}}",
      message_type = "error"
    )
  }
  counts <- SeuratObject::LayerData(
    srt,
    assay = assay,
    layer = layer_names[[1]]
  )
  values <- if (inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
  if (
    any(!is.finite(values)) ||
      any(values < 0) ||
      any(abs(values - round(values)) > sqrt(.Machine$double.eps))
  ) {
    log_message(
      "{.arg {label}_layer} must contain finite, non-negative raw integer counts",
      message_type = "error"
    )
  }
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    log_message(
      "{.arg {label}_layer} must have feature and cell names",
      message_type = "error"
    )
  }
  counts
}

cell2fate_write_input <- function(prepared, path, verbose = TRUE) {
  input <- Seurat::CreateSeuratObject(
    counts = prepared$spliced,
    assay = "spliced",
    meta.data = data.frame(
      cell2fate_cluster = unname(prepared$clusters),
      row.names = names(prepared$clusters),
      check.names = FALSE
    )
  )
  input[["unspliced"]] <- SeuratObject::CreateAssayObject(
    counts = prepared$unspliced
  )

  # Reuse the isolated runtime selected above without re-entering PrepareEnv.
  old_options <- options(scop_skip_python_prepare = TRUE)
  on.exit(options(old_options), add = TRUE)
  srt_to_h5ad(
    input,
    path = path,
    features = prepared$features,
    assay_x = "spliced",
    layer_x = "counts",
    assay_y = "unspliced",
    layer_y = "counts",
    reductions = character(),
    graphs = character(),
    neighbors = character(),
    overwrite = TRUE,
    verbose = verbose
  )
}

cell2fate_check_python <- function(envname = NULL, verbose = TRUE) {
  target_env <- envname %||% "cell2fate_env"
  PrepareEnv(
    envname = target_env,
    version = "3.9-1",
    modules = "cell2fate",
    verbose = verbose
  )
  base_runtime_ok <- check_python(
    c("setuptools<81", "jaxlib==0.4.10"),
    envname = target_env,
    pip = FALSE,
    verbose = FALSE
  )
  runtime_ok <- check_python(
    c(
      paste0(
        "git+https://github.com/",
        .cell2fate_repository,
        ".git@",
        .cell2fate_commit
      ),
      "scvi-tools==0.16.1",
      "anndata==0.8.0",
      "scanpy==1.9.1",
      "scvelo==0.2.4",
      "torch==1.11.0",
      "jax==0.4.10"
    ),
    envname = target_env,
    pip = TRUE,
    verbose = FALSE
  )
  if (!isTRUE(base_runtime_ok) || !isTRUE(runtime_ok)) {
    log_message(
      "The isolated Cell2fate environment is incomplete",
      message_type = "error"
    )
  }
  cache <- getOption("scop_env_cache", default = NULL)
  python <- cache[["python"]] %||% tryCatch(
    conda_python(
      envname = target_env,
      conda = resolve_conda("auto")
    ),
    error = function(...) NULL
  )
  if (is.null(python) || !file.exists(python)) {
    log_message(
      "Unable to resolve the Cell2fate Python executable",
      message_type = "error"
    )
  }
  normalizePath(python, winslash = "/", mustWork = TRUE)
}

cell2fate_result_files <- function(result_dir) {
  list(
    input = file.path(result_dir, "inputs", "input.h5ad"),
    cell_metadata = file.path(result_dir, "tables", "cell_metadata.csv"),
    velocity = file.path(result_dir, "tables", "velocity.csv"),
    posterior = file.path(result_dir, "posterior", "cell2fate_posterior.h5ad"),
    model = file.path(result_dir, "model"),
    manifest = file.path(result_dir, "manifest.json")
  )
}

cell2fate_align_cells <- function(x, cells, label) {
  if (is.null(rownames(x)) || anyDuplicated(rownames(x))) {
    log_message(
      "{.pkg Cell2fate} {.val {label}} must have unique cell names",
      message_type = "error"
    )
  }
  unknown <- setdiff(rownames(x), cells)
  if (length(unknown) > 0L || nrow(x) == 0L) {
    log_message(
      "{.pkg Cell2fate} {.val {label}} contains invalid cells",
      message_type = "error"
    )
  }
  x
}

cell2fate_expand_metadata <- function(metadata, cells, prefix) {
  output <- data.frame(row.names = cells)
  for (column in colnames(metadata)) {
    if (is.numeric(metadata[[column]])) {
      values <- rep(NA_real_, length(cells))
    } else {
      values <- rep(NA_character_, length(cells))
    }
    names(values) <- cells
    values[rownames(metadata)] <- metadata[[column]]
    output[[column]] <- unname(values)
  }
  output[[paste0(prefix, "_selected")]] <- cells %in% rownames(metadata)
  output
}
