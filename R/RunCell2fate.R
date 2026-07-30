#' @title Run Cell2fate RNA velocity analysis
#'
#' @description
#' Run the official Python Cell2fate model on raw spliced and unspliced counts
#' from a `Seurat` object. Cell2fate uses an isolated Python 3.9 environment
#' because its upstream dependency stack is not compatible with SCOP's default
#' scvi-tools environment. Inputs, model files, posterior output, logs, and a
#' reproducible manifest are persisted under `result_dir`.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt A `Seurat` object containing raw spliced and unspliced counts.
#' @param result_dir Directory used to persist inputs, model files, posterior
#' output, tables, logs, and the run manifest.
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
#' @param resume Reuse a completed run only when its input fingerprint and
#' parameters match.
#' @param overwrite Permit replacement of incompatible existing artifacts.
#' @param prefix Prefix used for Cell2fate metadata columns.
#' @param tool_name Name of the `srt@tools` result entry.
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

  result_dir <- normalizePath(
    path.expand(result_dir),
    winslash = "/",
    mustWork = FALSE
  )
  if (!dir.exists(result_dir) && !dir.create(result_dir, recursive = TRUE)) {
    log_message(
      "Unable to create {.arg result_dir}: {.file {result_dir}}",
      message_type = "error"
    )
  }
  logs_dir <- file.path(result_dir, "logs")
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

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
    store_velocity = store_velocity
  )
  config_path <- file.path(workdir, "config.json")
  runner_write_json(config, config_path)

  runner <- runner_script_path("cell2fate_runner.py", "Cell2fate")
  stdout_path <- file.path(logs_dir, "cell2fate_stdout.log")
  stderr_path <- file.path(logs_dir, "cell2fate_stderr.log")
  cache_dir <- file.path(result_dir, ".cache")
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

  files <- cell2fate_result_files(result_dir)
  required_names <- c("cell_metadata", "posterior", "manifest")
  if (isTRUE(store_velocity)) {
    required_names <- c(required_names, "velocity")
  }
  required_files <- unlist(files[required_names])
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0L || !dir.exists(files$model)) {
    log_message(
      "{.pkg Cell2fate} did not produce all required output artifacts",
      message_type = "error"
    )
  }

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
  manifest <- runner_read_json(files$manifest)
  srt_out@tools[[tool_name]] <- list(
    method = "Cell2fate",
    result_type = "rna_velocity",
    cells = rownames(metadata),
    features = if (is.null(velocity)) {
      manifest$features %||% character()
    } else {
      colnames(velocity)
    },
    cell_metadata = metadata,
    velocity = velocity,
    manifest = manifest,
    paths = lapply(files, normalizePath, winslash = "/", mustWork = FALSE),
    parameters = config[setdiff(names(config), c("input_path", "result_dir"))],
    provenance = list(
      producer = "RunCell2fate",
      backend_id = "cell2fate",
      backend_repository = cell2fate_repository(),
      backend_commit = manifest$backend_commit %||% cell2fate_backend_commit()
    )
  )

  log_message(
    "{.pkg Cell2fate} analysis completed",
    message_type = "success",
    verbose = verbose
  )
  srt_out
}

cell2fate_repository <- function() {
  "BayraktarLab/cell2fate"
}

cell2fate_backend_commit <- function() {
  "c03d1ca0bb963f550001c6070d4986a61ec8456a"
}

cell2fate_requirement <- function() {
  paste0(
    "git+https://github.com/",
    cell2fate_repository(),
    ".git@",
    cell2fate_backend_commit()
  )
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
  jaxlib_ok <- check_python(
    "jaxlib==0.4.10",
    envname = target_env,
    pip = FALSE,
    verbose = FALSE
  )
  runtime_ok <- check_python(
    c(
      cell2fate_requirement(),
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
  if (!isTRUE(jaxlib_ok) || !isTRUE(runtime_ok)) {
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
