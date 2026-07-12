#' @title Run cell2location spatial deconvolution
#'
#' @description
#' Estimate spot-level absolute cell abundance and proportions with the official
#' Python `cell2location` backend. The Python model runs in an isolated
#' subprocess while inputs, models, posterior outputs, logs, and a reproducible
#' manifest are persisted under `result_dir`.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt Spatial `Seurat` object containing raw counts.
#' @param result_dir Directory used to persist inputs, models, posterior results,
#' tables, logs, and the run manifest.
#' @param reference Optional single-cell `Seurat` reference. Required when
#' `reference_signatures` is `NULL`.
#' @param reference_label Metadata column in `reference` containing cell types.
#' @param reference_signatures Optional gene-by-cell-type numeric matrix,
#' data.frame, or CSV file. When supplied, reference-model training is skipped.
#' @param assay,reference_assay Assays containing spatial and reference counts.
#' @param layer,reference_layer Raw-count layers.
#' @param features Optional features to consider before shared non-zero genes are
#' selected.
#' @param spatial_batch,reference_batch Optional metadata columns describing
#' spatial and reference batches.
#' @param reference_covariates Optional categorical reference metadata columns
#' used to model technical effects.
#' @param min_cells Minimum number of reference cells retained per cell type.
#' @param N_cells_per_location Expected average number of cells per spatial
#' location.
#' @param detection_alpha Prior controlling within-batch variation in RNA
#' detection sensitivity.
#' @param gene_filter_params Named arguments passed to
#' `cell2location.utils.filtering.filter_genes()`.
#' @param reference_train_params,spatial_train_params Named arguments passed to
#' the reference and spatial model `train()` methods.
#' @param reference_posterior_params Named sampling arguments used when exporting
#' reference signatures.
#' @param spatial_posterior_params Named sampling arguments used when exporting
#' the spatial q05 posterior.
#' @param envname Name of the shared SCOP Python environment.
#' @param resume Reuse completed stages only when their manifest matches the
#' current inputs and parameters.
#' @param overwrite Permit replacement of incompatible existing artifacts.
#' @param prefix Prefix used for abundance/proportion metadata columns.
#' @param tool_name Name of the `srt@tools` result entry.
#' @param store_results Whether to store detailed result matrices and paths in
#' `srt@tools`.
#'
#' @return A `Seurat` object with cell2location abundance, proportion, dominant
#' cell type, and maximum-proportion metadata.
#'
#' @section Official human lymph node result:
#' The figure below was generated from the official cell2location Human Lymph
#' Node tutorial data using 600 genuine Visium locations, 1,600 genuine
#' reference cells across 20 annotated cell types, and the complete two-stage
#' model implemented by this function.
#'
#' \if{html}{\figure{cell2location_human_lymph_node.png}{options: width=100\% alt="cell2location cell-type proportions in the official Human Lymph Node Visium dataset"}}
#' @concept spatial-producer
#' @export
#'
#' @examples
#' \dontrun{
#' # Official cell2location Human Lymph Node tutorial data:
#' # https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html
#' reference <- h5ad_to_srt("reference_subset.h5ad")
#' spatial <- h5ad_to_srt("spatial_subset.h5ad")
#' spatial <- RunCell2location(
#'   srt = spatial,
#'   result_dir = "human_lymph_node_cell2location",
#'   reference = reference,
#'   reference_label = "Subset",
#'   reference_batch = "Sample",
#'   spatial_batch = "sample",
#'   N_cells_per_location = 30,
#'   detection_alpha = 20
#' )
#' Cell2locationPlot(
#'   spatial,
#'   plot_type = "proportion",
#'   cell_types = c("B_naive", "T_CD4+_naive", "FDC"),
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#' }
RunCell2location <- function(
  srt,
  result_dir,
  reference = NULL,
  reference_label = "celltype",
  reference_signatures = NULL,
  assay = NULL,
  reference_assay = NULL,
  layer = "counts",
  reference_layer = "counts",
  features = NULL,
  spatial_batch = NULL,
  reference_batch = NULL,
  reference_covariates = NULL,
  min_cells = 5L,
  N_cells_per_location = 30,
  detection_alpha = 20,
  gene_filter_params = list(
    cell_count_cutoff = 5,
    cell_percentage_cutoff2 = 0.03,
    nonz_mean_cutoff = 1.12
  ),
  reference_train_params = list(
    max_epochs = 250L,
    batch_size = 2500L,
    train_size = 1,
    lr = 0.002,
    accelerator = "auto",
    device = "auto"
  ),
  spatial_train_params = list(
    max_epochs = 30000L,
    batch_size = NULL,
    train_size = 1,
    accelerator = "auto",
    device = "auto"
  ),
  reference_posterior_params = list(
    num_samples = 1000L,
    batch_size = 2500L
  ),
  spatial_posterior_params = list(
    batch_size = 2500L
  ),
  envname = NULL,
  resume = TRUE,
  overwrite = FALSE,
  prefix = "Cell2location",
  tool_name = "Cell2location",
  store_results = TRUE,
  verbose = TRUE
) {
  cell2location_assert_flag(resume, "resume")
  cell2location_assert_flag(overwrite, "overwrite")
  cell2location_assert_flag(store_results, "store_results")
  cell2location_assert_string(prefix, "prefix")
  cell2location_assert_string(tool_name, "tool_name")
  cell2location_validate_param_list(gene_filter_params, "gene_filter_params")
  cell2location_validate_param_list(reference_train_params, "reference_train_params")
  cell2location_validate_param_list(spatial_train_params, "spatial_train_params")
  cell2location_validate_param_list(reference_posterior_params, "reference_posterior_params")
  cell2location_validate_param_list(spatial_posterior_params, "spatial_posterior_params")
  min_cells <- cell2location_positive_integer(min_cells, "min_cells")
  cell2location_positive_number(N_cells_per_location, "N_cells_per_location")
  cell2location_positive_number(detection_alpha, "detection_alpha")

  prepared <- cell2location_prepare_inputs(
    srt = srt,
    reference = reference,
    reference_label = reference_label,
    reference_signatures = reference_signatures,
    assay = assay,
    reference_assay = reference_assay,
    layer = layer,
    reference_layer = reference_layer,
    features = features,
    spatial_batch = spatial_batch,
    reference_batch = reference_batch,
    reference_covariates = reference_covariates,
    min_cells = min_cells,
    verbose = verbose
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

  python <- cell2location_check_python(
    envname = envname,
    verbose = verbose
  )
  workdir <- tempfile("cell2location_inputs_")
  dir.create(workdir, recursive = TRUE)
  on.exit(unlink(workdir, recursive = TRUE, force = TRUE), add = TRUE)

  # The cell2location environment was prepared above. Keep both h5ad exports on
  # that exact interpreter instead of letting srt_to_h5ad prepare/switch to a
  # second scanpy-only environment; restore the existing option on every exit.
  old_skip_prepare <- getOption("scop_skip_python_prepare", FALSE)
  options(scop_skip_python_prepare = TRUE)
  on.exit(
    options(scop_skip_python_prepare = old_skip_prepare),
    add = TRUE
  )

  spatial_path <- file.path(workdir, "spatial.h5ad")
  srt_to_h5ad(
    prepared$spatial,
    path = spatial_path,
    features = prepared$features,
    assay_x = prepared$assay,
    layer_x = layer,
    assay_y = character(),
    reductions = character(),
    graphs = character(),
    neighbors = character(),
    overwrite = TRUE,
    verbose = FALSE
  )

  reference_path <- NULL
  signatures_path <- NULL
  if (is.null(prepared$signatures)) {
    reference_path <- file.path(workdir, "reference.h5ad")
    srt_to_h5ad(
      prepared$reference,
      path = reference_path,
      features = prepared$features,
      assay_x = prepared$reference_assay,
      layer_x = reference_layer,
      assay_y = character(),
      reductions = character(),
      graphs = character(),
      neighbors = character(),
      overwrite = TRUE,
      verbose = FALSE
    )
  } else {
    signatures_path <- file.path(workdir, "reference_signatures.csv")
    utils::write.csv(
      prepared$signatures,
      file = signatures_path,
      quote = TRUE,
      row.names = TRUE
    )
  }

  config <- list(
    result_dir = result_dir,
    spatial_path = normalizePath(spatial_path, winslash = "/", mustWork = TRUE),
    reference_path = if (is.null(reference_path)) NULL else normalizePath(reference_path, winslash = "/", mustWork = TRUE),
    signatures_path = if (is.null(signatures_path)) NULL else normalizePath(signatures_path, winslash = "/", mustWork = TRUE),
    reference_label = if (is.null(prepared$signatures)) reference_label else NULL,
    spatial_batch = spatial_batch,
    reference_batch = reference_batch,
    reference_covariates = reference_covariates,
    N_cells_per_location = N_cells_per_location,
    detection_alpha = detection_alpha,
    gene_filter_params = gene_filter_params,
    reference_train_params = reference_train_params,
    spatial_train_params = spatial_train_params,
    reference_posterior_params = reference_posterior_params,
    spatial_posterior_params = spatial_posterior_params,
    resume = resume,
    overwrite = overwrite
  )
  config_path <- file.path(workdir, "config.json")
  cell2location_write_json(config, config_path)
  runner <- cell2location_runner_path()
  stdout_path <- file.path(logs_dir, "cell2location_stdout.log")
  stderr_path <- file.path(logs_dir, "cell2location_stderr.log")
  cache_dir <- file.path(result_dir, ".cache")
  dir.create(file.path(cache_dir, "numba"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cache_dir, "matplotlib"), recursive = TRUE, showWarnings = FALSE)

  log_message(
    "Run {.pkg cell2location} with {.val {length(prepared$features)}} shared genes and {.val {ncol(prepared$spatial)}} spatial locations",
    verbose = verbose
  )
  status <- cell2location_run_system2(
    command = python,
    args = c(shQuote(runner), "--config", shQuote(config_path)),
    env = c(
      "PYTHONNOUSERSITE=1",
      paste0("NUMBA_CACHE_DIR=", file.path(cache_dir, "numba")),
      paste0("MPLCONFIGDIR=", file.path(cache_dir, "matplotlib"))
    ),
    stdout = stdout_path,
    stderr = stderr_path
  )
  if (!identical(status, 0L)) {
    cell2location_runner_error(status, stdout_path, stderr_path)
  }

  files <- cell2location_result_files(result_dir)
  missing_files <- unlist(files[c("abundance", "proportions", "signatures", "manifest")])
  missing_files <- missing_files[!file.exists(missing_files)]
  if (length(missing_files) > 0L) {
    log_message(
      "{.pkg cell2location} did not produce required file{?s}: {.file {missing_files}}",
      message_type = "error"
    )
  }
  abundance <- cell2location_read_numeric_csv(files$abundance, "abundance")
  proportions <- cell2location_read_numeric_csv(files$proportions, "proportions")
  signatures <- cell2location_read_numeric_csv(files$signatures, "reference signatures")
  abundance <- cell2location_align_result(abundance, colnames(prepared$spatial), "abundance")
  proportions <- cell2location_align_result(proportions, colnames(prepared$spatial), "proportions")
  if (!identical(colnames(abundance), colnames(proportions))) {
    log_message(
      "cell2location abundance and proportion cell types do not match",
      message_type = "error"
    )
  }

  full_abundance <- matrix(
    NA_real_,
    nrow = ncol(srt),
    ncol = ncol(abundance),
    dimnames = list(colnames(srt), colnames(abundance))
  )
  full_abundance[rownames(abundance), ] <- as.matrix(abundance)
  abundance_meta <- as.data.frame(full_abundance, check.names = FALSE)
  colnames(abundance_meta) <- paste0(
    prefix,
    "_abundance_",
    make.unique(make.names(colnames(full_abundance)), sep = "_")
  )
  srt_out <- Seurat::AddMetaData(srt, metadata = abundance_meta)
  weight_summary <- scop_spatial_finalize_weights(
    proportions,
    all_spots = colnames(srt_out)
  )
  srt_out <- scop_spatial_add_deconv_metadata(
    srt_out,
    weights = proportions,
    prefix = prefix,
    metadata = weight_summary
  )

  if (isTRUE(store_results)) {
    manifest <- cell2location_read_json(files$manifest)
    result_parameters <- list(
      assay = prepared$assay,
      reference_assay = prepared$reference_assay,
      layer = layer,
      reference_layer = reference_layer,
      reference_label = reference_label,
      spatial_batch = spatial_batch,
      reference_batch = reference_batch,
      reference_covariates = reference_covariates,
      min_cells = min_cells,
      N_cells_per_location = N_cells_per_location,
      detection_alpha = detection_alpha,
      gene_filter_params = gene_filter_params,
      reference_train_params = reference_train_params,
      spatial_train_params = spatial_train_params,
      reference_posterior_params = reference_posterior_params,
      spatial_posterior_params = spatial_posterior_params,
      envname = get_envname(envname),
      resume = resume,
      overwrite = overwrite,
      prefix = prefix
    )
    backend_versions <- unlist(manifest$versions %||% character(), use.names = TRUE)
    backend_versions <- as.character(backend_versions[!is.na(backend_versions)])
    srt_out@tools[[tool_name]] <- spatial_result_build(
      bundle = list(
        abundance = abundance,
        proportions = proportions,
        reference_signatures = signatures,
        input_summary = prepared$summary,
        manifest = manifest,
        files = files,
        cells = rownames(proportions)
      ),
      method = "Cell2location",
      result_type = "deconvolution",
      source = list(
        image = character(),
        coordinate_space = "none",
        transform = NULL,
        assay = prepared$assay,
        layer = layer
      ),
      provenance = list(
        producer = "RunCell2location",
        backend_id = "cell2location",
        backend_versions = backend_versions
      ),
      parameters = result_parameters,
      summary = scop_spatial_weight_summary(proportions)
    )
  }

  log_message(
    "{.pkg cell2location} q05 abundance and proportions stored with prefix {.val {prefix}}",
    message_type = "success",
    verbose = verbose
  )
  srt_out
}

#' @title Plot cell2location spatial results
#'
#' @md
#' @inheritParams SpatialSpotPlot
#' @param srt A `Seurat` object returned by [RunCell2location()].
#' @param plot_type Result to draw: normalized proportion, q05 absolute
#' abundance, dominant cell type, or a spot-level proportion pie.
#' @param cell_types Optional cell types to display for abundance, proportion,
#' or pie plots.
#' @param prefix Metadata prefix used by [RunCell2location()].
#' @param tool_name Name of the `srt@tools` result entry.
#' @param ... Additional arguments passed to [SpatialSpotPlot()].
#'
#' @return A `ggplot`, `patchwork`, or list of plots.
#' @export
#'
#' @examples
#' \dontrun{
#' # Result from the official Human Lymph Node example in RunCell2location().
#' spatial <- readRDS("human_lymph_node_cell2location/official_human_lymph_node.rds")
#' selected <- names(sort(
#'   colMeans(spatial@tools$Cell2location$proportions),
#'   decreasing = TRUE
#' ))[1:6]
#' Cell2locationPlot(
#'   spatial,
#'   plot_type = "proportion",
#'   cell_types = selected,
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y"),
#'   ncol = 3
#' )
#' Cell2locationPlot(
#'   spatial,
#'   plot_type = "dominant",
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#' }
Cell2locationPlot <- function(
  srt,
  plot_type = c("proportion", "abundance", "dominant", "pie"),
  cell_types = NULL,
  prefix = "Cell2location",
  tool_name = "Cell2location",
  image = NULL,
  overlay_image = TRUE,
  coord.cols = c("col", "row"),
  ...
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  plot_type <- match.arg(plot_type)
  coords <- spatial_dim_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    overlay_image = overlay_image
  )$data
  if (!all(c("x", "y") %in% colnames(coords)) || any(!is.finite(as.matrix(coords[, c("x", "y"), drop = FALSE])))) {
    log_message("Spatial coordinates must exist and contain finite numeric values", message_type = "error")
  }

  tool <- srt@tools[[tool_name]]
  values <- NULL
  if (plot_type %in% c("proportion", "pie")) {
    values <- tool$proportions %||% cell2location_metadata_matrix(srt, paste0(prefix, "_prop_"))
  } else if (identical(plot_type, "abundance")) {
    values <- tool$abundance %||% cell2location_metadata_matrix(srt, paste0(prefix, "_abundance_"))
  }
  if (!is.null(values)) {
    values <- as.data.frame(values, check.names = FALSE)
    if (!is.null(cell_types)) {
      missing <- setdiff(cell_types, colnames(values))
      if (length(missing) > 0L) {
        log_message("Unknown {.arg cell_types}: {.val {missing}}", message_type = "error")
      }
      values <- values[, cell_types, drop = FALSE]
    }
  }

  defaults <- list(
    srt = srt,
    image = image,
    overlay_image = overlay_image,
    coord.cols = coord.cols
  )
  if (identical(plot_type, "dominant")) {
    defaults$group.by <- paste0(prefix, "_dominant_type")
  } else {
    defaults$values <- values
    defaults$plot_type <- if (identical(plot_type, "pie")) "pie" else "point"
  }
  do.call(SpatialSpotPlot, standard_scop_merge_args(defaults, list(...)))
}

cell2location_run_system2 <- function(command, args, env, stdout, stderr) {
  env_names <- sub("=.*$", "", env)
  env_values <- sub("^[^=]*=", "", env)
  old_values <- Sys.getenv(env_names, unset = NA_character_)
  names(old_values) <- env_names
  do.call(Sys.setenv, stats::setNames(as.list(env_values), env_names))
  on.exit({
    restore <- old_values[!is.na(old_values)]
    remove <- names(old_values)[is.na(old_values)]
    if (length(restore) > 0L) {
      do.call(Sys.setenv, as.list(restore))
    }
    if (length(remove) > 0L) {
      Sys.unsetenv(remove)
    }
  }, add = TRUE)
  system2(
    command = command,
    args = args,
    stdout = stdout,
    stderr = stderr
  )
}

cell2location_prepare_inputs <- function(
  srt,
  reference,
  reference_label,
  reference_signatures,
  assay,
  reference_assay,
  layer,
  reference_layer,
  features,
  spatial_batch,
  reference_batch,
  reference_covariates,
  min_cells,
  verbose
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  if (is.null(reference) && is.null(reference_signatures)) {
    log_message(
      "Provide {.arg reference} or {.arg reference_signatures}",
      message_type = "error"
    )
  }
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  cell2location_validate_assay(srt, assay, "assay")
  cell2location_validate_columns(srt, spatial_batch, "spatial_batch")
  cell2location_validate_unique_names(srt, "srt")

  spatial_counts <- cell2location_counts(srt, assay, layer, "Spatial")
  keep_spots <- Matrix::colSums(spatial_counts) > 0
  if (!any(keep_spots)) {
    log_message("No non-zero spatial locations remain", message_type = "error")
  }
  dropped_spots <- colnames(spatial_counts)[!keep_spots]
  spatial_counts <- spatial_counts[, keep_spots, drop = FALSE]
  srt_use <- srt[, colnames(spatial_counts)]

  signatures <- NULL
  reference_use <- NULL
  dropped_types <- data.frame(cell_type = character(), n_cells = integer())
  if (!is.null(reference_signatures)) {
    signatures <- cell2location_read_signatures(reference_signatures)
    common <- intersect(rownames(spatial_counts), rownames(signatures))
    if (!is.null(features)) {
      common <- intersect(features, common)
    }
    keep_features <- Matrix::rowSums(spatial_counts[common, , drop = FALSE]) > 0 &
      rowSums(signatures[common, , drop = FALSE]) > 0
    common <- common[keep_features]
    signatures <- signatures[common, , drop = FALSE]
    reference_assay <- NULL
  } else {
    if (!inherits(reference, "Seurat")) {
      log_message("{.arg reference} must be a {.cls Seurat} object", message_type = "error")
    }
    reference_assay <- reference_assay %||% SeuratObject::DefaultAssay(reference)
    cell2location_validate_assay(reference, reference_assay, "reference_assay")
    cell2location_validate_unique_names(reference, "reference")
    cell2location_validate_columns(reference, reference_label, "reference_label")
    cell2location_validate_columns(reference, reference_batch, "reference_batch")
    cell2location_validate_columns(reference, reference_covariates, "reference_covariates")
    labels <- as.character(reference[[reference_label, drop = TRUE]])
    names(labels) <- colnames(reference)
    valid <- !is.na(labels) & nzchar(labels)
    label_counts <- table(labels[valid])
    keep_types <- names(label_counts)[label_counts >= min_cells]
    drop_types <- names(label_counts)[label_counts < min_cells]
    dropped_types <- data.frame(
      cell_type = drop_types,
      n_cells = as.integer(label_counts[drop_types]),
      stringsAsFactors = FALSE
    )
    keep_cells <- valid & labels %in% keep_types
    if (!any(keep_cells)) {
      log_message("No reference cell types remain after {.arg min_cells} filtering", message_type = "error")
    }
    reference_use <- reference[, names(labels)[keep_cells]]
    reference_counts <- cell2location_counts(
      reference_use,
      reference_assay,
      reference_layer,
      "Reference"
    )
    keep_ref_cells <- Matrix::colSums(reference_counts) > 0
    reference_counts <- reference_counts[, keep_ref_cells, drop = FALSE]
    reference_use <- reference_use[, colnames(reference_counts)]
    labels_after <- as.character(reference_use[[reference_label, drop = TRUE]])
    counts_after <- table(labels_after)
    keep_types_after <- names(counts_after)[counts_after >= min_cells]
    reference_use <- reference_use[, labels_after %in% keep_types_after]
    if (ncol(reference_use) == 0L) {
      log_message("No reference cells remain after zero-count filtering", message_type = "error")
    }
    reference_counts <- reference_counts[, colnames(reference_use), drop = FALSE]
    common <- intersect(rownames(spatial_counts), rownames(reference_counts))
    if (!is.null(features)) {
      common <- intersect(features, common)
    }
    keep_features <- Matrix::rowSums(spatial_counts[common, , drop = FALSE]) > 0 &
      Matrix::rowSums(reference_counts[common, , drop = FALSE]) > 0
    common <- common[keep_features]
  }
  if (length(common) == 0L) {
    log_message("No shared non-zero features are available for {.fn RunCell2location}", message_type = "error")
  }

  log_message(
    "Use {.val {length(common)}} shared non-zero features for {.pkg cell2location}",
    verbose = verbose
  )
  list(
    spatial = srt_use[common, ],
    reference = if (is.null(reference_use)) NULL else reference_use[common, ],
    signatures = signatures,
    features = common,
    assay = assay,
    reference_assay = reference_assay,
    summary = list(
      n_spots = ncol(srt_use),
      n_reference_cells = if (is.null(reference_use)) NA_integer_ else ncol(reference_use),
      n_features = length(common),
      dropped_spots = dropped_spots,
      dropped_cell_types = dropped_types
    )
  )
}

cell2location_counts <- function(srt, assay, layer, label) {
  counts <- GetAssayData5(srt, assay = assay, layer = layer)
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    log_message("{.val {label}} counts must have feature and observation names", message_type = "error")
  }
  values <- if (inherits(counts, "sparseMatrix")) counts@x else as.numeric(counts)
  if (any(!is.finite(values)) || any(values < 0)) {
    log_message("{.val {label}} counts must be finite and non-negative", message_type = "error")
  }
  if (any(abs(values - round(values)) > sqrt(.Machine$double.eps))) {
    log_message(
      "{.val {label}} layer {.val {layer}} must contain raw integer counts; normalized values are not accepted",
      message_type = "error"
    )
  }
  counts
}

cell2location_read_signatures <- function(x) {
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) {
      log_message("{.arg reference_signatures} file does not exist: {.file {x}}", message_type = "error")
    }
    x <- utils::read.csv(x, row.names = 1, check.names = FALSE)
  }
  if (!is.matrix(x) && !is.data.frame(x)) {
    log_message("{.arg reference_signatures} must be a matrix, data.frame, or CSV file", message_type = "error")
  }
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (is.null(rownames(x)) || is.null(colnames(x)) || anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
    log_message("{.arg reference_signatures} must have unique gene and cell-type names", message_type = "error")
  }
  if (any(!is.finite(x)) || any(x < 0)) {
    log_message("{.arg reference_signatures} must be finite and non-negative", message_type = "error")
  }
  x
}

cell2location_check_python <- function(envname = NULL, verbose = TRUE) {
  PrepareEnv(envname = envname, modules = "cell2location", verbose = verbose)
  ok <- check_python(
    c("cell2location==0.1.5", "scvi-tools==1.3.3", "scanpy", "anndata", "numpy", "pandas", "scipy", "torch"),
    envname = envname,
    verbose = FALSE
  )
  if (!isTRUE(ok)) {
    log_message("The shared SCOP Python environment is missing cell2location dependencies", message_type = "error")
  }
  cache <- getOption("scop_env_cache", default = NULL)
  python <- cache[["python"]] %||% tryCatch(
    conda_python(envname = get_envname(envname), conda = resolve_conda("auto")),
    error = function(...) NULL
  )
  if (is.null(python) || !file.exists(python)) {
    log_message("Unable to resolve the cell2location Python executable", message_type = "error")
  }
  normalizePath(python, winslash = "/", mustWork = TRUE)
}

cell2location_runner_path <- function() {
  candidates <- c(
    system.file("python", "cell2location_runner.py", package = "scop", mustWork = FALSE),
    file.path("inst", "python", "cell2location_runner.py")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (length(candidates) == 0L) {
    log_message("Bundled cell2location Python runner was not found", message_type = "error")
  }
  normalizePath(candidates[1L], winslash = "/", mustWork = TRUE)
}

cell2location_write_json <- function(x, path) {
  check_r("jsonlite", verbose = FALSE)
  to_json <- get_namespace_fun("jsonlite", "toJSON")
  writeLines(
    as.character(to_json(x, auto_unbox = TRUE, null = "null", digits = NA, pretty = TRUE)),
    con = path,
    useBytes = TRUE
  )
}

cell2location_read_json <- function(path) {
  check_r("jsonlite", verbose = FALSE)
  get_namespace_fun("jsonlite", "fromJSON")(path, simplifyVector = FALSE)
}

cell2location_result_files <- function(result_dir) {
  list(
    abundance = file.path(result_dir, "tables", "abundance_q05.csv"),
    proportions = file.path(result_dir, "tables", "proportions.csv"),
    signatures = file.path(result_dir, "reference", "signatures.csv"),
    manifest = file.path(result_dir, "manifest.json"),
    spatial_posterior = file.path(result_dir, "spatial", "spatial_posterior.h5ad"),
    reference_posterior = file.path(result_dir, "reference", "reference_posterior.h5ad"),
    reference_model = file.path(result_dir, "reference", "model"),
    spatial_model = file.path(result_dir, "spatial", "model"),
    stdout = file.path(result_dir, "logs", "cell2location_stdout.log"),
    stderr = file.path(result_dir, "logs", "cell2location_stderr.log")
  )
}

cell2location_read_numeric_csv <- function(path, label) {
  x <- utils::read.csv(path, row.names = 1, check.names = FALSE)
  if (nrow(x) == 0L || ncol(x) == 0L || any(!vapply(x, is.numeric, logical(1)))) {
    log_message("cell2location {.val {label}} output must be a non-empty numeric table", message_type = "error")
  }
  if (any(!is.finite(as.matrix(x))) || is.null(rownames(x)) || is.null(colnames(x))) {
    log_message("cell2location {.val {label}} output contains invalid values or names", message_type = "error")
  }
  x
}

cell2location_align_result <- function(x, spot_ids, label) {
  if (anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
    log_message("cell2location {.val {label}} output names must be unique", message_type = "error")
  }
  missing <- setdiff(spot_ids, rownames(x))
  extra <- setdiff(rownames(x), spot_ids)
  if (length(missing) > 0L || length(extra) > 0L) {
    log_message(
      "cell2location {.val {label}} spot names do not match the prepared spatial input",
      message_type = "error"
    )
  }
  x[spot_ids, , drop = FALSE]
}

cell2location_runner_error <- function(status, stdout_path, stderr_path) {
  read_tail <- function(path) {
    if (!file.exists(path)) return(character())
    utils::tail(readLines(path, warn = FALSE), 30L)
  }
  lines <- c(read_tail(stderr_path), read_tail(stdout_path))
  if (length(lines) == 0L) lines <- "<no Python output captured>"
  log_message(
    "{.pkg cell2location} Python runner failed with status {.val {status}}:\n{.code {paste(lines, collapse = '\n')}}",
    message_type = "error"
  )
}

cell2location_metadata_matrix <- function(srt, prefix) {
  cols <- colnames(srt@meta.data)[startsWith(colnames(srt@meta.data), prefix)]
  if (length(cols) == 0L) {
    log_message("No cell2location result columns with prefix {.val {prefix}} were found", message_type = "error")
  }
  out <- srt@meta.data[, cols, drop = FALSE]
  colnames(out) <- sub(paste0("^", prefix), "", cols)
  out
}

cell2location_validate_assay <- function(srt, assay, arg) {
  if (length(assay) != 1L || is.na(assay) || !assay %in% SeuratObject::Assays(srt)) {
    log_message("{.arg {arg}} must identify one assay in the Seurat object", message_type = "error")
  }
}

cell2location_validate_columns <- function(srt, columns, arg) {
  if (is.null(columns) || length(columns) == 0L) return(invisible(TRUE))
  if (!is.character(columns) || any(!columns %in% colnames(srt@meta.data))) {
    log_message("{.arg {arg}} contains missing Seurat metadata columns", message_type = "error")
  }
  invisible(TRUE)
}

cell2location_validate_unique_names <- function(srt, arg) {
  if (anyDuplicated(colnames(srt)) || anyDuplicated(rownames(srt))) {
    log_message("{.arg {arg}} must have unique observation and feature names", message_type = "error")
  }
}

cell2location_validate_param_list <- function(x, arg) {
  if (
    !is.list(x) ||
      (
        length(x) > 0L &&
          (is.null(names(x)) || any(is.na(names(x)) | !nzchar(names(x))))
      )
  ) {
    log_message("{.arg {arg}} must be a named list", message_type = "error")
  }
}

cell2location_assert_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    log_message("{.arg {arg}} must be TRUE or FALSE", message_type = "error")
  }
}

cell2location_assert_string <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    log_message("{.arg {arg}} must be one non-empty string", message_type = "error")
  }
}

cell2location_positive_integer <- function(x, arg) {
  cell2location_positive_number(x, arg)
  if (x != as.integer(x)) {
    log_message("{.arg {arg}} must be an integer", message_type = "error")
  }
  as.integer(x)
}

cell2location_positive_number <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
    log_message("{.arg {arg}} must be one positive finite number", message_type = "error")
  }
  invisible(x)
}
