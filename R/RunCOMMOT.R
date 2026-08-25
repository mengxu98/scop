# Official COMMOT producer and stored-result plotting --------------------------

.commot_backend_commit <- "d117445bc07eaa19109c7609e97d9e35b26e99ca"

commot_validate_param_list <- function(x, name) {
  if (!is.list(x) || (length(x) > 0L && (is.null(names(x)) || any(!nzchar(names(x)))))) {
    log_message("{.arg {name}} must be a named list", message_type = "error")
  }
  invisible(TRUE)
}


commot_validate_matrix <- function(x) {
  if ((!is.matrix(x) && !methods::is(x, "Matrix")) ||
      is.null(rownames(x)) || is.null(colnames(x)) ||
      anyNA(rownames(x)) || anyNA(colnames(x)) ||
      any(!nzchar(rownames(x))) || any(!nzchar(colnames(x))) ||
      anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
    log_message("COMMOT expression must be a named matrix with unique feature and cell IDs", message_type = "error")
  }
  values <- matrix_values(x)
  if (any(!is.finite(values)) || any(values < 0)) {
    log_message("COMMOT expression must contain finite non-negative values", message_type = "error")
  }
  invisible(TRUE)
}

commot_input <- function(srt, group.by, assay, layer, image, coord.cols) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  if (!is.character(group.by) || length(group.by) != 1L || !group.by %in% colnames(srt[[]])) {
    log_message("{.arg group.by} must be one metadata column", message_type = "error")
  }
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!assay %in% SeuratObject::Assays(srt)) {
    log_message("Assay {.val {assay}} is absent from {.arg srt}", message_type = "error")
  }
  coords <- SpatialCoordinates(
    object = srt, image = image, coord.cols = coord.cols,
    space = "raw", image_policy = "strict"
  )
  cells <- as.character(coords$data$cell_id)
  if (length(cells) == 0L || anyDuplicated(cells) || any(!cells %in% colnames(srt))) {
    log_message("Spatial coordinates do not identify a valid unique subset of {.arg srt}", message_type = "error")
  }
  coords$data <- coords$data[match(cells, coords$data$cell_id), , drop = FALSE]
  if (any(!is.finite(coords$data$x)) || any(!is.finite(coords$data$y))) {
    log_message("COMMOT requires finite raw x and y coordinates", message_type = "error")
  }
  expression <- GetAssayData5(srt, assay = assay, layer = layer)[, cells, drop = FALSE]
  commot_validate_matrix(expression)
  groups <- as.character(srt[[]][cells, group.by])
  if (anyNA(groups) || any(!nzchar(groups))) {
    log_message("{.arg group.by} contains missing or empty labels", message_type = "error")
  }
  list(
    expression = expression, coordinates = coords$data, cells = cells,
    groups = groups, source = coords$source, assay = assay
  )
}

commot_runtime <- function(envname, verbose) {
  PrepareEnv(envname = envname, modules = "commot", verbose = verbose)
  requirements <- env_requirements(modules = "commot", verbose = FALSE)$packages
  available <- check_python(
    requirements,
    envname = envname, verbose = FALSE
  )
  if (!isTRUE(available)) {
    log_message("The isolated {.file {envname}} environment is missing official COMMOT dependencies", message_type = "error")
  }
  conda <- resolve_conda("auto")
  python <- resolve_python_executable(envname = envname, conda = conda)
  normalizePath(python, winslash = "/", mustWork = TRUE)
}


commot_write_inputs <- function(input, workdir) {
  matrix_path <- file.path(workdir, "expression.mtx")
  features_path <- file.path(workdir, "features.tsv")
  metadata_path <- file.path(workdir, "metadata.tsv")
  coordinates_path <- file.path(workdir, "coordinates.tsv")
  Matrix::writeMM(input$expression, matrix_path)
  write_tsv(data.frame(feature_id = rownames(input$expression)), features_path)
  write_tsv(data.frame(cell_id = input$cells, group = input$groups), metadata_path)
  write_tsv(data.frame(
    cell_id = input$cells,
    x = as.numeric(input$coordinates$x),
    y = as.numeric(input$coordinates$y)
  ), coordinates_path)
  list(
    matrix_path = normalizePath(matrix_path, winslash = "/", mustWork = TRUE),
    features_path = normalizePath(features_path, winslash = "/", mustWork = TRUE),
    metadata_path = normalizePath(metadata_path, winslash = "/", mustWork = TRUE),
    coordinates_path = normalizePath(coordinates_path, winslash = "/", mustWork = TRUE)
  )
}

commot_read_table <- function(path, required, numeric, label, allow_empty = FALSE) {
  if (!file.exists(path)) {
    log_message("COMMOT did not produce {.val {label}} output", message_type = "error")
  }
  value <- tryCatch(
    utils::read.csv(
      path, stringsAsFactors = FALSE, check.names = FALSE,
      na.strings = c("", "NA", "NaN", "nan")
    ),
    error = function(e) {
      log_message("Unable to read COMMOT {.val {label}} output: {conditionMessage(e)}", message_type = "error")
    }
  )
  if (!all(required %in% colnames(value)) || (!isTRUE(allow_empty) && nrow(value) == 0L)) {
    log_message("COMMOT returned invalid {.val {label}} columns or no rows", message_type = "error")
  }
  for (column in intersect(numeric, colnames(value))) {
    value[[column]] <- suppressWarnings(as.numeric(value[[column]]))
    finite_or_missing <- is.finite(value[[column]]) | is.na(value[[column]])
    if (!all(finite_or_missing)) {
      log_message("COMMOT {.val {label}} has invalid numeric values in {.field {column}}", message_type = "error")
    }
  }
  value
}

commot_validate_manifest <- function(manifest, input, communication, cluster_table, direction_table) {
  cell_ids <- as.character(unlist(manifest$cell_ids %||% character(), use.names = FALSE))
  feature_ids <- as.character(unlist(manifest$feature_ids %||% character(), use.names = FALSE))
  valid <- identical(as.integer(manifest$n_obs), ncol(input$expression)) &&
    identical(as.integer(manifest$n_vars), nrow(input$expression)) &&
    identical(cell_ids, input$cells) &&
    identical(feature_ids, rownames(input$expression)) &&
    identical(as.integer(manifest$n_communication_rows), nrow(communication)) &&
    identical(as.integer(manifest$n_cluster_rows), nrow(cluster_table)) &&
    identical(as.integer(manifest$n_direction_rows), nrow(direction_table)) &&
    identical(as.character(manifest$coordinate_space), "raw")
  if (!isTRUE(valid)) {
    log_message("COMMOT manifest does not match the submitted matrix, IDs, coordinates, or outputs", message_type = "error")
  }
  invisible(TRUE)
}

commot_artifact_path <- function(output.dir, result.name) {
  output.dir <- output.dir %||% file.path(getwd(), "commot_results")
  output.dir <- normalizePath(path.expand(output.dir), winslash = "/", mustWork = FALSE)
  filename <- paste0("commot_", gsub("[^A-Za-z0-9_.-]", "_", result.name), ".h5ad")
  file.path(output.dir, filename)
}

commot_copy_h5ad <- function(source, destination, overwrite) {
  if (!file.exists(source)) {
    log_message("COMMOT requested H5AD output is absent", message_type = "error")
  }
  if (file.exists(destination) && !isTRUE(overwrite)) {
    log_message("COMMOT H5AD already exists at {.file {destination}}; set {.arg overwrite = TRUE}", message_type = "error")
  }
  directory <- dirname(destination)
  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE)) {
    log_message("Unable to create COMMOT artifact directory {.file {directory}}", message_type = "error")
  }
  temporary <- tempfile(pattern = paste0(".", basename(destination), "."), tmpdir = directory)
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  if (!file.copy(source, temporary, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)) {
    log_message("Unable to stage COMMOT H5AD artifact", message_type = "error")
  }
  backup <- NULL
  published <- FALSE
  if (file.exists(destination) && isTRUE(overwrite)) {
    backup <- tempfile(pattern = paste0(".", basename(destination), ".backup."), tmpdir = directory)
    if (!file.rename(destination, backup)) {
      log_message("Unable to stage the existing COMMOT H5AD for replacement", message_type = "error")
    }
  }
  on.exit({
    if (!published && !is.null(backup) && file.exists(backup) && !file.exists(destination)) {
      file.rename(backup, destination)
    }
  }, add = TRUE)
  if (!file.rename(temporary, destination)) {
    log_message("Unable to atomically publish COMMOT H5AD artifact", message_type = "error")
  }
  published <- TRUE
  if (!is.null(backup) && file.exists(backup)) unlink(backup, force = TRUE)
  normalizePath(destination, winslash = "/", mustWork = TRUE)
}

commot_execute <- function(
  input, species, database, distance.threshold, cluster, direction,
  communication.args, cluster.args, direction.args, result.name,
  envname, output.dir, store.h5ad, overwrite, verbose
) {
  python <- commot_runtime(envname, verbose)
  workdir <- tempfile("commot_run_")
  dir.create(workdir, recursive = TRUE)
  on.exit(unlink(workdir, recursive = TRUE, force = TRUE), add = TRUE)
  output <- file.path(workdir, "output")
  dir.create(output)
  paths <- commot_write_inputs(input, workdir)
  config <- c(paths, list(
    output_dir = normalizePath(output, winslash = "/", mustWork = TRUE),
    species = species, database = database,
    distance_threshold = distance.threshold,
    cluster = cluster, direction = direction,
    communication_args = communication.args,
    cluster_args = cluster.args, direction_args = direction.args,
    store_h5ad = store.h5ad
  ))
  config_path <- file.path(workdir, "config.json")
  runner_write_json(config, config_path)
  stdout_path <- file.path(workdir, "stdout.log")
  stderr_path <- file.path(workdir, "stderr.log")
  runner <- runner_script_path("commot_runner.py", "COMMOT")
  status <- runner_system2(
    command = python,
    args = c(shQuote(runner), "--config", shQuote(config_path)),
    env = c(PYTHONNOUSERSITE = "1"),
    stdout = stdout_path,
    stderr = stderr_path
  )
  if (!identical(status, 0L)) {
    runner_error(status, stdout_path, stderr_path, backend = "COMMOT")
  }
  communication <- commot_read_table(
    file.path(output, "communication.csv"),
    required = c("sender", "receiver", "ligand", "receptor", "interaction_name", "pathway_name", "score", "pvalue", "method", "score_type"),
    numeric = c("score", "pvalue"), label = "communication"
  )
  cluster_table <- commot_read_table(
    file.path(output, "cluster.csv"),
    required = c("sender", "receiver", "score", "pvalue", "key"),
    numeric = c("score", "pvalue"), label = "cluster", allow_empty = !cluster
  )
  direction_table <- commot_read_table(
    file.path(output, "direction.csv"),
    required = c("cell_id", "group", "key", "perspective", "x", "y", "dx", "dy", "magnitude"),
    numeric = c("x", "y", "dx", "dy", "magnitude"), label = "direction", allow_empty = !direction
  )
  manifest <- runner_read_json(file.path(output, "manifest.json"))
  commot_validate_manifest(manifest, input, communication, cluster_table, direction_table)
  if (!all(communication$method == "COMMOT") ||
      !all(communication$score_type == "transport_mass") ||
      any(!is.na(communication$pvalue))) {
    log_message("COMMOT transport output has invalid method, score type, or p-value semantics", message_type = "error")
  }
  artifact <- NULL
  if (isTRUE(store.h5ad)) {
    destination <- commot_artifact_path(output.dir, result.name)
    destination <- commot_copy_h5ad(file.path(output, "commot_result.h5ad"), destination, overwrite)
    artifact <- list(
      path = destination, exists = file.exists(destination),
      size = unname(file.info(destination)$size), manifest = manifest
    )
  }
  list(
    communication = communication, cluster_table = cluster_table,
    direction_table = direction_table, manifest = manifest, h5ad = artifact
  )
}


#' @title Run official COMMOT spatial communication analysis
#'
#' @description
#' Run official Python COMMOT in an isolated subprocess. Raw expression,
#' metadata, and raw coordinates are exchanged through files; AnnData is not
#' placed inside the Seurat object.
#'
#' @param srt A `Seurat` spatial object.
#' @param group.by Metadata column used to aggregate communication by group.
#' @param species COMMOT ligand-receptor database species.
#' @param database Official COMMOT ligand-receptor database name.
#' @param assay,layer Assay and raw-count layer.
#' @param image Explicit spatial image. Required for multi-image objects.
#' @param coord.cols Metadata coordinate columns used when no image is present.
#' @param distance.threshold Maximum signaling distance in raw-coordinate units.
#' @param cluster Whether to run official cluster communication permutations.
#' @param direction Whether to run official communication direction analysis.
#' @param communication.args Named arguments for official spatial communication.
#'   Runner-only fields are `normalize`, `target_sum`, `signaling_type`, and
#'   nested `filter_args`.
#' @param cluster.args Named arguments for official cluster communication.
#' @param direction.args Named arguments for official communication direction.
#' @param result.name Stored result name.
#' @param envname Isolated Python environment name.
#' @param output.dir External directory for a requested H5AD artifact.
#' @param store.h5ad Whether to retain the complete AnnData as an external H5AD.
#' @param overwrite Whether to replace an existing named result or H5AD artifact.
#' @param backend Backend used only for SCOP CCC result aggregation.
#' @param verbose Whether to print progress messages.
#'
#' @return The input `Seurat` object with `srt@tools$COMMOT` and unified CCC
#'   results updated.
#' @references <https://github.com/zcang/COMMOT>
#' @export
RunCOMMOT <- function(
  srt,
  group.by,
  species = c("human", "mouse"),
  database = "CellChat",
  assay = NULL,
  layer = "counts",
  image = NULL,
  coord.cols = c("col", "row"),
  distance.threshold = NULL,
  cluster = FALSE,
  direction = FALSE,
  communication.args = list(),
  cluster.args = list(),
  direction.args = list(),
  result.name = "default",
  envname = "commot_env",
  output.dir = NULL,
  store.h5ad = FALSE,
  overwrite = FALSE,
  backend = c("cpp", "r"),
  verbose = TRUE
) {
  species <- match.arg(species)
  backend <- match.arg(backend)
  validate_scalar_flag(cluster, "cluster")
  validate_scalar_flag(direction, "direction")
  validate_scalar_flag(store.h5ad, "store.h5ad")
  validate_scalar_flag(overwrite, "overwrite")
  commot_validate_param_list(communication.args, "communication.args")
  commot_validate_param_list(cluster.args, "cluster.args")
  commot_validate_param_list(direction.args, "direction.args")
  if (!is.character(database) || length(database) != 1L || is.na(database) || !nzchar(database)) {
    log_message("{.arg database} must be one non-empty string", message_type = "error")
  }
  if (!is.character(result.name) || length(result.name) != 1L || is.na(result.name) || !nzchar(result.name)) {
    log_message("{.arg result.name} must be one non-empty string", message_type = "error")
  }
  if (!is.null(distance.threshold) &&
      (!is.numeric(distance.threshold) || length(distance.threshold) != 1L ||
       is.na(distance.threshold) || !is.finite(distance.threshold) || distance.threshold <= 0)) {
    log_message("{.arg distance.threshold} must be NULL or one positive finite raw-coordinate distance", message_type = "error")
  }
  input <- commot_input(srt, group.by, assay, layer, image, coord.cols)
  existing <- srt@tools[["COMMOT"]]
  if (!is.null(existing$results[[result.name]]) && !isTRUE(overwrite)) {
    log_message("COMMOT result {.val {result.name}} already exists; set {.arg overwrite = TRUE}", message_type = "error")
  }
  executed <- commot_execute(
    input = input, species = species, database = database,
    distance.threshold = distance.threshold, cluster = cluster,
    direction = direction, communication.args = communication.args,
    cluster.args = cluster.args, direction.args = direction.args,
    result.name = result.name, envname = envname, output.dir = output.dir,
    store.h5ad = store.h5ad, overwrite = overwrite, verbose = verbose
  )
  long_table <- standardize_long_df(executed$communication)
  result <- list(
    lr_table = executed$communication,
    pathway_table = data.frame(),
    tf_table = data.frame(),
    deconvolution = NULL,
    cluster_table = executed$cluster_table,
    direction_table = executed$direction_table,
    long_table = long_table,
    h5ad = executed$h5ad,
    native_object = NULL,
    manifest = executed$manifest
  )
  result <- spatial_tag_coordinate_contract(result)
  results <- existing$results %||% list()
  results[[result.name]] <- result
  parameters <- list(
    group.by = group.by, species = species, database = database,
    assay = input$assay, layer = layer,
    image = input$source$image %||% image, coord.cols = coord.cols,
    coordinate_space = "raw", distance.threshold = distance.threshold,
    distance_units = "raw coordinate units", cluster = cluster,
    direction = direction, communication.args = communication.args,
    cluster.args = cluster.args, direction.args = direction.args,
    result.name = result.name, envname = envname, output.dir = output.dir,
    store.h5ad = store.h5ad, backend = backend,
    backend_scope = "unified CCC result aggregation"
  )
  versions <- unlist(executed$manifest$versions %||% character(), use.names = TRUE)
  bundle <- list(
    method = "COMMOT", results = results, active_result = result.name,
    long_table = long_table, primary_table = long_table,
    cells = input$cells, coordinates = input$coordinates,
    parameters = parameters,
    summary = list(
      result_name = result.name, n_interactions = nrow(long_table),
      n_cluster_rows = nrow(executed$cluster_table),
      n_direction_rows = nrow(executed$direction_table),
      h5ad_stored = isTRUE(store.h5ad)
    ),
    provenance = list(
      producer = "RunCOMMOT", backend_id = "commot",
      backend_version = as.character(versions[["commot"]] %||% NA_character_),
      repository = "zcang/COMMOT", commit = .commot_backend_commit,
      manifest = executed$manifest
    )
  )
  bundle <- spatial_tag_coordinate_contract(bundle)
  validate_result_bundle(bundle, label = "COMMOT")
  srt@tools[["COMMOT"]] <- bundle
  srt <- ccc_update_unified_bundle(srt, method = "COMMOT", bundle = bundle, backend = backend)
  log_message("COMMOT analysis completed", message_type = "success", verbose = verbose)
  srt
}


commot_plot_object <- function(object, stored) {
  bundle <- stored$bundle
  spatial_require_coordinate_contract(stored$result, "RunCOMMOT()")
  bundle$active_result <- stored$result.name
  bundle$long_table <- stored$result$long_table
  bundle$primary_table <- stored$result$long_table
  object@tools[["COMMOT"]] <- spatial_tag_coordinate_contract(bundle)
  ccc_update_unified_bundle(object, method = "COMMOT", bundle = bundle, backend = "r")
}

commot_select_key <- function(table, key, label) {
  keys <- unique(as.character(table$key))
  keys <- keys[!is.na(keys) & nzchar(keys)]
  if (length(keys) == 0L) {
    log_message("Stored COMMOT {.val {label}} has no selectable keys", message_type = "error")
  }
  if (is.null(key) && length(keys) > 1L) {
    log_message("Stored COMMOT {.val {label}} has multiple keys; select {.arg key}", message_type = "error")
  }
  key <- key %||% keys[[1L]]
  if (!key %in% keys) log_message("Unknown COMMOT {.arg key}: {.val {key}}", message_type = "error")
  table[as.character(table$key) == key, , drop = FALSE]
}

#' @title Plot stored COMMOT results
#'
#' @description Plot stored COMMOT group communication, cluster matrices, or
#' raw-coordinate direction vectors without rerunning Python.
#'
#' @param object A `Seurat` object with COMMOT results.
#' @param result.name Stored COMMOT result name.
#' @param plot_type Network, cluster matrix, or direction-vector view.
#' @param key Stored cluster or direction selection key.
#' @param ... Arguments passed to the SCOP network plot for `"network"`.
#'
#' @return A `ggplot` or compatible plot object.
#' @export
COMMOTPlot <- function(
  object,
  result.name = NULL,
  plot_type = c("network", "matrix", "direction"),
  key = NULL,
  ...
) {
  plot_type <- match.arg(plot_type)
  stored <- tool_bundle_get_result(object, "COMMOT", result.name)
  spatial_require_coordinate_contract(stored$result, "RunCOMMOT()")
  if (identical(plot_type, "network")) {
    plot_object <- commot_plot_object(object, stored)
    return(do.call(CCCNetworkPlot, c(list(srt = plot_object, method = "COMMOT", plot_type = "circle"), list(...))))
  }
  if (identical(plot_type, "matrix")) {
    table <- stored$result$cluster_table
    if (!is.data.frame(table) || nrow(table) == 0L) {
      log_message("The stored COMMOT result has no cluster communication", message_type = "error")
    }
    table <- commot_select_key(table, key, "cluster output")
    return(
      ggplot2::ggplot(table, ggplot2::aes(
        x = .data[["receiver"]], y = .data[["sender"]], fill = .data[["score"]]
      )) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::labs(x = "Receiver", y = "Sender", fill = "Score") +
        ggplot2::theme_bw()
    )
  }
  table <- stored$result$direction_table
  if (!is.data.frame(table) || nrow(table) == 0L) {
    log_message("The stored COMMOT result has no communication directions", message_type = "error")
  }
  table <- commot_select_key(table, key, "direction output")
  ggplot2::ggplot(table, ggplot2::aes(x = .data[["x"]], y = .data[["y"]])) +
    ggplot2::geom_point(ggplot2::aes(color = .data[["group"]]), alpha = 0.45) +
    ggplot2::geom_segment(ggplot2::aes(
      xend = .data[["x"]] + .data[["dx"]],
      yend = .data[["y"]] + .data[["dy"]]
    ), arrow = grid::arrow(length = grid::unit(0.08, "inches")), alpha = 0.65) +
    ggplot2::facet_wrap(~perspective) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "Raw x", y = "Raw y", color = "Group") +
    ggplot2::theme_bw()
}
