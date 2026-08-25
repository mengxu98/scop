# Official SpatialDM producer, result accessor, and SCOP-style plots -----------

.spatialdm_backend_commit <- "9b0f559dfd152c361fdd311b129cda84692349f3"
.spatialdm_repository <- "StatBiomed/SpatialDM"


spatialdm_validate_matrix <- function(x, label, nonnegative = TRUE) {
  if ((!is.matrix(x) && !methods::is(x, "sparseMatrix")) ||
      is.null(rownames(x)) || is.null(colnames(x)) ||
      anyNA(rownames(x)) || anyNA(colnames(x)) ||
      any(!nzchar(rownames(x))) || any(!nzchar(colnames(x))) ||
      anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
    log_message("{.val {label}} must be a named matrix with unique feature and spot IDs", message_type = "error")
  }
  values <- matrix_values(x)
  if (any(!is.finite(values)) || (isTRUE(nonnegative) && any(values < 0))) {
    log_message("{.val {label}} must contain finite non-negative values", message_type = "error")
  }
  invisible(TRUE)
}

spatialdm_input <- function(
  srt, assay, layer, counts.layer, image, coord.cols, run.local
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!assay %in% SeuratObject::Assays(srt)) {
    log_message("Assay {.val {assay}} is absent from {.arg srt}", message_type = "error")
  }
  if (identical(layer, "scale.data") || identical(counts.layer, "scale.data")) {
    log_message("SpatialDM expects normalized/log-transformed expression, not {.val scale.data}", message_type = "error")
  }
  coords <- SpatialCoordinates(
    object = srt, image = image, coord.cols = coord.cols,
    space = "raw", image_policy = "strict"
  )
  cells <- as.character(coords$data$cell_id)
  if (length(cells) < 3L || anyDuplicated(cells) || any(!cells %in% colnames(srt))) {
    log_message("Spatial coordinates must identify at least three unique Seurat spots", message_type = "error")
  }
  coords$data <- coords$data[match(cells, coords$data$cell_id), , drop = FALSE]
  expression <- GetAssayData5(srt, assay = assay, layer = layer)[, cells, drop = FALSE]
  spatialdm_validate_matrix(expression, "normalized expression")
  counts <- NULL
  if (isTRUE(run.local)) {
    counts <- GetAssayData5(srt, assay = assay, layer = counts.layer)[, cells, drop = FALSE]
    spatialdm_validate_matrix(counts, "count expression")
  }
  list(
    expression = expression, counts = counts, coordinates = coords$data,
    cells = cells, source = coords$source, assay = assay, layer = layer,
    counts.layer = if (isTRUE(run.local)) counts.layer else NULL,
    transform = coords$transform
  )
}


spatialdm_write_inputs <- function(input, workdir, lr.database = NULL) {
  expression.path <- file.path(workdir, "expression.mtx")
  counts.path <- file.path(workdir, "counts.mtx")
  features.path <- file.path(workdir, "features.tsv")
  cells.path <- file.path(workdir, "cells.tsv")
  coordinates.path <- file.path(workdir, "coordinates.tsv")
  Matrix::writeMM(input$expression, expression.path)
  Matrix::writeMM(input$counts %||% input$expression, counts.path)
  write_tsv(data.frame(feature_id = rownames(input$expression)), features.path)
  write_tsv(data.frame(cell_id = input$cells), cells.path)
  write_tsv(data.frame(
    cell_id = input$cells,
    x = as.numeric(input$coordinates$x), y = as.numeric(input$coordinates$y)
  ), coordinates.path)
  lr.path <- NULL
  if (!is.null(lr.database)) {
    lr <- as.data.frame(lr.database, stringsAsFactors = FALSE, check.names = FALSE)
    required <- c("ligand", "receptor", "annotation")
    if (!all(required %in% colnames(lr)) || nrow(lr) == 0L) {
      log_message("{.arg lr.database} must contain ligand, receptor, and annotation columns", message_type = "error")
    }
    for (column in c("ligand", "receptor")) {
      lr[[column]] <- vapply(lr[[column]], function(value) {
        if (is.list(value)) paste(as.character(value), collapse = ";") else as.character(value)
      }, character(1))
    }
    lr.path <- file.path(workdir, "custom_lr.tsv")
    write_tsv(lr, lr.path)
  }
  list(
    expression_path = normalizePath(expression.path, winslash = "/", mustWork = TRUE),
    counts_path = normalizePath(counts.path, winslash = "/", mustWork = TRUE),
    features_path = normalizePath(features.path, winslash = "/", mustWork = TRUE),
    cells_path = normalizePath(cells.path, winslash = "/", mustWork = TRUE),
    coordinates_path = normalizePath(coordinates.path, winslash = "/", mustWork = TRUE),
    custom_lr_path = if (is.null(lr.path)) NULL else normalizePath(lr.path, winslash = "/", mustWork = TRUE)
  )
}

spatialdm_runtime <- function(envname, verbose) {
  PrepareEnv(envname = envname, modules = "spatialdm", verbose = verbose)
  requirements <- env_requirements(modules = "spatialdm", verbose = FALSE)$packages
  if (!isTRUE(check_python(requirements, envname = envname, verbose = FALSE))) {
    log_message("The isolated {.file {envname}} environment is missing official SpatialDM dependencies", message_type = "error")
  }
  python <- resolve_python_executable(envname = envname, conda = resolve_conda("auto"))
  normalizePath(python, winslash = "/", mustWork = TRUE)
}

spatialdm_read_table <- function(path, required = character(), label, allow.empty = FALSE, sep = "\t") {
  if (!file.exists(path)) log_message("SpatialDM did not produce {.val {label}} output", message_type = "error")
  value <- tryCatch(utils::read.table(path, header = TRUE, sep = sep, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA", "NaN", "nan")),
    error = function(e) log_message("Unable to read SpatialDM {.val {label}} output: {conditionMessage(e)}", message_type = "error"))
  if (!all(required %in% colnames(value)) || (!allow.empty && nrow(value) == 0L)) {
    log_message("SpatialDM returned invalid {.val {label}} output", message_type = "error")
  }
  value
}

spatialdm_read_matrix <- function(path, cells = NULL, label = "matrix", logical = FALSE) {
  table <- spatialdm_read_table(path, required = character(), label = label, allow.empty = TRUE)
  if (ncol(table) == 0L) return(matrix(numeric(), nrow = 0L, ncol = 0L))
  ids <- as.character(table[[1L]])
  table[[1L]] <- NULL
  out <- as.matrix(table)
  suppressWarnings(storage.mode(out) <- "numeric")
  rownames(out) <- ids
  if (!is.null(cells)) {
    if (!identical(colnames(out), as.character(cells))) log_message("SpatialDM {.val {label}} spot IDs are misaligned", message_type = "error")
  }
  if (isTRUE(logical)) out <- out != 0
  out
}

spatialdm_read_sparse <- function(path, cells, label) {
  if (!file.exists(path)) log_message("SpatialDM did not produce {.val {label}} weights", message_type = "error")
  out <- Matrix::readMM(path)
  if (!all(dim(out) == c(length(cells), length(cells)))) log_message("SpatialDM {.val {label}} dimensions do not match spots", message_type = "error")
  dimnames(out) <- list(as.character(cells), as.character(cells))
  out
}

spatialdm_validate_manifest <- function(manifest, input, global) {
  cell_ids <- as.character(unlist(manifest$cell_ids %||% character(), use.names = FALSE))
  feature_ids <- as.character(unlist(manifest$feature_ids %||% character(), use.names = FALSE))
  valid <- identical(as.character(manifest$backend), "SpatialDM") &&
    identical(as.character(manifest$coordinate_space), "raw") &&
    identical(as.integer(manifest$n_obs), length(input$cells)) &&
    identical(as.integer(manifest$n_vars), nrow(input$expression)) &&
    identical(cell_ids, input$cells) &&
    identical(feature_ids, rownames(input$expression)) &&
    identical(as.integer(manifest$n_candidate_pairs), nrow(global))
  if (!isTRUE(valid)) {
    log_message(
      "SpatialDM manifest does not match the submitted matrix, IDs, coordinates, or global output",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

spatialdm_execute <- function(input, species, lr.database, parameters, result.name, envname, verbose) {
  python <- spatialdm_runtime(envname, verbose)
  workdir <- tempfile("spatialdm_run_")
  dir.create(workdir, recursive = TRUE)
  on.exit(unlink(workdir, recursive = TRUE, force = TRUE), add = TRUE)
  output <- file.path(workdir, "output")
  dir.create(output)
  paths <- spatialdm_write_inputs(input, workdir, lr.database = lr.database)
  config <- c(paths, list(
    output_dir = normalizePath(output, winslash = "/", mustWork = TRUE),
    species = species, method = parameters$method, run.local = parameters$run.local,
    l = parameters$l, eff_dist = parameters$eff_dist, cutoff = parameters$cutoff,
    n_neighbors = parameters$n_neighbors, n_neighbor_layers = parameters$n_neighbor_layers,
    single_cell = parameters$single_cell, complex_mean = parameters$complex.mean,
    min_cell = parameters$min_cell, n_perm = parameters$n_perm, nproc = parameters$nproc,
    seed = parameters$seed, global_fdr = parameters$global.fdr,
    global_threshold = parameters$global.threshold, local_fdr = parameters$local.fdr,
    local_threshold = parameters$local.threshold
  ))
  config.path <- file.path(workdir, "config.json")
  runner_write_json(config, config.path)
  stdout.path <- file.path(workdir, "stdout.log")
  stderr.path <- file.path(workdir, "stderr.log")
  runner <- runner_script_path("spatialdm_runner.py", "SpatialDM")
  status <- runner_system2(
    command = python, args = c(shQuote(runner), "--config", shQuote(config.path)),
    env = c(PYTHONNOUSERSITE = "1"), stdout = stdout.path, stderr = stderr.path
  )
  if (!identical(status, 0L)) runner_error(status, stdout.path, stderr.path, backend = "SpatialDM")
  manifest <- runner_read_json(file.path(output, "manifest.json"))
  global <- spatialdm_read_table(file.path(output, "global.csv"),
    required = c("interaction", "moran_r", "z_score", "pvalue", "selected"), label = "global", allow.empty = FALSE, sep = ",")
  numeric.cols <- intersect(c("moran_r", "z_score", "pvalue", "z_pval", "perm_pval", "fdr", "n_spots"), colnames(global))
  for (column in numeric.cols) global[[column]] <- suppressWarnings(as.numeric(global[[column]]))
  global$selected <- as.logical(global$selected)
  spatialdm_validate_manifest(manifest, input, global)
  cells <- input$cells
  local <- list(
    local_i = spatialdm_read_matrix(file.path(output, "local_i.tsv"), cells, "local ligand statistic"),
    local_i_r = spatialdm_read_matrix(file.path(output, "local_i_r.tsv"), cells, "local receptor statistic"),
    local_z = spatialdm_read_matrix(file.path(output, "local_z.tsv"), cells, "local z statistic"),
    local_p = spatialdm_read_matrix(file.path(output, "local_p.tsv"), cells, "local p value"),
    selected_spots = spatialdm_read_matrix(file.path(output, "selected_spots.tsv"), cells, "selected spots", logical = TRUE)
  )
  coordinates <- spatialdm_read_table(file.path(output, "coordinates.tsv"), required = c("cell_id", "x", "y"), label = "coordinates")
  if (!identical(as.character(coordinates$cell_id), cells)) log_message("SpatialDM coordinate IDs are misaligned", message_type = "error")
  local$local_i <- local$local_i[global$interaction[global$selected], , drop = FALSE]
  local$local_i_r <- local$local_i_r[global$interaction[global$selected], , drop = FALSE]
  local$local_z <- local$local_z[global$interaction[global$selected], , drop = FALSE]
  local$local_p <- local$local_p[global$interaction[global$selected], , drop = FALSE]
  local$selected_spots <- local$selected_spots[global$interaction[global$selected], , drop = FALSE]
  list(global = global, local = local, coordinates = coordinates, weights = list(
    secreted = spatialdm_read_sparse(file.path(output, "weight.mtx"), cells, "secreted"),
    contact = spatialdm_read_sparse(file.path(output, "nearest_neighbors.mtx"), cells, "contact")
  ), manifest = manifest)
}

spatialdm_bundle <- function(input, executed, parameters, result.name) {
  global <- executed$global
  local <- executed$local
  result <- list(
    schema_version = 1L, method = "SpatialDM",
    result_type = "spatial_ligand_receptor_association",
    analysis_levels = c("lr_pair_global", "lr_spot_local"),
    source = list(spatial = input$source, transform = input$transform),
    cells = input$cells, features = rownames(input$expression),
    global = global, local = local, coordinates = executed$coordinates,
    weights = executed$weights, parameters = parameters,
    lr_database = list(
      name = executed$manifest$database_name %||% "CellChatDB",
      source = executed$manifest$database_source %||% "SpatialDM package",
      species = parameters$species,
      n_candidate_pairs = executed$manifest$n_candidate_pairs %||% nrow(global),
      n_selected_pairs = executed$manifest$n_selected_pairs %||% sum(global$selected),
      min_cell = parameters$min_cell, complex_mean = parameters$complex.mean
    ),
    summary = list(
      n_spots = length(input$cells), n_features = nrow(input$expression),
      n_candidate_pairs = nrow(global), n_selected_pairs = sum(global$selected),
      n_selected_spots = if (length(local$selected_spots)) sum(local$selected_spots) else 0L
    ),
    provenance = list(
      producer = "RunSpatialDM", backend_id = "SpatialDM",
      backend_version = executed$manifest$versions$SpatialDM %||% NA_character_,
      repository = .spatialdm_repository, commit = .spatialdm_backend_commit,
      manifest = executed$manifest
    )
  )
  result <- spatial_tag_coordinate_contract(result)
  spatial_tag_coordinate_contract(list(
    method = "SpatialDM", result_type = result$result_type,
    results = stats::setNames(list(result), result.name),
    active_result = result.name, result = result, global = global, local = local,
    source = result$source, cells = result$cells, features = result$features,
    coordinates = executed$coordinates, weights = executed$weights,
    lr_database = result$lr_database,
    parameters = parameters, summary = result$summary, provenance = result$provenance
  ))
}

#' @title Run official SpatialDM spatial ligand-receptor association analysis
#' @md
#' @inheritParams scop-params
#' @param srt A `Seurat` spatial object.
#' @param species Built-in SpatialDM species (`"human"` or `"mouse"`).
#' @param lr.database Optional custom LR table with `ligand`, `receptor`, and
#'   `annotation` columns.
#' @param assay,layer,counts.layer Assay and normalized/count expression layers.
#' @param image,coord.cols Spatial image or metadata coordinate columns.
#' @param method Official global/local test mode.
#' @param l,eff_dist SpatialDM RBF scale; one must be supplied.
#' @return A Seurat object with a schema-v1 `SpatialDM` result bundle.
#' @export
RunSpatialDM <- function(
  srt, species = c("human", "mouse"), lr.database = NULL, assay = NULL,
  layer = "data", counts.layer = "counts", image = NULL,
  coord.cols = c("col", "row"), coordinate.unit = NULL,
  method = c("z-score", "permutation"), run.local = TRUE, l = NULL,
  eff_dist = NULL, cutoff = 0.1, n_neighbors = NULL, n_neighbor_layers = 6,
  single_cell = FALSE, complex.mean = c("algebra", "geometric"), min_cell = 3,
  n_perm = 1000, nproc = 1, seed = 1, global.fdr = TRUE,
  global.threshold = 0.1, local.fdr = FALSE, local.threshold = 0.1,
  result.name = "default", envname = "spatialdm_env",
  overwrite = FALSE, verbose = TRUE
) {
  species <- match.arg(species)
  method <- match.arg(method)
  complex.mean <- match.arg(complex.mean)
  validate_scalar_flag(run.local, "run.local")
  validate_scalar_flag(single_cell, "single_cell")
  validate_scalar_flag(global.fdr, "global.fdr")
  validate_scalar_flag(local.fdr, "local.fdr")
  validate_scalar_flag(overwrite, "overwrite")
  if (is.null(l) && is.null(eff_dist)) log_message("One of {.arg l} and {.arg eff_dist} must be supplied", message_type = "error")
  if (!is.null(l) && (!is.numeric(l) || length(l) != 1L || !is.finite(l) || l <= 0)) log_message("{.arg l} must be one positive finite number", message_type = "error")
  if (!is.null(eff_dist) && (!is.numeric(eff_dist) || length(eff_dist) != 1L || !is.finite(eff_dist) || eff_dist <= 0)) log_message("{.arg eff_dist} must be one positive finite number", message_type = "error")
  if (!is.numeric(cutoff) || length(cutoff) != 1L || !is.finite(cutoff) || cutoff <= 0 || cutoff >= 1) log_message("{.arg cutoff} must be between zero and one", message_type = "error")
  if (!is.numeric(min_cell) || length(min_cell) != 1L || min_cell < 0) log_message("{.arg min_cell} must be non-negative", message_type = "error")
  if (!is.numeric(n_perm) || length(n_perm) != 1L || !is.finite(n_perm) || n_perm < 1) log_message("{.arg n_perm} must be one positive number", message_type = "error")
  if (!is.numeric(nproc) || length(nproc) != 1L || !is.finite(nproc) || nproc < 1) log_message("{.arg nproc} must be one positive number", message_type = "error")
  if (!is.numeric(n_neighbor_layers) || length(n_neighbor_layers) != 1L || !is.finite(n_neighbor_layers) || n_neighbor_layers < 1) log_message("{.arg n_neighbor_layers} must be one positive number", message_type = "error")
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) log_message("{.arg seed} must be one finite number", message_type = "error")
  if (!is.null(n_neighbors) && (!is.numeric(n_neighbors) || length(n_neighbors) != 1L || !is.finite(n_neighbors) || n_neighbors < 2)) log_message("{.arg n_neighbors} must be at least two", message_type = "error")
  for (threshold in c(global.threshold, local.threshold)) if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold) || threshold <= 0 || threshold >= 1) log_message("Selection thresholds must be between zero and one", message_type = "error")
  input <- spatialdm_input(srt, assay, layer, counts.layer, image, coord.cols, run.local)
  existing <- srt@tools[["SpatialDM"]]
  if (!is.null(existing) && !isTRUE(overwrite) && result.name %in% names(existing$results %||% list())) log_message("SpatialDM result {.val {result.name}} already exists; set {.arg overwrite = TRUE}", message_type = "error")
  parameters <- list(
    species = species, assay = input$assay, layer = input$layer,
    counts.layer = input$counts.layer, image = input$source$image,
    coordinate_space = "raw", coordinate_unit = coordinate.unit %||% "raw coordinate units",
    coord.cols = input$source$coord.cols, method = method, run.local = run.local,
    l = l, eff_dist = eff_dist, cutoff = cutoff, n_neighbors = n_neighbors,
    n_neighbor_layers = n_neighbor_layers, single_cell = single_cell,
    complex.mean = complex.mean, min_cell = min_cell, n_perm = n_perm,
    nproc = nproc, seed = seed, global.fdr = global.fdr,
    global.threshold = global.threshold, local.fdr = local.fdr,
    local.threshold = local.threshold
  )
  executed <- spatialdm_execute(input, species, lr.database, parameters, result.name, envname, verbose)
  bundle <- spatialdm_bundle(input, executed, parameters, result.name)
  if (is.null(existing) || isTRUE(overwrite)) {
    if (!is.null(existing) && isTRUE(overwrite)) bundle$results <- c(existing$results[setdiff(names(existing$results), result.name)], bundle$results)
    srt@tools[["SpatialDM"]] <- bundle
  } else {
    existing$results <- c(existing$results, bundle$results)
    existing$active_result <- result.name
    existing$result <- bundle$result
    existing$global <- bundle$global
    existing$local <- bundle$local
    existing$source <- bundle$source
    existing$cells <- bundle$cells
    existing$features <- bundle$features
    existing$coordinates <- bundle$coordinates
    existing$weights <- bundle$weights
    existing$lr_database <- bundle$lr_database
    existing$parameters <- bundle$parameters
    existing$summary <- bundle$summary
    existing$provenance <- bundle$provenance
    srt@tools[["SpatialDM"]] <- existing
  }
  log_message("SpatialDM analysis completed", message_type = "success", verbose = verbose)
  srt
}

spatialdm_get_result <- function(object, result.name = NULL) {
  if (!inherits(object, "Seurat")) log_message("{.arg object} must be a {.cls Seurat} object", message_type = "error")
  bundle <- object@tools[["SpatialDM"]]
  if (is.null(bundle)) log_message("SpatialDM results are absent", message_type = "error")
  if (is.null(result.name) && length(bundle$results) > 1L) log_message("Multiple SpatialDM results are stored; select {.arg result.name}", message_type = "error")
  result.name <- result.name %||% bundle$active_result
  result <- bundle$results[[result.name]]
  if (is.null(result)) log_message("Unknown SpatialDM result {.val {result.name}}", message_type = "error")
  spatial_require_coordinate_contract(result, "RunSpatialDM()")
  result
}


#' @title Access stored SpatialDM results
#' @param object A Seurat object with SpatialDM results.
#' @param result.name Stored result name.
#' @param type Result type to return (`"global"`, `"local"`, or `"weights"`).
#' @param pair Optional LR interaction name for local results.
#' @return The requested stored SpatialDM result payload.
#' @export
GetSpatialDMResult <- function(object, result.name = NULL, type = c("global", "local", "weights"), pair = NULL) {
  type <- match.arg(type)
  result <- spatialdm_get_result(object, result.name)
  if (identical(type, "global")) return(result$global)
  if (identical(type, "weights")) return(result$weights)
  if (is.null(pair)) return(result$local)
  if (!pair %in% rownames(result$local$local_p)) log_message("Unknown SpatialDM {.arg pair}: {.val {pair}}", message_type = "error")
  lapply(result$local, function(x) if (is.matrix(x)) x[pair, , drop = TRUE] else x)
}


#' @title Plot stored SpatialDM results
#' @param object A Seurat object with SpatialDM results.
#' @param result.name Stored result name.
#' @param plot_type Plot type (`"weights"`, `"global"`, or `"local"`).
#' @param pair LR interaction name for local plotting.
#' @param spot Spot identifier for weight plotting.
#' @param signaling Which stored weight matrix to use.
#' @param image.scale Image scale factor matching the selected raster.
#' @param highlight LR interactions to label in the global plot.
#' @param palette,palcolor,theme_use,theme_args SCOP plotting controls.
#' @param ... Additional arguments passed to SCOP spatial plotting functions.
#' @return A `ggplot` or `patchwork` object.
#' @export
SpatialDMPlot <- function(
  object, result.name = NULL, plot_type = c("weights", "global", "local"),
  pair = NULL, spot = NULL, signaling = c("secreted", "contact"),
  highlight = NULL, palette = "RdBu", palcolor = NULL,
  theme_use = "theme_scop", theme_args = list(),
  ..., image.scale = c("lowres", "hires")
) {
  plot_type <- match.arg(plot_type)
  image.scale <- match.arg(image.scale)
  signaling <- match.arg(signaling)
  result <- spatialdm_get_result(object, result.name)
  spatial_require_coordinate_contract(result, "RunSpatialDM()")
  if (identical(plot_type, "weights")) {
    if (is.null(spot) || length(spot) != 1L || !spot %in% as.character(result$coordinates$cell_id)) log_message("{.arg spot} must identify one stored SpatialDM spot", message_type = "error")
    values <- as.numeric(result$weights[[signaling]][spot, , drop = TRUE])
    names(values) <- colnames(result$weights[[signaling]])
    return(SpatialSpotPlot(object, values = values, image = result$parameters$image,
      image.scale = image.scale,
      coord.cols = result$parameters$coord.cols, palette = palette, palcolor = palcolor,
      theme_use = theme_use, theme_args = theme_args, legend.title = "Spatial weight", ...))
  }
  if (identical(plot_type, "global")) {
    df <- result$global
    df$selected <- as.logical(df$selected)
    global_p_label <- if ("fdr" %in% colnames(df)) "-log10(FDR)" else "-log10(p-value)"
    df$minus_log10_fdr <- -log10(pmax(as.numeric(df$fdr %||% df$pvalue), .Machine$double.xmin))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$moran_r, y = .data$minus_log10_fdr)) +
      ggplot2::geom_point(ggplot2::aes(color = .data$selected), alpha = 0.75) +
      ggplot2::scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#D55E00"), name = "Selected") +
      ggplot2::labs(x = "Global Moran's R", y = global_p_label, title = "Spatial association (Moran's R)") +
      apply_plot_theme(theme_use = theme_use, theme_args = theme_args)
    if (!is.null(highlight)) p <- p + ggrepel::geom_text_repel(data = df[df$interaction %in% highlight, , drop = FALSE], ggplot2::aes(label = .data$interaction), size = 3, max.overlaps = Inf)
    return(p)
  }
  if (is.null(pair) || length(pair) != 1L || !pair %in% rownames(result$local$local_p)) log_message("{.arg pair} must identify one stored local LR result", message_type = "error")
  pvalue <- as.numeric(result$local$local_p[pair, ])
  names(pvalue) <- colnames(result$local$local_p)
  interaction_row <- result$global[result$global$interaction == pair, , drop = FALSE]
  ligands <- as.character(interaction_row[1, grep("^Ligand", colnames(interaction_row), value = TRUE)])
  receptors <- as.character(interaction_row[1, grep("^Receptor", colnames(interaction_row), value = TRUE)])
  ligands <- ligands[!is.na(ligands) & nzchar(ligands)]
  receptors <- receptors[!is.na(receptors) & nzchar(receptors)]
  local_label <- if (isTRUE(result$parameters$local.fdr)) "1 - local FDR" else "1 - local p"
  plots <- list(SpatialSpotPlot(object, values = 1 - pvalue, image = result$parameters$image,
    image.scale = image.scale,
    coord.cols = result$parameters$coord.cols, palette = "Reds", theme_use = theme_use,
    theme_args = theme_args, legend.title = local_label, ...))
  for (gene in c(ligands, receptors)) if (gene %in% rownames(object)) plots[[length(plots) + 1L]] <- SpatialSpotPlot(object, features = gene,
    assay = result$parameters$assay, layer = result$parameters$layer, image = result$parameters$image,
    image.scale = image.scale,
    coord.cols = result$parameters$coord.cols, palette = palette, palcolor = palcolor,
    theme_use = theme_use, theme_args = theme_args, legend.title = gene, ...)
  patchwork::wrap_plots(plots, ncol = length(plots)) + patchwork::plot_annotation(title = pair)
}
