#' @title Run CHOIR clustering
#'
#' @description
#' Runs the optional CHOIR backend on a single-modality `Seurat` object.
#' CHOIR builds and prunes a hierarchical clustering tree using random-forest
#' classifiers and permutation tests to identify statistically distinct cell
#' populations.
#'
#' The upstream CHOIR records remain available in `srt@misc[[key]]`. This
#' wrapper also writes a stable cluster column to cell metadata and a lightweight
#' summary to `srt@tools[[tool_name]]`.
#'
#' @md
#' @param srt A `Seurat` object.
#' @param assay Assay used by CHOIR. If `NULL`, the default assay is used.
#' @param layer Assay layer used by CHOIR. If `NULL`, `"data"` is used for
#'   `RNA` and `sketch` assays, and `"scale.data"` is used for `SCT` and
#'   `integrated` assays. Other assays require an explicit layer.
#' @param key Name used by CHOIR to store its complete results in `srt@misc`.
#' @param cluster_colname Metadata column used for the final CHOIR clusters.
#' @param tool_name Name used to store the SCOP result summary in `srt@tools`.
#' @param alpha Significance level used for CHOIR permutation tests.
#' @param p_adjust Multiple-testing correction used by CHOIR.
#' @param feature_set Whether CHOIR random forests use variable or all features.
#' @param exclude_features Features excluded from CHOIR random forests.
#' @param n_iterations Number of bootstrap iterations for each permutation test.
#' @param n_trees Number of trees in each random forest.
#' @param min_accuracy Minimum classifier accuracy required to keep clusters
#'   separate.
#' @param max_clusters Must be `"auto"`. Numeric limits are rejected because
#'   the pinned upstream backend can fail to terminate when the number of
#'   clusters plateaus.
#' @param normalization_method Normalization performed inside CHOIR. Use
#'   `"none"` for previously normalized data or `"SCTransform"` with a counts
#'   layer. The pinned backend does not support `"SCTransform"` for Seurat v5
#'   `Assay5` objects.
#' @param batch.by Optional metadata column containing batch labels.
#' @param batch_correction_method Batch correction performed by CHOIR. If
#'   `NULL`, `"Harmony"` is selected when `batch.by` is supplied and `"none"`
#'   otherwise.
#' @param reduction Optional existing dimensional reduction supplied to CHOIR.
#'   This can be a reduction name in `srt` or a cell-by-dimension matrix.
#' @param var_features Features associated with `reduction`. If `NULL` and a
#'   reduction is supplied, variable features from `assay` are used.
#' @param atac Whether the selected assay contains ATAC-seq data.
#' @param n_cores Number of cores used by CHOIR. The pinned backend supports
#'   macOS and Linux; Windows execution is rejected before installation.
#' @param seed Random seed passed to CHOIR.
#' @param store_tool Whether to store a lightweight result summary in
#'   `srt@tools[[tool_name]]`.
#' @param verbose Whether to print progress messages.
#' @param overwrite Whether to replace existing CHOIR metadata, reduction, and
#'   `misc` entries, plus the `tools` entry when `store_tool = TRUE`.
#' @param ... Additional named arguments passed to the installed CHOIR entry
#'   point. Unsupported arguments produce an error rather than being silently
#'   ignored.
#'
#' @return A `Seurat` object containing CHOIR clusters in `cluster_colname`,
#' complete upstream records in `srt@misc[[key]]`, and, when
#' `store_tool = TRUE`, a lightweight summary in `srt@tools[[tool_name]]`.
#'
#' @references
#' Sant, C. et al. CHOIR improves significance-based detection of cell types
#' and states from single-cell data. \emph{Nature Genetics} 57, 1309-1319
#' (2025). \doi{10.1038/s41588-025-02148-8}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' if (check_r("corceslab/CHOIR", verbose = FALSE)) {
#'   data(pancreas_sub)
#'   pancreas_sub <- Seurat::NormalizeData(pancreas_sub, verbose = FALSE)
#'   pancreas_sub <- RunCHOIR(
#'     pancreas_sub,
#'     assay = "RNA",
#'     n_cores = 2,
#'     verbose = FALSE
#'   )
#'   CellDimPlot(pancreas_sub, group.by = "CHOIR_cluster")
#' }
#' }
RunCHOIR <- function(
  srt,
  assay = NULL,
  layer = NULL,
  key = "CHOIR",
  cluster_colname = "CHOIR_cluster",
  tool_name = "CHOIR",
  alpha = 0.05,
  p_adjust = c("bonferroni", "fdr", "none"),
  feature_set = c("var", "all"),
  exclude_features = NULL,
  n_iterations = 100,
  n_trees = 50,
  min_accuracy = 0.5,
  max_clusters = "auto",
  normalization_method = c("none", "SCTransform"),
  batch.by = NULL,
  batch_correction_method = NULL,
  reduction = NULL,
  var_features = NULL,
  atac = FALSE,
  n_cores = 1,
  seed = 1,
  store_tool = TRUE,
  verbose = TRUE,
  overwrite = FALSE,
  ...
) {
  validate_seurat_object(srt)
  validate_scalar_string(key, "key")
  validate_scalar_string(cluster_colname, "cluster_colname")
  validate_scalar_string(tool_name, "tool_name")
  validate_scalar_flag(atac, "atac")
  validate_scalar_flag(store_tool, "store_tool")
  validate_scalar_flag(verbose, "verbose")
  validate_scalar_flag(overwrite, "overwrite")

  p_adjust <- match.arg(p_adjust)
  feature_set <- match.arg(feature_set)
  normalization_method <- match.arg(normalization_method)
  max_clusters <- choir_validate_max_clusters(max_clusters)
  alpha <- choir_assert_probability(alpha, "alpha", open = TRUE)
  min_accuracy <- choir_assert_probability(min_accuracy, "min_accuracy")
  n_iterations <- validate_scalar_integer(n_iterations, "n_iterations")
  n_trees <- validate_scalar_integer(n_trees, "n_trees")
  n_cores <- validate_scalar_integer(n_cores, "n_cores")
  seed <- validate_scalar_integer(
    seed,
    "seed",
    minimum = 1L,
    message = "must be a positive integer"
  )

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  layer <- choir_resolve_layer(
    srt = srt,
    assay = assay,
    layer = layer,
    normalization_method = normalization_method
  )
  batch <- choir_resolve_batch(
    srt = srt,
    batch.by = batch.by,
    batch_correction_method = batch_correction_method,
    verbose = verbose
  )
  reduction_input <- choir_resolve_reduction(
    srt = srt,
    assay = assay,
    reduction = reduction,
    var_features = var_features
  )

  dots <- list(...)
  choir_validate_dots(dots)
  expected_cluster_col <- paste0("CHOIR_clusters_", alpha)
  if (identical(batch$labels, expected_cluster_col)) {
    log_message(
      "{.arg batch.by} cannot use CHOIR's managed cluster column {.val {expected_cluster_col}}",
      message_type = "error"
    )
  }
  backend_input <- choir_prepare_output(
    srt = srt,
    key = key,
    expected_cluster_col = expected_cluster_col,
    cluster_colname = cluster_colname,
    tool_name = tool_name,
    store_tool = store_tool,
    overwrite = overwrite
  )

  log_message(
    "Running {.pkg CHOIR} clustering",
    message_type = "running",
    verbose = verbose
  )
  backend_commit <- choir_check_r(verbose = FALSE)
  choir_fun <- choir_get_fun("CHOIR")

  args <- c(
    list(
      object = backend_input,
      key = key,
      alpha = alpha,
      p_adjust = p_adjust,
      feature_set = feature_set,
      exclude_features = exclude_features,
      n_iterations = n_iterations,
      n_trees = n_trees,
      min_accuracy = min_accuracy,
      max_clusters = max_clusters,
      normalization_method = normalization_method,
      batch_correction_method = batch$method,
      batch_labels = batch$labels,
      use_assay = assay,
      use_slot = layer,
      reduction = reduction_input$reduction,
      var_features = reduction_input$var_features,
      atac = atac,
      n_cores = n_cores,
      random_seed = seed,
      verbose = verbose
    ),
    dots
  )
  args <- choir_match_backend_args(choir_fun, args, dot_names = names(dots))
  result <- do.call(choir_fun, args)
  result <- choir_validate_result(
    result = result,
    input_cells = colnames(srt),
    input_metadata = colnames(backend_input@meta.data),
    alpha = alpha,
    key = key
  )

  backend_cluster_col <- attr(result, "choir_cluster_column")
  clusters <- factor(
    as.character(result@meta.data[[backend_cluster_col]]),
    levels = unique(as.character(result@meta.data[[backend_cluster_col]]))
  )
  names(clusters) <- colnames(result)
  result@meta.data[[cluster_colname]] <- clusters

  if (isTRUE(store_tool)) {
    result@tools[[tool_name]] <- list(
      clusters = data.frame(
        cell = colnames(result),
        cluster = as.character(clusters),
        row.names = colnames(result),
        stringsAsFactors = FALSE
      ),
      backend = list(
        package = .choir_package,
        repository = .choir_repository,
        commit = backend_commit,
        key = key,
        cluster_column = backend_cluster_col
      ),
      parameters = list(
        assay = assay,
        layer = layer,
        key = key,
        cluster_colname = cluster_colname,
        alpha = alpha,
        p_adjust = p_adjust,
        feature_set = feature_set,
        exclude_features = exclude_features,
        n_iterations = n_iterations,
        n_trees = n_trees,
        min_accuracy = min_accuracy,
        max_clusters = max_clusters,
        normalization_method = normalization_method,
        batch.by = batch$labels,
        batch_correction_method = batch$method,
        reduction = reduction_input$name,
        atac = atac,
        n_cores = n_cores,
        seed = seed,
        overwrite = overwrite
      )
    )
  }

  attr(result, "choir_cluster_column") <- NULL
  log_message(
    "{.pkg CHOIR} clusters stored in metadata column {.val {cluster_colname}}",
    message_type = "success",
    verbose = verbose
  )
  result
}

.choir_repository <- "corceslab/CHOIR"
.choir_package <- "CHOIR"
.choir_commit <- "e9ebfbc9089beeaf4ca088c7b81b18f39758b0bc"
.choir_reduction <- "CHOIR_P0_reduction"
.choir_dependencies <- c(
  "BiocGenerics",
  "bluster",
  "dplyr",
  "ggplot2",
  "ggtree",
  "harmony",
  "magrittr",
  "Matrix",
  "pengminshi/mrtree",
  "plyr",
  "progress",
  "ranger",
  "Seurat",
  "spatstat.univar",
  "stringr",
  "tidyr"
)

choir_check_r <- function(verbose = FALSE) {
  choir_validate_platform()
  if (choir_namespace_loaded()) {
    observed_commit <- choir_loaded_commit()
    if (identical(observed_commit, .choir_commit)) {
      return(invisible(.choir_commit))
    }
    log_message(
      "A different or unverifiable {.pkg CHOIR} backend is already loaded. Restart R before using the pinned backend.",
      message_type = "error"
    )
  }
  status <- check_r(
    .choir_package,
    dependencies = NA,
    install = FALSE,
    verbose = FALSE
  )
  available <- isTRUE(status) || isTRUE(unlist(status)[.choir_package])
  observed_commit <- choir_installed_commit()
  if (available && identical(observed_commit, .choir_commit)) {
    return(invisible(.choir_commit))
  }

  log_message(
    "Installing the pinned {.pkg CHOIR} backend without its upstream configure telemetry",
    verbose = verbose
  )
  choir_check_dependencies(verbose = verbose)
  choir_install_without_configure(verbose = verbose)
  status <- check_r(
    .choir_package,
    dependencies = NA,
    install = FALSE,
    verbose = FALSE
  )
  available <- isTRUE(status) || isTRUE(unlist(status)[.choir_package])
  observed_commit <- choir_installed_commit()
  if (!available || !identical(observed_commit, .choir_commit)) {
    log_message(
      "Unable to install the pinned optional {.pkg CHOIR} backend",
      message_type = "error"
    )
  }
  invisible(observed_commit)
}

choir_namespace_loaded <- function() {
  isNamespaceLoaded(.choir_package)
}

choir_loaded_commit <- function() {
  if (!choir_namespace_loaded()) {
    return(NULL)
  }
  namespace <- tryCatch(
    asNamespace(.choir_package),
    error = function(...) NULL
  )
  package_path <- tryCatch(
    getNamespaceInfo(namespace, "path"),
    error = function(...) NULL
  )
  if (is.null(package_path) || !nzchar(package_path)) {
    return(NULL)
  }
  description <- tryCatch(
    base::read.dcf(
      file.path(package_path, "DESCRIPTION"),
      fields = "RemoteSha"
    ),
    error = function(...) NULL
  )
  if (is.null(description) || length(description) != 1L) {
    return(NULL)
  }
  commit <- unname(description[[1L]])
  if (is.na(commit) || !nzchar(commit)) {
    return(NULL)
  }
  commit
}

choir_check_dependencies <- function(verbose = FALSE) {
  status <- check_r(
    .choir_dependencies,
    dependencies = NA,
    verbose = verbose
  )
  available <- unlist(status, use.names = FALSE)
  if (
    !isTRUE(status) &&
      (
        length(available) != length(.choir_dependencies) ||
          !all(available %in% TRUE)
      )
  ) {
    log_message(
      "Unable to prepare dependencies for the optional {.pkg CHOIR} backend",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

choir_installed_commit <- function() {
  description <- tryCatch(
    utils::packageDescription(.choir_package),
    error = function(...) NULL
  )
  if (is.null(description)) {
    return(NULL)
  }
  commit <- description[["RemoteSha"]]
  if (is.null(commit) || is.na(commit) || !nzchar(commit)) {
    return(NULL)
  }
  as.character(commit)
}

choir_install_without_configure <- function(verbose = FALSE) {
  workdir <- tempfile("choir_source_")
  dir.create(workdir, recursive = TRUE)
  on.exit(unlink(workdir, recursive = TRUE, force = TRUE), add = TRUE)

  archive <- file.path(workdir, "CHOIR.tar.gz")
  url <- paste0(
    "https://codeload.github.com/corceslab/CHOIR/tar.gz/",
    .choir_commit
  )
  current_timeout <- getOption("timeout", 60)
  old_options <- options(timeout = max(300, current_timeout))
  on.exit(options(old_options), add = TRUE)
  utils::download.file(url, archive, mode = "wb", quiet = !isTRUE(verbose))
  utils::untar(archive, exdir = workdir)
  descriptions <- list.files(
    workdir,
    pattern = "^DESCRIPTION$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(descriptions) != 1L) {
    log_message(
      "Unable to locate the downloaded {.pkg CHOIR} source",
      message_type = "error"
    )
  }
  description <- base::read.dcf(descriptions[[1]])
  if (!"RemoteSha" %in% colnames(description)) {
    description <- cbind(description, RemoteSha = .choir_commit)
  } else {
    description[, "RemoteSha"] <- .choir_commit
  }
  base::write.dcf(description, file = descriptions[[1]])
  source_dir <- dirname(descriptions[[1]])
  unlink(
    file.path(source_dir, c("configure", "configure.win")),
    force = TRUE
  )

  lib <- .libPaths()[[1]]
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  log <- file.path(workdir, "install.log")
  status <- system2(
    command = file.path(R.home("bin"), "R"),
    args = c(
      "CMD",
      "INSTALL",
      "--no-configure",
      paste0("--library=", shQuote(lib)),
      shQuote(source_dir)
    ),
    stdout = log,
    stderr = log
  )
  if (!identical(status, 0L)) {
    details <- runner_tail_lines(log, max_lines = 20L)
    log_message(
      paste(
        c("CHOIR installation failed", details),
        collapse = "\n"
      ),
      message_type = "error"
    )
  }
  invisible(TRUE)
}

choir_get_fun <- function(fun) {
  value <- get_namespace_fun(.choir_package, fun)
  if (!identical(choir_loaded_commit(), .choir_commit)) {
    log_message(
      "The loaded {.pkg CHOIR} namespace could not be matched to the pinned backend. Restart R and try again.",
      message_type = "error"
    )
  }
  value
}

choir_validate_platform <- function(os_type = .Platform$OS.type) {
  if (identical(os_type, "windows")) {
    log_message(
      "The pinned {.pkg CHOIR} backend supports macOS and Linux, but not Windows",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

choir_prepare_output <- function(
  srt,
  key,
  expected_cluster_col,
  cluster_colname,
  tool_name,
  store_tool,
  overwrite
) {
  metadata_names <- unique(c(expected_cluster_col, cluster_colname))
  metadata_conflicts <- intersect(metadata_names, colnames(srt@meta.data))
  conflicts <- c(
    if (length(metadata_conflicts) > 0L) {
      paste0("metadata:", metadata_conflicts)
    },
    if (.choir_reduction %in% names(srt@reductions)) {
      paste0("reduction:", .choir_reduction)
    },
    if (key %in% names(srt@misc)) paste0("misc:", key),
    if (isTRUE(store_tool) && tool_name %in% names(srt@tools)) {
      paste0("tools:", tool_name)
    }
  )
  if (length(conflicts) > 0L && !isTRUE(overwrite)) {
    log_message(
      "Existing CHOIR output would be replaced ({.val {conflicts}}). Set {.arg overwrite = TRUE} to replace it.",
      message_type = "error"
    )
  }
  if (isTRUE(overwrite)) {
    srt@meta.data[[expected_cluster_col]] <- NULL
    srt@reductions[[.choir_reduction]] <- NULL
    srt@misc[[key]] <- NULL
    if (isTRUE(store_tool)) {
      srt@tools[[tool_name]] <- NULL
    }
  }
  srt
}

choir_resolve_layer <- function(
  srt,
  assay,
  layer = NULL,
  normalization_method = "none"
) {
  if (!is.character(assay) || length(assay) != 1L || is.na(assay) ||
    !nzchar(assay)) {
    log_message(
      "{.arg assay} must be a single non-empty assay name",
      message_type = "error"
    )
  }
  if (!assay %in% SeuratObject::Assays(srt)) {
    log_message(
      "{.arg assay} {.val {assay}} is not present in {.arg srt}",
      message_type = "error"
    )
  }

  if (is.null(layer)) {
    layer <- if (assay %in% c("RNA", "sketch")) {
      "data"
    } else if (assay %in% c("SCT", "integrated")) {
      "scale.data"
    } else {
      log_message(
        "{.arg layer} must be supplied for assay {.val {assay}}",
        message_type = "error"
      )
    }
  }
  validate_scalar_string(layer, "layer")

  layers <- SeuratObject::Layers(srt[[assay]])
  if (!layer %in% layers) {
    split_layers <- grep(
      paste0("^", choir_escape_regex(layer), "\\."),
      layers,
      value = TRUE
    )
    if (length(split_layers) > 0L) {
      log_message(
        "Assay {.val {assay}} contains split {.arg {layer}} layers. Run {.code SeuratObject::JoinLayers()} before {.fn RunCHOIR}.",
        message_type = "error"
      )
    }
    log_message(
      "Layer {.val {layer}} is not present in assay {.val {assay}}",
      message_type = "error"
    )
  }

  expression <- tryCatch(
    GetAssayData5(srt, assay = assay, layer = layer),
    error = function(e) {
      log_message(
        "Unable to read layer {.val {layer}} from assay {.val {assay}}: {conditionMessage(e)}",
        message_type = "error"
      )
    }
  )
  if (!identical(colnames(expression), colnames(srt))) {
    log_message(
      "Layer {.val {layer}} must contain every cell in {.arg srt} in the same order",
      message_type = "error"
    )
  }
  if (nrow(expression) < 2L || ncol(expression) < 2L) {
    log_message(
      "CHOIR requires at least two features and two cells",
      message_type = "error"
    )
  }

  if (identical(normalization_method, "SCTransform") &&
    !identical(layer, "counts")) {
    log_message(
      "{.arg normalization_method = 'SCTransform'} requires {.arg layer = 'counts'}",
      message_type = "error"
    )
  }
  if (
    identical(normalization_method, "SCTransform") &&
      inherits(srt[[assay]], "Assay5")
  ) {
    log_message(
      "The pinned {.pkg CHOIR} backend does not support {.arg normalization_method = 'SCTransform'} with Seurat v5 {.cls Assay5} objects. Normalize the assay first and use {.arg normalization_method = 'none'}.",
      message_type = "error"
    )
  }
  if (identical(layer, "counts") &&
    identical(normalization_method, "none")) {
    log_message(
      "CHOIR expects normalized input when {.arg normalization_method = 'none'}. Consider normalizing first or use {.arg normalization_method = 'SCTransform'}.",
      message_type = "warning"
    )
  }

  layer
}

choir_resolve_batch <- function(
  srt,
  batch.by = NULL,
  batch_correction_method = NULL,
  verbose = TRUE
) {
  if (!is.null(batch.by)) {
    validate_scalar_string(batch.by, "batch.by")
    if (!batch.by %in% colnames(srt@meta.data)) {
      log_message(
        "{.arg batch.by} {.val {batch.by}} is not present in cell metadata",
        message_type = "error"
      )
    }
    labels <- srt@meta.data[[batch.by]]
    if (anyNA(labels) || any(!nzchar(as.character(labels)))) {
      log_message(
        "{.arg batch.by} must not contain missing or empty values",
        message_type = "error"
      )
    }
    if (length(unique(labels)) < 2L) {
      log_message(
        "{.arg batch.by} must contain at least two batches",
        message_type = "error"
      )
    }
  }

  if (is.null(batch_correction_method)) {
    batch_correction_method <- if (is.null(batch.by)) "none" else "Harmony"
  }
  if (!is.character(batch_correction_method) ||
    length(batch_correction_method) != 1L ||
    !batch_correction_method %in% c("none", "Harmony")) {
    log_message(
      "{.arg batch_correction_method} must be {.val 'none'} or {.val 'Harmony'}",
      message_type = "error"
    )
  }
  if (identical(batch_correction_method, "Harmony") && is.null(batch.by)) {
    log_message(
      "{.arg batch.by} is required when {.arg batch_correction_method = 'Harmony'}",
      message_type = "error"
    )
  }
  if (identical(batch_correction_method, "none") && !is.null(batch.by)) {
    log_message(
      "{.arg batch.by} is ignored when {.arg batch_correction_method = 'none'}",
      message_type = "warning",
      verbose = verbose
    )
    batch.by <- NULL
  }

  list(method = batch_correction_method, labels = batch.by)
}

choir_resolve_reduction <- function(
  srt,
  assay,
  reduction = NULL,
  var_features = NULL
) {
  reduction_name <- NULL
  if (is.character(reduction)) {
    validate_scalar_string(reduction, "reduction")
    if (!reduction %in% names(srt@reductions)) {
      log_message(
        "Reduction {.val {reduction}} is not present in {.arg srt}",
        message_type = "error"
      )
    }
    reduction_name <- reduction
    reduction <- SeuratObject::Embeddings(srt[[reduction]])
  }

  if (!is.null(reduction)) {
    if (!is.matrix(reduction) && !inherits(reduction, "Matrix")) {
      log_message(
        "{.arg reduction} must be a reduction name or matrix-like object",
        message_type = "error"
      )
    }
    if (is.null(rownames(reduction)) ||
      !setequal(rownames(reduction), colnames(srt))) {
      log_message(
        "{.arg reduction} row names must match all cells in {.arg srt}",
        message_type = "error"
      )
    }
    reduction <- as.matrix(reduction[colnames(srt), , drop = FALSE])
    if (any(!is.finite(reduction))) {
      log_message(
        "{.arg reduction} must contain only finite numeric values",
        message_type = "error"
      )
    }
    if (is.null(var_features)) {
      var_features <- SeuratObject::VariableFeatures(srt, assay = assay)
    }
    var_features <- unique(as.character(var_features))
    var_features <- intersect(var_features, rownames(srt[[assay]]))
    if (length(var_features) == 0L) {
      log_message(
        "{.arg var_features} must be supplied when using an existing reduction and no variable features are set",
        message_type = "error"
      )
    }
  } else if (!is.null(var_features)) {
    log_message(
      "{.arg var_features} can only be supplied together with {.arg reduction}",
      message_type = "error"
    )
  }

  list(
    reduction = reduction,
    var_features = var_features,
    name = reduction_name
  )
}

choir_validate_dots <- function(dots) {
  if (length(dots) == 0L) {
    return(invisible(TRUE))
  }
  dot_names <- names(dots)
  if (is.null(dot_names) || anyNA(dot_names) || any(!nzchar(dot_names))) {
    log_message(
      "All arguments passed through {.arg ...} must be named",
      message_type = "error"
    )
  }
  duplicates <- unique(dot_names[duplicated(dot_names)])
  if (length(duplicates) > 0L) {
    log_message(
      "Arguments passed through {.arg ...} must be unique: {.val {duplicates}}",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

choir_match_backend_args <- function(fun, args, dot_names = character()) {
  if (!is.function(fun)) {
    log_message(
      "The installed {.pkg CHOIR} package does not expose a callable {.fn CHOIR} entry point",
      message_type = "error"
    )
  }
  formal_names <- tryCatch(names(formals(fun)), error = function(e) NULL)
  required <- c("object", "key", "use_assay", "use_slot")
  if (is.null(formal_names) || !all(required %in% formal_names)) {
    log_message(
      "The installed {.pkg CHOIR} entry point is incompatible with {.fn RunCHOIR}",
      message_type = "error"
    )
  }
  unsupported <- setdiff(dot_names, formal_names)
  if (length(unsupported) > 0L && !"..." %in% formal_names) {
    log_message(
      "Unsupported CHOIR argument(s): {.val {unsupported}}",
      message_type = "error"
    )
  }
  if ("..." %in% formal_names) {
    return(args)
  }
  args[names(args) %in% formal_names]
}

choir_validate_result <- function(
  result,
  input_cells,
  input_metadata,
  alpha,
  key
) {
  if (!inherits(result, "Seurat")) {
    log_message(
      "{.pkg CHOIR} returned an object that is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  if (!identical(colnames(result), input_cells)) {
    log_message(
      "{.pkg CHOIR} changed cell identities or cell order",
      message_type = "error"
    )
  }

  expected <- paste0("CHOIR_clusters_", alpha)
  candidates <- setdiff(colnames(result@meta.data), input_metadata)
  candidates <- candidates[startsWith(candidates, "CHOIR_clusters_")]
  cluster_col <- if (expected %in% colnames(result@meta.data)) {
    expected
  } else if (length(candidates) == 1L) {
    candidates
  } else {
    log_message(
      "{.pkg CHOIR} did not return an unambiguous final cluster column",
      message_type = "error"
    )
  }
  clusters <- result@meta.data[[cluster_col]]
  labels <- if (is.atomic(clusters)) as.character(clusters) else character()
  invalid_numeric <- is.numeric(clusters) && any(!is.finite(clusters))
  if (
    length(clusters) != length(input_cells) ||
      length(labels) != length(input_cells) ||
      anyNA(clusters) ||
      any(!nzchar(labels)) ||
      invalid_numeric
  ) {
    log_message(
      "{.pkg CHOIR} returned invalid final cluster labels",
      message_type = "error"
    )
  }
  records <- result@misc[[key]]
  if (
    !is.list(records) ||
      !all(c("clusters", "records") %in% names(records))
  ) {
    log_message(
      "{.pkg CHOIR} did not return complete records under {.arg key} {.val {key}}",
      message_type = "error"
    )
  }

  attr(result, "choir_cluster_column") <- cluster_col
  result
}

choir_assert_probability <- function(x, arg, open = FALSE) {
  valid <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x >= 0 && x <= 1
  if (isTRUE(open)) {
    valid <- valid && x > 0 && x < 1
  }
  if (!valid) {
    interval <- if (isTRUE(open)) "(0, 1)" else "[0, 1]"
    log_message(
      "{.arg {arg}} must be a finite number in {.val {interval}}",
      message_type = "error"
    )
  }
  as.numeric(x)
}

choir_validate_max_clusters <- function(x) {
  if (is.character(x) && length(x) == 1L && !is.na(x) &&
    identical(x, "auto")) {
    return(x)
  }
  if (is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x >= 1 && x == floor(x)) {
    log_message(
      "Numeric {.arg max_clusters} is unsafe with the pinned {.pkg CHOIR} backend because its full-tree loop may not terminate when the cluster count plateaus. Use {.val 'auto'}.",
      message_type = "error"
    )
  }
  log_message(
    "{.arg max_clusters} must be {.val 'auto'}",
    message_type = "error"
  )
}

choir_escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}
