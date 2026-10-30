#' @title Run official CellOracle perturbation analysis
#'
#' @description
#' Fit cluster-specific regulatory models and simulate gene perturbations using
#' the unmodified CellOracle 0.20.0 Python implementation.
#' @md
#'
#' @param object A Seurat object with non-log RNA counts, clusters and a two-dimensional
#' embedding. Perform cell/gene selection before calling this function.
#' @param perturbation Named numeric vector of non-negative expression values.
#' Zero simulates knockout. Multiple genes specify one joint perturbation.
#' @param group.by Metadata column defining the clusters used for GRN fitting.
#' @param base_grn Official TF-info data.frame with `peak_id`, `gene_short_name`
#' and TF columns, a path to its parquet file, or a named list mapping target
#' genes to character vectors of candidate TFs (the official `TFdict` format).
#' Required for a new fit; do not pass SCOP `TF/target/importance` edge tables.
#' No expression-derived GRN is generated automatically.
#' @param assay Assay containing raw counts. Default is `"RNA"`.
#' @param layer Non-log count layer. Default is `"counts"`. Fractional counts are
#' accepted, matching the official raw-count import (including the official demo
#' input). Log-transformed, scaled, integrated and SCT residual data should not
#' be supplied; their provenance cannot be inferred from numeric values alone.
#' @param features Genes to retain, in the requested order. `NULL` uses all genes
#' in the assay. CellOracle recommends selecting approximately 1,000--3,000 genes
#' before import; this wrapper does not silently select or normalize genes.
#' @param reduction Reduction containing the two-dimensional embedding. `NULL`
#' uses [DefaultReduction()]. The first two columns are used without recomputing.
#' @param stage_args Named list of named argument lists for official methods:
#' `pca` (`perform_PCA`), `imputation` (`knn_imputation`), `links` (`get_links`),
#' `filter` (`filter_links`), `fit` (`fit_GRN_for_simulation`), `simulation`
#' (`simulate_shift`), `transition` (`estimate_transition_prob`), `embedding`
#' (`calculate_embedding_shift`), `grid` (`calculate_p_mass`), and `mass_filter`
#' (`calculate_mass_filter`). Unspecified arguments retain upstream defaults.
#' SCOP controls the perturbation, cluster-specific fitting, randomization
#' reference, worker arguments and seed. In CellOracle 0.20.0,
#' `calculate_p_mass` internally fixes `smooth=0.8` and `n_jobs=-1`; the wrapper
#' preserves that upstream behavior. Do not supply these two ineffective grid
#' arguments. For small inputs explicitly choose suitable `k`, `n_pca_dims`,
#' transition `n_neighbors` and grid `n_neighbors`; no silent clipping occurs.
#' @param model_dir Directory containing `model.celloracle.oracle` and
#' `network.celloracle.links` saved by a previous run. Reuse requires identical
#' cells, genes, raw counts, cluster labels and embedding. Supply `base_grn=NULL`
#' and no fitting-stage arguments when reusing a model.
#' @param output_dir Optional new directory for official Oracle and Links HDF5
#' files. Existing directories are never overwritten. These files retain the
#' full model and allow subsequent perturbations without refitting.
#' @param store_expression Store cell-by-gene `imputed_count`, `simulated_count`
#' and `delta_X` matrices in the returned object. Default `FALSE` limits size.
#' @param tool_name Unused name for the result in `object@tools`. Existing entries
#' are rejected; use a different name for another perturbation.
#' @param envname Isolated environment prepared with
#' `PrepareEnv(modules = "celloracle")`. Default `"celloracle_env"`.
#' @param conda Conda-compatible executable, passed to [PrepareEnv()].
#' @param docker_container Optional name or ID of an already running Linux Docker
#' container whose `python` has CellOracle 0.20.0 installed. This allows Windows R
#' to use the same official runner without a local Python installation. SCOP
#' transfers inputs/results through a temporary directory; `envname` and `conda`
#' are unused. Model directories always refer to the R host. The function does
#' not pull images, install packages or start containers.
#' @param cores Positive worker count for official methods that honor it. The
#' official grid routine may use all cores; see `stage_args`.
#' @param seed Non-negative integer seed. NumPy and standard Python RNG states
#' are restored on exit; the upstream transition routine also seeds Numba.
#' Docker launches its child with `PYTHONHASHSEED=seed` to make upstream
#' set-derived predictor ordering reproducible across processes. Local Python
#' inherits its interpreter's hash seed; separate reference runs must use the
#' same hash seed to compare bagging coefficients exactly.
#' @param verbose Whether to display SCOP progress messages. Upstream diagnostic
#' output and warnings remain visible. Docker output is also retained in the
#' result's `diagnostics` character vector.
#'
#' @details
#' Requires official `celloracle==0.20.0`, either in a Linux/macOS Python
#' environment or in `docker_container`. Native Windows Python is not supported.
#' Prepare the environment or container explicitly before analysis.
#' The raw import performs CellOracle's own log transformation; no additional
#' library-size normalization is introduced. PCA, imputation, network inference,
#' edge filtering, coefficient fitting, simulation and projection all run in
#' CellOracle. Network edge p-values are not perturbation-effect p-values.
#' Results are predictions, not experimental knockout measurements.
#'
#' @return A Seurat object with `object@tools[[tool_name]]` containing `cell_shift`
#' (one row per input cell, with original coordinates, displacement and randomized
#' displacement), `grid_flow` (official grid vectors, mass and exclusion mask),
#' `gene_effect` (mean signed and mean absolute `delta_X` across cells), `links`
#' and `filtered_links` (upstream edge tables per cluster), `parameters`, `cells`,
#' `genes`, and optionally `expression` and `model_dir`. Expression matrices have
#' input cell and selected gene dimnames. No live Python object is stored.
#' @seealso [CellOraclePlot()], [PrepareEnv()]
#' @references Kamimoto et al. (2023). Dissecting cell identity via network
#' inference and in silico gene perturbation. Nature. \doi{10.1038/s41586-022-05688-9}.
#' @export
#' @examples
#' \dontrun{
#' # On Linux/macOS; prepare dependencies explicitly before analysis.
#' PrepareEnv(modules = "celloracle")
#' # On Windows use docker_container = "celloracle" with an already prepared,
#' # running Linux container instead of preparing local Python.
#' data(pancreas_sub)
#' pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#' # Supply an official TF-info parquet prepared for the matching species.
#' pancreas_sub <- RunCellOracle(
#'   pancreas_sub, perturbation = c(Pdx1 = 0), group.by = "CellType",
#'   base_grn = "mouse_promoter_base_GRN.parquet",
#'   stage_args = list(
#'     imputation = list(k = 20L, n_pca_dims = 20L),
#'     transition = list(n_neighbors = 50L),
#'     grid = list(n_neighbors = 50L)
#'   )
#' )
#' CellOraclePlot(pancreas_sub, plot_type = "quiver")
#' }
RunCellOracle <- function(
  object, perturbation, group.by, base_grn = NULL,
  assay = "RNA", layer = "counts", features = NULL, reduction = NULL,
  stage_args = list(), model_dir = NULL, output_dir = NULL,
  store_expression = FALSE, tool_name = "CellOracle",
  envname = "celloracle_env", conda = "auto", docker_container = NULL, cores = 1L,
  seed = 1234L, verbose = TRUE
) {
  if (!inherits(object, "Seurat")) {
    log_message("{.arg object} must be a Seurat object", message_type = "error")
  }
  if (!is.character(tool_name) || length(tool_name) != 1L ||
      is.na(tool_name) || !nzchar(tool_name) || tool_name %in% names(object@tools)) {
    log_message("{.arg tool_name} must be a new, non-empty result name", message_type = "error")
  }
  if (!is.numeric(perturbation) || !length(perturbation) ||
      any(!is.finite(perturbation)) || any(perturbation < 0) ||
      is.null(names(perturbation)) || anyNA(names(perturbation)) ||
      any(!nzchar(names(perturbation))) || anyDuplicated(names(perturbation))) {
    log_message("{.arg perturbation} must have unique gene names and finite non-negative values", message_type = "error")
  }
  if (length(group.by) != 1L || is.na(group.by) || !group.by %in% colnames(object[[]])) {
    log_message("{.arg group.by} must name one metadata column", message_type = "error")
  }
  groups <- object[[]][[group.by]]
  if (anyNA(groups) || any(!nzchar(as.character(groups)))) {
    log_message("Cluster labels must not be missing or empty", message_type = "error")
  }
  if (!assay %in% SeuratObject::Assays(object)) {
    log_message("Cannot find assay {.val {assay}}", message_type = "error")
  }
  counts <- GetAssayData5(object, assay = assay, layer = layer)
  features <- features %||% rownames(counts)
  if (!is.character(features) || length(features) < 2L || anyNA(features) ||
      anyDuplicated(features) || !all(features %in% rownames(counts)) ||
      !all(names(perturbation) %in% features)) {
    log_message("{.arg features} must contain unique assay genes including all perturbed genes", message_type = "error")
  }
  if (!identical(colnames(counts), colnames(object))) {
    log_message("Counts must include all object cells in the same order", message_type = "error")
  }
  counts <- counts[features, , drop = FALSE]
  values <- if (inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
  if (any(!is.finite(values)) || any(values < 0)) {
    log_message("CellOracle count input must contain finite non-negative values", message_type = "error")
  }
  if (any(Matrix::colSums(counts) == 0) || any(Matrix::rowSums(counts) == 0)) {
    log_message("Remove zero-count cells and genes before CellOracle import", message_type = "error")
  }
  reduction <- reduction %||% DefaultReduction(object)
  if (length(reduction) != 1L || !reduction %in% SeuratObject::Reductions(object)) {
    log_message("Cannot find the requested reduction", message_type = "error")
  }
  embedding <- SeuratObject::Embeddings(object[[reduction]])
  if (ncol(embedding) < 2L || !identical(rownames(embedding), colnames(object)) ||
      any(!is.finite(embedding[, 1:2, drop = FALSE]))) {
    log_message("Reduction must contain finite coordinates for all cells in order", message_type = "error")
  }
  stages <- c("pca", "imputation", "links", "filter", "fit", "simulation",
              "transition", "embedding", "grid", "mass_filter")
  if (!is.list(stage_args) || (length(stage_args) &&
      (is.null(names(stage_args)) || anyNA(names(stage_args)) ||
       anyDuplicated(names(stage_args)) || !all(names(stage_args) %in% stages)))) {
    log_message("{.arg stage_args} must be named by documented CellOracle stages", message_type = "error")
  }
  reserved <- c("perturb_condition", "GRN_unit", "use_cluster_specific_TFdict",
                "calculate_randomized", "random_seed", "threads", "n_jobs",
                "plot", "cluster_name_for_GRN_unit", "test_mode", "use_randomized_GRN")
  for (stage in names(stage_args)) {
    args <- stage_args[[stage]]
    if (!is.list(args) || (length(args) &&
        (is.null(names(args)) || anyNA(names(args)) || any(!nzchar(names(args))) ||
         anyDuplicated(names(args)) || any(names(args) %in% reserved))) ||
        (stage == "grid" && "smooth" %in% names(args))) {
      log_message("Invalid or SCOP-controlled arguments in {.val {stage}}", message_type = "error")
    }
    integer_args <- switch(stage,
      pca = "n_components", imputation = c("k", "n_pca_dims", "b_sight", "b_maxl"),
      links = c("bagging_number", "verbose_level"), filter = "threshold_number",
      fit = "verbose_level", simulation = "n_propagation",
      transition = "n_neighbors", grid = c("n_grid", "n_neighbors"), character())
    for (name in intersect(names(args), integer_args)) {
      value <- args[[name]]
      if (is.null(value)) next
      if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
          value < 0 || value > .Machine$integer.max || value != floor(value)) {
        log_message("{.val {stage}} argument {.val {name}} must be an integer", message_type = "error")
      }
      stage_args[[stage]][[name]] <- as.integer(value)
    }
  }
  for (value in list(cores = cores, seed = seed)) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0 || value > .Machine$integer.max || value != floor(value)) {
      log_message("{.arg cores} and {.arg seed} must be non-negative integers", message_type = "error")
    }
  }
  if (cores < 1L) log_message("{.arg cores} must be positive", message_type = "error")
  if (!is.logical(store_expression) || length(store_expression) != 1L || is.na(store_expression)) {
    log_message("{.arg store_expression} must be TRUE or FALSE", message_type = "error")
  }
  if (is.null(model_dir)) {
    valid_grn <- is.data.frame(base_grn) && all(c("peak_id", "gene_short_name") %in% names(base_grn))
    if (is.character(base_grn) && length(base_grn) == 1L && !is.na(base_grn) && file.exists(base_grn)) {
      base_grn <- normalizePath(base_grn, winslash = "/", mustWork = TRUE)
      valid_grn <- TRUE
    } else if (is.list(base_grn) && !is.data.frame(base_grn) && length(base_grn)) {
      valid_grn <- !is.null(names(base_grn)) && !anyNA(names(base_grn)) &&
        all(nzchar(names(base_grn))) && !anyDuplicated(names(base_grn)) &&
        all(vapply(base_grn, function(x) is.character(x) && !anyNA(x) && all(nzchar(x)), logical(1)))
    }
    if (!valid_grn) {
      log_message("Supply official TF-info data/parquet or a target-to-TF named list as {.arg base_grn}", message_type = "error")
    }
  } else {
    if (!is.null(base_grn) || any(names(stage_args) %in% stages[1:5]) ||
        length(model_dir) != 1L || !dir.exists(model_dir) ||
        !all(file.exists(file.path(model_dir, c("model.celloracle.oracle", "network.celloracle.links"))))) {
      log_message("Model reuse requires saved Oracle/Links files, no {.arg base_grn} and no fitting arguments", message_type = "error")
    }
    model_dir <- normalizePath(model_dir, winslash = "/", mustWork = TRUE)
  }
  if (!is.null(output_dir) && (length(output_dir) != 1L || is.na(output_dir) ||
      !nzchar(output_dir) || file.exists(output_dir))) {
    log_message("{.arg output_dir} must be a new directory", message_type = "error")
  }
  if (!is.null(docker_container)) {
    result <- celloracle_run_docker(
      counts = counts, embedding = embedding[, 1:2, drop = FALSE], groups = groups,
      group.by = group.by, reduction = reduction, base_grn = base_grn,
      perturbation = perturbation, stage_args = stage_args,
      model_dir = model_dir, output_dir = output_dir,
      store_expression = store_expression, cores = cores, seed = seed,
      container = docker_container, verbose = verbose)
  } else {
    if (is_windows()) {
      log_message("Official CellOracle requires Linux/macOS. Run this analysis in a supported environment; native Windows is unsupported", message_type = "error")
    }
    conda <- resolve_conda(conda)
    if (is.null(conda) || !env_exist(conda = conda, envname = envname)) {
      log_message("Prepare the CellOracle environment explicitly with {.code PrepareEnv(modules = 'celloracle')} before analysis", message_type = "error")
    }
    requirements <- env_requirements(modules = "celloracle")
    ready <- exist_python_pkgs(requirements$packages, envname = envname, conda = conda, verbose = FALSE)
    if (!all(ready)) {
      log_message("CellOracle environment dependencies are missing or incompatible; run {.code PrepareEnv(modules = 'celloracle')} explicitly", message_type = "error")
    }
    PrepareEnv(envname = envname, conda = conda, modules = "celloracle", verbose = verbose)
    check_python("celloracle==0.20.0", envname = envname, conda = conda, verbose = verbose)
    adata <- srt_to_adata(object, features = features, assay_x = assay,
      layer_x = layer, assay_y = character(), reductions = reduction,
      graphs = character(), neighbors = character(), prepare_env = FALSE, verbose = verbose)
    # Use a private copy and explicit ordered coordinates, leaving Seurat unchanged.
    pd <- reticulate::import("pandas", convert = FALSE)
    reticulate::py_set_item(reticulate::py_get_attr(adata, "obs"), group.by,
      pd$Categorical(as.character(groups),
        categories = if (is.factor(groups)) levels(droplevels(groups)) else sort(unique(as.character(groups))),
        ordered = is.ordered(groups)))
    reticulate::py_set_item(reticulate::py_get_attr(adata, "obsm"), paste0("X_", reduction),
      reticulate::r_to_py(embedding[, 1:2, drop = FALSE]))
    runner <- scop_python_import("celloracle_runner", convert = TRUE)
    log_message("Running official CellOracle 0.20.0", verbose = verbose)
    result <- runner$run_celloracle(
      adata = adata, group_by = group.by, embedding_name = paste0("X_", reduction),
      base_grn = base_grn, perturbation = as.list(perturbation),
      stage_args = stage_args, cores = as.integer(cores), seed = as.integer(seed),
      model_dir = model_dir, output_dir = output_dir, store_expression = store_expression)
  }
  result <- celloracle_validate_result(result, colnames(object), features)
  result$parameters$input <- list(assay = assay, layer = layer, features = features,
    group.by = group.by, reduction = reduction, perturbation = perturbation,
    docker_container = docker_container,
    stage_args = stage_args,
    base_grn_source = if (!is.null(model_dir)) "saved_model" else if (is.character(base_grn)) base_grn else if (is.data.frame(base_grn)) "TF_info_matrix" else "TFdict",
    base_grn_md5 = if (is.character(base_grn)) unname(tools::md5sum(base_grn)) else NULL)
  object@tools[[tool_name]] <- result
  log_message("CellOracle results stored in {.code object@tools[[tool_name]]}", verbose = verbose)
  object
}

celloracle_run_docker <- function(
  counts, embedding, groups, group.by, reduction, base_grn, perturbation,
  stage_args, model_dir, output_dir, store_expression, cores, seed, container, verbose
) {
  if (!is.character(container) || length(container) != 1L || is.na(container) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", container)) {
    log_message("{.arg docker_container} must be one container name or ID", message_type = "error")
  }
  docker <- unname(Sys.which("docker"))
  if (!nzchar(docker)) log_message("Docker CLI is not available", message_type = "error")
  directory <- tempfile("scop-celloracle-")
  dir.create(directory)
  directory <- normalizePath(directory, winslash = "/", mustWork = TRUE)
  if (!identical(dirname(directory), normalizePath(tempdir(), winslash = "/", mustWork = TRUE))) {
    log_message("Unexpected CellOracle temporary directory", message_type = "error")
  }
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  remote <- paste0("/tmp/", basename(directory))
  logfile <- file.path(directory, "docker.log")
  run <- function(args) {
    status <- system2(docker, shQuote(args), stdout = logfile, stderr = logfile)
    if (status != 0L) {
      detail <- paste(readLines(logfile, warn = FALSE), collapse = "\n")
      log_message("CellOracle Docker execution failed: {detail}", message_type = "error")
    }
  }
  prior_type <- if (is.data.frame(base_grn)) "matrix" else if (is.character(base_grn)) "parquet" else "dict"
  Matrix::writeMM(Matrix::Matrix(counts, sparse = TRUE), file.path(directory, "counts.mtx"))
  utils::write.table(matrix(sprintf("%.17g", embedding), ncol = 2L),
    file.path(directory, "embedding.csv"), sep = ",", quote = FALSE,
    row.names = FALSE, col.names = FALSE)
  if (prior_type == "matrix") {
    utils::write.csv(base_grn, file.path(directory, "prior.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  } else if (prior_type == "parquet") {
    file.copy(base_grn, file.path(directory, "prior.parquet"))
  }
  if (!is.null(model_dir)) {
    dir.create(file.path(directory, "model"))
    files <- c("model.celloracle.oracle", "network.celloracle.links")
    if (!all(file.copy(file.path(model_dir, files), file.path(directory, "model")))) {
      log_message("Unable to copy saved CellOracle model files", message_type = "error")
    }
  }
  runner_write_json(list(
    cells = as.list(colnames(counts)), genes = as.list(rownames(counts)),
    groups = as.list(as.character(groups)),
    levels = as.list(if (is.factor(groups)) levels(droplevels(groups)) else sort(unique(as.character(groups)))),
    ordered = is.ordered(groups), group_by = group.by,
    embedding_name = paste0("X_", reduction), prior_type = prior_type,
    prior = if (prior_type == "dict") base_grn else NULL,
    perturbation = as.list(perturbation), stage_args = stage_args,
    cores = as.integer(cores), seed = as.integer(seed),
    reuse_model = !is.null(model_dir), save_model = !is.null(output_dir),
    store_expression = store_expression
  ), file.path(directory, "request.json"))
  file.copy(runner_script_path("celloracle_runner.py", "CellOracle"), directory)
  run(c("exec", container, "mkdir", remote))
  on.exit(system2(docker, shQuote(c("exec", container, "rm", "-rf", "--", remote)),
    stdout = FALSE, stderr = FALSE), add = TRUE)
  run(c("cp", paste0(directory, "/."), paste0(container, ":", remote)))
  log_message("Running official CellOracle in container {.val {container}}", verbose = verbose)
  run(c("exec", "-e", paste0("PYTHONHASHSEED=", seed), container,
    "python", paste0(remote, "/celloracle_runner.py"), remote))
  diagnostics <- readLines(logfile, warn = FALSE)
  if (length(diagnostics)) {
    diagnostic_text <- paste(diagnostics, collapse = "\n")
    log_message("{diagnostic_text}", verbose = verbose)
  }
  run(c("cp", paste0(container, ":", remote, "/output"), directory))
  output <- file.path(directory, "output")
  result <- runner_read_json(file.path(output, "manifest.json"))
  result$diagnostics <- diagnostics
  result$cells <- unlist(result$cells, use.names = FALSE)
  result$genes <- unlist(result$genes, use.names = FALSE)
  read_table <- function(filename, character_columns = character()) {
    utils::read.csv(file.path(output, filename), check.names = FALSE, fileEncoding = "UTF-8",
      stringsAsFactors = FALSE, na.strings = "",
      colClasses = if (length(character_columns)) stats::setNames(rep("character", length(character_columns)), character_columns) else NA)
  }
  result$cell_shift <- read_table("cell_shift.csv", c("cell", "group"))
  result$grid_flow <- read_table("grid_flow.csv")
  result$grid_flow$filtered <- as.logical(result$grid_flow$filtered)
  result$gene_effect <- read_table("gene_effect.csv", "gene")
  for (name in c("links", "filtered_links")) {
    result[[name]] <- lapply(result[[name]], read_table, character_columns = c("source", "target"))
  }
  result$expression <- lapply(result$expression, function(filename) {
    x <- as.matrix(utils::read.csv(file.path(output, filename), header = FALSE))
    storage.mode(x) <- "double"
    x
  })
  if (!length(result$expression)) result$expression <- NULL
  result$parameters["model_dir"] <- list(model_dir)
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (!all(file.copy(list.files(file.path(output, "model"), full.names = TRUE), output_dir))) {
      log_message("Unable to save CellOracle model files on the R host", message_type = "error")
    }
    result$model_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  }
  result
}

celloracle_validate_result <- function(result, cells, genes) {
  if (!identical(as.character(result$cells), cells) ||
      !identical(as.character(result$genes), genes) ||
      !identical(as.character(result$cell_shift$cell), cells) ||
      !identical(as.character(result$gene_effect$gene), genes)) {
    log_message("CellOracle returned a different cell or gene order; results were not stored", message_type = "error")
  }
  numeric_columns <- list(
    cell_shift = c("x", "y", "dx", "dy", "dx_random", "dy_random"),
    grid_flow = c("x", "y", "dx", "dy", "dx_random", "dy_random", "mass"),
    gene_effect = c("mean_delta", "mean_abs_delta")
  )
  for (name in names(numeric_columns)) {
    data <- result[[name]]
    columns <- numeric_columns[[name]]
    if (!is.data.frame(data) || !all(columns %in% names(data)) ||
        any(!is.finite(as.matrix(data[, columns, drop = FALSE])))) {
      log_message("CellOracle returned invalid numerical results in {.val {name}}", message_type = "error")
    }
  }
  if (!is.logical(result$grid_flow$filtered) || anyNA(result$grid_flow$filtered)) {
    log_message("CellOracle returned an invalid grid exclusion mask", message_type = "error")
  }
  for (key in names(result$expression)) {
    x <- result$expression[[key]]
    if (!identical(dim(x), c(length(cells), length(genes)))) {
      log_message("CellOracle expression dimensions do not match input", message_type = "error")
    }
    dimnames(x) <- list(cells, genes)
    result$expression[[key]] <- x
  }
  result
}
