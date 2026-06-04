#' @title Run SCENIC gene regulatory network analysis
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams PrepareEnv
#' @param srt A Seurat object.
#' @param layer Assay layer used as the count matrix.
#' @param ranking_dbs Character vector of cisTarget ranking feather files. If
#' `NULL`, the gene-based v10 cisTarget ranking databases are prepared from
#' `species`.
#' @param motif_annotations Motif annotation table used by `scenic ctx`. If
#' `NULL`, the v10 motif2tf table is prepared from `species`.
#' @param regulators Transcription factors used as candidate regulators in
#' GRNBoost2. This can be a character vector of gene names or one text file. If
#' `NULL`, the cisTarget TF list is prepared from `species`.
#' @param targets Optional target genes used to restrict the GRN output. This
#' can be a character vector of gene names or one text file. Regulator
#' expression is still kept as predictor input for GRNBoost2, and the adjacency
#' table passed to `scenic ctx` is filtered to these targets.
#' @param work_dir Directory used for SCENIC input and output files.
#' @param species Species used to select cisTarget reference files when
#' `ranking_dbs`, `motif_annotations`, or `regulators` is `NULL`. Supported
#' values are `"Homo_sapiens"`, `"Mus_musculus"`, and
#' `"Drosophila_melanogaster"`.
#' @param genome Genome build used to select cisTarget reference files when
#' automatic references are prepared. Human supports `"hg38"` (default) and
#' `"hg19"`. Mouse and fly currently use `"mm10"` and `"dm6"`, respectively.
#' @param data_dir Directory used to cache automatically prepared SCENIC
#' reference files. If `NULL`, files are stored under
#' `tools::R_user_dir("scop", "data")/SCENIC/<species>`.
#' @param prefix Prefix for SCENIC output files.
#' @param group.by Optional metadata column used to aggregate single-cell
#' counts before GRNBoost2. When `NULL` (default), GRNBoost2 runs on the
#' original single-cell count matrix. To use pre-built metacells, pass the
#' metacell membership column (e.g. `"Metacell_id"` from [RunMetaCell()]).
#' All cells sharing the same `group.by` value are summed into one metacell
#' profile.
#' @param min_expr_cells Minimum number of cells or metacells where a gene must
#' be detected before GRNBoost2.
#' @param min_regulon_size Minimum regulon size kept after `scenic ctx`.
#' @param backend SCENIC backend. `"cpp"` uses the native R/C++ path and
#' `"python"` uses the Python `scenicplus` path.
#' @param grn_method GRN inference method for the native backend.
#' @param cistarget_method cisTarget implementation for the native backend.
#' @param max_regulon_targets Maximum number of target genes kept per regulon in
#' the native backend.
#' @param n_rounds Number of boosting rounds used by the native GRN backend.
#' @param learning_rate Learning rate used by the native GRN backend.
#' @param max_depth Maximum tree depth used by the native GRN backend.
#' @param max_features Fraction of features sampled by the native GRN backend.
#' @param subsample Row subsampling fraction used by the native GRN backend.
#' @param early_stop_window_length Early-stopping window used by the native GRN
#' backend.
#' @param cores Number of workers used by GRNBoost2, `scenic ctx`, and
#' AUCell batch scoring. If multicore execution is not supported, this is
#' automatically reduced to one core.
#' @param aucell_batch_size Number of cells scored in each AUCell batch.
#' @param aucell_backend Backend used for AUCell regulon activity scoring.
#' `"r"` uses the AUCell package. `"cpp"` uses the package C++ gene-set scoring
#' implementation.
#' @param aucell_cpp_strategy C++ AUCell ranking strategy passed to the package
#' gene-set scoring backend.
#' @param seed Random seed used by GRNBoost2 and Seurat overclustering.
#' @param force Whether to rebuild existing SCENIC outputs.
#' @param assay_name Name of the assay used to store regulon activity scores.
#' @param tool_name Name of the `srt@tools` entry.
#' @param return_seurat Whether to return the modified Seurat object. If
#' `FALSE`, a result list is returned.
#' @param envname Python environment used for SCENIC. If `NULL`, the isolated
#' `"scenic_env"` environment is used.
#' @param conda The path or command name of a conda-compatible executable.
#' @param prepare_env Whether to prepare and configure the SCENIC Python
#' environment before running.
#'
#' @return A Seurat object with SCENIC results, or a result list when
#' `return_seurat = FALSE`.
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunSCENIC(
#'   pancreas_sub,
#'   species = "Mus_musculus"
#' )
#' }
RunSCENIC <- function(
  srt,
  assay = NULL,
  layer = "counts",
  ranking_dbs = NULL,
  motif_annotations = NULL,
  regulators = NULL,
  targets = NULL,
  work_dir = "scenic_output/",
  species = c("Homo_sapiens", "Mus_musculus", "Drosophila_melanogaster"),
  genome = NULL,
  data_dir = NULL,
  prefix = "scenic",
  group.by = NULL,
  min_expr_cells = 3,
  min_regulon_size = 10,
  backend = c("cpp", "python"),
  grn_method = c("grnboost2", "genie3"),
  cistarget_method = c("native_motif", "native_approx"),
  max_regulon_targets = 50,
  n_rounds = 5000,
  learning_rate = 0.01,
  max_depth = 3,
  max_features = 0.1,
  subsample = 0.9,
  early_stop_window_length = 25,
  cores = 1,
  aucell_batch_size = 500,
  aucell_backend = c("r", "cpp"),
  aucell_cpp_strategy = c("full", "sparse", "topk"),
  seed = 1234,
  force = FALSE,
  assay_name = "scenic",
  tool_name = "SCENIC",
  return_seurat = TRUE,
  envname = NULL,
  conda = "auto",
  prepare_env = TRUE,
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (missing(work_dir) || length(work_dir) != 1) {
    log_message("{.arg work_dir} must be one directory", message_type = "error")
  }
  species <- match.arg(species)
  backend <- match.arg(backend)
  grn_method <- match.arg(grn_method)
  cistarget_method <- match.arg(cistarget_method)
  aucell_backend <- match.arg(aucell_backend)
  aucell_cpp_strategy <- match.arg(aucell_cpp_strategy)

  if (identical(backend, "cpp")) {
    return(scenic_cpp(
      srt = srt,
      assay = assay,
      layer = layer,
      regulators = regulators,
      targets = targets,
      work_dir = work_dir,
      species = species,
      genome = genome,
      data_dir = data_dir,
      prefix = prefix,
      group.by = group.by,
      min_expr_cells = min_expr_cells,
      min_regulon_size = min_regulon_size,
      grn_method = grn_method,
      cistarget_method = cistarget_method,
      max_regulon_targets = max_regulon_targets,
      ranking_dbs = ranking_dbs,
      motif_annotations = motif_annotations,
      n_rounds = n_rounds,
      learning_rate = learning_rate,
      max_depth = max_depth,
      max_features = max_features,
      subsample = subsample,
      early_stop_window_length = early_stop_window_length,
      cores = cores,
      aucell_batch_size = aucell_batch_size,
      aucell_backend = aucell_backend,
      aucell_cpp_strategy = aucell_cpp_strategy,
      seed = seed,
      force = force,
      assay_name = assay_name,
      tool_name = tool_name,
      return_seurat = return_seurat,
      verbose = verbose
    ))
  }

  # Python path below.

  reference_data <- scenic_reference(
    species = species,
    genome = genome,
    data_dir = data_dir,
    ranking_dbs = ranking_dbs,
    motif_annotations = motif_annotations,
    regulators = regulators,
    verbose = verbose
  )
  ranking_dbs <- reference_data[["ranking_dbs"]]
  motif_annotations <- reference_data[["motif_annotations"]]
  regulators <- reference_data[["regulators"]]
  species <- reference_data[["species"]]
  genome <- reference_data[["genome"]]
  data_dir <- reference_data[["data_dir"]]

  if (length(ranking_dbs) == 0) {
    log_message("{.arg ranking_dbs} must be provided", message_type = "error")
  }
  if (length(motif_annotations) != 1) {
    log_message(
      "{.arg motif_annotations} must be one file",
      message_type = "error"
    )
  }
  if (is.null(regulators) || length(regulators) == 0) {
    log_message(
      "{.arg regulators} must contain at least one gene",
      message_type = "error"
    )
  }

  ranking_dbs <- normalizePath(ranking_dbs, mustWork = TRUE)
  motif_annotations <- normalizePath(motif_annotations, mustWork = TRUE)
  work_dir <- normalizePath(work_dir, mustWork = FALSE)
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  regulators_info <- scenic_prep_gene_arg(
    x = regulators,
    arg = "regulators",
    out_file = file.path(work_dir, paste0(prefix, "_regulators.txt")),
    required = TRUE
  )
  regulators <- regulators_info[["genes"]]
  regulators_file <- regulators_info[["file"]]
  targets <- scenic_prep_gene_arg(
    x = targets,
    arg = "targets",
    out_file = file.path(work_dir, paste0(prefix, "_targets.txt")),
    required = FALSE
  )[["genes"]]
  progress_state <- scenic_progress_init(verbose = verbose)
  on.exit(scenic_progress_close(progress_state), add = TRUE)
  scenic_progress_step(
    progress_state,
    value = 5,
    label = "Validated SCENIC inputs",
    verbose = verbose
  )

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  cores_requested <- suppressWarnings(as.integer(cores))
  if (
    length(cores_requested) != 1 ||
      is.na(cores_requested) ||
      cores_requested < 1L
  ) {
    cores_requested <- 1L
  }
  cores <- max(1L, cores_requested)
  detected_cores <- tryCatch(
    parallel::detectCores(logical = TRUE),
    error = function(...) NA_integer_
  )
  if (cores_requested > 1L && (is.na(detected_cores) || detected_cores < 2L)) {
    log_message(
      "{.arg cores} was set to {.val {cores_requested}}, but multicore support was not detected; using one core",
      message_type = "warning",
      verbose = verbose
    )
    cores <- 1L
  }
  if (cores_requested > 1L && .Platform$OS.type != "unix") {
    log_message(
      "{.arg cores} was set to {.val {cores_requested}}, but multicore AUCell scoring requires a Unix-like OS; using one core",
      message_type = "warning",
      verbose = verbose
    )
    cores <- 1L
  }
  aucell_batch_size <- max(1L, as.integer(aucell_batch_size))

  envname <- envname %||% "scenic_env"
  scenic_python_packages <- c(
    scenic_backend_requirement(),
    "arboreto==0.1.6",
    "ctxcore==0.2.0",
    "numpy==1.23.5",
    "dask==2024.2.1",
    "distributed==2024.2.1",
    "pyarrow"
  )
  if (isTRUE(prepare_env)) {
    scenic_progress_step(
      progress_state,
      value = 10,
      label = "Checking SCENIC Python environment",
      verbose = verbose
    )
    env_ready <- FALSE
    env_exists <- isTRUE(tryCatch(
      env_exist(conda = conda, envname = envname),
      error = function(...) FALSE
    ))
    if (env_exists) {
      pkg_installed <- tryCatch(
        exist_python_pkgs(
          packages = scenic_python_packages,
          envname = envname,
          conda = conda,
          verbose = FALSE
        ),
        error = function(...) {
          stats::setNames(
            rep(FALSE, length(scenic_python_packages)),
            scenic_python_packages
          )
        }
      )
      env_ready <- all(pkg_installed)
    }
    if (isTRUE(env_ready)) {
      log_message(
        "{.pkg SCENIC} Python environment already has the required packages; skipping {.fn PrepareEnv}",
        verbose = verbose
      )
    } else {
      if (isTRUE(env_exists)) {
        missing_packages <- names(pkg_installed)[!pkg_installed]
        log_message(
          "{.pkg SCENIC} Python environment is missing {.val {length(missing_packages)}} package{?s}; preparing environment",
          verbose = verbose
        )
      }
      scenic_progress_step(
        progress_state,
        value = 10,
        label = "Preparing SCENIC Python environment",
        verbose = verbose
      )
      PrepareEnv(
        envname = envname,
        conda = conda,
        version = "3.10-1",
        modules = "scenic"
      )
    }
  }
  scenic_progress_step(
    progress_state,
    value = 15,
    label = "Checking SCENIC Python packages",
    verbose = verbose
  )
  check_python(
    packages = scenic_python_packages,
    envname = envname,
    conda = conda,
    verbose = verbose
  )
  conda_resolved <- resolve_conda(conda)
  scenic_python <- conda_python(
    conda = conda_resolved,
    envname = envname
  )
  assert_python_runtime_switchable(
    scenic_python,
    restart_hint = scenic_runtime_restart_hint(envname = envname)
  )
  configure_python_runtime(scenic_python)

  functions <- reticulate::import_from_path(
    "functions",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )

  log_message(
    "Preparing {.pkg SCENIC} input matrix",
    verbose = verbose
  )
  scenic_progress_step(
    progress_state,
    value = 20,
    label = "Loading expression counts",
    verbose = verbose
  )
  counts <- GetAssayData5(srt, assay = assay, layer = layer)
  counts <- counts[, colnames(srt), drop = FALSE]
  if (nrow(counts) == 0 || ncol(counts) == 0) {
    log_message(
      "No expression values available for {.pkg SCENIC}",
      message_type = "error"
    )
  }
  log_message(
    "{.pkg SCENIC} expression input: assay {.val {assay}}, layer {.val {layer}}, {.val {nrow(counts)}} genes x {.val {ncol(counts)}} cells",
    verbose = verbose
  )

  if (!is.null(group.by)) {
    scenic_progress_step(
      progress_state,
      value = 30,
      label = "Aggregating counts by group.by for GRNBoost2",
      verbose = verbose
    )
    group_col <- srt[[group.by]][, 1]
    if (any(is.na(group_col))) {
      log_message(
        "{.arg group.by} column {.val {group.by}} contains NA values",
        message_type = "error"
      )
    }
    cell_groups <- split(colnames(srt), group_col)
    grn_counts <- do.call(
      cbind,
      lapply(cell_groups, function(cells) {
        Matrix::rowSums(counts[, cells, drop = FALSE])
      })
    )
    grn_counts <- Matrix::Matrix(grn_counts, sparse = TRUE)
    colnames(grn_counts) <- names(cell_groups)
    rownames(grn_counts) <- rownames(counts)
    grn_count_source <- "metacell matrix"
    grn_count_unit <- "groups"
  } else {
    grn_counts <- counts
    grn_count_source <- "single-cell matrix"
    grn_count_unit <- "cells"
    scenic_progress_step(
      progress_state,
      value = 30,
      label = "Using single-cell counts for GRNBoost2",
      verbose = verbose
    )
  }
  log_message(
    "{.pkg SCENIC} GRN count source: {.val {grn_count_source}}, {.val {nrow(grn_counts)}} genes x {.val {ncol(grn_counts)}} {grn_count_unit}",
    verbose = verbose
  )

  scenic_progress_step(
    progress_state,
    value = 40,
    label = "Reading cisTarget ranking database genes",
    verbose = verbose
  )
  ranking_gene_lists <- lapply(ranking_dbs, function(db) {
    as.character(unlist(functions$SCENICRankingGenes(db), use.names = FALSE))
  })
  ranking_genes <- Reduce(intersect, ranking_gene_lists)
  grn_gene_filter <- NULL
  if (!is.null(targets)) {
    grn_gene_filter <- union(targets, regulators)
  }
  grn_matrix <- scenic_prepare_grn_matrix(
    counts = grn_counts,
    ranking_genes = ranking_genes,
    genes = grn_gene_filter,
    min_expr_cells = min_expr_cells,
    verbose = verbose
  )
  regulators_detected <- intersect(regulators, colnames(grn_matrix))
  if (length(regulators_detected) == 0) {
    log_message(
      "No {.arg regulators} remain after matching expression features and cisTarget databases",
      message_type = "error"
    )
  }
  if (length(regulators_detected) < length(regulators)) {
    log_message(
      "{.pkg SCENIC} GRNBoost2 will use {.val {length(regulators_detected)}} of {.val {length(regulators)}} requested regulator{?s} present in the GRN matrix",
      message_type = "warning",
      verbose = verbose
    )
    regulators <- regulators_detected
    regulators_file <- scenic_write_gene_list(regulators, regulators_file)
  }
  if (!is.null(targets)) {
    targets_detected <- intersect(targets, colnames(grn_matrix))
    if (length(targets_detected) == 0) {
      log_message(
        "No {.arg targets} remain after matching expression features and cisTarget databases",
        message_type = "error"
      )
    }
    if (length(targets_detected) < length(targets)) {
      log_message(
        "{.pkg SCENIC} target constraint keeps {.val {length(targets_detected)}} of {.val {length(targets)}} requested target{?s} present in the GRN matrix",
        message_type = "warning",
        verbose = verbose
      )
    }
    targets <- targets_detected
  }

  scenic_progress_step(
    progress_state,
    value = 48,
    label = "Writing GRNBoost2 expression matrix",
    verbose = verbose
  )
  expr_csv <- file.path(work_dir, paste0(prefix, "_expression_mtx.csv"))
  grn_input_params_file <- file.path(
    work_dir,
    paste0(prefix, "_grn_input_params.rds")
  )
  grn_input_params <- list(
    regulators = regulators,
    targets = targets,
    genes = colnames(grn_matrix)
  )
  grn_inputs_changed <- scenic_grn_inputs_changed(
    grn_input_params_file,
    grn_input_params
  )
  grn_force <- isTRUE(force) || isTRUE(grn_inputs_changed)
  if (isTRUE(grn_inputs_changed) && file.exists(expr_csv) && isFALSE(force)) {
    log_message(
      "Rebuilding {.pkg SCENIC} GRNBoost2 inputs because {.arg regulators}, {.arg targets}, or GRN genes changed",
      verbose = verbose
    )
  }
  if (!file.exists(expr_csv) || isTRUE(grn_force)) {
    utils::write.csv(
      as.matrix(grn_matrix),
      file = expr_csv,
      row.names = TRUE
    )
    saveRDS(grn_input_params, grn_input_params_file)
  } else {
    log_message(
      "Reusing existing {.pkg SCENIC} expression matrix: {.file {expr_csv}}",
      verbose = verbose
    )
  }

  adj_file <- file.path(work_dir, paste0(prefix, "_step1_adj.tsv"))
  grn_adj_file <- if (is.null(targets)) {
    adj_file
  } else {
    file.path(work_dir, paste0(prefix, "_step1_adj_raw.tsv"))
  }
  ctx_file <- file.path(work_dir, paste0(prefix, "_step2_reg.tsv"))
  gmt_file <- file.path(work_dir, paste0(prefix, "_step2_regulons.gmt"))
  txt_file <- file.path(work_dir, paste0(prefix, "_step2_regulons.txt"))
  ras_file <- file.path(work_dir, paste0(prefix, "_regulon_activity_score.rds"))
  regulon_list_file <- file.path(work_dir, paste0(prefix, "_regulon_list.rds"))

  scenic_progress_step(
    progress_state,
    value = 55,
    label = "Running GRNBoost2",
    verbose = verbose
  )
  functions$RunSCENICGrn(
    expression_mtx = expr_csv,
    regulators = regulators_file,
    adj_output = grn_adj_file,
    cores = as.integer(cores),
    seed = as.integer(seed),
    force = isTRUE(grn_force),
    verbose = isTRUE(verbose)
  )
  if (!is.null(targets)) {
    scenic_flt_adj(
      input_file = grn_adj_file,
      output_file = adj_file,
      targets = targets,
      verbose = verbose
    )
  }

  scenic_progress_step(
    progress_state,
    value = 70,
    label = "Running SCENIC cisTarget pruning",
    verbose = verbose
  )
  functions$RunSCENICCtx(
    expression_mtx = expr_csv,
    ranking_dbs = as.list(ranking_dbs),
    motif_annotations = motif_annotations,
    adj_output = adj_file,
    ctx_output = ctx_file,
    cores = as.integer(cores),
    force = isTRUE(grn_force),
    verbose = isTRUE(verbose)
  )

  scenic_progress_step(
    progress_state,
    value = 82,
    label = "Converting SCENIC regulons",
    verbose = verbose
  )
  if (isTRUE(grn_force) || !file.exists(gmt_file) || !file.exists(txt_file)) {
    functions$SCENICRegulonsToFiles(
      regulon_file = ctx_file,
      gmt_file = gmt_file,
      txt_file = txt_file,
      min_regulon_size = as.integer(min_regulon_size)
    )
  } else {
    log_message(
      "Reusing existing SCENIC regulon files",
      verbose = verbose
    )
  }

  scenic_progress_step(
    progress_state,
    value = 88,
    label = "Reading regulon target lists",
    verbose = verbose
  )
  regulon_tbl <- scenic_read_regulon_txt(txt_file)
  regulon_list <- scenic_regulon_list(regulon_tbl)
  if (!file.exists(regulon_list_file) || isTRUE(grn_force)) {
    saveRDS(regulon_list, regulon_list_file)
  }

  if (file.exists(ras_file) && isFALSE(grn_force)) {
    log_message(
      "Reusing existing regulon activity matrix: {.file {ras_file}}",
      verbose = verbose
    )
    ras_mat <- readRDS(ras_file)
  } else {
    scenic_progress_step(
      progress_state,
      value = 92,
      label = "Calculating AUCell regulon activity scores",
      verbose = verbose
    )
    ras_mat <- scenic_compute_aucell_score(
      counts = counts,
      regulon_list = regulon_list,
      min_regulon_size = min_regulon_size,
      batch_size = aucell_batch_size,
      cores = cores,
      backend = aucell_backend,
      cpp_strategy = aucell_cpp_strategy,
      verbose = verbose
    )
    saveRDS(ras_mat, ras_file)
  }
  ras_mat <- ras_mat[colnames(srt), , drop = FALSE]

  scenic_progress_step(
    progress_state,
    value = 96,
    label = "Collecting SCENIC result files",
    verbose = verbose
  )
  adjacency <- if (file.exists(adj_file)) {
    utils::read.delim(adj_file, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    NULL
  }
  result <- list(
    scores = Matrix::t(Matrix::Matrix(as.matrix(ras_mat), sparse = TRUE)),
    scores_cells_by_regulon = ras_mat,
    regulons = regulon_tbl,
    regulon_list = regulon_list,
    adjacency = adjacency,
    files = list(
      expression_mtx = expr_csv,
      adj = adj_file,
      adj_raw = if (!identical(grn_adj_file, adj_file)) grn_adj_file else NULL,
      grn_input_params = grn_input_params_file,
      ctx = ctx_file,
      regulons_gmt = gmt_file,
      regulons_txt = txt_file,
      regulon_activity_score = ras_file,
      regulon_list = regulon_list_file
    ),
    group.by = group.by,
    parameters = list(
      assay = assay,
      layer = layer,
      species = species,
      genome = genome,
      data_dir = data_dir,
      ranking_dbs = ranking_dbs,
      motif_annotations = motif_annotations,
      regulators = regulators,
      regulators_file = regulators_file,
      targets = targets,
      group.by = group.by,
      min_expr_cells = min_expr_cells,
      min_regulon_size = min_regulon_size,
      cores = cores,
      aucell_batch_size = aucell_batch_size,
      aucell_backend = aucell_backend,
      aucell_cpp_strategy = aucell_cpp_strategy,
      seed = seed,
      assay_name = assay_name,
      tool_name = tool_name,
      envname = envname,
      progress = verbose
    )
  )

  if (isFALSE(return_seurat)) {
    scenic_progress_step(
      progress_state,
      value = 100,
      label = "Finished SCENIC",
      verbose = verbose
    )
    return(result)
  }

  srt[[assay_name]] <- Seurat::CreateAssayObject(data = result[["scores"]])
  srt[[assay_name]] <- Seurat::AddMetaData(
    object = srt[[assay_name]],
    metadata = data.frame(
      termnames = rownames(result[["scores"]]),
      target_count = lengths(regulon_list)[rownames(result[["scores"]])],
      row.names = rownames(result[["scores"]]),
      stringsAsFactors = FALSE
    )
  )
  srt@tools[[tool_name]] <- result
  scenic_progress_step(
    progress_state,
    value = 100,
    label = "Stored SCENIC results",
    verbose = verbose
  )
  log_message(
    "{.pkg SCENIC} results stored in assay {.val {assay_name}} and tools slot {.val {tool_name}}",
    verbose = verbose
  )
  srt
}

scenic_cpp <- function(
  srt,
  assay,
  layer,
  regulators,
  targets,
  work_dir,
  species,
  genome,
  data_dir,
  prefix,
  group.by,
  min_expr_cells,
  min_regulon_size,
  grn_method,
  cistarget_method,
  max_regulon_targets,
  ranking_dbs,
  motif_annotations,
  n_rounds,
  learning_rate,
  max_depth,
  max_features,
  subsample,
  early_stop_window_length,
  cores,
  aucell_batch_size,
  aucell_backend,
  aucell_cpp_strategy,
  seed,
  force,
  assay_name,
  tool_name,
  return_seurat,
  verbose
) {
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  species_config <- scenic_species_config(species, genome = genome)
  species <- species_config[["label"]]
  genome <- species_config[["genome"]]
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

  work_dir <- normalizePath(work_dir, mustWork = FALSE)
  ras_file <- file.path(
    work_dir,
    paste0(prefix, "_regulon_activity_score_cpp.rds")
  )
  adj_file <- file.path(work_dir, paste0(prefix, "_adj_cpp.tsv"))
  regulon_file <- file.path(work_dir, paste0(prefix, "_regulon_list_cpp.rds"))

  if (is.null(regulators) || length(regulators) == 0) {
    reference_data <- scenic_reference(
      species = species,
      genome = genome,
      data_dir = data_dir,
      verbose = verbose
    )
    regulators_raw <- reference_data[["regulators"]]
    if (is.null(regulators_raw) || length(regulators_raw) == 0) {
      log_message(
        "{.arg regulators} must be provided or resolvable from {.arg species}",
        message_type = "error"
      )
    }

    if (
      is.character(regulators_raw) &&
        length(regulators_raw) == 1 &&
        file.exists(regulators_raw)
    ) {
      regulators <- unique(readLines(regulators_raw, warn = FALSE))
      regulators <- regulators[nzchar(regulators) & !grepl("^#", regulators)]
    } else {
      regulators <- regulators_raw
    }
    if (length(regulators) == 0) {
      log_message(
        "No regulators found in reference data",
        message_type = "error"
      )
    }
  }

  counts <- GetAssayData5(srt, assay = assay, layer = layer)

  if (!is.null(group.by)) {
    if (!group.by %in% colnames(srt@meta.data)) {
      log_message(
        "{.arg group.by} column {.val {group.by}} not found",
        message_type = "error"
      )
    }
    groups <- srt@meta.data[[group.by]]
    grn_matrix <- do.call(
      cbind,
      lapply(split(seq_len(ncol(counts)), groups), function(idx) {
        Matrix::rowSums(counts[, idx, drop = FALSE])
      })
    )
    colnames(grn_matrix) <- levels(factor(groups))
  } else {
    grn_matrix <- counts
  }

  gene_expr_n <- Matrix::rowSums(grn_matrix > 0)
  keep_genes <- gene_expr_n >= min_expr_cells
  grn_matrix <- grn_matrix[keep_genes, , drop = FALSE]

  if (nrow(grn_matrix) == 0) {
    log_message(
      "No genes remain after filtering by expression",
      message_type = "error"
    )
  }

  regulators_detected <- intersect(regulators, rownames(grn_matrix))
  if (length(regulators_detected) == 0) {
    log_message(
      "No {.arg regulators} overlap with expression matrix",
      message_type = "error"
    )
  }
  regulators <- regulators_detected

  if (!is.null(targets)) {
    targets <- intersect(targets, rownames(grn_matrix))
    if (length(targets) == 0) {
      log_message(
        "No {.arg targets} overlap with expression matrix",
        message_type = "error"
      )
    }
  }

  log_message(
    "Running SCENIC ({.val {grn_method}}) on {.val {ncol(grn_matrix)}} cells with {.val {length(regulators)}} TFs",
    verbose = verbose
  )

  grn_force <- isTRUE(force) || !file.exists(adj_file)
  if (grn_force) {
    adjacency <- scenic_run_grn_method(
      grn_matrix = Matrix::t(grn_matrix),
      regulators = regulators,
      targets = targets,
      grn_method = grn_method,
      backend = "cpp",
      output_file = adj_file,
      max_edges_per_target = max_regulon_targets,
      n_rounds = n_rounds,
      learning_rate = learning_rate,
      max_depth = max_depth,
      max_features = max_features,
      subsample = subsample,
      early_stop_window_length = early_stop_window_length,
      seed = seed,
      cores = cores,
      force = TRUE,
      verbose = verbose
    )
  } else {
    log_message(
      "Reusing existing GRN output: {.file {adj_file}}",
      verbose = verbose
    )
    adjacency <- utils::read.delim(
      adj_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  if (nrow(adjacency) == 0) {
    log_message(
      "GRN inference produced no edges",
      message_type = "error"
    )
  }

  regulon_force <- grn_force || !file.exists(regulon_file)
  if (regulon_force) {
    if (identical(cistarget_method, "native_motif")) {
      ref_data <- tryCatch(
        scenic_reference(
          species = species,
          genome = genome,
          data_dir = data_dir,
          ranking_dbs = ranking_dbs,
          motif_annotations = motif_annotations,
          verbose = verbose
        ),
        error = function(e) NULL
      )
      if (
        !is.null(ref_data) &&
          length(ref_data[["ranking_dbs"]]) > 0 &&
          !is.null(ref_data[["motif_annotations"]])
      ) {
        regulon_list <- cistarget2(
          adjacency = adjacency,
          ranking_dbs = ref_data[["ranking_dbs"]],
          motif_annotations = ref_data[["motif_annotations"]],
          max_targets = max_regulon_targets,
          min_regulon_size = min_regulon_size,
          verbose = verbose
        )
      } else {
        log_message(
          "cisTarget reference data not available; falling back to native_approx",
          message_type = "warning",
          verbose = verbose
        )
        regulon_list <- build_regulons(
          adjacency = adjacency,
          max_targets = max_regulon_targets,
          min_regulon_size = min_regulon_size
        )
      }
    } else {
      regulon_list <- build_regulons(
        adjacency = adjacency,
        max_targets = max_regulon_targets,
        min_regulon_size = min_regulon_size
      )
    }
    saveRDS(regulon_list, regulon_file)
  } else {
    log_message(
      "Reusing existing regulon list: {.file {regulon_file}}",
      verbose = verbose
    )
    regulon_list <- readRDS(regulon_file)
  }

  if (length(regulon_list) == 0) {
    log_message(
      "No regulons remain after filtering",
      message_type = "error"
    )
  }

  log_message(
    "Native regulon builder produced {.val {length(regulon_list)}} regulons",
    verbose = verbose
  )

  if (file.exists(ras_file) && isFALSE(grn_force)) {
    log_message(
      "Reusing existing AUCell scores: {.file {ras_file}}",
      verbose = verbose
    )
    ras_mat <- readRDS(ras_file)
  } else {
    ras_mat <- scenic_compute_aucell_score(
      counts = counts,
      regulon_list = regulon_list,
      min_regulon_size = min_regulon_size,
      batch_size = aucell_batch_size,
      cores = cores,
      backend = aucell_backend,
      cpp_strategy = aucell_cpp_strategy,
      verbose = verbose
    )
    saveRDS(ras_mat, ras_file)
  }
  ras_mat <- ras_mat[colnames(srt), , drop = FALSE]

  regulon_tbl <- do.call(
    rbind,
    lapply(names(regulon_list), function(tf) {
      data.frame(
        regulon = tf,
        target = regulon_list[[tf]],
        target_count = length(regulon_list[[tf]]),
        stringsAsFactors = FALSE
      )
    })
  )

  result <- list(
    scores = Matrix::t(Matrix::Matrix(as.matrix(ras_mat), sparse = TRUE)),
    scores_cells_by_regulon = ras_mat,
    regulons = regulon_tbl,
    regulon_list = regulon_list,
    adjacency = adjacency,
    files = list(
      adj = adj_file,
      regulon_list = regulon_file,
      regulon_activity_score = ras_file
    ),
    group.by = group.by,
    parameters = list(
      backend = "cpp",
      grn_method = grn_method,
      cistarget_method = cistarget_method,
      method = if (identical(cistarget_method, "native_motif")) {
        "native_cistarget"
      } else {
        "native_approx"
      },
      assay = assay,
      layer = layer,
      species = species,
      genome = genome,
      regulators = regulators,
      targets = targets,
      group.by = group.by,
      min_expr_cells = min_expr_cells,
      min_regulon_size = min_regulon_size,
      max_regulon_targets = max_regulon_targets,
      cores = cores,
      aucell_batch_size = aucell_batch_size,
      aucell_backend = aucell_backend,
      aucell_cpp_strategy = aucell_cpp_strategy,
      seed = seed,
      assay_name = assay_name,
      tool_name = tool_name
    )
  )

  if (isFALSE(return_seurat)) {
    return(result)
  }

  srt[[assay_name]] <- Seurat::CreateAssayObject(data = result[["scores"]])
  srt[[assay_name]] <- Seurat::AddMetaData(
    object = srt[[assay_name]],
    metadata = data.frame(
      termnames = rownames(result[["scores"]]),
      target_count = lengths(regulon_list)[rownames(result[["scores"]])],
      row.names = rownames(result[["scores"]]),
      stringsAsFactors = FALSE
    )
  )
  srt@tools[[tool_name]] <- result
  log_message(
    "{.pkg SCENIC} results stored in assay {.val {assay_name}} and tools slot {.val {tool_name}}",
    message_type = "success",
    verbose = verbose
  )
  srt
}

build_regulons <- function(
  adjacency,
  max_targets = 50,
  min_regulon_size = 10
) {
  if (!all(c("TF", "target", "importance") %in% colnames(adjacency))) {
    log_message(
      "{.arg adjacency} must have columns TF, target, importance",
      message_type = "error"
    )
  }

  adj_sorted <- adjacency[
    order(
      adjacency[["TF"]],
      -as.numeric(adjacency[["importance"]])
    ),
  ]

  adj_split <- split(adj_sorted, adj_sorted[["TF"]])
  regulons <- lapply(adj_split, function(df) {
    targets <- unique(df[["target"]])
    utils::head(targets, max_targets)
  })

  regulons <- regulons[lengths(regulons) >= min_regulon_size]
  regulons
}

cistarget2 <- function(
  adjacency,
  ranking_dbs,
  motif_annotations,
  max_targets = 50,
  min_regulon_size = 10,
  nes_threshold = 1.5,
  auc_threshold = 0.001,
  cores = 1,
  verbose = TRUE
) {
  check_r("arrow", verbose = FALSE)

  motif_tbl <- tryCatch(
    utils::read.table(
      motif_annotations,
      header = FALSE,
      sep = "\t",
      stringsAsFactors = FALSE,
      comment.char = ""
    ),
    error = function(e) {
      log_message(
        "Failed to read motif annotations: {.val {conditionMessage(e)}}",
        message_type = "error"
      )
    }
  )
  colnames(motif_tbl) <- c(
    "motif",
    "tf",
    "direct_annotation",
    "inferred_annotation"
  )

  motif_to_tf <- split(motif_tbl[["tf"]], motif_tbl[["motif"]])
  motif_names <- names(motif_to_tf)

  log_message(
    "Native cisTarget: loaded {.val {length(motif_names)}} motifs from annotations",
    verbose = verbose
  )

  rank_matrices <- list()
  all_genes <- character(0)
  cluster_names_all <- character(0)

  for (db_path in ranking_dbs) {
    log_message(
      "  Reading ranking database: {.file {basename(db_path)}}",
      verbose = verbose
    )
    db <- tryCatch(
      arrow::read_feather(db_path),
      error = function(e) {
        log_message(
          "Failed to read {.file {db_path}}: {.val {conditionMessage(e)}}",
          message_type = "warning",
          verbose = verbose
        )
        NULL
      }
    )
    if (is.null(db)) {
      next
    }

    if (!"motifs" %in% colnames(db)) {
      log_message(
        "Ranking database missing 'motifs' column. Skipping.",
        message_type = "warning",
        verbose = verbose
      )
      next
    }
    clusters <- as.character(db[["motifs"]])
    gene_cols <- setdiff(colnames(db), "motifs")
    all_genes <- union(all_genes, gene_cols)
    cluster_names_all <- union(cluster_names_all, clusters)

    rank_matrices[[length(rank_matrices) + 1]] <- list(
      clusters = clusters,
      genes = gene_cols,
      data = db
    )
  }

  if (length(rank_matrices) == 0) {
    log_message(
      "Failed to read any ranking databases. Falling back to native_approx.",
      message_type = "warning",
      verbose = verbose
    )
    return(build_regulons(
      adjacency,
      max_targets,
      min_regulon_size
    ))
  }

  gene_to_clusters <- list()
  for (gene in unique(adjacency[["TF"]])) {
    gene_pattern <- paste0(
      "(^|__|_)",
      gene,
      "($|_)"
    )
    matching_rows <- grepl(
      gene_pattern,
      motif_tbl[["tf"]],
      ignore.case = TRUE
    )
    if (any(matching_rows)) {
      gene_to_clusters[[gene]] <- unique(
        motif_tbl[["motif"]][matching_rows]
      )
    }
  }

  log_message(
    "  Indexed {.val {length(all_genes)}} genes, {.val {length(cluster_names_all)}} clusters, {.val {length(gene_to_clusters)}} gene-to-cluster mappings",
    verbose = verbose
  )

  tfs <- unique(adjacency[["TF"]])
  regulons <- list()

  for (tf in tfs) {
    tf_adj <- adjacency[adjacency[["TF"]] == tf, , drop = FALSE]
    tf_adj <- tf_adj[order(-as.numeric(tf_adj[["importance"]])), , drop = FALSE]
    target_genes <- unique(tf_adj[["target"]])
    target_genes <- utils::head(target_genes, max_targets * 2)
    target_set <- intersect(target_genes, all_genes)

    if (length(target_set) < min_regulon_size) {
      next
    }

    tf_clusters <- gene_to_clusters[[tf]]
    if (is.null(tf_clusters) || length(tf_clusters) == 0) {
      next
    }

    cluster_scores <- data.frame(
      cluster = character(0),
      auc = numeric(0),
      nes = numeric(0),
      leading_edge_genes = I(list()),
      stringsAsFactors = FALSE
    )

    for (rm in rank_matrices) {
      gene_idx <- match(target_set, rm[["genes"]])
      gene_idx <- gene_idx[!is.na(gene_idx)]
      if (length(gene_idx) < min_regulon_size) {
        next
      }

      cluster_idx <- which(rm[["clusters"]] %in% tf_clusters)
      if (length(cluster_idx) == 0) {
        next
      }

      for (ci in cluster_idx) {
        cluster_name <- rm[["clusters"]][ci]
        rankings <- as.numeric(rm[["data"]][ci, rm[["genes"]]])
        names(rankings) <- rm[["genes"]]

        n_total <- length(rankings)
        target_ranks_sub <- rankings[target_set]
        target_ranks_sub <- sort(target_ranks_sub[!is.na(target_ranks_sub)])

        if (length(target_ranks_sub) < min_regulon_size) {
          next
        }

        max_r <- max(rankings, na.rm = TRUE)
        if (max_r <= 0) {
          next
        }
        target_ranks_sub <- target_ranks_sub / max_r

        n_targets <- length(target_ranks_sub)

        x <- c(0, target_ranks_sub, 1)
        y <- c(0, seq_len(n_targets) / n_targets, 1)
        auc <- sum(diff(x) * (y[-1] + y[-length(y)]) / 2)

        expected_auc <- 0.5
        nes <- (auc - expected_auc) /
          max(0.01, sqrt(expected_auc * (1 - expected_auc) / n_total))

        if (auc >= auc_threshold && nes >= nes_threshold) {
          n_lead <- max(2, as.integer(n_targets / 3))
          leading <- names(target_ranks_sub)[seq_len(min(n_lead, n_targets))]
          cluster_scores <- rbind(
            cluster_scores,
            data.frame(
              cluster = cluster_name,
              auc = auc,
              nes = nes,
              leading_edge_genes = I(list(unique(c(target_set, leading)))),
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }

    if (nrow(cluster_scores) == 0) {
      next
    }

    cluster_scores <- cluster_scores[
      order(-cluster_scores[["nes"]]),
      ,
      drop = FALSE
    ]
    top_clusters <- utils::head(cluster_scores, 5)

    regulon_genes <- unique(unlist(top_clusters[["leading_edge_genes"]]))
    regulon_genes <- intersect(regulon_genes, target_set)
    regulon_genes <- utils::head(regulon_genes, max_targets)

    if (length(regulon_genes) >= min_regulon_size) {
      regulons[[tf]] <- regulon_genes
    }
  }

  missing_tfs <- setdiff(tfs, names(regulons))
  if (length(missing_tfs) > 0) {
    log_message(
      "  {.val {length(missing_tfs)}} TFs had no enriched motifs; using top-N importance as fallback",
      message_type = "info",
      verbose = verbose
    )
    fallback <- build_regulons(
      adjacency[adjacency[["TF"]] %in% missing_tfs, , drop = FALSE],
      max_targets,
      min_regulon_size
    )
    regulons <- c(regulons, fallback)
  }

  log_message(
    "{.pkg cisTarget} produced {.val {length(regulons)}} regulons (motif-enriched: {.val {length(regulons) - length(missing_tfs)}}, fallback: {.val {length(missing_tfs)}})",
    verbose = verbose
  )

  regulons[order(names(regulons))]
}

scenic_reference <- function(
  species,
  genome = NULL,
  data_dir = NULL,
  ranking_dbs = NULL,
  motif_annotations = NULL,
  regulators = NULL,
  verbose = TRUE
) {
  species_config <- scenic_species_config(species, genome = genome)
  missing_ranking_dbs <- is.null(ranking_dbs) || length(ranking_dbs) == 0
  missing_motif_annotations <- is.null(motif_annotations) ||
    length(motif_annotations) == 0
  missing_regulators <- is.null(regulators) || length(regulators) == 0
  needs_download <- missing_ranking_dbs ||
    missing_motif_annotations ||
    missing_regulators

  reference_dir <- NULL
  if (isTRUE(needs_download)) {
    reference_dir <- scenic_reference_dir(
      data_dir = data_dir,
      species_key = species_config[["key"]]
    )
    dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
    log_message(
      "Preparing {.pkg SCENIC} reference data for {.val {species_config[['label']]}} in {.path {reference_dir}}",
      verbose = verbose
    )

    reference_files <- species_config[["files"]]
    if (isTRUE(missing_ranking_dbs)) {
      ranking_rows <- reference_files[
        reference_files[["role"]] == "ranking_dbs",
        ,
        drop = FALSE
      ]
      ranking_dbs <- scenic_dl_refs(
        ranking_rows,
        reference_dir,
        verbose = verbose
      )
    }
    if (isTRUE(missing_motif_annotations)) {
      motif_rows <- reference_files[
        reference_files[["role"]] == "motif_annotations",
        ,
        drop = FALSE
      ]
      motif_annotations <- scenic_dl_refs(
        motif_rows,
        reference_dir,
        verbose = verbose
      )
    }
    if (isTRUE(missing_regulators)) {
      tf_rows <- reference_files[
        reference_files[["role"]] == "tf_list",
        ,
        drop = FALSE
      ]
      regulators <- scenic_dl_refs(
        tf_rows,
        reference_dir,
        verbose = verbose
      )
    }
  } else if (!is.null(data_dir)) {
    reference_dir <- scenic_reference_dir(
      data_dir = data_dir,
      species_key = species_config[["key"]]
    )
  }

  list(
    species = species_config[["label"]],
    genome = species_config[["genome"]],
    data_dir = reference_dir,
    ranking_dbs = ranking_dbs,
    motif_annotations = motif_annotations,
    regulators = regulators
  )
}

scenic_species_config <- function(species, genome = NULL) {
  if (length(species) > 1L) {
    species <- species[[1L]]
  }
  if (
    !is.character(species) ||
      length(species) != 1L ||
      is.na(species) ||
      !nzchar(species)
  ) {
    log_message(
      "{.arg species} must be one supported species name",
      message_type = "error"
    )
  }
  species_key <- tolower(gsub("[ .-]+", "_", species))
  genome_key <- if (is.null(genome)) NULL else tolower(gsub("[ .-]+", "_", genome))
  species_key <- switch(
    species_key,
    homo_sapiens = "human",
    mus_musculus = "mouse",
    drosophila_melanogaster = "fly",
    NULL
  )
  if (is.null(species_key)) {
    log_message(
      "{.arg species} must be one of {.val Homo_sapiens}, {.val Mus_musculus}, or {.val Drosophila_melanogaster}",
      message_type = "error"
    )
  }

  cistarget_url <- "https://resources.aertslab.org/cistarget"
  genome_key <- switch(
    species_key,
    human = genome_key %||% "hg38",
    mouse = genome_key %||% "mm10",
    fly = genome_key %||% "dm6"
  )
  genome_key <- switch(
    genome_key,
    hg38 = "hg38",
    grch38 = "hg38",
    hg19 = "hg19",
    grch37 = "hg19",
    mm10 = "mm10",
    grcm38 = "mm10",
    dm6 = "dm6",
    genome_key
  )
  expected_genomes <- switch(
    species_key,
    human = c("hg38", "hg19"),
    mouse = "mm10",
    fly = "dm6"
  )
  if (!genome_key %in% expected_genomes) {
    log_message(
      "{.arg genome} must be one of {.val {expected_genomes}} for {.arg species} {.val {species}}",
      message_type = "error"
    )
  }
  switch(
    species_key,
    human = if (identical(genome_key, "hg19")) {
      list(
        key = "human_hg19",
        label = "Homo_sapiens",
        genome = "hg19",
        files = data.frame(
          role = c("ranking_dbs", "ranking_dbs", "motif_annotations", "tf_list"),
          filename = c(
            "hg19-500bp-upstream-10species.mc9nr.genes_vs_motifs.rankings.feather",
            "hg19-tss-centered-10kb-10species.mc9nr.genes_vs_motifs.rankings.feather",
            "motifs-v9-nr.hgnc-m0.001-o0.0.tbl",
            "allTFs_hgnc.txt"
          ),
          url = c(
            paste0(
              cistarget_url,
              "/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-10species.mc9nr.genes_vs_motifs.rankings.feather"
            ),
            paste0(
              cistarget_url,
              "/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-10species.mc9nr.genes_vs_motifs.rankings.feather"
            ),
            paste0(
              cistarget_url,
              "/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
            ),
            paste0(cistarget_url, "/tf_lists/allTFs_hg38.txt")
          ),
          stringsAsFactors = FALSE
        )
      )
    } else {
      list(
      key = "human",
      label = "Homo_sapiens",
      genome = "hg38",
      files = data.frame(
        role = c("ranking_dbs", "ranking_dbs", "motif_annotations", "tf_list"),
        filename = c(
          "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
          "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
          "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
          "allTFs_hg38.txt"
        ),
        url = c(
          paste0(
            cistarget_url,
            "/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
          ),
          paste0(
            cistarget_url,
            "/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
          ),
          paste0(
            cistarget_url,
            "/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
          ),
          paste0(cistarget_url, "/tf_lists/allTFs_hg38.txt")
        ),
        stringsAsFactors = FALSE
      )
    )
    },
    mouse = list(
      key = "mouse",
      label = "Mus_musculus",
      genome = "mm10",
      files = data.frame(
        role = c("ranking_dbs", "ranking_dbs", "motif_annotations", "tf_list"),
        filename = c(
          "mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
          "mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
          "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl",
          "allTFs_mm.txt"
        ),
        url = c(
          paste0(
            cistarget_url,
            "/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
          ),
          paste0(
            cistarget_url,
            "/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
          ),
          paste0(
            cistarget_url,
            "/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"
          ),
          paste0(cistarget_url, "/tf_lists/allTFs_mm.txt")
        ),
        stringsAsFactors = FALSE
      )
    ),
    fly = list(
      key = "fly",
      label = "Drosophila_melanogaster",
      genome = "dm6",
      files = data.frame(
        role = c("ranking_dbs", "motif_annotations", "tf_list"),
        filename = c(
          "dm6_v10_clust.genes_vs_motifs.rankings.feather",
          "motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl",
          "allTFs_dmel.txt"
        ),
        url = c(
          paste0(
            cistarget_url,
            "/databases/drosophila_melanogaster/dm6/flybase_r6.02/mc_v10_clust/gene_based/dm6_v10_clust.genes_vs_motifs.rankings.feather"
          ),
          paste0(
            cistarget_url,
            "/motif2tf/motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl"
          ),
          paste0(cistarget_url, "/tf_lists/allTFs_dmel.txt")
        ),
        stringsAsFactors = FALSE
      )
    )
  )
}

scenic_reference_dir <- function(data_dir = NULL, species_key) {
  if (is.null(data_dir)) {
    return(file.path(
      tools::R_user_dir("scop", "data"),
      "SCENIC",
      species_key
    ))
  }
  if (
    !is.character(data_dir) ||
      length(data_dir) != 1L ||
      is.na(data_dir) ||
      !nzchar(data_dir)
  ) {
    log_message(
      "{.arg data_dir} must be one directory path",
      message_type = "error"
    )
  }
  normalizePath(data_dir, mustWork = FALSE)
}

scenic_dl_refs <- function(
  reference_files,
  reference_dir,
  verbose = TRUE
) {
  vapply(
    seq_len(nrow(reference_files)),
    function(idx) {
      scenic_download_reference_file(
        url = reference_files[["url"]][[idx]],
        dest = file.path(reference_dir, reference_files[["filename"]][[idx]]),
        verbose = verbose
      )
    },
    character(1)
  )
}

scenic_download_reference_file <- function(url, dest, verbose = TRUE) {
  if (file.exists(dest)) {
    log_message(
      "Using cached {.pkg SCENIC} reference file: {.path {dest}}",
      verbose = verbose
    )
    return(normalizePath(dest, mustWork = TRUE))
  }

  log_message(
    "Downloading {.pkg SCENIC} reference file: {.file {basename(dest)}}",
    verbose = verbose
  )
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  temp_dest <- paste0(dest, ".download")
  if (file.exists(temp_dest)) {
    unlink(temp_dest)
  }
  old_timeout <- getOption("timeout")
  options(timeout = max(3600, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  download_ok <- tryCatch(
    {
      utils::download.file(
        url = url,
        destfile = temp_dest,
        mode = "wb",
        quiet = !isTRUE(verbose)
      )
      TRUE
    },
    error = function(err) {
      log_message(
        "Failed to download {.pkg SCENIC} reference file from {.val {url}}: {conditionMessage(err)}",
        message_type = "error"
      )
    }
  )
  if (!isTRUE(download_ok) || !file.exists(temp_dest)) {
    log_message(
      "Failed to download {.pkg SCENIC} reference file from {.val {url}}",
      message_type = "error"
    )
  }
  if (!file.rename(temp_dest, dest)) {
    log_message(
      "Failed to move downloaded {.pkg SCENIC} reference file to {.path {dest}}",
      message_type = "error"
    )
  }
  normalizePath(dest, mustWork = TRUE)
}

scenic_prep_gene_arg <- function(
  x,
  arg,
  out_file,
  required = FALSE
) {
  genes <- scenic_read_gene_list_argument(x, arg = arg, required = required)
  if (is.null(genes)) {
    return(list(genes = NULL, file = NULL))
  }

  list(
    genes = genes,
    file = scenic_write_gene_list(genes, out_file)
  )
}

scenic_read_gene_list_argument <- function(x, arg, required = FALSE) {
  if (is.null(x) || length(x) == 0) {
    if (isTRUE(required)) {
      log_message(
        "{.val {arg}} must contain at least one gene",
        message_type = "error"
      )
    }
    return(NULL)
  }

  x <- as.character(x)
  from_file <- length(x) == 1L && file.exists(x)
  path_like <- length(x) == 1L &&
    (grepl("[/\\\\]", x) ||
      grepl("\\.(txt|csv|tsv|list)$", x, ignore.case = TRUE))
  if (!isTRUE(from_file) && isTRUE(path_like)) {
    log_message(
      "{.val {arg}} file does not exist: {.path {x}}",
      message_type = "error"
    )
  }

  if (isTRUE(from_file)) {
    x <- readLines(x, warn = FALSE)
    x <- x[!grepl("^\\s*#", x)]
    x <- unlist(strsplit(x, "[,\\t ]+"), use.names = FALSE)
  }

  genes <- trimws(x)
  genes <- genes[nzchar(genes)]
  genes <- unique(genes)
  if (length(genes) == 0 && isTRUE(required)) {
    log_message(
      "{.val {arg}} must contain at least one gene",
      message_type = "error"
    )
  }
  genes
}

scenic_write_gene_list <- function(genes, out_file) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(genes, out_file, useBytes = TRUE)
  normalizePath(out_file, mustWork = TRUE)
}

scenic_grn_inputs_changed <- function(params_file, params) {
  if (!file.exists(params_file)) {
    return(TRUE)
  }
  old_params <- tryCatch(
    readRDS(params_file),
    error = function(...) NULL
  )
  !identical(old_params, params)
}

scenic_flt_adj <- function(
  input_file,
  output_file,
  targets,
  verbose = TRUE
) {
  if (!file.exists(input_file)) {
    log_message(
      "Cannot find GRNBoost2 adjacency file: {.file {input_file}}",
      message_type = "error"
    )
  }
  adjacency <- utils::read.delim(
    input_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  target_col <- if ("target" %in% colnames(adjacency)) {
    "target"
  } else if ("Target" %in% colnames(adjacency)) {
    "Target"
  } else if (ncol(adjacency) >= 2L) {
    colnames(adjacency)[[2L]]
  } else {
    log_message(
      "GRNBoost2 adjacency table must contain a target column",
      message_type = "error"
    )
  }

  n_before <- nrow(adjacency)
  adjacency <- adjacency[adjacency[[target_col]] %in% targets, , drop = FALSE]
  if (nrow(adjacency) == 0) {
    log_message(
      "No GRNBoost2 edges remain after filtering by {.arg targets}",
      message_type = "error"
    )
  }
  utils::write.table(
    adjacency,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  log_message(
    "{.pkg SCENIC} target constraint kept {.val {nrow(adjacency)}} of {.val {n_before}} GRNBoost2 edge{?s}",
    verbose = verbose
  )
  invisible(output_file)
}

scenic_progress_init <- function(verbose = TRUE) {
  if (!isTRUE(verbose)) {
    return(NULL)
  }
  progress_state <- new.env(parent = emptyenv())
  progress_state$pb <- utils::txtProgressBar(
    min = 0,
    max = 100,
    initial = 0,
    style = 3
  )
  progress_state$value <- 0L
  progress_state
}

scenic_progress_step <- function(progress_state, value, label, verbose = TRUE) {
  value <- max(0L, min(100L, as.integer(value)))
  if (!is.null(progress_state)) {
    value <- max(progress_state$value, value)
    utils::setTxtProgressBar(progress_state$pb, value)
    progress_state$value <- value
    if (isTRUE(verbose)) cat("\n")
  }
  log_message(
    sprintf("[%d%%] %s", value, label),
    verbose = verbose
  )
  invisible(progress_state)
}

scenic_progress_close <- function(progress_state) {
  if (!is.null(progress_state)) {
    close(progress_state$pb)
  }
  invisible(NULL)
}

scenic_prepare_grn_matrix <- function(
  counts,
  ranking_genes,
  genes = NULL,
  min_expr_cells = 3,
  verbose = TRUE
) {
  grn_matrix <- Matrix::t(counts)
  expr_in_cells <- Matrix::colSums(grn_matrix > 0)
  grn_matrix <- grn_matrix[, expr_in_cells >= min_expr_cells, drop = FALSE]
  if (!is.null(genes)) {
    grn_matrix <- grn_matrix[, colnames(grn_matrix) %in% genes, drop = FALSE]
  }
  grn_matrix <- grn_matrix[,
    colnames(grn_matrix) %in% ranking_genes,
    drop = FALSE
  ]
  if (ncol(grn_matrix) == 0) {
    log_message(
      "No genes remain after filtering by expression and cisTarget database",
      message_type = "error"
    )
  }
  log_message(
    "{.pkg SCENIC} GRN input matrix: {.val {nrow(grn_matrix)}} cells x {.val {ncol(grn_matrix)}} genes",
    verbose = verbose
  )
  grn_matrix
}

scenic_read_regulon_txt <- function(txt_file) {
  if (!file.exists(txt_file)) {
    log_message(
      "Cannot find SCENIC regulon table: {.file {txt_file}}",
      message_type = "error"
    )
  }
  regulons <- utils::read.table(
    txt_file,
    sep = "\t",
    header = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )
  if (ncol(regulons) < 3) {
    log_message(
      "SCENIC regulon table must contain regulon, motif, and target columns",
      message_type = "error"
    )
  }
  colnames(regulons)[1:3] <- c("regulon", "motif", "target")
  regulons[, 1:3, drop = FALSE]
}

scenic_regulon_list <- function(regulon_tbl) {
  rgnames <- unique(regulon_tbl[["regulon"]])
  regulon_list <- lapply(rgnames, function(rg) {
    target <- regulon_tbl[regulon_tbl[["regulon"]] == rg, "target"]
    unique(unlist(strsplit(target, ",", fixed = TRUE), use.names = FALSE))
  })
  names(regulon_list) <- sub("[0-9]+g", "+", rgnames)
  regulon_list
}

scenic_compute_aucell_score <- function(
  counts,
  regulon_list,
  min_regulon_size = 10,
  batch_size = 500,
  cores = 1,
  backend = c("r", "cpp"),
  cpp_strategy = c("full", "sparse", "topk"),
  verbose = TRUE
) {
  backend <- match.arg(backend)
  cpp_strategy <- match.arg(cpp_strategy)
  regulon_list <- lapply(regulon_list, intersect, rownames(counts))
  regulon_list <- regulon_list[lengths(regulon_list) >= min_regulon_size]
  if (length(regulon_list) == 0) {
    log_message(
      "No regulons remain after matching targets to expression matrix",
      message_type = "error"
    )
  }

  if (identical(backend, "cpp")) {
    log_message(
      "Calculating AUCell regulon activity scores with {.arg backend = 'cpp'} and {.arg cpp_strategy} = {.val {cpp_strategy}}",
      verbose = verbose
    )
    scores <- run_aucell_scores(
      expr_counts = counts,
      gene_sets = regulon_list,
      strategy = cpp_strategy
    )
    return(as.data.frame(scores, check.names = FALSE)[
      colnames(counts),
      ,
      drop = FALSE
    ])
  }

  check_r("AUCell", verbose = FALSE)
  batches <- split(
    seq_len(ncol(counts)),
    ceiling(seq_along(colnames(counts)) / batch_size)
  )
  log_message(
    "Calculating AUCell regulon activity scores with {.val {length(batches)}} batch{?es} and {.val {cores}} core{?s}",
    verbose = verbose
  )
  calc_auc <- function(cells_idx) {
    scenic_calc_auc_batch(
      counts = counts[, cells_idx, drop = FALSE],
      regulon_list = regulon_list
    )
  }

  scores_list <- thisutils::parallelize_fun(
    x = batches,
    fun = calc_auc,
    cores = cores,
    progress_bar_width = min(50L, length(batches)),
    verbose = verbose
  )

  scores <- do.call(rbind, scores_list)
  scores <- as.data.frame(scores, check.names = FALSE)
  scores[colnames(counts), , drop = FALSE]
}

scenic_calc_auc_batch <- function(counts, regulon_list) {
  rankings <- AUCell::AUCell_buildRankings(
    counts,
    nCores = 1,
    plotStats = FALSE,
    verbose = FALSE
  )
  auc <- AUCell::AUCell_calcAUC(
    regulon_list,
    rankings,
    nCores = 1,
    verbose = FALSE
  )
  auc_mat <- as.matrix(AUCell::getAUC(auc))
  auc_mat <- t(auc_mat)
  auc_mat
}
