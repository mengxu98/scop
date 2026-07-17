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
#' @param min_expr_cells Minimum number of cells where a gene must be detected
#' before GRNBoost2. To run SCENIC on metacells, first create a metacell-level
#' object with [RunMetaCell()] and pass that object to `RunSCENIC()`.
#' @param min_regulon_size Minimum regulon size kept after `scenic ctx`.
#' @param max_regulon_targets Maximum target genes retained per C++ regulon.
#' The default preserves the native fast path's bounded regulon size.
#' @param include_negative_regulons Whether the C++ backend should also build
#' negatively correlated regulons and label them as `TF(-)`. The default
#' matches pySCENIC's positive-regulon workflow and labels C++ regulons as
#' `TF(+)`.
#' @param backend SCENIC backend. `"cpp"` uses the R/C++ path and
#' `"python"` uses the Python pySCENIC path. The selected backend controls GRN,
#' cisTarget pruning, and AUCell scoring together.
#' @param n_rounds Number of boosting rounds used by GRNBoost2.
#' @param learning_rate Learning rate used by GRNBoost2.
#' @param max_depth Maximum tree depth used by GRNBoost2.
#' @param max_features Fraction of features sampled by GRNBoost2.
#' @param subsample Row subsampling fraction used by GRNBoost2.
#' @param early_stop_window_length Early-stopping window used by GRNBoost2.
#' @param cores Number of workers used by GRNBoost2, `scenic ctx`, and
#' AUCell scoring. If multicore execution is not supported, this is
#' automatically reduced to one core.
#' @param seed Random seed used by GRNBoost2 and Seurat overclustering.
#' @param force Whether to rebuild existing SCENIC outputs.
#' @param assay_name Name of the assay used to store regulon activity scores.
#' @param tool_name Name of the `srt@tools` entry.
#' @param return_seurat Whether to return the modified Seurat object. If
#' `FALSE`, a result list is returned.
#' @param envname Python environment used for SCENIC. If `NULL`, the isolated
#' `"scenic_env"` environment is used.
#' @param conda The path or command name of a conda-compatible executable.
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
  min_expr_cells = 3,
  min_regulon_size = 10,
  max_regulon_targets = 50,
  include_negative_regulons = FALSE,
  backend = c("cpp", "python"),
  n_rounds = 5000,
  learning_rate = 0.01,
  max_depth = 3,
  max_features = 0.1,
  subsample = 0.9,
  early_stop_window_length = 25,
  cores = 1,
  seed = 1234,
  force = FALSE,
  assay_name = "scenic",
  tool_name = "SCENIC",
  return_seurat = TRUE,
  envname = NULL,
  conda = "auto",
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  species <- match.arg(species)
  backend <- match.arg(backend)

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
      min_expr_cells = min_expr_cells,
      min_regulon_size = min_regulon_size,
      max_regulon_targets = max_regulon_targets,
      include_negative_regulons = include_negative_regulons,
      ranking_dbs = ranking_dbs,
      motif_annotations = motif_annotations,
      n_rounds = n_rounds,
      learning_rate = learning_rate,
      max_depth = max_depth,
      max_features = max_features,
      subsample = subsample,
      early_stop_window_length = early_stop_window_length,
      cores = cores,
      seed = seed,
      force = force,
      assay_name = assay_name,
      tool_name = tool_name,
      return_seurat = return_seurat,
      verbose = verbose
    ))
  }

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
  scenic_progress_step(
    progress_state,
    value = 15,
    label = "Checking SCENIC Python packages",
    verbose = verbose
  )
  PrepareEnv(
    envname = envname,
    conda = conda,
    modules = "scenic",
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

  grn_counts <- counts
  grn_count_source <- "input count matrix"
  grn_count_unit <- "cells"
  scenic_progress_step(
    progress_state,
    value = 30,
    label = "Using input counts for GRNBoost2",
    verbose = verbose
  )
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
    scenic_ranking_gene_columns(
      as.character(unlist(functions$SCENICRankingGenes(db), use.names = FALSE))
    )
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
  grn_gene_order <- if (is.null(targets)) {
    union(regulators, colnames(grn_matrix))
  } else {
    union(regulators, targets)
  }
  grn_gene_order <- intersect(grn_gene_order, colnames(grn_matrix))
  grn_matrix <- grn_matrix[, grn_gene_order, drop = FALSE]

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
  auc_file <- file.path(work_dir, paste0(prefix, "_regulon_activity_score.csv"))
  ras_file <- file.path(work_dir, paste0(prefix, "_regulon_activity_score.rds"))
  regulon_list_file <- file.path(work_dir, paste0(prefix, "_regulon_list.rds"))
  stage_timing_sec <- c(
    grn = NA_real_,
    cistarget = NA_real_,
    regulon_conversion = NA_real_,
    aucell = NA_real_
  )
  stage_time <- function(expr) {
    start <- proc.time()[["elapsed"]]
    value <- force(expr)
    list(
      value = value,
      elapsed = proc.time()[["elapsed"]] - start
    )
  }

  scenic_progress_step(
    progress_state,
    value = 55,
    label = "Running GRNBoost2",
    verbose = verbose
  )
  timed <- stage_time({
    functions$RunSCENICGrn(
      expression_mtx = expr_csv,
      regulators = regulators_file,
      adj_output = grn_adj_file,
      n_rounds = as.integer(n_rounds),
      learning_rate = as.numeric(learning_rate),
      max_depth = as.integer(max_depth),
      max_features = as.numeric(max_features),
      subsample = as.numeric(subsample),
      early_stop_window_length = as.integer(early_stop_window_length),
      cores = as.integer(cores),
      seed = as.integer(seed),
      force = isTRUE(grn_force),
      verbose = isTRUE(verbose)
    )
  })
  stage_timing_sec[["grn"]] <- timed[["elapsed"]]
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
  timed <- stage_time({
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
  })
  stage_timing_sec[["cistarget"]] <- timed[["elapsed"]]

  scenic_progress_step(
    progress_state,
    value = 82,
    label = "Converting SCENIC regulons",
    verbose = verbose
  )
  if (isTRUE(grn_force) || !file.exists(gmt_file) || !file.exists(txt_file)) {
    timed <- stage_time({
      functions$SCENICRegulonsToFiles(
        regulon_file = ctx_file,
        gmt_file = gmt_file,
        txt_file = txt_file,
        min_regulon_size = as.integer(min_regulon_size)
      )
    })
    stage_timing_sec[["regulon_conversion"]] <- timed[["elapsed"]]
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
    timed <- stage_time({
      functions$RunSCENICAUCell(
        expression_mtx = expr_csv,
        regulons_gmt = gmt_file,
        auc_output = auc_file,
        cores = as.integer(cores),
        seed = as.integer(seed),
        force = isTRUE(grn_force),
        verbose = isTRUE(verbose)
      )
    })
    stage_timing_sec[["aucell"]] <- timed[["elapsed"]]
    ras_mat <- scenic_read_python_aucell(auc_file)
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
  score_mat <- Matrix::t(Matrix::Matrix(
    if (is.data.frame(ras_mat)) as.matrix(ras_mat) else ras_mat,
    sparse = TRUE
  ))
  dimnames(score_mat) <- list(colnames(ras_mat), rownames(ras_mat))
  result <- list(
    scores = score_mat,
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
      regulon_activity_score_csv = auc_file,
      regulon_activity_score = ras_file,
      regulon_list = regulon_list_file
    ),
    parameters = list(
      backend = "python",
      grn_method = "grnboost2",
      cistarget_method = "pySCENIC_ctx",
      method = "pySCENIC",
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
      min_expr_cells = min_expr_cells,
      min_regulon_size = min_regulon_size,
      n_rounds = n_rounds,
      learning_rate = learning_rate,
      max_depth = max_depth,
      max_features = max_features,
      subsample = subsample,
      early_stop_window_length = early_stop_window_length,
      cores = cores,
      aucell_method = "pySCENIC",
      seed = seed,
      assay_name = assay_name,
      tool_name = tool_name,
      envname = envname,
      progress = verbose
    ),
    details = list(
      stage_timing_sec = stage_timing_sec
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
  min_expr_cells,
  min_regulon_size,
  max_regulon_targets,
  include_negative_regulons,
  ranking_dbs,
  motif_annotations,
  n_rounds,
  learning_rate,
  max_depth,
  max_features,
  subsample,
  early_stop_window_length,
  cores,
  seed,
  force,
  assay_name,
  tool_name,
  return_seurat,
  verbose
) {
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  max_regulon_targets <- scenic_normalize_max_regulon_targets(max_regulon_targets)
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

  work_dir <- normalizePath(work_dir, mustWork = FALSE)
  ras_file <- file.path(
    work_dir,
    paste0(prefix, "_regulon_activity_score_cpp.rds")
  )
  adj_file <- file.path(work_dir, paste0(prefix, "_adj_cpp.tsv"))
  grn_input_params_file <- file.path(
    work_dir,
    paste0(prefix, "_grn_input_params_cpp.rds")
  )
  regulon_file <- file.path(work_dir, paste0(prefix, "_regulon_list_cpp.rds"))
  regulon_params_file <- file.path(work_dir, paste0(prefix, "_regulon_params_cpp.rds"))
  regulators_file <- file.path(work_dir, paste0(prefix, "_regulators.txt"))

  reference_data <- scenic_reference(
    species = species,
    genome = genome,
    data_dir = data_dir,
    ranking_dbs = ranking_dbs,
    motif_annotations = motif_annotations,
    regulators = regulators,
    verbose = verbose
  )
  species <- reference_data[["species"]]
  genome <- reference_data[["genome"]]
  data_dir <- reference_data[["data_dir"]]
  ranking_dbs <- reference_data[["ranking_dbs"]]
  motif_annotations <- reference_data[["motif_annotations"]]
  regulators_raw <- reference_data[["regulators"]]
  if (is.null(regulators_raw) || length(regulators_raw) == 0) {
    log_message(
      "{.arg regulators} must be provided or resolvable from {.arg species}",
      message_type = "error"
    )
  }
  regulators <- scenic_read_gene_list_argument(
    regulators_raw,
    arg = "regulators",
    required = TRUE
  )
  if (length(regulators) == 0) {
    log_message(
      "No regulators found in reference data",
      message_type = "error"
    )
  }

  counts <- GetAssayData5(srt, assay = assay, layer = layer)
  grn_matrix <- counts

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

  ref_data <- reference_data
  if (!is.null(ref_data) && length(ref_data[["ranking_dbs"]]) > 0) {
    check_r("arrow", verbose = FALSE)
    ranking_gene_lists <- lapply(ref_data[["ranking_dbs"]], function(db) {
      scenic_ranking_gene_columns(scenic_feather_column_names(db))
    })
    ranking_genes <- Reduce(intersect, ranking_gene_lists)
    keep_ranking <- rownames(grn_matrix) %in% ranking_genes
    if (any(keep_ranking)) {
      grn_matrix <- grn_matrix[keep_ranking, , drop = FALSE]
      regulators <- intersect(regulators, rownames(grn_matrix))
      if (!is.null(targets)) {
        targets <- intersect(targets, rownames(grn_matrix))
      }
    }
    if (length(regulators) == 0) {
      log_message(
        "No {.arg regulators} remain after matching expression features and cisTarget databases",
        message_type = "error"
      )
    }
    if (!is.null(targets) && length(targets) == 0) {
      log_message(
        "No {.arg targets} remain after matching expression features and cisTarget databases",
        message_type = "error"
      )
    }
  }
  regulators_file <- scenic_write_gene_list(regulators, regulators_file)

  grn_input_params <- list(
    regulators = regulators,
    targets = targets,
    genes = rownames(grn_matrix),
    n_rounds = as.integer(max(1L, n_rounds)),
    learning_rate = as.numeric(learning_rate),
    max_depth = as.integer(max(1L, max_depth)),
    max_features = as.numeric(max_features),
    subsample = as.numeric(subsample),
    early_stop_window_length = as.integer(max(0L, early_stop_window_length)),
    seed = as.integer(seed %||% 1234L)
  )
  grn_inputs_changed <- scenic_grn_inputs_changed(
    grn_input_params_file,
    grn_input_params
  )

  log_message(
    "Running SCENIC (grnboost2) on {.val {ncol(grn_matrix)}} cells with {.val {length(regulators)}} TFs",
    verbose = verbose
  )

  grn_force <- isTRUE(force) ||
    !file.exists(adj_file) ||
    isTRUE(grn_inputs_changed)
  stage_timing_sec <- c(
    grn = NA_real_,
    cistarget = NA_real_,
    aucell = NA_real_
  )
  stage_time <- function(expr) {
    start <- proc.time()[["elapsed"]]
    value <- force(expr)
    list(
      value = value,
      elapsed = proc.time()[["elapsed"]] - start
    )
  }
  if (isTRUE(grn_inputs_changed) && file.exists(adj_file) && isFALSE(force)) {
    log_message(
      "Rebuilding C++ {.pkg SCENIC} GRNBoost2 adjacency because input genes or GRN parameters changed",
      verbose = verbose
    )
  }
  if (grn_force) {
    timed <- stage_time({
      RunGRN(
        object = Matrix::t(grn_matrix),
        regulators = regulators,
        targets = targets,
        genes_in = "columns",
        grn_method = "grnboost2",
        backend = "cpp",
        output_file = adj_file,
        max_edges_per_target = Inf,
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
    })
    adjacency <- timed[["value"]]
    stage_timing_sec[["grn"]] <- timed[["elapsed"]]
    saveRDS(grn_input_params, grn_input_params_file)
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

  regulon_params <- list(
    include_negative_regulons = isTRUE(include_negative_regulons),
    min_regulon_size = as.integer(min_regulon_size),
    max_regulon_targets = max_regulon_targets,
    signed_regulon_names = TRUE
  )
  regulon_force <- grn_force ||
    !file.exists(regulon_file) ||
    scenic_grn_inputs_changed(regulon_params_file, regulon_params)
  if (regulon_force) {
    if (
      !is.null(ref_data) &&
        length(ref_data[["ranking_dbs"]]) > 0 &&
        !is.null(ref_data[["motif_annotations"]])
    ) {
      timed <- stage_time({
        cistarget2(
          adjacency = adjacency,
          expr_mtx = Matrix::t(grn_matrix),
          ranking_dbs = ref_data[["ranking_dbs"]],
          motif_annotations = ref_data[["motif_annotations"]],
          max_targets = max_regulon_targets,
          min_regulon_size = min_regulon_size,
          include_negative_regulons = include_negative_regulons,
          fallback = FALSE,
          cores = cores,
          verbose = verbose
        )
      })
      regulon_list <- timed[["value"]]
      stage_timing_sec[["cistarget"]] <- timed[["elapsed"]]
    } else {
      log_message(
        "C++ cisTarget reference data are required for {.arg backend = 'cpp'}",
        message_type = "error"
      )
    }
    saveRDS(regulon_list, regulon_file)
    saveRDS(regulon_params, regulon_params_file)
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
    "C++ regulon builder produced {.val {length(regulon_list)}} regulons",
    verbose = verbose
  )

  if (file.exists(ras_file) && isFALSE(regulon_force)) {
    log_message(
      "Reusing existing AUCell scores: {.file {ras_file}}",
      verbose = verbose
    )
    ras_mat <- readRDS(ras_file)
  } else {
    timed <- stage_time({
      scenic_compute_aucell_score(
        counts = grn_matrix,
        regulon_list = regulon_list,
        min_regulon_size = min_regulon_size,
        cores = cores,
        backend = "cpp",
        cpp_algorithm = "ctxcore",
        seed = seed,
        verbose = verbose
      )
    })
    ras_mat <- timed[["value"]]
    stage_timing_sec[["aucell"]] <- timed[["elapsed"]]
    saveRDS(ras_mat, ras_file)
  }
  ras_mat <- ras_mat[colnames(srt), , drop = FALSE]

  regulon_tbl <- scenic_regulon_table(regulon_list)
  score_mat <- Matrix::t(Matrix::Matrix(
    if (is.data.frame(ras_mat)) as.matrix(ras_mat) else ras_mat,
    sparse = TRUE
  ))
  dimnames(score_mat) <- list(colnames(ras_mat), rownames(ras_mat))

  result <- list(
    scores = score_mat,
    scores_cells_by_regulon = ras_mat,
    regulons = regulon_tbl,
    regulon_list = regulon_list,
    adjacency = adjacency,
    files = list(
      regulators = regulators_file,
      adj = adj_file,
      grn_input_params = grn_input_params_file,
      regulon_params = regulon_params_file,
      regulon_list = regulon_file,
      regulon_activity_score = ras_file
    ),
    parameters = list(
      backend = "cpp",
      grn_method = "grnboost2",
      cistarget_method = "cpp_motif",
      method = "cpp_cistarget",
      assay = assay,
      layer = layer,
      species = species,
      genome = genome,
      regulators = regulators,
      regulators_file = regulators_file,
      targets = targets,
      min_expr_cells = min_expr_cells,
      min_regulon_size = min_regulon_size,
      include_negative_regulons = include_negative_regulons,
      n_rounds = n_rounds,
      learning_rate = learning_rate,
      max_depth = max_depth,
      max_features = max_features,
      subsample = subsample,
      early_stop_window_length = early_stop_window_length,
      max_regulon_targets = max_regulon_targets,
      cores = cores,
      aucell_method = "cpp_ctxcore_grn_space",
      seed = seed,
      assay_name = assay_name,
      tool_name = tool_name
    ),
    details = list(
      stage_timing_sec = stage_timing_sec,
      cistarget_profile_sec = attr(regulon_list, "profile_sec", exact = TRUE)
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
  max_targets = Inf,
  min_regulon_size = 10,
  suffix = "(+)"
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
    if (is.finite(max_targets)) {
      targets <- utils::head(targets, max_targets)
    }
    targets
  })

  regulons <- regulons[lengths(regulons) >= min_regulon_size]
  names(regulons) <- scenic_regulon_name(names(regulons), suffix = suffix)
  regulons
}

scenic_feather_column_names <- function(path) {
  check_r("arrow", verbose = FALSE)
  tryCatch(
    names(arrow::open_dataset(path, format = "feather")),
    error = function(...) colnames(arrow::read_feather(path))
  )
}

scenic_ranking_index_columns <- function(columns) {
  intersect(c("motifs", "tracks", "regions", "genes"), columns)
}

scenic_ranking_gene_columns <- function(columns) {
  setdiff(columns, scenic_ranking_index_columns(columns))
}

cistarget2 <- function(
  adjacency,
  expr_mtx = NULL,
  ranking_dbs,
  motif_annotations,
  max_targets = Inf,
  min_regulon_size = 10,
  rank_threshold = 5000,
  nes_threshold = 3.0,
  auc_threshold = 0.05,
  include_negative_regulons = FALSE,
  fallback = FALSE,
  cores = 1,
  verbose = TRUE
) {
  check_r("arrow", verbose = FALSE)
  profile_time <- function(expr) {
    start <- proc.time()[["elapsed"]]
    value <- force(expr)
    list(
      value = value,
      elapsed = proc.time()[["elapsed"]] - start
    )
  }
  profile <- c(
    read_annotations = 0,
    read_rankings = 0,
    index_annotations = 0,
    module_build = 0,
    matrix_subset = 0,
    cistarget_stats = 0,
    annotation_filter = 0,
    leading_edge = 0
  )

  profiled <- profile_time({
    tryCatch(
      utils::read.table(
        motif_annotations,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE,
        comment.char = ""
      ),
      error = function(e) {
        log_message(
          "Failed to read motif annotations: {.val {conditionMessage(e)}}",
          message_type = "error"
        )
      }
    )
  })
  motif_tbl <- profiled[["value"]]
  profile[["read_annotations"]] <- profile[["read_annotations"]] + profiled[["elapsed"]]
  motif_col <- intersect(c("#motif_id", "motif_id", "motif"), colnames(motif_tbl))[1]
  gene_col <- intersect(c("gene_name", "tf", "TF"), colnames(motif_tbl))[1]
  if (is.na(motif_col) || is.na(gene_col)) {
    log_message(
      "{.arg motif_annotations} must contain motif and gene annotation columns",
      message_type = "error"
    )
  }
  motif_tbl[["motif"]] <- as.character(motif_tbl[[motif_col]])
  motif_tbl[["tf"]] <- as.character(motif_tbl[[gene_col]])
  qvalue_col <- intersect(
    c("motif_similarity_qvalue", "MotifSimilarityQvalue"),
    colnames(motif_tbl)
  )[1]
  orthology_col <- intersect(
    c("orthologous_identity", "OrthologousIdentity"),
    colnames(motif_tbl)
  )[1]
  if (!is.na(qvalue_col)) {
    motif_tbl[[".motif_similarity_qvalue"]] <- as.numeric(motif_tbl[[qvalue_col]])
  } else {
    motif_tbl[[".motif_similarity_qvalue"]] <- 0
  }
  if (!is.na(orthology_col)) {
    motif_tbl[[".orthologous_identity"]] <- as.numeric(motif_tbl[[orthology_col]])
  } else {
    motif_tbl[[".orthologous_identity"]] <- 1
  }
  motif_tbl <- motif_tbl[
    motif_tbl[[".motif_similarity_qvalue"]] <= 0.001 &
      motif_tbl[[".orthologous_identity"]] >= 0,
    ,
    drop = FALSE
  ]

  motif_to_tf <- split(motif_tbl[["tf"]], motif_tbl[["motif"]])
  motif_names <- names(motif_to_tf)
  profiled <- profile_time({
    motif_tbl <- motif_tbl[
      order(
        motif_tbl[["tf"]],
        -motif_tbl[[".motif_similarity_qvalue"]],
        motif_tbl[[".orthologous_identity"]]
      ),
      ,
      drop = FALSE
    ]
    motif_key <- paste(motif_tbl[["tf"]], motif_tbl[["motif"]], sep = "\r")
    motif_tbl <- motif_tbl[!duplicated(motif_key, fromLast = TRUE), , drop = FALSE]
    split(motif_tbl, motif_tbl[["tf"]], drop = TRUE)
  })
  motif_tbl_by_tf <- profiled[["value"]]
  profile[["index_annotations"]] <- profile[["index_annotations"]] + profiled[["elapsed"]]

  log_message(
    "C++ cisTarget: loaded {.val {length(motif_names)}} motifs from annotations",
    verbose = verbose
  )

  rank_matrices <- list()
  all_genes <- character(0)
  cluster_names_all <- character(0)
  ranking_required_genes <- unique(c(adjacency[["TF"]], adjacency[["target"]]))

  for (db_path in ranking_dbs) {
    log_message(
      "  Reading ranking database: {.file {basename(db_path)}}",
      verbose = verbose
    )
    db_columns <- scenic_feather_column_names(db_path)
    index_cols <- scenic_ranking_index_columns(db_columns)
    index_col <- rev(index_cols)[1]
    if (is.na(index_col)) {
      log_message(
        "Ranking database missing motif/index column. Skipping.",
        message_type = "warning",
        verbose = verbose
      )
      next
    }
    gene_cols_all <- scenic_ranking_gene_columns(db_columns)
    gene_cols <- intersect(ranking_required_genes, gene_cols_all)
    if (length(gene_cols) < min_regulon_size) {
      log_message(
        "Ranking database has fewer than {.val {min_regulon_size}} usable genes for the current adjacency. Skipping.",
        message_type = "warning",
        verbose = verbose
      )
      next
    }
    read_cols <- c(gene_cols, index_col)
    profiled <- profile_time({
      tryCatch(
        arrow::read_feather(db_path, col_select = read_cols),
        error = function(e) {
          log_message(
            "Failed to read {.file {db_path}}: {.val {conditionMessage(e)}}",
            message_type = "warning",
            verbose = verbose
          )
          NULL
        }
      )
    })
    db <- profiled[["value"]]
    profile[["read_rankings"]] <- profile[["read_rankings"]] + profiled[["elapsed"]]
    if (is.null(db)) {
      next
    }

    clusters <- as.character(db[[index_col]])
    gene_cols <- setdiff(colnames(db), index_col)
    all_genes <- union(all_genes, gene_cols)
    cluster_names_all <- union(cluster_names_all, clusters)
    profiled <- profile_time({
      data <- as.matrix(db[, gene_cols])
      storage.mode(data) <- "integer"
      data
    })
    ranking_data <- profiled[["value"]]
    profile[["read_rankings"]] <- profile[["read_rankings"]] + profiled[["elapsed"]]

    rank_matrices[[length(rank_matrices) + 1]] <- list(
      clusters = clusters,
      genes = gene_cols,
      gene_index = stats::setNames(seq_along(gene_cols), gene_cols),
      data = ranking_data,
      total_genes = length(gene_cols_all),
      rank_cutoff = as.integer(round(auc_threshold * length(gene_cols_all))) - 1L
    )
  }
  rm(db)
  gc()

  if (length(rank_matrices) == 0) {
    log_message(
      "Failed to read any ranking databases. Falling back to rank-based approximation.",
      message_type = "warning",
      verbose = verbose
    )
    return(build_regulons(
      adjacency,
      max_targets,
      min_regulon_size,
      suffix = "(+)"
    ))
  }

  gene_to_clusters <- list()
  for (gene in unique(adjacency[["TF"]])) {
    matching_rows <- motif_tbl[["tf"]] == gene
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
  profiled <- profile_time({
    scenic_modules_from_adjacencies(
      adjacency = adjacency,
      expr_mtx = expr_mtx,
      min_genes = 20,
      top_n_targets = max_targets,
      keep_only_activating = !isTRUE(include_negative_regulons),
      verbose = verbose
    )
  })
  modules <- profiled[["value"]]
  profile[["module_build"]] <- profile[["module_build"]] + profiled[["elapsed"]]

  tf_importance_map <- split(adjacency[, c("target", "importance")], adjacency[["TF"]])
  tf_importance_map <- lapply(tf_importance_map, function(tf_df) {
    imp_vec <- stats::setNames(tf_df[["importance"]], tf_df[["target"]])
    imp_list <- split(imp_vec, names(imp_vec))
    lapply(imp_list, max)
  })

  process_module <- function(module) {
    local_profile <- profile * 0
    tf <- module[["tf"]]
    suffix <- module[["suffix"]] %||% "(+)"
    regulon_name <- scenic_regulon_name(tf, suffix = suffix)
    target_set <- intersect(module[["genes"]], all_genes)

    if (length(target_set) < min_regulon_size) {
      return(NULL)
    }

    tf_clusters <- gene_to_clusters[[tf]]
    if (is.null(tf_clusters) || length(tf_clusters) == 0) {
      return(NULL)
    }

    cluster_score_list <- list()

    for (rm in rank_matrices) {
      gene_idx <- unname(rm[["gene_index"]][target_set])
      gene_idx <- gene_idx[!is.na(gene_idx)]
      if (length(gene_idx) < min_regulon_size) {
        next
      }

      profiled <- profile_time({
        rm[["data"]][, gene_idx, drop = FALSE]
      })
      ranking_sub <- profiled[["value"]]
      local_profile[["matrix_subset"]] <- local_profile[["matrix_subset"]] + profiled[["elapsed"]]
      profiled <- profile_time({
        scenic_ctx_auc_avg2sd(
          ranks = ranking_sub,
          total_genes = rm[["total_genes"]],
          rank_threshold = as.integer(rank_threshold),
          rank_cutoff = rm[["rank_cutoff"]]
        )
      })
      cistarget_stats <- profiled[["value"]]
      local_profile[["cistarget_stats"]] <- local_profile[["cistarget_stats"]] + profiled[["elapsed"]]
      aucs <- cistarget_stats[["auc"]]
      auc_sd <- scenic_population_sd(aucs)
      if (!is.finite(auc_sd) || auc_sd <= 0) {
        next
      }
      ness <- (aucs - mean(aucs)) / auc_sd
      enriched_idx <- which(ness >= nes_threshold)
      if (length(enriched_idx) == 0) {
        next
      }

      profiled <- profile_time({
        annotated_tf <- motif_tbl_by_tf[[tf]]
        if (is.null(annotated_tf)) {
          motif_tbl[FALSE, , drop = FALSE]
        } else {
          annotated_tf[
            annotated_tf[["motif"]] %in% rm[["clusters"]][enriched_idx],
            ,
            drop = FALSE
          ]
        }
      })
      annotated_features <- profiled[["value"]]
      local_profile[["annotation_filter"]] <- local_profile[["annotation_filter"]] + profiled[["elapsed"]]
      if (nrow(annotated_features) == 0) {
        next
      }
      annotated_feature_idx <- match(annotated_features[["motif"]], rm[["clusters"]])
      annotated_feature_idx <- annotated_feature_idx[!is.na(annotated_feature_idx)]
      if (length(annotated_feature_idx) == 0) {
        next
      }
      annotated_features <- annotated_features[
        match(rm[["clusters"]][annotated_feature_idx], annotated_features[["motif"]]),
        ,
        drop = FALSE
      ]
      py_row_idx <- enriched_idx[seq_len(min(length(enriched_idx), nrow(annotated_features)))]
      avg_recovery <- cistarget_stats[["avg2sd_recovery"]]

      for (row_i in seq_along(py_row_idx)) {
        ci <- py_row_idx[[row_i]]
        cluster_name <- annotated_features[["motif"]][[row_i]]
        target_ranks_sub <- as.integer(ranking_sub[ci, ])
        names(target_ranks_sub) <- rm[["genes"]][gene_idx]
        profiled <- profile_time({
          scenic_cistarget_leading_edge(
            ranks = target_ranks_sub,
            avg2sd_recovery = avg_recovery
          )
        })
        leading <- profiled[["value"]]
        local_profile[["leading_edge"]] <- local_profile[["leading_edge"]] + profiled[["elapsed"]]

        if (length(leading) == 0) {
          next
        }

        cluster_score_list[[length(cluster_score_list) + 1L]] <- data.frame(
          cluster = cluster_name,
          auc = aucs[[annotated_feature_idx[[row_i]]]],
          nes = ness[[annotated_feature_idx[[row_i]]]],
          leading_edge_genes = I(list(unique(leading))),
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(cluster_score_list) == 0L) {
      return(list(tf = NULL, genes = NULL, profile = local_profile))
    }
    cluster_scores <- do.call(rbind, cluster_score_list)

    cluster_scores <- cluster_scores[
      order(-cluster_scores[["nes"]]),
      ,
      drop = FALSE
    ]

    # ── Merge leading-edge genes across enriched motifs (ctxcore-compatible) ──
    # ctxcore groups motifs by shared TF annotation (already done: all motifs
    # here are annotated to this TF).  It merges via Regulon.union() which takes
    # the *maximum GRN importance weight* per gene across all enriched motifs.
    # Reference: ctxcore/genesig.py:206-220 (Regulon.union),
    #            ctxcore/transform.py:262-307 (_regulon4group).
    n_motifs <- nrow(cluster_scores)
    tf_adj <- tf_importance_map[[tf]]
    if (!is.null(tf_adj) && length(tf_adj) > 0L) {
      gene_scores <- list()
      for (ri in seq_len(n_motifs)) {
        genes_i <- cluster_scores[["leading_edge_genes"]][[ri]]
        for (g in genes_i) {
          imp <- tf_adj[[g]]
          if (is.null(imp)) imp <- 0.0
          gene_scores[[g]] <- max(gene_scores[[g]] %||% 0.0, imp)
        }
      }
      regulon_genes <- names(sort(unlist(gene_scores), decreasing = TRUE))
    } else {
      regulon_genes <- unique(unlist(cluster_scores[["leading_edge_genes"]]))
    }
    regulon_genes <- intersect(regulon_genes, target_set)

    if (length(regulon_genes) > 0) {
      return(list(
        tf = tf,
        regulon = regulon_name,
        genes = regulon_genes,
        profile = local_profile
      ))
    }
    list(tf = NULL, genes = NULL, profile = local_profile)
  }

  module_cores <- max(1L, as.integer(cores))
  if (.Platform[["OS.type"]] == "unix" && module_cores > 1L && length(modules) > 1L) {
    module_results <- parallel::mclapply(
      modules,
      process_module,
      mc.cores = module_cores,
      mc.preschedule = TRUE
    )
  } else {
    module_results <- lapply(modules, process_module)
  }
  for (module_result in module_results) {
    if (is.null(module_result)) {
      next
    }
    profile <- profile + module_result[["profile"]]
    if (is.null(module_result[["tf"]])) {
      next
    }
    regulon_name <- module_result[["regulon"]] %||%
      scenic_regulon_name(module_result[["tf"]], suffix = "(+)")
    regulons[[regulon_name]] <- unique(c(
      regulons[[regulon_name]],
      module_result[["genes"]]
    ))
  }

  regulons <- regulons[lengths(regulons) >= min_regulon_size]
  positive_regulon_names <- scenic_regulon_name(tfs, suffix = "(+)")
  missing_tfs <- tfs[!positive_regulon_names %in% names(regulons)]
  if (length(missing_tfs) > 0 && isTRUE(fallback)) {
    log_message(
      "  {.val {length(missing_tfs)}} TFs had no enriched motifs; using top-N importance as fallback",
      verbose = verbose
    )
    fallback_regulons <- build_regulons(
      adjacency[adjacency[["TF"]] %in% missing_tfs, , drop = FALSE],
      max_targets,
      min_regulon_size,
      suffix = "(+)"
    )
    regulons <- c(regulons, fallback_regulons)
  }

  motif_enriched_n <- length(setdiff(tfs, missing_tfs))
  fallback_n <- if (isTRUE(fallback)) length(missing_tfs) else 0L
  log_message(
    "{.pkg cisTarget} produced {.val {length(regulons)}} regulons (motif-enriched: {.val {motif_enriched_n}}, fallback: {.val {fallback_n}})",
    verbose = verbose
  )
  log_message(
    "C++ cisTarget timing: {.val {paste(names(profile), sprintf('%.3fs', profile), sep = '=', collapse = ', ')}}",
    verbose = verbose
  )

  if (length(regulons) == 0L) {
    attr(regulons, "profile_sec") <- profile
    return(regulons)
  }
  regulons <- regulons[order(names(regulons))]
  attr(regulons, "profile_sec") <- profile
  regulons
}

scenic_normalize_max_regulon_targets <- function(max_regulon_targets) {
  if (is.null(max_regulon_targets) || length(max_regulon_targets) != 1L) {
    return(Inf)
  }
  max_regulon_targets <- suppressWarnings(as.numeric(max_regulon_targets))
  if (!is.finite(max_regulon_targets)) {
    return(Inf)
  }
  max(1L, as.integer(max_regulon_targets))
}

scenic_cistarget_auc <- function(
  rankings,
  total_genes,
  rank_threshold = 5000,
  auc_threshold = 0.05
) {
  rank_threshold <- min(as.integer(rank_threshold), as.integer(total_genes) - 1L)
  rank_cutoff <- as.integer(round(auc_threshold * total_genes)) - 1L
  if (rank_cutoff < 1L || rank_cutoff > rank_threshold) {
    log_message(
      "{.arg auc_threshold} and {.arg rank_threshold} produce an invalid cisTarget rank cutoff",
      message_type = "error"
    )
  }
  n_genes <- ncol(rankings)
  max_auc <- (rank_cutoff + 1L) * n_genes
  apply(rankings, 1, function(ranks) {
    ranks <- ranks[is.finite(ranks) & ranks <= rank_cutoff]
    if (length(ranks) == 0) {
      return(0)
    }
    x <- c(sort(as.integer(ranks)), rank_cutoff)
    y <- seq_len(length(x) - 1L)
    sum(diff(x) * y) / max_auc
  })
}

scenic_cistarget_recovery <- function(
  ranks,
  rank_threshold = 5000
) {
  ranks <- ranks[is.finite(ranks) & ranks < rank_threshold]
  if (length(ranks) == 0) {
    return(rep(0, rank_threshold))
  }
  recovered <- tabulate(as.integer(ranks) + 1L, nbins = rank_threshold)
  cumsum(recovered)
}

scenic_cistarget_avg2sd_recovery <- function(
  rankings,
  total_genes,
  rank_threshold = 5000
) {
  rank_threshold <- min(as.integer(rank_threshold), as.integer(total_genes) - 1L)
  n_features <- nrow(rankings)
  if (n_features == 0L) {
    return(rep(0, rank_threshold))
  }
  recovery_sum <- numeric(rank_threshold)
  recovery_sumsq <- numeric(rank_threshold)
  for (idx in seq_len(n_features)) {
    recovery <- scenic_cistarget_recovery(rankings[idx, ], rank_threshold = rank_threshold)
    recovery_sum <- recovery_sum + recovery
    recovery_sumsq <- recovery_sumsq + recovery * recovery
  }
  recovery_mean <- recovery_sum / n_features
  recovery_var <- pmax(recovery_sumsq / n_features - recovery_mean * recovery_mean, 0)
  recovery_mean + 2 * sqrt(recovery_var)
}

scenic_population_sd <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  sqrt(mean((x - mean(x))^2))
}

scenic_cistarget_leading_edge <- function(
  ranks,
  avg2sd_recovery
) {
  rank_threshold <- length(avg2sd_recovery)
  recovery <- scenic_cistarget_recovery(ranks, rank_threshold = rank_threshold)
  rank_at_max <- which.max(recovery - avg2sd_recovery) - 1L
  genes <- names(ranks)[is.finite(ranks) & ranks <= rank_at_max]
  genes[order(ranks[genes])]
}

scenic_modules_from_adjacencies <- function(
  adjacency,
  expr_mtx = NULL,
  thresholds = c(0.75, 0.9),
  top_n_targets = 50,
  top_n_regulators = c(5, 10, 50),
  min_genes = 20,
  rho_threshold = 0.03,
  keep_only_activating = TRUE,
  verbose = TRUE
) {
  if (!all(c("TF", "target", "importance") %in% colnames(adjacency))) {
    log_message(
      "{.arg adjacency} must have columns TF, target, importance",
      message_type = "error"
    )
  }
  adjacency <- adjacency[!is.na(adjacency[["importance"]]), , drop = FALSE]
  adjacency[["importance"]] <- as.numeric(adjacency[["importance"]])
  adjacency <- adjacency[order(-adjacency[["importance"]]), , drop = FALSE]
  threshold_values <- stats::quantile(
    adjacency[["importance"]],
    probs = thresholds,
    names = FALSE,
    na.rm = TRUE
  )

  if (!is.null(expr_mtx)) {
    adjacency <- scenic_add_correlation(
      adjacency = adjacency,
      expr_mtx = expr_mtx,
      rho_threshold = rho_threshold
    )
    if (isTRUE(keep_only_activating)) {
      adjacency <- adjacency[adjacency[["regulation"]] == 1L, , drop = FALSE]
    }
  }

  modules <- list()
  add_module <- function(tf, genes, context, regulation = 1L) {
    genes <- unique(c(tf, genes))
    if (length(genes) < min_genes) {
      return()
    }
    modules[[length(modules) + 1L]] <<- list(
      tf = tf,
      genes = genes,
      context = context,
      regulation = as.integer(regulation),
      suffix = if (as.integer(regulation) < 0L) "(-)" else "(+)"
    )
  }

  add_modules_from_adjacency <- function(adjacency_sign, regulation = 1L) {
    if (nrow(adjacency_sign) == 0) {
      return()
    }
    adjacency_sign <- adjacency_sign[order(-adjacency_sign[["importance"]]), , drop = FALSE]
    for (idx in seq_along(threshold_values)) {
      adj_thr <- adjacency_sign[adjacency_sign[["importance"]] > threshold_values[[idx]], , drop = FALSE]
      for (tf in unique(adj_thr[["TF"]])) {
        tf_adj <- adj_thr[adj_thr[["TF"]] == tf, , drop = FALSE]
        add_module(
          tf,
          tf_adj[["target"]],
          paste0("weight>", thresholds[[idx]] * 100, "%"),
          regulation = regulation
        )
      }
    }

    for (tf in unique(adjacency_sign[["TF"]])) {
      tf_adj <- adjacency_sign[adjacency_sign[["TF"]] == tf, , drop = FALSE]
      tf_adj <- tf_adj[order(-tf_adj[["importance"]]), , drop = FALSE]
      for (n_targets in top_n_targets) {
        n_targets_label <- if (is.finite(n_targets)) as.character(as.integer(n_targets)) else "all"
        n_targets_use <- if (is.finite(n_targets)) as.integer(n_targets) else nrow(tf_adj)
        add_module(
          tf,
          utils::head(tf_adj[["target"]], n_targets_use),
          paste0("top", n_targets_label),
          regulation = regulation
        )
      }
    }

    for (n in top_n_regulators) {
      top_by_target <- do.call(
        rbind,
        lapply(split(adjacency_sign, adjacency_sign[["target"]]), function(df) {
          df <- df[order(-df[["importance"]]), , drop = FALSE]
          utils::head(df, n)
        })
      )
      if (!is.null(top_by_target) && nrow(top_by_target) > 0) {
        for (tf in unique(top_by_target[["TF"]])) {
          tf_adj <- top_by_target[top_by_target[["TF"]] == tf, , drop = FALSE]
          add_module(
            tf,
            tf_adj[["target"]],
            paste0("top", n, "perTarget"),
            regulation = regulation
          )
        }
      }
    }
  }

  if (nrow(adjacency) > 0) {
    for (regulation in scenic_module_regulations(adjacency)) {
      add_modules_from_adjacency(
        scenic_filter_module_regulation(adjacency, regulation),
        regulation = regulation
      )
    }
  }

  log_message(
    "C++ pySCENIC module builder produced {.val {length(modules)}} candidate modules",
    verbose = verbose
  )
  modules
}

scenic_module_regulations <- function(adjacency) {
  if (!"regulation" %in% colnames(adjacency)) {
    return(1L)
  }
  regulations <- sort(unique(as.integer(adjacency[["regulation"]])))
  regulations[regulations %in% c(-1L, 1L)]
}

scenic_filter_module_regulation <- function(adjacency, regulation) {
  if (!"regulation" %in% colnames(adjacency)) {
    return(adjacency)
  }
  adjacency[
    as.integer(adjacency[["regulation"]]) == as.integer(regulation),
    ,
    drop = FALSE
  ]
}

scenic_add_correlation <- function(
  adjacency,
  expr_mtx,
  rho_threshold = 0.03
) {
  tfs <- intersect(unique(adjacency[["TF"]]), colnames(expr_mtx))
  targets <- intersect(unique(adjacency[["target"]]), colnames(expr_mtx))
  adjacency <- adjacency[
    adjacency[["TF"]] %in% tfs & adjacency[["target"]] %in% targets,
    ,
    drop = FALSE
  ]
  if (nrow(adjacency) == 0) {
    adjacency[["rho"]] <- numeric(0)
    adjacency[["regulation"]] <- integer(0)
    return(adjacency)
  }

  genes <- unique(c(tfs, targets))
  expr_mtx <- as.matrix(expr_mtx[, genes, drop = FALSE])
  tf_index <- match(adjacency[["TF"]], genes)
  target_index <- match(adjacency[["target"]], genes)
  rho <- tryCatch(
    scenic_edge_correlation_cpp(
      expr = expr_mtx,
      tf_index = tf_index,
      target_index = target_index
    ),
    error = function(e) NULL
  )
  if (is.null(rho)) {
    corr_mtx <- stats::cor(
      expr_mtx[, tfs, drop = FALSE],
      expr_mtx[, targets, drop = FALSE]
    )
    rho <- mapply(
      function(tf, target) corr_mtx[tf, target],
      adjacency[["TF"]],
      adjacency[["target"]],
      USE.NAMES = FALSE
    )
  }
  rho[is.na(rho)] <- 0
  adjacency[["rho"]] <- rho
  adjacency[["regulation"]] <- as.integer(rho > rho_threshold) -
    as.integer(rho < -rho_threshold)
  adjacency
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
  regulons[["regulon"]] <- scenic_regulon_tf_name(regulons[["regulon"]])
  regulons[, 1:3, drop = FALSE]
}

scenic_regulon_tf_name <- function(x) {
  sub("\\([0-9]+g\\)$", "", x)
}

scenic_regulon_name <- function(tf, suffix = "(+)") {
  tf <- sub("\\([+-]\\)$", "", as.character(tf))
  paste0(tf, suffix)
}

scenic_regulon_list <- function(regulon_tbl) {
  rgnames <- unique(regulon_tbl[["regulon"]])
  regulon_list <- lapply(rgnames, function(rg) {
    target <- regulon_tbl[regulon_tbl[["regulon"]] == rg, "target"]
    unique(unlist(strsplit(target, ",", fixed = TRUE), use.names = FALSE))
  })
  names(regulon_list) <- scenic_regulon_tf_name(rgnames)
  regulon_list
}

scenic_regulon_table <- function(regulon_list) {
  do.call(
    rbind,
    lapply(names(regulon_list), function(regulon) {
      targets <- unique(as.character(regulon_list[[regulon]]))
      targets <- targets[nzchar(targets)]
      data.frame(
        regulon = regulon,
        target = paste(targets, collapse = ","),
        target_count = length(targets),
        stringsAsFactors = FALSE
      )
    })
  )
}

scenic_read_python_aucell <- function(auc_file) {
  auc <- utils::read.csv(
    auc_file,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  auc <- as.data.frame(auc, check.names = FALSE)
  colnames(auc) <- scenic_regulon_tf_name(colnames(auc))
  auc
}

scenic_compute_aucell_score <- function(
  counts,
  regulon_list,
  min_regulon_size = 10,
  batch_size = 500,
  cores = 1,
  backend = c("r", "cpp"),
  cpp_algorithm = c("aucell", "ctxcore"),
  seed = NULL,
  verbose = TRUE
) {
  backend <- match.arg(backend)
  cpp_algorithm <- match.arg(cpp_algorithm)
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
      "Calculating AUCell regulon activity scores with {.arg backend = 'cpp'}",
      verbose = verbose
    )
    if (!is.null(seed)) {
      set.seed(as.integer(seed))
      counts <- counts[sample(rownames(counts)), , drop = FALSE]
      regulon_list <- lapply(regulon_list, intersect, rownames(counts))
    }
    scores <- run_aucell_scores(
      expr_counts = counts,
      gene_sets = regulon_list,
      strategy = "full",
      algorithm = cpp_algorithm,
      seed = if (!is.null(seed)) -1L else 0L
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
