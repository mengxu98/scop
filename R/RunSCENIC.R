#' @title Run SCENIC gene regulatory network analysis
#'
#' @description
#' Run SCENIC from a Seurat object. The GRN and cisTarget steps use a
#' metacell count matrix by default, while AUCell scores are calculated on the
#' original single-cell count matrix.
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
#' values include `"Homo_sapiens"`, `"Mus_musculus"`,
#' `"Drosophila_melanogaster"` and aliases such as `"human"`, `"mouse"`, and
#' `"fly"`.
#' @param data_dir Directory used to cache automatically prepared SCENIC
#' reference files. If `NULL`, files are stored under
#' `tools::R_user_dir("scop", "data")/SCENIC/<species>`.
#' @param prefix Prefix for SCENIC output files.
#' @param metacell Whether to build a metacell count matrix for GRNBoost2.
#' @param metacell.by Optional metadata column(s) used to keep metacells within
#' groups such as samples or cell types.
#' @param metacell_resolution Resolution passed to [Seurat::FindClusters()] for
#' overclustering. If `NULL`, candidate resolutions are scanned and the one
#' closest to `metacell_target` is used.
#' @param metacell_target Target number of metacells used when
#' `metacell_resolution = NULL`. If `NULL`, a default target is chosen from the
#' number of cells.
#' @param metacell_resolution_candidates Candidate resolutions scanned when
#' `metacell_resolution = NULL`.
#' @param metacell_reduction Reduction used to build the metacell neighbor graph.
#' The default `"pca"` keeps the original behavior and recomputes PCA from the
#' selected assay. To use an already batch-corrected embedding such as Harmony,
#' run it before `RunSCENIC()` and set `metacell_reduction = "Harmony"`.
#' Only the metacell grouping uses this reduction; GRNBoost2 still uses raw
#' count sums per metacell.
#' @param metacell_dims Dimensions used for metacell overclustering.
#' @param min_expr_cells Minimum number of cells or metacells where a gene must
#' be detected before GRNBoost2.
#' @param min_regulon_size Minimum regulon size kept after `scenic ctx`.
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
#' @param progress Whether to show a stage-level progress bar.
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
  data_dir = NULL,
  prefix = "scenic",
  metacell = TRUE,
  metacell.by = NULL,
  metacell_resolution = NULL,
  metacell_target = NULL,
  metacell_resolution_candidates = c(0.5, 1, 2, 5, 10, 20, 30, 40, 50, 75, 100),
  metacell_reduction = "pca",
  metacell_dims = 1:30,
  min_expr_cells = 3,
  min_regulon_size = 10,
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
  progress = verbose,
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

  reference_data <- scenic_resolve_reference_data(
    species = species,
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
  regulators_info <- scenic_prepare_gene_list_argument(
    x = regulators,
    arg = "regulators",
    out_file = file.path(work_dir, paste0(prefix, "_regulators.txt")),
    required = TRUE
  )
  regulators <- regulators_info[["genes"]]
  regulators_file <- regulators_info[["file"]]
  targets <- scenic_prepare_gene_list_argument(
    x = targets,
    arg = "targets",
    out_file = file.path(work_dir, paste0(prefix, "_targets.txt")),
    required = FALSE
  )[["genes"]]
  progress_state <- scenic_progress_init(progress = progress, verbose = verbose)
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
  aucell_backend <- match.arg(aucell_backend)
  aucell_cpp_strategy <- match.arg(aucell_cpp_strategy)

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

  if (isTRUE(metacell)) {
    scenic_progress_step(
      progress_state,
      value = 30,
      label = "Building metacell counts for GRNBoost2",
      verbose = verbose
    )
    metacell_result <- scenic_build_metacell_counts(
      srt = srt,
      counts = counts,
      assay = assay,
      metacell.by = metacell.by,
      metacell_resolution = metacell_resolution,
      metacell_target = metacell_target,
      metacell_resolution_candidates = metacell_resolution_candidates,
      metacell_reduction = metacell_reduction,
      metacell_dims = metacell_dims,
      seed = seed,
      verbose = verbose
    )
    grn_counts <- metacell_result[["counts"]]
    metacell_info <- metacell_result[["info"]]
  } else {
    grn_counts <- counts
    metacell_info <- NULL
    scenic_progress_step(
      progress_state,
      value = 30,
      label = "Using single-cell counts for GRNBoost2",
      verbose = verbose
    )
  }
  grn_count_source <- if (isTRUE(metacell)) {
    "metacell matrix"
  } else {
    "single-cell matrix"
  }
  grn_count_unit <- if (isTRUE(metacell)) "cells/metacells" else "cells"
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
    scenic_filter_adjacency_targets(
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
    metacell = metacell_info,
    parameters = list(
      assay = assay,
      layer = layer,
      species = species,
      data_dir = data_dir,
      ranking_dbs = ranking_dbs,
      motif_annotations = motif_annotations,
      regulators = regulators,
      regulators_file = regulators_file,
      targets = targets,
      metacell = metacell,
      metacell.by = metacell.by,
      metacell_resolution = metacell_resolution,
      metacell_target = metacell_target,
      metacell_resolution_candidates = metacell_resolution_candidates,
      metacell_reduction = metacell_reduction,
      metacell_dims = metacell_dims,
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
      progress = progress
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

scenic_resolve_reference_data <- function(
  species,
  data_dir = NULL,
  ranking_dbs = NULL,
  motif_annotations = NULL,
  regulators = NULL,
  verbose = TRUE
) {
  species_config <- scenic_species_config(species)
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
      ranking_dbs <- scenic_download_reference_files(
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
      motif_annotations <- scenic_download_reference_files(
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
      regulators <- scenic_download_reference_files(
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
    data_dir = reference_dir,
    ranking_dbs = ranking_dbs,
    motif_annotations = motif_annotations,
    regulators = regulators
  )
}

scenic_species_config <- function(species) {
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
  species_key <- switch(
    species_key,
    homo_sapiens = "human",
    human = "human",
    hsa = "human",
    hs = "human",
    hg38 = "human",
    mus_musculus = "mouse",
    mouse = "mouse",
    mm = "mouse",
    mm10 = "mouse",
    drosophila_melanogaster = "fly",
    drosophila = "fly",
    fly = "fly",
    dmel = "fly",
    dm6 = "fly",
    NULL
  )
  if (is.null(species_key)) {
    log_message(
      "{.arg species} must be one of {.val Homo_sapiens}, {.val Mus_musculus}, or {.val Drosophila_melanogaster}",
      message_type = "error"
    )
  }

  cistarget_url <- "https://resources.aertslab.org/cistarget"
  switch(
    species_key,
    human = list(
      key = "human",
      label = "Homo_sapiens",
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
    ),
    mouse = list(
      key = "mouse",
      label = "Mus_musculus",
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

scenic_download_reference_files <- function(
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

scenic_prepare_gene_list_argument <- function(
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

scenic_filter_adjacency_targets <- function(
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

scenic_progress_init <- function(progress = TRUE, verbose = TRUE) {
  if (!isTRUE(progress) || !isTRUE(verbose)) {
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
    cat("\n")
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

scenic_default_target_metacells <- function(n_cells) {
  n_cells <- as.integer(n_cells)
  if (length(n_cells) != 1 || is.na(n_cells) || n_cells <= 0) {
    log_message(
      "{.arg n_cells} must be a positive integer",
      message_type = "error"
    )
  }

  target <- if (n_cells <= 500) {
    min(n_cells, 80L)
  } else if (n_cells <= 1000) {
    100L
  } else if (n_cells <= 2000) {
    150L
  } else if (n_cells <= 5000) {
    300L
  } else if (n_cells <= 10000) {
    500L
  } else if (n_cells <= 20000) {
    700L
  } else if (n_cells <= 50000) {
    1000L
  } else if (n_cells <= 100000) {
    1500L
  } else if (n_cells <= 200000) {
    2000L
  } else {
    min(round(n_cells / 100), 10000L)
  }

  as.integer(target)
}

scenic_metacell_labels <- function(cluster_vec, group_df = NULL) {
  cluster_vec <- as.character(cluster_vec)
  if (!is.null(group_df)) {
    group_df[] <- lapply(group_df, as.character)
    return(do.call(
      interaction,
      c(group_df, list(cluster = cluster_vec, drop = TRUE, sep = "_"))
    ))
  }

  factor(cluster_vec)
}

scenic_resolution_summary <- function(
  srt,
  cluster_cols,
  resolution_candidates,
  group_df = NULL
) {
  cluster_resolutions <- suppressWarnings(
    as.numeric(sub("^.*_snn_res\\.", "", cluster_cols))
  )

  summary_list <- lapply(resolution_candidates, function(resolution) {
    resolution_idx <- which(
      !is.na(cluster_resolutions) &
        abs(cluster_resolutions - resolution) < sqrt(.Machine$double.eps)
    )
    if (length(resolution_idx) == 0) {
      return(NULL)
    }

    cluster_col <- cluster_cols[[resolution_idx[[length(resolution_idx)]]]]
    cluster_vec <- as.character(srt[[cluster_col]][, 1])
    metacell_labels <- scenic_metacell_labels(cluster_vec, group_df = group_df)
    cell_counts <- lengths(split(colnames(srt), metacell_labels))

    data.frame(
      resolution = resolution,
      cluster_col = cluster_col,
      n_metacells = length(cell_counts),
      min_cells = min(cell_counts),
      median_cells = stats::median(cell_counts),
      mean_cells = mean(cell_counts),
      max_cells = max(cell_counts),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, summary_list)
}

scenic_select_metacell_resolution <- function(
  resolution_summary,
  target_metacells
) {
  if (is.null(resolution_summary) || nrow(resolution_summary) == 0) {
    log_message(
      "Cannot summarize metacell candidates from {.fn Seurat::FindClusters}",
      message_type = "error"
    )
  }

  resolution_summary[["target_distance"]] <- abs(
    resolution_summary[["n_metacells"]] - target_metacells
  )
  resolution_summary <- resolution_summary[
    order(
      resolution_summary[["target_distance"]],
      -resolution_summary[["n_metacells"]],
      resolution_summary[["resolution"]]
    ),
    ,
    drop = FALSE
  ]

  list(
    selected_resolution = resolution_summary[["resolution"]][[1]],
    selected_cluster_col = resolution_summary[["cluster_col"]][[1]],
    resolution_summary = resolution_summary
  )
}

scenic_build_metacell_counts <- function(
  srt,
  counts,
  assay,
  metacell.by = NULL,
  metacell_resolution = NULL,
  metacell_target = NULL,
  metacell_resolution_candidates = c(0.5, 1, 2, 5, 10, 20, 30, 40, 50, 75, 100),
  metacell_reduction = "pca",
  metacell_dims = 1:30,
  seed = 1234,
  verbose = TRUE
) {
  set.seed(seed)
  log_message(
    "Building metacells for {.pkg SCENIC} GRN input",
    verbose = verbose
  )
  metacell_reduction <- as.character(metacell_reduction)
  if (
    length(metacell_reduction) != 1 ||
      is.na(metacell_reduction) ||
      !nzchar(metacell_reduction)
  ) {
    log_message(
      "{.arg metacell_reduction} must be one reduction name",
      message_type = "error"
    )
  }
  metacell_dims <- unique(as.integer(metacell_dims))
  metacell_dims <- metacell_dims[!is.na(metacell_dims) & metacell_dims > 0L]
  if (length(metacell_dims) == 0) {
    log_message(
      "{.arg metacell_dims} must contain at least one positive integer",
      message_type = "error"
    )
  }
  if (identical(metacell_reduction, "pca")) {
    srt <- Seurat::NormalizeData(srt, assay = assay, verbose = FALSE)
    srt <- Seurat::FindVariableFeatures(
      srt,
      assay = assay,
      selection.method = "vst",
      nfeatures = 2000,
      verbose = FALSE
    )
    srt <- Seurat::ScaleData(srt, assay = assay, verbose = FALSE)
    pca_features <- SeuratObject::VariableFeatures(srt, assay = assay)
    max_pca_dims <- min(length(pca_features), ncol(srt) - 1L)
    if (max_pca_dims < 1L) {
      log_message(
        "At least two cells and one variable feature are required to build SCENIC metacells with PCA",
        message_type = "error"
      )
    }
    if (max(metacell_dims) > max_pca_dims) {
      log_message(
        "{.arg metacell_dims} requests dimension {.val {max(metacell_dims)}}, but PCA can use at most {.val {max_pca_dims}} dimensions for this object; using available dimensions only",
        message_type = "warning",
        verbose = verbose
      )
      metacell_dims <- metacell_dims[metacell_dims <= max_pca_dims]
    }
    if (length(metacell_dims) == 0) {
      log_message(
        "No valid {.arg metacell_dims} remain after checking PCA limits",
        message_type = "error"
      )
    }
    srt <- Seurat::RunPCA(
      srt,
      assay = assay,
      npcs = max(metacell_dims),
      verbose = FALSE
    )
  } else {
    available_reductions <- SeuratObject::Reductions(srt)
    if (!metacell_reduction %in% available_reductions) {
      available_reductions_text <- if (length(available_reductions) > 0) {
        paste(available_reductions, collapse = ", ")
      } else {
        "none"
      }
      log_message(
        "{.arg metacell_reduction} {.val {metacell_reduction}} is not present in {.arg srt}. Run the batch-corrected reduction before {.fn RunSCENIC} or use {.val {'pca'}}. Available reductions: {.val {available_reductions_text}}",
        message_type = "error"
      )
    }
    reduction_dims <- ncol(Seurat::Embeddings(
      srt,
      reduction = metacell_reduction
    ))
    if (max(metacell_dims) > reduction_dims) {
      log_message(
        "{.arg metacell_dims} requests dimension {.val {max(metacell_dims)}}, but {.arg metacell_reduction} {.val {metacell_reduction}} has only {.val {reduction_dims}} dimensions",
        message_type = "error"
      )
    }
    log_message(
      "Using existing reduction {.val {metacell_reduction}} for {.pkg SCENIC} metacell overclustering",
      verbose = verbose
    )
  }
  srt <- Seurat::FindNeighbors(
    srt,
    reduction = metacell_reduction,
    dims = metacell_dims,
    k.param = 20,
    verbose = FALSE
  )

  group_df <- NULL
  if (!is.null(metacell.by)) {
    missing_cols <- setdiff(metacell.by, colnames(srt@meta.data))
    if (length(missing_cols) > 0) {
      log_message(
        "{.arg metacell.by} columns not found in {.arg srt}: {.val {missing_cols}}",
        message_type = "error"
      )
    }
    group_df <- srt@meta.data[, metacell.by, drop = FALSE]
    missing_group <- !stats::complete.cases(group_df)
    if (any(missing_group)) {
      missing_group_cols <- names(group_df)[
        vapply(group_df, function(x) any(is.na(x)), logical(1))
      ]
      log_message(
        "{.arg metacell.by} contains missing values in {.val {sum(missing_group)}} cell{?s} across column{?s} {.val {missing_group_cols}}; fill or remove missing annotations before building SCENIC metacells.",
        message_type = "error"
      )
    }
  }

  auto_resolution <- is.null(metacell_resolution)
  if (isTRUE(auto_resolution)) {
    metacell_target <- metacell_target %||%
      scenic_default_target_metacells(ncol(srt))
    metacell_target <- max(1L, as.integer(metacell_target))
    metacell_resolution_candidates <- unique(as.numeric(
      metacell_resolution_candidates
    ))
    metacell_resolution_candidates <- metacell_resolution_candidates[
      !is.na(metacell_resolution_candidates) &
        metacell_resolution_candidates > 0
    ]
    if (length(metacell_resolution_candidates) == 0) {
      log_message(
        "{.arg metacell_resolution_candidates} must contain at least one positive number",
        message_type = "error"
      )
    }
    log_message(
      "Auto {.pkg SCENIC} target metacells: {.val {metacell_target}} for {.val {ncol(srt)}} cells",
      verbose = verbose
    )
  } else {
    metacell_target <- NULL
    metacell_resolution_candidates <- as.numeric(metacell_resolution)
  }

  srt <- Seurat::FindClusters(
    srt,
    resolution = metacell_resolution_candidates,
    random.seed = seed,
    verbose = FALSE
  )
  cluster_cols <- colnames(srt@meta.data)
  cluster_cols <- cluster_cols[startsWith(
    cluster_cols,
    paste0(assay, "_snn_res.")
  )]
  cluster_resolutions <- suppressWarnings(
    as.numeric(sub("^.*_snn_res\\.", "", cluster_cols))
  )
  cluster_cols <- cluster_cols[
    vapply(
      cluster_resolutions,
      function(x) {
        !is.na(x) &&
          any(
            abs(x - metacell_resolution_candidates) < sqrt(.Machine$double.eps)
          )
      },
      logical(1)
    )
  ]
  if (length(cluster_cols) == 0) {
    log_message(
      "Cannot find metacell overclustering column after {.fn Seurat::FindClusters}",
      message_type = "error"
    )
  }
  resolution_summary <- scenic_resolution_summary(
    srt = srt,
    cluster_cols = cluster_cols,
    resolution_candidates = metacell_resolution_candidates,
    group_df = group_df
  )

  if (isTRUE(auto_resolution)) {
    selection <- scenic_select_metacell_resolution(
      resolution_summary = resolution_summary,
      target_metacells = metacell_target
    )
    resolution_summary <- selection[["resolution_summary"]]
    cluster_col <- selection[["selected_cluster_col"]]
    selected_resolution <- selection[["selected_resolution"]]
    metacell_scan <- paste(
      sprintf(
        "%s:%s",
        format(
          resolution_summary[["resolution"]],
          trim = TRUE,
          scientific = FALSE
        ),
        resolution_summary[["n_metacells"]]
      ),
      collapse = ", "
    )
    log_message(
      "{.pkg SCENIC} metacell resolution scan (resolution:n_metacells): {metacell_scan}",
      verbose = verbose
    )
  } else {
    cluster_col <- resolution_summary[["cluster_col"]][[1]]
    selected_resolution <- resolution_summary[["resolution"]][[1]]
  }

  cluster_vec <- as.character(srt[[cluster_col]][, 1])
  metacell_labels <- scenic_metacell_labels(cluster_vec, group_df = group_df)

  metacells <- split(colnames(srt), metacell_labels)
  metacell_names <- paste0("metacell_", seq_along(metacells))
  original_labels <- names(metacells)
  names(metacells) <- metacell_names
  metacell_sizes <- lengths(metacells)
  log_message(
    "{.pkg SCENIC} selected metacell resolution {.val {selected_resolution}} with {.val {length(metacells)}} metacells",
    verbose = verbose
  )
  log_message(
    "{.pkg SCENIC} metacell size summary: min {.val {min(metacell_sizes)}}, median {.val {stats::median(metacell_sizes)}}, mean {.val {round(mean(metacell_sizes), 2)}}, max {.val {max(metacell_sizes)}} cells",
    verbose = verbose
  )

  meta_counts <- do.call(
    cbind,
    lapply(metacells, function(cells) {
      Matrix::rowSums(counts[, cells, drop = FALSE])
    })
  )
  meta_counts <- Matrix::Matrix(meta_counts, sparse = TRUE)
  rownames(meta_counts) <- rownames(counts)
  colnames(meta_counts) <- metacell_names

  cell_map <- data.frame(
    cell = unlist(metacells, use.names = FALSE),
    metacell = rep(names(metacells), lengths(metacells)),
    metacell_label = rep(original_labels, lengths(metacells)),
    stringsAsFactors = FALSE
  )
  info <- list(
    cluster_col = cluster_col,
    resolution = selected_resolution,
    reduction = metacell_reduction,
    auto_resolution = auto_resolution,
    target_metacells = metacell_target,
    selected_resolution = selected_resolution,
    resolution_summary = resolution_summary,
    dims = metacell_dims,
    n_metacells = length(metacells),
    cell_counts = data.frame(
      metacell = names(metacells),
      metacell_label = original_labels,
      n_cells = lengths(metacells),
      stringsAsFactors = FALSE
    ),
    cell_map = cell_map
  )

  list(counts = meta_counts, info = info)
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
    "{.pkg SCENIC} GRN input matrix: {.val {nrow(grn_matrix)}} cells/metacells x {.val {ncol(grn_matrix)}} genes",
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

  if (cores > 1 && .Platform$OS.type == "unix") {
    scores_list <- parallel::mclapply(
      batches,
      calc_auc,
      mc.cores = cores
    )
  } else {
    if (cores > 1) {
      log_message(
        "{.arg cores} > 1 requires a Unix-like OS for AUCell batch scoring; using one core",
        message_type = "warning",
        verbose = verbose
      )
    }
    scores_list <- lapply(batches, calc_auc)
  }

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
