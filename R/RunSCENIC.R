#' @title Run pySCENIC gene regulatory network analysis
#'
#' @description
#' Run pySCENIC from a Seurat object. The GRN and cisTarget steps use a
#' metacell count matrix by default, while AUCell scores are calculated on the
#' original single-cell count matrix.
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams PrepareEnv
#' @param srt A Seurat object.
#' @param layer Assay layer used as the count matrix.
#' @param ranking_dbs Character vector of cisTarget ranking feather files.
#' @param motif_annotations Motif annotation table used by `pyscenic ctx`.
#' @param tf_list Transcription-factor list used by GRNBoost2.
#' @param work_dir Directory used for pySCENIC input and output files.
#' @param prefix Prefix for pySCENIC output files.
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
#' @param min_regulon_size Minimum regulon size kept after `pyscenic ctx`.
#' @param cores Number of workers used by GRNBoost2, `pyscenic ctx`, and
#' AUCell batch scoring. If multicore execution is not supported, this is
#' automatically reduced to one core.
#' @param aucell_batch_size Number of cells scored in each AUCell batch.
#' @param seed Random seed used by GRNBoost2 and Seurat overclustering.
#' @param force Whether to rebuild existing pySCENIC outputs.
#' @param assay_name Name of the assay used to store regulon activity scores.
#' @param tool_name Name of the `srt@tools` entry.
#' @param return_seurat Whether to return the modified Seurat object. If
#' `FALSE`, a result list is returned.
#' @param envname Python environment used for pySCENIC. If `NULL`, the isolated
#' `"pyscenic_env"` environment is used.
#' @param conda The path or command name of a conda-compatible executable.
#' @param prepare_env Whether to prepare and configure the pySCENIC Python
#' environment before running.
#' @param progress Whether to show a stage-level progress bar.
#'
#' @return A Seurat object with pySCENIC results, or a result list when
#' `return_seurat = FALSE`.
#' @export
#'
#' @examples
#' \dontrun{
#' srt <- RunSCENIC(
#'   srt,
#'   ranking_dbs = c(
#'     "/path/to/hg38_500bp_up_100bp_down.genes_vs_motifs.rankings.feather",
#'     "/path/to/hg38_10kbp_up_10kbp_down.genes_vs_motifs.rankings.feather"
#'   ),
#'   motif_annotations = "/path/to/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
#'   tf_list = "/path/to/hsa_hgnc_tfs.motifs-v10.txt",
#'   work_dir = "./pyscenic",
#'   cores = 8
#' )
#' }
RunSCENIC <- function(
  srt,
  assay = NULL,
  layer = "counts",
  ranking_dbs,
  motif_annotations,
  tf_list,
  work_dir,
  prefix = "pyscenic",
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
  seed = 1234,
  force = FALSE,
  assay_name = "pyscenic",
  tool_name = "Pyscenic",
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
  if (missing(ranking_dbs) || length(ranking_dbs) == 0) {
    log_message("{.arg ranking_dbs} must be provided", message_type = "error")
  }
  if (missing(motif_annotations) || length(motif_annotations) != 1) {
    log_message("{.arg motif_annotations} must be one file", message_type = "error")
  }
  if (missing(tf_list) || length(tf_list) != 1) {
    log_message("{.arg tf_list} must be one file", message_type = "error")
  }
  if (missing(work_dir) || length(work_dir) != 1) {
    log_message("{.arg work_dir} must be one directory", message_type = "error")
  }

  ranking_dbs <- normalizePath(ranking_dbs, mustWork = TRUE)
  motif_annotations <- normalizePath(motif_annotations, mustWork = TRUE)
  tf_list <- normalizePath(tf_list, mustWork = TRUE)
  work_dir <- normalizePath(work_dir, mustWork = FALSE)
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  progress_state <- pyscenic_progress_init(progress = progress, verbose = verbose)
  on.exit(pyscenic_progress_close(progress_state), add = TRUE)
  pyscenic_progress_step(
    progress_state,
    value = 5,
    label = "Validated pySCENIC inputs",
    verbose = verbose
  )

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  cores_requested <- suppressWarnings(as.integer(cores))
  if (length(cores_requested) != 1 || is.na(cores_requested) || cores_requested < 1L) {
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

  envname <- envname %||% "pyscenic_env"
  pyscenic_python_packages <- c(
    "pyscenic==0.12.1",
    "arboreto==0.1.6",
    "ctxcore==0.2.0",
    "numpy==1.23.5",
    "dask==2024.2.1",
    "distributed==2024.2.1",
    "pyarrow"
  )
  if (isTRUE(prepare_env)) {
    pyscenic_progress_step(
      progress_state,
      value = 10,
      label = "Checking pySCENIC Python environment",
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
          packages = pyscenic_python_packages,
          envname = envname,
          conda = conda,
          verbose = FALSE
        ),
        error = function(...) {
          stats::setNames(rep(FALSE, length(pyscenic_python_packages)), pyscenic_python_packages)
        }
      )
      env_ready <- all(pkg_installed)
    }
    if (isTRUE(env_ready)) {
      log_message(
        "{.pkg pySCENIC} Python environment already has the required packages; skipping {.fn PrepareEnv}",
        verbose = verbose
      )
    } else {
      if (isTRUE(env_exists)) {
        missing_packages <- names(pkg_installed)[!pkg_installed]
        log_message(
          "{.pkg pySCENIC} Python environment is missing {.val {length(missing_packages)}} package{?s}; preparing environment",
          verbose = verbose
        )
      }
      pyscenic_progress_step(
        progress_state,
        value = 10,
        label = "Preparing pySCENIC Python environment",
        verbose = verbose
      )
      PrepareEnv(
        envname = envname,
        conda = conda,
        version = "3.10-1",
        modules = "pyscenic"
      )
    }
  }
  pyscenic_progress_step(
    progress_state,
    value = 15,
    label = "Checking pySCENIC Python packages",
    verbose = verbose
  )
  check_python(
    packages = pyscenic_python_packages,
    envname = envname,
    conda = conda,
    verbose = verbose
  )
  conda_resolved <- resolve_conda(conda)
  pyscenic_python <- conda_python(
    conda = conda_resolved,
    envname = envname
  )
  assert_python_runtime_switchable(
    pyscenic_python,
    restart_hint = pyscenic_runtime_restart_hint(envname = envname)
  )
  configure_python_runtime(pyscenic_python)

  functions <- reticulate::import_from_path(
    "functions",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )

  log_message(
    "Preparing {.pkg pySCENIC} input matrix",
    verbose = verbose
  )
  pyscenic_progress_step(
    progress_state,
    value = 20,
    label = "Loading expression counts",
    verbose = verbose
  )
  counts <- GetAssayData5(srt, assay = assay, layer = layer)
  counts <- counts[, colnames(srt), drop = FALSE]
  if (nrow(counts) == 0 || ncol(counts) == 0) {
    log_message(
      "No expression values available for {.pkg pySCENIC}",
      message_type = "error"
    )
  }
  log_message(
    "{.pkg pySCENIC} expression input: assay {.val {assay}}, layer {.val {layer}}, {.val {nrow(counts)}} genes x {.val {ncol(counts)}} cells",
    verbose = verbose
  )

  if (isTRUE(metacell)) {
    pyscenic_progress_step(
      progress_state,
      value = 30,
      label = "Building metacell counts for GRNBoost2",
      verbose = verbose
    )
    metacell_result <- pyscenic_build_metacell_counts(
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
    pyscenic_progress_step(
      progress_state,
      value = 30,
      label = "Using single-cell counts for GRNBoost2",
      verbose = verbose
    )
  }
  grn_count_source <- if (isTRUE(metacell)) "metacell matrix" else "single-cell matrix"
  grn_count_unit <- if (isTRUE(metacell)) "cells/metacells" else "cells"
  log_message(
    "{.pkg pySCENIC} GRN count source: {.val {grn_count_source}}, {.val {nrow(grn_counts)}} genes x {.val {ncol(grn_counts)}} {grn_count_unit}",
    verbose = verbose
  )

  pyscenic_progress_step(
    progress_state,
    value = 40,
    label = "Reading cisTarget ranking database genes",
    verbose = verbose
  )
  ranking_gene_lists <- lapply(ranking_dbs, function(db) {
    as.character(unlist(functions$PyscenicRankingGenes(db), use.names = FALSE))
  })
  ranking_genes <- Reduce(intersect, ranking_gene_lists)
  grn_matrix <- pyscenic_prepare_grn_matrix(
    counts = grn_counts,
    ranking_genes = ranking_genes,
    min_expr_cells = min_expr_cells,
    verbose = verbose
  )

  pyscenic_progress_step(
    progress_state,
    value = 48,
    label = "Writing GRNBoost2 expression matrix",
    verbose = verbose
  )
  expr_csv <- file.path(work_dir, paste0(prefix, "_expression_mtx.csv"))
  if (!file.exists(expr_csv) || isTRUE(force)) {
    utils::write.csv(
      as.matrix(grn_matrix),
      file = expr_csv,
      row.names = TRUE
    )
  } else {
    log_message(
      "Reusing existing {.pkg pySCENIC} expression matrix: {.file {expr_csv}}",
      verbose = verbose
    )
  }

  adj_file <- file.path(work_dir, paste0(prefix, "_step1_adj.tsv"))
  ctx_file <- file.path(work_dir, paste0(prefix, "_step2_reg.tsv"))
  gmt_file <- file.path(work_dir, paste0(prefix, "_step2_regulons.gmt"))
  txt_file <- file.path(work_dir, paste0(prefix, "_step2_regulons.txt"))
  ras_file <- file.path(work_dir, paste0(prefix, "_regulon_activity_score.rds"))
  regulon_list_file <- file.path(work_dir, paste0(prefix, "_regulon_list.rds"))

  pyscenic_progress_step(
    progress_state,
    value = 55,
    label = "Running GRNBoost2",
    verbose = verbose
  )
  functions$RunSCENICGrn(
    expression_mtx = expr_csv,
    tf_list = tf_list,
    adj_output = adj_file,
    cores = as.integer(cores),
    seed = as.integer(seed),
    force = isTRUE(force),
    verbose = isTRUE(verbose)
  )

  pyscenic_progress_step(
    progress_state,
    value = 70,
    label = "Running pySCENIC cisTarget pruning",
    verbose = verbose
  )
  functions$RunSCENICCtx(
    expression_mtx = expr_csv,
    ranking_dbs = as.list(ranking_dbs),
    motif_annotations = motif_annotations,
    adj_output = adj_file,
    ctx_output = ctx_file,
    cores = as.integer(cores),
    force = isTRUE(force),
    verbose = isTRUE(verbose)
  )

  pyscenic_progress_step(
    progress_state,
    value = 82,
    label = "Converting pySCENIC regulons",
    verbose = verbose
  )
  if (isTRUE(force) || !file.exists(gmt_file) || !file.exists(txt_file)) {
    functions$PyscenicRegulonsToFiles(
      regulon_file = ctx_file,
      gmt_file = gmt_file,
      txt_file = txt_file,
      min_regulon_size = as.integer(min_regulon_size)
    )
  } else {
    log_message(
      "Reusing existing pySCENIC regulon files",
      verbose = verbose
    )
  }

  pyscenic_progress_step(
    progress_state,
    value = 88,
    label = "Reading regulon target lists",
    verbose = verbose
  )
  regulon_tbl <- pyscenic_read_regulon_txt(txt_file)
  regulon_list <- pyscenic_regulon_list(regulon_tbl)
  if (!file.exists(regulon_list_file) || isTRUE(force)) {
    saveRDS(regulon_list, regulon_list_file)
  }

  if (file.exists(ras_file) && isFALSE(force)) {
    log_message(
      "Reusing existing regulon activity matrix: {.file {ras_file}}",
      verbose = verbose
    )
    ras_mat <- readRDS(ras_file)
  } else {
    pyscenic_progress_step(
      progress_state,
      value = 92,
      label = "Calculating AUCell regulon activity scores",
      verbose = verbose
    )
    ras_mat <- pyscenic_compute_aucell_score(
      counts = counts,
      regulon_list = regulon_list,
      min_regulon_size = min_regulon_size,
      batch_size = aucell_batch_size,
      cores = cores,
      verbose = verbose
    )
    saveRDS(ras_mat, ras_file)
  }
  ras_mat <- ras_mat[colnames(srt), , drop = FALSE]

  pyscenic_progress_step(
    progress_state,
    value = 96,
    label = "Collecting pySCENIC result files",
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
      ranking_dbs = ranking_dbs,
      motif_annotations = motif_annotations,
      tf_list = tf_list,
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
      seed = seed,
      assay_name = assay_name,
      tool_name = tool_name,
      envname = envname,
      progress = progress
    )
  )

  if (isFALSE(return_seurat)) {
    pyscenic_progress_step(
      progress_state,
      value = 100,
      label = "Finished pySCENIC",
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
  pyscenic_progress_step(
    progress_state,
    value = 100,
    label = "Stored pySCENIC results",
    verbose = verbose
  )
  log_message(
    "{.pkg pySCENIC} results stored in assay {.val {assay_name}} and tools slot {.val {tool_name}}",
    verbose = verbose
  )
  srt
}

pyscenic_progress_init <- function(progress = TRUE, verbose = TRUE) {
  if (!isTRUE(progress) || !isTRUE(verbose)) {
    return(NULL)
  }
  progress_state <- new.env(parent = emptyenv())
  progress_state$pb <- utils::txtProgressBar(min = 0, max = 100, initial = 0, style = 3)
  progress_state$value <- 0L
  progress_state
}

pyscenic_progress_step <- function(progress_state, value, label, verbose = TRUE) {
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

pyscenic_progress_close <- function(progress_state) {
  if (!is.null(progress_state)) {
    close(progress_state$pb)
  }
  invisible(NULL)
}

pyscenic_default_target_metacells <- function(n_cells) {
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

pyscenic_metacell_labels <- function(cluster_vec, group_df = NULL) {
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

pyscenic_resolution_summary <- function(
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
    metacell_labels <- pyscenic_metacell_labels(cluster_vec, group_df = group_df)
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

pyscenic_select_metacell_resolution <- function(resolution_summary, target_metacells) {
  if (is.null(resolution_summary) || nrow(resolution_summary) == 0) {
    log_message(
      "Cannot summarize metacell candidates from {.fn Seurat::FindClusters}",
      message_type = "error"
    )
  }

  resolution_summary[["target_distance"]] <- abs(
    resolution_summary[["n_metacells"]] - target_metacells
  )
  resolution_summary <- resolution_summary[order(
    resolution_summary[["target_distance"]],
    -resolution_summary[["n_metacells"]],
    resolution_summary[["resolution"]]
  ), , drop = FALSE]

  list(
    selected_resolution = resolution_summary[["resolution"]][[1]],
    selected_cluster_col = resolution_summary[["cluster_col"]][[1]],
    resolution_summary = resolution_summary
  )
}

pyscenic_build_metacell_counts <- function(
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
    "Building metacells for {.pkg pySCENIC} GRN input",
    verbose = verbose
  )
  metacell_reduction <- as.character(metacell_reduction)
  if (length(metacell_reduction) != 1 || is.na(metacell_reduction) || !nzchar(metacell_reduction)) {
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
        "At least two cells and one variable feature are required to build pySCENIC metacells with PCA",
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
    reduction_dims <- ncol(Seurat::Embeddings(srt, reduction = metacell_reduction))
    if (max(metacell_dims) > reduction_dims) {
      log_message(
        "{.arg metacell_dims} requests dimension {.val {max(metacell_dims)}}, but {.arg metacell_reduction} {.val {metacell_reduction}} has only {.val {reduction_dims}} dimensions",
        message_type = "error"
      )
    }
    log_message(
      "Using existing reduction {.val {metacell_reduction}} for {.pkg pySCENIC} metacell overclustering",
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
        "{.arg metacell.by} contains missing values in {.val {sum(missing_group)}} cell{?s} across column{?s} {.val {missing_group_cols}}; fill or remove missing annotations before building pySCENIC metacells.",
        message_type = "error"
      )
    }
  }

  auto_resolution <- is.null(metacell_resolution)
  if (isTRUE(auto_resolution)) {
    metacell_target <- metacell_target %||% pyscenic_default_target_metacells(ncol(srt))
    metacell_target <- max(1L, as.integer(metacell_target))
    metacell_resolution_candidates <- unique(as.numeric(metacell_resolution_candidates))
    metacell_resolution_candidates <- metacell_resolution_candidates[
      !is.na(metacell_resolution_candidates) & metacell_resolution_candidates > 0
    ]
    if (length(metacell_resolution_candidates) == 0) {
      log_message(
        "{.arg metacell_resolution_candidates} must contain at least one positive number",
        message_type = "error"
      )
    }
    log_message(
      "Auto {.pkg pySCENIC} target metacells: {.val {metacell_target}} for {.val {ncol(srt)}} cells",
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
  cluster_cols <- cluster_cols[startsWith(cluster_cols, paste0(assay, "_snn_res."))]
  cluster_resolutions <- suppressWarnings(
    as.numeric(sub("^.*_snn_res\\.", "", cluster_cols))
  )
  cluster_cols <- cluster_cols[
    vapply(cluster_resolutions, function(x) {
      !is.na(x) && any(abs(x - metacell_resolution_candidates) < sqrt(.Machine$double.eps))
    }, logical(1))
  ]
  if (length(cluster_cols) == 0) {
    log_message(
      "Cannot find metacell overclustering column after {.fn Seurat::FindClusters}",
      message_type = "error"
    )
  }
  resolution_summary <- pyscenic_resolution_summary(
    srt = srt,
    cluster_cols = cluster_cols,
    resolution_candidates = metacell_resolution_candidates,
    group_df = group_df
  )

  if (isTRUE(auto_resolution)) {
    selection <- pyscenic_select_metacell_resolution(
      resolution_summary = resolution_summary,
      target_metacells = metacell_target
    )
    resolution_summary <- selection[["resolution_summary"]]
    cluster_col <- selection[["selected_cluster_col"]]
    selected_resolution <- selection[["selected_resolution"]]
    metacell_scan <- paste(
      sprintf(
        "%s:%s",
        format(resolution_summary[["resolution"]], trim = TRUE, scientific = FALSE),
        resolution_summary[["n_metacells"]]
      ),
      collapse = ", "
    )
    log_message(
      "{.pkg pySCENIC} metacell resolution scan (resolution:n_metacells): {metacell_scan}",
      verbose = verbose
    )
  } else {
    cluster_col <- resolution_summary[["cluster_col"]][[1]]
    selected_resolution <- resolution_summary[["resolution"]][[1]]
  }

  cluster_vec <- as.character(srt[[cluster_col]][, 1])
  metacell_labels <- pyscenic_metacell_labels(cluster_vec, group_df = group_df)

  metacells <- split(colnames(srt), metacell_labels)
  metacell_names <- paste0("metacell_", seq_along(metacells))
  original_labels <- names(metacells)
  names(metacells) <- metacell_names
  metacell_sizes <- lengths(metacells)
  log_message(
    "{.pkg pySCENIC} selected metacell resolution {.val {selected_resolution}} with {.val {length(metacells)}} metacells",
    verbose = verbose
  )
  log_message(
    "{.pkg pySCENIC} metacell size summary: min {.val {min(metacell_sizes)}}, median {.val {stats::median(metacell_sizes)}}, mean {.val {round(mean(metacell_sizes), 2)}}, max {.val {max(metacell_sizes)}} cells",
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

pyscenic_prepare_grn_matrix <- function(
  counts,
  ranking_genes,
  min_expr_cells = 3,
  verbose = TRUE
) {
  grn_matrix <- Matrix::t(counts)
  expr_in_cells <- Matrix::colSums(grn_matrix > 0)
  grn_matrix <- grn_matrix[, expr_in_cells >= min_expr_cells, drop = FALSE]
  grn_matrix <- grn_matrix[, colnames(grn_matrix) %in% ranking_genes, drop = FALSE]
  if (ncol(grn_matrix) == 0) {
    log_message(
      "No genes remain after filtering by expression and cisTarget database",
      message_type = "error"
    )
  }
  log_message(
    "{.pkg pySCENIC} GRN input matrix: {.val {nrow(grn_matrix)}} cells/metacells x {.val {ncol(grn_matrix)}} genes",
    verbose = verbose
  )
  grn_matrix
}

pyscenic_read_regulon_txt <- function(txt_file) {
  if (!file.exists(txt_file)) {
    log_message(
      "Cannot find pySCENIC regulon table: {.file {txt_file}}",
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
      "pySCENIC regulon table must contain regulon, motif, and target columns",
      message_type = "error"
    )
  }
  colnames(regulons)[1:3] <- c("regulon", "motif", "target")
  regulons[, 1:3, drop = FALSE]
}

pyscenic_regulon_list <- function(regulon_tbl) {
  rgnames <- unique(regulon_tbl[["regulon"]])
  regulon_list <- lapply(rgnames, function(rg) {
    target <- regulon_tbl[regulon_tbl[["regulon"]] == rg, "target"]
    unique(unlist(strsplit(target, ",", fixed = TRUE), use.names = FALSE))
  })
  names(regulon_list) <- sub("[0-9]+g", "+", rgnames)
  regulon_list
}

pyscenic_compute_aucell_score <- function(
  counts,
  regulon_list,
  min_regulon_size = 10,
  batch_size = 500,
  cores = 1,
  verbose = TRUE
) {
  check_r("AUCell", verbose = FALSE)
  regulon_list <- lapply(regulon_list, intersect, rownames(counts))
  regulon_list <- regulon_list[lengths(regulon_list) >= min_regulon_size]
  if (length(regulon_list) == 0) {
    log_message(
      "No regulons remain after matching targets to expression matrix",
      message_type = "error"
    )
  }

  batches <- split(
    seq_len(ncol(counts)),
    ceiling(seq_along(colnames(counts)) / batch_size)
  )
  log_message(
    "Calculating AUCell regulon activity scores with {.val {length(batches)}} batch{?es} and {.val {cores}} core{?s}",
    verbose = verbose
  )
  calc_auc <- function(cells_idx) {
    pyscenic_calc_auc_batch(
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

pyscenic_calc_auc_batch <- function(counts, regulon_list) {
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
