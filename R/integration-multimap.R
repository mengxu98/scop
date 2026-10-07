#' @title MultiMAP integration function
#'
#' @inheritParams RunIntegration
#' @param gene_activity_assay Name of the gene activity assay used to provide
#' a shared feature space for RNA-ATAC integration.
#' @param MultiMAP_params A list of parameters passed to `MultiMAP::Integration`.
#' The following keys are managed internally and should not be supplied:
#' `"adatas"`, `"use_reps"`, `"embedding"`, and `"seed"`.
#'
#' @export
#' @examples
#' \dontrun{
#' data("pbmcmultiome_sub", package = "scop")
#' pbmcmultiome_sub$batch <- rep(c("batch1", "batch2"), length.out = ncol(pbmcmultiome_sub))
#' pbmcmultiome_sub <- MultiMAP_integrate(
#'   srt_merge = pbmcmultiome_sub,
#'   batch = "batch",
#'   linear_reduction_dims = 20,
#'   linear_reduction_dims_use = 1:10
#' )
#' }
MultiMAP_integrate <- function(
  srt_merge = NULL,
  batch = NULL,
  append = TRUE,
  srt_list = NULL,
  assay = NULL,
  do_normalization = NULL,
  normalization_method = "LogNormalize",
  do_HVF_finding = TRUE,
  HVF_source = "separate",
  HVF_method = "vst",
  nHVF = 2000,
  HVF_min_intersection = 1,
  HVF = NULL,
  do_scaling = TRUE,
  vars_to_regress = NULL,
  regression_model = "linear",
  scale_within_batch = FALSE,
  linear_reduction = "pca",
  linear_reduction_dims = 50,
  linear_reduction_dims_use = NULL,
  linear_reduction_params = list(),
  force_linear_reduction = FALSE,
  nonlinear_reduction = "umap",
  nonlinear_reduction_dims = c(2, 3),
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  gene_activity_assay = "ACTIVITY",
  MultiMAP_params = list(),
  verbose = TRUE,
  seed = 11
) {
  if (!is.list(MultiMAP_params)) {
    log_message(
      "{.arg MultiMAP_params} must be a list",
      message_type = "error"
    )
  }
  reserved_multimap_params <- c("adatas", "use_reps", "embedding", "seed")
  invalid_multimap_params <- intersect(
    names(MultiMAP_params),
    reserved_multimap_params
  )
  if (length(invalid_multimap_params) > 0) {
    log_message(
      "{.arg MultiMAP_params} contains reserved keys managed by {.fn MultiMAP_integrate}: {.val {invalid_multimap_params}}",
      message_type = "error"
    )
  }

  cluster_algorithm_index <- resolve_cluster_algorithm_index(cluster_algorithm)

  set.seed(seed)
  if (is.null(srt_merge) && is.null(srt_list)) {
    log_message(
      "{.arg srt_list} or {.arg srt_merge} must be provided",
      message_type = "error"
    )
  }
  if (!is.null(srt_list)) {
    srt_merge <- Reduce(merge, srt_list)
  }
  srt_merge_raw <- srt_merge

  assay_pair <- wnn_assays(
    srt = srt_merge,
    assay = assay
  )
  rna_assay <- assay_pair[["rna"]]
  atac_assay <- assay_pair[["atac"]]
  rna_prefix <- resolve_assay_prefix(srt = srt_merge, assay = rna_assay)
  atac_prefix <- resolve_assay_prefix(srt = srt_merge, assay = atac_assay)

  PrepareEnv(modules = "multimap")
  check_python(c("multimap", "scanpy"))

  srt_merge <- atac_add_activity(
    srt = srt_merge,
    assay = atac_assay,
    gene_activity_assay = gene_activity_assay,
    verbose = verbose
  )

  srt_merge <- RunStandardWorkflow(
    object = srt_merge,
    prefix = "Standard",
    assay = c(rna_assay, atac_assay),
    do_normalization = do_normalization,
    normalization_method = normalization_method,
    do_HVF_finding = do_HVF_finding,
    HVF_method = HVF_method,
    nHVF = nHVF,
    HVF = HVF,
    do_scaling = do_scaling,
    vars_to_regress = vars_to_regress,
    regression_model = regression_model,
    linear_reduction = linear_reduction,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_dims_use = linear_reduction_dims_use,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = force_linear_reduction,
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    neighbor_metric = neighbor_metric,
    neighbor_k = neighbor_k,
    cluster_algorithm = cluster_algorithm,
    cluster_resolution = cluster_resolution,
    verbose = verbose,
    seed = seed
  )

  rna_reduction <- paste0(rna_prefix, "pca")
  atac_reduction <- paste0(atac_prefix, "lsi")
  if (!all(c(rna_reduction, atac_reduction) %in% SeuratObject::Reductions(srt_merge))) {
    log_message(
      "MultiMAP requires reductions {.val {c(rna_reduction, atac_reduction)}}",
      message_type = "error"
    )
  }

  rna_hvf <- SeuratObject::VariableFeatures(srt_merge, assay = rna_assay)
  if (length(rna_hvf) == 0) {
    rna_hvf <- SeuratObject::VariableFeatures(srt_merge[[rna_assay]])
  }
  shared_features <- intersect(
    rna_hvf,
    rownames(srt_merge[[gene_activity_assay]])
  )
  if (length(shared_features) < 50) {
    shared_features <- intersect(
      rownames(srt_merge[[rna_assay]]),
      rownames(srt_merge[[gene_activity_assay]])
    )
  }
  if (length(shared_features) < 50) {
    log_message(
      "Need at least 50 shared RNA/gene-activity features for {.pkg MultiMAP}",
      message_type = "error"
    )
  }

  rna_adata <- srt_to_adata(
    object = srt_merge,
    features = shared_features,
    assay_x = rna_assay,
    layer_x = "counts",
    assay_y = NULL,
    reductions = rna_reduction,
    graphs = character(0),
    neighbors = character(0),
    verbose = FALSE
  )
  atac_adata <- srt_to_adata(
    object = srt_merge,
    features = shared_features,
    assay_x = gene_activity_assay,
    layer_x = "counts",
    assay_y = NULL,
    reductions = atac_reduction,
    graphs = character(0),
    neighbors = character(0),
    verbose = FALSE
  )

  rna_names <- paste0(colnames(srt_merge), "__RNA")
  atac_names <- paste0(colnames(srt_merge), "__ATAC")
  rna_adata$obs_names <- rna_names
  atac_adata$obs_names <- atac_names
  rna_adata$obs[["orig_cell"]] <- colnames(srt_merge)
  atac_adata$obs[["orig_cell"]] <- colnames(srt_merge)
  rna_adata$obs[["modality"]] <- "RNA"
  atac_adata$obs[["modality"]] <- "ATAC"

  multimap_params <- MultiMAP_params
  multimap_params[["adatas"]] <- list(rna_adata, atac_adata)
  multimap_params[["use_reps"]] <- c(rna_reduction, atac_reduction)
  multimap_params[["embedding"]] <- multimap_params[["embedding"]] %||% TRUE
  multimap_params[["seed"]] <- multimap_params[["seed"]] %||% as.integer(seed)
  multimap_params[["n_components"]] <- multimap_params[["n_components"]] %||%
    as.integer(max(10L, max(nonlinear_reduction_dims)))

  multimap_result <- run_multimap_python(
    rna_adata = rna_adata,
    atac_adata = atac_adata,
    MultiMAP_params = multimap_params,
    verbose = verbose
  )
  embed <- multimap_result[["embedding"]]
  obs_joint <- multimap_result[["obs"]]
  obs_names_joint <- multimap_result[["obs_names"]]
  if ("orig_cell" %in% colnames(obs_joint)) {
    cell_order <- as.character(obs_joint[["orig_cell"]])
  } else {
    cell_order <- sub(
      pattern = "__(RNA|ATAC)(-[0-9]+)?$",
      replacement = "",
      x = obs_names_joint,
      perl = TRUE
    )
  }
  cell_count <- rowsum(
    matrix(1, nrow = nrow(embed), ncol = 1),
    group = cell_order,
    reorder = FALSE
  )
  if (!all(as.vector(cell_count[, 1]) == 2L)) {
    log_message(
      "Current {.pkg MultiMAP} integration supports paired RNA-ATAC inputs with exactly two modality observations per cell",
      message_type = "error"
    )
  }
  embed_mean <- rowsum(
    embed,
    group = cell_order,
    reorder = FALSE
  )
  embed_mean <- embed_mean / as.vector(cell_count[, 1])
  embed_mean <- embed_mean[colnames(srt_merge), , drop = FALSE]
  colnames(embed_mean) <- paste0("MultiMAP_", seq_len(ncol(embed_mean)))

  srt_merge[["MultiMAP"]] <- CreateDimReducObject(
    embeddings = embed_mean,
    key = "MultiMAP_",
    assay = rna_assay
  )
  dims_use <- seq_len(ncol(embed_mean))
  SeuratObject::DefaultAssay(srt_merge) <- rna_assay

  hvf_use <- SeuratObject::VariableFeatures(srt_merge, assay = rna_assay)
  if (length(hvf_use) == 0) {
    hvf_use <- SeuratObject::VariableFeatures(srt_merge[[rna_assay]])
  }
  if (length(hvf_use) == 0) {
    hvf_use <- shared_features
  }
  srt_merge <- find_neighbors_and_clusters(
    srt = srt_merge,
    reduction = "MultiMAP",
    dims_use = dims_use,
    graph_prefix = "MultiMAP_",
    graph_snn = "MultiMAP_SNN",
    cluster_colname = "MultiMAPclusters",
    HVF = hvf_use,
    neighbor_metric = neighbor_metric,
    neighbor_k = neighbor_k,
    cluster_algorithm = cluster_algorithm,
    cluster_algorithm_index = cluster_algorithm_index,
    cluster_resolution = cluster_resolution,
    verbose = verbose
  )

  srt_merge <- run_nonlinear_reduction(
    srt = srt_merge,
    prefix = "MultiMAP",
    reduction_use = "MultiMAP",
    reduction_dims = dims_use,
    graph_use = "MultiMAP_SNN",
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed,
    verbose = verbose
  )

  srt_merge@misc[["Default_reduction"]] <- if ("MultiMAPUMAP2D" %in% names(srt_merge@reductions)) {
    "MultiMAPUMAP"
  } else {
    "MultiMAP"
  }
  srt_merge@misc[["MultiMAP_reduction_list"]] <- c(rna_reduction, atac_reduction)
  srt_merge@misc[["MultiMAP_shared_features"]] <- shared_features
  SeuratObject::DefaultAssay(srt_merge) <- rna_assay

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_merge,
      pattern = paste0(
        rna_assay,
        "|",
        atac_assay,
        "|",
        gene_activity_assay,
        "|MultiMAP|Default_reduction"
      ),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  }

  srt_merge
}


run_multimap_python <- function(
  rna_adata,
  atac_adata,
  MultiMAP_params = list(),
  verbose = TRUE
) {
  env_cache <- getOption("scop_env_cache", default = NULL)
  python <- env_cache[["python"]] %||%
    tryCatch(
      conda_python(envname = get_envname(), conda = resolve_conda("auto")),
      error = function(...) NULL
    )
  if (is.null(python) || !file.exists(python)) {
    log_message(
      "Unable to resolve python executable for {.pkg MultiMAP}",
      message_type = "error"
    )
  }

  workdir <- tempfile(pattern = "multimap_run_")
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  numba_cache_dir <- file.path(workdir, "numba_cache")
  mpl_config_dir <- file.path(workdir, "matplotlib")
  dir.create(numba_cache_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(mpl_config_dir, recursive = TRUE, showWarnings = FALSE)

  rna_path <- file.path(workdir, "rna.h5ad")
  atac_path <- file.path(workdir, "atac.h5ad")
  embed_out <- file.path(workdir, "multimap_embedding.csv")
  obs_out <- file.path(workdir, "multimap_obs.csv")
  script_path <- file.path(workdir, "run_multimap.py")
  stdout_path <- file.path(workdir, "multimap_stdout.log")
  stderr_path <- file.path(workdir, "multimap_stderr.log")

  rna_adata$write_h5ad(rna_path, compression = "gzip")
  atac_adata$write_h5ad(atac_path, compression = "gzip")

  params <- MultiMAP_params
  params[["adatas"]] <- NULL
  script_body <- c(
    "import anndata as ad",
    "import pandas as pd",
    "",
    "try:",
    "    import MultiMAP as multimap",
    "except ImportError:",
    "    import multimap",
    "",
    sprintf("rna = ad.read_h5ad(%s)", glue_python_literal(rna_path)),
    sprintf("atac = ad.read_h5ad(%s)", glue_python_literal(atac_path)),
    sprintf("params = %s", glue_python_literal(params)),
    "params['adatas'] = [rna, atac]",
    "adata_joint = multimap.Integration(**params)",
    "if 'X_multimap' not in adata_joint.obsm:",
    "    raise ValueError('MultiMAP did not produce obsm[\"X_multimap\"]')",
    "embed = pd.DataFrame(adata_joint.obsm['X_multimap'], index=adata_joint.obs_names)",
    "obs = adata_joint.obs.copy()",
    "obs.insert(0, 'obs_name', adata_joint.obs_names.astype(str))",
    sprintf("embed.to_csv(%s)", glue_python_literal(embed_out)),
    sprintf("obs.to_csv(%s, index=False)", glue_python_literal(obs_out))
  )
  script_lines <- c(
    "def main():",
    paste0("    ", script_body),
    "",
    "if __name__ == '__main__':",
    "    main()"
  )
  writeLines(script_lines, con = script_path, useBytes = TRUE)

  status <- system2(
    command = python,
    args = script_path,
    env = c(
      "PYTHONNOUSERSITE=1",
      "KMP_DUPLICATE_LIB_OK=TRUE",
      "KMP_WARNINGS=0",
      "OMP_NUM_THREADS=1",
      "OPENBLAS_NUM_THREADS=1",
      "MKL_NUM_THREADS=1",
      "VECLIB_MAXIMUM_THREADS=1",
      "NUMEXPR_NUM_THREADS=1",
      sprintf("NUMBA_CACHE_DIR=%s", numba_cache_dir),
      sprintf("MPLCONFIGDIR=%s", mpl_config_dir)
    ),
    stdout = stdout_path,
    stderr = stderr_path
  )
  if (!identical(status, 0L)) {
    stderr_lines <- if (file.exists(stderr_path)) {
      readLines(stderr_path, warn = FALSE)
    } else {
      character(0)
    }
    stdout_lines <- if (file.exists(stdout_path)) {
      readLines(stdout_path, warn = FALSE)
    } else {
      character(0)
    }
    error_lines <- c(utils::tail(stderr_lines, 20), utils::tail(stdout_lines, 20))
    if (length(error_lines) == 0) {
      error_lines <- sprintf("<no output captured; workdir: %s>", workdir)
    }
    log_message(
      "{.pkg MultiMAP} python runner failed:\n{.code {paste(error_lines, collapse = '\n')}}",
      message_type = "error"
    )
  }
  if (!file.exists(embed_out) || !file.exists(obs_out)) {
    log_message(
      "{.pkg MultiMAP} python runner did not produce embedding files. Logs are in {.file {workdir}}",
      message_type = "error"
    )
  }

  log_message(
    "MultiMAP python runner completed",
    verbose = verbose
  )

  obs <- utils::read.csv(obs_out, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"obs_name" %in% colnames(obs)) {
    log_message(
      "{.pkg MultiMAP} python runner output is missing {.field obs_name}",
      message_type = "error"
    )
  }
  list(
    embedding = as.matrix(utils::read.csv(embed_out, row.names = 1, check.names = FALSE)),
    obs = obs,
    obs_names = as.character(obs[["obs_name"]])
  )
}
