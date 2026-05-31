#' @title The GLUE integration function
#'
#' @description
#' Integrate paired RNA and ATAC data using `scglue`.
#' The current implementation supports one `Seurat` object containing one RNA
#' assay and one `ChromatinAssay`, and folds the modality embeddings back to the
#' original cells by averaging the paired RNA/ATAC embeddings.
#'
#' @md
#' @inheritParams integration_scop
#' @param gene_annotation_by How RNA feature names are matched to gene
#' annotation. One of `"auto"`, `"gene_name"`, or `"gene_id"`.
#' Default is `"auto"`.
#' @param GLUE_params A list of parameters passed to
#' `scglue.models.fit_SCGLUE()`.
#' Reserved keys are `"adatas"` and `"graph"`.
#'
#' @return A `Seurat` object.
#'
#' @seealso [integration_scop]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("pbmcmultiome_sub", package = "scop")
#' pbmcmultiome_sub$batch <- rep(c("batch1", "batch2"), length.out = ncol(pbmcmultiome_sub))
#' pbmcmultiome_sub <- GLUE_integrate(
#'   srt_merge = pbmcmultiome_sub,
#'   batch = "batch",
#'   GLUE_params = list(
#'     skip_balance = TRUE,
#'     init_kws = list(latent_dim = 20L),
#'     fit_kws = list(max_epochs = 10L, patience = 2L, reduce_lr_patience = 1L)
#'   )
#' )
#' }
GLUE_integrate <- function(
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
  gene_annotation_by = c("auto", "gene_name", "gene_id"),
  GLUE_params = list(),
  verbose = TRUE,
  seed = 11
) {
  if (!requireNamespace("Signac", quietly = TRUE)) {
    log_message(
      "{.pkg GLUE} integration requires {.pkg Signac}",
      message_type = "error"
    )
  }
  if (!is.list(GLUE_params)) {
    log_message(
      "{.arg GLUE_params} must be a list",
      message_type = "error"
    )
  }
  invalid_glue_params <- intersect(names(GLUE_params), c("adatas", "graph"))
  if (length(invalid_glue_params) > 0) {
    log_message(
      "{.arg GLUE_params} contains reserved keys managed by {.fn GLUE_integrate}: {.val {invalid_glue_params}}",
      message_type = "error"
    )
  }

  gene_annotation_by <- match.arg(gene_annotation_by)
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
  rna_prefix <- standard_scop_assay_prefix(srt = srt_merge, assay = rna_assay)
  atac_prefix <- standard_scop_assay_prefix(srt = srt_merge, assay = atac_assay)

  PrepareEnv(modules = "glue")
  check_python(c("scglue", "scanpy"))

  t_standard <- Sys.time()
  srt_merge <- standard_scop(
    srt = srt_merge,
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
  log_message(
    "GLUE standard_scop completed in {round(as.numeric(difftime(Sys.time(), t_standard, units = 'secs')), 1)}s",
    verbose = verbose
  )

  rna_reduction <- paste0(rna_prefix, "pca")
  atac_reduction <- paste0(atac_prefix, "lsi")
  if (!all(c(rna_reduction, atac_reduction) %in% SeuratObject::Reductions(srt_merge))) {
    log_message(
      "GLUE requires reductions {.val {c(rna_reduction, atac_reduction)}}",
      message_type = "error"
    )
  }

  t_coords <- Sys.time()
  srt_glue <- glue_add_feature_coords(
    srt = srt_merge,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    gene_annotation_by = gene_annotation_by
  )
  rna_feature_data <- GetFeaturesData(srt_glue, assay = rna_assay)
  rna_features_use <- rownames(rna_feature_data)[
    !is.na(rna_feature_data[["chrom"]]) &
      nzchar(as.character(rna_feature_data[["chrom"]])) &
      !is.na(rna_feature_data[["chromStart"]]) &
      !is.na(rna_feature_data[["chromEnd"]]) &
      !is.na(rna_feature_data[["strand"]]) &
      nzchar(as.character(rna_feature_data[["strand"]]))
  ]
  if (length(rna_features_use) < 100L) {
    log_message(
      "Need at least 100 RNA features with genomic coordinates for {.pkg GLUE}",
      message_type = "error"
    )
  }
  log_message(
    "GLUE feature coordinates prepared in {round(as.numeric(difftime(Sys.time(), t_coords, units = 'secs')), 1)}s; RNA features with coords: {length(rna_features_use)}",
    verbose = verbose
  )

  t_adata <- Sys.time()
  old_skip_prepare <- getOption("scop_skip_python_prepare", FALSE)
  options(scop_skip_python_prepare = TRUE)
  on.exit(
    options(scop_skip_python_prepare = old_skip_prepare),
    add = TRUE
  )
  rna_adata <- srt_to_adata(
    srt = srt_glue,
    features = rna_features_use,
    assay_x = rna_assay,
    layer_x = "counts",
    reductions = rna_reduction,
    graphs = character(0),
    neighbors = character(0),
    verbose = FALSE
  )
  atac_adata <- srt_to_adata(
    srt = srt_glue,
    assay_x = atac_assay,
    layer_x = "counts",
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
  rna_adata$layers[["counts"]] <- rna_adata$X$copy()
  rna_adata$obsm[["X_pca"]] <- get_adata_element(rna_adata$obsm, rna_reduction)
  atac_adata$obsm[["X_lsi"]] <- get_adata_element(atac_adata$obsm, atac_reduction)
  log_message(
    "GLUE AnnData export completed in {round(as.numeric(difftime(Sys.time(), t_adata, units = 'secs')), 1)}s; RNA cells/features: {nrow(rna_adata$obs)}/{nrow(rna_adata$var)}, ATAC cells/features: {nrow(atac_adata$obs)}/{nrow(atac_adata$var)}",
    verbose = verbose
  )

  use_batch <- length(batch) == 1 && batch %in% colnames(srt_merge@meta.data)
  glue_params <- GLUE_params
  glue_params[["init_kws"]] <- glue_params[["init_kws"]] %||% list()
  glue_params[["fit_kws"]] <- glue_params[["fit_kws"]] %||% list()
  if (!is.list(glue_params[["init_kws"]])) {
    log_message(
      "{.arg GLUE_params[['init_kws']]} must be a list",
      message_type = "error"
    )
  }
  if (!is.list(glue_params[["fit_kws"]])) {
    log_message(
      "{.arg GLUE_params[['fit_kws']]} must be a list",
      message_type = "error"
    )
  }
  glue_params[["init_kws"]][["random_seed"]] <- glue_params[["init_kws"]][["random_seed"]] %||%
    as.integer(seed)
  glue_params[["fit_kws"]][["directory"]] <- glue_params[["fit_kws"]][["directory"]] %||%
    tempfile(pattern = "glue_")

  t_python <- Sys.time()
  log_message(
    "Starting GLUE python runner",
    verbose = verbose
  )
  glue_embed <- run_glue_python(
    rna_adata = rna_adata,
    atac_adata = atac_adata,
    GLUE_params = glue_params,
    batch = if (use_batch) batch else NULL,
    rna_rep = rna_reduction,
    atac_rep = atac_reduction,
    verbose = verbose
  )
  log_message(
    "GLUE python runner returned in {round(as.numeric(difftime(Sys.time(), t_python, units = 'secs')), 1)}s",
    verbose = verbose
  )
  rna_embed <- glue_embed[["rna"]]
  atac_embed <- glue_embed[["atac"]]
  hvf_nodes <- glue_embed[["hvf_nodes"]]
  rownames(rna_embed) <- rna_names
  rownames(atac_embed) <- atac_names
  embed <- rbind(rna_embed, atac_embed)
  cell_order <- sub(
    pattern = "__(RNA|ATAC)(-[0-9]+)?$",
    replacement = "",
    x = rownames(embed),
    perl = TRUE
  )
  cell_count <- rowsum(
    matrix(1, nrow = nrow(embed), ncol = 1),
    group = cell_order,
    reorder = FALSE
  )
  if (!all(as.vector(cell_count[, 1]) == 2L)) {
    log_message(
      "Current {.pkg GLUE} integration supports paired RNA-ATAC inputs with exactly two modality observations per cell",
      message_type = "error"
    )
  }
  t_fold <- Sys.time()
  embed_mean <- rowsum(
    embed,
    group = cell_order,
    reorder = FALSE
  )
  embed_mean <- embed_mean / as.vector(cell_count[, 1])
  embed_mean <- embed_mean[colnames(srt_merge), , drop = FALSE]
  colnames(embed_mean) <- paste0("GLUE_", seq_len(ncol(embed_mean)))
  log_message(
    "GLUE paired embedding fold-back completed in {round(as.numeric(difftime(Sys.time(), t_fold, units = 'secs')), 1)}s",
    verbose = verbose
  )

  srt_merge[["GLUE"]] <- CreateDimReducObject(
    embeddings = embed_mean,
    key = "GLUE_",
    assay = rna_assay
  )
  dims_use <- seq_len(ncol(embed_mean))
  SeuratObject::DefaultAssay(srt_merge) <- rna_assay

  hvf_use <- SeuratObject::VariableFeatures(srt_merge, assay = rna_assay)
  if (length(hvf_use) == 0) {
    hvf_use <- SeuratObject::VariableFeatures(srt_merge[[rna_assay]])
  }
  if (length(hvf_use) == 0) {
    hvf_use <- rownames(srt_merge[[rna_assay]])
  }

  srt_merge <- find_neighbors_and_clusters(
    srt = srt_merge,
    reduction = "GLUE",
    dims_use = dims_use,
    graph_prefix = "GLUE_",
    graph_snn = "GLUE_SNN",
    cluster_colname = "GLUEclusters",
    HVF = hvf_use,
    neighbor_metric = neighbor_metric,
    neighbor_k = neighbor_k,
    cluster_algorithm = cluster_algorithm,
    cluster_algorithm_index = switch(
      EXPR = tolower(cluster_algorithm),
      "louvain" = 1,
      "louvain_refined" = 2,
      "slm" = 3,
      "leiden" = 4
    ),
    cluster_resolution = cluster_resolution,
    verbose = verbose
  )

  srt_merge <- run_nonlinear_reduction(
    srt = srt_merge,
    prefix = "GLUE",
    reduction_use = "GLUE",
    reduction_dims = dims_use,
    graph_use = "GLUE_SNN",
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed,
    verbose = verbose
  )

  srt_merge@misc[["Default_reduction"]] <- if ("GLUEUMAP2D" %in% names(srt_merge@reductions)) {
    "GLUEUMAP"
  } else {
    "GLUE"
  }
  srt_merge@misc[["GLUE_reduction_list"]] <- c(rna_reduction, atac_reduction)
  srt_merge@misc[["GLUE_graph_vertices"]] <- length(hvf_nodes)
  SeuratObject::DefaultAssay(srt_merge) <- rna_assay

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_merge,
      pattern = "GLUE|Default_reduction",
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  }

  srt_merge
}

glue_add_feature_coords <- function(
  srt,
  rna_assay,
  atac_assay,
  gene_annotation_by = c("auto", "gene_name", "gene_id")
) {
  gene_annotation_by <- match.arg(gene_annotation_by)

  rna_features <- GetFeaturesData(srt, assay = rna_assay)
  rna_features <- glue_add_gene_coords(
    features = rna_features,
    feature_names = glue_feature_names(srt, assay = rna_assay),
    annotation = Signac::Annotation(srt[[atac_assay]]),
    gene_annotation_by = gene_annotation_by
  )
  srt <- AddFeaturesData(
    srt,
    features = rna_features,
    assay = rna_assay
  )

  atac_features <- GetFeaturesData(srt, assay = atac_assay)
  atac_features <- glue_add_peak_coords(
    features = atac_features,
    feature_names = glue_feature_names(srt, assay = atac_assay)
  )
  srt <- AddFeaturesData(
    srt,
    features = atac_features,
    assay = atac_assay
  )

  srt
}

glue_add_gene_coords <- function(
  features,
  feature_names,
  annotation,
  gene_annotation_by = c("auto", "gene_name", "gene_id")
) {
  gene_annotation_by <- match.arg(gene_annotation_by)
  features <- as.data.frame(features)
  if (nrow(features) == 0 && length(feature_names) > 0) {
    features <- data.frame(row.names = feature_names)
  } else if (is.null(rownames(features))) {
    rownames(features) <- feature_names
  }
  annotation_df <- as.data.frame(annotation)
  required_cols <- c("seqnames", "start", "end")
  if (!all(required_cols %in% colnames(annotation_df))) {
    log_message(
      "ATAC gene annotation must contain {.val {required_cols}}",
      message_type = "error"
    )
  }

  coords_existing <- glue_existing_coords(features)
  coords_missing <- is.na(coords_existing[["chrom"]]) |
    is.na(coords_existing[["chromStart"]]) |
    is.na(coords_existing[["chromEnd"]]) |
    is.na(coords_existing[["strand"]])

  gene_tables <- list()
  gene_matches <- c(gene_name = 0L, gene_id = 0L)
  for (field in c("gene_name", "gene_id")) {
    if (!field %in% colnames(annotation_df)) {
      next
    }
    coord_table <- glue_gene_coord_table(
      annotation_df = annotation_df,
      field = field
    )
    gene_tables[[field]] <- coord_table
    gene_matches[[field]] <- sum(rownames(features) %in% rownames(coord_table))
  }

  field_order <- switch(gene_annotation_by,
    auto = names(sort(gene_matches, decreasing = TRUE)),
    gene_name = c("gene_name", "gene_id"),
    gene_id = c("gene_id", "gene_name")
  )
  field_order <- field_order[field_order %in% names(gene_tables)]

  coords_fill <- coords_existing
  for (field in field_order) {
    coord_table <- gene_tables[[field]]
    matched <- coord_table[rownames(features), c("chrom", "chromStart", "chromEnd", "strand"), drop = FALSE]
    keep <- coords_missing & !is.na(matched[["chrom"]])
    if (any(keep)) {
      coords_fill[keep, c("chrom", "chromStart", "chromEnd", "strand")] <- matched[keep, , drop = FALSE]
      coords_missing <- is.na(coords_fill[["chrom"]]) |
        is.na(coords_fill[["chromStart"]]) |
        is.na(coords_fill[["chromEnd"]]) |
        is.na(coords_fill[["strand"]])
    }
    if (!any(coords_missing)) {
      break
    }
  }

  if (sum(!coords_missing) < max(50L, ceiling(length(feature_names) * 0.1))) {
    log_message(
      "Unable to derive enough RNA gene coordinates for {.pkg GLUE} from {.cls ChromatinAssay} annotation",
      message_type = "error"
    )
  }

  features[, c("chrom", "chromStart", "chromEnd", "strand")] <- coords_fill[, c("chrom", "chromStart", "chromEnd", "strand"), drop = FALSE]
  features
}

glue_add_peak_coords <- function(features, feature_names) {
  features <- as.data.frame(features)
  if (nrow(features) == 0 && length(feature_names) > 0) {
    features <- data.frame(row.names = feature_names)
  } else if (is.null(rownames(features))) {
    rownames(features) <- feature_names
  }
  coords_existing <- glue_existing_coords(features)
  coords_missing <- is.na(coords_existing[["chrom"]]) |
    is.na(coords_existing[["chromStart"]]) |
    is.na(coords_existing[["chromEnd"]])
  coords_parsed <- utils::strcapture(
    pattern = "^([^:]+):([0-9]+)-([0-9]+)$",
    x = rownames(features),
    proto = list(
      chrom = character(),
      chromStart = integer(),
      chromEnd = integer()
    )
  )
  coords_parsed_alt <- utils::strcapture(
    pattern = "^([^-]+)-([0-9]+)-([0-9]+)$",
    x = rownames(features),
    proto = list(
      chrom = character(),
      chromStart = integer(),
      chromEnd = integer()
    )
  )
  use_alt <- (is.na(coords_parsed[["chrom"]]) | !nzchar(coords_parsed[["chrom"]])) &
    !is.na(coords_parsed_alt[["chrom"]]) &
    nzchar(coords_parsed_alt[["chrom"]])
  coords_parsed[use_alt, ] <- coords_parsed_alt[use_alt, , drop = FALSE]
  if (any(coords_missing & (is.na(coords_parsed[["chrom"]]) | !nzchar(coords_parsed[["chrom"]])))) {
    log_message(
      "Unable to parse ATAC peak coordinates for {.pkg GLUE}",
      message_type = "error"
    )
  }
  coords_existing[coords_missing, c("chrom", "chromStart", "chromEnd")] <- coords_parsed[coords_missing, , drop = FALSE]
  features[, c("chrom", "chromStart", "chromEnd")] <- coords_existing[, c("chrom", "chromStart", "chromEnd"), drop = FALSE]
  features
}

glue_existing_coords <- function(features) {
  feature_names <- rownames(features)
  if (is.null(feature_names)) {
    feature_names <- seq_len(nrow(features))
  }
  n_features <- length(feature_names)
  coords <- data.frame(
    chrom = rep(NA_character_, n_features),
    chromStart = rep(NA_integer_, n_features),
    chromEnd = rep(NA_integer_, n_features),
    strand = rep(NA_character_, n_features),
    row.names = feature_names
  )
  if ("chrom" %in% colnames(features)) {
    coords[["chrom"]] <- as.character(features[["chrom"]])
  }
  if ("chromStart" %in% colnames(features)) {
    coords[["chromStart"]] <- suppressWarnings(as.integer(features[["chromStart"]]))
  }
  if ("chromEnd" %in% colnames(features)) {
    coords[["chromEnd"]] <- suppressWarnings(as.integer(features[["chromEnd"]]))
  }
  if ("strand" %in% colnames(features)) {
    coords[["strand"]] <- as.character(features[["strand"]])
  }
  coords
}

glue_gene_coord_table <- function(annotation_df, field = c("gene_name", "gene_id")) {
  field <- match.arg(field)
  annotation_df <- annotation_df[
    !is.na(annotation_df[[field]]) &
      nzchar(annotation_df[[field]]) &
      !is.na(annotation_df[["seqnames"]]) &
      !is.na(annotation_df[["start"]]) &
      !is.na(annotation_df[["end"]]), ,
    drop = FALSE
  ]
  keys <- as.character(annotation_df[[field]])
  keep <- !duplicated(keys)
  chrom <- as.character(annotation_df[["seqnames"]][keep])
  strand <- as.character(annotation_df[["strand"]][keep])
  chrom_start <- as.numeric(stats::ave(annotation_df[["start"]], keys, FUN = min))[keep]
  chrom_end <- as.numeric(stats::ave(annotation_df[["end"]], keys, FUN = max))[keep]
  data.frame(
    chrom = chrom,
    chromStart = as.integer(chrom_start),
    chromEnd = as.integer(chrom_end),
    strand = strand,
    row.names = keys[keep]
  )
}

glue_logical_col <- function(x) {
  if (is.logical(x)) {
    return(x %in% TRUE)
  }
  if (is.factor(x)) {
    x <- as.character(x)
  }
  if (is.character(x)) {
    return(toupper(x) %in% "TRUE")
  }
  as.logical(x) %in% TRUE
}

glue_feature_names <- function(srt, assay = NULL) {
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  feature_names <- tryCatch(
    rownames(GetAssayData5(srt, assay = assay, layer = "counts")),
    error = function(...) NULL
  )
  if (is.null(feature_names) || length(feature_names) == 0) {
    feature_names <- tryCatch(
      rownames(GetAssayData5(srt, assay = assay, layer = "data")),
      error = function(...) NULL
    )
  }
  if (is.null(feature_names) || length(feature_names) == 0) {
    feature_names <- rownames(GetFeaturesData(srt, assay = assay))
  }
  if (is.null(feature_names) || length(feature_names) == 0) {
    log_message(
      "Unable to resolve feature names for assay {.val {assay}}",
      message_type = "error"
    )
  }
  feature_names
}

run_glue_python <- function(
  rna_adata,
  atac_adata,
  GLUE_params = list(),
  batch = NULL,
  rna_rep = "X_pca",
  atac_rep = "X_lsi",
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
      "Unable to resolve python executable for {.pkg GLUE}",
      message_type = "error"
    )
  }

  workdir <- tempfile(pattern = "glue_run_")
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  numba_cache_dir <- file.path(workdir, "numba_cache")
  mpl_config_dir <- file.path(workdir, "matplotlib")
  dir.create(numba_cache_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(mpl_config_dir, recursive = TRUE, showWarnings = FALSE)
  rna_path <- file.path(workdir, "rna.h5ad")
  atac_path <- file.path(workdir, "atac.h5ad")
  rna_out <- file.path(workdir, "rna_glue.csv")
  atac_out <- file.path(workdir, "atac_glue.csv")
  meta_out <- file.path(workdir, "glue_meta.txt")
  script_path <- file.path(workdir, "run_glue.py")
  stdout_path <- file.path(workdir, "glue_stdout.log")
  stderr_path <- file.path(workdir, "glue_stderr.log")

  rna_adata <- glue_sanitize_adata_var(rna_adata)
  atac_adata <- glue_sanitize_adata_var(atac_adata)
  rna_adata$write_h5ad(rna_path, compression = "gzip")
  atac_adata$write_h5ad(atac_path, compression = "gzip")

  fit_params <- GLUE_params
  fit_params[["skip_balance"]] <- fit_params[["skip_balance"]] %||% FALSE
  script_body <- c(
    "import anndata as ad",
    "import pandas as pd",
    "import scglue",
    "import sys",
    "",
    sprintf("rna = ad.read_h5ad(%s)", glue_python_literal(rna_path)),
    sprintf("atac = ad.read_h5ad(%s)", glue_python_literal(atac_path)),
    sprintf("batch_col = %s", glue_python_literal(batch)),
    sprintf("rna_rep = %s", glue_python_literal(rna_rep)),
    sprintf("atac_rep = %s", glue_python_literal(atac_rep)),
    sprintf("fit_params = %s", glue_python_literal(fit_params)),
    "if 'counts' not in rna.layers:",
    "    rna.layers['counts'] = rna.X.copy()",
    "rna_mask = (rna.var['chrom'].astype(str) != '') & (rna.var['chromStart'] >= 0) & (rna.var['chromEnd'] >= 0) & (rna.var['strand'].astype(str) != '')",
    "atac_mask = (atac.var['chrom'].astype(str) != '') & (atac.var['chromStart'] >= 0) & (atac.var['chromEnd'] >= 0)",
    "rna = rna[:, rna_mask].copy()",
    "atac = atac[:, atac_mask].copy()",
    "if rna.n_vars == 0 or atac.n_vars == 0:",
    "    raise ValueError(f'GLUE input has zero valid features after coordinate filtering: rna={rna.n_vars}, atac={atac.n_vars}')",
    "rna_cfg = {'use_highly_variable': False, 'use_layer': 'counts', 'use_rep': rna_rep}",
    "atac_cfg = {'use_highly_variable': False, 'use_rep': atac_rep}",
    "if batch_col is not None:",
    "    rna_cfg['use_batch'] = batch_col",
    "    atac_cfg['use_batch'] = batch_col",
    "scglue.models.configure_dataset(rna, 'NB', **rna_cfg)",
    "scglue.models.configure_dataset(atac, 'NB', **atac_cfg)",
    "guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)",
    "active_nodes = list(dict.fromkeys(rna.var_names.tolist() + atac.var_names.tolist()))",
    "guidance_fit = guidance.subgraph(active_nodes).copy()",
    "hvf_nodes = list(guidance_fit.nodes)",
    "glue = scglue.models.fit_SCGLUE({'rna': rna, 'atac': atac}, guidance_fit, **fit_params)",
    "rna_glue = pd.DataFrame(glue.encode_data('rna', rna), index=rna.obs_names)",
    "atac_glue = pd.DataFrame(glue.encode_data('atac', atac), index=atac.obs_names)",
    sprintf("rna_glue.to_csv(%s)", glue_python_literal(rna_out)),
    sprintf("atac_glue.to_csv(%s)", glue_python_literal(atac_out)),
    sprintf("with open(%s, 'w', encoding='utf-8') as fh:", glue_python_literal(meta_out)),
    "    fh.write('\\n'.join(hvf_nodes))"
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
      sprintf("PYTHONNOUSERSITE=%s", "1"),
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
      error_lines <- "<no output captured>"
    }
    log_message(
      "{.pkg GLUE} python runner failed:\n{.code {paste(error_lines, collapse = '\n')}}",
      message_type = "error"
    )
  }
  if (!file.exists(rna_out) || !file.exists(atac_out)) {
    log_message(
      "{.pkg GLUE} python runner did not produce embedding files",
      message_type = "error"
    )
  }

  log_message(
    "GLUE python runner completed",
    verbose = verbose
  )

  list(
    rna = as.matrix(utils::read.csv(rna_out, row.names = 1, check.names = FALSE)),
    atac = as.matrix(utils::read.csv(atac_out, row.names = 1, check.names = FALSE)),
    hvf_nodes = if (file.exists(meta_out)) readLines(meta_out, warn = FALSE) else character(0)
  )
}

glue_python_literal <- function(x) {
  if (is.null(x)) {
    return("None")
  }
  if (is.list(x)) {
    if (length(x) == 0) {
      return("{}")
    }
    if (is.null(names(x)) || any(names(x) == "")) {
      return(sprintf(
        "[%s]",
        paste(vapply(x, glue_python_literal, character(1)), collapse = ", ")
      ))
    }
    return(sprintf(
      "{%s}",
      paste(
        sprintf(
          "%s: %s",
          vapply(names(x), glue_python_literal, character(1)),
          vapply(x, glue_python_literal, character(1))
        ),
        collapse = ", "
      )
    ))
  }
  if (length(x) > 1) {
    return(sprintf(
      "[%s]",
      paste(vapply(as.list(x), glue_python_literal, character(1)), collapse = ", ")
    ))
  }
  if (is.logical(x)) {
    return(if (isTRUE(x)) "True" else "False")
  }
  if (is.numeric(x)) {
    if (is.na(x)) {
      return("None")
    }
    return(as.character(x))
  }
  if (is.character(x)) {
    if (is.na(x) || !nzchar(x)) {
      return("None")
    }
    return(sprintf("'%s'", gsub("'", "\\\\'", x, fixed = TRUE)))
  }
  sprintf("'%s'", gsub("'", "\\\\'", as.character(x), fixed = TRUE))
}

glue_sanitize_adata_var <- function(adata) {
  var <- as.data.frame(py_to_r2(adata$var))
  if ("chrom" %in% colnames(var)) {
    var[["chrom"]] <- as.character(var[["chrom"]])
    var[["chrom"]][is.na(var[["chrom"]])] <- ""
  }
  if ("chromStart" %in% colnames(var)) {
    var[["chromStart"]] <- suppressWarnings(as.integer(var[["chromStart"]]))
    var[["chromStart"]][is.na(var[["chromStart"]])] <- -1L
  }
  if ("chromEnd" %in% colnames(var)) {
    var[["chromEnd"]] <- suppressWarnings(as.integer(var[["chromEnd"]]))
    var[["chromEnd"]][is.na(var[["chromEnd"]])] <- -1L
  }
  if ("strand" %in% colnames(var)) {
    var[["strand"]] <- as.character(var[["strand"]])
    var[["strand"]][is.na(var[["strand"]])] <- ""
  }
  adata$var <- var
  adata$var_names <- rownames(var)
  adata
}
