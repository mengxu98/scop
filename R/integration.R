#' @title The Uncorrected integration function
#'
#' @inheritParams integration_scop
#'
#' @export
Uncorrected_integrate <- function(
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
  verbose = TRUE,
  seed = 11
) {
  if (length(linear_reduction) > 1) {
    log_message(
      "Only the first method in the {.arg linear_reduction} will be used",
      message_type = "warning",
      verbose = verbose
    )
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srt_merge)) {
    reduc_test <- c(reduc_test, SeuratObject::Reductions(srt_merge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    log_message(
      "{.arg linear_reduction} must be one of {.val {reduc_test}}",
      message_type = "error"
    )
  }
  if (
    !is.null(linear_reduction_dims_use) &&
      max(linear_reduction_dims_use) > linear_reduction_dims
  ) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }

  nonlinear_reductions <- c(
    "umap",
    "umap-naive",
    "tsne",
    "dm",
    "phate",
    "pacmap",
    "trimap",
    "largevis",
    "fr"
  )
  if (any(!nonlinear_reduction %in% nonlinear_reductions)) {
    log_message(
      "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
      message_type = "error"
    )
  }
  cluster_algorithms <- c("louvain", "slm", "leiden")
  if (!cluster_algorithm %in% cluster_algorithms) {
    log_message(
      "{.arg cluster_algorithm} must be one of {.val {cluster_algorithms}}",
      message_type = "error"
    )
  }
  cluster_algorithm_index <- switch(
    EXPR = tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "{.arg srt_list} and {.arg srt_merge} were all empty",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "{.arg srt_list} and {.arg srt_merge} have different cells",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- CheckDataList(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      verbose = verbose,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srt_merge <- Reduce(merge, srt_list)
    SeuratObject::VariableFeatures(srt_merge) <- HVF
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    checked <- CheckDataMerge(
      srt_merge = srt_merge,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      verbose = verbose,
      seed = seed
    )
    srt_merge <- checked[["srt_merge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    log_message(
      "{.arg normalization_method} is {.val TFIDF}. Use {.pkg lsi} workflow..."
    )
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  log_message(
    "Perform {.pkg Uncorrected} integration"
  )
  scale_features <- rownames(
    GetAssayData5(
      srt_merge,
      layer = "scale.data",
      assay = SeuratObject::DefaultAssay(srt_merge)
    )
  )
  if (
    isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))
  ) {
    if (normalization_method != "SCT") {
      log_message(
        "Perform {.fn Seurat::ScaleData}",
        verbose = verbose
      )
      srt_merge <- Seurat::ScaleData(
        object = srt_merge,
        split.by = if (isTRUE(scale_within_batch)) batch else NULL,
        assay = SeuratObject::DefaultAssay(srt_merge),
        features = HVF,
        vars.to.regress = vars_to_regress,
        model.use = regression_model,
        verbose = FALSE
      )
    }
  }

  log_message(
    "Perform {.val {linear_reduction}} linear dimension reduction",
    verbose = verbose
  )
  srt_merge <- RunDimsReduction(
    srt_merge,
    prefix = "Uncorrected",
    features = HVF,
    assay = SeuratObject::DefaultAssay(srt_merge),
    linear_reduction = linear_reduction,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = force_linear_reduction,
    verbose = verbose,
    seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- resolve_linear_dims_use(
      srt = srt_merge,
      reduction = paste0("Uncorrected", linear_reduction),
      normalization_method = normalization_method,
      reduction_method = linear_reduction
    )
  }

  srt_merge <- find_neighbors_and_clusters(
    srt = srt_merge,
    reduction = paste0("Uncorrected", linear_reduction),
    dims_use = linear_reduction_dims_use,
    graph_prefix = "Uncorrected_",
    graph_snn = "Uncorrected_SNN",
    cluster_colname = paste0("Uncorrected", linear_reduction, "clusters"),
    HVF = HVF,
    neighbor_metric = neighbor_metric,
    neighbor_k = neighbor_k,
    cluster_algorithm = cluster_algorithm,
    cluster_algorithm_index = cluster_algorithm_index,
    cluster_resolution = cluster_resolution,
    verbose = verbose
  )

  srt_merge <- run_nonlinear_reduction(
    srt = srt_merge,
    prefix = "Uncorrected",
    reduction_use = paste0("Uncorrected", linear_reduction),
    reduction_dims = linear_reduction_dims_use,
    graph_use = "Uncorrected_SNN",
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed,
    verbose = verbose
  )

  SeuratObject::DefaultAssay(srt_merge) <- assay
  SeuratObject::VariableFeatures(srt_merge) <- srt_merge@misc[[
    "Uncorrected_HVF"
  ]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_merge,
      pattern = paste0(assay, "|Uncorrected|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_merge)
  }
}

#' @title The WNN integration function
#'
#' @inheritParams integration_scop
#'
#' @export
#' @examples
#' data("pbmcmultiome_sub", package = "scop")
#' pbmcmultiome_sub$batch <- rep(c("batch1", "batch2"), length.out = ncol(pbmcmultiome_sub))
#' pbmcmultiome_sub <- WNN_integrate(
#'   srt_merge = pbmcmultiome_sub,
#'   batch = "batch",
#'   linear_reduction_dims = 20,
#'   linear_reduction_dims_use = 1:10
#' )
WNN_integrate <- function(
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
  verbose = TRUE,
  seed = 11
) {
  cluster_algorithms <- c("louvain", "slm", "leiden")
  if (!cluster_algorithm %in% cluster_algorithms) {
    log_message(
      "{.arg cluster_algorithm} must be one of {.val {cluster_algorithms}}",
      message_type = "error"
    )
  }
  cluster_algorithm_index <- switch(
    EXPR = tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

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

  rna_reduction <- paste0(rna_prefix, "pca")
  atac_reduction <- paste0(atac_prefix, "lsi")
  if (!all(c(rna_reduction, atac_reduction) %in% SeuratObject::Reductions(srt_merge))) {
    log_message(
      "WNN requires reductions {.val {c(rna_reduction, atac_reduction)}}",
      message_type = "error"
    )
  }

  rna_dims_use <- wnn_dims(
    srt = srt_merge,
    reduction = rna_reduction,
    dims_use = linear_reduction_dims_use,
    reduction_method = "pca",
    normalization_method = "LogNormalize"
  )
  atac_dims_use <- wnn_dims(
    srt = srt_merge,
    reduction = atac_reduction,
    dims_use = if (is.null(linear_reduction_dims_use)) NULL else linear_reduction_dims_use,
    reduction_method = "svd",
    normalization_method = "TFIDF"
  )

  neighbor_k_use <- min(as.integer(neighbor_k), max(1L, ncol(srt_merge) - 1L))
  knn_range_use <- min(
    max(neighbor_k_use + 1L, neighbor_k_use * 4L),
    max(1L, ncol(srt_merge) - 1L)
  )
  if (!identical(neighbor_k_use, as.integer(neighbor_k))) {
    log_message(
      "Adjust neighbor k from {.val {neighbor_k}} to {.val {neighbor_k_use}} for small-sample WNN graph construction",
      verbose = verbose
    )
  }
  if (knn_range_use < 200L) {
    log_message(
      "Adjust WNN knn.range to {.val {knn_range_use}} for small-sample graph construction",
      verbose = verbose
    )
  }

  log_message(
    "Perform {.pkg WNN} integration using {.pkg {rna_reduction}} and {.pkg {atac_reduction}}",
    verbose = verbose
  )
  SeuratObject::DefaultAssay(srt_merge) <- rna_assay
  srt_merge <- Seurat::FindMultiModalNeighbors(
    object = srt_merge,
    reduction.list = list(rna_reduction, atac_reduction),
    dims.list = list(rna_dims_use, atac_dims_use),
    k.nn = neighbor_k_use,
    knn.range = knn_range_use,
    knn.graph.name = "WNNKNN",
    snn.graph.name = "WNNSNN",
    weighted.nn.name = "WNN",
    modality.weight.name = c(
      paste0(rna_prefix, ".weight"),
      paste0(atac_prefix, ".weight")
    ),
    verbose = verbose
  )

  hvf_use <- SeuratObject::VariableFeatures(srt_merge, assay = rna_assay)
  if (length(hvf_use) == 0) {
    hvf_use <- SeuratObject::VariableFeatures(srt_merge[[rna_assay]])
  }
  if (length(hvf_use) == 0) {
    hvf_use <- utils::head(rownames(srt_merge[[rna_assay]]), 2000L)
  }
  srt_merge <- find_neighbors_and_clusters(
    srt = srt_merge,
    reduction = rna_reduction,
    dims_use = rna_dims_use,
    graph_prefix = "WNN_",
    graph_snn = "WNNSNN",
    cluster_colname = "WNNclusters",
    HVF = hvf_use,
    neighbor_metric = neighbor_metric,
    neighbor_k = neighbor_k_use,
    cluster_algorithm = cluster_algorithm,
    cluster_algorithm_index = cluster_algorithm_index,
    cluster_resolution = cluster_resolution,
    run_find_neighbors = FALSE,
    verbose = verbose
  )

  srt_merge <- run_wnn_reduction(
    srt = srt_merge,
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    verbose = verbose,
    seed = seed
  )

  wnn_reductions <- grep(
    "^WNN(UMAP|FR)",
    names(srt_merge@reductions),
    value = TRUE
  )
  srt_merge@misc[["Default_reduction"]] <- if ("WNNUMAP2D" %in% names(srt_merge@reductions)) {
    "WNNUMAP"
  } else if (length(wnn_reductions) > 0) {
    sub("(2D|3D)$", "", wnn_reductions[[1]])
  } else {
    srt_merge@misc[["Default_reduction"]] %||% NULL
  }
  srt_merge@misc[["WNN_reduction_list"]] <- c(rna_reduction, atac_reduction)
  srt_merge@misc[["WNN_dims_list"]] <- list(
    rna = rna_dims_use,
    atac = atac_dims_use
  )
  SeuratObject::DefaultAssay(srt_merge) <- rna_assay

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_merge,
      pattern = paste0(rna_assay, "|", atac_assay, "|WNN|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  }

  srt_merge
}

#' @title The MultiMAP integration function
#'
#' @inheritParams integration_scop
#' @param gene_activity_assay Name of the gene activity assay used to provide
#' a shared feature space for RNA-ATAC integration. Default is `"ACTIVITY"`.
#' @param MultiMAP_params A list of parameters passed to `MultiMAP::Integration`.
#' The following keys are managed internally and should not be supplied:
#' `"adatas"`, `"use_reps"`, `"embedding"`, and `"seed"`.
#' Default is `list()`.
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

  cluster_algorithms <- c("louvain", "slm", "leiden")
  if (!cluster_algorithm %in% cluster_algorithms) {
    log_message(
      "{.arg cluster_algorithm} must be one of {.val {cluster_algorithms}}",
      message_type = "error"
    )
  }
  cluster_algorithm_index <- switch(
    EXPR = tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

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

  PrepareEnv(modules = "multimap")
  check_python(c("multimap", "scanpy"))

  srt_merge <- atac_add_activity(
    srt = srt_merge,
    assay = atac_assay,
    gene_activity_assay = gene_activity_assay,
    verbose = verbose
  )

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
    srt = srt_merge,
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
    srt = srt_merge,
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

wnn_assays <- function(srt, assay = NULL) {
  assays_available <- SeuratObject::Assays(srt)
  chrom_assays <- assays_available[vapply(
    assays_available,
    function(x) inherits(srt[[x]], "ChromatinAssay"),
    logical(1)
  )]
  rna_assays <- setdiff(assays_available, chrom_assays)
  if (length(chrom_assays) == 0 || length(rna_assays) == 0) {
    log_message(
      "WNN requires at least one RNA assay and one {.cls ChromatinAssay}",
      message_type = "error"
    )
  }
  if (is.null(assay)) {
    assay_default <- SeuratObject::DefaultAssay(srt)
    rna_assay <- if (assay_default %in% rna_assays) assay_default else rna_assays[[1]]
    atac_assay <- if ("peaks" %in% chrom_assays) "peaks" else chrom_assays[[1]]
    return(list(rna = rna_assay, atac = atac_assay))
  }

  assay <- unique(as.character(assay))
  if (length(assay) == 1) {
    if (assay %in% chrom_assays) {
      return(list(
        rna = if ("RNA" %in% rna_assays) "RNA" else rna_assays[[1]],
        atac = assay
      ))
    }
    if (assay %in% rna_assays) {
      return(list(
        rna = assay,
        atac = if ("peaks" %in% chrom_assays) "peaks" else chrom_assays[[1]]
      ))
    }
  }

  rna_assay <- assay[assay %in% rna_assays][[1]] %||% NULL
  atac_assay <- assay[assay %in% chrom_assays][[1]] %||% NULL
  if (is.null(rna_assay) || is.null(atac_assay)) {
    log_message(
      "{.arg assay} for WNN must include one RNA assay and one {.cls ChromatinAssay}",
      message_type = "error"
    )
  }
  list(rna = rna_assay, atac = atac_assay)
}

wnn_dims <- function(
  srt,
  reduction,
  dims_use = NULL,
  reduction_method,
  normalization_method
) {
  dims_use <- dims_use %||% resolve_linear_dims_use(
    srt = srt,
    reduction = reduction,
    linear_reduction_dims_use = NULL,
    normalization_method = normalization_method,
    reduction_method = reduction_method,
    verbose = FALSE
  )
  available_dims <- seq_len(ncol(Seurat::Embeddings(srt, reduction = reduction)))
  dims_use <- intersect(as.integer(dims_use), available_dims)
  if (length(dims_use) == 0) {
    log_message(
      "No valid dimensions remain for {.pkg {reduction}}",
      message_type = "error"
    )
  }
  dims_use
}

run_wnn_reduction <- function(
  srt,
  nonlinear_reduction,
  nonlinear_reduction_dims,
  nonlinear_reduction_params,
  force_nonlinear_reduction,
  verbose,
  seed
) {
  supported <- c("umap", "umap-naive", "fr")
  unsupported <- setdiff(nonlinear_reduction, supported)
  if (length(unsupported) > 0) {
    log_message(
      "WNN currently supports only {.val {supported}} nonlinear reductions. Skip {.val {unsupported}}",
      message_type = "warning",
      verbose = verbose
    )
  }
  nonlinear_use <- intersect(nonlinear_reduction, supported)
  if (length(nonlinear_use) == 0) {
    log_message(
      "No supported WNN nonlinear reduction was requested. Fall back to {.val umap}",
      message_type = "warning",
      verbose = verbose
    )
    nonlinear_use <- "umap"
  }
  for (nr in nonlinear_use) {
    for (n in nonlinear_reduction_dims) {
      if (identical(nr, "fr")) {
        srt <- RunDimsReduction(
          srt = srt,
          prefix = "WNN",
          graph_use = "WNNSNN",
          nonlinear_reduction = nr,
          nonlinear_reduction_dims = n,
          nonlinear_reduction_params = nonlinear_reduction_params,
          force_nonlinear_reduction = force_nonlinear_reduction,
          verbose = verbose,
          seed = seed
        )
      } else {
        srt <- RunDimsReduction(
          srt = srt,
          prefix = "WNN",
          neighbor_use = "WNN",
          nonlinear_reduction = nr,
          nonlinear_reduction_dims = n,
          nonlinear_reduction_params = nonlinear_reduction_params,
          force_nonlinear_reduction = force_nonlinear_reduction,
          verbose = verbose,
          seed = seed
        )
      }
    }
  }
  srt
}

#' @title The Seurat integration function
#'
#' @inheritParams integration_scop
#' @param FindIntegrationAnchors_params A list of parameters for the Seurat::FindIntegrationAnchors function.
#' Default is `list()`.
#' @param IntegrateData_params A list of parameters for the Seurat::IntegrateData function.
#' Default is `list()`.
#' @param IntegrateEmbeddings_params A list of parameters for the Seurat::IntegrateEmbeddings function.
#' Default is `list()`.
#' @export
Seurat_integrate <- function(
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
  FindIntegrationAnchors_params = list(),
  IntegrateData_params = list(),
  IntegrateEmbeddings_params = list(),
  verbose = TRUE,
  seed = 11
) {
  if (length(linear_reduction) > 1) {
    log_message(
      "Only the first method in the {.arg linear_reduction} will be used",
      message_type = "warning",
      verbose = verbose
    )
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "ica", "svd", "nmf", "mds", "glmpca")
  if (!is.null(srt_merge)) {
    reduc_test <- c(reduc_test, SeuratObject::Reductions(srt_merge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    log_message(
      "{.arg linear_reduction} must be one of {.val {reduc_test}}",
      message_type = "error"
    )
  }
  if (
    !is.null(linear_reduction_dims_use) &&
      max(linear_reduction_dims_use) > linear_reduction_dims
  ) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  nonlinear_reductions <- c(
    "umap",
    "umap-naive",
    "tsne",
    "dm",
    "phate",
    "pacmap",
    "trimap",
    "largevis",
    "fr"
  )
  if (any(!nonlinear_reduction %in% nonlinear_reductions)) {
    log_message(
      "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
      message_type = "error"
    )
  }
  cluster_algorithms <- c("louvain", "slm", "leiden")
  if (!cluster_algorithm %in% cluster_algorithms) {
    log_message(
      "{.arg cluster_algorithm} must be one of {.val {cluster_algorithms}}",
      message_type = "error"
    )
  }
  cluster_algorithm_index <- switch(
    EXPR = tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "{.arg srt_list} and {.arg srt_merge} were all empty",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "{.arg srt_list} and {.arg srt_merge} have different cells",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- CheckDataList(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    if (normalization_method == "TFIDF") {
      srt_merge <- Reduce(merge, srt_list)
      SeuratObject::VariableFeatures(srt_merge) <- HVF
    }
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    checked <- CheckDataMerge(
      srt_merge = srt_merge,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      verbose = verbose,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    srt_merge <- checked[["srt_merge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  reduction_key <- FindIntegrationAnchors_params[["reduction"]]
  if (!is.null(reduction_key) && normalization_method != "TFIDF") {
    reduction_key <- tolower(
      gsub("[^a-z]", "", as.character(reduction_key)[1L], perl = TRUE)
    )
    if (reduction_key %in% c("cca", "ccaintegration")) {
      FindIntegrationAnchors_params[["reduction"]] <- "cca"
    } else if (reduction_key %in% c("rpca", "rpcaintegration")) {
      FindIntegrationAnchors_params[["reduction"]] <- "rpca"
    }
  }

  if (min(sapply(srt_list, ncol)) < 50) {
    log_message(
      "The cell count in some batches is lower than 50, which may not be suitable for the current integration method",
      message_type = "warning"
    )
    if (interactive()) {
      answer <- utils::askYesNo("Are you sure to continue?", default = FALSE)
      if (isFALSE(answer)) {
        return(srt_merge)
      }
    } else {
      log_message(
        "Non-interactive session detected. Continue integration after warning",
        message_type = "warning"
      )
    }
  }

  if (normalization_method == "TFIDF") {
    log_message(
      "{.arg normalization_method} is {.val TFIDF}. Use {.pkg rlsi} integration workflow..."
    )
    do_scaling <- FALSE
    linear_reduction <- "svd"
    FindIntegrationAnchors_params[["reduction"]] <- "rlsi"
    if (is.null(FindIntegrationAnchors_params[["dims"]])) {
      max_anchor_dim <- min(
        linear_reduction_dims,
        30,
        min(sapply(srt_list, ncol)) - 1L
      )
      FindIntegrationAnchors_params[["dims"]] <- if (max_anchor_dim >= 2) {
        2:max_anchor_dim
      } else {
        1L
      }
    }
    srt_merge <- Signac::RunTFIDF(
      object = srt_merge,
      assay = SeuratObject::DefaultAssay(srt_merge),
      verbose = FALSE
    )
    srt_merge <- RunDimsReduction(
      srt_merge,
      prefix = "",
      features = HVF,
      assay = SeuratObject::DefaultAssay(srt_merge),
      linear_reduction = "svd",
      linear_reduction_dims = linear_reduction_dims,
      linear_reduction_params = linear_reduction_params,
      force_linear_reduction = force_linear_reduction,
      verbose = verbose,
      seed = seed
    )
    srt_merge[["lsi"]] <- srt_merge[["svd"]]
    for (i in seq_along(srt_list)) {
      srt <- srt_list[[i]]
      log_message(
        "Perform {.pkg svd} linear dimension reduction on {.val {i}} of {.arg srt_list}"
      )
      srt <- RunDimsReduction(
        srt,
        prefix = "",
        features = HVF,
        assay = SeuratObject::DefaultAssay(srt),
        linear_reduction = "svd",
        linear_reduction_dims = linear_reduction_dims,
        linear_reduction_params = linear_reduction_params,
        force_linear_reduction = force_linear_reduction,
        verbose = verbose,
        seed = seed
      )
      srt[["lsi"]] <- srt[["svd"]]
      srt_list[[i]] <- srt
    }
  }

  if (isTRUE(FindIntegrationAnchors_params[["reduction"]] == "rpca")) {
    log_message("Use {.pkg rpca} integration workflow...")
    for (i in seq_along(srt_list)) {
      srt <- srt_list[[i]]
      scale_features <- rownames(
        GetAssayData5(
          srt,
          layer = "scale.data",
          assay = SeuratObject::DefaultAssay(srt)
        )
      )
      if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))) {
        log_message(
          "Perform {.fn Seurat::ScaleData} on {.arg srt}"
        )
        srt <- Seurat::ScaleData(
          object = srt,
          assay = SeuratObject::DefaultAssay(srt),
          features = HVF,
          vars.to.regress = vars_to_regress,
          model.use = regression_model,
          verbose = FALSE
        )
      }
      log_message(
        "Perform {.pkg pca} linear dimension reduction on {.val {i}} of {.arg srt_list}"
      )
      srt <- RunDimsReduction(
        srt,
        prefix = "",
        features = HVF,
        assay = SeuratObject::DefaultAssay(srt),
        linear_reduction = "pca",
        linear_reduction_dims = linear_reduction_dims,
        linear_reduction_params = linear_reduction_params,
        force_linear_reduction = force_linear_reduction,
        verbose = verbose,
        seed = seed
      )
      srt_list[[i]] <- srt
    }
  }

  if (is.null(names(srt_list))) {
    names(srt_list) <- paste0("srt_", seq_along(srt_list))
  }

  if (normalization_method %in% c("LogNormalize", "SCT")) {
    log_message("Perform FindIntegrationAnchors")
    params1 <- list(
      object.list = srt_list,
      normalization.method = normalization_method,
      anchor.features = HVF,
      verbose = FALSE
    )
    for (nm in names(FindIntegrationAnchors_params)) {
      params1[[nm]] <- FindIntegrationAnchors_params[[nm]]
    }
    srt_anchors <- invoke_fun(
      Seurat::FindIntegrationAnchors,
      params1
    )

    log_message("Perform {.pkg Seurat} integration")
    params2 <- list(
      anchorset = srt_anchors,
      new.assay.name = "Seuratcorrected",
      normalization.method = normalization_method,
      features.to.integrate = HVF,
      verbose = FALSE
    )
    for (nm in names(IntegrateData_params)) {
      params2[[nm]] <- IntegrateData_params[[nm]]
    }
    srt_integrated <- invoke_fun(Seurat::IntegrateData, params2)

    SeuratObject::DefaultAssay(srt_integrated) <- "Seuratcorrected"
    SeuratObject::VariableFeatures(srt_integrated[["Seuratcorrected"]]) <- HVF

    scale_features <- rownames(
      GetAssayData5(
        srt_integrated,
        layer = "scale.data",
        assay = SeuratObject::DefaultAssay(srt_integrated)
      )
    )
    if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))) {
      log_message("Perform ScaleData on {.arg srt_integrated}")
      srt_integrated <- Seurat::ScaleData(
        object = srt_integrated,
        split.by = if (isTRUE(scale_within_batch)) batch else NULL,
        assay = SeuratObject::DefaultAssay(srt_integrated),
        features = HVF,
        vars.to.regress = vars_to_regress,
        model.use = regression_model,
        verbose = FALSE
      )
    }

    log_message(
      "Perform {.val {linear_reduction}} linear dimension reduction"
    )
    srt_integrated <- RunDimsReduction(
      srt_integrated,
      prefix = "Seurat",
      features = HVF,
      assay = SeuratObject::DefaultAssay(srt_integrated),
      linear_reduction = linear_reduction,
      linear_reduction_dims = linear_reduction_dims,
      linear_reduction_params = linear_reduction_params,
      force_linear_reduction = force_linear_reduction,
      verbose = verbose,
      seed = seed
    )
    if (is.null(linear_reduction_dims_use)) {
      linear_reduction_dims_use <- srt_integrated@reductions[[paste0(
        "Seurat",
        linear_reduction
      )]]@misc[["dims_estimate"]] %||%
        1:linear_reduction_dims
    }
  } else if (normalization_method == "TFIDF") {
    log_message(
      "Perform {.fn FindIntegrationAnchors} with {.arg reduction = rlsi}"
    )
    params1 <- list(
      object.list = srt_list,
      normalization.method = "LogNormalize",
      anchor.features = HVF,
      reduction = "rlsi",
      verbose = FALSE
    )
    for (nm in names(FindIntegrationAnchors_params)) {
      params1[[nm]] <- FindIntegrationAnchors_params[[nm]]
    }
    srt_anchors <- invoke_fun(Seurat::FindIntegrationAnchors, params1)

    log_message("Perform {.pkg Seurat} integration")
    params2 <- list(
      anchorset = srt_anchors,
      reductions = srt_merge[["lsi"]],
      new.reduction.name = "Seuratlsi",
      verbose = FALSE
    )
    for (nm in names(IntegrateEmbeddings_params)) {
      params2[[nm]] <- IntegrateEmbeddings_params[[nm]]
    }
    srt_integrated <- invoke_fun(IntegrateEmbeddings, params2)

    if (is.null(linear_reduction_dims_use)) {
      linear_reduction_dims_use <- 2:max(srt_integrated@reductions[[paste0(
        "Seurat",
        linear_reduction
      )]]@misc[["dims_estimate"]]) %||%
        2:linear_reduction_dims
    }
    linear_reduction <- "lsi"
  }

  srt_integrated <- find_neighbors_and_clusters(
    srt = srt_integrated,
    reduction = paste0("Seurat", linear_reduction),
    dims_use = linear_reduction_dims_use,
    graph_prefix = "Seurat_",
    graph_snn = "Seurat_SNN",
    cluster_colname = "Seuratclusters",
    HVF = HVF,
    neighbor_metric = neighbor_metric,
    neighbor_k = neighbor_k,
    cluster_algorithm = cluster_algorithm,
    cluster_algorithm_index = cluster_algorithm_index,
    cluster_resolution = cluster_resolution,
    verbose = verbose
  )

  srt_integrated <- run_nonlinear_reduction(
    srt = srt_integrated,
    prefix = "Seurat",
    reduction_use = paste0("Seurat", linear_reduction),
    reduction_dims = linear_reduction_dims_use,
    graph_use = "Seurat_SNN",
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed,
    verbose = verbose
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[["Seurat_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_integrated,
      pattern = paste0(assay, "|Seurat|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_integrated)
  }
}

#' @title The scVI integration function
#'
#' @inheritParams integration_scop
#' @param scVI_dims_use A vector specifying the dimensions returned by scVI that will be utilized for downstream cell cluster finding and nonlinear reduction.
#' If set to NULL, all the returned dimensions will be used by default.
#' @param model A string indicating the scVI model to be used.
#' Options are "SCVI", "PEAKVI", and "POISSONVI".
#' Default is `"SCVI"`.
#' @param SCVI_params A list of parameters for the SCVI model.
#' Default is `list()`.
#' @param PEAKVI_params A list of parameters for the PEAKVI model.
#' Default is `list()`.
#' @param POISSONVI_params A list of parameters for the POISSONVI model.
#' Default is `list()`.
#' @param train_params A list of parameters passed to the model `train()` method.
#' Default is `list()`.
#' @param cores An integer setting the number of threads for `scVI`.
#' Default is `1`.
#' @export
#' @examples
#' \dontrun{
#' data("pbmcmultiome_sub", package = "scop")
#' pbmcmultiome_sub$batch <- rep(c("batch1", "batch2"), length.out = ncol(pbmcmultiome_sub))
#' pbmcmultiome_sub <- scVI_integrate(
#'   srt_merge = pbmcmultiome_sub,
#'   batch = "batch",
#'   assay = "peaks",
#'   model = "PEAKVI",
#'   train_params = list(max_epochs = 2L)
#' )
#' pbmcmultiome_sub <- scVI_integrate(
#'   srt_merge = pbmcmultiome_sub,
#'   batch = "batch",
#'   assay = "peaks",
#'   model = "POISSONVI",
#'   train_params = list(max_epochs = 2L)
#' )
#' }
scVI_integrate <- function(
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
  scVI_dims_use = NULL,
  nonlinear_reduction = "umap",
  nonlinear_reduction_dims = c(2, 3),
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  model = "SCVI",
  SCVI_params = list(),
  PEAKVI_params = list(),
  POISSONVI_params = list(),
  train_params = list(),
  cores = 1,
  verbose = TRUE,
  seed = 11
) {
  model <- toupper(model)
  if (!model %in% c("SCVI", "PEAKVI", "POISSONVI")) {
    log_message(
      "{.arg model} must be one of {.val {c('SCVI', 'PEAKVI', 'POISSONVI')}}",
      message_type = "error"
    )
  }
  nonlinear_reductions <- c(
    "umap",
    "umap-naive",
    "tsne",
    "dm",
    "phate",
    "pacmap",
    "trimap",
    "largevis",
    "fr"
  )
  if (any(!nonlinear_reduction %in% nonlinear_reductions)) {
    log_message(
      "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
      message_type = "error"
    )
  }
  cluster_algorithms <- c("louvain", "slm", "leiden")
  if (!cluster_algorithm %in% cluster_algorithms) {
    log_message(
      "{.arg cluster_algorithm} must be one of {.val {cluster_algorithms}}",
      message_type = "error"
    )
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  PrepareEnv(modules = "scvi")
  check_python("scvi-tools")
  scvi <- reticulate::import("scvi")
  scipy <- reticulate::import("scipy")
  set.seed(seed)

  scvi$settings$num_threads <- as.integer(cores)

  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "{.arg srt_list} and {.arg srt_merge} were all empty",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "{.arg srt_list} and {.arg srt_merge} have different cells",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- CheckDataList(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srt_merge <- Reduce(merge, srt_list)
    SeuratObject::VariableFeatures(srt_merge) <- HVF
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    checked <- CheckDataMerge(
      srt_merge = srt_merge,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      verbose = verbose,
      seed = seed
    )
    srt_merge <- checked[["srt_merge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }
  if (
    identical(model, "PEAKVI") &&
      !inherits(srt_merge[[assay]], "ChromatinAssay")
  ) {
    log_message(
      "{.arg model = 'PEAKVI'} requires {.cls ChromatinAssay}",
      message_type = "error"
    )
  }
  if (
    identical(model, "POISSONVI") &&
      !inherits(srt_merge[[assay]], "ChromatinAssay")
  ) {
    log_message(
      "{.arg model = 'POISSONVI'} requires {.cls ChromatinAssay}",
      message_type = "error"
    )
  }

  reduction_name <- switch(model,
    SCVI = "scVI",
    PEAKVI = "PeakVI",
    POISSONVI = "PoissonVI"
  )
  reduction_key <- paste0(reduction_name, "_")
  graph_prefix <- reduction_key
  graph_snn <- paste0(reduction_name, "_SNN")
  cluster_colname <- paste0(reduction_name, "clusters")
  hvf_key <- paste0(reduction_name, "_HVF")
  append_tag <- reduction_name

  adata <- srt_to_adata(
    srt_merge,
    features = HVF,
    assay_x = SeuratObject::DefaultAssay(srt_merge),
    assay_y = NULL,
    verbose = FALSE
  )
  adata[["X"]] <- scipy$sparse$csr_matrix(adata[["X"]])

  if (model == "SCVI") {
    scvi$model$SCVI$setup_anndata(adata, batch_key = batch)
    model_params <- list(
      adata = adata
    )
    for (nm in names(SCVI_params)) {
      model_params[[nm]] <- SCVI_params[[nm]]
    }
    model <- invoke_fun(scvi$model$SCVI, model_params)
    invoke_fun(model$train, train_params)
    srt_integrated <- srt_merge
    srt_merge <- NULL
    corrected <- Matrix::t(
      as_matrix(
        model$get_normalized_expression()
      )
    )
    srt_integrated[["scVIcorrected"]] <- SeuratObject::CreateAssayObject(
      counts = corrected
    )
    SeuratObject::DefaultAssay(srt_integrated) <- "scVIcorrected"
    SeuratObject::VariableFeatures(srt_integrated[["scVIcorrected"]]) <- HVF
  } else if (model == "PEAKVI") {
    log_message("Assay is ChromatinAssay. Using PeakVI workflow.")
    scvi$model$PEAKVI$setup_anndata(adata, batch_key = batch)
    model_params <- list(
      adata = adata
    )
    for (nm in names(PEAKVI_params)) {
      model_params[[nm]] <- PEAKVI_params[[nm]]
    }
    model <- invoke_fun(scvi$model$PEAKVI, model_params)
    invoke_fun(model$train, train_params)
    srt_integrated <- srt_merge
    srt_merge <- NULL
  } else if (model == "POISSONVI") {
    log_message("Assay is ChromatinAssay. Using PoissonVI workflow.")
    scvi$external$POISSONVI$setup_anndata(adata, batch_key = batch)
    model_params <- list(
      adata = adata
    )
    for (nm in names(POISSONVI_params)) {
      model_params[[nm]] <- POISSONVI_params[[nm]]
    }
    model <- invoke_fun(scvi$external$POISSONVI, model_params)
    invoke_fun(model$train, train_params)
    srt_integrated <- srt_merge
    srt_merge <- NULL
  }

  latent <- as_matrix(model$get_latent_representation())
  rownames(latent) <- colnames(srt_integrated)
  colnames(latent) <- paste0(reduction_key, seq_len(ncol(latent)))
  srt_integrated[[reduction_name]] <- CreateDimReducObject(
    embeddings = latent,
    key = reduction_key,
    assay = SeuratObject::DefaultAssay(srt_integrated)
  )
  if (is.null(scVI_dims_use)) {
    scVI_dims_use <- 1:ncol(srt_integrated[[reduction_name]]@cell.embeddings)
  }

  srt_integrated <- find_neighbors_and_clusters(
    srt = srt_integrated,
    reduction = reduction_name,
    dims_use = scVI_dims_use,
    graph_prefix = graph_prefix,
    graph_snn = graph_snn,
    cluster_colname = cluster_colname,
    HVF = HVF,
    neighbor_metric = neighbor_metric,
    neighbor_k = neighbor_k,
    cluster_algorithm = cluster_algorithm,
    cluster_algorithm_index = cluster_algorithm_index,
    cluster_resolution = cluster_resolution,
    verbose = verbose
  )

  srt_integrated <- run_nonlinear_reduction(
    srt = srt_integrated,
    prefix = reduction_name,
    reduction_use = reduction_name,
    reduction_dims = scVI_dims_use,
    graph_use = graph_snn,
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed,
    verbose = verbose
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[[hvf_key]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_integrated,
      pattern = paste0(assay, "|", append_tag, "|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_integrated)
  }
}

#' @title The MNN integration function
#'
#' @inheritParams integration_scop
#' @param mnnCorrect_params A list of parameters for the batchelor::mnnCorrect function,
#' default is an empty list.
#' @export
MNN_integrate <- function(
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
  mnnCorrect_params = list(),
  verbose = TRUE,
  seed = 11
) {
  if (length(linear_reduction) > 1) {
    log_message(
      "Only the first method in the 'linear_reduction' will be used.",
      message_type = "warning"
    )
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srt_merge)) {
    reduc_test <- c(reduc_test, SeuratObject::Reductions(srt_merge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    log_message(
      "{.arg linear_reduction} must be one of {.val {reduc_test}}",
      message_type = "error"
    )
  }
  if (
    !is.null(linear_reduction_dims_use) &&
      max(linear_reduction_dims_use) > linear_reduction_dims
  ) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  nonlinear_reductions <- c(
    "umap",
    "umap-naive",
    "tsne",
    "dm",
    "phate",
    "pacmap",
    "trimap",
    "largevis",
    "fr"
  )
  if (any(!nonlinear_reduction %in% nonlinear_reductions)) {
    log_message(
      "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
      message_type = "error"
    )
  }
  cluster_algorithms <- c("louvain", "slm", "leiden")
  if (!cluster_algorithm %in% cluster_algorithms) {
    log_message(
      "{.arg cluster_algorithm} must be one of {.val {cluster_algorithms}}",
      message_type = "error"
    )
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_r("batchelor", verbose = FALSE)
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "{.arg srt_list} and {.arg srt_merge} were all empty",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "srt_list and srt_merge have different cells.",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- CheckDataList(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    srt_list <- Seurat::SplitObject(
      object = srt_merge,
      split.by = batch
    )
    checked <- CheckDataList(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (is.null(srt_merge) && !is.null(srt_list)) {
    srt_merge <- Reduce(merge, srt_list)
  }

  if (normalization_method == "TFIDF") {
    log_message(
      "{.arg normalization_method} is {.val TFIDF}. Use {.pkg lsi} workflow..."
    )
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }
  mnn_fallback_warned <- FALSE
  sce_list <- lapply(
    srt_list,
    function(srt) {
      data_matrix <- GetAssayData5(
        srt,
        layer = "data",
        assay = SeuratObject::DefaultAssay(srt)
      )
      if (
        is.null(dim(data_matrix)) ||
          nrow(data_matrix) == 0 ||
          ncol(data_matrix) == 0
      ) {
        if (!mnn_fallback_warned) {
          log_message(
            "Layer {.val data} is empty for MNN input. Fallback to {.val counts} with {.fn log1p} transform.",
            message_type = "warning",
            verbose = verbose
          )
          mnn_fallback_warned <- TRUE
        }
        data_matrix <- GetAssayData5(
          srt,
          layer = "counts",
          assay = SeuratObject::DefaultAssay(srt)
        )
        if (inherits(data_matrix, "dgCMatrix")) {
          data_matrix <- as_matrix(data_matrix)
        }
        data_matrix <- log1p(data_matrix)
      }
      data_matrix <- data_matrix[HVF, , drop = FALSE]
      if (inherits(data_matrix, "dgCMatrix")) {
        data_matrix <- as_matrix(data_matrix)
      }
      if (nrow(data_matrix) == 0 || ncol(data_matrix) == 0) {
        log_message(
          "No available features/cells for MNN after preparing {.val logcounts} matrix.",
          message_type = "error"
        )
      }
      sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(logcounts = data_matrix)
      )
      return(sce)
    }
  )
  if (is.null(names(sce_list))) {
    names(sce_list) <- paste0("sce_", seq_along(sce_list))
  }

  log_message("Perform {.pkg MNN} integration")
  params <- c(
    sce_list,
    list(cos.norm.out = FALSE)
  )
  for (nm in names(mnnCorrect_params)) {
    params[[nm]] <- mnnCorrect_params[[nm]]
  }
  out <- invoke_fun(batchelor::mnnCorrect, params)

  srt_integrated <- srt_merge
  srt_merge <- NULL
  srt_integrated[["MNNcorrected"]] <- CreateAssayObject(
    counts = out@assays@data$corrected
  )
  SeuratObject::VariableFeatures(srt_integrated[["MNNcorrected"]]) <- HVF
  SeuratObject::DefaultAssay(srt_integrated) <- "MNNcorrected"
  scale_features <- rownames(
    GetAssayData5(
      srt_integrated,
      layer = "scale.data",
      assay = SeuratObject::DefaultAssay(srt_integrated)
    )
  )
  if (
    isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))
  ) {
    log_message("Perform ScaleData")
    srt_integrated <- Seurat::ScaleData(
      object = srt_integrated,
      split.by = if (isTRUE(scale_within_batch)) batch else NULL,
      assay = SeuratObject::DefaultAssay(srt_integrated),
      features = HVF,
      vars.to.regress = vars_to_regress,
      model.use = regression_model,
      verbose = FALSE
    )
  }

  log_message(
    "Perform {.val {linear_reduction}} linear dimension reduction",
    verbose = verbose
  )
  srt_integrated <- RunDimsReduction(
    srt_integrated,
    prefix = "MNN",
    features = HVF,
    assay = SeuratObject::DefaultAssay(srt_integrated),
    linear_reduction = linear_reduction,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = force_linear_reduction,
    verbose = verbose,
    seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- resolve_linear_dims_use(
      srt = srt_integrated,
      reduction = paste0("MNN", linear_reduction),
      normalization_method = normalization_method,
      reduction_method = linear_reduction
    )
  }

  srt_integrated <- find_neighbors_and_clusters(
    srt = srt_integrated,
    reduction = paste0("MNN", linear_reduction),
    dims_use = linear_reduction_dims_use,
    graph_prefix = "MNN_",
    graph_snn = "MNN_SNN",
    cluster_colname = "MNNclusters",
    HVF = HVF,
    neighbor_metric = neighbor_metric,
    neighbor_k = neighbor_k,
    cluster_algorithm = cluster_algorithm,
    cluster_algorithm_index = cluster_algorithm_index,
    cluster_resolution = cluster_resolution,
    verbose = verbose
  )

  srt_integrated <- run_nonlinear_reduction(
    srt = srt_integrated,
    prefix = "MNN",
    reduction_use = paste0("MNN", linear_reduction),
    reduction_dims = linear_reduction_dims_use,
    graph_use = "MNN_SNN",
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed,
    verbose = verbose
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[[
    "MNN_HVF"
  ]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_integrated,
      pattern = paste0(assay, "|MNN|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_integrated)
  }
}

#' @title The fastMNN integration function
#'
#' @inheritParams integration_scop
#' @param fastMNN_dims_use A vector specifying the dimensions returned by fastMNN that will be utilized for downstream cell cluster finding and nonlinear reduction.
#' If set to NULL, all the returned dimensions will be used by default.
#' @param fastMNN_params A list of parameters for the batchelor::fastMNN function, default is an empty list.
#'
#' @export
fastMNN_integrate <- function(
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
  fastMNN_dims_use = NULL,
  nonlinear_reduction = "umap",
  nonlinear_reduction_dims = c(2, 3),
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  fastMNN_params = list(),
  verbose = TRUE,
  seed = 11
) {
  if (
    any(
      !nonlinear_reduction %in%
        c(
          "umap",
          "umap-naive",
          "tsne",
          "dm",
          "phate",
          "pacmap",
          "trimap",
          "largevis",
          "fr"
        )
    )
  ) {
    log_message(
      "'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'.",
      message_type = "error"
    )
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    log_message(
      "'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.",
      message_type = "error"
    )
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_r("batchelor", verbose = FALSE)
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "{.arg srt_list} and {.arg srt_merge} were all empty",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "srt_list and srt_merge have different cells.",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- CheckDataList(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    srt_list <- Seurat::SplitObject(object = srt_merge, split.by = batch)
    checked <- CheckDataList(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (is.null(srt_merge) && !is.null(srt_list)) {
    srt_merge <- Reduce(merge, srt_list)
  }

  fastmnn_fallback_warned <- FALSE
  sce_list <- lapply(srt_list, function(srt) {
    data_matrix <- GetAssayData5(
      srt,
      layer = "data",
      assay = SeuratObject::DefaultAssay(srt)
    )
    if (
      is.null(dim(data_matrix)) ||
        nrow(data_matrix) == 0 ||
        ncol(data_matrix) == 0
    ) {
      if (!fastmnn_fallback_warned) {
        log_message(
          "Layer {.val {'data'}} is empty for fastMNN input. Fallback to {.val {'counts'}} with {.fn log1p} transform.",
          message_type = "warning",
          verbose = verbose
        )
        fastmnn_fallback_warned <- TRUE
      }
      data_matrix <- GetAssayData5(
        srt,
        layer = "counts",
        assay = SeuratObject::DefaultAssay(srt)
      )
      data_matrix <- log1p(data_matrix)
    }
    data_matrix <- data_matrix[HVF, , drop = FALSE]
    if (nrow(data_matrix) == 0 || ncol(data_matrix) == 0) {
      log_message(
        "No available features/cells for fastMNN after preparing {.val {'logcounts'}} matrix.",
        message_type = "error"
      )
    }
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(logcounts = data_matrix)
    )
    return(sce)
  })
  if (is.null(names(sce_list))) {
    names(sce_list) <- paste0("sce_", seq_along(sce_list))
  }

  log_message("Perform {.pkg fastMNN} integration")
  params <- c(
    sce_list,
    list()
  )
  for (nm in names(fastMNN_params)) {
    params[[nm]] <- fastMNN_params[[nm]]
  }
  out <- invoke_fun(batchelor::fastMNN, params)

  srt_integrated <- srt_merge
  srt_merge <- NULL
  srt_integrated[["fastMNNcorrected"]] <- CreateAssayObject(
    counts = out@assays@data$reconstructed
  )
  SeuratObject::DefaultAssay(srt_integrated) <- "fastMNNcorrected"
  SeuratObject::VariableFeatures(srt_integrated[["fastMNNcorrected"]]) <- HVF
  reduction <- out@int_colData$reducedDims$corrected
  colnames(reduction) <- paste0("fastMNN_", seq_len(ncol(reduction)))
  srt_integrated[["fastMNN"]] <- CreateDimReducObject(
    embeddings = reduction,
    key = "fastMNN_",
    assay = "fastMNNcorrected"
  )

  if (is.null(fastMNN_dims_use)) {
    fastMNN_dims_use <- 1:ncol(srt_integrated[["fastMNN"]]@cell.embeddings)
  }

  srt_integrated <- find_neighbors_and_clusters(
    srt = srt_integrated,
    reduction = "fastMNN",
    dims_use = fastMNN_dims_use,
    graph_prefix = "fastMNN_",
    graph_snn = "fastMNN_SNN",
    cluster_colname = "fastMNNclusters",
    HVF = HVF,
    neighbor_metric = neighbor_metric,
    neighbor_k = neighbor_k,
    cluster_algorithm = cluster_algorithm,
    cluster_algorithm_index = cluster_algorithm_index,
    cluster_resolution = cluster_resolution,
    verbose = verbose
  )

  srt_integrated <- run_nonlinear_reduction(
    srt = srt_integrated,
    prefix = "fastMNN",
    reduction_use = "fastMNN",
    reduction_dims = fastMNN_dims_use,
    graph_use = "fastMNN_SNN",
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed,
    verbose = verbose
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[[
    "fastMNN_HVF"
  ]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_integrated,
      pattern = paste0(assay, "|fastMNN|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_integrated)
  }
}

#' @title The Harmony integration function
#'
#' @inheritParams integration_scop
#' @param harmony_dims_use A vector specifying the dimensions returned by RunHarmony that will be utilized for downstream cell cluster finding and nonlinear reduction.
#' If set to NULL, all the returned dimensions will be used by default.
#' @param RunHarmony_params A list of parameters for [harmony::RunHarmony], default is an empty list.
#'
#' @export
Harmony_integrate <- function(
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
  harmony_dims_use = NULL,
  nonlinear_reduction = "umap",
  nonlinear_reduction_dims = c(2, 3),
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  RunHarmony_params = list(),
  verbose = TRUE,
  seed = 11
) {
  if (length(linear_reduction) > 1) {
    log_message(
      "Only the first of {.val {linear_reduction}} will be used",
      message_type = "warning",
      verbose = verbose
    )
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srt_merge)) {
    reduc_test <- c(reduc_test, SeuratObject::Reductions(srt_merge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    log_message(
      "{.arg linear_reduction} must be one of {.val {reduc_test}}",
      message_type = "error"
    )
  }
  if (
    !is.null(linear_reduction_dims_use) &&
      max(linear_reduction_dims_use) > linear_reduction_dims
  ) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  nonlinear_reductions <- c(
    "umap",
    "umap-naive",
    "tsne",
    "dm",
    "phate",
    "pacmap",
    "trimap",
    "largevis",
    "fr"
  )
  if (any(!nonlinear_reduction %in% nonlinear_reductions)) {
    log_message(
      "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
      message_type = "error"
    )
  }
  cluster_algorithms <- c("louvain", "slm", "leiden")
  if (!cluster_algorithm %in% cluster_algorithms) {
    log_message(
      "{.arg cluster_algorithm} must be one of {.val {cluster_algorithms}}",
      message_type = "error"
    )
  }
  cluster_algorithm_index <- switch(
    EXPR = tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "{.arg srt_list} and {.arg srt_merge} were all empty",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "{.arg srt_list} and {.arg srt_merge} have different cells",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- CheckDataList(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srt_merge <- Reduce(merge, srt_list)
    SeuratObject::VariableFeatures(srt_merge) <- HVF
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    checked <- CheckDataMerge(
      srt_merge = srt_merge,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      verbose = verbose,
      seed = seed
    )
    srt_merge <- checked[["srt_merge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    log_message(
      "{.arg normalization_method} is {.val TFIDF}. Use {.pkg lsi} workflow..."
    )
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }
  scale_features <- rownames(
    GetAssayData5(
      srt_merge,
      layer = "scale.data",
      assay = SeuratObject::DefaultAssay(srt_merge)
    )
  )
  if (
    isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))
  ) {
    log_message("Perform {.fn Seurat::ScaleData}")
    srt_merge <- Seurat::ScaleData(
      object = srt_merge,
      split.by = if (isTRUE(scale_within_batch)) batch else NULL,
      assay = SeuratObject::DefaultAssay(srt_merge),
      features = HVF,
      vars.to.regress = vars_to_regress,
      model.use = regression_model,
      verbose = FALSE
    )
  }

  log_message(
    "Perform linear dimension reduction({.val {linear_reduction}})"
  )
  srt_merge <- RunDimsReduction(
    srt_merge,
    prefix = "Harmony",
    features = HVF,
    assay = SeuratObject::DefaultAssay(srt_merge),
    linear_reduction = linear_reduction,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = force_linear_reduction,
    verbose = verbose,
    seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- resolve_linear_dims_use(
      srt = srt_merge,
      reduction = paste0("Harmony", linear_reduction),
      normalization_method = normalization_method,
      reduction_method = linear_reduction
    )
  }

  log_message(
    "Perform {.pkg Harmony} integration",
    verbose = verbose
  )
  log_message(
    "Using {.val {paste0('Harmony', linear_reduction)}} ({.val {min(linear_reduction_dims_use)}}:{.val {max(linear_reduction_dims_use)}}) as input",
    verbose = verbose
  )
  params <- list(
    object = srt_merge,
    group.by.vars = batch,
    assay = SeuratObject::DefaultAssay(srt_merge),
    reduction = paste0("Harmony", linear_reduction),
    dims.use = linear_reduction_dims_use,
    reduction.name = "Harmony",
    reduction.key = "Harmony_",
    verbose = FALSE
  )
  feature_num <- nrow(
    GetAssayData5(
      srt_merge,
      layer = "scale.data",
      assay = SeuratObject::DefaultAssay(srt_merge)
    )
  )
  if (feature_num == 0) {
    params[["project.dim"]] <- FALSE
  }
  if (!is.null(RunHarmony_params[["reduction.save"]])) {
    RunHarmony_params[["reduction.name"]] <- RunHarmony_params[["reduction.name"]] %||%
      RunHarmony_params[["reduction.save"]]
    RunHarmony_params[["reduction.save"]] <- NULL
  }
  for (nm in names(RunHarmony_params)) {
    params[[nm]] <- RunHarmony_params[[nm]]
  }
  srt_integrated <- invoke_fun(RunHarmony2, params)
  harmony_reduction <- params[["reduction.name"]] %||% "Harmony"

  if (is.null(harmony_dims_use)) {
    harmony_dims_use <- seq_len(
      ncol(
        srt_integrated[[harmony_reduction]]@cell.embeddings
      )
    )
  }

  srt_integrated <- find_neighbors_and_clusters(
    srt = srt_integrated,
    reduction = harmony_reduction,
    dims_use = harmony_dims_use,
    graph_prefix = "Harmony_",
    graph_snn = "Harmony_SNN",
    cluster_colname = "Harmonyclusters",
    HVF = HVF,
    neighbor_metric = neighbor_metric,
    neighbor_k = neighbor_k,
    cluster_algorithm = cluster_algorithm,
    cluster_algorithm_index = cluster_algorithm_index,
    cluster_resolution = cluster_resolution,
    verbose = verbose
  )

  srt_integrated <- run_nonlinear_reduction(
    srt = srt_integrated,
    prefix = "Harmony",
    reduction_use = harmony_reduction,
    reduction_dims = harmony_dims_use,
    graph_use = "Harmony_SNN",
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed,
    verbose = verbose
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[[
    "Harmony_HVF"
  ]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_integrated,
      pattern = paste0(assay, "|Harmony|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_integrated)
  }
}

#' @title The Scanorama integration function
#'
#' @inheritParams integration_scop
#' @param Scanorama_dims_use  A vector specifying the dimensions returned by Scanorama that will be utilized for downstream cell cluster finding and nonlinear reduction.
#' If set to NULL, all the returned dimensions will be used by default.
#' @param return_corrected Logical indicating whether to return the corrected data.
#' Default is `FALSE`.
#' @param Scanorama_params A list of parameters for the scanorama.correct function.
#' Default is `list()`.
#'
#' @export
Scanorama_integrate <- function(
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
  Scanorama_dims_use = NULL,
  nonlinear_reduction = "umap",
  nonlinear_reduction_dims = c(2, 3),
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  return_corrected = FALSE,
  Scanorama_params = list(),
  verbose = TRUE,
  seed = 11
) {
  PrepareEnv(modules = "scanorama")

  nonlinear_reductions <- c(
    "umap",
    "umap-naive",
    "tsne",
    "dm",
    "phate",
    "pacmap",
    "trimap",
    "largevis",
    "fr"
  )
  if (any(!nonlinear_reduction %in% nonlinear_reductions)) {
    log_message(
      "'nonlinear_reduction' must be one of ",
      paste(nonlinear_reductions, collapse = ", "),
      message_type = "error"
    )
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    log_message(
      "'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.",
      message_type = "error"
    )
  }
  cluster_algorithm_index <- switch(
    EXPR = tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_python("scanorama")
  scanorama <- reticulate::import("scanorama")
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "srt_list and srt_merge were all empty.",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "srt_list and srt_merge have different cells.",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- CheckDataList(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    srt_list <- Seurat::SplitObject(object = srt_merge, split.by = batch)
    checked <- CheckDataList(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }
  srt_integrated <- Reduce(merge, srt_list)

  log_message("Perform {.pkg Scanorama} integration")
  assaylist <- list()
  genelist <- list()
  for (i in seq_along(srt_list)) {
    assaylist[[i]] <- Matrix::t(
      as_matrix(
        GetAssayData5(
          object = srt_list[[i]],
          layer = "data",
          assay = SeuratObject::DefaultAssay(srt_list[[i]])
        )[HVF, , drop = FALSE]
      )
    )
    genelist[[i]] <- HVF
  }
  if (isTRUE(return_corrected)) {
    params <- list(
      datasets_full = assaylist,
      genes_list = genelist,
      return_dimred = TRUE,
      return_dense = TRUE,
      verbose = FALSE
    )
    for (nm in names(Scanorama_params)) {
      params[[nm]] <- Scanorama_params[[nm]]
    }
    corrected <- invoke_fun(scanorama$correct, params)

    cor_value <- Matrix::t(invoke_fun(rbind, corrected[[2]]))
    rownames(cor_value) <- corrected[[3]]
    colnames(cor_value) <- unlist(sapply(assaylist, rownames))
    srt_integrated[["Scanoramacorrected"]] <- CreateAssayObject(
      data = cor_value
    )
    SeuratObject::VariableFeatures(srt_integrated[[
      "Scanoramacorrected"
    ]]) <- HVF

    dim_reduction <- invoke_fun(rbind, corrected[[1]])
    rownames(dim_reduction) <- unlist(sapply(assaylist, rownames))
    colnames(dim_reduction) <- paste0(
      "Scanorama_",
      seq_len(ncol(dim_reduction))
    )
  } else {
    params <- list(
      datasets_full = assaylist,
      genes_list = genelist,
      verbose = FALSE
    )
    for (nm in names(Scanorama_params)) {
      params[[nm]] <- Scanorama_params[[nm]]
    }
    integrated <- invoke_fun(scanorama$integrate, params)

    dim_reduction <- invoke_fun(rbind, integrated[[1]])
    rownames(dim_reduction) <- unlist(sapply(assaylist, rownames))
    colnames(dim_reduction) <- paste0(
      "Scanorama_",
      seq_len(ncol(dim_reduction))
    )
  }
  srt_integrated[["Scanorama"]] <- CreateDimReducObject(
    embeddings = dim_reduction,
    key = "Scanorama_",
    assay = SeuratObject::DefaultAssay(srt_integrated)
  )

  if (is.null(Scanorama_dims_use)) {
    Scanorama_dims_use <- 1:ncol(srt_integrated[["Scanorama"]]@cell.embeddings)
  }

  srt_integrated <- find_neighbors_and_clusters(
    srt = srt_integrated,
    reduction = "Scanorama",
    dims_use = Scanorama_dims_use,
    graph_prefix = "Scanorama_",
    graph_snn = "Scanorama_SNN",
    cluster_colname = "Scanoramaclusters",
    HVF = HVF,
    neighbor_metric = neighbor_metric,
    neighbor_k = neighbor_k,
    cluster_algorithm = cluster_algorithm,
    cluster_algorithm_index = cluster_algorithm_index,
    cluster_resolution = cluster_resolution,
    verbose = verbose
  )

  srt_integrated <- run_nonlinear_reduction(
    srt = srt_integrated,
    prefix = "Scanorama",
    reduction_use = "Scanorama",
    reduction_dims = Scanorama_dims_use,
    graph_use = "Scanorama_SNN",
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed,
    verbose = verbose
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[[
    "Scanorama_HVF"
  ]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_integrated,
      pattern = paste0(assay, "|Scanorama|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_integrated)
  }
}

#' @title The BBKNN integration function
#'
#' @inheritParams integration_scop
#' @param bbknn_params A list of parameters for the bbknn.matrix.bbknn function, default is an empty list.
#' @param bbknn_backend Backend used for BBKNN graph construction. `"python"`
#' uses the original `bbknn.matrix.bbknn` implementation. `"r"` uses a native
#' R/C++ balanced KNN graph with approximate connectivity weights.
#'
#' @export
BBKNN_integrate <- function(
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
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  bbknn_params = list(),
  bbknn_backend = c("python", "r"),
  verbose = TRUE,
  seed = 11
) {
  bbknn_backend <- match.arg(bbknn_backend)
  if (identical(bbknn_backend, "python")) {
    PrepareEnv(modules = "bbknn")
  }

  if (length(linear_reduction) > 1) {
    log_message(
      "Only the first method in the {.arg linear_reduction} will be used",
      message_type = "warning"
    )
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c(
    "pca",
    "svd",
    "ica",
    "nmf",
    "mds",
    "glmpca"
  )
  if (!is.null(srt_merge)) {
    reduc_test <- c(reduc_test, SeuratObject::Reductions(srt_merge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    log_message(
      "{.arg linear_reduction} must be one of {.val {reduc_test}}",
      message_type = "error"
    )
  }
  if (
    !is.null(linear_reduction_dims_use) &&
      max(linear_reduction_dims_use) > linear_reduction_dims
  ) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  nonlinear_reductions <- c("umap", "umap-naive", "fr")
  if (any(!nonlinear_reduction %in% nonlinear_reductions)) {
    log_message(
      "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
      message_type = "error"
    )
  }
  cluster_algorithms <- c("louvain", "slm", "leiden")
  if (!cluster_algorithm %in% cluster_algorithms) {
    log_message(
      "{.arg cluster_algorithm} must be one of {.val {cluster_algorithms}}",
      message_type = "error"
    )
  }
  cluster_algorithm_index <- switch(
    EXPR = tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  if (identical(bbknn_backend, "python")) {
    check_python("bbknn")
    bbknn <- reticulate::import("bbknn")
  }
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "{.arg srt_list} and {.arg srt_merge} were all empty",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "{.arg srt_list} and {.arg srt_merge} have different cells",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- CheckDataList(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srt_merge <- Reduce(merge, srt_list)
    SeuratObject::VariableFeatures(srt_merge) <- HVF
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    checked <- CheckDataMerge(
      srt_merge = srt_merge,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      verbose = verbose,
      seed = seed
    )
    srt_merge <- checked[["srt_merge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    log_message(
      "{.arg normalization_method} is {.pkg TFIDF}. Use {.pkg lsi} workflow..."
    )
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }
  scale_features <- rownames(
    GetAssayData5(
      srt_merge,
      layer = "scale.data",
      assay = SeuratObject::DefaultAssay(srt_merge)
    )
  )
  if (
    isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))
  ) {
    log_message("Perform {.fn Seurat::ScaleData}")
    srt_merge <- Seurat::ScaleData(
      object = srt_merge,
      split.by = if (isTRUE(scale_within_batch)) batch else NULL,
      assay = SeuratObject::DefaultAssay(srt_merge),
      features = HVF,
      vars.to.regress = vars_to_regress,
      model.use = regression_model,
      verbose = FALSE
    )
  }

  log_message(
    "Perform {.val {linear_reduction}} linear dimension reduction"
  )
  srt_merge <- RunDimsReduction(
    srt_merge,
    prefix = "BBKNN",
    features = HVF,
    assay = SeuratObject::DefaultAssay(srt_merge),
    linear_reduction = linear_reduction,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = force_linear_reduction,
    verbose = verbose,
    seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- resolve_linear_dims_use(
      srt = srt_merge,
      reduction = paste0("BBKNN", linear_reduction),
      normalization_method = normalization_method,
      reduction_method = linear_reduction
    )
  }

  log_message(
    "Perform {.pkg BBKNN} integration with {.arg bbknn_backend = {bbknn_backend}}",
    verbose = verbose
  )
  log_message(
    "Using {.val {paste0('BBKNN', linear_reduction)}} ({.val {min(linear_reduction_dims_use)}}:{.val {max(linear_reduction_dims_use)}}) as input",
    verbose = verbose
  )
  emb <- Embeddings(srt_merge, reduction = paste0("BBKNN", linear_reduction))[,
    linear_reduction_dims_use,
    drop = FALSE
  ]
  params <- list(
    pca = emb,
    batch_list = srt_merge[[batch, drop = TRUE]]
  )
  for (nm in names(bbknn_params)) {
    params[[nm]] <- bbknn_params[[nm]]
  }
  bem <- if (identical(bbknn_backend, "python")) {
    invoke_fun(bbknn$matrix$bbknn, params)
  } else {
    params[["verbose"]] <- verbose
    invoke_fun(bbknn_native_graph, params)
  }
  n.neighbors <- bem[[3]]$n_neighbors
  srt_integrated <- srt_merge

  bbknn_graph <- SeuratObject::as.sparse(
    bem[[2]][1:nrow(bem[[2]]), , drop = FALSE]
  )
  rownames(bbknn_graph) <- colnames(bbknn_graph) <- rownames(emb)
  bbknn_graph <- SeuratObject::as.Graph(bbknn_graph)
  bbknn_graph@assay.used <- SeuratObject::DefaultAssay(srt_integrated)
  srt_integrated@graphs[["BBKNN"]] <- bbknn_graph

  bbknn_dist <- Matrix::t(
    SeuratObject::as.sparse(
      bem[[1]][1:nrow(bem[[1]]), , drop = FALSE]
    )
  )
  rownames(bbknn_dist) <- colnames(bbknn_dist) <- rownames(emb)
  bbknn_dist <- SeuratObject::as.Graph(bbknn_dist)
  bbknn_dist@assay.used <- SeuratObject::DefaultAssay(srt_integrated)
  srt_integrated@graphs[["BBKNN_dist"]] <- bbknn_dist

  val <- split(
    bbknn_dist@x,
    rep(seq_len(ncol(bbknn_dist)), diff(bbknn_dist@p))
  )
  pos <- split(
    bbknn_dist@i + 1,
    rep(seq_len(ncol(bbknn_dist)), diff(bbknn_dist@p))
  )
  idx <- Matrix::t(
    mapply(
      function(x, y) {
        out <- y[utils::head(order(x, decreasing = FALSE), n.neighbors - 1)]
        length(out) <- n.neighbors - 1
        out
      },
      x = val,
      y = pos
    )
  )
  idx[is.na(idx)] <- sample(
    seq_len(nrow(idx)),
    size = sum(is.na(idx)),
    replace = TRUE
  )
  idx <- cbind(seq_len(nrow(idx)), idx)
  dist <- Matrix::t(mapply(
    function(x, y) {
      out <- y[utils::head(order(x, decreasing = FALSE), n.neighbors - 1)]
      length(out) <- n.neighbors - 1
      out[is.na(out)] <- 0
      out
    },
    x = val,
    y = val
  ))
  dist <- cbind(0, dist)
  srt_integrated[["BBKNN_neighbors"]] <- methods::new(
    Class = "Neighbor",
    nn.idx = idx,
    nn.dist = dist,
    alg.info = list(),
    cell.names = rownames(emb)
  )
  nonlinear_reduction_params[["n.neighbors"]] <- n.neighbors

  srt_integrated <- find_neighbors_and_clusters(
    srt = srt_integrated,
    reduction = NULL,
    dims_use = NULL,
    graph_prefix = "BBKNN_",
    graph_snn = "BBKNN",
    cluster_colname = "BBKNNclusters",
    HVF = HVF,
    neighbor_metric = "euclidean",
    neighbor_k = 20L,
    cluster_algorithm = cluster_algorithm,
    cluster_algorithm_index = cluster_algorithm_index,
    cluster_resolution = cluster_resolution,
    run_find_neighbors = FALSE,
    verbose = verbose
  )

  srt_integrated <- run_nonlinear_reduction(
    srt = srt_integrated,
    prefix = "BBKNN",
    graph_use = "BBKNN",
    neighbor_use = "BBKNN_neighbors",
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed,
    verbose = verbose
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[[
    "BBKNN_HVF"
  ]] <- HVF
  srt_integrated@misc[["BBKNN_backend"]] <- bbknn_backend
  srt_integrated@misc[["BBKNN_parameters"]] <- bem[[3]]

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_integrated,
      pattern = paste0(assay, "|BBKNN|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_integrated)
  }
}

bbknn_native_graph <- function(
  pca,
  batch_list,
  neighbors_within_batch = 3,
  n_pcs = 50,
  trim = NULL,
  metric = "euclidean",
  ...,
  verbose = TRUE
) {
  unsupported <- names(list(...))
  unsupported <- unsupported[nzchar(unsupported)]
  if (length(unsupported) > 0L) {
    log_message(
      "{.arg bbknn_backend = 'r'} ignores unsupported {.arg bbknn_params}: {.val {unsupported}}",
      message_type = "warning",
      verbose = verbose
    )
  }
  pca <- as.matrix(pca)
  storage.mode(pca) <- "double"
  batch_list <- as.character(batch_list)
  if (nrow(pca) != length(batch_list)) {
    log_message(
      "{.arg pca} and {.arg batch_list} must contain the same number of cells",
      message_type = "error"
    )
  }
  n_pcs <- min(as.integer(n_pcs), ncol(pca))
  pca <- pca[, seq_len(n_pcs), drop = FALSE]
  neighbors_within_batch <- as.integer(neighbors_within_batch)
  batches <- sort(unique(batch_list))
  batch_counts <- table(batch_list)
  if (min(batch_counts[batches]) < neighbors_within_batch) {
    log_message(
      "Not all batches have at least {.arg neighbors_within_batch} cells",
      message_type = "error"
    )
  }
  metric_use <- switch(metric,
    euclidean = "euclidean",
    cosine = "cosine",
    angular = "cosine",
    {
      log_message(
        "{.arg bbknn_backend = 'r'} supports {.val euclidean}, {.val cosine}, and {.val angular}; using {.val euclidean} instead of {.val {metric}}",
        message_type = "warning",
        verbose = verbose
      )
      "euclidean"
    }
  )

  n_cells <- nrow(pca)
  n_neighbors <- neighbors_within_batch * length(batches)
  idx <- matrix(NA_integer_, nrow = n_cells, ncol = n_neighbors)
  dist <- matrix(NA_real_, nrow = n_cells, ncol = n_neighbors)
  col_start <- 1L
  for (batch in batches) {
    batch_idx <- which(batch_list == batch)
    knn <- run_cpp_knn(
      reference = pca[batch_idx, , drop = FALSE],
      query = pca,
      k = neighbors_within_batch,
      metric = metric_use,
      exclude_self = FALSE,
      n_threads = 0L
    )
    cols <- col_start:(col_start + neighbors_within_batch - 1L)
    mapped_idx <- knn[["idx"]]
    mapped_idx[] <- batch_idx[mapped_idx]
    idx[, cols] <- mapped_idx
    dist[, cols] <- knn[["dist"]]
    col_start <- col_start + neighbors_within_batch
  }

  order_idx <- t(apply(dist, 1L, order))
  gather <- cbind(
    rep(seq_len(n_cells), each = n_neighbors),
    as.vector(t(order_idx))
  )
  idx <- matrix(idx[gather], nrow = n_cells, byrow = TRUE)
  dist <- matrix(dist[gather], nrow = n_cells, byrow = TRUE)

  distances <- Matrix::sparseMatrix(
    i = rep(seq_len(n_cells), each = n_neighbors),
    j = as.vector(t(idx)),
    x = as.vector(t(dist)),
    dims = c(n_cells, n_cells)
  )
  connectivities_directed <- Matrix::sparseMatrix(
    i = rep(seq_len(n_cells), each = n_neighbors),
    j = as.vector(t(idx)),
    x = 1 / (1 + as.vector(t(dist))),
    dims = c(n_cells, n_cells)
  )
  connectivities <- pmax(connectivities_directed, Matrix::t(connectivities_directed))
  if (is.null(trim)) {
    trim <- 10L * n_neighbors
  }
  trim <- as.integer(trim)
  if (trim > 0L) {
    connectivities <- trim_sparse_rows(connectivities, trim = trim)
    connectivities <- pmax(connectivities, Matrix::t(connectivities))
  }
  list(
    distances,
    connectivities,
    list(
      n_neighbors = n_neighbors,
      method = "native_approx",
      metric = metric_use,
      n_pcs = n_pcs,
      bbknn = list(
        trim = trim,
        computation = "r_cpp_exact_knn_approx_connectivity",
        neighbors_within_batch = neighbors_within_batch
      )
    )
  )
}

trim_sparse_rows <- function(x, trim) {
  x <- methods::as(x, "dgCMatrix")
  keep <- vector("list", nrow(x))
  x_row <- Matrix::t(x)
  for (i in seq_len(nrow(x))) {
    start <- x_row@p[[i]] + 1L
    end <- x_row@p[[i + 1L]]
    if (start > end) {
      keep[[i]] <- NULL
      next
    }
    vals <- x_row@x[start:end]
    rows <- x_row@i[start:end] + 1L
    ord <- order(vals, decreasing = TRUE)
    ord <- utils::head(ord, trim)
    keep[[i]] <- data.frame(
      i = rep(i, length(ord)),
      j = rows[ord],
      x = vals[ord]
    )
  }
  keep <- do.call(rbind, keep)
  if (is.null(keep) || nrow(keep) == 0L) {
    return(Matrix::sparseMatrix(dims = dim(x)))
  }
  Matrix::sparseMatrix(
    i = keep$i,
    j = keep$j,
    x = keep$x,
    dims = dim(x)
  )
}

#' @title The CSS integration function
#'
#' @inheritParams integration_scop
#' @param CSS_dims_use A vector specifying the dimensions returned by CSS that will be utilized for downstream cell cluster finding and nonlinear reduction.
#' If set to NULL, all the returned dimensions will be used by default.
#' @param CSS_params A list of parameters for the [simspec::cluster_sim_spectrum](https://github.com/quadbio/simspec) function.
#' Default is `list()`.
#'
#' @export
CSS_integrate <- function(
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
  CSS_dims_use = NULL,
  nonlinear_reduction = "umap",
  nonlinear_reduction_dims = c(2, 3),
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  CSS_params = list(),
  verbose = TRUE,
  seed = 11
) {
  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "{.arg srt_list} and {.arg srt_merge} were all empty",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "{.arg srt_list} and {.arg srt_merge} have different cells",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- CheckDataList(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srt_merge <- Reduce(merge, srt_list)
    SeuratObject::VariableFeatures(srt_merge) <- HVF
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    checked <- CheckDataMerge(
      srt_merge = srt_merge,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      verbose = verbose,
      seed = seed
    )
    srt_merge <- checked[["srt_merge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (length(linear_reduction) > 1) {
    log_message(
      "Only the first method in the {.arg linear_reduction} will be used",
      message_type = "warning"
    )
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srt_merge)) {
    reduc_test <- c(reduc_test, SeuratObject::Reductions(srt_merge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    log_message(
      "{.arg linear_reduction} must be one of {.val {reduc_test}}",
      message_type = "error"
    )
  }
  if (
    !is.null(linear_reduction_dims_use) &&
      max(linear_reduction_dims_use) > linear_reduction_dims
  ) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  nonlinear_reductions <- c(
    "umap",
    "umap-naive",
    "tsne",
    "dm",
    "phate",
    "pacmap",
    "trimap",
    "largevis",
    "fr"
  )
  if (any(!nonlinear_reduction %in% nonlinear_reductions)) {
    log_message(
      "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
      message_type = "error"
    )
  }
  cluster_algorithms <- c("louvain", "slm", "leiden")
  if (!cluster_algorithm %in% cluster_algorithms) {
    log_message(
      "{.arg cluster_algorithm} must be one of {.val {cluster_algorithms}}",
      message_type = "error"
    )
  }
  cluster_algorithm_index <- switch(
    EXPR = tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_r(c("quadbio/simspec", "qlcMatrix"), verbose = FALSE)
  set.seed(seed)

  if (normalization_method == "TFIDF") {
    log_message(
      "{.arg normalization_method} is {.val TFIDF}. Use {.pkg lsi} workflow..."
    )
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }
  scale_features <- rownames(
    GetAssayData5(
      srt_merge,
      layer = "scale.data",
      assay = SeuratObject::DefaultAssay(srt_merge)
    )
  )
  if (
    isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))
  ) {
    log_message("Perform ScaleData")
    srt_merge <- Seurat::ScaleData(
      object = srt_merge,
      split.by = if (isTRUE(scale_within_batch)) batch else NULL,
      assay = SeuratObject::DefaultAssay(srt_merge),
      features = HVF,
      vars.to.regress = vars_to_regress,
      model.use = regression_model,
      verbose = FALSE
    )
  }

  log_message(
    "Perform {.val {linear_reduction}} linear dimension reduction"
  )
  srt_merge <- RunDimsReduction(
    srt_merge,
    prefix = "CSS",
    features = HVF,
    assay = SeuratObject::DefaultAssay(srt_merge),
    linear_reduction = linear_reduction,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = force_linear_reduction,
    verbose = verbose,
    seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- resolve_linear_dims_use(
      srt = srt_merge,
      reduction = paste0("CSS", linear_reduction),
      normalization_method = normalization_method,
      reduction_method = linear_reduction
    )
  }

  log_message(
    "Perform {.pkg CSS} integration",
    verbose = verbose
  )
  log_message(
    "Using {.val {paste0('CSS', linear_reduction)}} ({.val {min(linear_reduction_dims_use)}}:{.val {max(linear_reduction_dims_use)}}) as input",
    verbose = verbose
  )
  srt_css <- srt_merge
  if (inherits(Seurat::GetAssay(srt_css, assay = SeuratObject::DefaultAssay(srt_css)), "Assay5")) {
    srt_css <- SeuratObject::JoinLayers(srt_css)
  }
  params <- list(
    object = srt_css,
    use_dr = paste0("CSS", linear_reduction),
    dims_use = linear_reduction_dims_use,
    var_genes = HVF,
    label_tag = batch,
    reduction.name = "CSS",
    reduction.key = "CSS_",
    verbose = FALSE
  )
  for (nm in names(CSS_params)) {
    params[[nm]] <- CSS_params[[nm]]
  }
  srt_integrated <- invoke_fun(
    get_namespace_fun("simspec", "cluster_sim_spectrum"),
    params
  )

  if (any(is.na(srt_integrated@reductions[["CSS"]]@cell.embeddings))) {
    log_message(
      "NA detected in the CSS embeddings. You can try to use a lower resolution value in the {.arg CSS_params}.",
      message_type = "error"
    )
  }
  if (is.null(CSS_dims_use)) {
    CSS_dims_use <- seq_len(
      ncol(srt_integrated[["CSS"]]@cell.embeddings)
    )
  }

  srt_integrated <- find_neighbors_and_clusters(
    srt = srt_integrated,
    reduction = "CSS",
    dims_use = CSS_dims_use,
    graph_prefix = "CSS_",
    graph_snn = "CSS_SNN",
    cluster_colname = "CSSclusters",
    HVF = HVF,
    neighbor_metric = neighbor_metric,
    neighbor_k = neighbor_k,
    cluster_algorithm = cluster_algorithm,
    cluster_algorithm_index = cluster_algorithm_index,
    cluster_resolution = cluster_resolution,
    verbose = verbose
  )

  srt_integrated <- run_nonlinear_reduction(
    srt = srt_integrated,
    prefix = "CSS",
    reduction_use = "CSS",
    reduction_dims = CSS_dims_use,
    graph_use = "CSS_SNN",
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed,
    verbose = verbose
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[[
    "CSS_HVF"
  ]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_integrated,
      pattern = paste0(assay, "|CSS|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_integrated)
  }
}

#' @title The LIGER integration function
#'
#' @md
#' @inheritParams integration_scop
#' @param liger_dims_use A vector specifying the dimensions returned by LIGER that will be utilized for downstream cell cluster finding and nonlinear reduction.
#' If set to NULL, all the returned dimensions will be used by default.
#' @param optimizeALS_params A list of parameters for the [rliger::runIntegration] function.
#' Default is `list()`.
#' @param quantilenorm_params A list of parameters for the [rliger::quantileNorm] function.
#' Default is `list()`.
#'
#' @export
#'
#' @examples
#' data(panc8_sub)
#' panc8_sub <- LIGER_integrate(
#'   panc8_sub,
#'   batch = "tech"
#' )
#' CellDimPlot(
#'   panc8_sub,
#'   group.by = c("tech", "celltype")
#' )
LIGER_integrate <- function(
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
  liger_dims_use = NULL,
  nonlinear_reduction = "umap",
  nonlinear_reduction_dims = c(2, 3),
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  optimizeALS_params = list(),
  quantilenorm_params = list(),
  verbose = TRUE,
  seed = 11
) {
  nonlinear_reductions <- c(
    "umap",
    "umap-naive",
    "tsne",
    "dm",
    "phate",
    "pacmap",
    "trimap",
    "largevis",
    "fr"
  )
  if (any(!nonlinear_reduction %in% nonlinear_reductions)) {
    log_message(
      "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
      message_type = "error"
    )
  }
  cluster_algorithms <- c("louvain", "slm", "leiden")
  if (!cluster_algorithm %in% cluster_algorithms) {
    log_message(
      "{.arg cluster_algorithm} must be one of {.val {cluster_algorithms}}",
      message_type = "error"
    )
  }
  cluster_algorithm_index <- switch(
    EXPR = tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_r("rliger", verbose = FALSE)
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "{.arg srt_list} and {.arg srt_merge} were all empty",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "{.arg srt_list} and {.arg srt_merge} have different cells",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- CheckDataList(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srt_merge <- Reduce(merge, srt_list)
    SeuratObject::VariableFeatures(srt_merge) <- HVF
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    checked <- CheckDataMerge(
      srt_merge = srt_merge,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      verbose = verbose,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    srt_merge <- checked[["srt_merge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (min(sapply(srt_list, ncol)) < 30) {
    log_message(
      "The cell count in some batches is lower than 30, which may not be suitable for the current integration method",
      message_type = "warning",
      verbose = verbose
    )
    answer <- log_message(
      "Are you sure to continue?",
      message_type = "ask"
    )
    if (isFALSE(answer)) {
      return(srt_merge)
    }
  }

  SeuratObject::VariableFeatures(srt_merge) <- HVF
  liger_scale_features <- tryCatch(
    rownames(
      GetAssayData5(
        object = srt_merge,
        layer = "ligerScaleData",
        assay = SeuratObject::DefaultAssay(srt_merge)
      )
    ),
    error = function(e) character(0)
  )
  if (isFALSE(do_scaling) && length(liger_scale_features) == 0) {
    log_message(
      "When {.arg do_scaling} is FALSE, the layer {.val ligerScaleData} must already exist",
      message_type = "error"
    )
  }
  if (
    isTRUE(do_scaling) ||
      (is.null(do_scaling) && any(!HVF %in% liger_scale_features))
  ) {
    log_message(
      "Prepare {.pkg rliger} layer {.val ligerScaleData} ...",
      verbose = verbose
    )
    srt_merge <- invoke_fun(
      rliger::scaleNotCenter,
      list(
        object = srt_merge,
        assay = SeuratObject::DefaultAssay(srt_merge),
        layer = "data",
        save = "ligerScaleData",
        datasetVar = batch,
        features = HVF
      )
    )
  }

  log_message(
    "Perform {.pkg LIGER} integration",
    verbose = verbose
  )
  params1 <- list(
    object = srt_merge,
    k = 20,
    method = "iNMF",
    datasetVar = batch,
    useLayer = "ligerScaleData",
    assay = SeuratObject::DefaultAssay(srt_merge),
    seed = seed,
    verbose = FALSE
  )
  for (nm in names(optimizeALS_params)) {
    params1[[nm]] <- optimizeALS_params[[nm]]
  }
  srt_merge <- invoke_fun(rliger::runIntegration, params1)

  reduction1 <- Embeddings(object = srt_merge[["inmf"]])
  colnames(reduction1) <- paste0("riNMF_", seq_len(ncol(reduction1)))
  loadings1 <- SeuratObject::Loadings(object = srt_merge[["inmf"]])
  if (ncol(loadings1) == ncol(reduction1)) {
    colnames(loadings1) <- colnames(reduction1)
  }
  srt_merge[["iNMF_raw"]] <- CreateDimReducObject(
    embeddings = reduction1,
    loadings = loadings1,
    assay = SeuratObject::DefaultAssay(srt_merge),
    key = "riNMF_"
  )

  ref_dataset <- names(
    sort(table(srt_merge[[batch]][, 1]), decreasing = TRUE)
  )[1]
  params2 <- list(
    object = srt_merge,
    reduction = "inmf",
    reference = ref_dataset,
    useDims = seq_len(ncol(reduction1)),
    verbose = FALSE
  )
  for (nm in names(quantilenorm_params)) {
    params2[[nm]] <- quantilenorm_params[[nm]]
  }
  srt_merge <- invoke_fun(rliger::quantileNorm, params2)
  srt_merge[["LIGER"]] <- CreateDimReducObject(
    embeddings = Embeddings(object = srt_merge[["inmfNorm"]]),
    assay = SeuratObject::DefaultAssay(srt_merge),
    key = "LIGER_"
  )
  srt_integrated <- srt_merge
  srt_merge <- NULL
  if (is.null(liger_dims_use)) {
    liger_dims_use <- seq_len(
      ncol(srt_integrated[["LIGER"]]@cell.embeddings)
    )
  }

  srt_integrated <- find_neighbors_and_clusters(
    srt = srt_integrated,
    reduction = "LIGER",
    dims_use = liger_dims_use,
    graph_prefix = "LIGER_",
    graph_snn = "LIGER_SNN",
    cluster_colname = "LIGERclusters",
    HVF = HVF,
    neighbor_metric = neighbor_metric,
    neighbor_k = neighbor_k,
    cluster_algorithm = cluster_algorithm,
    cluster_algorithm_index = cluster_algorithm_index,
    cluster_resolution = cluster_resolution,
    verbose = verbose
  )

  srt_integrated <- run_nonlinear_reduction(
    srt = srt_integrated,
    prefix = "LIGER",
    reduction_use = "LIGER",
    reduction_dims = liger_dims_use,
    graph_use = "LIGER_SNN",
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed,
    verbose = verbose
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[[
    "LIGER_HVF"
  ]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_integrated,
      pattern = paste0(assay, "|LIGER|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_integrated)
  }
}

#' @title The Conos integration function
#'
#' @inheritParams integration_scop
#' @param buildGraph_params A list of parameters for the buildGraph function.
#' Default is `list()`.
#' @param cores  An integer setting the number of threads for `Conos`.
#' Default is `2`.
#'
#' @export
Conos_integrate <- function(
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
  linear_reduction = "pca",
  linear_reduction_dims = 50,
  linear_reduction_dims_use = NULL,
  linear_reduction_params = list(),
  force_linear_reduction = FALSE,
  nonlinear_reduction = "umap",
  nonlinear_reduction_dims = c(2, 3),
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  buildGraph_params = list(),
  cores = 2,
  verbose = TRUE,
  seed = 11
) {
  if (length(linear_reduction) > 1) {
    log_message(
      "Only the first method in the {.arg linear_reduction} will be used",
      message_type = "warning",
      verbose = verbose
    )
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srt_merge)) {
    reduc_test <- c(reduc_test, SeuratObject::Reductions(srt_merge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    log_message(
      "{.arg linear_reduction} must be one of {.val {reduc_test}}",
      message_type = "error"
    )
  }
  if (
    !is.null(linear_reduction_dims_use) &&
      max(linear_reduction_dims_use) > linear_reduction_dims
  ) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  nonlinear_reductions <- c("umap", "umap-naive", "fr")
  if (any(!nonlinear_reduction %in% nonlinear_reductions)) {
    log_message(
      "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
      message_type = "error"
    )
  }
  cluster_algorithms <- c("louvain", "slm", "leiden")
  if (!cluster_algorithm %in% cluster_algorithms) {
    log_message(
      "{.arg cluster_algorithm} must be one of {.val {cluster_algorithms}}",
      message_type = "error"
    )
  }
  cluster_algorithm_index <- switch(
    EXPR = tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_r("conos", verbose = FALSE)
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "{.arg srt_list} and {.arg srt_merge} were all empty",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "{.arg srt_list} and {.arg srt_merge} have different cells",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- CheckDataList(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      verbose = verbose,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srt_merge <- Reduce(merge, srt_list)
    SeuratObject::VariableFeatures(srt_merge) <- HVF
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    checked <- CheckDataMerge(
      srt_merge = srt_merge,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      verbose = verbose,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    srt_merge <- checked[["srt_merge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (min(sapply(srt_list, ncol)) < 30) {
    log_message(
      "The cell count in some batches is lower than {.val 30}, which may not be suitable for the {.pkg Conos} integration method",
      message_type = "warning"
    )
    answer <- utils::askYesNo(
      "Are you sure to continue?",
      default = FALSE
    )
    if (isFALSE(answer)) {
      return(srt_merge)
    }
  }

  srt_integrated <- srt_merge
  srt_merge <- NULL

  if (normalization_method == "TFIDF") {
    log_message(
      "{.arg normalization_method} is {.val TFIDF}. Use {.pkg lsi} workflow..."
    )
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  for (i in seq_along(srt_list)) {
    srt <- srt_list[[i]]
    scale_features <- rownames(
      GetAssayData5(
        srt,
        layer = "scale.data",
        assay = SeuratObject::DefaultAssay(srt)
      )
    )
    if (
      isTRUE(do_scaling) ||
        (is.null(do_scaling) && any(!HVF %in% scale_features))
    ) {
      log_message(
        "Perform ScaleData on the data {.val {i}} ..."
      )
      srt <- Seurat::ScaleData(
        object = srt,
        assay = SeuratObject::DefaultAssay(srt),
        features = HVF,
        vars.to.regress = vars_to_regress,
        model.use = regression_model,
        verbose = FALSE
      )
    }
    log_message(
      "Perform {.val {linear_reduction}} linear dimension reduction"
    )
    srt <- RunDimsReduction(
      srt,
      prefix = "Conos",
      features = HVF,
      assay = SeuratObject::DefaultAssay(srt),
      linear_reduction = linear_reduction,
      linear_reduction_dims = linear_reduction_dims,
      linear_reduction_params = linear_reduction_params,
      force_linear_reduction = force_linear_reduction,
      verbose = verbose,
      seed = seed
    )
    srt[["pca"]] <- srt[[paste0("Conos", linear_reduction)]]
    srt_list[[i]] <- srt
  }
  if (is.null(names(srt_list))) {
    names(srt_list) <- paste0("srt_", seq_along(srt_list))
  }

  if (is.null(linear_reduction_dims_use)) {
    maxdims <- max(unlist(sapply(
      srt_list,
      function(srt) {
        max(RunDimsEstimate(
          srt = srt,
          reduction = paste0("Conos", linear_reduction),
          reduction_method = linear_reduction,
          skip_first = normalization_method == "TFIDF",
          use_stored = TRUE,
          verbose = FALSE
        ))
      }
    )))
  } else {
    maxdims <- max(linear_reduction_dims_use)
  }

  log_message(
    " Perform {.pkg Conos} integration"
  )
  log_message(
    "{.pkg Conos} integration using {.pkg {linear_reduction}} ({.val {1}}:{.val {maxdims}}) as input",
    verbose = verbose
  )
  srt_list_con <- NULL
  conos_fun <- get_namespace_fun("conos", "Conos")$new
  invisible(
    utils::capture.output(
      srt_list_con <- suppressWarnings(
        suppressMessages(
          conos_fun(
            srt_list,
            n.cores = cores,
            verbose = FALSE
          )
        )
      ),
      type = "output"
    )
  )
  params <- list(
    ncomps = maxdims,
    verbose = FALSE
  )
  for (nm in names(buildGraph_params)) {
    params[[nm]] <- buildGraph_params[[nm]]
  }
  invisible(
    utils::capture.output(
      suppressWarnings(
        suppressMessages(
          invoke_fun(srt_list_con[["buildGraph"]], params)
        )
      ),
      type = "output"
    )
  )
  conos_graph <- igraph::as_adjacency_matrix(
    srt_list_con$graph,
    type = "both",
    attr = "weight",
    names = TRUE,
    sparse = TRUE
  )
  graph_cells <- colnames(conos_graph)
  object_cells <- colnames(srt_integrated)
  if (!setequal(graph_cells, object_cells)) {
    log_message(
      "Cell names in {.pkg Conos} graph do not match {.arg srt_integrated}",
      message_type = "error"
    )
  }
  conos_graph <- conos_graph[object_cells, object_cells, drop = FALSE]
  conos_graph <- SeuratObject::as.Graph(conos_graph)
  conos_graph@assay.used <- SeuratObject::DefaultAssay(srt_integrated)
  srt_integrated@graphs[["Conos"]] <- conos_graph
  nonlinear_reduction_params[["n.neighbors"]] <- params[["k"]]

  srt_integrated <- find_neighbors_and_clusters(
    srt = srt_integrated,
    reduction = NULL,
    dims_use = NULL,
    graph_prefix = "Conos_",
    graph_snn = "Conos",
    cluster_colname = "Conosclusters",
    HVF = HVF,
    neighbor_metric = "euclidean",
    neighbor_k = 20L,
    cluster_algorithm = cluster_algorithm,
    cluster_algorithm_index = cluster_algorithm_index,
    cluster_resolution = cluster_resolution,
    run_find_neighbors = FALSE,
    verbose = verbose
  )

  srt_integrated <- run_nonlinear_reduction(
    srt = srt_integrated,
    prefix = "Conos",
    graph_use = "Conos",
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed,
    verbose = verbose
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[[
    "Conos_HVF"
  ]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_integrated,
      pattern = paste0(assay, "|Conos|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_integrated)
  }
}

#' @title The ComBat integration function
#'
#' @inheritParams integration_scop
#' @param ComBat_params A list of parameters for the sva::ComBat function.
#' Default is `list()`.
#'
#' @export
ComBat_integrate <- function(
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
  ComBat_params = list(),
  verbose = TRUE,
  seed = 11
) {
  if (length(linear_reduction) > 1) {
    log_message(
      "Only the first method in the {.arg linear_reduction} will be used",
      message_type = "warning",
      verbose = verbose,
    )
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srt_merge)) {
    reduc_test <- c(reduc_test, SeuratObject::Reductions(srt_merge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    log_message(
      "{.arg linear_reduction} must be one of {.val {reduc_test}}",
      message_type = "error"
    )
  }
  if (
    !is.null(linear_reduction_dims_use) &&
      max(linear_reduction_dims_use) > linear_reduction_dims
  ) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  nonlinear_reductions <- c(
    "umap",
    "umap-naive",
    "tsne",
    "dm",
    "phate",
    "pacmap",
    "trimap",
    "largevis",
    "fr"
  )
  if (any(!nonlinear_reduction %in% nonlinear_reductions)) {
    log_message(
      "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
      message_type = "error"
    )
  }
  cluster_algorithms <- c("louvain", "slm", "leiden")
  if (!cluster_algorithm %in% cluster_algorithms) {
    log_message(
      "{.arg cluster_algorithm} must be one of {.val {cluster_algorithms}}",
      message_type = "error"
    )
  }
  cluster_algorithm_index <- switch(
    EXPR = tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_r("sva", verbose = FALSE)
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "{.arg srt_list} and {.arg srt_merge} were all empty",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "{.arg srt_list} and {.arg srt_merge} have different cells",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- CheckDataList(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      verbose = verbose,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srt_merge <- Reduce(merge, srt_list)
    SeuratObject::VariableFeatures(srt_merge) <- HVF
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    checked <- CheckDataMerge(
      srt_merge = srt_merge,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      verbose = verbose,
      seed = seed
    )
    srt_merge <- checked[["srt_merge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    log_message(
      "{.arg normalization_method} is {.pkg TFIDF}. Use {.pkg lsi} workflow..."
    )
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  log_message("Perform {.pkg ComBat} integration")
  dat <- GetAssayData5(
    srt_merge,
    layer = "data",
    assay = SeuratObject::DefaultAssay(srt_merge)
  )[HVF, , drop = FALSE]
  batch <- srt_merge[[batch, drop = TRUE]]
  params <- list(
    dat = dat,
    batch = batch
  )
  for (nm in names(ComBat_params)) {
    params[[nm]] <- ComBat_params[[nm]]
  }
  corrected <- suppressMessages(
    invoke_fun(
      sva::ComBat,
      params
    )
  )

  srt_integrated <- srt_merge
  srt_merge <- NULL
  srt_integrated[["ComBatcorrected"]] <- CreateAssayObject(data = corrected)
  SeuratObject::DefaultAssay(srt_integrated) <- "ComBatcorrected"
  SeuratObject::VariableFeatures(srt_integrated[["ComBatcorrected"]]) <- HVF
  scale_features <- rownames(
    GetAssayData5(
      srt_integrated,
      layer = "scale.data",
      assay = SeuratObject::DefaultAssay(srt_integrated)
    )
  )
  if (
    isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))
  ) {
    log_message("Perform {.fn Seurat::ScaleData}")
    srt_integrated <- Seurat::ScaleData(
      srt_integrated,
      split.by = if (isTRUE(scale_within_batch)) batch else NULL,
      assay = SeuratObject::DefaultAssay(srt_integrated),
      features = HVF,
      vars.to.regress = vars_to_regress,
      model.use = regression_model,
      verbose = FALSE
    )
  }

  log_message(
    "Perform linear dimension reduction ({.val {linear_reduction}})"
  )
  srt_integrated <- RunDimsReduction(
    srt_integrated,
    prefix = "ComBat",
    features = HVF,
    assay = SeuratObject::DefaultAssay(srt_integrated),
    linear_reduction = linear_reduction,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = force_linear_reduction,
    verbose = verbose,
    seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- resolve_linear_dims_use(
      srt = srt_integrated,
      reduction = paste0("ComBat", linear_reduction),
      normalization_method = normalization_method,
      reduction_method = linear_reduction
    )
  }

  srt_integrated <- find_neighbors_and_clusters(
    srt = srt_integrated,
    reduction = paste0("ComBat", linear_reduction),
    dims_use = linear_reduction_dims_use,
    graph_prefix = "ComBat_",
    graph_snn = "ComBat_SNN",
    cluster_colname = "ComBatclusters",
    HVF = HVF,
    neighbor_metric = neighbor_metric,
    neighbor_k = neighbor_k,
    cluster_algorithm = cluster_algorithm,
    cluster_algorithm_index = cluster_algorithm_index,
    cluster_resolution = cluster_resolution,
    verbose = verbose
  )

  srt_integrated <- run_nonlinear_reduction(
    srt = srt_integrated,
    prefix = "ComBat",
    reduction_use = paste0("ComBat", linear_reduction),
    reduction_dims = linear_reduction_dims_use,
    graph_use = "ComBat_SNN",
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed,
    verbose = verbose
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[[
    "ComBat_HVF"
  ]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_integrated,
      pattern = paste0(assay, "|ComBat|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_integrated)
  }
}

#' @title The Coralysis integration function
#'
#' @inheritParams integration_scop
#' @inheritParams CCA_integrate
#' @param coralysis_dims_use A vector specifying the dimensions returned by Coralysis PCA
#' that will be utilized for downstream cell cluster finding and nonlinear reduction.
#' If `NULL`, all available Coralysis PCA dimensions will be used by default.
#' @param cores Number of threads passed to [Coralysis::RunParallelDivisiveICP].
#' @param PrepareData_params A list of parameters for [Coralysis::PrepareData].
#' Default is `list()`.
#' @param RunParallelDivisiveICP_params A list of parameters for
#' [Coralysis::RunParallelDivisiveICP]. Default is `list()`.
#' @param RunPCA_params A list of parameters for [Coralysis::RunPCA].
#' Default is `list()`.
#'
#' @export
Coralysis_integrate <- function(
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
  coralysis_dims_use = NULL,
  nonlinear_reduction = "umap",
  nonlinear_reduction_dims = c(2, 3),
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  cores = NULL,
  PrepareData_params = list(),
  RunParallelDivisiveICP_params = list(),
  RunPCA_params = list(),
  verbose = TRUE,
  seed = 11
) {
  check_r("Coralysis", verbose = FALSE)

  if (!identical(normalization_method, "LogNormalize")) {
    log_message(
      "{.pkg Coralysis} requires {.arg normalization_method = 'LogNormalize'}",
      message_type = "error"
    )
  }

  nonlinear_reductions <- c(
    "umap",
    "umap-naive",
    "tsne",
    "dm",
    "phate",
    "pacmap",
    "trimap",
    "largevis",
    "fr"
  )
  if (any(!nonlinear_reduction %in% nonlinear_reductions)) {
    log_message(
      "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
      message_type = "error"
    )
  }
  cluster_algorithms <- c("louvain", "slm", "leiden")
  if (!cluster_algorithm %in% cluster_algorithms) {
    log_message(
      "{.arg cluster_algorithm} must be one of {.val {cluster_algorithms}}",
      message_type = "error"
    )
  }
  cluster_algorithm_index <- switch(
    EXPR = tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "{.arg srt_list} and {.arg srt_merge} were all empty",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "{.arg srt_list} and {.arg srt_merge} have different cells",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- CheckDataList(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = NULL,
      verbose = verbose,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    srt_merge <- Reduce(merge, srt_list)
    SeuratObject::VariableFeatures(srt_merge) <- HVF
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    checked <- CheckDataMerge(
      srt_merge = srt_merge,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = NULL,
      verbose = verbose,
      seed = seed
    )
    srt_merge <- checked[["srt_merge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
  }

  assay_use <- SeuratObject::DefaultAssay(srt_merge)
  if (!is.null(assay) && assay %in% SeuratObject::Assays(srt_merge)) {
    SeuratObject::DefaultAssay(srt_merge) <- assay
    assay_use <- assay
  }

  log_message(
    "Perform {.pkg Coralysis} integration",
    verbose = verbose
  )
  srt_sce <- srt_merge
  if (inherits(Seurat::GetAssay(srt_sce, assay = assay_use), "Assay5")) {
    srt_sce[[assay_use]] <- SeuratObject::JoinLayers(srt_sce[[assay_use]])
  }
  sce <- Seurat::as.SingleCellExperiment(srt_sce, assay = assay_use)
  if (!is.null(HVF)) {
    HVF <- intersect(HVF, rownames(sce))
  }
  if (length(HVF) == 0L) {
    log_message(
      "No highly variable features were available for {.pkg Coralysis}",
      message_type = "error"
    )
  }
  sce <- sce[HVF, , drop = FALSE]

  assay_names <- SummarizedExperiment::assayNames(sce)
  if (!"logcounts" %in% assay_names) {
    log_message(
      "{.pkg Coralysis} requires a {.val logcounts} assay after conversion from {.cls Seurat}",
      message_type = "error"
    )
  }

  prep_params <- utils::modifyList(
    list(object = sce),
    PrepareData_params
  )
  sce <- invoke_fun(Coralysis::PrepareData, prep_params)

  run_params <- utils::modifyList(
    list(
      object = sce,
      batch.label = batch
    ),
    RunParallelDivisiveICP_params
  )
  if (is.null(run_params[["threads"]]) && !is.null(cores)) {
    run_params[["threads"]] <- as.integer(cores)
  }
  set.seed(seed)
  sce <- invoke_fun(Coralysis::RunParallelDivisiveICP, run_params)

  pca_params <- utils::modifyList(
    list(
      object = sce,
      assay.name = "joint.probability",
      dimred.name = "Coralysis"
    ),
    RunPCA_params
  )
  set.seed(seed)
  sce <- invoke_fun(Coralysis::RunPCA, pca_params)

  if (!"Coralysis" %in% SingleCellExperiment::reducedDimNames(sce)) {
    log_message(
      "{.pkg Coralysis} did not return the expected {.val Coralysis} reduced dimension",
      message_type = "error"
    )
  }

  coralysis_emb <- SingleCellExperiment::reducedDim(sce, "Coralysis")
  coralysis_emb <- as.matrix(coralysis_emb)
  coralysis_emb <- coralysis_emb[colnames(srt_merge), , drop = FALSE]
  colnames(coralysis_emb) <- paste0("Coralysis_", seq_len(ncol(coralysis_emb)))

  srt_integrated <- srt_merge
  srt_integrated[["Coralysis"]] <- Seurat::CreateDimReducObject(
    embeddings = coralysis_emb,
    assay = assay_use,
    key = "Coralysis_"
  )

  dims_use <- coralysis_dims_use
  if (is.null(dims_use)) {
    dims_use <- seq_len(ncol(coralysis_emb))
  }
  dims_use <- sort(unique(as.integer(dims_use)))
  dims_use <- dims_use[!is.na(dims_use) & dims_use >= 1L]
  if (length(dims_use) == 0L) {
    log_message(
      "{.arg coralysis_dims_use} must contain positive integers",
      message_type = "error"
    )
  }
  if (max(dims_use) > ncol(coralysis_emb)) {
    log_message(
      "{.arg coralysis_dims_use} exceeds the number of available Coralysis dimensions",
      message_type = "error"
    )
  }

  srt_integrated <- find_neighbors_and_clusters(
    srt = srt_integrated,
    reduction = "Coralysis",
    dims_use = dims_use,
    graph_prefix = "Coralysis_",
    graph_snn = "Coralysis_SNN",
    cluster_colname = "Coralysisclusters",
    HVF = HVF,
    neighbor_metric = neighbor_metric,
    neighbor_k = neighbor_k,
    cluster_algorithm = cluster_algorithm,
    cluster_algorithm_index = cluster_algorithm_index,
    cluster_resolution = cluster_resolution,
    verbose = verbose
  )

  srt_integrated <- run_nonlinear_reduction(
    srt = srt_integrated,
    prefix = "Coralysis",
    reduction_use = "Coralysis",
    reduction_dims = dims_use,
    graph_use = "Coralysis_SNN",
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed,
    verbose = verbose
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay_use
  SeuratObject::VariableFeatures(srt_integrated) <- HVF
  srt_integrated@misc[["Coralysis_HVF"]] <- HVF
  srt_integrated@misc[["Coralysis_clusters"]] <- as.vector(
    SummarizedExperiment::colData(sce)[[run_params[["label.name"]] %||% "cluster"]]
  )

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_output <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_integrated,
      pattern = paste0(assay_use, "|Coralysis|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    SeuratObject::DefaultAssay(srt_output) <- assay_use
    SeuratObject::VariableFeatures(srt_output) <- HVF
    srt_output@misc[["Coralysis_HVF"]] <- HVF
    srt_output@misc[["Coralysis_clusters"]] <- srt_integrated@misc[[
      "Coralysis_clusters"
    ]]
    return(srt_output)
  }

  srt_integrated
}
