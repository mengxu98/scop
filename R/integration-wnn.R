#' @title WNN integration function
#'
#' @inheritParams RunIntegration
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
  srt_merge <- FindMultiModalNeighbors(
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
          object = srt,
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
          object = srt,
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
