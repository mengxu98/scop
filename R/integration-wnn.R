#' @title WNN integration function
#'
#' @md
#' @inheritParams RunIntegration
#'
#' @description Combine paired RNA and ATAC representations with weighted
#' nearest neighbors, retaining both modality reductions and the WNN graph.
#' @param assay One RNA assay and one chromatin assay, or one assay whose
#' partner can be inferred. Automatic selection prefers an active RNA/SCT
#' assay, then RNA, then SCT; custom RNA assay names must be explicit.
#' Known SCOP gene-activity and chromVAR assays cannot substitute for RNA.
#' @param batch Accepted for integration-interface compatibility. WNN does not
#' perform batch correction; supplying this argument emits a warning.
#' @param linear_reduction RNA reduction, which must be `"pca"`. ATAC uses SVD.
#' @param nonlinear_reduction One or more of `"umap"`, `"umap-naive"` and `"fr"`,
#'   or `character()` for graphs and clusters without an embedding.
#' @param HVF Optional RNA feature names. ATAC features are selected independently.
#' @details
#' Cells must represent paired measurements. RNA and ATAC preprocessing is
#' followed by native multimodal neighbor construction. This function does not
#' use batch labels to correct either modality. Specify assays explicitly after
#' adding derived assays, for example `assay = c("RNA", "peaks")`.
#' Both assays must contain the same cells as the object. Subset to verified
#' paired cells explicitly before running WNN. Nonlinear reductions supported
#' here are `"umap"`, `"umap-naive"` and `"fr"`; use `character()` to retain
#' graphs and clusters without computing an embedding.
#' RNA uses the requested normalization and PCA. ATAC always uses TF-IDF and
#' SVD, including when RNA uses SCT. `HVF` supplies RNA features; ATAC features
#' are selected independently from accessibility counts.
#' @return A Seurat object with RNA PCA and ATAC LSI reductions, `WNN` neighbors,
#' `WNNKNN` and `WNNSNN` graphs, `WNNclusters` metadata and requested WNN
#' embeddings. `Misc(object, "WNN_reduction_list")` names the retained source
#' reductions; `WNN_dims_list` records their selected dimensions. Both return
#' modes retain complete processed modality assays (including scaled RNA), ATAC
#' reductions, metadata and representation receipts, with RNA as the active assay. The
#' caller's random-number state is restored.
#' @seealso [ReadMultiomeData()], [RunStandardWorkflow()], [RunLinkPeaks()]
#'
#' @export
#' @examples
#' # Small paired subset; no external fragments or reference downloads.
#' data("pbmcmultiome_sub", package = "scop")
#' paired <- subset(pbmcmultiome_sub, cells = colnames(pbmcmultiome_sub)[1:80])
#' paired <- WNN_integrate(
#'   srt_merge = paired, assay = c("RNA", "peaks"),
#'   linear_reduction_dims = 10, linear_reduction_dims_use = 2:8,
#'   nHVF = 100, neighbor_k = 10, nonlinear_reduction_dims = 2,
#'   verbose = FALSE
#' )
#' SeuratObject::Reductions(paired)
#' stopifnot(all(paired@misc$WNN_reduction_list %in% SeuratObject::Reductions(paired)))
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
  if (!is.character(nonlinear_reduction) || anyNA(nonlinear_reduction) ||
    any(!nonlinear_reduction %in% c("umap", "umap-naive", "fr"))) {
    log_message("WNN supports only umap, umap-naive and fr; use character() for graphs and clusters without an embedding.", message_type = "error")
  }
  if (!identical(linear_reduction, "pca")) {
    log_message("WNN requires linear_reduction = 'pca' for RNA and computes ATAC LSI separately.", message_type = "error")
  }

  atac_backend_fun("withr", "local_preserve_seed")()
  set.seed(seed)
  if (!is.null(batch)) {
    # Supported thisutils releases silence warnings when verbose is FALSE.
    log_message(
      "WNN combines modalities but does not perform batch correction; batch is not used to correct the input representations.",
      message_type = "warning", verbose = TRUE
    )
  }
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
  if (!setequal(colnames(srt_merge[[rna_assay]]), colnames(srt_merge)) ||
    !setequal(colnames(srt_merge[[atac_assay]]), colnames(srt_merge))) {
    log_message("WNN requires the same paired cells in the RNA assay, ATAC assay and object. Subset explicitly to verified paired cells before analysis.", message_type = "error")
  }
  rna_prefix <- resolve_assay_prefix(srt = srt_merge, assay = rna_assay)
  atac_prefix <- resolve_assay_prefix(srt = srt_merge, assay = atac_assay)

  for (modality in c("rna", "atac")) {
    srt_merge <- RunStandardWorkflow(
      object = srt_merge,
      prefix = if (modality == "rna") rna_prefix else atac_prefix,
      assay = assay_pair[[modality]],
      do_normalization = do_normalization,
      normalization_method = if (modality == "rna") normalization_method else "TFIDF",
      do_HVF_finding = do_HVF_finding,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF = if (modality == "rna") HVF else NULL,
      do_scaling = do_scaling,
      vars_to_regress = vars_to_regress,
      regression_model = regression_model,
      linear_reduction = if (modality == "rna") "pca" else "svd",
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
    if (modality == "rna") rna_assay <- SeuratObject::DefaultAssay(srt_merge)
  }

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
    verbose = verbose,
    stop_on_error = TRUE
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
    # Copy whole modality assays: srt_append only synchronizes selected layers.
    for (name in c(rna_assay, atac_assay)) srt_merge_raw[[name]] <- srt_merge[[name]]
    SeuratObject::DefaultAssay(srt_merge_raw) <- rna_assay
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_merge,
      pattern = paste0(
        "^(", paste(gsub(".", "[.]", unique(c(rna_assay, atac_assay, rna_prefix, atac_prefix)), fixed = TRUE), collapse = "|"),
        "|WNN|Default_reduction)"
      ),
      slots = setdiff(methods::slotNames(srt_merge), "assays"),
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
  derived <- vapply(rna_assays, function(name) {
    name %in% c("ACTIVITY", "chromvar") ||
      !is.null(SeuratObject::Misc(srt[[name]], "scop_gene_activity")) ||
      !is.null(SeuratObject::Misc(srt[[name]], "scop_chromvar"))
  }, logical(1))
  rna_assays <- rna_assays[!derived]
  if (length(chrom_assays) == 0 || length(rna_assays) == 0) {
    log_message(
      "WNN requires at least one RNA assay and one {.cls ChromatinAssay}",
      message_type = "error"
    )
  }
  if (!is.null(assay) && (!is.character(assay) || !length(assay) ||
    length(assay) > 2L || anyNA(assay) || anyDuplicated(assay) ||
    !all(assay %in% c(rna_assays, chrom_assays)))) {
    log_message(
      "assay for WNN must select measured RNA and chromatin assays; derived scores and unknown or duplicate assay names are not supported.",
      message_type = "error"
    )
  }
  rna_assay <- assay[assay %in% rna_assays]
  atac_assay <- assay[assay %in% chrom_assays]
  if (!length(rna_assay)) {
    preferred <- unique(c(SeuratObject::DefaultAssay(srt), "RNA", "SCT"))
    preferred <- preferred[preferred %in% intersect(rna_assays, c("RNA", "SCT"))]
    rna_assay <- utils::head(preferred, 1L)
  }
  if (!length(atac_assay)) {
    atac_assay <- if ("peaks" %in% chrom_assays) "peaks" else chrom_assays
  }
  if (length(rna_assay) != 1L || length(atac_assay) != 1L) {
    log_message(
      "WNN requires one RNA assay and one ChromatinAssay; specify assay = c('RNA_assay', 'ATAC_assay') when selection is ambiguous or RNA uses a custom name.",
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
  for (nr in nonlinear_reduction) {
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
