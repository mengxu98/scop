#' @title Coralysis integration function
#'
#' @inheritParams RunIntegration
#' @inheritParams CCA_integrate
#' @param coralysis_dims_use Dimensions returned by Coralysis PCA
#' that will be utilized for downstream cell cluster finding and nonlinear reduction.
#' If `NULL`, all available Coralysis PCA dimensions will be used by default.
#' @param cores Number of threads passed to [Coralysis::RunParallelDivisiveICP].
#' @param PrepareData_params A list of parameters for [Coralysis::PrepareData].
#' @param RunParallelDivisiveICP_params A list of parameters for
#' [Coralysis::RunParallelDivisiveICP].
#' @param RunPCA_params A list of parameters for [Coralysis::RunPCA].
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

  validate_nonlinear_reductions(nonlinear_reduction)
  cluster_algorithm_index <- resolve_cluster_algorithm_index(cluster_algorithm)

  set.seed(seed)
  validate_integration_input_cells(srt_list, srt_merge)
  srt_merge_raw <- srt_merge
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
