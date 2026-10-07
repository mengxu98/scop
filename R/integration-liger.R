#' @title LIGER integration function
#'
#' @md
#' @inheritParams RunIntegration
#' @param liger_dims_use Dimensions returned by LIGER that will be utilized for downstream cell cluster finding and nonlinear reduction.
#' If set to NULL, all the returned dimensions will be used by default.
#' @param optimizeALS_params A list of parameters for the [rliger::runIntegration] function.
#' @param quantilenorm_params A list of parameters for the [rliger::quantileNorm] function.
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
  validate_nonlinear_reductions(nonlinear_reduction)
  cluster_algorithm_index <- resolve_cluster_algorithm_index(cluster_algorithm)

  check_r("rliger", verbose = FALSE)
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
