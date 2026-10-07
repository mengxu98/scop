#' @title MNN integration function
#'
#' @inheritParams RunIntegration
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
  validate_nonlinear_reductions(nonlinear_reduction)
  cluster_algorithm_index <- resolve_cluster_algorithm_index(cluster_algorithm)

  check_r("batchelor", verbose = FALSE)
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
    srt_integrated <- ScaleData(
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


#' @title FastMNN integration function
#'
#' @inheritParams RunIntegration
#' @param fastMNN_dims_use Dimensions returned by fastMNN that will be utilized for downstream cell cluster finding and nonlinear reduction.
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
  cluster_algorithm_index <- resolve_cluster_algorithm_index(cluster_algorithm)

  check_r("batchelor", verbose = FALSE)
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
  corrected_matrix <- Matrix::Matrix(
    as.matrix(out@assays@data$reconstructed),
    sparse = TRUE
  )
  srt_integrated[["fastMNNcorrected"]] <- CreateAssayObject(
    counts = corrected_matrix
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
