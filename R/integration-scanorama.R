#' @title Scanorama integration function
#'
#' @inheritParams RunIntegration
#' @param Scanorama_dims_use Dimensions returned by Scanorama that will be utilized for downstream cell cluster finding and nonlinear reduction.
#' If set to NULL, all the returned dimensions will be used by default.
#' @param return_corrected Whether to return the corrected data.
#' @param Scanorama_params A list of parameters for the scanorama.correct function.
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
  cluster_algorithm_index <- resolve_cluster_algorithm_index(cluster_algorithm)

  check_python("scanorama")
  scanorama <- reticulate::import("scanorama")
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
    assaylist[[i]] <- python_cells_by_features(
      GetAssayData5(
        object = srt_list[[i]],
        layer = "data",
        assay = SeuratObject::DefaultAssay(srt_list[[i]])
      )[HVF, , drop = FALSE]
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
