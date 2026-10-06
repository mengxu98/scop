#' @title ScVI integration function
#'
#' @inheritParams RunIntegration
#' @param scVI_dims_use Dimensions returned by scVI that will be utilized for downstream cell cluster finding and nonlinear reduction.
#' If set to NULL, all the returned dimensions will be used by default.
#' @param model A string indicating the scVI model to be used.
#' Options are "SCVI", "PEAKVI", and "POISSONVI".
#' @param SCVI_params A list of parameters for the SCVI model.
#' @param PEAKVI_params A list of parameters for the PEAKVI model.
#' @param POISSONVI_params A list of parameters for the POISSONVI model.
#' @param train_params A list of parameters passed to the model `train()` method.
#' @param cores An integer setting the number of threads for `scVI`.
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
  validate_nonlinear_reductions(nonlinear_reduction)
  cluster_algorithm_index <- resolve_cluster_algorithm_index(cluster_algorithm)

  PrepareEnv(modules = "scvi")
  check_python("scvi-tools")
  scvi <- reticulate::import("scvi")
  scipy <- reticulate::import("scipy")
  set.seed(seed)

  scvi$settings$num_threads <- as.integer(cores)

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
