#' @title CSS integration function
#'
#' @inheritParams RunIntegration
#' @param CSS_dims_use Dimensions returned by CSS that will be utilized for downstream cell cluster finding and nonlinear reduction.
#' If set to NULL, all the returned dimensions will be used by default.
#' @param CSS_params A list of parameters for the [simspec::cluster_sim_spectrum](https://github.com/quadbio/simspec) function.
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
  validate_nonlinear_reductions(nonlinear_reduction)
  cluster_algorithm_index <- resolve_cluster_algorithm_index(cluster_algorithm)

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
    assay_merge <- SeuratObject::DefaultAssay(srt_merge)
    if (inherits(srt_merge[[assay_merge]], "Assay5")) {
      srt_merge[[assay_merge]] <- SeuratObject::JoinLayers(
        srt_merge[[assay_merge]]
      )
    }
    srt_merge <- ScaleData(
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
