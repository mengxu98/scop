#' @title Standard single-cell and spatial workflow
#'
#' @description
#' Normalize, find variable features, scale, reduce dimensions, and cluster a
#' `Seurat` object. `workflow = "spatial"` adds spot QC, spatial variable
#' features, optional BayesSpace clustering, and optional deconvolution. It is
#' not a multi-slice integration orchestrator.
#'
#' Objects that also contain a `ChromatinAssay` are preprocessed sequentially
#' (default assay, then chromatin).
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams scop-params
#' @param prefix Prefix for intermediate object names.
#' @param workflow `"single_cell"` or `"spatial"` (basic single-image Visium-style).
#' @param do_spot_qc,spot_qc_params Run [RunSpotQC()] and extra arguments.
#' @param do_spatial_variable_features,spatial_variable_features_params
#' Run [RunSpatialVariableFeatures()]. The workflow defaults
#' `set_variable_features = FALSE` so expression HVFs are kept.
#' @param do_spatial_cluster,spatial_cluster_method,spatial_q,bayesspace_params
#' Spatial clustering. Only `"BayesSpace"` is supported. `spatial_q = NULL`
#' uses the number of ordinary spot clusters.
#' @param reference,reference_label,reference_assay Optional single-cell reference
#' for deconvolution.
#' @param do_deconvolution,deconvolution_method,deconvolution_params
#' Deconvolution (`"RCTD"`, `"SPOTlight"`, or `"Cell2location"`). `NULL` runs
#' when `reference` or cell2location signatures are provided.
#' @param do_normalization,normalization_method Normalize if the assay has no
#' scaled data (`NULL`). One of `"LogNormalize"`, `"SCT"`, `"TFIDF"`, `"scran"`.
#' `"SCT"` switches downstream steps to the `"SCT"` assay.
#' @param do_HVF_finding,HVF_method,nHVF,HVF Highly variable features.
#' `HVF_method` is `"vst"`, `"mvp"`, `"disp"`, or `"scran"`.
#' @param do_scaling Force scaling via ScaleData.
#' @param vars_to_regress Regressors used in scaling. `NULL` uses none.
#' @param regression_model `"linear"`, `"poisson"`, or `"negativebinomial"`.
#' @param linear_reduction,linear_reduction_dims,linear_reduction_dims_use,linear_reduction_params,force_linear_reduction
#' Linear reduction (`"pca"`, `"svd"`, `"ica"`, `"nmf"`, `"mds"`, `"glmpca"`).
#' `linear_reduction_dims_use = NULL` uses estimated dimensions, else the first 50.
#' @param nonlinear_reduction,nonlinear_reduction_dims,nonlinear_reduction_params,force_nonlinear_reduction
#' Nonlinear reduction (`"umap"`, `"umap-naive"`, `"tsne"`, `"dm"`, `"phate"`,
#' `"pacmap"`, `"trimap"`, `"largevis"`, `"fr"`).
#' @param neighbor_metric,neighbor_k Neighbor graph (`"euclidean"`, `"cosine"`,
#' `"manhattan"`, `"hamming"`).
#' @param cluster_algorithm,cluster_resolution Clustering (`"louvain"`, `"slm"`,
#' `"leiden"`). Larger `cluster_resolution` yields fewer clusters.
#' @param ... Additional arguments passed to reduction methods.
#'
#' @return A `Seurat` object. Spatial workflows store per-stage execution and
#' result-storage state in `srt@tools[["run_standard_spatial_workflow"]]`,
#' including the effective `set_variable_features` value, `result_stored`, and
#' `result_location`. A failed stage signals an error with the same table in
#' attribute `standard_spatial_stages`.
#'
#' @export
#'
#' @examples
#' library(Matrix)
#' data(pancreas_sub)
#' pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType"
#' )
#'
#' # Use a combination of different linear
#' # or nonlinear dimension reduction methods
#' linear_reductions <- c(
#'   "pca", "nmf", "mds"
#' )
#' pancreas_sub <- RunStandardWorkflow(
#'   pancreas_sub,
#'   linear_reduction = linear_reductions,
#'   nonlinear_reduction = "umap"
#' )
#' plist1 <- lapply(
#'   linear_reductions, function(lr) {
#'     CellDimPlot(
#'       pancreas_sub,
#'       group.by = "SubCellType",
#'       reduction = paste0(
#'         "Standard", lr, "UMAP2D"
#'       ),
#'       xlab = "", ylab = "",
#'       title = paste0(lr, "_umap"),
#'       legend.position = "none",
#'       theme_use = "theme_blank"
#'     )
#'   }
#' )
#' patchwork::wrap_plots(plist1)
#'
#' nonlinear_reductions <- c(
#'   "umap", "tsne", "fr"
#' )
#' pancreas_sub <- RunStandardWorkflow(
#'   pancreas_sub,
#'   linear_reduction = "pca",
#'   nonlinear_reduction = nonlinear_reductions
#' )
#' plist2 <- lapply(
#'   nonlinear_reductions, function(nr) {
#'     CellDimPlot(
#'       pancreas_sub,
#'       group.by = "SubCellType",
#'       reduction = paste0(
#'         "Standardpca", nr, "2D"
#'       ),
#'       xlab = "", ylab = "",
#'       title = paste0("pca_", nr),
#'       legend.position = "none",
#'       theme_use = "theme_blank"
#'     )
#'   }
#' )
#' patchwork::wrap_plots(plist2)
#'
#' data(visium_human_pancreas_sub)
#' spatial <- RunStandardWorkflow(
#'   visium_human_pancreas_sub,
#'   workflow = "spatial",
#'   assay = "Spatial",
#'   do_spatial_cluster = FALSE,
#'   spatial_cluster_method = "BayesSpace",
#'   do_deconvolution = FALSE,
#'   deconvolution_method = "RCTD",
#'   linear_reduction_dims = 10,
#'   linear_reduction_dims_use = 1:5,
#'   nonlinear_reduction_dims = 2,
#'   spatial_variable_features_params = list(nfeatures = 50)
#' )
#' SpatialSpotPlot(spatial, group.by = "SpotQC")
#' SpatialSpotPlot(
#'   spatial,
#'   features = spatial@tools[["SpatialVariableFeatures"]]$summary$top_features[1:2]
#' )
#'
#' spatial_bayes <- RunStandardWorkflow(
#'   visium_human_pancreas_sub,
#'   workflow = "spatial",
#'   assay = "Spatial",
#'   do_spatial_cluster = TRUE,
#'   spatial_cluster_method = "BayesSpace",
#'   spatial_q = 3,
#'   do_deconvolution = FALSE,
#'   deconvolution_method = "RCTD",
#'   bayesspace_params = list(
#'     n.PCs = 5,
#'     n.HVGs = 200,
#'     store_sce = FALSE,
#'     spatial_cluster_params = list(
#'       nrep = 200,
#'       burn.in = 50,
#'       thin = 10,
#'       save.chain = FALSE
#'     )
#'   )
#' )
#' SpatialSpotPlot(spatial_bayes, group.by = "BayesSpace_cluster")
RunStandardWorkflow <- function(
  srt,
  prefix = "Standard",
  workflow = c("single_cell", "spatial"),
  assay = NULL,
  image = NULL,
  coord.cols = c("x", "y"),
  do_spot_qc = TRUE,
  spot_qc_params = list(),
  do_spatial_variable_features = TRUE,
  spatial_variable_features_params = list(),
  do_spatial_cluster = FALSE,
  spatial_cluster_method = "BayesSpace",
  spatial_q = NULL,
  bayesspace_params = list(),
  reference = NULL,
  reference_label = NULL,
  reference_assay = NULL,
  do_deconvolution = !is.null(reference),
  deconvolution_method = "RCTD",
  deconvolution_params = list(),
  do_normalization = NULL,
  normalization_method = "LogNormalize",
  do_HVF_finding = TRUE,
  HVF_method = "vst",
  nHVF = 2000,
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
  nonlinear_reduction_dims = 2,
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  cores = 1L,
  verbose = TRUE,
  seed = 11,
  ...
) {
  workflow <- match.arg(workflow)
  if (identical(workflow, "spatial")) {
    return(run_standard_spatial_workflow(
      srt = srt,
      prefix = prefix,
      assay = assay,
      image = image,
      coord.cols = coord.cols,
      do_spot_qc = do_spot_qc,
      spot_qc_params = spot_qc_params,
      do_spatial_variable_features = do_spatial_variable_features,
      spatial_variable_features_params = spatial_variable_features_params,
      do_spatial_cluster = do_spatial_cluster,
      spatial_cluster_method = spatial_cluster_method,
      spatial_q = spatial_q,
      bayesspace_params = bayesspace_params,
      reference = reference,
      reference_label = reference_label,
      reference_assay = reference_assay,
      do_deconvolution = do_deconvolution,
      deconvolution_method = deconvolution_method,
      deconvolution_params = deconvolution_params,
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
      cores = cores,
      verbose = verbose,
      seed = seed,
      ...
    ))
  }

  log_message(
    "Start standard processing workflow...",
    text_color = "blue",
    verbose = verbose
  )

  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }

  assays_to_run <- resolve_standard_assays(
    srt = srt,
    assay = assay
  )
  if (length(assays_to_run) > 1) {
    assay_default <- SeuratObject::DefaultAssay(srt)
    log_message(
      "Auto preprocess assays: {.val {assays_to_run}}",
      verbose = verbose
    )
    for (i in seq_along(assays_to_run)) {
      assay_i <- assays_to_run[[i]]
      prefix_i <- resolve_standard_prefix(
        srt = srt,
        assay = assay_i,
        prefix = prefix,
        multi_assay = TRUE,
        primary_assay = assays_to_run[[1]]
      )
      srt <- RunStandardWorkflow(
        srt = srt,
        prefix = prefix_i,
        assay = assay_i,
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
    }
    SeuratObject::DefaultAssay(srt) <- assay_default
    return(srt)
  }
  assay <- assays_to_run[[1]]
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  is_atac_assay <- inherits(srt[[assay]], "ChromatinAssay")

  if (is_atac_assay) {
    atac_args <- atac_defaults(
      prefix = prefix,
      do_normalization = do_normalization,
      normalization_method = normalization_method,
      do_HVF_finding = do_HVF_finding,
      nHVF = nHVF,
      do_scaling = do_scaling,
      linear_reduction = linear_reduction,
      linear_reduction_dims_use = linear_reduction_dims_use,
      neighbor_metric = neighbor_metric
    )
    prefix <- atac_args$prefix
    do_normalization <- atac_args$do_normalization
    normalization_method <- atac_args$normalization_method
    do_HVF_finding <- atac_args$do_HVF_finding
    nHVF <- atac_args$nHVF
    do_scaling <- atac_args$do_scaling
    linear_reduction <- atac_args$linear_reduction
    linear_reduction_dims_use <- atac_args$linear_reduction_dims_use
    neighbor_metric <- atac_args$neighbor_metric
  }

  linear_reductions <- c(
    "pca",
    "svd",
    "ica",
    "nmf",
    "mds",
    "glmpca"
  )
  if (
    any(
      !linear_reduction %in% c(linear_reductions, SeuratObject::Reductions(srt))
    )
  ) {
    log_message(
      "{.arg linear_reduction} must be one of: {.val {linear_reductions}}",
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
      "{.arg nonlinear_reduction} must be one of: {.val {nonlinear_reductions}}",
      message_type = "error"
    )
  }
  cluster_algorithm_index <- resolve_cluster_algorithm_index(cluster_algorithm)

  set.seed(seed)

  checked <- CheckDataList(
    srt_list = list(srt),
    batch = "",
    assay = assay,
    do_normalization = do_normalization,
    do_HVF_finding = do_HVF_finding,
    normalization_method = normalization_method,
    HVF_source = "separate",
    HVF_method = HVF_method,
    nHVF = nHVF,
    HVF = HVF,
    vars_to_regress = vars_to_regress,
    cores = cores,
    seed = seed,
    verbose = verbose
  )
  srt <- checked[["srt_list"]][[1]]
  HVF <- checked[["HVF"]]
  assay <- checked[["assay"]]
  type <- checked[["type"]]
  rm(checked)

  if (normalization_method == "SCT" && type == "RNA") {
    assay <- "SCT"
    SeuratObject::DefaultAssay(srt) <- assay
  }

  if (normalization_method == "TFIDF") {
    log_message(
      "{.arg normalization_method} is {.pkg TFIDF}. Use {.pkg lsi} workflow",
      verbose = verbose
    )
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }
  assay_obj <- srt[[assay]]
  scale_features <- if (inherits(assay_obj, "Assay5")) {
    if ("scale.data" %in% names(assay_obj@layers)) {
      rownames(SeuratObject::GetAssayData(assay_obj, layer = "scale.data"))
    } else {
      character(0)
    }
  } else {
    sc <- assay_obj@scale.data
    if (!is.null(sc) && nrow(sc) > 0) rownames(sc) else character(0)
  }
  if (
    isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))
  ) {
    if (normalization_method != "SCT") {
      log_message(
        "Perform {.fn ScaleData}",
        verbose = verbose
      )
      srt <- suppressWarnings(
        ScaleData(
          object = srt,
          assay = assay,
          features = HVF,
          vars.to.regress = vars_to_regress,
          model.use = regression_model,
          verbose = FALSE
        )
      )
    }
  }

  for (lr in linear_reduction) {
    log_message(
      "Perform {.pkg {lr}} linear dimension reduction",
      verbose = verbose
    )
    if (
      identical(lr, "pca") &&
        isTRUE(force_linear_reduction) &&
        length(linear_reduction_params) == 0L
    ) {
      srt <- RunPCA(
        object = srt,
        assay = assay,
        features = HVF,
        npcs = linear_reduction_dims,
        reduction.name = paste0(prefix, lr),
        reduction.key = paste0(prefix, "PC_"),
        verbose = FALSE,
        seed.use = seed
      )
      srt@misc[["Default_reduction"]] <- paste0(prefix, lr)
    } else {
      srt <- RunDimsReduction(
        srt,
        prefix = prefix,
        features = HVF,
        assay = assay,
        linear_reduction = lr,
        linear_reduction_dims = linear_reduction_dims,
        linear_reduction_params = linear_reduction_params,
        force_linear_reduction = force_linear_reduction,
        verbose = verbose,
        seed = seed
      )
    }

    if (is.null(linear_reduction_dims_use)) {
      linear_reduction_dims_use_current <- RunDimsEstimate(
        srt = srt,
        reduction = paste0(prefix, lr),
        reduction_method = lr,
        skip_first = normalization_method == "TFIDF",
        use_stored = TRUE,
        verbose = verbose,
        ...
      )
    } else {
      linear_reduction_dims_use_current <- linear_reduction_dims_use
    }

    srt <- tryCatch(
      {
        srt <- FindNeighbors(
          object = srt,
          reduction = paste0(prefix, lr),
          dims = linear_reduction_dims_use_current,
          annoy.metric = neighbor_metric,
          k.param = neighbor_k,
          graph.name = paste0(prefix, lr, "_", c("KNN", "SNN")),
          verbose = FALSE
        )

        log_message(
          "Perform {.fn Seurat::FindClusters} with {.arg cluster_algorithm = '{cluster_algorithm}'} and {.arg cluster_resolution = {cluster_resolution}}",
          verbose = verbose
        )
        find_clusters_method <- utils::getS3method(
          "FindClusters",
          "Seurat",
          envir = asNamespace("Seurat")
        )
        if (
          identical(tolower(cluster_algorithm), "leiden") &&
            "leiden_method" %in% names(formals(find_clusters_method))
        ) {
          srt <- FindClusters(
            object = srt,
            resolution = cluster_resolution,
            algorithm = cluster_algorithm_index,
            graph.name = paste0(prefix, lr, "_SNN"),
            verbose = FALSE,
            leiden_method = "igraph"
          )
        } else {
          srt <- FindClusters(
            object = srt,
            resolution = cluster_resolution,
            algorithm = cluster_algorithm_index,
            graph.name = paste0(prefix, lr, "_SNN"),
            verbose = FALSE
          )
        }
        log_message("Reorder clusters...", verbose = verbose)
        srt <- srt_reorder(
          srt,
          features = HVF,
          reorder_by = "seurat_clusters",
          layer = "data"
        )
        srt[["seurat_clusters"]] <- NULL
        srt[[paste0(prefix, lr, "clusters")]] <- Idents(srt)
        srt
      },
      error = function(error) {
        err_msg <- conditionMessage(error)
        err_msg <- gsub("{", "{{", err_msg, fixed = TRUE)
        err_msg <- gsub("}", "}}", err_msg, fixed = TRUE)
        log_message(err_msg, message_type = "warning", verbose = verbose)
        log_message(
          "Error when performing {.fn Seurat::FindClusters}. Skip it",
          message_type = "warning",
          verbose = verbose
        )
        srt
      }
    )

    srt <- tryCatch(
      {
        for (nr in nonlinear_reduction) {
          log_message(
            "Perform {.pkg {nr}} nonlinear dimension reduction",
            verbose = verbose
          )
          for (n in nonlinear_reduction_dims) {
            if (
              identical(nr, "umap") &&
                isTRUE(force_nonlinear_reduction) &&
                !is.null(linear_reduction_dims_use_current) &&
                identical(as.integer(n), 2L)
            ) {
              params <- c(
                list(
                  object = srt,
                  reduction = paste0(prefix, lr),
                  dims = linear_reduction_dims_use_current,
                  n.components = n,
                  reduction.name = paste0(prefix, lr, "UMAP", n, "D"),
                  reduction.key = paste0(prefix, lr, "UMAP", n, "D_"),
                  verbose = FALSE,
                  seed.use = seed
                ),
                nonlinear_reduction_params
              )
              if (is.null(params[["cores"]])) {
                params[["cores"]] <- cores
              }
              srt <- do.call(RunUMAP, params)
              srt@misc[["Default_reduction"]] <- paste0(prefix, lr, "UMAP")
            } else {
              srt <- RunDimsReduction(
                srt,
                prefix = paste0(prefix, lr),
                reduction_use = paste0(prefix, lr),
                reduction_dims = linear_reduction_dims_use_current,
                graph_use = paste0(prefix, lr, "_SNN"),
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
      },
      error = function(error) {
        err_msg <- conditionMessage(error)
        err_msg <- gsub("{", "{{", err_msg, fixed = TRUE)
        err_msg <- gsub("}", "}}", err_msg, fixed = TRUE)
        log_message(err_msg, message_type = "warning", verbose = verbose)
        log_message(
          "Error when performing {.pkg {nr}} nonlinear dimension reduction. Skip it",
          message_type = "warning",
          verbose = verbose
        )
        srt
      }
    )
  }

  cluster_name <- paste0(prefix, linear_reduction[1], "clusters")
  if (cluster_name %in% colnames(srt@meta.data)) {
    srt[[paste0(prefix, "clusters")]] <- srt[[cluster_name]]
  }
  for (nr in nonlinear_reduction) {
    for (n in nonlinear_reduction_dims) {
      reductions_name <- paste0(
        prefix,
        linear_reduction[1],
        toupper(nr),
        n,
        "D"
      )
      if (reductions_name %in% names(srt@reductions)) {
        reduc <- srt@reductions[[reductions_name]]
        srt@reductions[[paste0(prefix, toupper(nr), n, "D")]] <- reduc
      }
    }
    srt@misc[["Default_reduction"]] <- paste0(prefix, toupper(nr))
  }

  SeuratObject::DefaultAssay(srt) <- assay
  SeuratObject::VariableFeatures(srt) <- srt@misc[["Standard_HVF"]] <- HVF

  if (is_atac_assay) {
    srt <- standardize_atac(
      srt = srt,
      prefix = prefix
    )
  }

  log_message(
    "Standard processing workflow completed",
    message_type = "success",
    text_color = "green",
    verbose = verbose
  )

  return(srt)
}

run_standard_spatial_workflow <- function(
  srt,
  prefix = "Standard",
  assay = NULL,
  image = NULL,
  coord.cols = c("x", "y"),
  do_spot_qc = TRUE,
  spot_qc_params = list(),
  do_spatial_variable_features = TRUE,
  spatial_variable_features_params = list(),
  do_spatial_cluster = FALSE,
  spatial_cluster_method = "BayesSpace",
  spatial_q = NULL,
  bayesspace_params = list(),
  reference = NULL,
  reference_label = NULL,
  reference_assay = NULL,
  do_deconvolution = !is.null(reference),
  deconvolution_method = "RCTD",
  deconvolution_params = list(),
  do_normalization = NULL,
  normalization_method = "LogNormalize",
  do_HVF_finding = TRUE,
  HVF_method = "vst",
  nHVF = 2000,
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
  nonlinear_reduction_dims = 2,
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  cores = 1L,
  verbose = TRUE,
  seed = 11,
  ...
) {
  log_message(
    "Start standard spot-level spatial workflow...",
    text_color = "blue",
    verbose = verbose
  )

  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  if (is.null(assay) || length(assay) != 1L || is.na(assay)) {
    log_message(
      "{.arg assay} must be specified for {.arg workflow = 'spatial'}",
      message_type = "error"
    )
  }
  if (!assay %in% SeuratObject::Assays(srt)) {
    log_message(
      "{.arg assay} {.val {assay}} is not present in {.cls Seurat}",
      message_type = "error"
    )
  }
  is_atac_assay <- inherits(srt[[assay]], "ChromatinAssay")
  if (is_atac_assay) {
    atac_args <- atac_defaults(
      prefix = prefix,
      do_normalization = do_normalization,
      normalization_method = normalization_method,
      do_HVF_finding = do_HVF_finding,
      nHVF = nHVF,
      do_scaling = do_scaling,
      linear_reduction = linear_reduction,
      linear_reduction_dims_use = linear_reduction_dims_use,
      neighbor_metric = neighbor_metric
    )
    prefix <- atac_args$prefix
    do_normalization <- atac_args$do_normalization
    normalization_method <- atac_args$normalization_method
    do_HVF_finding <- atac_args$do_HVF_finding
    nHVF <- atac_args$nHVF
    do_scaling <- atac_args$do_scaling
    linear_reduction <- atac_args$linear_reduction
    linear_reduction_dims_use <- atac_args$linear_reduction_dims_use
    neighbor_metric <- atac_args$neighbor_metric
  }
  cores <- validate_scalar_integer(cores, "cores")
  validate_scalar_flag(do_spot_qc, "do_spot_qc")
  validate_scalar_flag(
    do_spatial_variable_features,
    "do_spatial_variable_features"
  )
  validate_scalar_flag(do_spatial_cluster, "do_spatial_cluster")
  if (!is.null(do_deconvolution)) {
    validate_scalar_flag(do_deconvolution, "do_deconvolution")
  }

  spatial_cluster_method_label <- if (
    is.character(spatial_cluster_method) &&
      length(spatial_cluster_method) == 1L &&
      !is.na(spatial_cluster_method) &&
      nzchar(spatial_cluster_method)
  ) {
    spatial_cluster_method
  } else {
    "BayesSpace"
  }
  deconvolution_method_label <- if (
    is.character(deconvolution_method) &&
      length(deconvolution_method) == 1L &&
      !is.na(deconvolution_method) &&
      nzchar(deconvolution_method)
  ) {
    deconvolution_method
  } else {
    "RCTD"
  }
  deconvolution_requested <- isTRUE(do_deconvolution) ||
    is.null(do_deconvolution)
  stages <- data.frame(
    stage = c("quality_control", "spatial_variable_features", "spatial_clustering", "deconvolution"),
    requested = c(
      isTRUE(do_spot_qc), isTRUE(do_spatial_variable_features),
      isTRUE(do_spatial_cluster), deconvolution_requested
    ),
    status = ifelse(
      c(
        isTRUE(do_spot_qc), isTRUE(do_spatial_variable_features),
        isTRUE(do_spatial_cluster), deconvolution_requested
      ),
      "requested",
      "skipped"
    ),
    requested_method = c(
      "RunSpotQC", "RunSpatialVariableFeatures",
      spatial_cluster_method_label, deconvolution_method_label
    ),
    actual_method = NA_character_,
    result_tool_key = NA_character_,
    result_metadata_key = NA_character_,
    result_stored = FALSE,
    result_location = NA_character_,
    variable_features_before = NA_integer_,
    variable_features_after = NA_integer_,
    set_variable_features = NA,
    reason = ifelse(
      c(
        isTRUE(do_spot_qc), isTRUE(do_spatial_variable_features),
        isTRUE(do_spatial_cluster), deconvolution_requested
      ),
      NA_character_,
      "not requested"
    ),
    stringsAsFactors = FALSE
  )
  update_stage <- function(
    stage,
    status = NULL,
    actual_method = NULL,
    result_tool_key = NULL,
    result_metadata_key = NULL,
    result_stored = NULL,
    result_location = NULL,
    variable_features_before = NULL,
    variable_features_after = NULL,
    set_variable_features = NULL,
    reason = NULL
  ) {
    i <- match(stage, stages$stage)
    if (!is.null(status)) {
      stages$status[[i]] <<- status
    }
    if (!is.null(actual_method)) {
      stages$actual_method[[i]] <<- actual_method
    }
    if (!is.null(result_tool_key)) {
      stages$result_tool_key[[i]] <<- result_tool_key
    }
    if (!is.null(result_metadata_key)) {
      stages$result_metadata_key[[i]] <<- result_metadata_key
    }
    if (!is.null(result_stored)) {
      stages$result_stored[[i]] <<- isTRUE(result_stored)
    }
    if (!is.null(result_location)) {
      stages$result_location[[i]] <<- result_location
    }
    if (!is.null(variable_features_before)) {
      stages$variable_features_before[[i]] <<- variable_features_before
    }
    if (!is.null(variable_features_after)) {
      stages$variable_features_after[[i]] <<- variable_features_after
    }
    if (!is.null(set_variable_features)) {
      stages$set_variable_features[[i]] <<- set_variable_features
    }
    if (!is.null(reason)) {
      stages$reason[[i]] <<- reason
    }
    invisible(NULL)
  }
  fail_stage <- function(stage, actual_method, error) {
    i <- match(stage, stages$stage)
    if (!identical(stages$status[[i]], "failed")) {
      update_stage(
        stage = stage,
        status = "failed",
        actual_method = actual_method,
        result_tool_key = NA_character_,
        result_metadata_key = NA_character_,
        result_stored = FALSE,
        result_location = NA_character_,
        reason = conditionMessage(error)
      )
    }
    attr(error, "standard_spatial_stages") <- stages
    stop(error)
  }
  run_stage_setup <- function(stage, expr, actual_method) {
    tryCatch(
      force(expr),
      error = function(e) fail_stage(stage, actual_method, e)
    )
  }
  run_stage <- function(stage, expr, actual_method, result_probe) {
    tryCatch(
      {
        value <- force(expr)
        probe <- result_probe(value)
        if (!is.list(probe)) {
          log_message(
            "Spatial stage result probe must return a list",
            message_type = "error"
          )
        }
        update_stage(
          stage = stage,
          status = if (isTRUE(probe$result_complete)) "completed" else "failed",
          actual_method = actual_method,
          result_tool_key = probe$result_tool_key %||% NA_character_,
          result_metadata_key = probe$result_metadata_key %||% NA_character_,
          result_stored = isTRUE(probe$result_stored),
          result_location = probe$result_location %||% NA_character_,
          reason = if (isTRUE(probe$result_complete)) {
            NA_character_
          } else {
            probe$reason %||% "expected result is missing"
          }
        )
        if (!isTRUE(probe$result_complete)) {
          log_message(
            paste0(
              "Spatial stage {.val {stage}} completed without ",
              "the expected result"
            ),
            message_type = "error"
          )
        }
        value
      },
      error = function(e) fail_stage(stage, actual_method, e)
    )
  }

  spatial_cluster_method <- if (do_spatial_cluster) {
    run_stage_setup(
      stage = "spatial_clustering",
      actual_method = "RunBayesSpace",
      expr = match.arg(spatial_cluster_method, "BayesSpace")
    )
  } else {
    "BayesSpace"
  }
  spatial_clustering_row <- match("spatial_clustering", stages$stage)
  stages$requested_method[[spatial_clustering_row]] <- spatial_cluster_method

  deconvolution_method <- if (deconvolution_requested) {
    run_stage_setup(
      stage = "deconvolution",
      actual_method = "deconvolution_method validation",
      expr = match.arg(
        deconvolution_method,
        c("RCTD", "SPOTlight", "Cell2location")
      )
    )
  } else {
    "RCTD"
  }
  deconvolution_row <- match("deconvolution", stages$stage)
  stages$requested_method[[deconvolution_row]] <- deconvolution_method
  deconv_default_name <- switch(deconvolution_method,
    RCTD = "RCTD",
    SPOTlight = "SPOTlight",
    Cell2location = "Cell2location"
  )
  deconv_producer <- switch(deconvolution_method,
    RCTD = "RunRCTD",
    SPOTlight = "RunSPOTlight",
    Cell2location = "RunCell2location"
  )
  deconv_fun <- switch(deconvolution_method,
    RCTD = RunRCTD,
    SPOTlight = RunSPOTlight,
    Cell2location = RunCell2location
  )

  spot_qc_params <- if (do_spot_qc) {
    run_stage_setup(
      stage = "quality_control",
      actual_method = "RunSpotQC",
      expr = {
        validate_named_list(spot_qc_params, "spot_qc_params")
        spot_qc_params
      }
    )
  } else {
    list()
  }
  spatial_variable_features_params <- if (do_spatial_variable_features) {
    run_stage_setup(
      stage = "spatial_variable_features",
      actual_method = "RunSpatialVariableFeatures",
      expr = {
        validate_named_list(
          spatial_variable_features_params,
          "spatial_variable_features_params"
        )
        spatial_variable_features_params
      }
    )
  } else {
    list()
  }
  bayesspace_params <- if (do_spatial_cluster) {
    run_stage_setup(
      stage = "spatial_clustering",
      actual_method = "RunBayesSpace",
      expr = {
        validate_named_list(bayesspace_params, "bayesspace_params")
        bayesspace_params
      }
    )
  } else {
    list()
  }
  deconvolution_params <- if (deconvolution_requested) {
    run_stage_setup(
      stage = "deconvolution",
      actual_method = deconv_producer,
      expr = {
        validate_named_list(
          deconvolution_params,
          "deconvolution_params"
        )
        deconvolution_params
      }
    )
  } else {
    list()
  }
  if (is.null(do_deconvolution)) {
    do_deconvolution <- !is.null(reference) ||
      (
        identical(deconvolution_method, "Cell2location") &&
          !is.null(deconvolution_params[["reference_signatures"]])
      )
  }
  stages$requested[[deconvolution_row]] <- isTRUE(do_deconvolution)
  stages$status[[deconvolution_row]] <- if (isTRUE(do_deconvolution)) {
    "requested"
  } else {
    "skipped"
  }
  stages$reason[[deconvolution_row]] <- if (isTRUE(do_deconvolution)) {
    NA_character_
  } else {
    "not requested"
  }
  has_cell2location_signatures <- identical(
    deconvolution_method,
    "Cell2location"
  ) && !is.null(deconvolution_params[["reference_signatures"]])
  deconv_will_run <- isTRUE(do_deconvolution) &&
    (!is.null(reference) || has_cell2location_signatures)
  planned_deconv_prefix <- deconvolution_params[["prefix"]] %||%
    deconv_default_name
  planned_deconv_tool_name <- deconvolution_params[["tool_name"]] %||%
    deconv_default_name
  planned_deconv_store_results <-
    deconvolution_params[["store_results"]] %||% TRUE
  planned_bayesspace_cluster_colname <-
    bayesspace_params[["cluster_colname"]] %||% "BayesSpace_cluster"
  planned_bayesspace_init_colname <- if (
    "init_colname" %in% names(bayesspace_params)
  ) {
    bayesspace_params[["init_colname"]]
  } else {
    "BayesSpace_init"
  }
  cluster_col <- paste0(prefix, "clusters")
  preprocessing_cluster_cols <- standard_spatial_preprocessing_metadata_targets(
    prefix = prefix,
    linear_reduction = linear_reduction,
    normalization_method = normalization_method,
    is_atac_assay = is_atac_assay,
    cluster_resolution = cluster_resolution
  )
  spot_qc_args <- NULL
  spot_qc_metadata_targets <- character()
  if (isTRUE(do_spot_qc)) {
    spot_qc_args <- merge_call_args(
      list(
        srt = srt,
        assay = assay,
        qc_metrics = c("outlier", "umi", "gene", "mito"),
        outlier_threshold = c(
          "log10_nCount:lower:3",
          "log10_nFeature:lower:3",
          "spot_featurecount_dist:lower:3"
        ),
        verbose = verbose
      ),
      spot_qc_params
    )
    spot_qc_args$assay <- spot_qc_args$assay %||%
      SeuratObject::DefaultAssay(srt)
    spot_qc_args$qc_metrics <- spot_qc_args$qc_metrics %||%
      c("outlier", "umi", "gene", "mito")
    spot_qc_args$outlier_threshold <- spot_qc_args$outlier_threshold %||%
      c(
        "log10_nCount:lower:3",
        "log10_nFeature:lower:3",
        "spot_featurecount_dist:lower:3"
      )
    spot_qc_metadata_targets <- standard_spatial_spot_qc_metadata_targets(
      srt = srt,
      assay = spot_qc_args$assay,
      qc_metrics = spot_qc_args$qc_metrics,
      outlier_threshold = spot_qc_args$outlier_threshold
    )
  }

  metadata_output_plan <- if (isTRUE(do_spot_qc)) {
    run_stage_setup(
      stage = "quality_control",
      actual_method = "RunSpotQC",
      expr = standard_spatial_validate_metadata_output_plan(
        preprocessing_cluster_cols = preprocessing_cluster_cols,
        do_spot_qc = TRUE,
        spot_qc_metadata_targets = spot_qc_metadata_targets
      )
    )
  } else {
    standard_spatial_validate_metadata_output_plan(
      preprocessing_cluster_cols = preprocessing_cluster_cols,
      do_spot_qc = FALSE,
      spot_qc_metadata_targets = spot_qc_metadata_targets
    )
  }
  if (isTRUE(do_spatial_cluster)) {
    metadata_output_plan <- run_stage_setup(
      stage = "spatial_clustering",
      actual_method = "RunBayesSpace",
      expr = standard_spatial_add_clustering_output_plan(
        metadata_targets = metadata_output_plan$targets,
        metadata_owners = metadata_output_plan$owners,
        cluster_colname = planned_bayesspace_cluster_colname,
        init_colname = planned_bayesspace_init_colname
      )
    )
  }
  deconv_plan <- NULL
  if (deconv_will_run) {
    deconv_plan <- run_stage_setup(
      stage = "deconvolution",
      actual_method = deconv_producer,
      expr = {
        validate_scalar_string(
          planned_deconv_prefix,
          "deconvolution_params$prefix"
        )
        validate_scalar_flag(
          planned_deconv_store_results,
          "deconvolution_params$store_results"
        )
        standard_spatial_validate_deconvolution_tool_name(
          planned_deconv_tool_name,
          check_reserved = planned_deconv_store_results
        )
        standard_spatial_validate_deconvolution_output_plan(
          metadata_targets = metadata_output_plan$targets,
          metadata_owners = metadata_output_plan$owners,
          prefix = planned_deconv_prefix,
          method = deconvolution_method
        )
        if (is.null(reference_label) && !has_cell2location_signatures) {
          log_message(
            "{.arg reference_label} must be provided when running deconvolution",
            message_type = "error"
          )
        }
        params <- if (identical(deconvolution_method, "RCTD")) {
          standard_spatial_prepare_rctd_params(
            deconvolution_params = deconvolution_params,
            verbose = verbose
          )
        } else {
          deconvolution_params
        }
        defaults <- list(
          srt = srt,
          reference = reference,
          reference_label = reference_label,
          assay = assay,
          reference_assay = reference_assay,
          prefix = deconv_default_name,
          tool_name = deconv_default_name,
          store_results = TRUE,
          verbose = verbose
        )
        if (identical(deconvolution_method, "RCTD")) {
          defaults$image <- image
          defaults$coord.cols <- coord.cols
          defaults$coordinate_space <- "raw"
        }
        args <- merge_call_args(defaults, params)
        args$prefix <- planned_deconv_prefix
        args$tool_name <- planned_deconv_tool_name
        args$store_results <- planned_deconv_store_results
        list(
          args = args,
          prefix = planned_deconv_prefix,
          tool_name = planned_deconv_tool_name,
          store_results = planned_deconv_store_results
        )
      }
    )
  }

  if (isTRUE(do_spot_qc)) {
    spot_qc_args <- run_stage_setup(
      stage = "quality_control",
      actual_method = "RunSpotQC",
      expr = {
        standard_spatial_validate_spot_qc_effective_args(
          srt = srt,
          assay = spot_qc_args$assay,
          qc_metrics = spot_qc_args$qc_metrics,
          outlier_threshold = spot_qc_args$outlier_threshold
        )
        spot_qc_args
      }
    )
    spot_qc_args$srt <- standard_spatial_clear_outputs(
      srt,
      metadata_keys = "SpotQC"
    )
    srt <- run_stage(
      stage = "quality_control",
      expr = do.call(RunSpotQC, spot_qc_args),
      actual_method = "RunSpotQC",
      result_probe = function(result) {
        standard_spatial_result_probe(
          srt = result,
          metadata_keys = "SpotQC",
          metadata_required = TRUE,
          metadata_expected = 'meta.data[["SpotQC"]]'
        )
      }
    )
  }

  srt <- RunStandardWorkflow(
    srt = srt,
    prefix = prefix,
    workflow = "single_cell",
    assay = assay,
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
    cores = cores,
    verbose = verbose,
    seed = seed,
    ...
  )

  if (isTRUE(do_spatial_variable_features)) {
    svf_setup <- run_stage_setup(
      stage = "spatial_variable_features",
      actual_method = "RunSpatialVariableFeatures",
      expr = {
        svf_args <- merge_call_args(
          list(
            srt = srt,
            assay = assay,
            image = image,
            coord.cols = coord.cols,
            coordinate_space = "raw",
            set_variable_features = FALSE,
            store_results = TRUE,
            verbose = verbose,
            seed = seed
          ),
          spatial_variable_features_params
        )
        svf_assay <- svf_args$assay %||% SeuratObject::DefaultAssay(srt)
        validate_scalar_string(
          svf_assay,
          "spatial_variable_features_params$assay"
        )
        if (!svf_assay %in% SeuratObject::Assays(srt)) {
          log_message(
            "Effective spatial variable feature assay {.val {svf_assay}} is not present in {.cls Seurat}",
            message_type = "error"
          )
        }
        svf_args$assay <- svf_assay
        svf_args$set_variable_features <-
          svf_args$set_variable_features %||% TRUE
        svf_args$store_results <- svf_args$store_results %||% TRUE
        validate_scalar_flag(
          svf_args$set_variable_features,
          "spatial_variable_features_params$set_variable_features"
        )
        validate_scalar_flag(
          svf_args$store_results,
          "spatial_variable_features_params$store_results"
        )
        svf_store_results <- isTRUE(svf_args$store_results)
        svf_set_variable_features <- isTRUE(svf_args$set_variable_features)
        variable_features_before <- length(
          standard_spatial_variable_features(srt, assay = svf_assay)
        )
        if (is.null(svf_args$features)) {
          input_features <- standard_spatial_variable_features(
            srt,
            assay = svf_assay
          )
          if (length(input_features) > 0L) {
            svf_args$features <- input_features
          }
        }
        svf_input <- if (svf_set_variable_features) {
          standard_spatial_suspend_variable_features(srt, assay = svf_assay)
        } else {
          list(srt = srt, restore_metadata = NULL)
        }
        svf_args$srt <- standard_spatial_clear_outputs(
          svf_input$srt,
          tool_keys = if (svf_store_results) {
            "SpatialVariableFeatures"
          } else {
            character()
          }
        )
        list(
          args = svf_args,
          assay = svf_assay,
          store_results = svf_store_results,
          set_variable_features = svf_set_variable_features,
          variable_features_before = variable_features_before,
          restore_metadata = svf_input$restore_metadata,
          selection_marker = svf_input$selection_marker %||% NULL
        )
      }
    )
    svf_args <- svf_setup$args
    svf_assay <- svf_setup$assay
    svf_store_results <- svf_setup$store_results
    svf_set_variable_features <- svf_setup$set_variable_features
    variable_features_before <- svf_setup$variable_features_before
    svf_restore_metadata <- svf_setup$restore_metadata
    svf_selection_marker <- svf_setup$selection_marker
    srt <- run_stage(
      stage = "spatial_variable_features",
      expr = do.call(RunSpatialVariableFeatures, svf_args),
      actual_method = "RunSpatialVariableFeatures",
      result_probe = function(result) {
        variable_features <- if (svf_set_variable_features) {
          standard_spatial_variable_features(result, assay = svf_assay)
        } else {
          character()
        }
        explicit_empty_selection <- svf_set_variable_features &&
          length(variable_features) == 0L &&
          spatial_has_explicit_empty_variable_features(
            result,
            assay = svf_assay,
            expected_token = svf_selection_marker$token %||% NULL
          )
        valid_empty_selection <- explicit_empty_selection &&
          (
            !svf_store_results ||
              standard_spatial_svf_has_valid_empty_selection(result)
          )
        variable_features_stored <- length(variable_features) > 0L ||
          valid_empty_selection
        standard_spatial_result_probe(
          srt = result,
          tool_key = if (svf_store_results) {
            "SpatialVariableFeatures"
          } else {
            NULL
          },
          tool_required = svf_store_results,
          extra_locations = if (
            svf_set_variable_features && variable_features_stored
          ) {
            sprintf('VariableFeatures("%s")', svf_assay)
          } else {
            character()
          },
          extra_required = svf_set_variable_features &&
            !valid_empty_selection,
          extra_expected = sprintf('VariableFeatures("%s")', svf_assay)
        )
      }
    )
    if (svf_set_variable_features) {
      preserve_empty_selection <-
        length(standard_spatial_variable_features(srt, assay = svf_assay)) ==
          0L &&
        spatial_has_explicit_empty_variable_features(
          srt,
          assay = svf_assay,
          expected_token = svf_selection_marker$token %||% NULL
        ) &&
        (
          !svf_store_results ||
            standard_spatial_svf_has_valid_empty_selection(srt)
        )
      srt <- standard_spatial_restore_variable_feature_info(
        srt,
        assay = svf_assay,
        metadata = svf_restore_metadata,
        preserve_empty_selection = preserve_empty_selection,
        selection_marker = svf_selection_marker
      )
    }
    update_stage(
      stage = "spatial_variable_features",
      variable_features_before = variable_features_before,
      variable_features_after = length(
        standard_spatial_variable_features(srt, assay = svf_assay)
      ),
      set_variable_features = svf_set_variable_features
    )
  }

  if (isTRUE(do_spatial_cluster)) {
    bayesspace_setup <- run_stage_setup(
      stage = "spatial_clustering",
      actual_method = "RunBayesSpace",
      expr = {
        if (!identical(spatial_cluster_method, "BayesSpace")) {
          log_message(
            "{.arg spatial_cluster_method} only supports {.val BayesSpace}",
            message_type = "error"
          )
        }
        spatial_q_use <- spatial_q
        if (is.null(spatial_q_use)) {
          if (!cluster_col %in% colnames(srt@meta.data)) {
            log_message(
              "Unable to infer {.arg spatial_q}; metadata column {.val {cluster_col}} was not found",
              message_type = "error"
            )
          }
          clusters <- as.character(srt[[cluster_col, drop = TRUE]])
          spatial_q_use <- length(unique(stats::na.omit(clusters)))
          if (spatial_q_use < 2L) {
            log_message(
              "Unable to infer {.arg spatial_q}; {.val {cluster_col}} contains fewer than 2 clusters",
              message_type = "error"
            )
          }
        }
        bayesspace_args <- list(
          srt = srt,
          q = spatial_q_use,
          assay = assay,
          image = image,
          verbose = verbose
        )
        linear_reduction_use <- if (identical(
          normalization_method,
          "TFIDF"
        )) {
          "svd"
        } else {
          linear_reduction[[1L]]
        }
        reduction_use <- paste0(prefix, linear_reduction_use)
        if (
          reduction_use %in% SeuratObject::Reductions(srt) &&
            is.null(bayesspace_params[["use_reduction"]])
        ) {
          bayesspace_args[["use_reduction"]] <- reduction_use
          if (is.null(bayesspace_params[["dims"]])) {
            emb <- SeuratObject::Embeddings(srt, reduction = reduction_use)
            if (!is.null(linear_reduction_dims_use)) {
              dims_use <- linear_reduction_dims_use
            } else {
              dims_use <- seq_len(min(15L, ncol(emb)))
            }
            bayesspace_args[["dims"]] <- dims_use[dims_use <= ncol(emb)]
          }
        }
        bayesspace_args <- merge_call_args(
          bayesspace_args,
          bayesspace_params
        )
        bayesspace_args$cluster_colname <-
          planned_bayesspace_cluster_colname
        bayesspace_args["init_colname"] <- list(
          planned_bayesspace_init_colname
        )
        bayesspace_cluster_colname <- planned_bayesspace_cluster_colname
        bayesspace_init_colname <- planned_bayesspace_init_colname
        bayesspace_metadata_keys <- c(
          bayesspace_cluster_colname,
          if (!is.null(bayesspace_init_colname)) bayesspace_init_colname
        )
        bayesspace_args$srt <- standard_spatial_clear_outputs(
          srt,
          tool_keys = "BayesSpace",
          metadata_keys = bayesspace_metadata_keys
        )
        list(
          args = bayesspace_args,
          cluster_colname = bayesspace_cluster_colname,
          metadata_keys = bayesspace_metadata_keys
        )
      }
    )
    bayesspace_args <- bayesspace_setup$args
    bayesspace_cluster_colname <- bayesspace_setup$cluster_colname
    bayesspace_metadata_keys <- bayesspace_setup$metadata_keys
    srt <- run_stage(
      stage = "spatial_clustering",
      expr = do.call(RunBayesSpace, bayesspace_args),
      actual_method = "RunBayesSpace",
      result_probe = function(result) {
        standard_spatial_result_probe(
          srt = result,
          tool_key = "BayesSpace",
          tool_required = TRUE,
          metadata_keys = bayesspace_metadata_keys,
          metadata_required = TRUE,
          metadata_complete = bayesspace_cluster_colname %in%
            colnames(result@meta.data),
          metadata_expected = sprintf(
            'meta.data[["%s"]]',
            bayesspace_cluster_colname
          )
        )
      }
    )
  }

  if (isTRUE(do_deconvolution)) {
    if (is.null(reference) && !has_cell2location_signatures) {
      update_stage(
        stage = "deconvolution",
        status = "skipped",
        reason = "reference and reference_signatures are unavailable"
      )
      log_message(
        "Skip deconvolution because {.arg reference} is {.val NULL}",
        message_type = "warning",
        verbose = verbose
      )
    } else {
      deconv_runtime <- run_stage_setup(
        stage = "deconvolution",
        actual_method = deconv_producer,
        expr = {
          deconv_args <- deconv_plan$args
          deconv_args$srt <- standard_spatial_clear_outputs(
            srt,
            tool_keys = if (deconv_plan$store_results) {
              deconv_plan$tool_name
            } else {
              character()
            },
            metadata_keys = standard_spatial_deconv_metadata_keys(
              srt,
              prefix = deconv_plan$prefix,
              method = deconvolution_method
            )
          )
          list(args = deconv_args)
        }
      )
      deconv_args <- deconv_runtime$args
      deconv_prefix <- deconv_plan$prefix
      deconv_tool_name <- deconv_plan$tool_name
      deconv_store_results <- deconv_plan$store_results
      srt <- run_stage(
        stage = "deconvolution",
        actual_method = deconv_producer,
        expr = do.call(deconv_fun, deconv_args),
        result_probe = function(result) {
          metadata_keys <- standard_spatial_deconv_metadata_keys(
            result,
            prefix = deconv_prefix,
            method = deconvolution_method
          )
          standard_spatial_result_probe(
            srt = result,
            tool_key = if (deconv_store_results) {
              deconv_tool_name
            } else {
              NULL
            },
            tool_required = deconv_store_results,
            metadata_keys = metadata_keys,
            metadata_required = TRUE,
            metadata_complete = standard_spatial_deconv_metadata_complete(
              metadata_keys,
              prefix = deconv_prefix,
              method = deconvolution_method
            ),
            metadata_expected = if (identical(
              deconvolution_method,
              "Cell2location"
            )) {
              sprintf(
                paste0(
                  'meta.data columns "%s_abundance_*", "%s_prop_*", ',
                  '"%s_dominant_type", and "%s_max_prop"'
                ),
                deconv_prefix,
                deconv_prefix,
                deconv_prefix,
                deconv_prefix
              )
            } else {
              sprintf(
                paste0(
                  'meta.data columns "%s_prop_*", ',
                  '"%s_dominant_type", and "%s_max_prop"'
                ),
                deconv_prefix,
                deconv_prefix,
                deconv_prefix
              )
            }
          )
        }
      )
    }
  }

  incomplete_requested <- stages$requested & stages$status != "completed"
  workflow_status <- if (any(stages$status == "failed")) {
    "failed"
  } else if (any(incomplete_requested)) {
    "partial"
  } else {
    "completed"
  }
  srt@tools[["run_standard_spatial_workflow"]] <- list(
    status = workflow_status,
    stages = stages,
    parameters = list(
      prefix = prefix,
      assay = assay,
      image = image,
      coord.cols = coord.cols,
      do_spot_qc = do_spot_qc,
      do_spatial_variable_features = do_spatial_variable_features,
      spatial_variable_features_params = spatial_variable_features_params,
      do_spatial_cluster = do_spatial_cluster,
      spatial_cluster_method = spatial_cluster_method,
      spatial_q = spatial_q,
      do_deconvolution = do_deconvolution,
      deconvolution_method = deconvolution_method,
      cores = cores
    ),
    cluster_col = if (cluster_col %in% colnames(srt@meta.data)) {
      cluster_col
    } else {
      NULL
    }
  )

  if (identical(workflow_status, "completed")) {
    log_message(
      "Standard spot-level spatial workflow completed",
      message_type = "success",
      text_color = "green",
      verbose = verbose
    )
  } else {
    log_message(
      "Standard spot-level spatial workflow finished with skipped requested stages",
      message_type = "warning",
      verbose = verbose
    )
  }
  srt
}

standard_spatial_variable_features <- function(srt, assay) {
  variable_features <- suppressWarnings(
    SeuratObject::VariableFeatures(srt, assay = assay)
  )
  variable_features <- as.character(variable_features)
  variable_features[
    !is.na(variable_features) & nzchar(variable_features)
  ]
}

standard_spatial_suspend_variable_features <- function(srt, assay) {
  assay_object <- srt[[assay]]
  if (!inherits(assay_object, "StdAssay")) {
    SeuratObject::VariableFeatures(assay_object) <- character()
    srt[[assay]] <- assay_object
    marker <- spatial_begin_svf_selection(srt, assay = assay)
    return(list(
      srt = marker$srt,
      restore_metadata = NULL,
      selection_marker = marker$state
    ))
  }

  feature_metadata <- assay_object[[]]
  clear_columns <- unique(c(
    intersect(
      c("var.features", "var.features.rank"),
      colnames(feature_metadata)
    ),
    grep("^vf_", colnames(feature_metadata), value = TRUE)
  ))
  restore_columns <- setdiff(
    clear_columns,
    c("var.features", "var.features.rank")
  )
  restore_metadata <- feature_metadata[, restore_columns, drop = FALSE]
  for (column in clear_columns) {
    assay_object[[column]] <- NULL
  }
  srt[[assay]] <- assay_object
  list(
    srt = srt,
    restore_metadata = restore_metadata,
    selection_marker = NULL
  )
}

standard_spatial_restore_variable_feature_info <- function(
  srt,
  assay,
  metadata,
  preserve_empty_selection = FALSE,
  selection_marker = NULL
) {
  if (
    (is.null(metadata) || ncol(metadata) == 0L) &&
      !isTRUE(preserve_empty_selection) &&
      is.null(selection_marker)
  ) {
    return(srt)
  }
  assay_object <- srt[[assay]]
  if (!is.null(metadata) && ncol(metadata) > 0L) {
    for (column in colnames(metadata)) {
      values <- metadata[[column]]
      names(values) <- rownames(metadata)
      assay_object[[column]] <- values[rownames(assay_object)]
    }
  }
  srt[[assay]] <- assay_object
  if (isTRUE(preserve_empty_selection)) {
    # SeuratObject may expose an NA sentinel for the explicit StdAssay
    # zero-feature state; standard_spatial_variable_features() removes it.
    srt <- spatial_set_active_variable_features(
      srt,
      assay = assay,
      features = character()
    )
  }
  srt <- spatial_restore_svf_selection_marker(
    srt,
    assay = assay,
    state = selection_marker
  )
  srt
}

standard_spatial_svf_has_valid_empty_selection <- function(srt) {
  tool <- srt@tools[["SpatialVariableFeatures"]]
  required_result_columns <- c(
    "feature", "rank", "method", "statistic", "score", "p_value",
    "q_value", "mean", "variance", "n_spots"
  )
  required_summary_fields <- c(
    "n_features", "top_features", "top_feature_summary"
  )
  if (
    !is.list(tool) ||
      !is.data.frame(tool$result) ||
      nrow(tool$result) == 0L ||
      !all(required_result_columns %in% colnames(tool$result)) ||
      !is.character(tool$result$feature) ||
      !is.null(dim(tool$result$feature)) ||
      length(tool$result$feature) != nrow(tool$result) ||
      anyNA(tool$result$feature) ||
      any(!nzchar(tool$result$feature)) ||
      anyDuplicated(tool$result$feature) > 0L ||
      !is.numeric(tool$result$score) ||
      !is.null(dim(tool$result$score)) ||
      length(tool$result$score) != nrow(tool$result) ||
      !is.list(tool$summary) ||
      !all(required_summary_fields %in% names(tool$summary)) ||
      !is.numeric(tool$summary$n_features) ||
      length(tool$summary$n_features) != 1L ||
      is.na(tool$summary$n_features) ||
      tool$summary$n_features != nrow(tool$result) ||
      !identical(tool$summary$top_features, character()) ||
      !is.data.frame(tool$summary$top_feature_summary) ||
      nrow(tool$summary$top_feature_summary) != 0L ||
      !all(
        c("feature", "rank", "score") %in%
          colnames(tool$summary$top_feature_summary)
      ) ||
      !is.character(tool$summary$top_feature_summary$feature) ||
      !is.integer(tool$summary$top_feature_summary$rank) ||
      !is.numeric(tool$summary$top_feature_summary$score) ||
      !is.list(tool$parameters) ||
      !isTRUE(tool$parameters$set_variable_features)
  ) {
    return(FALSE)
  }

  all(!is.finite(tool$result$score))
}

standard_spatial_validate_deconvolution_tool_name <- function(
  tool_name,
  check_reserved = TRUE
) {
  validate_scalar_string(tool_name, "deconvolution_params$tool_name")
  workflow_owned_keys <- c(
    "run_standard_spatial_workflow",
    "SpotQC",
    "SpatialVariableFeatures",
    "BayesSpace"
  )
  if (isTRUE(check_reserved) && tool_name %in% workflow_owned_keys) {
    log_message(
      paste0(
        "{.arg deconvolution_params$tool_name} {.val {tool_name}} is ",
        "reserved for output owned by the spatial workflow"
      ),
      message_type = "error"
    )
  }
  invisible(TRUE)
}

standard_spatial_preprocessing_metadata_targets <- function(
  prefix,
  linear_reduction,
  normalization_method,
  is_atac_assay = FALSE,
  cluster_resolution = 0.6
) {
  prefix <- prefix %||% ""
  validate_scalar_string(prefix, "prefix", require_character = FALSE)
  valid_linear_reduction <- is.character(linear_reduction) &&
    length(linear_reduction) > 0L &&
    !anyNA(linear_reduction) &&
    all(nzchar(linear_reduction))
  if (!isTRUE(valid_linear_reduction)) {
    log_message(
      "{.arg linear_reduction} must be a non-empty character vector",
      message_type = "error"
    )
  }
  validate_scalar_string(normalization_method, "normalization_method")
  validate_scalar_flag(is_atac_assay, "is_atac_assay")
  effective_linear_reduction <- if (identical(
    normalization_method,
    "TFIDF"
  )) {
    "svd"
  } else {
    linear_reduction
  }
  resolution_cluster_cols <- unlist(
    lapply(
      effective_linear_reduction,
      function(reduction) {
        paste0(
          prefix,
          reduction,
          "_SNN_res.",
          cluster_resolution
        )
      }
    ),
    use.names = FALSE
  )
  unique(c(
    paste0(prefix, effective_linear_reduction, "clusters"),
    paste0(prefix, "clusters"),
    if (is_atac_assay) paste0(prefix, "lsiclusters"),
    resolution_cluster_cols
  ))
}

standard_spatial_spot_qc_metadata_targets <- function(
  srt,
  assay,
  qc_metrics,
  outlier_threshold
) {
  qc_metrics_available <- c("outlier", "umi", "gene", "mito")
  valid_assay <- is.character(assay) &&
    length(assay) == 1L &&
    !is.na(assay) &&
    nzchar(assay) &&
    assay %in% SeuratObject::Assays(srt)
  valid_qc_metrics <- is.character(qc_metrics) &&
    !anyNA(qc_metrics) &&
    all(nzchar(qc_metrics)) &&
    all(qc_metrics %in% qc_metrics_available)
  targets <- c(
    "percent.mito",
    "spot_featurecount_dist",
    "SpotQC",
    "spot_umi_qc",
    "spot_gene_qc",
    "spot_mito_qc",
    "spot_outlier_qc"
  )
  if (isTRUE(valid_assay)) {
    targets <- c(
      paste0("nCount_", assay),
      paste0("nFeature_", assay),
      paste0("log10_nCount_", assay),
      paste0("log10_nFeature_", assay),
      targets
    )
  }
  valid_outlier_threshold <- is.character(outlier_threshold) &&
    length(outlier_threshold) > 0L &&
    !anyNA(outlier_threshold) &&
    all(nzchar(outlier_threshold))
  if (
    isTRUE(valid_qc_metrics) &&
      "outlier" %in% qc_metrics &&
      isTRUE(valid_outlier_threshold)
  ) {
    targets <- c(
      targets,
      make.names(paste0("spot_", outlier_threshold))
    )
  }
  targets
}

standard_spatial_validate_spot_qc_effective_args <- function(
  srt,
  assay,
  qc_metrics,
  outlier_threshold
) {
  validate_scalar_string(assay, "spot_qc_params$assay")
  if (!assay %in% SeuratObject::Assays(srt)) {
    log_message(
      "Effective spot QC assay {.val {assay}} is not present in {.cls Seurat}",
      message_type = "error"
    )
  }
  qc_metrics_available <- c("outlier", "umi", "gene", "mito")
  valid_qc_metrics <- is.character(qc_metrics) &&
    !anyNA(qc_metrics) &&
    all(nzchar(qc_metrics)) &&
    all(qc_metrics %in% qc_metrics_available)
  if (!isTRUE(valid_qc_metrics)) {
    log_message(
      "{.arg spot_qc_params$qc_metrics} must contain only {.val {qc_metrics_available}}",
      message_type = "error"
    )
  }
  if ("outlier" %in% qc_metrics) {
    valid_outlier_threshold <- is.character(outlier_threshold) &&
      length(outlier_threshold) > 0L &&
      !anyNA(outlier_threshold) &&
      all(nzchar(outlier_threshold))
    if (!isTRUE(valid_outlier_threshold)) {
      log_message(
        paste0(
          "{.arg spot_qc_params$outlier_threshold} must be a non-empty ",
          "character vector of {.val metric:direction:nmads} rules"
        ),
        message_type = "error"
      )
    }
  }
  invisible(TRUE)
}

standard_spatial_validate_metadata_output_plan <- function(
  preprocessing_cluster_cols,
  do_spot_qc,
  spot_qc_metadata_targets
) {
  valid_preprocessing_targets <- is.character(preprocessing_cluster_cols) &&
    length(preprocessing_cluster_cols) > 0L &&
    !anyNA(preprocessing_cluster_cols) &&
    all(nzchar(preprocessing_cluster_cols))
  if (!isTRUE(valid_preprocessing_targets)) {
    log_message(
      "Planned preprocessing cluster outputs must be non-empty strings",
      message_type = "error"
    )
  }
  preprocessing_cluster_cols <- unique(preprocessing_cluster_cols)
  metadata_targets <- preprocessing_cluster_cols
  metadata_owners <- rep(
    "single_cell_preprocessing cluster_col",
    length(preprocessing_cluster_cols)
  )
  if (isTRUE(do_spot_qc)) {
    metadata_targets <- c(metadata_targets, spot_qc_metadata_targets)
    metadata_owners <- c(
      metadata_owners,
      rep("quality_control", length(spot_qc_metadata_targets))
    )
  }
  duplicated_targets <- unique(metadata_targets[
    duplicated(metadata_targets) |
      duplicated(metadata_targets, fromLast = TRUE)
  ])
  if (length(duplicated_targets) > 0L) {
    log_message(
      paste0(
        "Planned spatial workflow metadata outputs collide at {.val ",
        "{duplicated_targets}}; choose distinct output names"
      ),
      message_type = "error"
    )
  }

  list(targets = metadata_targets, owners = metadata_owners)
}

standard_spatial_add_clustering_output_plan <- function(
  metadata_targets,
  metadata_owners,
  cluster_colname,
  init_colname
) {
  validate_scalar_string(
    cluster_colname,
    "bayesspace_params$cluster_colname"
  )
  metadata_targets <- c(metadata_targets, as.character(cluster_colname))
  metadata_owners <- c(metadata_owners, "spatial_clustering cluster_colname")
  if (!is.null(init_colname)) {
    validate_scalar_string(
      init_colname,
      "bayesspace_params$init_colname"
    )
    metadata_targets <- c(metadata_targets, as.character(init_colname))
    metadata_owners <- c(metadata_owners, "spatial_clustering init_colname")
  }

  duplicated_targets <- unique(metadata_targets[
    duplicated(metadata_targets) |
      duplicated(metadata_targets, fromLast = TRUE)
  ])
  if (length(duplicated_targets) > 0L) {
    log_message(
      paste0(
        "Planned spatial workflow metadata outputs collide at {.val ",
        "{duplicated_targets}}; choose distinct output names"
      ),
      message_type = "error"
    )
  }

  list(targets = metadata_targets, owners = metadata_owners)
}

standard_spatial_validate_deconvolution_output_plan <- function(
  metadata_targets,
  metadata_owners,
  prefix,
  method
) {
  deconvolution_collisions <- metadata_targets[
    standard_spatial_is_deconv_metadata_key(
      metadata_targets,
      prefix = as.character(prefix),
      method = method
    )
  ]
  if (length(deconvolution_collisions) > 0L) {
    collision_owners <- metadata_owners[
      metadata_targets %in% deconvolution_collisions
    ]
    log_message(
      paste0(
        "Planned deconvolution metadata with prefix {.val ",
        "{prefix}} collides with {.val ",
        "{deconvolution_collisions}} from spatial workflow stage(s) ",
        "{.val {collision_owners}}; choose a distinct prefix or output name"
      ),
      message_type = "error"
    )
  }
  invisible(TRUE)
}

standard_spatial_result_probe <- function(
  srt,
  tool_key = NULL,
  tool_required = FALSE,
  metadata_keys = character(),
  metadata_required = FALSE,
  metadata_complete = NULL,
  metadata_expected = NULL,
  extra_locations = character(),
  extra_required = FALSE,
  extra_expected = NULL
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "Spatial stage result probe requires a {.cls Seurat} object",
      message_type = "error"
    )
  }

  tool_key <- if (
    is.character(tool_key) &&
      length(tool_key) == 1L &&
      !is.na(tool_key) &&
      nzchar(tool_key)
  ) {
    tool_key
  } else {
    NULL
  }
  tool_stored <- !is.null(tool_key) && !is.null(srt@tools[[tool_key]])

  metadata_keys <- unique(as.character(metadata_keys))
  metadata_keys <- metadata_keys[
    !is.na(metadata_keys) &
      nzchar(metadata_keys) &
      metadata_keys %in% colnames(srt@meta.data)
  ]
  metadata_stored <- length(metadata_keys) > 0L
  metadata_complete <- if (is.null(metadata_complete)) {
    metadata_stored
  } else {
    isTRUE(metadata_complete)
  }

  extra_locations <- unique(as.character(extra_locations))
  extra_locations <- extra_locations[
    !is.na(extra_locations) & nzchar(extra_locations)
  ]
  extra_stored <- length(extra_locations) > 0L

  locations <- character()
  if (tool_stored) {
    locations <- c(locations, sprintf('tools[["%s"]]', tool_key))
  }
  if (metadata_stored) {
    locations <- c(
      locations,
      paste0("meta.data: ", paste(metadata_keys, collapse = ", "))
    )
  }
  if (extra_stored) {
    locations <- c(locations, extra_locations)
  }

  missing <- character()
  if (isTRUE(tool_required) && !tool_stored) {
    tool_expected <- if (is.null(tool_key)) {
      "required tools output"
    } else {
      sprintf('tools[["%s"]]', tool_key)
    }
    missing <- c(
      missing,
      tool_expected
    )
  }
  if (isTRUE(metadata_required) && !metadata_complete) {
    missing <- c(
      missing,
      metadata_expected %||% "required metadata output"
    )
  }
  if (isTRUE(extra_required) && !extra_stored) {
    missing <- c(
      missing,
      extra_expected %||% "required additional output"
    )
  }

  list(
    result_stored = length(locations) > 0L,
    result_complete = length(missing) == 0L,
    result_location = if (length(locations) > 0L) {
      paste(locations, collapse = "; ")
    } else {
      NA_character_
    },
    result_tool_key = if (tool_stored) tool_key else NA_character_,
    result_metadata_key = if (metadata_stored) {
      paste(metadata_keys, collapse = "; ")
    } else {
      NA_character_
    },
    reason = if (length(missing) > 0L) {
      paste0("missing expected output: ", paste(missing, collapse = ", "))
    } else {
      NA_character_
    }
  )
}

standard_spatial_clear_outputs <- function(
  srt,
  tool_keys = character(),
  metadata_keys = character()
) {
  tool_keys <- unique(as.character(tool_keys))
  tool_keys <- tool_keys[!is.na(tool_keys) & nzchar(tool_keys)]
  for (tool_key in tool_keys) {
    srt@tools[[tool_key]] <- NULL
  }

  metadata_keys <- unique(as.character(metadata_keys))
  metadata_keys <- metadata_keys[
    !is.na(metadata_keys) &
      nzchar(metadata_keys) &
      metadata_keys %in% colnames(srt@meta.data)
  ]
  for (metadata_key in metadata_keys) {
    srt[[metadata_key]] <- NULL
  }
  srt
}

standard_spatial_deconv_metadata_keys <- function(srt, prefix, method) {
  metadata_keys <- colnames(srt@meta.data)
  metadata_keys[
    standard_spatial_is_deconv_metadata_key(metadata_keys, prefix, method)
  ]
}

standard_spatial_is_deconv_metadata_key <- function(
  metadata_keys,
  prefix,
  method
) {
  summary_key <- metadata_keys %in% paste0(
    prefix,
    c("_dominant_type", "_max_prop")
  )
  proportion_key <- startsWith(metadata_keys, paste0(prefix, "_prop_"))
  abundance_key <- if (identical(method, "Cell2location")) {
    startsWith(metadata_keys, paste0(prefix, "_abundance_"))
  } else {
    rep(FALSE, length(metadata_keys))
  }
  proportion_key | summary_key | abundance_key
}

standard_spatial_deconv_metadata_complete <- function(
  metadata_keys,
  prefix,
  method
) {
  core_complete <- any(
    startsWith(metadata_keys, paste0(prefix, "_prop_"))
  ) && all(
    paste0(prefix, c("_dominant_type", "_max_prop")) %in% metadata_keys
  )
  if (identical(method, "Cell2location")) {
    return(
      core_complete &&
        any(startsWith(metadata_keys, paste0(prefix, "_abundance_")))
    )
  }
  core_complete
}

merge_call_args <- function(defaults, extra) {
  if (length(extra) == 0L) {
    return(defaults)
  }
  for (nm in names(extra)) {
    defaults[[nm]] <- extra[[nm]]
  }
  defaults
}

standard_spatial_prepare_rctd_params <- function(
  deconvolution_params,
  verbose = TRUE
) {
  min_cells <- deconvolution_params[["min_cells"]] %||% 25
  if (
    is.numeric(min_cells) &&
      length(min_cells) == 1L &&
      !is.na(min_cells) &&
      min_cells < 25
  ) {
    log_message(
      "{.pkg spacexr} RCTD requires at least 25 reference cells per cell type; set {.arg min_cells} from {.val {min_cells}} to {.val 25}",
      message_type = "warning",
      verbose = verbose
    )
    deconvolution_params[["min_cells"]] <- 25
  }
  deconvolution_params
}

resolve_standard_assays <- function(srt, assay = NULL) {
  assays_available <- SeuratObject::Assays(srt)
  if (is.null(assay)) {
    assay_default <- SeuratObject::DefaultAssay(srt)
    chrom_assays <- assays_available[vapply(
      assays_available,
      function(x) inherits(srt[[x]], "ChromatinAssay"),
      logical(1)
    )]
    return(unique(c(assay_default, setdiff(chrom_assays, assay_default))))
  }
  assay <- unique(assay)
  if (any(!assay %in% assays_available)) {
    log_message(
      "{.arg assay} must be present in {.cls Seurat}: {.val {setdiff(assay, assays_available)}}",
      message_type = "error"
    )
  }
  assay
}

resolve_assay_prefix <- function(srt, assay) {
  if (inherits(srt[[assay]], "ChromatinAssay")) {
    if (identical(tolower(assay), "peaks")) {
      return("ATAC")
    }
    return(assay)
  }
  assay
}

resolve_standard_prefix <- function(
  srt,
  assay,
  prefix = "Standard",
  multi_assay = FALSE,
  primary_assay = NULL
) {
  if (!multi_assay) {
    return(prefix)
  }
  if (identical(prefix, "Standard")) {
    return(resolve_assay_prefix(srt = srt, assay = assay))
  }
  if (!is.null(primary_assay) && identical(assay, primary_assay)) {
    return(prefix)
  }
  resolve_assay_prefix(srt = srt, assay = assay)
}

atac_defaults <- function(
  prefix,
  do_normalization,
  normalization_method,
  do_HVF_finding,
  nHVF,
  do_scaling,
  linear_reduction,
  linear_reduction_dims_use,
  neighbor_metric
) {
  list(
    prefix = if (identical(prefix, "Standard")) "ATAC" else prefix,
    do_normalization = if (is.null(do_normalization)) {
      TRUE
    } else {
      do_normalization
    },
    normalization_method = if (
      identical(normalization_method, "LogNormalize")
    ) {
      "TFIDF"
    } else {
      normalization_method
    },
    do_HVF_finding = if (is.null(do_HVF_finding)) TRUE else do_HVF_finding,
    nHVF = if (identical(nHVF, 2000)) 20000 else nHVF,
    do_scaling = if (isTRUE(do_scaling)) FALSE else do_scaling,
    linear_reduction = if (identical(linear_reduction, "pca")) {
      "svd"
    } else {
      linear_reduction
    },
    linear_reduction_dims_use = if (is.null(linear_reduction_dims_use)) {
      2:30
    } else {
      linear_reduction_dims_use
    },
    neighbor_metric = if (identical(neighbor_metric, "euclidean")) {
      "cosine"
    } else {
      neighbor_metric
    }
  )
}

standardize_atac <- function(srt, prefix = "ATAC") {
  if (!inherits(srt, "Seurat")) {
    return(srt)
  }
  assay <- SeuratObject::DefaultAssay(srt)
  if (!inherits(srt[[assay]], "ChromatinAssay")) {
    return(srt)
  }

  prefix <- prefix %||% ""
  svd_name <- paste0(prefix, "svd")
  lsi_name <- paste0(prefix, "lsi")
  if (svd_name %in% names(srt@reductions)) {
    reduc <- srt@reductions[[svd_name]]
    SeuratObject::Key(reduc) <- paste0(lsi_name, "_")
    srt@reductions[[lsi_name]] <- reduc
  }

  svd_cluster <- paste0(prefix, "svdclusters")
  lsi_cluster <- paste0(prefix, "lsiclusters")
  if (svd_cluster %in% colnames(srt@meta.data)) {
    srt[[lsi_cluster]] <- srt[[svd_cluster]]
  } else {
    prefix_cluster <- paste0(prefix, "clusters")
    if (prefix_cluster %in% colnames(srt@meta.data)) {
      srt[[lsi_cluster]] <- srt[[prefix_cluster]]
    }
  }

  default_reduction <- srt@misc[["Default_reduction"]] %||% NULL
  if (is.character(default_reduction) && length(default_reduction) == 1) {
    default_reduction <- sub(
      paste0("^", prefix, "svd"),
      paste0(prefix, "lsi"),
      default_reduction
    )
    srt@misc[["Default_reduction"]] <- default_reduction
  }

  srt@misc[["ATAC_default_linear_reduction"]] <- lsi_name
  srt@misc[["ATAC_default_cluster_col"]] <- if (
    lsi_cluster %in% colnames(srt@meta.data)
  ) {
    lsi_cluster
  } else {
    NULL
  }
  srt
}
