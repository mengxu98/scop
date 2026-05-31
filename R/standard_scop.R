#' @title Standard workflow for scop
#'
#' @description
#' This function performs a standard single-cell or spot-level spatial analysis
#' workflow.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt A Seurat object.
#' @param prefix A prefix to add to the names of intermediate objects created by the function.
#' Default is `"Standard"`.
#' @param workflow Workflow to run. `"single_cell"` keeps the original standard
#' workflow. `"spatial"` runs a spot-level spatial workflow that wraps the
#' original workflow with spot QC, spatial variable features, optional spatial
#' clustering, and optional RCTD deconvolution.
#' @param assay Which assay to use.
#' If `NULL`, the default assay of the Seurat object will be used.
#' When the object also contains `ChromatinAssay`, the default assay and
#' additional `ChromatinAssay` will be preprocessed sequentially.
#' @param image Name of the Seurat spatial image used by the spatial workflow.
#' If `NULL`, the first image is used when present.
#' @param coord.cols Metadata coordinate columns used by the spatial workflow
#' when no image is available.
#' @param do_spot_qc Whether to run [RunSpotQC()] in the spatial workflow.
#' @param spot_qc_params Named list of additional arguments passed to
#' [RunSpotQC()].
#' @param do_spatial_variable_features Whether to run
#' [RunSpatialVariableFeatures()] in the spatial workflow.
#' @param spatial_variable_features_params Named list of additional arguments
#' passed to [RunSpatialVariableFeatures()].
#' @param do_spatial_cluster Whether to run spatial-aware clustering in the
#' spatial workflow.
#' @param spatial_cluster_method Spatial clustering method. Only
#' `"BayesSpace"` is supported in this workflow.
#' @param spatial_q Number of spatial clusters for [RunBayesSpace()]. If
#' `NULL`, the number of ordinary spot clusters is used.
#' @param bayesspace_params Named list of additional arguments passed to
#' [RunBayesSpace()].
#' @param reference Optional single-cell reference used for spatial
#' deconvolution.
#' @param reference_label Metadata column in `reference` containing cell type
#' labels.
#' @param reference_assay Assay used in `reference` for deconvolution.
#' @param do_deconvolution Whether to run deconvolution in the spatial workflow.
#' If `NULL`, deconvolution is run only when `reference` is provided.
#' @param deconvolution_method Deconvolution method. Only `"RCTD"` is supported
#' in this workflow.
#' @param deconvolution_params Named list of additional arguments passed to
#' [RunRCTD()].
#' @param do_normalization Whether to perform normalization.
#' If `NULL`, normalization will be performed if the specified assay does not have scaled data.
#' @param normalization_method The method to use for normalization.
#' Options are `"LogNormalize"`, `"SCT"`, `"TFIDF"`, or `"scran"`.
#' When `"SCT"` is used on an RNA assay, downstream reductions and clustering
#' are run on the generated `"SCT"` assay.
#' Default is `"LogNormalize"`.
#' @param do_HVF_finding Whether to perform high variable feature finding.
#' If `TRUE`, the function will force to find the highly variable features (HVF) using the specified HVF method.
#' @param HVF_method The method to use for finding highly variable features.
#' Options are `"vst"`, `"mvp"`, `"disp"`, or `"scran"`.
#' Default is `"vst"`.
#' @param nHVF The number of highly variable features to select.
#' If NULL, all highly variable features will be used.
#' Default is `2000`.
#' @param HVF A vector of feature names to use as highly variable features.
#' If NULL, the function will use the highly variable features identified by the HVF method.
#' @param do_scaling Whether to perform scaling.
#' If `TRUE`, the function will force to scale the data using the [Seurat::ScaleData] function.
#' @param vars_to_regress A vector of feature names to use as regressors in the scaling step.
#' If NULL, no regressors will be used.
#' @param regression_model The regression model to use for scaling.
#' Options are `"linear"`, `"poisson"`, or `"negativebinomial"`.
#' Default is `"linear"`.
#' @param linear_reduction The linear dimensionality reduction method to use.
#' Options are `"pca"`, `"svd"`, `"ica"`, `"nmf"`, `"mds"`, or `"glmpca"`.
#' Default is `"pca"`.
#' @param linear_reduction_dims
#' The number of dimensions to keep after linear dimensionality reduction.
#' Default is `50`.
#' @param linear_reduction_dims_use The dimensions to use for downstream analysis.
#' If `NULL`, estimated dimensions stored in the linear reduction will be used when available;
#' otherwise, the first up to `50` dimensions will be used as a fallback.
#' @param linear_reduction_params A list of parameters to pass to the linear dimensionality reduction method.
#' @param force_linear_reduction Whether to force linear dimensionality reduction even if the specified reduction is already present in the Seurat object.
#' @param nonlinear_reduction The nonlinear dimensionality reduction method to use.
#' Options are `"umap"`, `"umap-naive"`, `"tsne"`, `"dm"`,
#' `"phate"`, `"pacmap"`, `"trimap"`, `"largevis"`, or `"fr"`.
#' Default is `"umap"`.
#' @param nonlinear_reduction_dims The number of dimensions to keep after nonlinear dimensionality reduction.
#' If a vector is provided, different numbers of dimensions can be specified for each method.
#' Default is `c(2, 3)`.
#' @param nonlinear_reduction_params A list of parameters to pass to the nonlinear dimensionality reduction method.
#' @param force_nonlinear_reduction Whether to force nonlinear dimensionality reduction even if the specified reduction is already present in the Seurat object.
#' Default is `TRUE`.
#' @param neighbor_metric The distance metric to use for finding neighbors.
#' Options are `"euclidean"`, `"cosine"`, `"manhattan"`, or `"hamming"`.
#' Default is `"euclidean"`.
#' @param neighbor_k The number of nearest neighbors to use for finding neighbors.
#' Default is `20`.
#' @param cluster_algorithm The clustering algorithm to use.
#' Options are `"louvain"`, `"slm"`, or `"leiden"`.
#' Default is `"louvain"`.
#' @param cluster_resolution The resolution parameter to use for clustering.
#' Larger values result in fewer clusters.
#' Default is `0.6`.
#' @param seed Random seed for reproducibility.
#' Default is `11`.
#' @param ... Additional parameters to pass to the dimensionality reduction methods.
#'
#' @return A `Seurat` object.
#'
#' @export
#'
#' @examples
#' library(Matrix)
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
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
#' pancreas_sub <- standard_scop(
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
#' pancreas_sub <- standard_scop(
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
#' if (requireNamespace("scran", quietly = TRUE)) {
#'   data(pancreas_sub)
#'   pancreas_scran <- pancreas_sub[, 1:80]
#'   pancreas_scran <- standard_scop(
#'     pancreas_scran,
#'     assay = "RNA",
#'     do_normalization = TRUE,
#'     normalization_method = "scran",
#'     do_HVF_finding = TRUE,
#'     HVF_method = "scran",
#'     nHVF = 100,
#'     linear_reduction_dims = 10,
#'     linear_reduction_dims_use = 1:5,
#'     nonlinear_reduction = "umap",
#'     nonlinear_reduction_dims = 2
#'   )
#'   CellDimPlot(
#'     pancreas_scran,
#'     reduction = "StandardUMAP2D",
#'     group.by = "Standardclusters"
#'   )
#' }
#'
#' \dontrun{
#' data(visium_human_pancreas_sub)
#' spatial <- standard_scop(
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
#' SpatialSpotPlot(spatial, group.by = "Standardclusters")
#' SpatialSpotPlot(
#'   spatial,
#'   features = spatial@misc[["SpatialVariableFeatures"]][1:2]
#' )
#'
#' spatial_bayes <- standard_scop(
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
#'
#' data(panc8_sub)
#' spatial_rctd <- standard_scop(
#'   visium_human_pancreas_sub,
#'   workflow = "spatial",
#'   assay = "Spatial",
#'   do_spatial_cluster = FALSE,
#'   spatial_cluster_method = "BayesSpace",
#'   do_deconvolution = TRUE,
#'   deconvolution_method = "RCTD",
#'   reference = panc8_sub,
#'   reference_assay = "RNA",
#'   reference_label = "celltype",
#'   deconvolution_params = list(
#'     max_cores = 1,
#'     min_cells = 25
#'   )
#' )
#' SpatialSpotPlot(spatial_rctd, group.by = "RCTD_dominant_type")
#'
#' rctd_cols <- grep(
#'   "^RCTD_prop_",
#'   colnames(spatial_rctd@meta.data),
#'   value = TRUE
#' )
#' SpatialSpotPlot(
#'   spatial_rctd,
#'   group.by = rctd_cols[1:min(4, length(rctd_cols))]
#' )
#' SpatialSpotPlot(
#'   spatial_rctd,
#'   group.by = "RCTD_dominant_type",
#'   plot_type = "pie"
#' )
#' }
standard_scop <- function(
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
  nonlinear_reduction_dims = c(2, 3),
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  verbose = TRUE,
  seed = 11,
  ...
) {
  workflow <- match.arg(workflow)
  if (identical(workflow, "spatial")) {
    return(standard_spatial_scop(
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

  assays_to_run <- standard_scop_resolve_assays(
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
      prefix_i <- standard_scop_resolve_prefix(
        srt = srt,
        assay = assay_i,
        prefix = prefix,
        multi_assay = TRUE,
        primary_assay = assays_to_run[[1]]
      )
      srt <- standard_scop(
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
        "Perform {.fn Seurat::ScaleData}",
        verbose = verbose
      )
      srt <- suppressWarnings(
        Seurat::ScaleData(
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
        srt <- Seurat::FindNeighbors(
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
        srt <- Seurat::FindClusters(
          object = srt,
          resolution = cluster_resolution,
          algorithm = cluster_algorithm_index,
          leiden_method = "igraph",
          graph.name = paste0(prefix, lr, "_SNN"),
          verbose = FALSE
        )
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
        srt
      },
      error = function(error) {
        err_msg <- conditionMessage(error)
        err_msg <- gsub("{", "{{", err_msg, fixed = TRUE)
        err_msg <- gsub("}", "}}", err_msg, fixed = TRUE)
        log_message(err_msg, message_type = "warning", verbose = verbose)
        log_message(
          "Error when performing {.pkg {nr}} nonlinear dimension reduction. Skip it",
          message_type = "error"
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

standard_spatial_scop <- function(
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
  nonlinear_reduction_dims = c(2, 3),
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
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

  standard_scop_validate_named_list(spot_qc_params, "spot_qc_params")
  standard_scop_validate_named_list(
    spatial_variable_features_params,
    "spatial_variable_features_params"
  )
  standard_scop_validate_named_list(bayesspace_params, "bayesspace_params")
  standard_scop_validate_named_list(
    deconvolution_params,
    "deconvolution_params"
  )

  spatial_cluster_method <- match.arg(spatial_cluster_method, "BayesSpace")
  deconvolution_method <- match.arg(deconvolution_method, "RCTD")
  if (is.null(do_deconvolution)) {
    do_deconvolution <- !is.null(reference)
  }

  if (isTRUE(do_spot_qc)) {
    spot_qc_args <- standard_scop_merge_args(
      list(
        srt = srt,
        assay = assay,
        verbose = verbose
      ),
      spot_qc_params
    )
    srt <- do.call(RunSpotQC, spot_qc_args)
  }

  srt <- standard_scop(
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
    verbose = verbose,
    seed = seed,
    ...
  )

  cluster_col <- paste0(prefix, "clusters")
  if (isTRUE(do_spatial_variable_features)) {
    svf_args <- standard_scop_merge_args(
      list(
        srt = srt,
        assay = assay,
        image = image,
        coord.cols = coord.cols,
        verbose = verbose,
        seed = seed
      ),
      spatial_variable_features_params
    )
    srt <- do.call(RunSpatialVariableFeatures, svf_args)
  }

  if (isTRUE(do_spatial_cluster)) {
    if (!identical(spatial_cluster_method, "BayesSpace")) {
      log_message(
        "{.arg spatial_cluster_method} only supports {.val BayesSpace}",
        message_type = "error"
      )
    }
    spatial_q_use <- spatial_q %||% standard_spatial_infer_q(
      srt = srt,
      cluster_col = cluster_col
    )
    bayesspace_args <- list(
      srt = srt,
      q = spatial_q_use,
      assay = assay,
      image = image,
      verbose = verbose
    )
    linear_reduction_use <- linear_reduction[[1L]]
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
    bayesspace_args <- standard_scop_merge_args(
      bayesspace_args,
      bayesspace_params
    )
    srt <- do.call(RunBayesSpace, bayesspace_args)
  }

  if (isTRUE(do_deconvolution)) {
    if (is.null(reference)) {
      log_message(
        "Skip deconvolution because {.arg reference} is {.val NULL}",
        message_type = "warning",
        verbose = verbose
      )
    } else {
      if (is.null(reference_label)) {
        log_message(
          "{.arg reference_label} must be provided when running deconvolution",
          message_type = "error"
        )
      }
      if (!identical(deconvolution_method, "RCTD")) {
        log_message(
          "{.arg deconvolution_method} only supports {.val RCTD}",
          message_type = "error"
        )
      }
      deconvolution_params <- standard_spatial_prepare_rctd_params(
        deconvolution_params = deconvolution_params,
        verbose = verbose
      )
      rctd_args <- standard_scop_merge_args(
        list(
          srt = srt,
          reference = reference,
          reference_label = reference_label,
          assay = assay,
          reference_assay = reference_assay,
          image = image,
          coord.cols = coord.cols,
          verbose = verbose
        ),
        deconvolution_params
      )
      srt <- do.call(RunRCTD, rctd_args)
    }
  }

  srt@tools[["standard_spatial_scop"]] <- list(
    parameters = list(
      prefix = prefix,
      assay = assay,
      image = image,
      coord.cols = coord.cols,
      do_spot_qc = do_spot_qc,
      do_spatial_variable_features = do_spatial_variable_features,
      do_spatial_cluster = do_spatial_cluster,
      spatial_cluster_method = spatial_cluster_method,
      spatial_q = spatial_q,
      do_deconvolution = do_deconvolution,
      deconvolution_method = deconvolution_method
    ),
    cluster_col = if (cluster_col %in% colnames(srt@meta.data)) {
      cluster_col
    } else {
      NULL
    }
  )

  log_message(
    "Standard spot-level spatial workflow completed",
    message_type = "success",
    text_color = "green",
    verbose = verbose
  )
  srt
}

standard_scop_validate_named_list <- function(x, arg_name) {
  if (!is.list(x)) {
    log_message(
      "{.arg {arg_name}} must be a named list",
      message_type = "error"
    )
  }
  if (
    length(x) > 0L &&
      (is.null(names(x)) || any(is.na(names(x)) | !nzchar(names(x))))
  ) {
    log_message(
      "{.arg {arg_name}} must be a named list",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

standard_scop_merge_args <- function(defaults, extra) {
  if (length(extra) == 0L) {
    return(defaults)
  }
  for (nm in names(extra)) {
    defaults[[nm]] <- extra[[nm]]
  }
  defaults
}

standard_spatial_infer_q <- function(srt, cluster_col) {
  if (!cluster_col %in% colnames(srt@meta.data)) {
    log_message(
      "Unable to infer {.arg spatial_q}; metadata column {.val {cluster_col}} was not found",
      message_type = "error"
    )
  }
  clusters <- as.character(srt[[cluster_col, drop = TRUE]])
  q <- length(unique(stats::na.omit(clusters)))
  if (q < 2L) {
    log_message(
      "Unable to infer {.arg spatial_q}; {.val {cluster_col}} contains fewer than 2 clusters",
      message_type = "error"
    )
  }
  q
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

standard_scop_resolve_assays <- function(srt, assay = NULL) {
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

standard_scop_assay_prefix <- function(srt, assay) {
  if (inherits(srt[[assay]], "ChromatinAssay")) {
    if (identical(tolower(assay), "peaks")) {
      return("ATAC")
    }
    return(assay)
  }
  assay
}

standard_scop_resolve_prefix <- function(
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
    return(standard_scop_assay_prefix(srt = srt, assay = assay))
  }
  if (!is.null(primary_assay) && identical(assay, primary_assay)) {
    return(prefix)
  }
  standard_scop_assay_prefix(srt = srt, assay = assay)
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
  if (
    svd_name %in% names(srt@reductions) && !lsi_name %in% names(srt@reductions)
  ) {
    reduc <- srt@reductions[[svd_name]]
    SeuratObject::Key(reduc) <- paste0(lsi_name, "_")
    srt@reductions[[lsi_name]] <- reduc
  }

  svd_cluster <- paste0(prefix, "svdclusters")
  lsi_cluster <- paste0(prefix, "lsiclusters")
  if (
    svd_cluster %in%
      colnames(srt@meta.data) &&
      !lsi_cluster %in% colnames(srt@meta.data)
  ) {
    srt[[lsi_cluster]] <- srt[[svd_cluster]]
  }
  prefix_cluster <- paste0(prefix, "clusters")
  if (
    prefix_cluster %in%
      colnames(srt@meta.data) &&
      !lsi_cluster %in% colnames(srt@meta.data)
  ) {
    srt[[lsi_cluster]] <- srt[[prefix_cluster]]
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
