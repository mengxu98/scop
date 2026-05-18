#' @title Standard workflow for scop
#'
#' @description
#' This function performs a standard single-cell analysis workflow.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt A Seurat object.
#' @param prefix A prefix to add to the names of intermediate objects created by the function.
#' Default is `"Standard"`.
#' @param assay Which assay to use.
#' If `NULL`, the default assay of the Seurat object will be used.
#' When the object also contains `ChromatinAssay`, the default assay and
#' additional `ChromatinAssay` will be preprocessed sequentially.
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
standard_scop <- function(
  srt,
  prefix = "Standard",
  assay = NULL,
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
  if (cluster_algorithm == "leiden") {
    PrepareEnv(modules = "scanpy")
    check_python("leidenalg")
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
