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
#' @param assay The name of the assay to use for the analysis.
#' If `NULL`, the default assay of the Seurat object will be used.
#' @param do_normalization Whether to perform normalization.
#' If `NULL`, normalization will be performed if the specified assay does not have scaled data.
#' @param normalization_method The method to use for normalization.
#' Options are `"LogNormalize"`, `"SCT"`, or `"TFIDF"`.
#' Default is `"LogNormalize"`.
#' @param do_HVF_finding Whether to perform high variable feature finding.
#' If `TRUE`, the function will force to find the highly variable features (HVF) using the specified HVF method.
#' @param HVF_method The method to use for finding highly variable features.
#' Options are `"vst"`, `"mvp"`, or `"disp"`.
#' Default is `"vst"`.
#' @param nHVF The number of highly variable features to select.
#' If NULL, all highly variable features will be used.
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
#' If `NULL`, all dimensions will be used.
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
#' @param seed The random seed to use for reproducibility.
#' Default is `11`.
#'
#' @return A `Seurat` object.
#'
#' @seealso [integration_scop]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(Matrix)
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType"
#' )
#'
#' # Use a combination of different linear
#' # or non-linear dimension reduction methods
#' linear_reductions <- c(
#'   "pca", "nmf", "mds", "glmpca"
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
#'       xlab = "", ylab = "", title = lr,
#'       legend.position = "none",
#'       theme_use = "theme_blank"
#'     )
#'   }
#' )
#' patchwork::wrap_plots(plotlist = plist1)
#'
#' nonlinear_reductions <- c(
#'   "umap", "tsne", "dm", "phate",
#'   "pacmap", "trimap", "largevis", "fr"
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
#'         "Standardpca", toupper(nr), "2D"
#'       ),
#'       xlab = "", ylab = "", title = nr,
#'       legend.position = "none",
#'       theme_use = "theme_blank"
#'     )
#'   }
#' )
#' patchwork::wrap_plots(plotlist = plist2)
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
    seed = 11) {
  log_message(
    "Start standard scop workflow...",
    verbose = verbose
  )

  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat} object",
      message_type = "error"
    )
  }

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  linear_reductions <- c(
    "pca", "svd", "ica",
    "nmf", "mds", "glmpca"
  )
  if (any(!linear_reduction %in% c(linear_reductions, SeuratObject::Reductions(srt)))) {
    log_message(
      "{.arg linear_reduction} must be one of: {.val {linear_reductions}}",
      message_type = "error"
    )
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  nonlinear_reductions <- c(
    "umap", "umap-naive", "tsne",
    "dm", "phate", "pacmap",
    "trimap", "largevis", "fr"
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
    PrepareEnv()
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

  if (normalization_method == "TFIDF") {
    log_message(
      "{.arg normalization_method} is {.pkg TFIDF}. Use {.pkg lsi} workflow",
      verbose = verbose
    )
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }
  scale_features <- GetAssayData5(
    srt,
    layer = "scale.data",
    assay = assay
  ) |>
    rownames()
  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))) {
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
    srt <- RunDimReduction(
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
      linear_reduction_dims_use_current <- 1:ncol(
        srt@reductions[[paste0(
          prefix,
          lr
        )]]@cell.embeddings
      )
      if (normalization_method == "TFIDF") {
        linear_reduction_dims_use_current <- 2:max(
          linear_reduction_dims_use_current
        )
      }
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
        log_message(error, message_type = "warning", verbose = verbose)
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
            srt <- RunDimReduction(
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
        log_message(
          error,
          message_type = "warning",
          verbose = verbose
        )
        log_message(
          "Error when performing {.pkg {nr}} nonlinear dimension reduction. Skip it",
          message_type = "error"
        )
        srt
      }
    )
  }

  if (paste0(prefix, linear_reduction[1], "clusters") %in% colnames(srt@meta.data)) {
    srt[[paste0(prefix, "clusters")]] <- srt[[paste0(
      prefix,
      linear_reduction[1],
      "clusters"
    )]]
  }
  for (nr in nonlinear_reduction) {
    for (n in nonlinear_reduction_dims) {
      if (paste0(prefix, linear_reduction[1], toupper(nr), n, "D") %in% names(srt@reductions)) {
        reduc <- srt@reductions[[paste0(
          prefix,
          linear_reduction[1],
          toupper(nr),
          n,
          "D"
        )]]
        srt@reductions[[paste0(prefix, toupper(nr), n, "D")]] <- reduc
      }
    }
    srt@misc[["Default_reduction"]] <- paste0(prefix, toupper(nr))
  }

  SeuratObject::DefaultAssay(srt) <- assay
  SeuratObject::VariableFeatures(srt) <- srt@misc[["Standard_HVF"]] <- HVF

  log_message(
    "Run scop standard workflow done",
    message_type = "success",
    verbose = verbose
  )

  return(srt)
}
