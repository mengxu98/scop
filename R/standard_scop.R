#' Standard scop
#'
#' This function performs a standard single-cell analysis workflow.
#'
#' @param srt A Seurat object.
#' @param prefix A prefix to add to the names of intermediate objects created by the function (default is "Standard").
#' @param assay The name of the assay to use for the analysis. If NULL, the default assay of the Seurat object will be used.
#' @param do_normalization A logical value indicating whether to perform normalization. If NULL, normalization will be performed if the specified assay does not have scaled data.
#' @param normalization_method The method to use for normalization. Options are "LogNormalize", "SCT", or "TFIDF" (default is "LogNormalize").
#' @param do_HVF_finding A logical value indicating whether to perform high variable feature finding. If TRUE, the function will force to find the highly variable features (HVF) using the specified HVF method.
#' @param HVF_method The method to use for finding highly variable features. Options are "vst", "mvp" or "disp" (default is "vst").
#' @param nHVF The number of highly variable features to select. If NULL, all highly variable features will be used.
#' @param HVF A vector of feature names to use as highly variable features. If NULL, the function will use the highly variable features identified by the HVF method.
#' @param do_scaling A logical value indicating whether to perform scaling. If TRUE, the function will force to scale the data using the ScaleData function.
#' @param vars_to_regress A vector of feature names to use as regressors in the scaling step. If NULL, no regressors will be used.
#' @param regression_model The regression model to use for scaling. Options are "linear", "poisson", or "negativebinomial" (default is "linear").
#' @param linear_reduction The linear dimensionality reduction method to use. Options are "pca", "svd", "ica", "nmf", "mds", or "glmpca" (default is "pca").
#' @param linear_reduction_dims The number of dimensions to keep after linear dimensionality reduction (default is 50).
#' @param linear_reduction_dims_use The dimensions to use for downstream analysis. If NULL, all dimensions will be used.
#' @param linear_reduction_params A list of parameters to pass to the linear dimensionality reduction method.
#' @param force_linear_reduction A logical value indicating whether to force linear dimensionality reduction even if the specified reduction is already present in the Seurat object.
#' @param nonlinear_reduction The nonlinear dimensionality reduction method to use. Options are "umap","umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", or "fr" (default is "umap").
#' @param nonlinear_reduction_dims The number of dimensions to keep after nonlinear dimensionality reduction. If a vector is provided, different numbers of dimensions can be specified for each method (default is c(2, 3)).
#' @param nonlinear_reduction_params A list of parameters to pass to the nonlinear dimensionality reduction method.
#' @param force_nonlinear_reduction A logical value indicating whether to force nonlinear dimensionality reduction even if the specified reduction is already present in the Seurat object.
#' @param neighbor_metric The distance metric to use for finding neighbors. Options are "euclidean", "cosine", "manhattan", or "hamming" (default is "euclidean").
#' @param neighbor_k The number of nearest neighbors to use for finding neighbors (default is 20).
#' @param cluster_algorithm The clustering algorithm to use. Options are "louvain", "slm", or "leiden" (default is "louvain").
#' @param cluster_resolution The resolution parameter to use for clustering. Larger values result in fewer clusters (default is 0.6).
#' @param seed The random seed to use for reproducibility (default is 11).
#'
#' @return A \code{Seurat} object.
#'
#' @seealso \code{\link{integration_scop}}
#'
#' @export
#'
#' @examples
#' data("pancreas_sub")
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
    seed = 11) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "'srt' is not a Seurat object.",
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
      "'linear_reduction' must be one of: ", paste(linear_reductions, ""),
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
      "'nonlinear_reduction' must be one of: ", paste(nonlinear_reductions, ""),
      message_type = "error"
    )
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    log_message(
      "'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.",
      message_type = "error"
    )
  }
  if (cluster_algorithm == "leiden") {
    check_python("leidenalg")
  }
  cluster_algorithm_index <- switch(
    EXPR = tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  time_start <- Sys.time()
  set.seed(seed)

  log_message("Start standard_scop")

  checked <- check_srt_list(
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
    seed = seed
  )
  srt <- checked[["srt_list"]][[1]]
  HVF <- checked[["HVF"]]
  assay <- checked[["assay"]]
  type <- checked[["type"]]
  rm(checked)

  if (normalization_method == "TFIDF") {
    log_message(
      "normalization_method is 'TFIDF'. Use 'lsi' workflow..."
    )
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  if (
    isTRUE(do_scaling) ||
      (is.null(do_scaling) &&
        any(
          !HVF %in%
            rownames(GetAssayData5(
              srt,
              layer = "scale.data",
              assay = assay
            ))
        ))
  ) {
    if (normalization_method != "SCT") {
      log_message(
        "Perform ScaleData on the data..."
      )
      srt <- Seurat::ScaleData(
        object = srt,
        assay = assay,
        features = HVF,
        vars.to.regress = vars_to_regress,
        model.use = regression_model,
        verbose = FALSE
      )
    }
  }

  for (lr in linear_reduction) {
    log_message(
      paste0("Perform linear dimension reduction (", lr, ") on the data...")
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
      verbose = TRUE,
      seed = seed
    )

    if (is.null(linear_reduction_dims_use)) {
      # linear_reduction_dims_use_current <- srt@reductions[[paste0(
      #   prefix,
      #   lr
      # )]]@misc[["dims_estimate"]]
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
          # force.recalc = TRUE,
          graph.name = paste0(prefix, lr, "_", c("KNN", "SNN")),
          verbose = FALSE
        )

        log_message(
          paste0("Perform FindClusters (", cluster_algorithm, ") on the data...")
        )
        srt <- Seurat::FindClusters(
          object = srt,
          resolution = cluster_resolution,
          algorithm = cluster_algorithm_index,
          graph.name = paste0(prefix, lr, "_SNN"),
          verbose = FALSE
        )
        log_message("Reorder clusters...")
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
        log_message(error, message_type = "warning")
        log_message(
          "Error when performing FindClusters. Skip this step...",
          message_type = "warning"
        )
        return(srt)
      }
    )

    srt <- tryCatch(
      {
        for (nr in nonlinear_reduction) {
          log_message(
            "Perform nonlinear dimension reduction (",
            nr,
            ") on the data..."
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
              verbose = FALSE,
              seed = seed
            )
          }
        }
        srt
      },
      error = function(error) {
        log_message(
          error,
          message_type = "warning"
        )
        log_message(
          "Error when performing nonlinear dimension reduction. Skip this step...",
          message_type = "error"
        )
        return(srt)
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

  time_end <- Sys.time()
  log_message("Run standard_scop done", message_type = "success")
  log_message(
    "Elapsed time: ",
    format(
      round(difftime(time_end, time_start), 2),
      format = "%Y-%m-%d %H:%M:%S"
    )
  )

  return(srt)
}
