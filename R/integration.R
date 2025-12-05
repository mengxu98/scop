#' @title The integration_scop function
#'
#' @description
#' Integrate single-cell RNA-seq data using various integration methods.
#'
#' @md
#' @inheritParams CheckDataList
#' @inheritParams CheckDataMerge
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param scale_within_batch  Whether to scale data within each batch.
#' Only valid when the `integration_method` is one of `"Uncorrected"`,
#' `"Seurat"`, `"MNN"`, `"Harmony"`, `"BBKNN"`, `"CSS"`, `"ComBat"`.
#' @param integration_method A character vector specifying the integration method to use.
#' Supported methods are: `"Uncorrected"`, `"Seurat"`, `"scVI"`, `"MNN"`, `"fastMNN"`,
#' `"Harmony"`, `"Scanorama"`, `"BBKNN"`, `"CSS"`, `"LIGER"`, `"Conos"`, `"ComBat"`.
#' Default is `"Uncorrected"`.
#' @param append Whether the integrated data will be appended to the original Seurat object (`srt_merge`).
#' Default is `TRUE`.
#' @param ... Additional arguments to be passed to the integration method functions.
#'
#' @return A `Seurat` object.
#'
#' @seealso
#' [Seurat_integrate],
#' [scVI_integrate],
#' [MNN_integrate],
#' [fastMNN_integrate],
#' [Harmony_integrate],
#' [Scanorama_integrate],
#' [BBKNN_integrate],
#' [CSS_integrate],
#' [LIGER_integrate],
#' [Conos_integrate],
#' [ComBat_integrate],
#' [standard_scop]
#'
#' @export
#' @examples
#' data(panc8_sub)
#' panc8_sub <- integration_scop(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_method = "Uncorrected"
#' )
#' CellDimPlot(
#'   panc8_sub,
#'   group.by = c("tech", "celltype")
#' )
#'
#' panc8_sub <- integration_scop(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_method = "Uncorrected",
#'   HVF_min_intersection = 5
#' )
#' CellDimPlot(
#'   panc8_sub,
#'   group.by = c("tech", "celltype")
#' )
#'
#' panc8_sub <- integration_scop(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_method = "Uncorrected",
#'   HVF_min_intersection = 5,
#'   scale_within_batch = TRUE
#' )
#' CellDimPlot(
#'   panc8_sub,
#'   group.by = c("tech", "celltype")
#' )
#'
#' \dontrun{
#' panc8_sub <- integration_scop(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_method = "Seurat"
#' )
#' CellDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' panc8_sub <- integration_scop(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_method = "Seurat",
#'   FindIntegrationAnchors_params = list(reduction = "rpca")
#' )
#' CellDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' integration_methods <- c(
#'   "Uncorrected", "Seurat", "scVI", "MNN", "fastMNN", "Harmony",
#'   "Scanorama", "BBKNN", "CSS", "LIGER", "Conos", "ComBat"
#' )
#' for (method in integration_methods) {
#'   panc8_sub <- integration_scop(
#'     panc8_sub,
#'     batch = "tech",
#'     integration_method = method,
#'     linear_reduction_dims_use = 1:50,
#'     nonlinear_reduction = "umap"
#'   )
#'   print(
#'     CellDimPlot(panc8_sub,
#'       group.by = c("tech", "celltype"),
#'       reduction = paste0(method, "UMAP2D"),
#'       xlab = "", ylab = "", title = method,
#'       legend.position = "none", theme_use = "theme_blank"
#'     )
#'   )
#' }
#'
#' nonlinear_reductions <- c(
#'   "umap", "tsne", "dm", "phate",
#'   "pacmap", "trimap", "largevis", "fr"
#' )
#' panc8_sub <- integration_scop(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_method = "Seurat",
#'   linear_reduction_dims_use = 1:50,
#'   nonlinear_reduction = nonlinear_reductions
#' )
#' for (nr in nonlinear_reductions) {
#'   print(
#'     CellDimPlot(
#'       panc8_sub,
#'       group.by = c("tech", "celltype"),
#'       reduction = paste0("Seurat", nr, "2D"),
#'       xlab = "", ylab = "", title = nr,
#'       legend.position = "none", theme_use = "theme_blank"
#'     )
#'   )
#' }
#' }
integration_scop <- function(
    srt_merge = NULL,
    batch,
    append = TRUE,
    srt_list = NULL,
    assay = NULL,
    integration_method = c(
      "Uncorrected",
      "Seurat",
      "scVI",
      "MNN",
      "fastMNN",
      "Harmony",
      "Scanorama",
      "BBKNN",
      "CSS",
      "LIGER",
      "Conos",
      "ComBat"
    ),
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
    seed = 11,
    verbose = TRUE,
    ...) {
  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "Neither {.arg srt_list} nor {.arg srt_merge} was found",
      message_type = "error"
    )
  }
  integration_method <- match.arg(integration_method)

  args <- as.list(match.call())[-1]
  new_env <- new.env(parent = parent.frame())
  args <- lapply(args, function(x) eval(x, envir = new_env))

  formals <- mget(names(formals()))
  formals <- formals[names(formals) != "..."]

  args <- utils::modifyList(formals, args)

  log_message(
    "Run {.pkg {integration_method}} integration...",
    message_type = "running",
    verbose = verbose
  )
  integrate_fun <- paste0(integration_method, "_integrate")
  srt_integrated <- invoke_fun(
    integrate_fun,
    args[names(args) %in% methods::formalArgs(integrate_fun)]
  )
  log_message(
    "Run {.pkg {integration_method}} integration done",
    message_type = "success",
    verbose = verbose
  )

  return(srt_integrated)
}

#' @title The Uncorrected integration function
#'
#' @inheritParams integration_scop
#'
#' @export
Uncorrected_integrate <- function(
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
    verbose = TRUE,
    seed = 11) {
  if (length(linear_reduction) > 1) {
    log_message(
      "Only the first method in the {.arg linear_reduction} will be used",
      message_type = "warning",
      verbose = verbose
    )
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srt_merge)) {
    reduc_test <- c(reduc_test, SeuratObject::Reductions(srt_merge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    log_message(
      "{.arg linear_reduction} must be one of {.val {reduc_test}}",
      message_type = "error"
    )
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
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
      "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
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
  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "srt_list and srt_merge were all empty",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "srt_list and srt_merge have different cells",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
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
      verbose = verbose,
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

  if (normalization_method == "TFIDF") {
    log_message(
      "{.arg normalization_method} is {.val TFIDF}. Use {.pkg lsi} workflow..."
    )
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  log_message(
    "Perform {.pkg Uncorrected} integration"
  )
  scale_features <- rownames(
    GetAssayData5(
      srt_merge,
      layer = "scale.data",
      assay = SeuratObject::DefaultAssay(srt_merge)
    )
  )
  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))) {
    if (normalization_method != "SCT") {
      log_message(
        "Perform {.fn Seurat::ScaleData}",
        verbose = verbose
      )
      srt_merge <- Seurat::ScaleData(
        object = srt_merge,
        split.by = if (isTRUE(scale_within_batch)) batch else NULL,
        assay = SeuratObject::DefaultAssay(srt_merge),
        features = HVF,
        vars.to.regress = vars_to_regress,
        model.use = regression_model,
        verbose = FALSE
      )
    }
  }

  log_message(
    "Perform linear dimension reduction({.val {linear_reduction}})",
    verbose = verbose
  )
  srt_merge <- RunDimReduction(
    srt_merge,
    prefix = "Uncorrected",
    features = HVF,
    assay = SeuratObject::DefaultAssay(srt_merge),
    linear_reduction = linear_reduction,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = force_linear_reduction,
    verbose = FALSE,
    seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srt_merge@reductions[[paste0(
      "Uncorrected",
      linear_reduction
    )]]@misc[["dims_estimate"]]
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  srt_merge <- tryCatch(
    {
      srt_merge <- Seurat::FindNeighbors(
        object = srt_merge,
        reduction = paste0("Uncorrected", linear_reduction),
        dims = linear_reduction_dims_use,
        annoy.metric = neighbor_metric,
        k.param = neighbor_k,
        graph.name = paste0("Uncorrected_", c("KNN", "SNN")),
        verbose = FALSE
      )

      log_message(
        "Perform Seurat::FindClusters ({.val {cluster_algorithm}})"
      )
      srt_merge <- Seurat::FindClusters(
        object = srt_merge,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        leiden_method = "igraph",
        graph.name = "Uncorrected_SNN",
        verbose = FALSE
      )
      log_message("Reorder clusters...")
      srt_merge <- srt_reorder(
        srt_merge,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srt_merge[["seurat_clusters"]] <- NULL
      srt_merge[[paste0("Uncorrected", linear_reduction, "clusters")]] <- SeuratObject::Idents(
        srt_merge
      )
      srt_merge
    },
    error = function(error) {
      log_message(error, message_type = "warning")
      log_message(
        "Error when performing FindClusters. Skip this step",
        message_type = "warning"
      )
      srt_merge
    }
  )

  srt_merge <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        log_message(
          "Perform nonlinear dimension reduction ({.val {nr}})"
        )
        for (n in nonlinear_reduction_dims) {
          srt_merge <- RunDimReduction(
            srt_merge,
            prefix = "Uncorrected",
            reduction_use = paste0("Uncorrected", linear_reduction),
            reduction_dims = linear_reduction_dims_use,
            graph_use = "Uncorrected_SNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srt_merge
    },
    error = function(error) {
      log_message(error, message_type = "warning")
      log_message(
        "Error when performing nonlinear dimension reduction. Skip this step",
        message_type = "warning"
      )
      srt_merge
    }
  )

  SeuratObject::DefaultAssay(srt_merge) <- assay
  SeuratObject::VariableFeatures(srt_merge) <- srt_merge@misc[["Uncorrected_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_merge,
      pattern = paste0(assay, "|Uncorrected|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_merge)
  }
}

#' @title The Seurat integration function
#'
#' @inheritParams integration_scop
#' @param FindIntegrationAnchors_params A list of parameters for the Seurat::FindIntegrationAnchors function.
#' Default is `list()`.
#' @param IntegrateData_params A list of parameters for the Seurat::IntegrateData function.
#' Default is `list()`.
#' @param IntegrateEmbeddings_params A list of parameters for the Seurat::IntegrateEmbeddings function.
#' Default is `list()`.
#'
#' @export
Seurat_integrate <- function(
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
    FindIntegrationAnchors_params = list(),
    IntegrateData_params = list(),
    IntegrateEmbeddings_params = list(),
    verbose = TRUE,
    seed = 11) {
  if (length(linear_reduction) > 1) {
    log_message(
      "Only the first method in the {.arg linear_reduction} will be used",
      message_type = "warning",
      verbose = verbose
    )
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "ica", "svd", "nmf", "mds", "glmpca")
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
      "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
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
  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "srt_list and srt_merge were all empty",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "srt_list and srt_merge have different cells",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
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
    if (normalization_method == "TFIDF") {
      srt_merge <- Reduce(merge, srt_list)
      SeuratObject::VariableFeatures(srt_merge) <- HVF
    }
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

  if (min(sapply(srt_list, ncol)) < 50) {
    log_message(
      "The cell count in some batches is lower than 50, which may not be suitable for the current integration method",
      message_type = "warning"
    )
    answer <- utils::askYesNo("Are you sure to continue?", default = FALSE)
    if (isFALSE(answer)) {
      return(srt_merge)
    }
  }

  if (normalization_method == "TFIDF") {
    log_message(
      "{.arg normalization_method} is {.val TFIDF}. Use {.pkg rlsi} integration workflow..."
    )
    do_scaling <- FALSE
    linear_reduction <- "svd"
    FindIntegrationAnchors_params[["reduction"]] <- "rlsi"
    if (is.null(FindIntegrationAnchors_params[["dims"]])) {
      FindIntegrationAnchors_params[["dims"]] <- 2:min(
        linear_reduction_dims,
        30
      )
    }
    srt_merge <- Signac::RunTFIDF(
      object = srt_merge,
      assay = SeuratObject::DefaultAssay(srt_merge),
      verbose = FALSE
    )
    srt_merge <- RunDimReduction(
      srt_merge,
      prefix = "",
      features = HVF,
      assay = SeuratObject::DefaultAssay(srt_merge),
      linear_reduction = "svd",
      linear_reduction_dims = linear_reduction_dims,
      linear_reduction_params = linear_reduction_params,
      force_linear_reduction = force_linear_reduction,
      verbose = FALSE,
      seed = seed
    )
    srt_merge[["lsi"]] <- srt_merge[["svd"]]
    for (i in seq_along(srt_list)) {
      srt <- srt_list[[i]]
      log_message(
        "Perform linear dimension reduction (svd) on the data {.val {i}} ..."
      )
      srt <- RunDimReduction(
        srt,
        prefix = "",
        features = HVF,
        assay = SeuratObject::DefaultAssay(srt),
        linear_reduction = "svd",
        linear_reduction_dims = linear_reduction_dims,
        linear_reduction_params = linear_reduction_params,
        force_linear_reduction = force_linear_reduction,
        verbose = FALSE,
        seed = seed
      )
      srt[["lsi"]] <- srt[["svd"]]
      srt_list[[i]] <- srt
    }
  }

  if (isTRUE(FindIntegrationAnchors_params[["reduction"]] == "rpca")) {
    log_message("Use {.pkg rpca} integration workflow...")
    for (i in seq_along(srt_list)) {
      srt <- srt_list[[i]]
      scale_features <- rownames(
        GetAssayData5(
          srt,
          layer = "scale.data",
          assay = SeuratObject::DefaultAssay(srt)
        )
      )
      if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))) {
        log_message(
          "Perform ScaleData on the data {.val {i}} ..."
        )
        srt <- Seurat::ScaleData(
          object = srt,
          assay = SeuratObject::DefaultAssay(srt),
          features = HVF,
          vars.to.regress = vars_to_regress,
          model.use = regression_model,
          verbose = FALSE
        )
      }
      log_message(
        "Perform linear dimension reduction (pca) on the data {.val {i}} ..."
      )
      srt <- RunDimReduction(
        srt,
        prefix = "",
        features = HVF,
        assay = SeuratObject::DefaultAssay(srt),
        linear_reduction = "pca",
        linear_reduction_dims = linear_reduction_dims,
        linear_reduction_params = linear_reduction_params,
        force_linear_reduction = force_linear_reduction,
        verbose = FALSE,
        seed = seed
      )
      srt_list[[i]] <- srt
    }
  }

  if (is.null(names(srt_list))) {
    names(srt_list) <- paste0("srt_", seq_along(srt_list))
  }

  if (normalization_method %in% c("LogNormalize", "SCT")) {
    log_message("Perform FindIntegrationAnchors")
    params1 <- list(
      object.list = srt_list,
      normalization.method = normalization_method,
      anchor.features = HVF,
      verbose = FALSE
    )
    for (nm in names(FindIntegrationAnchors_params)) {
      params1[[nm]] <- FindIntegrationAnchors_params[[nm]]
    }
    srt_anchors <- invoke_fun(
      Seurat::FindIntegrationAnchors,
      params1
    )

    log_message("Perform integration(Seurat)")
    params2 <- list(
      anchorset = srt_anchors,
      new.assay.name = "Seuratcorrected",
      normalization.method = normalization_method,
      features.to.integrate = HVF,
      verbose = FALSE
    )
    for (nm in names(IntegrateData_params)) {
      params2[[nm]] <- IntegrateData_params[[nm]]
    }
    srt_integrated <- invoke_fun(Seurat::IntegrateData, params2)

    SeuratObject::DefaultAssay(srt_integrated) <- "Seuratcorrected"
    SeuratObject::VariableFeatures(srt_integrated[["Seuratcorrected"]]) <- HVF

    scale_features <- rownames(
      GetAssayData5(
        srt_integrated,
        layer = "scale.data",
        assay = SeuratObject::DefaultAssay(srt_integrated)
      )
    )
    if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))) {
      log_message("Perform ScaleData")
      srt_integrated <- Seurat::ScaleData(
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
      "Perform linear dimension reduction ({.val {linear_reduction}})"
    )
    srt_integrated <- RunDimReduction(
      srt_integrated,
      prefix = "Seurat",
      features = HVF,
      assay = SeuratObject::DefaultAssay(srt_integrated),
      linear_reduction = linear_reduction,
      linear_reduction_dims = linear_reduction_dims,
      linear_reduction_params = linear_reduction_params,
      force_linear_reduction = force_linear_reduction,
      verbose = FALSE,
      seed = seed
    )
    if (is.null(linear_reduction_dims_use)) {
      linear_reduction_dims_use <- srt_integrated@reductions[[paste0(
        "Seurat",
        linear_reduction
      )]]@misc[["dims_estimate"]] %||%
        1:linear_reduction_dims
    }
  } else if (normalization_method == "TFIDF") {
    log_message("Perform FindIntegrationAnchors")
    params1 <- list(
      object.list = srt_list,
      normalization.method = "LogNormalize",
      anchor.features = HVF,
      reduction = "rlsi",
      verbose = FALSE
    )
    for (nm in names(FindIntegrationAnchors_params)) {
      params1[[nm]] <- FindIntegrationAnchors_params[[nm]]
    }
    srt_anchors <- invoke_fun(Seurat::FindIntegrationAnchors, params1)

    log_message("Perform integration(Seurat)")
    params2 <- list(
      anchorset = srt_anchors,
      reductions = srt_merge[["lsi"]],
      new.reduction.name = "Seuratlsi",
      verbose = FALSE
    )
    for (nm in names(IntegrateEmbeddings_params)) {
      params2[[nm]] <- IntegrateEmbeddings_params[[nm]]
    }
    srt_integrated <- invoke_fun(IntegrateEmbeddings, params2)

    if (is.null(linear_reduction_dims_use)) {
      linear_reduction_dims_use <- 2:max(srt_integrated@reductions[[paste0(
        "Seurat",
        linear_reduction
      )]]@misc[["dims_estimate"]]) %||%
        2:linear_reduction_dims
    }
    linear_reduction <- "lsi"
  }

  srt_integrated <- tryCatch(
    {
      srt_integrated <- Seurat::FindNeighbors(
        object = srt_integrated,
        reduction = paste0("Seurat", linear_reduction),
        dims = linear_reduction_dims_use,
        annoy.metric = neighbor_metric,
        k.param = neighbor_k,
        graph.name = paste0("Seurat_", c("KNN", "SNN")),
        verbose = FALSE
      )

      log_message(
        "Perform Seurat::FindClusters ({.val {cluster_algorithm}})"
      )
      srt_integrated <- Seurat::FindClusters(
        object = srt_integrated,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        leiden_method = "igraph",
        graph.name = "Seurat_SNN",
        verbose = FALSE
      )
      log_message("Reorder clusters...")
      srt_integrated <- srt_reorder(
        srt_integrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srt_integrated[["seurat_clusters"]] <- NULL
      srt_integrated[["Seuratclusters"]] <- SeuratObject::Idents(srt_integrated)
      srt_integrated
    },
    error = function(error) {
      log_message(error, message_type = "warning")
      log_message(
        "Error when performing {.fn Seurat::FindClusters}. Skip this step",
        message_type = "warning"
      )
      srt_integrated
    }
  )

  srt_integrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        log_message(
          "Perform nonlinear dimension reduction ({.pkg {nr}})"
        )
        for (n in nonlinear_reduction_dims) {
          srt_integrated <- RunDimReduction(
            srt_integrated,
            prefix = "Seurat",
            reduction_use = paste0("Seurat", linear_reduction),
            reduction_dims = linear_reduction_dims_use,
            graph_use = "Seurat_SNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srt_integrated
    },
    error = function(error) {
      log_message(error, message_type = "warning")
      log_message(
        "Error when performing nonlinear dimension reduction. Skip this step",
        message_type = "warning"
      )
      return(srt_integrated)
    }
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[["Seurat_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_integrated,
      pattern = paste0(assay, "|Seurat|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_integrated)
  }
}

#' @title The scVI integration function
#'
#' @inheritParams integration_scop
#' @param scVI_dims_use A vector specifying the dimensions returned by scVI that will be utilized for downstream cell cluster finding and non-linear reduction.
#' If set to NULL, all the returned dimensions will be used by default.
#' @param model A string indicating the scVI model to be used.
#' Options are "SCVI" and "PEAKVI".
#' Default is `"SCVI"`.
#' @param SCVI_params A list of parameters for the SCVI model.
#' Default is `list()`.
#' @param PEAKVI_params A list of parameters for the PEAKVI model.
#' Default is `list()`.
#' @param num_threads An integer setting the number of threads for scVI.
#' Default is `8`.
#'
#' @export
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
    num_threads = 1,
    verbose = TRUE,
    seed = 11) {
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
      "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
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
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  if (
    .Platform$OS.type == "windows" &&
      !exist_python_pkgs(packages = "scvi-tools")
  ) {
    suppressWarnings(
      system2(
        command = conda_python(),
        args = "-m pip install jax[cpu]===0.3.20 -f https://whls.blob.core.windows.net/unstable/index.html --use-deprecated legacy-resolver",
        stdout = TRUE
      )
    )
  }

  PrepareEnv()
  check_python("scvi-tools")
  scvi <- reticulate::import("scvi")
  scipy <- reticulate::import("scipy")
  set.seed(seed)

  scvi$settings$num_threads <- as.integer(num_threads)

  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "srt_list and srt_merge were all empty.",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "srt_list and srt_merge have different cells.",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
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
    params <- list(
      adata = adata
    )
    for (nm in names(SCVI_params)) {
      params[[nm]] <- SCVI_params[[nm]]
    }
    model <- invoke_fun(scvi$model$SCVI, params)
    model$train()
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
    params <- list(
      adata = adata
    )
    for (nm in names(PEAKVI_params)) {
      params[[nm]] <- PEAKVI_params[[nm]]
    }
    model <- invoke_fun(scvi$model$PEAKVI, params)
    model$train()
    srt_integrated <- srt_merge
    srt_merge <- NULL
  }

  latent <- as_matrix(model$get_latent_representation())
  rownames(latent) <- colnames(srt_integrated)
  colnames(latent) <- paste0("scVI_", seq_len(ncol(latent)))
  srt_integrated[["scVI"]] <- CreateDimReducObject(
    embeddings = latent,
    key = "scVI_",
    assay = SeuratObject::DefaultAssay(srt_integrated)
  )
  if (is.null(scVI_dims_use)) {
    scVI_dims_use <- 1:ncol(srt_integrated[["scVI"]]@cell.embeddings)
  }

  srt_integrated <- tryCatch(
    {
      srt_integrated <- Seurat::FindNeighbors(
        object = srt_integrated,
        reduction = "scVI",
        dims = scVI_dims_use,
        annoy.metric = neighbor_metric,
        k.param = neighbor_k,
        graph.name = paste0("scVI_", c("KNN", "SNN")),
        verbose = FALSE
      )

      log_message(
        "Perform Seurat::FindClusters ({.val {cluster_algorithm}})"
      )
      srt_integrated <- Seurat::FindClusters(
        object = srt_integrated,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        leiden_method = "igraph",
        graph.name = "scVI_SNN",
        verbose = FALSE
      )
      log_message("Reorder clusters...")
      srt_integrated <- srt_reorder(
        srt_integrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srt_integrated[["seurat_clusters"]] <- NULL
      srt_integrated[["scVIclusters"]] <- SeuratObject::Idents(srt_integrated)
      srt_integrated
    },
    error = function(error) {
      log_message(error, message_type = "warning")
      log_message(
        "Error when performing {.fn Seurat::FindClusters}. Skip this step",
        message_type = "warning"
      )
      srt_integrated
    }
  )

  srt_integrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        log_message(
          paste0("Perform nonlinear dimension reduction (", nr, ")")
        )
        for (n in nonlinear_reduction_dims) {
          srt_integrated <- RunDimReduction(
            srt_integrated,
            prefix = "scVI",
            reduction_use = "scVI",
            reduction_dims = scVI_dims_use,
            graph_use = "scVI_SNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srt_integrated
    },
    error = function(error) {
      log_message(error, message_type = "warning")
      log_message(
        "Error when performing nonlinear dimension reduction. Skip this step",
        message_type = "warning"
      )
      return(srt_integrated)
    }
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[["scVI_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_integrated,
      pattern = paste0(assay, "|scVI|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_integrated)
  }
}

#' @title The MNN integration function
#'
#' @inheritParams integration_scop
#' @param mnnCorrect_params A list of parameters for the batchelor::mnnCorrect function,
#' default is an empty list.
#'
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
    seed = 11) {
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
      "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
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
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_r("batchelor")
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "srt_list and srt_merge were all empty.",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "srt_list and srt_merge have different cells.",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
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
      object = srt_merge, split.by = batch
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

  if (normalization_method == "TFIDF") {
    log_message(
      "{.arg normalization_method} is {.val TFIDF}. Use {.pkg lsi} workflow..."
    )
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }
  sce_list <- lapply(
    srt_list, function(srt) {
      data_matrix <- GetAssayData5(
        srt,
        layer = "data",
        assay = SeuratObject::DefaultAssay(srt)
      )
      sce <- Seurat::as.SingleCellExperiment(
        Seurat::CreateSeuratObject(
          counts = data_matrix[HVF, ]
        )
      )
      if (inherits(sce@assays@data$logcounts, "dgCMatrix")) {
        sce@assays@data$logcounts <- as_matrix(
          sce@assays@data$logcounts
        )
      }
      return(sce)
    }
  )
  if (is.null(names(sce_list))) {
    names(sce_list) <- paste0("sce_", seq_along(sce_list))
  }

  log_message("Perform {.pkg MNN} integration")
  params <- list(
    sce_list,
    cos.norm.out = FALSE
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
  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))) {
    log_message("Perform ScaleData")
    srt_integrated <- Seurat::ScaleData(
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
    paste0("Perform linear dimension reduction (", linear_reduction, ")")
  )
  srt_integrated <- RunDimReduction(
    srt_integrated,
    prefix = "MNN",
    features = HVF,
    assay = SeuratObject::DefaultAssay(srt_integrated),
    linear_reduction = linear_reduction,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = force_linear_reduction,
    verbose = FALSE,
    seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srt_integrated@reductions[[paste0(
      "MNN",
      linear_reduction
    )]]@misc[["dims_estimate"]] %||%
      1:linear_reduction_dims
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  srt_integrated <- tryCatch(
    {
      srt_integrated <- Seurat::FindNeighbors(
        object = srt_integrated,
        reduction = paste0("MNN", linear_reduction),
        dims = linear_reduction_dims_use,
        annoy.metric = neighbor_metric,
        k.param = neighbor_k,
        graph.name = paste0("MNN_", c("KNN", "SNN")),
        verbose = FALSE
      )

      log_message(
        paste0("Perform Seurat::FindClusters (", cluster_algorithm, ")")
      )
      srt_integrated <- Seurat::FindClusters(
        object = srt_integrated,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        leiden_method = "igraph",
        graph.name = "MNN_SNN",
        verbose = FALSE
      )
      log_message("Reorder clusters...")
      srt_integrated <- srt_reorder(
        srt_integrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srt_integrated[["seurat_clusters"]] <- NULL
      srt_integrated[["MNNclusters"]] <- SeuratObject::Idents(srt_integrated)
      srt_integrated
    },
    error = function(error) {
      log_message(error, message_type = "warning")
      log_message(
        "Error when performing FindClusters. Skip this step",
        message_type = "warning"
      )
      srt_integrated
    }
  )

  srt_integrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        log_message(
          paste0("Perform nonlinear dimension reduction (", nr, ")")
        )
        for (n in nonlinear_reduction_dims) {
          srt_integrated <- RunDimReduction(
            srt_integrated,
            prefix = "MNN",
            reduction_use = paste0("MNN", linear_reduction),
            reduction_dims = linear_reduction_dims_use,
            graph_use = "MNN_SNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srt_integrated
    },
    error = function(error) {
      log_message(error, message_type = "warning")
      log_message(
        "Error when performing nonlinear dimension reduction. Skip this step",
        message_type = "warning"
      )
      return(srt_integrated)
    }
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[["MNN_HVF"]] <- HVF

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

#' @title The fastMNN integration function
#'
#' @inheritParams integration_scop
#' @param fastMNN_dims_use A vector specifying the dimensions returned by fastMNN that will be utilized for downstream cell cluster finding and non-linear reduction.
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
    seed = 11) {
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
  if (cluster_algorithm == "leiden") {
    PrepareEnv()
    check_python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_r("batchelor")
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "srt_list and srt_merge were all empty.",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "srt_list and srt_merge have different cells.",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
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

  sce_list <- lapply(srt_list, function(srt) {
    sce <- Seurat::as.SingleCellExperiment(
      Seurat::CreateSeuratObject(
        counts = GetAssayData5(
          srt,
          layer = "data",
          assay = SeuratObject::DefaultAssay(srt)
        )[HVF, , drop = FALSE]
      )
    )
    if (inherits(sce@assays@data$logcounts, "dgCMatrix")) {
      sce@assays@data$logcounts <- as_matrix(sce@assays@data$logcounts)
    }
    return(sce)
  })
  if (is.null(names(sce_list))) {
    names(sce_list) <- paste0("sce_", seq_along(sce_list))
  }

  log_message("Perform integration(fastMNN)")
  params <- list(
    sce_list
  )
  for (nm in names(fastMNN_params)) {
    params[[nm]] <- fastMNN_params[[nm]]
  }
  out <- invoke_fun(batchelor::fastMNN, params)

  srt_integrated <- srt_merge
  srt_merge <- NULL
  srt_integrated[["fastMNNcorrected"]] <- CreateAssayObject(
    counts = as_matrix(out@assays@data$reconstructed)
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

  srt_integrated <- tryCatch(
    {
      srt_integrated <- Seurat::FindNeighbors(
        object = srt_integrated,
        reduction = "fastMNN",
        dims = fastMNN_dims_use,
        annoy.metric = neighbor_metric,
        k.param = neighbor_k,
        graph.name = paste0("fastMNN", "_", c("KNN", "SNN")),
        verbose = FALSE
      )

      log_message(
        paste0("Perform Seurat::FindClusters (", cluster_algorithm, ")")
      )
      srt_integrated <- Seurat::FindClusters(
        object = srt_integrated,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        leiden_method = "igraph",
        graph.name = "fastMNN_SNN",
        verbose = FALSE
      )
      log_message("Reorder clusters...")
      srt_integrated <- srt_reorder(
        srt_integrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srt_integrated[["seurat_clusters"]] <- NULL
      srt_integrated[["fastMNNclusters"]] <- SeuratObject::Idents(srt_integrated)
      srt_integrated
    },
    error = function(error) {
      log_message(error, message_type = "warning")
      log_message(
        "Error when performing {.fn Seurat::FindClusters}. Skip this step",
        message_type = "warning"
      )
      srt_integrated
    }
  )

  srt_integrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        log_message(
          paste0("Perform nonlinear dimension reduction (", nr, ")")
        )
        for (n in nonlinear_reduction_dims) {
          srt_integrated <- RunDimReduction(
            srt_integrated,
            prefix = "fastMNN",
            reduction_use = "fastMNN",
            reduction_dims = fastMNN_dims_use,
            graph_use = "fastMNN_SNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srt_integrated
    },
    error = function(error) {
      log_message(error, message_type = "warning")
      log_message(
        "Error when performing nonlinear dimension reduction. Skip this step",
        message_type = "warning"
      )
      srt_integrated
    }
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[["fastMNN_HVF"]] <- HVF

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

#' @title The Harmony integration function
#'
#' @inheritParams integration_scop
#' @param harmony_dims_use A vector specifying the dimensions returned by RunHarmony that will be utilized for downstream cell cluster finding and non-linear reduction.
#' If set to NULL, all the returned dimensions will be used by default.
#' @param RunHarmony_params A list of parameters for [harmony::RunHarmony], default is an empty list.
#'
#' @export
Harmony_integrate <- function(
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
    harmony_dims_use = NULL,
    nonlinear_reduction = "umap",
    nonlinear_reduction_dims = c(2, 3),
    nonlinear_reduction_params = list(),
    force_nonlinear_reduction = TRUE,
    neighbor_metric = "euclidean",
    neighbor_k = 20L,
    cluster_algorithm = "louvain",
    cluster_resolution = 0.6,
    RunHarmony_params = list(),
    verbose = TRUE,
    seed = 11) {
  if (length(linear_reduction) > 1) {
    log_message(
      "Only the first of {.val {linear_reduction}} will be used",
      message_type = "warning",
      verbose = verbose
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
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
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
      "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
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

  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "{.arg srt_list} and {.arg srt_merge} were all empty",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "{.arg srt_list} and {.arg srt_merge} have different cells",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
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
  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))) {
    log_message("Perform {.fn Seurat::ScaleData}")
    srt_merge <- Seurat::ScaleData(
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
    "Perform linear dimension reduction({.val {linear_reduction}})"
  )
  srt_merge <- RunDimReduction(
    srt_merge,
    prefix = "Harmony",
    features = HVF,
    assay = SeuratObject::DefaultAssay(srt_merge),
    linear_reduction = linear_reduction,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = force_linear_reduction,
    verbose = FALSE,
    seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srt_merge@reductions[[paste0(
      "Harmony",
      linear_reduction
    )]]@misc[["dims_estimate"]] %||%
      1:linear_reduction_dims
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  log_message("Perform {.pkg Harmony} integration")
  log_message(
    "Harmony integration using {.val {paste0('Harmony', linear_reduction)}}, dims:{.val {min(linear_reduction_dims_use)}}-{.val {max(linear_reduction_dims_use)}} as input"
  )
  params <- list(
    object = srt_merge,
    group.by.vars = batch,
    assay = SeuratObject::DefaultAssay(srt_merge),
    reduction = paste0("Harmony", linear_reduction),
    dims.use = linear_reduction_dims_use,
    reduction.save = "Harmony",
    verbose = FALSE
  )
  feature_num <- nrow(
    GetAssayData5(
      srt_merge,
      layer = "scale.data",
      assay = SeuratObject::DefaultAssay(srt_merge)
    )
  )
  if (feature_num == 0) {
    params[["project.dim"]] <- FALSE
  }
  for (nm in names(RunHarmony_params)) {
    params[[nm]] <- RunHarmony_params[[nm]]
  }
  srt_integrated <- invoke_fun(RunHarmony2, params)

  if (is.null(harmony_dims_use)) {
    harmony_dims_use <- seq_len(
      ncol(
        srt_integrated[["Harmony"]]@cell.embeddings
      )
    )
  }

  srt_integrated <- tryCatch(
    {
      srt_integrated <- Seurat::FindNeighbors(
        object = srt_integrated,
        reduction = "Harmony",
        dims = harmony_dims_use,
        annoy.metric = neighbor_metric,
        k.param = neighbor_k,
        graph.name = paste0("Harmony", "_", c("KNN", "SNN")),
        verbose = FALSE
      )

      log_message(
        paste0("Perform Seurat::FindClusters (", cluster_algorithm, ")")
      )
      srt_integrated <- Seurat::FindClusters(
        object = srt_integrated,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        leiden_method = "igraph",
        graph.name = "Harmony_SNN",
        verbose = FALSE
      )
      log_message("Reorder clusters...")
      srt_integrated <- srt_reorder(
        srt_integrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srt_integrated[["seurat_clusters"]] <- NULL
      srt_integrated[["Harmonyclusters"]] <- SeuratObject::Idents(srt_integrated)
      srt_integrated
    },
    error = function(error) {
      log_message(error, message_type = "warning")
      log_message(
        "Error when performing FindClusters. Skip this step",
        message_type = "warning"
      )
      return(srt_integrated)
    }
  )

  srt_integrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        log_message(
          paste0("Perform nonlinear dimension reduction (", nr, ")")
        )
        for (n in nonlinear_reduction_dims) {
          srt_integrated <- RunDimReduction(
            srt_integrated,
            prefix = "Harmony",
            reduction_use = "Harmony",
            reduction_dims = harmony_dims_use,
            graph_use = "Harmony_SNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            force_nonlinear_reduction = force_nonlinear_reduction,
            nonlinear_reduction_params = nonlinear_reduction_params,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srt_integrated
    },
    error = function(error) {
      log_message(error, message_type = "warning")
      log_message(
        "Error when performing nonlinear dimension reduction. Skip this step",
        message_type = "warning"
      )
      return(srt_integrated)
    }
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[["Harmony_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_integrated,
      pattern = paste0(assay, "|Harmony|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_integrated)
  }
}

#' @title The Scanorama integration function
#'
#' @inheritParams integration_scop
#' @param Scanorama_dims_use  A vector specifying the dimensions returned by Scanorama that will be utilized for downstream cell cluster finding and non-linear reduction.
#' If set to NULL, all the returned dimensions will be used by default.
#' @param return_corrected Logical indicating whether to return the corrected data.
#' Default is `FALSE`.
#' @param Scanorama_params A list of parameters for the scanorama.correct function.
#' Default is `list()`.
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
    seed = 11) {
  PrepareEnv()

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
      "'nonlinear_reduction' must be one of ", paste(nonlinear_reductions, collapse = ", "),
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

  check_python("scanorama")
  scanorama <- reticulate::import("scanorama")
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "srt_list and srt_merge were all empty.",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "srt_list and srt_merge have different cells.",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
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

  log_message("Perform integration(Scanorama)")
  assaylist <- list()
  genelist <- list()
  for (i in seq_along(srt_list)) {
    assaylist[[i]] <- Matrix::t(
      as_matrix(
        GetAssayData5(
          object = srt_list[[i]],
          layer = "data",
          assay = SeuratObject::DefaultAssay(srt_list[[i]])
        )[HVF, , drop = FALSE]
      )
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
    srt_integrated[["Scanoramacorrected"]] <- CreateAssayObject(data = cor_value)
    SeuratObject::VariableFeatures(srt_integrated[["Scanoramacorrected"]]) <- HVF

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

  srt_integrated <- tryCatch(
    {
      srt_integrated <- Seurat::FindNeighbors(
        object = srt_integrated,
        reduction = "Scanorama",
        dims = Scanorama_dims_use,
        annoy.metric = neighbor_metric,
        k.param = neighbor_k,
        graph.name = paste0("Scanorama_", c("KNN", "SNN")),
        verbose = FALSE
      )

      log_message(
        paste0("Perform Seurat::FindClusters (", cluster_algorithm, ")")
      )
      srt_integrated <- Seurat::FindClusters(
        object = srt_integrated,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        leiden_method = "igraph",
        graph.name = "Scanorama_SNN",
        verbose = FALSE
      )
      log_message("Reorder clusters...")
      srt_integrated <- srt_reorder(
        srt_integrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srt_integrated[["seurat_clusters"]] <- NULL
      srt_integrated[["Scanoramaclusters"]] <- SeuratObject::Idents(srt_integrated)
      srt_integrated
    },
    error = function(error) {
      log_message(error, message_type = "warning")
      log_message(
        "Error when performing {.fn Seurat::FindClusters}. Skip this step",
        message_type = "warning"
      )
      srt_integrated
    }
  )

  srt_integrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        log_message(
          paste0("Perform nonlinear dimension reduction (", nr, ")")
        )
        for (n in nonlinear_reduction_dims) {
          srt_integrated <- RunDimReduction(
            srt_integrated,
            prefix = "Scanorama",
            reduction_use = "Scanorama",
            reduction_dims = Scanorama_dims_use,
            graph_use = "Scanorama_SNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srt_integrated
    },
    error = function(error) {
      log_message(error, message_type = "warning")
      log_message(
        "Error when performing nonlinear dimension reduction. Skip this step",
        message_type = "warning"
      )
      return(srt_integrated)
    }
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

#' @title The BBKNN integration function
#'
#' @inheritParams integration_scop
#' @param bbknn_params A list of parameters for the bbknn.matrix.bbknn function, default is an empty list.
#'
#' @export
BBKNN_integrate <- function(
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
    cluster_algorithm = "louvain",
    cluster_resolution = 0.6,
    bbknn_params = list(),
    verbose = TRUE,
    seed = 11) {
  PrepareEnv()

  if (length(linear_reduction) > 1) {
    log_message(
      "Only the first method in the {.arg linear_reduction} will be used",
      message_type = "warning"
    )
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c(
    "pca", "svd", "ica", "nmf", "mds", "glmpca"
  )
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
  nonlinear_reductions <- c("umap", "umap-naive", "fr")
  if (any(!nonlinear_reduction %in% nonlinear_reductions)) {
    log_message(
      "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
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
    check_python("leidenalg")
  }
  cluster_algorithm_index <- switch(
    EXPR = tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_python("bbknn")
  bbknn <- reticulate::import("bbknn")
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "{.arg srt_list} and {.arg srt_merge} were all empty",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "{.arg srt_list} and {.arg srt_merge} have different cells",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
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

  if (normalization_method == "TFIDF") {
    log_message(
      "{.arg normalization_method} is {.pkg TFIDF}. Use {.pkg lsi} workflow..."
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
  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))) {
    log_message("Perform {.fn Seurat::ScaleData}")
    srt_merge <- Seurat::ScaleData(
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
    "Perform linear dimension reduction({.val {linear_reduction}})"
  )
  srt_merge <- RunDimReduction(
    srt_merge,
    prefix = "BBKNN",
    features = HVF,
    assay = SeuratObject::DefaultAssay(srt_merge),
    linear_reduction = linear_reduction,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = force_linear_reduction,
    verbose = FALSE,
    seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srt_merge@reductions[[paste0(
      "BBKNN",
      linear_reduction
    )]]@misc[["dims_estimate"]]
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  log_message(
    "Perform integration(BBKNN)"
  )
  log_message(
    "Using Reduction({.val {paste0('BBKNN', linear_reduction)}}, dims:{.val {min(linear_reduction_dims_use)}}-{.val {max(linear_reduction_dims_use)}}) as input"
  )
  emb <- Embeddings(srt_merge, reduction = paste0("BBKNN", linear_reduction))[,
    linear_reduction_dims_use,
    drop = FALSE
  ]
  params <- list(
    pca = emb,
    batch_list = srt_merge[[batch, drop = TRUE]]
  )
  for (nm in names(bbknn_params)) {
    params[[nm]] <- bbknn_params[[nm]]
  }
  bem <- invoke_fun(bbknn$matrix$bbknn, params)
  n.neighbors <- bem[[3]]$n_neighbors
  srt_integrated <- srt_merge

  bbknn_graph <- SeuratObject::as.sparse(
    bem[[2]][1:nrow(bem[[2]]), , drop = FALSE]
  )
  rownames(bbknn_graph) <- colnames(bbknn_graph) <- rownames(emb)
  bbknn_graph <- SeuratObject::as.Graph(bbknn_graph)
  bbknn_graph@assay.used <- SeuratObject::DefaultAssay(srt_integrated)
  srt_integrated@graphs[["BBKNN"]] <- bbknn_graph

  bbknn_dist <- Matrix::t(
    SeuratObject::as.sparse(
      bem[[1]][1:nrow(bem[[1]]), , drop = FALSE]
    )
  )
  rownames(bbknn_dist) <- colnames(bbknn_dist) <- rownames(emb)
  bbknn_dist <- SeuratObject::as.Graph(bbknn_dist)
  bbknn_dist@assay.used <- SeuratObject::DefaultAssay(srt_integrated)
  srt_integrated@graphs[["BBKNN_dist"]] <- bbknn_dist

  val <- split(
    bbknn_dist@x,
    rep(seq_len(ncol(bbknn_dist)), diff(bbknn_dist@p))
  )
  pos <- split(
    bbknn_dist@i + 1,
    rep(seq_len(ncol(bbknn_dist)), diff(bbknn_dist@p))
  )
  idx <- Matrix::t(
    mapply(
      function(x, y) {
        out <- y[utils::head(order(x, decreasing = FALSE), n.neighbors - 1)]
        length(out) <- n.neighbors - 1
        out
      },
      x = val,
      y = pos
    )
  )
  idx[is.na(idx)] <- sample(
    seq_len(nrow(idx)),
    size = sum(is.na(idx)),
    replace = TRUE
  )
  idx <- cbind(seq_len(nrow(idx)), idx)
  dist <- Matrix::t(mapply(
    function(x, y) {
      out <- y[utils::head(order(x, decreasing = FALSE), n.neighbors - 1)]
      length(out) <- n.neighbors - 1
      out[is.na(out)] <- 0
      out
    },
    x = val,
    y = val
  ))
  dist <- cbind(0, dist)
  srt_integrated[["BBKNN_neighbors"]] <- methods::new(
    Class = "Neighbor",
    nn.idx = idx,
    nn.dist = dist,
    alg.info = list(),
    cell.names = rownames(emb)
  )
  nonlinear_reduction_params[["n.neighbors"]] <- n.neighbors

  srt_integrated <- tryCatch(
    {
      log_message(
        "Perform Seurat::FindClusters ({.val {cluster_algorithm}})"
      )
      srt_integrated <- Seurat::FindClusters(
        object = srt_integrated,
        graph.name = "BBKNN",
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        leiden_method = "igraph",
        verbose = FALSE
      )
      log_message("Reorder clusters...")
      srt_integrated <- srt_reorder(
        srt_integrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srt_integrated[["seurat_clusters"]] <- NULL
      srt_integrated[["BBKNNclusters"]] <- SeuratObject::Idents(srt_integrated)
      srt_integrated
    },
    error = function(error) {
      log_message(error, message_type = "warning")
      log_message(
        "Error when performing {.fn Seurat::FindClusters}. Skip this step",
        message_type = "warning"
      )
      srt_integrated
    }
  )

  srt_integrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        log_message(
          "Perform nonlinear dimension reduction ({.val {nr}})"
        )
        if (nr %in% c("fr")) {
          nonlinear_reduction_params[["n.neighbors"]] <- NULL
        } else {
          nonlinear_reduction_params[["n.neighbors"]] <- n.neighbors
        }
        for (n in nonlinear_reduction_dims) {
          srt_integrated <- RunDimReduction(
            srt_integrated,
            prefix = "BBKNN",
            neighbor_use = "BBKNN_neighbors",
            graph_use = "BBKNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srt_integrated
    },
    error = function(error) {
      log_message(error, message_type = "warning")
      log_message(
        "Error when performing nonlinear dimension reduction. Skip this step",
        message_type = "warning"
      )
      srt_integrated
    }
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[["BBKNN_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_integrated,
      pattern = paste0(assay, "|BBKNN|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_integrated)
  }
}

#' @title The CSS integration function
#'
#' @inheritParams integration_scop
#' @param CSS_dims_use A vector specifying the dimensions returned by CSS that will be utilized for downstream cell cluster finding and non-linear reduction.
#' If set to NULL, all the returned dimensions will be used by default.
#' @param CSS_params A list of parameters for the [simspec::cluster_sim_spectrum](https://github.com/quadbio/simspec) function.
#' Default is `list()`.
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
    seed = 11) {
  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "{.arg srt_list} and {.arg srt_merge} were all empty",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "{.arg srt_list} and {.arg srt_merge} have different cells",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
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
      "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
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

  check_r(c("quadbiolab/simspec", "qlcMatrix"))
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
  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))) {
    log_message("Perform ScaleData")
    srt_merge <- Seurat::ScaleData(
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
    "Perform linear dimension reduction({.val {linear_reduction}})"
  )
  srt_merge <- RunDimReduction(
    srt_merge,
    prefix = "CSS",
    features = HVF,
    assay = SeuratObject::DefaultAssay(srt_merge),
    linear_reduction = linear_reduction,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = force_linear_reduction,
    verbose = FALSE,
    seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srt_merge@reductions[[paste0(
      "CSS",
      linear_reduction
    )]]@misc[["dims_estimate"]]
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  log_message("Perform {.pkg CSS} integration")
  log_message(
    "Using Reduction({.val {paste0('CSS', linear_reduction)}}, dims:{.val {min(linear_reduction_dims_use)}}-{.val {max(linear_reduction_dims_use)}}) as input"
  )
  params <- list(
    object = SeuratObject::JoinLayers(srt_merge),
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

  srt_integrated <- tryCatch(
    {
      srt_integrated <- Seurat::FindNeighbors(
        object = srt_integrated,
        reduction = "CSS",
        dims = CSS_dims_use,
        annoy.metric = neighbor_metric,
        k.param = neighbor_k,
        graph.name = paste0("CSS", "_", c("KNN", "SNN")),
        verbose = FALSE
      )

      log_message(
        "Perform Seurat::FindClusters ({.val {cluster_algorithm}})"
      )
      srt_integrated <- Seurat::FindClusters(
        object = srt_integrated,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        leiden_method = "igraph",
        graph.name = "CSS_SNN",
        verbose = FALSE
      )
      log_message("Reorder clusters...")
      srt_integrated <- srt_reorder(
        srt_integrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srt_integrated[["seurat_clusters"]] <- NULL
      srt_integrated[["CSSclusters"]] <- SeuratObject::Idents(srt_integrated)
      srt_integrated
    },
    error = function(error) {
      log_message(error, message_type = "warning")
      log_message(
        "Error when performing {.fn Seurat::FindClusters}. Skip this step",
        message_type = "warning"
      )
      srt_integrated
    }
  )

  srt_integrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        log_message(
          paste0("Perform nonlinear dimension reduction (", nr, ")")
        )
        for (n in nonlinear_reduction_dims) {
          srt_integrated <- RunDimReduction(
            srt_integrated,
            prefix = "CSS",
            reduction_use = "CSS",
            reduction_dims = CSS_dims_use,
            graph_use = "CSS_SNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srt_integrated
    },
    error = function(error) {
      log_message(error, message_type = "warning")
      log_message(
        "Error when performing nonlinear dimension reduction. Skip this step",
        message_type = "warning"
      )
      srt_integrated
    }
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[["CSS_HVF"]] <- HVF

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

#' @title The LIGER integration function
#'
#' @md
#' @inheritParams integration_scop
#' @param LIGER_dims_use A vector specifying the dimensions returned by LIGER that will be utilized for downstream cell cluster finding and non-linear reduction.
#' If set to NULL, all the returned dimensions will be used by default.
#' @param optimizeALS_params A list of parameters for the [rliger::optimizeALS] function.
#' Default is `list()`.
#' @param quantilenorm_params A list of parameters for the [rliger::quantile_norm] function.
#' Default is `list()`.
#'
#' @export
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
    LIGER_dims_use = NULL,
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
    seed = 11) {
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
      "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
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

  check_r("rliger")
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "srt_list and srt_merge were all empty",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "srt_list and srt_merge have different cells",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
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
      message_type = "warning"
    )
    answer <- utils::askYesNo("Are you sure to continue?", default = FALSE)
    if (isFALSE(answer)) {
      return(srt_merge)
    }
  }

  scale.data <- list()
  for (i in seq_along(srt_list)) {
    srt <- srt_list[[i]]
    scale_features <- rownames(
      GetAssayData5(
        srt,
        layer = "scale.data",
        assay = SeuratObject::DefaultAssay(srt)
      )
    )
    if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))) {
      log_message("Perform ScaleData on the data {.val {i}} ...")
      srt <- Seurat::ScaleData(
        object = srt,
        assay = SeuratObject::DefaultAssay(srt),
        features = HVF,
        do.center = FALSE,
        vars.to.regress = vars_to_regress,
        model.use = regression_model,
        verbose = FALSE
      )
    }
    scale.data[[i]] <- Matrix::t(
      x = GetAssayData5(
        object = srt,
        layer = "scale.data",
        assay = SeuratObject::DefaultAssay(srt)
      )
    )
  }

  log_message("Perform {.pkg LIGER} integration")
  params1 <- list(
    object = scale.data,
    k = 20,
    verbose = FALSE
  )
  for (nm in names(optimizeALS_params)) {
    params1[[nm]] <- optimizeALS_params[[nm]]
  }
  out1 <- invoke_fun(rliger::optimizeALS, params1)
  colnames(x = out1$W) <- colnames(scale.data[[1]])
  reduction1 <- invoke_fun(what = "rbind", args = out1$H)
  colnames(reduction1) <- paste0("riNMF_", seq_len(ncol(reduction1)))
  loadings1 <- Matrix::t(x = out1$W)
  rownames(loadings1) <- colnames(scale.data[[1]])
  colnames(loadings1) <- paste0("riNMF_", seq_len(ncol(loadings1)))
  srt_merge[["iNMF_raw"]] <- CreateDimReducObject(
    embeddings = reduction1,
    loadings = loadings1,
    assay = SeuratObject::DefaultAssay(srt_merge),
    key = "riNMF_"
  )

  embeddings <- sapply(
    X = Seurat::SplitObject(object = srt_merge, split.by = batch),
    FUN = function(x) {
      return(Embeddings(object = x[["iNMF_raw"]]))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  num.samples <- vapply(
    X = embeddings,
    FUN = nrow,
    FUN.VALUE = integer(length = 1L)
  )
  ref_dataset <- names(x = embeddings)[which.max(x = num.samples)]
  params2 <- list(
    object = embeddings,
    ref_dataset = ref_dataset
  )
  for (nm in names(quantilenorm_params)) {
    params2[[nm]] <- quantilenorm_params[[nm]]
  }
  out2 <- invoke_fun(rliger::quantile_norm, params2)
  srt_merge[["LIGER"]] <- CreateDimReducObject(
    embeddings = out2$H.norm,
    assay = SeuratObject::DefaultAssay(srt_merge),
    key = "LIGER_"
  )
  srt_integrated <- srt_merge
  srt_merge <- NULL
  if (is.null(LIGER_dims_use)) {
    LIGER_dims_use <- seq_len(
      ncol(srt_integrated[["LIGER"]]@cell.embeddings)
    )
  }

  srt_integrated <- tryCatch(
    {
      srt_integrated <- Seurat::FindNeighbors(
        object = srt_integrated,
        reduction = "LIGER",
        dims = LIGER_dims_use,
        annoy.metric = neighbor_metric,
        k.param = neighbor_k,
        graph.name = paste0("LIGER", "_", c("KNN", "SNN")),
        verbose = FALSE
      )

      log_message(
        paste0("Perform Seurat::FindClusters (", cluster_algorithm, ")")
      )
      srt_integrated <- Seurat::FindClusters(
        object = srt_integrated,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        leiden_method = "igraph",
        graph.name = "LIGER_SNN",
        verbose = FALSE
      )
      log_message("Reorder clusters...")
      srt_integrated <- srt_reorder(
        srt_integrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srt_integrated[["seurat_clusters"]] <- NULL
      srt_integrated[["LIGERclusters"]] <- SeuratObject::Idents(srt_integrated)
      srt_integrated
    },
    error = function(error) {
      log_message(error, message_type = "warning", verbose = verbose)
      log_message(
        "Error when performing {.fn Seurat::FindClusters}. Skip this step",
        message_type = "warning",
        verbose = verbose
      )
      srt_integrated
    }
  )

  srt_integrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        log_message(
          paste0("Perform nonlinear dimension reduction (", nr, ")")
        )
        for (n in nonlinear_reduction_dims) {
          srt_integrated <- RunDimReduction(
            srt_integrated,
            prefix = "LIGER",
            reduction_use = "LIGER",
            reduction_dims = LIGER_dims_use,
            graph_use = "LIGER_SNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srt_integrated
    },
    error = function(error) {
      log_message(error, message_type = "warning", verbose = verbose)
      log_message(
        "Error when performing nonlinear dimension reduction. Skip this step",
        message_type = "warning",
        verbose = verbose
      )
      srt_integrated
    }
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[["LIGER_HVF"]] <- HVF

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

#' @title The Conos integration function
#'
#' @inheritParams integration_scop
#' @param buildGraph_params A list of parameters for the buildGraph function.
#' Default is `list()`.
#' @param num_threads  An integer setting the number of threads for Conos.
#' Default is `2`.
#'
#' @export
Conos_integrate <- function(
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
    linear_reduction = "pca",
    linear_reduction_dims = 50,
    linear_reduction_dims_use = NULL,
    linear_reduction_params = list(),
    force_linear_reduction = FALSE,
    nonlinear_reduction = "umap",
    nonlinear_reduction_dims = c(2, 3),
    nonlinear_reduction_params = list(),
    force_nonlinear_reduction = TRUE,
    cluster_algorithm = "louvain",
    cluster_resolution = 0.6,
    buildGraph_params = list(),
    num_threads = 2,
    verbose = TRUE,
    seed = 11) {
  if (length(linear_reduction) > 1) {
    log_message(
      "Only the first method in the {.arg linear_reduction} will be used",
      message_type = "warning",
      verbose = verbose
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
  nonlinear_reductions <- c("umap", "umap-naive", "fr")
  if (any(!nonlinear_reduction %in% nonlinear_reductions)) {
    log_message(
      "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
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

  check_r("conos")
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "srt_list and srt_merge were all empty",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "srt_list and srt_merge have different cells",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
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
      verbose = verbose,
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
      "The cell count in some batches is lower than {.val 30}, which may not be suitable for the {.pkg Conos} integration method",
      message_type = "warning"
    )
    answer <- utils::askYesNo(
      "Are you sure to continue?",
      default = FALSE
    )
    if (isFALSE(answer)) {
      return(srt_merge)
    }
  }

  srt_integrated <- srt_merge
  srt_merge <- NULL

  if (normalization_method == "TFIDF") {
    log_message(
      "{.arg normalization_method} is {.val TFIDF}. Use {.pkg lsi} workflow..."
    )
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  for (i in seq_along(srt_list)) {
    srt <- srt_list[[i]]
    scale_features <- rownames(
      GetAssayData5(
        srt,
        layer = "scale.data",
        assay = SeuratObject::DefaultAssay(srt)
      )
    )
    if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))) {
      log_message(
        "Perform ScaleData on the data {.val {i}} ..."
      )
      srt <- Seurat::ScaleData(
        object = srt,
        assay = SeuratObject::DefaultAssay(srt),
        features = HVF,
        vars.to.regress = vars_to_regress,
        model.use = regression_model,
        verbose = FALSE
      )
    }
    log_message(
      "Perform linear dimension reduction ({.val {linear_reduction}}) on the data {.val {i}} ..."
    )
    srt <- RunDimReduction(
      srt,
      prefix = "Conos",
      features = HVF,
      assay = SeuratObject::DefaultAssay(srt),
      linear_reduction = linear_reduction,
      linear_reduction_dims = linear_reduction_dims,
      linear_reduction_params = linear_reduction_params,
      force_linear_reduction = force_linear_reduction,
      verbose = FALSE,
      seed = seed
    )
    srt[["pca"]] <- srt[[paste0("Conos", linear_reduction)]]
    srt_list[[i]] <- srt
  }
  if (is.null(names(srt_list))) {
    names(srt_list) <- paste0("srt_", seq_along(srt_list))
  }

  if (is.null(linear_reduction_dims_use)) {
    maxdims <- max(unlist(sapply(
      srt_list,
      function(srt) {
        max(srt@reductions[[paste0("Conos", linear_reduction)]]@misc[[
          "dims_estimate"
        ]])
      }
    )))
  } else {
    maxdims <- max(linear_reduction_dims_use)
  }

  log_message(
    " Perform {.pkg Conos} integration"
  )
  log_message(
    "Conos integration using Reduction(",
    linear_reduction,
    ", dims_max:",
    maxdims,
    ") as input"
  )
  srt_list_con <- conos::Conos$new(srt_list, n.cores = num_threads)
  params <- list(
    ncomps = maxdims,
    verbose = FALSE
  )
  for (nm in names(buildGraph_params)) {
    params[[nm]] <- buildGraph_params[[nm]]
  }
  invoke_fun(srt_list_con[["buildGraph"]], params)
  conos_graph <- igraph::as_adjacency_matrix(
    srt_list_con$graph,
    type = "both",
    attr = "weight",
    names = TRUE,
    sparse = TRUE
  )
  conos_graph <- SeuratObject::as.Graph(conos_graph)
  conos_graph@assay.used <- SeuratObject::DefaultAssay(srt_integrated)
  srt_integrated@graphs[["Conos"]] <- conos_graph
  nonlinear_reduction_params[["n.neighbors"]] <- params[["k"]]

  srt_integrated <- tryCatch(
    {
      log_message(
        "Perform {.fn Seurat::FindClusters} ({.val {cluster_algorithm}})"
      )
      srt_integrated <- Seurat::FindClusters(
        object = srt_integrated,
        graph.name = "Conos",
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        leiden_method = "igraph",
        verbose = FALSE
      )
      log_message("Reorder clusters...")
      srt_integrated <- srt_reorder(
        srt_integrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srt_integrated[["seurat_clusters"]] <- NULL
      srt_integrated[["Conosclusters"]] <- SeuratObject::Idents(srt_integrated)
      srt_integrated
    },
    error = function(error) {
      log_message(
        error,
        message_type = "warning"
      )
      log_message(
        "Error when performing {.fn Seurat::FindClusters}. Skip this step",
        message_type = "warning"
      )
      srt_integrated
    }
  )

  srt_integrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        log_message(
          "Perform nonlinear dimension reduction ({.val {nr}})"
        )
        if (nr %in% c("fr")) {
          nonlinear_reduction_params[["n.neighbors"]] <- NULL
        } else {
          nonlinear_reduction_params[["n.neighbors"]] <- params[["k"]]
        }
        for (n in nonlinear_reduction_dims) {
          srt_integrated <- RunDimReduction(
            srt_integrated,
            prefix = "Conos",
            graph_use = "Conos",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srt_integrated
    },
    error = function(error) {
      log_message(
        error,
        message_type = "warning"
      )
      log_message(
        "Error when performing nonlinear dimension reduction. Skip this step",
        message_type = "warning"
      )
      srt_integrated
    }
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[["Conos_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_integrated,
      pattern = paste0(assay, "|Conos|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_integrated)
  }
}

#' @title The ComBat integration function
#'
#' @inheritParams integration_scop
#' @param ComBat_params A list of parameters for the sva::ComBat function.
#' Default is `list()`.
#'
#' @export
ComBat_integrate <- function(
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
    ComBat_params = list(),
    verbose = TRUE,
    seed = 11) {
  if (length(linear_reduction) > 1) {
    log_message(
      "Only the first method in the {.arg linear_reduction} will be used",
      message_type = "warning",
      verbose = verbose,
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
      "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
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

  check_r("sva")
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "srt_list and srt_merge were all empty",
      message_type = "error"
    )
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      log_message(
        "srt_list and srt_merge have different cells",
        message_type = "error"
      )
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
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
      verbose = verbose,
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

  if (normalization_method == "TFIDF") {
    log_message(
      "{.arg normalization_method} is {.pkg TFIDF}. Use {.pkg lsi} workflow..."
    )
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  log_message("Perform {.pkg ComBat} integration")
  dat <- GetAssayData5(
    srt_merge,
    layer = "data",
    assay = SeuratObject::DefaultAssay(srt_merge)
  )[HVF, , drop = FALSE]
  batch <- srt_merge[[batch, drop = TRUE]]
  params <- list(
    dat = dat,
    batch = batch
  )
  for (nm in names(ComBat_params)) {
    params[[nm]] <- ComBat_params[[nm]]
  }
  corrected <- suppressMessages(
    invoke_fun(
      sva::ComBat,
      params
    )
  )

  srt_integrated <- srt_merge
  srt_merge <- NULL
  srt_integrated[["ComBatcorrected"]] <- CreateAssayObject(data = corrected)
  SeuratObject::DefaultAssay(srt_integrated) <- "ComBatcorrected"
  SeuratObject::VariableFeatures(srt_integrated[["ComBatcorrected"]]) <- HVF
  scale_features <- rownames(
    GetAssayData5(
      srt_integrated,
      layer = "scale.data",
      assay = SeuratObject::DefaultAssay(srt_integrated)
    )
  )
  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))) {
    log_message("Perform {.fn Seurat::ScaleData}")
    srt_integrated <- Seurat::ScaleData(
      srt_integrated,
      split.by = if (isTRUE(scale_within_batch)) batch else NULL,
      assay = SeuratObject::DefaultAssay(srt_integrated),
      features = HVF,
      vars.to.regress = vars_to_regress,
      model.use = regression_model,
      verbose = FALSE
    )
  }

  log_message(
    "Perform linear dimension reduction ({.val {linear_reduction}})"
  )
  srt_integrated <- RunDimReduction(
    srt_integrated,
    prefix = "ComBat",
    features = HVF,
    assay = SeuratObject::DefaultAssay(srt_integrated),
    linear_reduction = linear_reduction,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = force_linear_reduction,
    verbose = FALSE,
    seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srt_integrated@reductions[[paste0(
      "ComBat",
      linear_reduction
    )]]@misc[["dims_estimate"]]
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  srt_integrated <- tryCatch(
    {
      srt_integrated <- Seurat::FindNeighbors(
        object = srt_integrated,
        reduction = paste0("ComBat", linear_reduction),
        dims = linear_reduction_dims_use,
        annoy.metric = neighbor_metric,
        k.param = neighbor_k,
        graph.name = paste0("ComBat_", c("KNN", "SNN")),
        verbose = FALSE
      )

      log_message(
        "Perform Seurat::FindClusters ({.val {cluster_algorithm}})"
      )
      srt_integrated <- Seurat::FindClusters(
        object = srt_integrated,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        leiden_method = "igraph",
        graph.name = "ComBat_SNN",
        verbose = FALSE
      )
      log_message("Reorder clusters...")
      srt_integrated <- srt_reorder(
        srt_integrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srt_integrated[["seurat_clusters"]] <- NULL
      srt_integrated[["ComBatclusters"]] <- SeuratObject::Idents(srt_integrated)
      srt_integrated
    },
    error = function(error) {
      log_message(
        error,
        message_type = "warning"
      )
      log_message(
        "Error when performing {.fn Seurat::FindClusters}. Skip this step",
        message_type = "warning"
      )
      srt_integrated
    }
  )

  srt_integrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        log_message(
          paste0("Perform nonlinear dimension reduction (", nr, ")")
        )
        for (n in nonlinear_reduction_dims) {
          srt_integrated <- RunDimReduction(
            srt_integrated,
            prefix = "ComBat",
            reduction_use = paste0("ComBat", linear_reduction),
            reduction_dims = linear_reduction_dims_use,
            graph_use = "ComBat_SNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srt_integrated
    },
    error = function(error) {
      log_message(
        error,
        message_type = "warning",
        verbose = verbose
      )
      log_message(
        "Error when performing nonlinear dimension reduction. Skip this step",
        message_type = "warning",
        verbose = verbose
      )
      return(srt_integrated)
    }
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[["ComBat_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_integrated,
      pattern = paste0(assay, "|ComBat|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_integrated)
  }
}
