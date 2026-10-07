#' @title Seurat integration function
#'
#' @inheritParams RunIntegration
#' @param FindIntegrationAnchors_params A list of parameters for the Seurat::FindIntegrationAnchors function.
#' @param IntegrateData_params A list of parameters for the Seurat::IntegrateData function.
#' @param IntegrateEmbeddings_params A list of parameters for the Seurat::IntegrateEmbeddings function.
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
  seed = 11
) {
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
  validate_nonlinear_reductions(nonlinear_reduction)
  cluster_algorithm_index <- resolve_cluster_algorithm_index(cluster_algorithm)

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

  reduction_key <- FindIntegrationAnchors_params[["reduction"]]
  if (!is.null(reduction_key) && normalization_method != "TFIDF") {
    reduction_key <- tolower(
      gsub("[^a-z]", "", as.character(reduction_key)[1L], perl = TRUE)
    )
    if (reduction_key %in% c("cca", "ccaintegration")) {
      FindIntegrationAnchors_params[["reduction"]] <- "cca"
    } else if (reduction_key %in% c("rpca", "rpcaintegration")) {
      FindIntegrationAnchors_params[["reduction"]] <- "rpca"
    }
  }

  if (min(sapply(srt_list, ncol)) < 50) {
    log_message(
      "The cell count in some batches is lower than 50, which may not be suitable for the current integration method",
      message_type = "warning"
    )
    if (interactive()) {
      answer <- utils::askYesNo("Are you sure to continue?", default = FALSE)
      if (isFALSE(answer)) {
        return(srt_merge)
      }
    } else {
      log_message(
        "Non-interactive session detected. Continue integration after warning",
        message_type = "warning"
      )
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
      max_anchor_dim <- min(
        linear_reduction_dims,
        30,
        min(sapply(srt_list, ncol)) - 1L
      )
      FindIntegrationAnchors_params[["dims"]] <- if (max_anchor_dim >= 2) {
        2:max_anchor_dim
      } else {
        1L
      }
    }
    srt_merge <- Signac::RunTFIDF(
      object = srt_merge,
      assay = SeuratObject::DefaultAssay(srt_merge),
      verbose = FALSE
    )
    srt_merge <- RunDimsReduction(
      srt_merge,
      prefix = "",
      features = HVF,
      assay = SeuratObject::DefaultAssay(srt_merge),
      linear_reduction = "svd",
      linear_reduction_dims = linear_reduction_dims,
      linear_reduction_params = linear_reduction_params,
      force_linear_reduction = force_linear_reduction,
      verbose = verbose,
      seed = seed
    )
    srt_merge[["lsi"]] <- srt_merge[["svd"]]
    for (i in seq_along(srt_list)) {
      srt <- srt_list[[i]]
      log_message(
        "Perform {.pkg svd} linear dimension reduction on {.val {i}} of {.arg srt_list}"
      )
      srt <- RunDimsReduction(
        srt,
        prefix = "",
        features = HVF,
        assay = SeuratObject::DefaultAssay(srt),
        linear_reduction = "svd",
        linear_reduction_dims = linear_reduction_dims,
        linear_reduction_params = linear_reduction_params,
        force_linear_reduction = force_linear_reduction,
        verbose = verbose,
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
          "Perform {.fn Seurat::ScaleData} on {.arg srt}"
        )
        srt <- ScaleData(
          object = srt,
          assay = SeuratObject::DefaultAssay(srt),
          features = HVF,
          vars.to.regress = vars_to_regress,
          model.use = regression_model,
          verbose = FALSE
        )
      }
      log_message(
        "Perform {.pkg pca} linear dimension reduction on {.val {i}} of {.arg srt_list}"
      )
      srt <- RunDimsReduction(
        srt,
        prefix = "",
        features = HVF,
        assay = SeuratObject::DefaultAssay(srt),
        linear_reduction = "pca",
        linear_reduction_dims = linear_reduction_dims,
        linear_reduction_params = linear_reduction_params,
        force_linear_reduction = force_linear_reduction,
        verbose = verbose,
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

    log_message("Perform {.pkg Seurat} integration")
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
      log_message("Perform ScaleData on {.arg srt_integrated}")
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
      "Perform {.val {linear_reduction}} linear dimension reduction"
    )
    srt_integrated <- RunDimsReduction(
      srt_integrated,
      prefix = "Seurat",
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
      linear_reduction_dims_use <- srt_integrated@reductions[[paste0(
        "Seurat",
        linear_reduction
      )]]@misc[["dims_estimate"]] %||%
        1:linear_reduction_dims
    }
  } else if (normalization_method == "TFIDF") {
    log_message(
      "Perform {.fn FindIntegrationAnchors} with {.arg reduction = rlsi}"
    )
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

    log_message("Perform {.pkg Seurat} integration")
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

  srt_integrated <- find_neighbors_and_clusters(
    srt = srt_integrated,
    reduction = paste0("Seurat", linear_reduction),
    dims_use = linear_reduction_dims_use,
    graph_prefix = "Seurat_",
    graph_snn = "Seurat_SNN",
    cluster_colname = "Seuratclusters",
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
    prefix = "Seurat",
    reduction_use = paste0("Seurat", linear_reduction),
    reduction_dims = linear_reduction_dims_use,
    graph_use = "Seurat_SNN",
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed,
    verbose = verbose
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
