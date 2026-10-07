#' @title Conos integration function
#'
#' @inheritParams RunIntegration
#' @param buildGraph_params A list of parameters for the buildGraph function.
#' @param cores An integer setting the number of threads for `Conos`.
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
  cores = 2,
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
  validate_nonlinear_reductions(
    nonlinear_reduction,
    allowed = c("umap", "umap-naive", "fr")
  )
  cluster_algorithm_index <- resolve_cluster_algorithm_index(cluster_algorithm)

  check_r("conos", verbose = FALSE)
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
    if (
      isTRUE(do_scaling) ||
        (is.null(do_scaling) && any(!HVF %in% scale_features))
    ) {
      log_message(
        "Perform ScaleData on the data {.val {i}} ..."
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
      "Perform {.val {linear_reduction}} linear dimension reduction"
    )
    srt <- RunDimsReduction(
      srt,
      prefix = "Conos",
      features = HVF,
      assay = SeuratObject::DefaultAssay(srt),
      linear_reduction = linear_reduction,
      linear_reduction_dims = linear_reduction_dims,
      linear_reduction_params = linear_reduction_params,
      force_linear_reduction = force_linear_reduction,
      verbose = verbose,
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
        max(RunDimsEstimate(
          object = srt,
          reduction = paste0("Conos", linear_reduction),
          reduction_method = linear_reduction,
          skip_first = normalization_method == "TFIDF",
          use_stored = TRUE,
          verbose = FALSE
        ))
      }
    )))
  } else {
    maxdims <- max(linear_reduction_dims_use)
  }

  log_message(
    " Perform {.pkg Conos} integration"
  )
  log_message(
    "{.pkg Conos} integration using {.pkg {linear_reduction}} ({.val {1}}:{.val {maxdims}}) as input",
    verbose = verbose
  )
  srt_list_con <- NULL
  conos_fun <- get_namespace_fun("conos", "Conos")$new
  invisible(
    utils::capture.output(
      srt_list_con <- suppressWarnings(
        suppressMessages(
          conos_fun(
            srt_list,
            n.cores = cores,
            verbose = FALSE
          )
        )
      ),
      type = "output"
    )
  )
  params <- list(
    ncomps = maxdims,
    verbose = FALSE
  )
  for (nm in names(buildGraph_params)) {
    params[[nm]] <- buildGraph_params[[nm]]
  }
  invisible(
    utils::capture.output(
      suppressWarnings(
        suppressMessages(
          invoke_fun(srt_list_con[["buildGraph"]], params)
        )
      ),
      type = "output"
    )
  )
  conos_graph <- igraph::as_adjacency_matrix(
    srt_list_con$graph,
    type = "both",
    attr = "weight",
    names = TRUE,
    sparse = TRUE
  )
  graph_cells <- colnames(conos_graph)
  object_cells <- colnames(srt_integrated)
  if (!setequal(graph_cells, object_cells)) {
    log_message(
      "Cell names in {.pkg Conos} graph do not match {.arg srt_integrated}",
      message_type = "error"
    )
  }
  conos_graph <- conos_graph[object_cells, object_cells, drop = FALSE]
  conos_graph <- SeuratObject::as.Graph(conos_graph)
  conos_graph@assay.used <- SeuratObject::DefaultAssay(srt_integrated)
  srt_integrated@graphs[["Conos"]] <- conos_graph
  nonlinear_reduction_params[["n.neighbors"]] <- params[["k"]]

  srt_integrated <- find_neighbors_and_clusters(
    srt = srt_integrated,
    reduction = NULL,
    dims_use = NULL,
    graph_prefix = "Conos_",
    graph_snn = "Conos",
    cluster_colname = "Conosclusters",
    HVF = HVF,
    neighbor_metric = "euclidean",
    neighbor_k = 20L,
    cluster_algorithm = cluster_algorithm,
    cluster_algorithm_index = cluster_algorithm_index,
    cluster_resolution = cluster_resolution,
    run_find_neighbors = FALSE,
    verbose = verbose
  )

  srt_integrated <- run_nonlinear_reduction(
    srt = srt_integrated,
    prefix = "Conos",
    graph_use = "Conos",
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed,
    verbose = verbose
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- srt_integrated@misc[[
    "Conos_HVF"
  ]] <- HVF

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
