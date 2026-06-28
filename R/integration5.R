#' @title Seurat v5 CCA integration
#'
#' @inheritParams integration_scop
#' @param IntegrateLayers_params A list of parameters passed to [Seurat::IntegrateLayers].
#'
#' @export
CCA_integrate <- function(
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
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  IntegrateLayers_params = list(),
  verbose = TRUE,
  seed = 11
) {
  run_integration5(
    method_fun = get_namespace_fun("Seurat", "CCAIntegration"),
    method_label = "CCAIntegration",
    reduction_name = "CCA",
    graph_prefix = "CCA_",
    graph_snn = "CCA_SNN",
    cluster_colname = "CCAclusters",
    misc_hvf_key = "CCA_HVF",
    append_tag = "CCA",
    srt_merge = srt_merge,
    batch = batch,
    append = append,
    srt_list = srt_list,
    assay = assay,
    do_normalization = do_normalization,
    normalization_method = normalization_method,
    do_HVF_finding = do_HVF_finding,
    HVF_source = HVF_source,
    HVF_method = HVF_method,
    nHVF = nHVF,
    HVF_min_intersection = HVF_min_intersection,
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
    IntegrateLayers_params = IntegrateLayers_params,
    verbose = verbose,
    seed = seed
  )
}

#' @title Seurat v5 RPCA integration
#'
#' @inheritParams integration_scop
#' @inheritParams CCA_integrate
#'
#' @export
RPCA_integrate <- function(
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
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  IntegrateLayers_params = list(),
  verbose = TRUE,
  seed = 11
) {
  run_integration5(
    method_fun = get_namespace_fun("Seurat", "RPCAIntegration"),
    method_label = "RPCAIntegration",
    reduction_name = "RPCA",
    graph_prefix = "RPCA_",
    graph_snn = "RPCA_SNN",
    cluster_colname = "RPCAclusters",
    misc_hvf_key = "RPCA_HVF",
    append_tag = "RPCA",
    srt_merge = srt_merge,
    batch = batch,
    append = append,
    srt_list = srt_list,
    assay = assay,
    do_normalization = do_normalization,
    normalization_method = normalization_method,
    do_HVF_finding = do_HVF_finding,
    HVF_source = HVF_source,
    HVF_method = HVF_method,
    nHVF = nHVF,
    HVF_min_intersection = HVF_min_intersection,
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
    IntegrateLayers_params = IntegrateLayers_params,
    verbose = verbose,
    seed = seed
  )
}

#' @title Seurat v5 fastMNN integration
#'
#' @inheritParams integration_scop
#' @inheritParams CCA_integrate
#' @param fastMNN_dims_use A vector specifying the integrated dimensions used for downstream clustering and nonlinear reduction.
#'
#' @export
fastMNN5_integrate <- function(
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
  IntegrateLayers_params = list(),
  verbose = TRUE,
  seed = 11
) {
  check_r("satijalab/seurat-wrappers", verbose = FALSE)
  run_integration5(
    method_fun = get_namespace_fun("SeuratWrappers", "FastMNNIntegration"),
    method_label = "FastMNNIntegration",
    reduction_name = "fastMNN5",
    graph_prefix = "fastMNN5_",
    graph_snn = "fastMNN5_SNN",
    cluster_colname = "fastMNN5clusters",
    misc_hvf_key = "fastMNN5_HVF",
    append_tag = "fastMNN5",
    srt_merge = srt_merge,
    batch = batch,
    append = append,
    srt_list = srt_list,
    assay = assay,
    do_normalization = do_normalization,
    normalization_method = normalization_method,
    do_HVF_finding = do_HVF_finding,
    HVF_source = HVF_source,
    HVF_method = HVF_method,
    nHVF = nHVF,
    HVF_min_intersection = HVF_min_intersection,
    HVF = HVF,
    do_scaling = TRUE,
    vars_to_regress = NULL,
    regression_model = "linear",
    linear_reduction = "pca",
    linear_reduction_dims = 50L,
    linear_reduction_dims_use = NULL,
    linear_reduction_params = list(),
    force_linear_reduction = TRUE,
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    neighbor_metric = neighbor_metric,
    neighbor_k = neighbor_k,
    cluster_algorithm = cluster_algorithm,
    cluster_resolution = cluster_resolution,
    IntegrateLayers_params = IntegrateLayers_params,
    dims_use_override = fastMNN_dims_use,
    verbose = verbose,
    seed = seed
  )
}

#' @title Seurat v5 Harmony integration
#'
#' @inheritParams integration_scop
#' @inheritParams CCA_integrate
#' @param harmony_dims_use A vector specifying the integrated dimensions used for downstream clustering and nonlinear reduction.
#'
#' @export
Harmony5_integrate <- function(
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
  harmony_dims_use = NULL,
  nonlinear_reduction = "umap",
  nonlinear_reduction_dims = c(2, 3),
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  IntegrateLayers_params = list(),
  verbose = TRUE,
  seed = 11
) {
  run_integration5(
    method_fun = get_namespace_fun("Seurat", "HarmonyIntegration"),
    method_label = "HarmonyIntegration",
    reduction_name = "Harmony5",
    graph_prefix = "Harmony5_",
    graph_snn = "Harmony5_SNN",
    cluster_colname = "Harmony5clusters",
    misc_hvf_key = "Harmony5_HVF",
    append_tag = "Harmony5",
    srt_merge = srt_merge,
    batch = batch,
    append = append,
    srt_list = srt_list,
    assay = assay,
    do_normalization = do_normalization,
    normalization_method = normalization_method,
    do_HVF_finding = do_HVF_finding,
    HVF_source = HVF_source,
    HVF_method = HVF_method,
    nHVF = nHVF,
    HVF_min_intersection = HVF_min_intersection,
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
    IntegrateLayers_params = IntegrateLayers_params,
    dims_use_override = harmony_dims_use,
    verbose = verbose,
    seed = seed
  )
}

#' @title Seurat v5 scVI integration
#'
#' @md
#' @inheritParams integration_scop
#' @inheritParams CCA_integrate
#' @param scVI_dims_use A vector specifying the integrated dimensions used for downstream clustering and nonlinear reduction.
#' @param cores Number of DataLoader worker processes for `scVI` training.
#' `NULL` (default) uses the PyTorch default (0, single-process).
#' Increase to speed up data loading on multi-core machines.
#'
#' @export
scVI5_integrate <- function(
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
  cores = NULL,
  IntegrateLayers_params = list(),
  verbose = TRUE,
  seed = 11
) {
  PrepareEnv(modules = "scvi")
  check_python("scvi-tools")
  check_r("satijalab/seurat-wrappers", verbose = FALSE)

  scvi_method_fun <- function(
    object,
    features = NULL,
    layers = "counts",
    conda_env = NULL,
    new.reduction = "integrated.dr",
    ndims = 30,
    nlayers = 2,
    gene_likelihood = "nb",
    max_epochs = NULL,
    ...
  ) {
    reticulate::use_condaenv(conda_env, required = TRUE)
    sc <- reticulate::import("scanpy", convert = FALSE)
    scvi <- reticulate::import("scvi", convert = FALSE)
    scipy <- reticulate::import("scipy", convert = FALSE)
    if (is.null(max_epochs)) {
      max_epochs <- reticulate::r_to_py(max_epochs)
    } else {
      max_epochs <- as.integer(max_epochs)
    }
    if (inherits(object, what = "SCTAssay")) {
      batches <- get_namespace_fun("SeuratWrappers", ".FindSCTBatches")(object)
    } else {
      batches <- get_namespace_fun("SeuratWrappers", ".FindBatches")(
        object,
        layers = layers
      )
      object <- SeuratObject::JoinLayers(object = object, layers = "counts")
    }
    adata <- sc$AnnData(
      X = scipy$sparse$csr_matrix(
        Matrix::t(SeuratObject::LayerData(object, layer = "counts")[features, ])
      ),
      obs = batches,
      var = object[[]][features, ]
    )
    scvi$model$SCVI$setup_anndata(adata, batch_key = "batch")
    model <- scvi$model$SCVI(
      adata = adata,
      n_latent = as.integer(ndims),
      n_layers = as.integer(nlayers),
      gene_likelihood = gene_likelihood
    )
    train_args <- list(max_epochs = max_epochs)
    if (!is.null(cores)) {
      train_args[["datasplitter_kwargs"]] <- list(
        num_workers = as.integer(cores),
        persistent_workers = TRUE
      )
    }
    py_warnings <- reticulate::import("warnings")
    py_warnings$filterwarnings("ignore")
    old_pw <- Sys.getenv("PYTHONWARNINGS", unset = NA_character_)
    Sys.setenv(PYTHONWARNINGS = "ignore")
    tryCatch(
      do.call(model$train, train_args),
      finally = {
        py_warnings$resetwarnings()
        if (is.na(old_pw)) {
          Sys.unsetenv("PYTHONWARNINGS")
        } else {
          Sys.setenv(PYTHONWARNINGS = old_pw)
        }
      }
    )
    latent <- model$get_latent_representation()
    latent <- as.matrix(latent)
    rownames(latent) <- reticulate::py_to_r(adata$obs$index$values)
    colnames(latent) <- paste0(new.reduction, "_", seq_len(ncol(latent)))
    suppressWarnings(
      latent.dr <- Seurat::CreateDimReducObject(
        embeddings = latent,
        key = new.reduction
      )
    )
    output.list <- list(latent.dr)
    names(output.list) <- new.reduction
    return(output.list)
  }

  ilp <- IntegrateLayers_params
  if (is.null(ilp[["conda_env"]])) {
    ilp[["conda_env"]] <- get_envname()
  }
  run_integration5(
    method_fun = scvi_method_fun,
    method_label = "scVIIntegration",
    reduction_name = "scVI5",
    graph_prefix = "scVI5_",
    graph_snn = "scVI5_SNN",
    cluster_colname = "scVI5clusters",
    misc_hvf_key = "scVI5_HVF",
    append_tag = "scVI5",
    srt_merge = srt_merge,
    batch = batch,
    append = append,
    srt_list = srt_list,
    assay = assay,
    do_normalization = do_normalization,
    normalization_method = normalization_method,
    do_HVF_finding = do_HVF_finding,
    HVF_source = HVF_source,
    HVF_method = HVF_method,
    nHVF = nHVF,
    HVF_min_intersection = HVF_min_intersection,
    HVF = HVF,
    do_scaling = TRUE,
    vars_to_regress = NULL,
    regression_model = "linear",
    linear_reduction = "pca",
    linear_reduction_dims = 50L,
    linear_reduction_dims_use = NULL,
    linear_reduction_params = list(),
    force_linear_reduction = TRUE,
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    neighbor_metric = neighbor_metric,
    neighbor_k = neighbor_k,
    cluster_algorithm = cluster_algorithm,
    cluster_resolution = cluster_resolution,
    IntegrateLayers_params = ilp,
    dims_use_override = scVI_dims_use,
    verbose = verbose,
    seed = seed
  )
}

run_integration5 <- function(
  method_fun,
  method_label,
  reduction_name,
  graph_prefix,
  graph_snn,
  cluster_colname,
  misc_hvf_key,
  append_tag,
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
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  IntegrateLayers_params = list(),
  dims_use_override = NULL,
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
  if (!identical(tolower(linear_reduction), "pca")) {
    log_message(
      "{.pkg {method_label}} require {.arg linear_reduction = 'pca'}",
      message_type = "error"
    )
  }
  if (!normalization_method %in% c("LogNormalize", "SCT")) {
    log_message(
      "{.pkg {method_label}} require {.arg normalization_method} to be one of {.val {c('LogNormalize', 'SCT')}}",
      message_type = "error"
    )
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

  srt_merge_raw <- srt_merge
  if (!is.null(srt_list)) {
    if (
      !inherits(srt_list, "list") ||
        any(sapply(srt_list, function(x) !inherits(x, "Seurat")))
    ) {
      log_message(
        "{.arg srt_list} is not a list of {.cls Seurat}",
        message_type = "error"
      )
    }
    which_less_2 <- which(sapply(srt_list, ncol) < 2)
    if (length(which_less_2) > 0) {
      log_message(
        "{.cls Seurat} in {.arg srt_list} contain less than 2 cells. Index: {.val {which_less_2}}",
        message_type = "error"
      )
    }
    default_assays <- unique(sapply(srt_list, SeuratObject::DefaultAssay))
    if (is.null(assay)) {
      if (length(default_assays) != 1) {
        log_message(
          "Default assay of {.cls Seurat} objects in {.arg srt_list} is inconsistent",
          message_type = "error"
        )
      }
      assay <- default_assays[[1]]
    }
    if (length(batch) == length(srt_list)) {
      srt_list_tmp <- list()
      for (bat in unique(batch)) {
        srt_list_tmp[[bat]] <- Reduce(merge, srt_list[batch == bat])
      }
      srt_list <- srt_list_tmp
      batch_colname <- ".seurat5_batch"
      for (i in seq_along(srt_list)) {
        srt_list[[i]][[batch_colname]] <- names(srt_list)[i]
      }
      batch <- batch_colname
    }
    log_message(
      "Merge {.arg srt_list} into a single {.cls Seurat} object",
      verbose = verbose
    )
    srt_merge <- Reduce(merge, srt_list)
  } else {
    assay <- assay %||% SeuratObject::DefaultAssay(srt_merge)
  }

  if (length(batch) != 1 || !batch %in% colnames(srt_merge@meta.data)) {
    log_message(
      "{.arg batch} must be a single column name present in meta.data",
      message_type = "error"
    )
  }
  if (!is.factor(srt_merge[[batch, drop = TRUE]])) {
    srt_merge[[batch]] <- factor(
      srt_merge[[batch, drop = TRUE]],
      levels = unique(srt_merge[[batch, drop = TRUE]])
    )
  }

  if (!inherits(Seurat::GetAssay(srt_merge, assay = assay), "Assay5")) {
    log_message(
      "{.pkg {method_label}} requires an {.cls Assay5} (RNA) assay",
      message_type = "error"
    )
  }

  srt_integrated <- srt_merge

  if (identical(normalization_method, "LogNormalize")) {
    SeuratObject::DefaultAssay(srt_integrated) <- assay
    assay_use <- assay

    srt_integrated[[assay]] <- SeuratObject::JoinLayers(srt_integrated[[assay]])
    if ("scale.data" %in% SeuratObject::Layers(srt_integrated[[assay]])) {
      srt_integrated[[assay]]$scale.data <- NULL
    }
    srt_integrated[[assay]] <- split(
      x = srt_integrated[[assay]],
      f = srt_integrated[[batch, drop = TRUE]]
    )

    data_layers <- SeuratObject::Layers(
      srt_integrated[[assay]],
      search = "data"
    )
    if (is.null(do_normalization)) {
      if (length(data_layers) == 0) {
        do_normalization <- TRUE
      } else {
        status <- suppressWarnings(
          CheckDataType(
            srt_integrated,
            layer = data_layers[[1]],
            assay = assay,
            verbose = FALSE
          )
        )
        if (identical(status, "log_normalized_counts")) {
          log_message(
            "Data appears already log-normalized. Skipping {.fn Seurat::NormalizeData}",
            verbose = verbose
          )
          do_normalization <- FALSE
        } else {
          log_message(
            "Data is {.val {status}}. Will perform {.fn Seurat::NormalizeData}",
            message_type = "warning",
            verbose = verbose
          )
          do_normalization <- TRUE
        }
      }
    }
    if (isTRUE(do_normalization)) {
      log_message(
        "Perform {.fn Seurat::NormalizeData} on split layers for Seurat v5 integration",
        verbose = verbose
      )
      srt_integrated <- Seurat::NormalizeData(
        object = srt_integrated,
        assay = assay,
        normalization.method = "LogNormalize",
        verbose = FALSE
      )
    }

    if (is.null(HVF)) {
      if (
        isTRUE(do_HVF_finding) ||
          is.null(do_HVF_finding) ||
          length(SeuratObject::VariableFeatures(
            srt_integrated,
            assay = assay
          )) ==
            0
      ) {
        if (identical(HVF_source, "separate")) {
          log_message(
            "Perform {.fn Seurat::FindVariableFeatures} per batch ({.arg HVF_source = 'separate'})",
            verbose = verbose
          )
          srt_integrated <- Seurat::FindVariableFeatures(
            object = srt_integrated,
            assay = assay,
            nfeatures = nHVF,
            selection.method = HVF_method,
            verbose = FALSE
          )
          HVF <- SeuratObject::VariableFeatures(
            srt_integrated,
            assay = assay,
            nfeatures = nHVF
          )
          if (HVF_min_intersection > 1) {
            hvf_layer_search <- if (identical(HVF_method, "vst")) {
              "counts"
            } else {
              "data"
            }
            per_layer_hvf <- lapply(
              SeuratObject::Layers(
                srt_integrated[[assay]],
                search = hvf_layer_search
              ),
              function(lyr) {
                SeuratObject::VariableFeatures(
                  srt_integrated[[assay]],
                  layer = lyr,
                  nfeatures = nHVF
                )
              }
            )
            hvf_freq <- table(unlist(per_layer_hvf))
            HVF <- HVF[
              HVF %in% names(hvf_freq[hvf_freq >= HVF_min_intersection])
            ]
          }
        } else {
          log_message(
            "Use global HVF from merged dataset ({.arg HVF_source = 'global'})",
            verbose = verbose
          )
          srt_global <- srt_integrated
          srt_global[[assay]] <- SeuratObject::JoinLayers(srt_global[[assay]])
          srt_global <- Seurat::FindVariableFeatures(
            object = srt_global,
            assay = assay,
            nfeatures = nHVF,
            selection.method = HVF_method,
            verbose = FALSE
          )
          HVF <- SeuratObject::VariableFeatures(srt_global, assay = assay)
          rm(srt_global)
        }
      } else {
        HVF <- SeuratObject::VariableFeatures(srt_integrated, assay = assay)
      }
    } else {
      cf <- rownames(Seurat::GetAssay(srt_integrated, assay = assay))
      HVF <- HVF[HVF %in% cf]
    }

    if (length(HVF) == 0) {
      log_message("No HVF available", message_type = "error")
    }
    log_message(
      "Number of available HVF: {.val {length(HVF)}}",
      verbose = verbose
    )

    hvf_counts <- tryCatch(
      Matrix::colSums(
        GetAssayData5(srt_integrated, layer = "counts", assay = assay)[
          HVF, ,
          drop = FALSE
        ]
      ),
      error = function(e) NULL
    )
    if (!is.null(hvf_counts)) {
      cell_abnormal <- names(hvf_counts)[hvf_counts == 0]
      if (length(cell_abnormal) > 0) {
        log_message(
          "Some cells do not express any HVF: {.val {cell_abnormal}}",
          message_type = "warning",
          verbose = verbose
        )
      }
    }

    scale_features <- tryCatch(
      rownames(GetAssayData5(
        srt_integrated,
        layer = "scale.data",
        assay = assay
      )),
      error = function(e) character(0)
    )
    if (
      isTRUE(do_scaling) ||
        (is.null(do_scaling) && any(!HVF %in% scale_features))
    ) {
      log_message(
        "Perform {.fn Seurat::ScaleData} on split layers for Seurat v5 integration",
        verbose = verbose
      )
      srt_integrated <- Seurat::ScaleData(
        object = srt_integrated,
        assay = assay,
        features = HVF,
        vars.to.regress = vars_to_regress,
        model.use = regression_model,
        verbose = FALSE
      )
    }
  } else {
    assay_use <- "SCT"
    if (
      isTRUE(do_normalization) ||
        is.null(do_normalization) ||
        !"SCT" %in% SeuratObject::Assays(srt_integrated)
    ) {
      check_r("glmGamPoi", verbose = FALSE)
      srt_integrated[[assay]] <- SeuratObject::JoinLayers(srt_integrated[[
        assay
      ]])
      srt_integrated[[assay]] <- split(
        x = srt_integrated[[assay]],
        f = srt_integrated[[batch, drop = TRUE]]
      )
      log_message(
        "Perform {.fn SCTransform} per batch for Seurat v5 integration",
        verbose = verbose
      )
      srt_integrated <- SCTransform(
        object = srt_integrated,
        assay = assay,
        variable.features.n = nHVF,
        vars.to.regress = vars_to_regress,
        new.assay.name = "SCT",
        cores = cores %||% 1L,
        verbose = FALSE
      )
    } else {
      SeuratObject::DefaultAssay(srt_integrated) <- "SCT"
    }

    if (is.null(HVF)) {
      HVF <- SeuratObject::VariableFeatures(
        srt_integrated,
        assay = "SCT",
        nfeatures = nHVF
      )
      if (HVF_min_intersection > 1) {
        sct_models <- srt_integrated[["SCT"]]@SCTModel.list
        if (length(sct_models) > 1) {
          per_batch_hvf <- lapply(sct_models, function(model) {
            fa <- Seurat::SCTResults(model, slot = "feature.attributes")
            nf <- min(nHVF, nrow(fa))
            rownames(fa)[utils::head(
              order(fa$residual_variance, decreasing = TRUE),
              nf
            )]
          })
          hvf_freq <- table(unlist(per_batch_hvf))
          HVF <- HVF[HVF %in% names(hvf_freq[hvf_freq >= HVF_min_intersection])]
        }
      }
    } else {
      cf <- rownames(Seurat::GetAssay(srt_integrated, assay = "SCT"))
      HVF <- HVF[HVF %in% cf]
    }

    scale_features <- tryCatch(
      rownames(GetAssayData5(
        srt_integrated,
        layer = "scale.data",
        assay = "SCT"
      )),
      error = function(e) character(0)
    )
    if (length(scale_features) > 0) {
      hvf_missing <- setdiff(HVF, scale_features)
      if (length(hvf_missing) > 0) {
        log_message(
          "Some HVF are absent from SCT scale.data and will be dropped: {.val {hvf_missing}}",
          message_type = "warning",
          verbose = verbose
        )
        HVF <- intersect(HVF, scale_features)
      }
    }

    if (length(HVF) == 0) {
      log_message("No HVF available", message_type = "error")
    }
    log_message(
      "Number of available HVF: {.val {length(HVF)}}",
      verbose = verbose
    )

    SeuratObject::DefaultAssay(srt_integrated) <- "SCT"
  }

  SeuratObject::VariableFeatures(
    object = srt_integrated,
    assay = assay_use
  ) <- HVF

  pca_layer <- if (
    isFALSE(do_scaling) && identical(normalization_method, "LogNormalize")
  ) {
    "data"
  } else {
    "scale.data"
  }
  log_message(
    "Perform {.pkg PCA} on split layers before {.fn Seurat::IntegrateLayers}",
    verbose = verbose
  )
  srt_integrated <- RunDimsReduction(
    srt_integrated,
    prefix = "",
    features = HVF,
    assay = assay_use,
    layer = pca_layer,
    linear_reduction = "pca",
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = TRUE,
    verbose = verbose,
    seed = seed
  )

  params <- list(
    object = srt_integrated,
    method = method_fun,
    orig.reduction = "pca",
    assay = assay_use,
    features = HVF,
    new.reduction = reduction_name,
    normalization.method = normalization_method,
    verbose = FALSE
  )
  for (nm in names(IntegrateLayers_params)) {
    params[[nm]] <- IntegrateLayers_params[[nm]]
  }

  if (method_label %in% c("CCAIntegration", "RPCAIntegration")) {
    assay_obj <- Seurat::GetAssay(
      object = params[["object"]],
      assay = params[["assay"]]
    )
    if (identical(normalization_method, "SCT")) {
      group_names <- levels(assay_obj)
      group_ncells <- vapply(
        group_names,
        function(group_name) {
          length(SeuratObject::Cells(x = assay_obj, layer = group_name))
        },
        integer(1)
      )
    } else {
      group_names <- SeuratObject::Layers(
        object = params[["object"]],
        assay = params[["assay"]],
        search = "data"
      )
      group_ncells <- vapply(
        group_names,
        function(layer_name) {
          ncol(assay_obj[layer_name])
        },
        integer(1)
      )
    }
    if (length(group_ncells) == 0L) {
      log_message(
        "No integration groups available for Seurat v5 workflow",
        message_type = "error"
      )
    }

    max_dims <- min(group_ncells)
    if (identical(method_label, "CCAIntegration")) {
      max_dims <- max_dims - 1L
    }
    if (max_dims < 2L) {
      log_message(
        "At least one batch has too few cells for the selected Seurat v5 integration method",
        message_type = "error"
      )
    }

    dims_use <- params[["dims"]] %||%
      seq_len(min(30L, linear_reduction_dims, max_dims))
    dims_use <- sort(unique(as.integer(dims_use)))
    dims_use <- dims_use[!is.na(dims_use) & dims_use >= 1L]
    if (length(dims_use) == 0L) {
      log_message(
        "{.arg IntegrateLayers_params$dims} must contain positive integers",
        message_type = "error"
      )
    }
    if (max(dims_use) > max_dims) {
      log_message(
        paste(
          "{.arg IntegrateLayers_params$dims} exceeds the supported maximum",
          "{.val {max_dims}} for {.fn {method_label}}"
        ),
        message_type = "error"
      )
    }
    params[["dims"]] <- dims_use

    dims_to_integrate <- params[["dims.to.integrate"]]
    if (!is.null(dims_to_integrate)) {
      dims_to_integrate <- sort(unique(as.integer(dims_to_integrate)))
      dims_to_integrate <- dims_to_integrate[!is.na(dims_to_integrate)]
      if (!all(dims_to_integrate %in% dims_use)) {
        log_message(
          "{.arg IntegrateLayers_params$dims.to.integrate} must be a subset of {.arg IntegrateLayers_params$dims}",
          message_type = "error"
        )
      }
      params[["dims.to.integrate"]] <- dims_to_integrate
    }

    min_group_cells <- min(group_ncells)
    if (
      !is.null(params[["k.weight"]]) &&
        as.integer(params[["k.weight"]]) > min_group_cells
    ) {
      log_message(
        "{.arg IntegrateLayers_params$k.weight} exceeds the size of the smallest batch",
        message_type = "error"
      )
    }
    if (
      !is.null(params[["k.filter"]]) &&
        as.integer(params[["k.filter"]]) >= min_group_cells
    ) {
      log_message(
        "{.arg IntegrateLayers_params$k.filter} must be smaller than the smallest batch size",
        message_type = "error"
      )
    }
  } else {
    params[["normalization.method"]] <- NULL
  }

  log_message(
    "Perform {.pkg Seurat v5} integration with {.fn {method_label}}",
    verbose = verbose
  )
  if (
    method_label %in% c("CCAIntegration", "RPCAIntegration") &&
      requireNamespace("future", quietly = TRUE)
  ) {
    old_future_plan <- future::plan()
    old_future_max_size <- getOption("future.globals.maxSize")
    on.exit(future::plan(old_future_plan), add = TRUE)
    on.exit(options(future.globals.maxSize = old_future_max_size), add = TRUE)
    future::plan(future::sequential)
    options(future.globals.maxSize = Inf)
  }
  srt_integrated <- invoke_fun(Seurat::IntegrateLayers, params)

  if (identical(normalization_method, "LogNormalize")) {
    srt_integrated[[assay_use]] <- SeuratObject::JoinLayers(srt_integrated[[
      assay_use
    ]])
  }

  dims_for_downstream <- if (!is.null(dims_use_override)) {
    dims_use_override
  } else {
    resolve_linear_dims_use(
      srt = srt_integrated,
      reduction = reduction_name,
      linear_reduction_dims_use = linear_reduction_dims_use,
      normalization_method = normalization_method,
      reduction_method = "pca",
      verbose = verbose
    )
  }

  srt_integrated <- find_neighbors_and_clusters(
    srt = srt_integrated,
    reduction = reduction_name,
    dims_use = dims_for_downstream,
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
    prefix = gsub("_$", "", graph_prefix),
    reduction_use = reduction_name,
    reduction_dims = dims_for_downstream,
    graph_use = graph_snn,
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed,
    verbose = verbose
  )

  SeuratObject::DefaultAssay(srt_integrated) <- assay
  SeuratObject::VariableFeatures(srt_integrated) <- HVF
  srt_integrated@misc[[misc_hvf_key]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_output <- srt_merge_raw
    if (
      assay %in%
        SeuratObject::Assays(srt_output) &&
        inherits(Seurat::GetAssay(srt_output, assay = assay), "Assay5")
    ) {
      srt_output[[assay]] <- SeuratObject::JoinLayers(srt_output[[assay]])
    }
    srt_output <- srt_append(
      srt_raw = srt_output,
      srt_append = srt_integrated,
      pattern = paste0(assay, "|", append_tag, "|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    SeuratObject::DefaultAssay(srt_output) <- assay
    SeuratObject::VariableFeatures(srt_output) <- HVF
    return(srt_output)
  }
  srt_integrated
}
