#' @title Run CellTypist cell type annotation
#'
#' @description
#' CellTypist is an automated cell type annotation tool for scRNA-seq datasets
#' based on logistic regression classifiers. This function runs CellTypist annotation
#' on a Seurat object or AnnData object.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams RunCellRank
#' @inheritParams FeatureDimPlot
#' @param adata Optional AnnData object used as input.
#' @param assay Which assay to use.
#' Default is `"RNA"`.
#' @param model Model name or path. Default is `"Immune_All_Low.pkl"`.
#' Supports three formats:
#' 1. Model name (e.g., `"Immune_All_Low.pkl"`): automatically searched in `~/.celltypist/data/models/`
#' 2. Full path (contains `/`): use the provided path directly
#' 3. `NULL`: use default model
#' 4. A summary list returned by [TrainCellTypist()]: use its `model_path`
#' @param mode Prediction mode: `"best match"` or `"prob match"`.
#' Default is `"best match"`.
#' @param p_thres Probability threshold for `"prob match"` mode.
#' Default is `0.5`.
#' @param majority_voting Whether to use majority voting.
#' Default is `FALSE`.
#' @param over_clustering Over-clustering result. Can be:
#' - String: column name in Seurat metadata or AnnData obs
#' - Vector: over-clustering labels
#' - `NULL`: use heuristic over-clustering
#' @param min_prop Minimum proportion for majority voting.
#' Default is `0`.
#' @param use_GPU Whether to use GPU for over-clustering.
#' Default is `FALSE`.
#' @param insert_labels Whether to insert predicted labels.
#' Default is `TRUE`.
#' @param insert_conf Whether to insert confidence scores.
#' Default is `TRUE`.
#' @param insert_conf_by Which prediction type to base confidence on.
#' Default is `"predicted_labels"`.
#' @param insert_prob Whether to insert probability matrix.
#' Default is `FALSE`.
#' @param insert_decision Whether to insert decision matrix.
#' Default is `FALSE`.
#' @param prefix Prefix for inserted columns.
#' Default is `"celltypist_"`.
#'
#' @return
#' An AnnData object or a Seurat object depending on the `return_seurat` argument.
#'
#' @seealso
#' [CellTypistModels], [TrainCellTypist], [RunSingleR], [RunScmap]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunCellTypist(
#'   pancreas_sub,
#'   model = "Developing_Mouse_Brain.pkl"
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "celltypist_predicted_labels",
#'   legend.position = "none"
#' )
#'
#' # Use prob match mode
#' pancreas_sub <- RunCellTypist(
#'   pancreas_sub,
#'   model = "Developing_Mouse_Brain.pkl",
#'   mode = "prob match",
#'   p_thres = 0.5
#' )
#'
#' # Use majority voting
#' pancreas_sub <- RunCellTypist(
#'   pancreas_sub,
#'   model = "Developing_Mouse_Brain.pkl",
#'   majority_voting = TRUE
#' )
#' }
RunCellTypist <- function(
  srt = NULL,
  adata = NULL,
  assay = "RNA",
  layer = "data",
  model = "Immune_All_Low.pkl",
  mode = "best match",
  p_thres = 0.5,
  majority_voting = FALSE,
  over_clustering = NULL,
  min_prop = 0,
  use_GPU = FALSE,
  insert_labels = TRUE,
  insert_conf = TRUE,
  insert_conf_by = "predicted_labels",
  insert_prob = FALSE,
  insert_decision = FALSE,
  prefix = "celltypist_",
  return_seurat = !is.null(srt),
  verbose = TRUE
) {
  log_message(
    "Running {.pkg CellTypist} annotation...",
    verbose = verbose
  )
  PrepareEnv()
  check_python("celltypist==1.7.1", verbose = verbose)
  if (all(is.null(srt), is.null(adata))) {
    log_message(
      "One of {.arg srt} or {.arg adata} must be provided",
      message_type = "error"
    )
  }

  if (is.list(model)) {
    if (is.null(model[["model_path"]]) || !nzchar(model[["model_path"]])) {
      log_message(
        "{.arg model} summary list must contain a valid {.field model_path}",
        message_type = "error"
      )
    }
    model <- model[["model_path"]]
  }

  args <- mget(names(formals()))
  args <- lapply(
    args, function(x) {
      if (is.numeric(x)) {
        y <- ifelse(
          grepl(
            "\\.",
            as.character(x)
          ),
          as.double(x),
          as.integer(x)
        )
      } else {
        y <- x
      }
      y
    }
  )
  call_envir <- parent.frame(1)
  args <- lapply(
    args, function(arg) {
      if (is.symbol(arg)) {
        eval(arg, envir = call_envir)
      } else if (is.call(arg)) {
        eval(arg, envir = call_envir)
      } else {
        arg
      }
    }
  )

  params <- c(
    "srt",
    "assay",
    "layer",
    "return_seurat"
  )
  args <- args[!names(args) %in% params]

  if (!is.null(over_clustering)) {
    if (is.character(over_clustering) && length(over_clustering) == 1) {
      if (!is.null(srt)) {
        if (over_clustering %in% colnames(srt[[]])) {
          args[["over_clustering"]] <- over_clustering
        } else {
          log_message(
            "{.val {over_clustering}} not found in Seurat metadata. Will use heuristic over-clustering",
            message_type = "warning",
            verbose = verbose
          )
          args[["over_clustering"]] <- NULL
        }
      } else {
        args[["over_clustering"]] <- over_clustering
      }
    } else {
      args[["over_clustering"]] <- over_clustering
    }
  }

  if (!is.null(srt)) {
    args[["adata"]] <- srt_to_adata(
      srt = srt,
      assay_x = assay,
      layer_x = layer
    )
    if (!is.null(over_clustering) &&
      is.character(over_clustering) &&
      length(over_clustering) == 1 &&
      over_clustering %in% colnames(srt[[]])) {
      args[["over_clustering"]] <- as.character(srt[[over_clustering]][, 1])
    }
  }

  functions <- reticulate::import_from_path(
    "functions",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )
  adata <- do.call(functions$CellTypist, args)

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      srt_out1 <- srt_append(
        srt_raw = srt,
        srt_append = srt_out
      )
      srt_out2 <- srt_append(
        srt_raw = srt_out1,
        srt_append = srt_out,
        pattern = paste0("^", prefix),
        overwrite = TRUE,
        verbose = FALSE
      )
      return(srt_out2)
    }
  } else {
    return(adata)
  }
}

#' @title Train a CellTypist model
#'
#' @description
#' Train a CellTypist model from a Seurat object, AnnData object, or h5ad file.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams RunCellTypist
#' @param h5ad Optional path to an input `.h5ad` file.
#' @param labels Cell labels used for training. Can be a metadata column name when
#' `srt`/`adata` is supplied, or a vector aligned with cells.
#' @param genes Optional gene names. Usually inferred from the input object.
#' @param transpose_input Whether to transpose the input matrix before training.
#' @param with_mean Whether to center features during scaling.
#' @param check_expression Whether to validate expected CellTypist input format.
#' @param C Inverse regularization strength for logistic regression.
#' @param solver Optional solver passed to CellTypist.
#' @param max_iter Optional maximum iterations.
#' @param n_jobs Number of CPUs used by CellTypist training.
#' @param use_SGD Whether to use SGD training.
#' @param alpha Regularization strength for SGD training.
#' @param mini_batch Whether to enable mini-batch training.
#' @param batch_number Number of batches per epoch for mini-batch training.
#' @param batch_size Batch size for mini-batch training.
#' @param epochs Number of epochs for mini-batch training.
#' @param balance_cell_type Whether to balance cell types during mini-batch training.
#' @param feature_selection Whether to run CellTypist feature selection.
#' @param top_genes Number of top genes used during feature selection.
#' @param date,details,url,source,version Free-text metadata stored in the model.
#' @param model_path Optional output path for the trained model.
#' @param return Return mode. One of `"summary"`, `"path"`, or `"model"`.
#' Default is `"summary"`.
#'
#' @return
#' Depends on `return`: a summary list, the model path, or a Python model object.
#'
#' @seealso
#' [RunCellTypist], [CellTypistModels]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#'
#' model_info <- TrainCellTypist(
#'   srt = pancreas_sub,
#'   labels = "SubCellType",
#'   model_path = tempfile(fileext = ".pkl")
#' )
#'
#' data(panc8_sub)
#' genenames <- make.unique(
#'   thisutils::capitalize(
#'     rownames(panc8_sub),
#'     force_tolower = TRUE
#'   )
#' )
#' names(genenames) <- rownames(panc8_sub)
#' panc8_sub <- RenameFeatures(
#'   panc8_sub,
#'   newnames = genenames
#' )
#' panc8_sub <- standard_scop(panc8_sub)
#'
#' pancreas_sub <- RunCellTypist(
#'   srt = pancreas_sub,
#'   model = model_info
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = c("SubCellType", "celltypist_predicted_labels")
#' )
#'
#' panc8_sub <- RunCellTypist(
#'   srt = panc8_sub,
#'   model = model_info
#' )
#' CellDimPlot(
#'   panc8_sub,
#'   group.by = c("celltype", "celltypist_predicted_labels")
#' )
#'
#' ht <- CellCorHeatmap(
#'   srt_query = panc8_sub,
#'   srt_ref = panc8_sub,
#'   query_group = "celltypist_predicted_labels",
#'   ref_group = "celltype",
#'   width = 4,
#'   height = 3
#' )
#' ht$plot
#' }
TrainCellTypist <- function(
  srt = NULL,
  adata = NULL,
  h5ad = NULL,
  assay = "RNA",
  layer = "data",
  labels,
  genes = NULL,
  transpose_input = FALSE,
  with_mean = TRUE,
  check_expression = TRUE,
  C = 1,
  solver = NULL,
  max_iter = NULL,
  n_jobs = 1,
  use_SGD = FALSE,
  alpha = 1e-04,
  use_GPU = FALSE,
  mini_batch = FALSE,
  batch_number = 100,
  batch_size = 1000,
  epochs = 10,
  balance_cell_type = FALSE,
  feature_selection = FALSE,
  top_genes = 300,
  date = "",
  details = "",
  url = "",
  source = "",
  version = "",
  model_path = NULL,
  return = c("summary", "path", "model"),
  verbose = TRUE
) {
  return <- match.arg(return)
  log_message(
    "Training {.pkg CellTypist} model...",
    verbose = verbose
  )
  PrepareEnv()
  check_python("celltypist==1.7.1", verbose = verbose)

  if (is.null(srt) && is.null(adata) && is.null(h5ad)) {
    log_message(
      "One of {.arg srt}, {.arg adata}, or {.arg h5ad} must be provided",
      message_type = "error"
    )
  }

  if (!is.null(srt) && !inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }

  if (!is.null(srt) && !is.character(labels)) {
    if (length(labels) != ncol(srt)) {
      log_message(
        "{.arg labels} length must match the number of cells in {.arg srt}",
        message_type = "error"
      )
    }
  }

  if (!is.null(adata) && !is.character(labels)) {
    n_obs <- tryCatch(adata$n_obs, error = function(e) NULL)
    if (!is.null(n_obs) && length(labels) != n_obs) {
      log_message(
        "{.arg labels} length must match the number of cells in {.arg adata}",
        message_type = "error"
      )
    }
  }

  if (!is.null(srt)) {
    adata <- srt_to_adata(
      srt = srt,
      assay_x = assay,
      layer_x = layer
    )
    if (is.character(labels) && length(labels) == 1) {
      if (!labels %in% colnames(srt[[]])) {
        log_message(
          "{.arg labels} must be a valid metadata column in {.cls Seurat}",
          message_type = "error"
        )
      }
    }
  }

  functions <- reticulate::import_from_path(
    "functions",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )

  model <- functions$TrainCellTypist(
    adata = adata,
    h5ad = h5ad,
    labels = labels,
    genes = genes,
    transpose_input = transpose_input,
    with_mean = with_mean,
    check_expression = check_expression,
    C = C,
    solver = solver,
    max_iter = max_iter,
    n_jobs = as.integer(n_jobs),
    use_SGD = use_SGD,
    alpha = alpha,
    use_GPU = use_GPU,
    mini_batch = mini_batch,
    batch_number = as.integer(batch_number),
    batch_size = as.integer(batch_size),
    epochs = as.integer(epochs),
    balance_cell_type = balance_cell_type,
    feature_selection = feature_selection,
    top_genes = as.integer(top_genes),
    date = date,
    details = details,
    url = url,
    source = source,
    version = version,
    model_path = model_path,
    return_model = identical(return, "model"),
    verbose = verbose
  )

  if (identical(return, "model")) {
    return(model)
  }

  model <- reticulate::py_to_r(model)
  if (identical(return, "path")) {
    return(model$model_path)
  }

  model
}

#' @title Get available CellTypist models
#'
#' @description
#' Unified CellTypist model management interface. Use it to list available models,
#' download models, inspect model metadata, extract model markers, subset models,
#' convert models, or delete local models.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param action Action to perform. One of `"list"`, `"download"`, `"info"`,
#' `"markers"`, `"subset"`, `"convert"`, or `"delete"`.
#' @param on_the_fly If `TRUE`, only show downloaded models.
#' If `FALSE`, show all available models (fetch list from server).
#' Used when `action = "list"`. Default is `FALSE`.
#' @param model Model name(s) used for actions other than `"list"`.
#' @param force_update Whether to refresh the model index before downloading.
#' Used when `action = "download"`.
#' @param cell_type Cell type used when `action = "markers"`.
#' @param top_n Number of markers to extract when `action = "markers"`.
#' @param only_positive Whether to return only positive markers.
#' @param keep_cell_types Cell types retained when `action = "subset"`.
#' @param exclude_cell_types Cell types removed when `action = "subset"`.
#' @param output_model_path Optional output path used by `"subset"` or `"convert"`.
#' If `NULL`, the original model file will be overwritten.
#' @param map_file Optional two-column gene mapping file used by `action = "convert"`.
#' @param sep Delimiter used by `map_file`.
#' @param convert_from,convert_to Column indices used by `map_file`.
#' @param unique_only Whether to use only one-to-one mappings in `action = "convert"`.
#' @param collapse How to collapse one-to-many mappings in `action = "convert"`.
#' @param random_state Random seed used when `collapse = "random"`.
#'
#' @return
#' Depends on `action`: a data frame, a summary list, or a character vector.
#'
#' @seealso
#' [RunCellTypist], [TrainCellTypist]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Get available models
#' models <- CellTypistModels()
#' print(models)
#'
#' # Show downloaded models only
#' local_models <- CellTypistModels(on_the_fly = TRUE)
#'
#' # Download a model
#' CellTypistModels(
#'   action = "download",
#'   model = "Immune_All_Low.pkl"
#' )
#'
#' # Inspect a model
#' CellTypistModels(
#'   action = "info",
#'   model = "Immune_All_Low.pkl"
#' )
#'
#' # Extract top markers for a cell type
#' CellTypistModels(
#'   action = "markers",
#'   model = "Immune_All_Low.pkl",
#'   cell_type = "B cells",
#'   top_n = 20
#' )
#' }
CellTypistModels <- function(
  action = c("list", "download", "info", "markers", "subset", "convert", "delete"),
  on_the_fly = FALSE,
  model = NULL,
  force_update = FALSE,
  cell_type = NULL,
  top_n = 10,
  only_positive = TRUE,
  keep_cell_types = NULL,
  exclude_cell_types = NULL,
  output_model_path = NULL,
  map_file = NULL,
  sep = ",",
  convert_from = NULL,
  convert_to = NULL,
  unique_only = TRUE,
  collapse = c("average", "random"),
  random_state = 0,
  verbose = TRUE
) {
  action <- match.arg(action)
  collapse <- match.arg(collapse)

  log_message(
    "Running {.pkg CellTypist} model action: {.val {action}}",
    verbose = verbose
  )
  PrepareEnv()
  check_python("celltypist==1.7.1", verbose = verbose)
  functions <- reticulate::import_from_path(
    "functions",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )

  result <- functions$CellTypistModels(
    action = action,
    on_the_fly = on_the_fly,
    model = model,
    force_update = force_update,
    cell_type = cell_type,
    top_n = as.integer(top_n),
    only_positive = only_positive,
    keep_cell_types = keep_cell_types,
    exclude_cell_types = exclude_cell_types,
    output_model_path = output_model_path,
    map_file = map_file,
    sep = sep,
    convert_from = convert_from,
    convert_to = convert_to,
    unique_only = unique_only,
    collapse = collapse,
    random_state = as.integer(random_state),
    verbose = verbose
  )

  if (is.null(result)) {
    return(NULL)
  }

  reticulate::py_to_r(result)
}
