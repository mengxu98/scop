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
#' @param assay Which assay to use.
#' Default is `"RNA"`.
#' @param model Model name or path. Default is `"Immune_All_Low.pkl"`.
#' Supports three formats:
#' 1. Model name (e.g., `"Immune_All_Low.pkl"`): automatically searched in `~/.celltypist/data/models/`
#' 2. Full path (contains `/`): use the provided path directly
#' 3. `NULL`: use default model
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
#' @seealso
#' [CellTypistModels], [RunSingleR], [RunScmap]
#'
#' @return
#' An AnnData object or a Seurat object depending on the `return_seurat` argument.
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
    verbose = TRUE) {
  log_message(
    "Running {.pkg CellTypist} annotation...",
    verbose = verbose
  )
  PrepareEnv()
  if (all(is.null(srt), is.null(adata))) {
    log_message(
      "One of {.arg srt} or {.arg adata} must be provided",
      message_type = "error"
    )
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

#' @title Get available CellTypist models
#'
#' @description
#' Get a list of all available CellTypist models,
#' either downloaded locally or available from the server.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param on_the_fly If `TRUE`, only show downloaded models.
#' If `FALSE`, show all available models (fetch list from server).
#' Default is `FALSE`.
#'
#' @return
#' A data frame containing model names and descriptions.
#'
#' @seealso
#' [RunCellTypist]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Get available models
#' models <- CellTypistModels()
#' print(models)
#' }
CellTypistModels <- function(
    on_the_fly = FALSE,
    verbose = TRUE) {
  log_message("Fetching CellTypist models...", verbose = verbose)
  PrepareEnv()
  functions <- reticulate::import_from_path(
    "functions",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )

  models_df <- functions$CellTypistModels(
    on_the_fly = on_the_fly,
    verbose = verbose
  )

  if (!is.null(models_df)) {
    models_r <- reticulate::py_to_r(models_df)
    return(models_r)
  } else {
    return(NULL)
  }
}
