#' @title Run scOMM label prediction
#'
#' @description
#' Run `scOMM` on shared features between a reference object and a query object,
#' write predicted labels and class scores into query metadata, and optionally
#' evaluate predictions against a truth label.
#'
#' @md
#' @inheritParams standard_scop
#' @param reference Reference `Seurat` object used for supervision.
#' @param reference_assay Assay used in the reference object.
#' @param query_assay Assay used in the query object.
#' @param reference_label Metadata column in the reference used as supervision labels.
#' @param features Shared features passed to `scOMM`. If `NULL`, reference variable
#' features are used.
#' @param prediction_prefix Prefix added to prediction metadata columns.
#' @param evaluate Whether to compute prediction metrics against a truth label.
#' @param truth_col Metadata column in `srt` used as the truth label when
#' `evaluate = TRUE`.
#' @param tool_name Name used to store detailed results in `srt@tools`.
#' @param rare_threshold Maximum class proportion used to define rare classes
#' when calculating `rare_recall`.
#' @param scomm_python Optional Python binary used by the `scOMM` backend.
#' If `NULL`, `SCOP_SCOMM_PYTHON` is consulted and reticulate defaults are used
#' otherwise.
#' @param scomm_hidden_nodes,scomm_epochs,scomm_batch_size,scomm_threshold,scomm_seed
#' Parameters passed to the `scOMM` backend.
#'
#' @return A `Seurat` object with `scOMM` predictions stored in metadata and `tools`.
#' @export
RunscOMM <- function(
  srt,
  reference,
  reference_assay = NULL,
  query_assay = NULL,
  reference_label = NULL,
  features = NULL,
  prediction_prefix = "predicted_",
  evaluate = FALSE,
  truth_col = NULL,
  tool_name = "scOMM",
  rare_threshold = 0.05,
  scomm_python = NULL,
  scomm_hidden_nodes = c(128, 64),
  scomm_epochs = 10,
  scomm_batch_size = 32,
  scomm_threshold = 0.5,
  scomm_seed = 11,
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  if (!inherits(reference, "Seurat")) {
    log_message(
      "{.arg reference} is not a {.cls Seurat}",
      message_type = "error"
    )
  }

  reference_assay <- reference_assay %||% SeuratObject::DefaultAssay(reference)
  query_assay <- query_assay %||% SeuratObject::DefaultAssay(srt)
  tool_name <- tool_name %||% "scOMM"
  if (!reference_assay %in% SeuratObject::Assays(reference)) {
    log_message(
      "{.arg reference_assay} not found in {.arg reference}",
      message_type = "error"
    )
  }
  if (!query_assay %in% SeuratObject::Assays(srt)) {
    log_message(
      "{.arg query_assay} not found in {.arg srt}",
      message_type = "error"
    )
  }

  scomm_res <- run_scomm(
    reference = reference,
    query = srt,
    reference_assay = reference_assay,
    query_assay = query_assay,
    reference_label = reference_label,
    features = features,
    python = scomm_python,
    hidden_nodes = scomm_hidden_nodes,
    epochs = scomm_epochs,
    batch_size = scomm_batch_size,
    threshold = scomm_threshold,
    seed = scomm_seed,
    verbose = verbose
  )
  pred_res <- add_prediction_meta(
    srt = srt,
    ids = scomm_res$ids,
    probabilities = scomm_res$probabilities,
    prediction_prefix = prediction_prefix
  )
  srt <- pred_res$srt

  if (isTRUE(evaluate)) {
    if (is.null(truth_col) || !truth_col %in% colnames(srt@meta.data)) {
      log_message(
        "{.arg truth_col} must be provided in {.arg srt@meta.data} when {.arg evaluate = TRUE}",
        message_type = "error"
      )
    }
    eval_res <- collect_mapping_metrics(
      srt = srt,
      predicted_col = pred_res$predicted_col,
      truth_col = truth_col,
      probability_col = pred_res$probability_col,
      rare_threshold = rare_threshold
    )
  } else {
    eval_res <- NULL
  }

  srt@tools[[tool_name]] <- list(
    method = "scOMM",
    predicted_col = pred_res$predicted_col,
    probability_col = pred_res$probability_col,
    probability_cols = pred_res$probability_cols,
    truth_col = truth_col,
    metrics = eval_res,
    reference_assay = reference_assay,
    query_assay = query_assay,
    reference_label = reference_label,
    features = scomm_res$features
  )
  srt
}

add_prediction_meta <- function(
  srt,
  ids,
  probabilities = NULL,
  prediction_prefix = "predicted_"
) {
  pred_col <- paste0(prediction_prefix, "prediction")
  srt[[pred_col]] <- factor(ids)
  prob_col <- NULL
  prob_cols <- character(0)
  if (!is.null(probabilities)) {
    prob_df <- as.data.frame(probabilities, check.names = FALSE)
    prob_cols <- paste0(prediction_prefix, "score.", make.names(colnames(prob_df)))
    colnames(prob_df) <- prob_cols
    rownames(prob_df) <- colnames(srt)
    srt <- SeuratObject::AddMetaData(srt, metadata = prob_df)
    prob_col <- paste0(prediction_prefix, "prediction.score.max")
    srt[[prob_col]] <- apply(prob_df, 1, max, na.rm = TRUE)
  }
  list(
    srt = srt,
    predicted_col = pred_col,
    probability_col = prob_col,
    probability_cols = prob_cols
  )
}
