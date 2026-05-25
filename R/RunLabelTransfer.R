#' @title Transfer reference labels to query cells
#'
#' @description
#' Standalone label-transfer workflow for query cells using a reference object.
#' The current implementation is optimized for scATAC query objects mapped to a
#' scRNA-seq reference via gene activity.
#'
#' @md
#' @inheritParams standard_scop
#' @param reference RNA reference `Seurat` object used for label transfer.
#' @param method Label-transfer backend. One of `"Seurat"` or `"scOMM"`.
#' @param prefix Prefix used to resolve ATAC reductions. Default is `"ATAC"`.
#' @param reference_assay Assay used in the reference object.
#' @param reference_reduction Reduction used in the reference object.
#' @param reference_dims Dimensions used from the reference reduction.
#' @param reference_label Metadata column in the reference used as transfer labels.
#' @param add_gene_activity Whether to calculate a gene activity assay for the query.
#' @param gene_activity_assay Name of the gene activity assay used for mapping.
#' @param weight_reduction Reduction in `srt` used to weight transferred labels.
#' If `NULL`, an ATAC linear reduction is resolved automatically from
#' `ATAC_default_linear_reduction`, `{prefix}lsi`, `{prefix}svd`, or the current
#' default reduction.
#' @param dims Query reduction dimensions used by `TransferData`.
#' @param features Features used by `FindTransferAnchors`. If `NULL`, reference
#' variable features are used.
#' @param prediction_prefix Prefix added to prediction metadata columns. If
#' `NULL`, `"predicted_"` is used for `method = "Seurat"` and `"scomm_"` is used
#' for `method = "scOMM"`.
#' @param k.weight Number of neighbors used when weighting transfer anchors.
#' @param evaluate Whether to compute mapping metrics against a truth label.
#' @param truth_col Metadata column in `srt` used as the truth label when
#' `evaluate = TRUE`.
#' @param tool_name Name used to store detailed results in `srt@tools`.
#' @param rare_threshold Maximum class proportion used to define rare classes
#' when calculating `rare_recall`.
#' @param scomm_python Optional Python binary used by the `scOMM` backend.
#' If `NULL`, `SCOP_SCOMM_PYTHON` is consulted and reticulate defaults are used
#' otherwise.
#' @param scomm_hidden_nodes,scomm_epochs,scomm_batch_size,scomm_threshold,scomm_seed
#' Parameters passed to the optional `scOMM` backend.
#'
#' @return A `Seurat` object with prediction metadata added.
#' @export
#' @examples
#' data("pbmcmultiome_sub", package = "scop")
#' pbmcmultiome_sub <- standard_scop(
#'   pbmcmultiome_sub,
#'   assay = "RNA",
#'   linear_reduction_dims = 20
#' )
#' reference <- subset(pbmcmultiome_sub, cells = colnames(pbmcmultiome_sub)[1:250])
#' query <- subset(pbmcmultiome_sub, cells = colnames(pbmcmultiome_sub)[251:350])
#' query <- standard_scop(
#'   query,
#'   assay = "peaks",
#'   normalization_method = "TFIDF",
#'   linear_reduction_dims = 20
#' )
#' query <- RunLabelTransfer(
#'   srt = query,
#'   reference = reference,
#'   assay = "peaks",
#'   reference_assay = "RNA",
#'   reference_reduction = "Standardpca",
#'   reference_label = "CellType",
#'   reference_dims = 1:10,
#'   dims = 2:10
#' )
RunLabelTransfer <- function(
  srt,
  reference,
  assay = NULL,
  method = c("Seurat", "scOMM"),
  prefix = "ATAC",
  reference_assay = NULL,
  reference_reduction = "pca",
  reference_dims = 1:30,
  reference_label = NULL,
  add_gene_activity = TRUE,
  gene_activity_assay = "ACTIVITY",
  weight_reduction = NULL,
  dims = 2:30,
  features = NULL,
  prediction_prefix = NULL,
  k.weight = 100,
  evaluate = FALSE,
  truth_col = NULL,
  tool_name = NULL,
  rare_threshold = 0.05,
  scomm_python = NULL,
  scomm_hidden_nodes = c(128, 64),
  scomm_epochs = 10,
  scomm_batch_size = 32,
  scomm_threshold = 0.5,
  scomm_seed = 11,
  verbose = TRUE
) {
  method <- match.arg(method)
  prediction_prefix <- prediction_prefix %||%
    if (identical(method, "scOMM")) {
      "scomm_"
    } else {
      "predicted_"
    }
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  reference_assay <- reference_assay %||% SeuratObject::DefaultAssay(reference)
  if (!inherits(srt[[assay]], "ChromatinAssay")) {
    log_message(
      "{.arg assay} must refer to a {.cls ChromatinAssay}",
      message_type = "error"
    )
  }
  if (
    isTRUE(add_gene_activity) &&
      !gene_activity_assay %in% SeuratObject::Assays(srt) &&
      reference_assay %in% SeuratObject::Assays(srt)
  ) {
    log_message(
      "Use existing query assay {.val {reference_assay}} as {.arg gene_activity_assay}",
      verbose = verbose
    )
    gene_activity_assay <- reference_assay
  }
  tool_name <- tool_name %||% paste0(prefix, "_", method, "_LabelTransfer")
  if (identical(method, "Seurat")) {
    srt <- atac_transfer_labels(
      srt = srt,
      assay = assay,
      prefix = prefix,
      linear_reduction_dims_use = dims,
      reference = reference,
      reference_assay = reference_assay,
      reference_reduction = reference_reduction,
      reference_dims = reference_dims,
      reference_label = reference_label,
      add_gene_activity = add_gene_activity,
      gene_activity_assay = gene_activity_assay,
      weight_reduction = weight_reduction,
      features = features,
      prediction_prefix = prediction_prefix,
      k.weight = k.weight,
      verbose = verbose
    )
    predicted_col <- resolve_pred_col(
      srt = srt,
      prediction_prefix = prediction_prefix
    )
    probability_col <- resolve_prob_col(
      srt = srt,
      prediction_prefix = prediction_prefix
    )
    probability_cols <- grep(
      paste0("^", prediction_prefix, "prediction\\.score\\."),
      colnames(srt@meta.data),
      value = TRUE
    )
    used_features <- features
  } else {
    if (isTRUE(add_gene_activity) && !gene_activity_assay %in% SeuratObject::Assays(srt)) {
      srt <- atac_add_activity(
        srt = srt,
        assay = assay,
        gene_activity_assay = gene_activity_assay,
        verbose = verbose
      )
    }
    query_assay <- if (gene_activity_assay %in% SeuratObject::Assays(srt)) {
      gene_activity_assay
    } else {
      assay
    }
    srt <- RunscOMM(
      srt = srt,
      reference = reference,
      reference_assay = reference_assay,
      query_assay = query_assay,
      reference_label = reference_label,
      features = features,
      prediction_prefix = prediction_prefix,
      evaluate = FALSE,
      tool_name = tool_name,
      scomm_python = scomm_python,
      scomm_hidden_nodes = scomm_hidden_nodes,
      scomm_epochs = scomm_epochs,
      scomm_batch_size = scomm_batch_size,
      scomm_threshold = scomm_threshold,
      scomm_seed = scomm_seed,
      verbose = verbose
    )
    tool_info <- srt@tools[[tool_name]]
    predicted_col <- tool_info$predicted_col
    probability_col <- tool_info$probability_col
    probability_cols <- tool_info$probability_cols
    used_features <- tool_info$features
  }
  if (isTRUE(evaluate)) {
    if (is.null(truth_col) || !truth_col %in% colnames(srt@meta.data)) {
      log_message(
        "{.arg truth_col} must be provided in {.arg srt@meta.data} when {.arg evaluate = TRUE}",
        message_type = "error"
      )
    }
    eval_res <- collect_mapping_metrics(
      srt = srt,
      predicted_col = predicted_col,
      truth_col = truth_col,
      probability_col = probability_col,
      rare_threshold = rare_threshold
    )
  } else {
    eval_res <- NULL
  }
  srt@tools[[tool_name]] <- list(
    method = method,
    predicted_col = predicted_col,
    probability_col = probability_col,
    probability_cols = probability_cols,
    truth_col = truth_col,
    metrics = eval_res,
    reference_assay = reference_assay,
    features = used_features
  )
  srt
}

resolve_pred_col <- function(srt, prediction_prefix = "predicted_") {
  candidates <- c(
    paste0(prediction_prefix, "predicted.id"),
    paste0(prediction_prefix, "prediction"),
    grep(
      paste0("^", prediction_prefix, "(predicted\\.id|prediction)$"),
      colnames(srt@meta.data),
      value = TRUE
    )
  )
  candidates <- unique(candidates)
  candidates[candidates %in% colnames(srt@meta.data)][1] %||% NULL
}

resolve_prob_col <- function(srt, prediction_prefix = "predicted_") {
  candidates <- c(
    paste0(prediction_prefix, "score.max"),
    paste0(prediction_prefix, "prediction.score.max"),
    grep(
      paste0("^", prediction_prefix, "(score\\.max|prediction\\.score\\.max)$"),
      colnames(srt@meta.data),
      value = TRUE
    )
  )
  candidates <- unique(candidates)
  candidates[candidates %in% colnames(srt@meta.data)][1] %||% NULL
}

atac_transfer_labels <- function(
  srt,
  assay,
  prefix,
  linear_reduction_dims_use,
  reference,
  reference_assay = NULL,
  reference_reduction = "pca",
  reference_dims = 1:30,
  reference_label = NULL,
  add_gene_activity = TRUE,
  gene_activity_assay = "ACTIVITY",
  weight_reduction = NULL,
  features = NULL,
  prediction_prefix = "predicted_",
  k.weight = 100,
  verbose = TRUE
) {
  if (!inherits(reference, "Seurat")) {
    log_message(
      "{.arg reference} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  if (is.null(reference_label) || !reference_label %in% colnames(reference@meta.data)) {
    log_message(
      "{.arg reference_label} must be a valid metadata column in {.arg reference}",
      message_type = "error"
    )
  }

  reference_assay <- reference_assay %||% SeuratObject::DefaultAssay(reference)
  features <- features %||% SeuratObject::VariableFeatures(reference, assay = reference_assay)
  weight_reduction <- atac_weight_reduction(
    srt = srt,
    prefix = prefix,
    weight_reduction = weight_reduction,
    verbose = verbose
  )
  available_dims <- seq_len(ncol(Seurat::Embeddings(srt, reduction = weight_reduction)))
  linear_reduction_dims_use <- intersect(linear_reduction_dims_use, available_dims)
  if (length(linear_reduction_dims_use) == 0) {
    log_message(
      "No valid dimensions remain after checking {.arg weight_reduction}",
      message_type = "error"
    )
  }

  if (isTRUE(add_gene_activity) && !gene_activity_assay %in% SeuratObject::Assays(srt)) {
    srt <- atac_add_activity(
      srt = srt,
      assay = assay,
      gene_activity_assay = gene_activity_assay,
      verbose = verbose
    )
  }

  query_assay <- if (gene_activity_assay %in% SeuratObject::Assays(srt)) {
    gene_activity_assay
  } else {
    assay
  }
  features <- Reduce(
    intersect,
    list(
      features,
      rownames(reference[[reference_assay]]),
      rownames(srt[[query_assay]])
    )
  )
  if (length(features) == 0) {
    log_message(
      "No shared features between {.arg reference_assay} and {.arg query_assay} for ATAC label transfer",
      message_type = "error"
    )
  }
  anchor_params <- atac_anchor_params(
    reference_cells = ncol(reference),
    query_cells = ncol(srt),
    verbose = verbose
  )

  anchor_reduction <- if (!is.null(reference_reduction)) "pcaproject" else "cca"
  log_message("Running RNA reference label transfer for ATAC cells...", verbose = verbose)
  ref_labels <- reference[[reference_label]]
  if (inherits(ref_labels, "data.frame")) {
    ref_labels <- ref_labels[[1]]
  }
  anchors <- Seurat::FindTransferAnchors(
    reference = reference,
    query = srt,
    reference.assay = reference_assay,
    query.assay = query_assay,
    reduction = anchor_reduction,
    reference.reduction = reference_reduction,
    dims = reference_dims,
    features = features,
    k.anchor = anchor_params$k.anchor,
    k.filter = anchor_params$k.filter,
    k.score = anchor_params$k.score,
    verbose = FALSE
  )
  k.weight_use <- atac_k_weight(
    anchors = anchors,
    query_cells = ncol(srt),
    k.weight = k.weight,
    verbose = verbose
  )
  predictions <- Seurat::TransferData(
    anchorset = anchors,
    refdata = ref_labels,
    weight.reduction = srt[[weight_reduction]],
    dims = linear_reduction_dims_use,
    k.weight = k.weight_use,
    verbose = FALSE
  )
  colnames(predictions) <- paste0(prediction_prefix, colnames(predictions))
  srt <- SeuratObject::AddMetaData(srt, metadata = predictions)
  srt@misc[["scop_mode"]] <- "atac"
  srt
}

atac_add_activity <- function(
  srt,
  assay,
  gene_activity_assay = "ACTIVITY",
  verbose = TRUE
) {
  if (gene_activity_assay %in% SeuratObject::Assays(srt)) {
    return(srt)
  }

  log_message("Calculating gene activity assay...", verbose = verbose)
  ga <- tryCatch(
    Signac::GeneActivity(object = srt, assay = assay),
    error = function(error) error
  )

  if (inherits(ga, "error")) {
    log_message(
      "Unable to calculate gene activity assay: {.val {conditionMessage(ga)}}. Please provide a valid ATAC annotation.",
      message_type = "error"
    )
  }

  srt[[gene_activity_assay]] <- Seurat::CreateAssayObject(counts = ga)
  srt <- Seurat::NormalizeData(
    object = srt,
    assay = gene_activity_assay,
    normalization.method = "LogNormalize",
    verbose = FALSE
  )
  srt
}

atac_anchor_params <- function(
  reference_cells,
  query_cells,
  k.anchor = 5,
  k.filter = 200,
  k.score = 30,
  verbose = TRUE
) {
  max_shared <- max(1L, min(as.integer(reference_cells), as.integer(query_cells)) - 1L)
  params <- list(
    k.anchor = min(as.integer(k.anchor), max_shared),
    k.filter = min(as.integer(k.filter), max_shared),
    k.score = min(as.integer(k.score), max_shared)
  )
  if (!identical(params$k.anchor, as.integer(k.anchor))) {
    log_message(
      "Adjust {.arg k.anchor} from {.val {k.anchor}} to {.val {params$k.anchor}} for small-sample ATAC mapping",
      verbose = verbose
    )
  }
  if (!identical(params$k.filter, as.integer(k.filter))) {
    log_message(
      "Adjust {.arg k.filter} from {.val {k.filter}} to {.val {params$k.filter}} for small-sample ATAC mapping",
      verbose = verbose
    )
  }
  if (!identical(params$k.score, as.integer(k.score))) {
    log_message(
      "Adjust {.arg k.score} from {.val {k.score}} to {.val {params$k.score}} for small-sample ATAC mapping",
      verbose = verbose
    )
  }
  params
}

atac_k_weight <- function(
  anchors,
  query_cells,
  k.weight = 100,
  verbose = TRUE
) {
  anchor_matrix <- tryCatch(
    anchors@anchors,
    error = function(...) NULL
  )
  anchor_n <- if (!is.null(anchor_matrix)) {
    nrow(anchor_matrix)
  } else {
    NA_integer_
  }
  anchor_ref_cells <- if (!is.null(anchor_matrix) && ncol(anchor_matrix) >= 1) {
    length(unique(anchor_matrix[, 1]))
  } else {
    NA_integer_
  }
  anchor_query_cells <- if (!is.null(anchor_matrix) && ncol(anchor_matrix) >= 2) {
    length(unique(anchor_matrix[, 2]))
  } else {
    NA_integer_
  }
  limits <- c(
    as.integer(k.weight),
    as.integer(query_cells) - 1L,
    as.integer(anchor_n) - 1L,
    as.integer(anchor_ref_cells) - 1L,
    as.integer(anchor_query_cells) - 1L
  )
  limits <- limits[is.finite(limits) & !is.na(limits)]
  limits <- limits[limits > 0]
  if (length(limits) == 0) {
    return(1L)
  }
  k.weight_use <- min(limits)
  if (!identical(k.weight_use, as.integer(k.weight))) {
    log_message(
      "Adjust {.arg k.weight} from {.val {k.weight}} to {.val {k.weight_use}} for small-sample ATAC mapping",
      verbose = verbose
    )
  }
  as.integer(k.weight_use)
}

atac_weight_reduction <- function(
  srt,
  prefix = NULL,
  weight_reduction = NULL,
  verbose = TRUE
) {
  available_reductions <- SeuratObject::Reductions(srt)
  if (!is.null(weight_reduction)) {
    if (!weight_reduction %in% available_reductions) {
      log_message(
        "{.arg weight_reduction} is not present in {.arg srt}",
        message_type = "error"
      )
    }
    return(weight_reduction)
  }

  default_reduction <- tryCatch(
    DefaultReduction(srt),
    error = function(...) NULL
  )
  candidate_reductions <- unique(c(
    srt@misc[["ATAC_default_linear_reduction"]] %||% NULL,
    if (!is.null(prefix)) {
      c(paste0(prefix, "lsi"), paste0(prefix, "svd"))
    },
    tryCatch(
      DefaultReduction(srt, pattern = "(lsi|svd)$"),
      error = function(...) NULL
    ),
    default_reduction
  ))
  candidate_reductions <- candidate_reductions[candidate_reductions %in% available_reductions]
  if (length(candidate_reductions) == 0) {
    log_message(
      "Unable to resolve an ATAC linear reduction for {.arg weight_reduction}",
      message_type = "error"
    )
  }
  log_message(
    "Use {.val {candidate_reductions[[1]]}} as the ATAC weight reduction",
    verbose = verbose
  )
  candidate_reductions[[1]]
}
