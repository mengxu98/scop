#' @title Transfer reference labels to query cells
#'
#' @description
#' Standalone label-transfer workflow for query cells using a reference object.
#' The current implementation is optimized for scATAC query objects mapped to a
#' scRNA-seq reference via gene activity.
#'
#' @md
#' @inheritParams RunStandardWorkflow
#' @param reference RNA reference `Seurat` object used for label transfer.
#' @param method Label-transfer backend. One of `"Seurat"` or `"scOMM"`.
#' @param prefix Prefix used to resolve ATAC reductions.
#' @param reference_assay Assay used in the reference object.
#' @param reference_reduction Reduction used in the reference object.
#' @param reference_dims Dimensions used from the reference reduction.
#' @param reference_label Metadata column in the reference used as transfer labels.
#' @param add_gene_activity Whether to calculate a gene activity assay for the query.
#' If the requested activity assay is absent, an existing query assay named
#' `reference_assay` is used before attempting gene activity calculation. Set
#' `FALSE` with an explicit existing `gene_activity_assay` to avoid this fallback.
#' @param gene_activity_assay Name of the gene activity assay used for mapping.
#' @param weight_reduction Reduction in `object` used to weight transferred labels.
#' If `NULL`, an ATAC linear reduction is resolved automatically from
#' `ATAC_default_linear_reduction`, `{prefix}lsi`, `{prefix}svd`, or the current
#' default reduction.
#' @param dims Query reduction dimensions used by `TransferData`.
#' @param features Features used by `FindTransferAnchors`. If `NULL`, reference
#' variable features are used.
#' @param prediction_prefix Prefix added to prediction metadata columns. If
#' `NULL`, `"predicted_"` is used for `method = "Seurat"` and `"scomm_"` is used
#' for `method = "scOMM"`. Seurat prefixes are literal strings. A rerun with
#' the same tool and prefix replaces only columns recorded by that prior result;
#' collisions with unrecorded metadata columns are rejected.
#' @param k.weight Number of neighbors used when weighting transfer anchors.
#' @param evaluate Whether to compute mapping metrics against a truth label.
#' @param truth_col Metadata column in `object` used as the truth label when
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
#' @return A `Seurat` object with prediction metadata added. For the Seurat
#' backend, the tool result records actual anchor features and current prediction
#' columns. This is a calculation-time snapshot, not automatically updated by
#' later subsetting or renaming. Caller RNG state and original cell IDs are
#' preserved; temporary query IDs isolate native suffix handling. Predictions
#' describe similarity to the supplied reference, not independently established
#' biological accuracy.
#' @export
#' @examples
#' # A small synthetic technical example with independent reference/query cells.
#' # ACTIVITY below is synthetic input, not measured RNA or inferred gene activity.
#' set.seed(91)
#' genes <- paste0("GENE", seq_len(150))
#' make_counts <- function(n, prefix) {
#'   mu <- matrix(8, 150, n)
#'   mu[1:30, seq_len(n / 2)] <- 24
#'   mu[31:60, (n / 2 + 1):n] <- 24
#'   Matrix::Matrix(matrix(rnbinom(length(mu), mu = mu, size = 10), 150,
#'     dimnames = list(genes, paste0(prefix, seq_len(n)))
#'   ), sparse = TRUE)
#' }
#' reference <- SeuratObject::CreateSeuratObject(make_counts(60, "ref"))
#' reference$label <- rep(c("classA", "classB"), each = 30)
#' reference <- Seurat::NormalizeData(reference, verbose = FALSE)
#' reference <- Seurat::FindVariableFeatures(reference, nfeatures = 150, verbose = FALSE)
#' reference <- Seurat::ScaleData(reference, features = genes, verbose = FALSE)
#' reference <- Seurat::RunPCA(reference, features = genes, npcs = 10, verbose = FALSE)
#' activity <- make_counts(40, "query")
#' peaks <- Matrix::Matrix(matrix(rpois(200 * 40, 4), 200,
#'   dimnames = list(paste0(
#'     "chr1-", seq(1, by = 200, length.out = 200), "-",
#'     seq(100, by = 200, length.out = 200)
#'   ), colnames(activity))
#' ), sparse = TRUE)
#' query <- SeuratObject::CreateSeuratObject(
#'   Signac::CreateChromatinAssay(counts = peaks),
#'   assay = "peaks"
#' )
#' query[["ACTIVITY"]] <- SeuratObject::CreateAssayObject(counts = activity)
#' query <- Seurat::NormalizeData(query, assay = "ACTIVITY", verbose = FALSE)
#' query <- Signac::RunTFIDF(query, assay = "peaks", verbose = FALSE)
#' query <- Signac::FindTopFeatures(query, assay = "peaks", min.cutoff = "q0", verbose = FALSE)
#' query <- Signac::RunSVD(query, assay = "peaks", n = 10, verbose = FALSE)
#' query <- RunLabelTransfer(query, reference,
#'   assay = "peaks",
#'   reference_assay = "RNA", reference_reduction = "pca", reference_dims = 1:5,
#'   reference_label = "label", add_gene_activity = FALSE,
#'   gene_activity_assay = "ACTIVITY", weight_reduction = "lsi", dims = 2:5,
#'   features = genes, k.weight = 5, verbose = FALSE
#' )
#' head(query$predicted_predicted.id)
RunLabelTransfer <- function(
  object,
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
  verbose = TRUE,
  srt = NULL
) {
  srt <- resolve_deprecated_srt(object, srt, missing(object))
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
      verbose = verbose,
      tool_name = tool_name
    )
    current <- srt@tools[[tool_name]]
    predicted_col <- current$predicted_col
    probability_col <- current$probability_col
    probability_cols <- current$probability_cols
    used_features <- current$features
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
      object = srt,
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
  verbose = TRUE,
  tool_name = NULL
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

  seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (seed_exists) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit(
    {
      if (seed_exists) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
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
  # TransferData removes every literal "_query" token, including user ID content.
  # Native temporary names avoid both lost labels and pre-writeback collisions.
  query_cells <- colnames(srt)
  reference_cells <- colnames(reference)
  safe_cells <- make.unique(c(reference_cells, paste0("scopQuery", seq_along(query_cells))))[
    length(reference_cells) + seq_along(query_cells)
  ]
  query <- SeuratObject::RenameCells(srt, new.names = stats::setNames(safe_cells, query_cells))
  anchors <- Seurat::FindTransferAnchors(
    reference = reference,
    query = query,
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
    weight.reduction = query[[weight_reduction]],
    dims = linear_reduction_dims_use,
    k.weight = k.weight_use,
    verbose = FALSE
  )
  prediction_cells <- rownames(predictions)
  if (length(prediction_cells) != length(safe_cells) || anyDuplicated(prediction_cells) ||
    !setequal(prediction_cells, safe_cells)) {
    log_message("Native label-transfer predictions must cover every query cell exactly once.", message_type = "error")
  }
  predictions <- predictions[match(safe_cells, prediction_cells), , drop = FALSE]
  rownames(predictions) <- query_cells
  colnames(predictions) <- paste0(prediction_prefix, colnames(predictions))
  predicted_col <- paste0(prediction_prefix, "predicted.id")
  probability_col <- paste0(prediction_prefix, "prediction.score.max")
  probability_cols <- colnames(predictions)[startsWith(colnames(predictions), paste0(prediction_prefix, "prediction.score."))]
  if (!is.null(tool_name)) {
    previous <- srt@tools[[tool_name]]
    previous_columns <- if (is.list(previous) && identical(previous$method, "Seurat") &&
      identical(previous$predicted_col, predicted_col)) {
      unique(c(previous$predicted_col, previous$probability_col, previous$probability_cols))
    } else {
      character()
    }
    previous_columns <- intersect(previous_columns, colnames(srt@meta.data))
    collisions <- intersect(colnames(predictions), setdiff(colnames(srt@meta.data), previous_columns))
    if (length(collisions)) {
      log_message("Prediction columns already exist outside this tool's recorded results; use an unused prediction_prefix.", message_type = "error")
    }
    srt@meta.data <- srt@meta.data[, setdiff(colnames(srt@meta.data), previous_columns), drop = FALSE]
  }
  srt <- SeuratObject::AddMetaData(srt, metadata = predictions)
  if (!is.null(tool_name)) {
    srt@tools[[tool_name]] <- list(
      predicted_col = predicted_col, probability_col = probability_col,
      probability_cols = probability_cols, features = anchors@anchor.features
    )
  }
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

  check_r("Signac", verbose = FALSE)

  log_message("Calculating gene activity assay...", verbose = verbose)
  ga <- tryCatch(
    get_namespace_fun("Signac", "GeneActivity")(object = srt, assay = assay),
    error = function(error) error
  )

  if (inherits(ga, "error")) {
    log_message(
      "Unable to calculate gene activity assay: {.val {conditionMessage(ga)}}. Please provide a valid ATAC annotation.",
      message_type = "error"
    )
  }

  srt[[gene_activity_assay]] <- Seurat::CreateAssayObject(counts = ga)
  srt <- NormalizeData(
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
