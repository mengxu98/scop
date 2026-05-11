#' @title Annotate single cells using native SciBet
#'
#' @description
#' Run a native `scop` implementation of the core SciBet classifier using a
#' labeled reference `Seurat` object.
#'
#' @md
#' @inheritParams RunSingleR
#' @inheritParams thisutils::log_message
#' @param features Candidate features used by SciBet. If `NULL`, common genes
#' between query and reference are used.
#' @param nfeatures Number of entropy-test features selected from `features`.
#' @param additional_features_per_class Additional high-expression features
#' selected per reference class.
#' @param query_layer,ref_layer Assay layers used for query and reference.
#' @param input_transform How to transform extracted values before SciBet's
#' internal `log1p` step. `"auto"` applies `expm1` to `"data"` layers and no
#' transform otherwise.
#' @param prefix Prefix for metadata columns.
#' @param store_model Whether to store the SciBet core and probabilities in
#' `srt_query@tools`.
#'
#' @return A `Seurat` object with SciBet annotations in metadata and results in
#' `srt_query@tools[["SciBet"]]`.
#' @export
#'
#' @examples
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
#'
#' data(pancreas_sub)
#' pancreas_sub <- RunSciBet(
#'   srt_query = pancreas_sub,
#'   srt_ref = panc8_sub,
#'   ref_group = "celltype",
#'   nfeatures = 200
#' )
#' pancreas_sub <- standard_scop(pancreas_sub, verbose = FALSE)
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = c("SubCellType", "scibet_annotation"),
#'   xlab = "UMAP_1",
#'   ylab = "UMAP_2"
#' )
#'
#' ht <- CellCorHeatmap(
#'   srt_query = pancreas_sub,
#'   srt_ref = pancreas_sub,
#'   query_group = "scibet_annotation",
#'   ref_group = "SubCellType",
#'   width = 3,
#'   height = 3
#' )
#' ht$plot
RunSciBet <- function(
  srt_query,
  srt_ref,
  ref_group,
  query_group = NULL,
  query_assay = NULL,
  ref_assay = NULL,
  query_layer = "counts",
  ref_layer = "counts",
  features = NULL,
  nfeatures = 1000,
  additional_features_per_class = 0,
  input_transform = c("auto", "none", "expm1"),
  prefix = "scibet",
  store_model = TRUE,
  verbose = TRUE
) {
  if (!inherits(srt_query, "Seurat")) {
    log_message(
      "{.arg srt_query} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (!inherits(srt_ref, "Seurat")) {
    log_message(
      "{.arg srt_ref} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  input_transform <- match.arg(input_transform)
  query_assay <- query_assay %||% SeuratObject::DefaultAssay(srt_query)
  ref_assay <- ref_assay %||% SeuratObject::DefaultAssay(srt_ref)

  if (is.null(ref_group)) {
    log_message(
      "{.arg ref_group} must be provided",
      message_type = "error"
    )
  }
  if (length(ref_group) == 1L && ref_group %in% colnames(srt_ref@meta.data)) {
    ref_labels <- srt_ref[[ref_group, drop = TRUE]]
  } else if (length(ref_group) == ncol(srt_ref)) {
    ref_labels <- ref_group
  } else {
    log_message(
      "{.arg ref_group} must be one metadata column or one value per reference cell",
      message_type = "error"
    )
  }
  ref_labels <- as.character(ref_labels)
  keep_ref <- !is.na(ref_labels)
  if (!all(keep_ref)) {
    log_message(
      "Drop {.val {sum(!keep_ref)}} reference cells with missing {.arg ref_group}",
      verbose = verbose
    )
    srt_ref <- srt_ref[, keep_ref]
    ref_labels <- ref_labels[keep_ref]
  }
  ref_labels <- factor(ref_labels)
  if (nlevels(ref_labels) < 2L) {
    log_message(
      "{.arg ref_group} must contain at least two non-missing classes",
      message_type = "error"
    )
  }

  if (!is.null(query_group)) {
    if (
      length(query_group) == 1L &&
        query_group %in% colnames(srt_query@meta.data)
    ) {
      query_labels <- srt_query[[query_group, drop = TRUE]]
    } else if (length(query_group) == ncol(srt_query)) {
      query_labels <- query_group
    } else {
      log_message(
        "{.arg query_group} must be one metadata column or one value per query cell",
        message_type = "error"
      )
    }
    query_labels <- as.character(query_labels)
  } else {
    query_labels <- NULL
  }

  query_features <- rownames(srt_query[[query_assay]])
  ref_features <- rownames(srt_ref[[ref_assay]])
  features <- features %||% query_features
  common_features <- unique(intersect(
    features,
    intersect(query_features, ref_features)
  ))
  if (length(common_features) == 0L) {
    log_message(
      "No common features are available for {.fn RunSciBet}",
      message_type = "error"
    )
  }
  log_message(
    "Run native SciBet with {.val {length(common_features)}} candidate features and {.val {nlevels(ref_labels)}} reference classes",
    verbose = verbose
  )

  ref_expr <- scibet_expr_matrix(
    srt = srt_ref,
    assay = ref_assay,
    layer = ref_layer,
    features = common_features,
    input_transform = input_transform
  )
  query_expr <- scibet_expr_matrix(
    srt = srt_query,
    assay = query_assay,
    layer = query_layer,
    features = common_features,
    input_transform = input_transform
  )

  query_expr_run <- query_expr
  query_run_index <- colnames(query_expr)
  if (!is.null(query_labels)) {
    keep_query <- !is.na(query_labels)
    if (!all(keep_query)) {
      log_message(
        "Cells with missing {.arg query_group} will receive {.val NA} SciBet annotation",
        message_type = "warning",
        verbose = verbose
      )
    }
    group <- factor(query_labels[keep_query])
    if (nlevels(group) == 0L) {
      log_message(
        "{.arg query_group} does not contain any non-missing groups",
        message_type = "error"
      )
    }
    design <- stats::model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    query_expr_run <- query_expr[, keep_query, drop = FALSE] %*% design
    query_expr_run <- sweep(query_expr_run, 2, as.numeric(table(group)), "/")
    query_expr_run <- as.matrix(query_expr_run)
    storage.mode(query_expr_run) <- "double"
    query_run_index <- colnames(query_expr_run)
  }

  result <- scibet_fit_predict(
    ref = ref_expr,
    query = query_expr_run,
    labels = as.integer(ref_labels),
    n_labels = nlevels(ref_labels),
    n_top = as.integer(nfeatures),
    additional_per_label = as.integer(additional_features_per_class)
  )

  classes <- levels(ref_labels)
  prob_run <- result[["probabilities"]]
  colnames(prob_run) <- classes
  rownames(prob_run) <- query_run_index
  predicted_run <- classes[result[["predicted_index"]]]
  names(predicted_run) <- query_run_index
  score_run <- apply(prob_run, 1, max)

  if (!is.null(query_labels)) {
    predicted <- rep(NA_character_, ncol(srt_query))
    names(predicted) <- colnames(srt_query)
    score <- rep(NA_real_, ncol(srt_query))
    names(score) <- colnames(srt_query)
    matched <- as.character(query_labels) %in% names(predicted_run)
    predicted[matched] <- predicted_run[as.character(query_labels[matched])]
    score[matched] <- score_run[as.character(query_labels[matched])]
    prob_store <- prob_run
  } else {
    predicted <- predicted_run[colnames(srt_query)]
    score <- score_run[colnames(srt_query)]
    prob_store <- prob_run[colnames(srt_query), , drop = FALSE]
  }

  annotation_col <- paste0(prefix, "_annotation")
  score_col <- paste0(prefix, "_score")
  srt_query[[annotation_col]] <- unname(predicted)
  srt_query[[score_col]] <- unname(score)

  core <- result[["core"]]
  feature_index <- result[["feature_index"]]
  selected_features <- common_features[feature_index]
  rownames(core) <- selected_features
  colnames(core) <- classes

  if (isTRUE(store_model)) {
    srt_query@tools[["SciBet"]] <- list(
      probabilities = prob_store,
      model = core,
      features = selected_features,
      candidate_features = common_features,
      classes = classes,
      parameters = list(
        query_assay = query_assay,
        ref_assay = ref_assay,
        query_layer = query_layer,
        ref_layer = ref_layer,
        query_group = query_group,
        ref_group = ref_group,
        nfeatures = nfeatures,
        additional_features_per_class = additional_features_per_class,
        input_transform = input_transform,
        prefix = prefix
      )
    )
  }

  log_message(
    "{.pkg SciBet} annotations stored in metadata column {.val {annotation_col}}",
    verbose = verbose
  )
  srt_query
}

scibet_expr_matrix <- function(
  srt,
  assay,
  layer,
  features,
  input_transform = c("auto", "none", "expm1")
) {
  input_transform <- match.arg(input_transform)
  x <- GetAssayData5(
    object = srt,
    assay = assay,
    layer = layer
  )[features, , drop = FALSE]
  x <- as.matrix(x)
  storage.mode(x) <- "double"

  transform_use <- input_transform
  if (transform_use == "auto") {
    transform_use <- if (identical(layer, "data")) "expm1" else "none"
  }
  if (transform_use == "expm1") {
    x <- expm1(x)
  }
  x[!is.finite(x) | x < 0] <- 0
  x
}
