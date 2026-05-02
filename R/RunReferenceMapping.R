#' @title Map query cells into a reference space
#'
#' @description
#' Reference mapping workflow modeled after Seurat's `MapQuery` pattern. The
#' current implementation computes gene activity for ATAC query cells, finds
#' anchors to an RNA reference, integrates the query into the reference embedding
#' space, transfers labels, and projects query cells to the reference UMAP.
#'
#' @md
#' @inheritParams RunLabelTransfer
#' @param reference RNA reference `Seurat` object used for mapping.
#' @param prefix Prefix used to resolve ATAC reductions. Default is `"ATAC"`.
#' @param reference_assay Assay used in the reference object.
#' @param reference_reduction Reduction used in the reference object.
#' @param reference_dims Dimensions used from the reference reduction.
#' @param reference_label Metadata column in the reference used as transfer labels.
#' @param add_gene_activity Whether to calculate a gene activity assay for the query.
#' @param gene_activity_assay Name of the gene activity assay used for mapping.
#' @param ref_umap UMAP reduction in the RNA reference used for query projection.
#' @param label_method Label-transfer backend used after anchor mapping. One of
#' `"Seurat"` or `"scOMM"`.
#' @param reduction_project_method Anchor projection method. Default is
#' `"pcaproject"` for RNA reference mapping.
#' @param k.anchor,k.filter,k.score,k.weight Parameters passed to Seurat anchor
#' finding/integration.
#' @param projection_method,nn_method,k,distance_metric,vote_fun Passed to
#' [RunKNNMap] for UMAP projection and neighbor voting.
#'
#' @return A query `Seurat` object mapped into the RNA reference space.
#' @export
#' @examples
#' \dontrun{
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
#' query <- RunReferenceMapping(
#'   srt = query,
#'   reference = reference,
#'   assay = "peaks",
#'   reference_assay = "RNA",
#'   reference_reduction = "Standardpca",
#'   ref_umap = "StandardUMAP2D",
#'   reference_label = "CellType",
#'   reference_dims = 1:10,
#'   dims = 2:10
#' )
#' }
RunReferenceMapping <- function(
  srt,
  reference,
  assay = NULL,
  prefix = "ATAC",
  reference_assay = NULL,
  reference_reduction = "pca",
  ref_umap = NULL,
  reference_dims = 1:30,
  reference_label = NULL,
  label_method = c("Seurat", "scOMM"),
  add_gene_activity = TRUE,
  gene_activity_assay = "ACTIVITY",
  dims = 2:30,
  features = NULL,
  reduction_project_method = "pcaproject",
  k.anchor = 5,
  k.filter = 200,
  k.score = 30,
  k.weight = 100,
  projection_method = c("model", "knn"),
  nn_method = NULL,
  k = 30,
  distance_metric = "cosine",
  vote_fun = "mean",
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
  label_method <- match.arg(label_method)
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  reference_assay <- reference_assay %||% SeuratObject::DefaultAssay(reference)
  tool_name <- tool_name %||% paste0(prefix, "_ReferenceMapping")
  if (!inherits(srt[[assay]], "ChromatinAssay")) {
    log_message(
      "{.arg assay} must refer to a {.cls ChromatinAssay}",
      message_type = "error"
    )
  }
  if (!inherits(reference, "Seurat")) {
    log_message(
      "{.arg reference} is not a {.cls Seurat}",
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
  features <- features %||% SeuratObject::VariableFeatures(reference, assay = reference_assay)
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
      "No shared features available for ATAC reference mapping",
      message_type = "error"
    )
  }

  if (is.null(ref_umap)) {
    ref_umap <- sort(
      SeuratObject::Reductions(reference)[grep(
        "umap",
        SeuratObject::Reductions(reference),
        ignore.case = TRUE
      )]
    )[1]
  }
  if (is.na(ref_umap) || length(ref_umap) == 0) {
    log_message(
      "{.arg ref_umap} must refer to a UMAP reduction in {.arg reference}",
      message_type = "error"
    )
  }
  dims_use <- dims
  ref_dims_use <- reference_dims
  ndims_use <- min(length(dims_use), length(ref_dims_use))
  dims_use <- dims_use[seq_len(ndims_use)]
  ref_dims_use <- ref_dims_use[seq_len(ndims_use)]
  anchor_params <- atac_anchor_params(
    reference_cells = ncol(reference),
    query_cells = ncol(srt),
    k.anchor = k.anchor,
    k.filter = k.filter,
    k.score = k.score,
    verbose = verbose
  )

  log_message("Finding RNA-to-ATAC anchors for query mapping...", verbose = verbose)
  anchors <- Seurat::FindTransferAnchors(
    reference = reference,
    query = srt,
    reference.assay = reference_assay,
    query.assay = query_assay,
    reduction = reduction_project_method,
    reference.reduction = reference_reduction,
    dims = ref_dims_use,
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

  srt <- Seurat::IntegrateEmbeddings(
    anchorset = anchors,
    reference = reference,
    query = srt,
    reductions = reduction_project_method,
    new.reduction.name = "ref.pca",
    weight.reduction = if (identical(reduction_project_method, "pcaproject")) {
      "pcaproject"
    } else {
      reduction_project_method
    },
    dims.to.integrate = dims_use,
    k.weight = k.weight_use
  )
  query_dims_use <- intersect(
    dims_use,
    seq_len(ncol(Seurat::Embeddings(srt, reduction = "ref.pca")))
  )
  ref_dims_use <- ref_dims_use[seq_len(min(length(ref_dims_use), length(query_dims_use)))]
  query_dims_use <- query_dims_use[seq_len(min(length(ref_dims_use), length(query_dims_use)))]

  predicted_col <- NULL
  probability_col <- NULL
  probability_cols <- character(0)
  used_features <- features
  if (!is.null(reference_label)) {
    if (identical(label_method, "Seurat")) {
      srt <- atac_transfer_labels(
        srt = srt,
        assay = assay,
        prefix = prefix,
        linear_reduction_dims_use = query_dims_use,
        reference = reference,
        reference_assay = reference_assay,
        reference_reduction = reference_reduction,
        reference_dims = ref_dims_use,
        reference_label = reference_label,
        add_gene_activity = FALSE,
        gene_activity_assay = gene_activity_assay,
        weight_reduction = "ref.pca",
        features = features,
        prediction_prefix = "predicted_",
        k.weight = k.weight_use,
        verbose = verbose
      )
      predicted_col <- resolve_pred_col(
        srt = srt,
        prediction_prefix = "predicted_"
      )
      probability_col <- resolve_prob_col(
        srt = srt,
        prediction_prefix = "predicted_"
      )
      probability_cols <- grep(
        "^predicted_prediction\\.score\\.",
        colnames(srt@meta.data),
        value = TRUE
      )
    } else {
      srt <- RunscOMM(
        srt = srt,
        reference = reference,
        reference_assay = reference_assay,
        query_assay = query_assay,
        reference_label = reference_label,
        features = features,
        prediction_prefix = "predicted_",
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
  }

  srt <- RunKNNMap(
    srt_query = srt,
    query_assay = query_assay,
    srt_ref = reference,
    ref_assay = reference_assay,
    ref_group = reference_label,
    ref_umap = ref_umap,
    query_reduction = "ref.pca",
    ref_reduction = reference_reduction,
    query_dims = query_dims_use,
    ref_dims = ref_dims_use,
    projection_method = projection_method,
    nn_method = nn_method,
    k = k,
    distance_metric = distance_metric,
    vote_fun = vote_fun
  )
  srt@misc[["scop_atac_mapquery"]] <- list(
    query_assay = assay,
    gene_activity_assay = query_assay,
    reference_reduction = reference_reduction,
    ref_umap = ref_umap
  )
  if (isTRUE(evaluate)) {
    if (is.null(predicted_col)) {
      log_message(
        "{.arg reference_label} must be provided to evaluate reference mapping results.",
        message_type = "error"
      )
    }
    if (is.null(truth_col) || !truth_col %in% colnames(srt@meta.data)) {
      log_message(
        "{.arg truth_col} must be present in {.arg srt@meta.data} when {.arg evaluate = TRUE}",
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
    method = label_method,
    predicted_col = predicted_col,
    probability_col = probability_col,
    probability_cols = probability_cols,
    truth_col = truth_col,
    metrics = eval_res,
    reference_assay = reference_assay,
    query_assay = query_assay,
    reference_reduction = reference_reduction,
    ref_umap = ref_umap,
    features = used_features
  )
  srt
}
