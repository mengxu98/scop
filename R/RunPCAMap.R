#' @title Single-cell reference mapping with PCA method
#'
#' @inheritParams RunKNNMap
#' @inheritParams thisutils::log_message
#' @param ref_pca Name of a PCA reduction in the reference object to use for calculating the distance metric. If NULL (default), it will be automatically detected as the first PCA reduction.
#'
#' @export
#'
#' @examples
#' data(panc8_sub)
#' panc8_sub <- RunStandardWorkflow(panc8_sub)
#' srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
#' srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
#' srt_ref <- RunIntegration(
#'   srt_ref,
#'   batch = "tech",
#'   integration_method = "Uncorrected"
#' )
#' CellDimPlot(srt_ref, group.by = c("celltype", "tech"))
#'
#' # Projection
#' srt_query <- RunPCAMap(
#'   srt_query = srt_query,
#'   srt_ref = srt_ref,
#'   ref_pca = "Uncorrectedpca",
#'   ref_umap = "UncorrectedUMAP2D"
#' )
#' ProjectionPlot(
#'   srt_query = srt_query,
#'   srt_ref = srt_ref,
#'   query_group = "celltype",
#'   ref_group = "celltype"
#' )
RunPCAMap <- function(
  srt_query,
  srt_ref,
  query_assay = NULL,
  ref_assay = srt_ref[[ref_pca]]@assay.used,
  ref_pca = NULL,
  ref_dims = 1:30,
  ref_umap = NULL,
  ref_group = NULL,
  projection_method = c("model", "knn"),
  nn_method = NULL,
  k = 30,
  distance_metric = "cosine",
  vote_fun = "mean",
  verbose = TRUE
) {
  query_assay <- query_assay %||% SeuratObject::DefaultAssay(srt_query)
  ref_assay <- ref_assay %||% SeuratObject::DefaultAssay(srt_ref)
  if (!is.null(ref_group)) {
    srt_ref <- prepare_mapping_reference(
      srt_ref, ref_group,
      column_error = "ref_group must be one of the column names in the meta.data",
      length_error = "Length of ref_group must be one or length of srt_ref."
    )
  }

  if (is.null(ref_pca)) {
    ref_pca <- sort(
      SeuratObject::Reductions(srt_ref)[grep(
        "pca",
        SeuratObject::Reductions(srt_ref),
        ignore.case = TRUE
      )]
    )[1]
    if (length(ref_pca) == 0) {
      log_message(
        "Cannot find PCA reduction in the srt_ref",
        message_type = "error"
      )
    } else {
      log_message("Set ref_pca to ", ref_pca, verbose = verbose)
    }
  }
  if (is.null(ref_umap)) {
    ref_umap <- mapping_default_umap(srt_ref)
    log_message("Set ref_umap to ", ref_umap, verbose = verbose)
  }
  projection_method <- mapping_projection_method(
    srt_ref, ref_umap, projection_method, distance_metric, verbose
  )

  pca.out <- srt_ref[[ref_pca]]
  query_data <- GetAssayData5(srt_query, layer = "data", assay = query_assay)
  ref_data <- GetAssayData5(srt_ref, layer = "data", assay = ref_assay)
  status_query <- CheckDataType(object = query_data)
  log_message("Detected srt_query data type: ", status_query, verbose = verbose)
  status_ref <- CheckDataType(object = ref_data)
  log_message("Detected srt_ref data type: ", status_ref, verbose = verbose)
  if (
    status_ref != status_query ||
      any(status_query == "unknown", status_ref == "unknown")
  ) {
    log_message(
      "Data type is unknown or different between srt_query and srt_ref.",
      message_type = "warning",
      verbose = verbose
    )
  }

  log_message("Run PCA projection", verbose = verbose)
  rotation <- pca.out@feature.loadings
  features_common <- Reduce(
    intersect,
    list(rownames(rotation), rownames(query_data), rownames(ref_data))
  )
  if (!length(features_common)) {
    log_message(
      "No shared PCA features in the query and reference data layers.",
      message_type = "error"
    )
  }
  ref_pca_data <- ref_data[features_common, , drop = FALSE]
  center <- Matrix::rowMeans(ref_pca_data)
  sds <- sqrt(fast_row_vars(ref_pca_data))
  if (length(features_common) == 1L && (!is.finite(sds) || sds <= 0)) {
    log_message(
      "The sole shared PCA feature must have positive finite reference variance.",
      message_type = "error"
    )
  }

  log_message("Use ", length(features_common), " features to calculate PC.", verbose = verbose)
  query_data <- Matrix::t(query_data[features_common, , drop = FALSE])
  query_pca <- scale(query_data, center, sds) %*%
    rotation[features_common, , drop = FALSE]

  srt_query[["ref.pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = query_pca,
    key = pca.out@key,
    assay = query_assay
  )

  log_message("Run UMAP projection", verbose = verbose)
  srt_query <- RunKNNMap(
    srt_query = srt_query,
    query_assay = query_assay,
    srt_ref = srt_ref,
    ref_assay = ref_assay,
    ref_group = ref_group,
    ref_umap = ref_umap,
    query_reduction = "ref.pca",
    ref_reduction = ref_pca,
    query_dims = ref_dims,
    ref_dims = ref_dims,
    projection_method = projection_method,
    nn_method = nn_method,
    k = k,
    distance_metric = distance_metric,
    vote_fun = vote_fun
  )
  return(srt_query)
}
