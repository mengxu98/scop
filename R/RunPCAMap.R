#' Single-cell reference mapping with PCA method
#'
#' @inheritParams RunKNNMap
#' @param ref_pca  A character string specifying the name of a PCA reduction in the reference object to use for calculating the distance metric. If NULL (default), it will be automatically detected as the first PCA reduction.
#'
#' @export
#'
#' @examples
#' data(panc8_sub)
#' srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
#' srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
#' srt_ref <- integration_scop(
#'   srt_ref,
#'   batch = "tech",
#'   integration_method = "Seurat"
#' )
#' CellDimPlot(srt_ref, group.by = c("celltype", "tech"))
#'
#' # Projection
#' srt_query <- RunPCAMap(
#'   srt_query = srt_query,
#'   srt_ref = srt_ref,
#'   ref_pca = "Seuratpca",
#'   ref_umap = "SeuratUMAP2D"
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
    vote_fun = "mean") {
  query_assay <- query_assay %||% SeuratObject::DefaultAssay(srt_query)
  ref_assay <- ref_assay %||% SeuratObject::DefaultAssay(srt_ref)
  if (!is.null(ref_group)) {
    if (length(ref_group) == ncol(srt_ref)) {
      srt_ref[["ref_group"]] <- ref_group
    } else if (length(ref_group) == 1) {
      if (!ref_group %in% colnames(srt_ref@meta.data)) {
        log_message(
          "ref_group must be one of the column names in the meta.data",
          message_type = "error"
        )
      } else {
        srt_ref[["ref_group"]] <- srt_ref[[ref_group]]
      }
    } else {
      log_message(
        "Length of ref_group must be one or length of srt_ref.",
        message_type = "error"
      )
    }
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
      log_message("Set ref_pca to ", ref_pca)
    }
  }
  if (is.null(ref_umap)) {
    ref_umap <- sort(
      SeuratObject::Reductions(srt_ref)[grep(
        "umap",
        SeuratObject::Reductions(srt_ref),
        ignore.case = TRUE
      )]
    )[1]
    if (length(ref_umap) == 0) {
      log_message(
        "Cannot find UMAP reduction in the srt_ref",
        message_type = "error"
      )
    } else {
      log_message("Set ref_umap to ", ref_umap)
    }
  }
  projection_method <- match.arg(projection_method)
  if (
    projection_method == "model" &&
      !"model" %in% names(srt_ref[[ref_umap]]@misc)
  ) {
    log_message(
      "No UMAP model detected. Set the projection_method to 'knn'",
      message_type = "warning"
    )
    projection_method <- "knn"
  }
  if (
    projection_method == "model" &&
      !distance_metric %in% c("euclidean", "cosine", "manhattan", "hamming")
  ) {
    log_message(
      "distance_metric must be one of euclidean, cosine, manhattan, and hamming when projection_method='model'",
      message_type = "error"
    )
  }

  pca.out <- srt_ref[[ref_pca]]
  status_query <- CheckDataType(
    object = GetAssayData5(
      srt_query,
      layer = "data",
      assay = query_assay
    )
  )
  log_message("Detected srt_query data type: ", status_query)
  status_ref <- CheckDataType(
    object = GetAssayData5(
      srt_ref,
      layer = "data",
      assay = ref_assay
    )
  )
  log_message("Detected srt_ref data type: ", status_ref)
  if (
    status_ref != status_query ||
      any(status_query == "unknown", status_ref == "unknown")
  ) {
    log_message(
      "Data type is unknown or different between srt_query and srt_ref.",
      message_type = "warning"
    )
  }

  log_message("Run PCA projection")
  features <- rownames(pca.out@feature.loadings)
  center <- apply(
    GetAssayData5(
      object = srt_ref,
      layer = "data",
      assay = ref_assay
    )[
      features,
    ],
    1,
    mean
  )
  names(center) <- features
  sds <- apply(
    GetAssayData5(
      object = srt_ref,
      layer = "data",
      assay = ref_assay
    )[
      features,
    ],
    1,
    stats::sd
  )
  names(sds) <- features
  rotation <- pca.out@feature.loadings

  features_common <- Reduce(
    intersect,
    list(
      features,
      rownames(srt_query[[query_assay]]),
      rownames(srt_ref[[ref_assay]])
    )
  )
  log_message("Use ", length(features_common), " features to calculate PC.")
  query_data <- Matrix::t(
    GetAssayData5(
      srt_query,
      layer = "data",
      assay = query_assay
    )[
      features_common,
    ]
  )
  query_pca <- scale(
    query_data[, features_common],
    center[features_common],
    sds[features_common]
  ) %*%
    rotation[features_common, ]

  srt_query[["ref.pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = query_pca,
    key = pca.out@key,
    assay = query_assay,
    verbose = FALSE
  )

  log_message("Run UMAP projection")
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
