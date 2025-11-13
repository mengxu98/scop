#' Single-cell reference mapping with Seurat method
#'
#' @inheritParams RunKNNMap
#' @param ref_pca A character string specifying the name of the PCA reduction in the reference object to use for calculating the distance metric.
#' @param normalization.method The normalization method to use.
#' Default is `"LogNormalize"`.
#' @param reduction_project_method Dimensional reduction to perform when finding anchors.
#' Default is `"pcaproject"`.
#' @param k.anchor How many neighbors (k) to use when finding anchors.
#' Default is `5`.
#' @param k.filter How many neighbors (k) to use when filtering anchors. Set to NA to turn off filtering.
#' Default is `200`.
#' @param k.score How many neighbors (k) to use when scoring anchors.
#' Default is `30`.
#' @param k.weight Number of neighbors to consider when weighting anchors.
#' Default is `100`.
#'
#' @export
#'
#' @examples
#' data(panc8_sub)
#' panc8_sub <- standard_scop(panc8_sub)
#' srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
#' srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
#' srt_ref <- integration_scop(
#'   srt_ref,
#'   batch = "tech",
#'   integration_method = "Uncorrected"
#' )
#' CellDimPlot(srt_ref, group.by = c("celltype", "tech"))
#'
#' # Projection
#' srt_query <- RunSeuratMap(
#'   srt_query = srt_query,
#'   srt_ref = srt_ref,
#'   ref_pca = "Uncorrectedpca",
#'   ref_umap = "UncorrectedUMAP2D",
#'   k.weight = 50
#' )
#' ProjectionPlot(
#'   srt_query = srt_query,
#'   srt_ref = srt_ref,
#'   query_group = "celltype",
#'   ref_group = "celltype"
#' )
RunSeuratMap <- function(
    srt_query,
    srt_ref,
    query_assay = NULL,
    ref_pca = NULL,
    ref_assay = srt_ref[[ref_pca]]@assay.used,
    ref_dims = 1:30,
    ref_umap = NULL,
    ref_group = NULL,
    normalization.method = "LogNormalize",
    reduction_project_method = "pcaproject",
    k.anchor = 5,
    k.filter = 200,
    k.score = 30,
    k.weight = 100,
    projection_method = c("model", "knn"),
    nn_method = NULL,
    k = 30,
    distance_metric = "cosine",
    vote_fun = "mean") {
  query_assay <- query_assay %||% SeuratObject::DefaultAssay(srt_query)
  ref_assay <- ref_assay %||% SeuratObject::DefaultAssay(srt_ref)
  weight_reduction <- switch(reduction_project_method,
    "pcaproject" = "pcaproject",
    "lsiproject" = "lsiproject",
    "rpca" = "pcaproject",
    "cca" = "cca"
  )
  if (is.null(ref_pca)) {
    if (any(grepl("pca", SeuratObject::Reductions(srt_ref), ignore.case = TRUE))) {
      ref_pca <- sort(
        SeuratObject::Reductions(srt_ref)[grep("pca", SeuratObject::Reductions(srt_ref))]
      )[1]
    } else {
      log_message(
        "'ref_pca' is NUll and no pca reduction detected. Run standard_scop first.\n"
      )
      srt_ref <- standard_scop(srt_ref)
      ref_pca <- "Standardpca"
    }
    log_message("Set the ref_pca to '", ref_pca, "'")
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

  status_query <- CheckDataType(
    GetAssayData5(
      srt_query,
      layer = "data",
      assay = query_assay
    )
  )
  log_message("Detected srt_query data type: ", status_query)
  status_ref <- CheckDataType(
    GetAssayData5(
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

  log_message("Run FindTransferAnchors")
  anchors <- Seurat::FindTransferAnchors(
    query = srt_query,
    query.assay = query_assay,
    reference = srt_ref,
    normalization.method = normalization.method,
    dims = ref_dims,
    reference.reduction = ref_pca,
    reduction = reduction_project_method,
    reference.assay = ref_assay,
    k.anchor = k.anchor,
    k.filter = k.filter,
    k.score = k.score
  )
  if (reduction_project_method != "cca") {
    srt_query <- Seurat::IntegrateEmbeddings(
      anchorset = anchors,
      reference = srt_ref,
      query = srt_query,
      reductions = reduction_project_method,
      new.reduction.name = "ref.pca",
      weight.reduction = weight_reduction,
      k.weight = k.weight
    )
  }

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
