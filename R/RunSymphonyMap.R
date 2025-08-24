#' @title Single-cell reference mapping with Symphony method
#'
#' @md
#' @inheritParams RunKNNMap
#' @param ref_pca A character string specifying the name of the PCA reduction in the reference object to use for calculating the distance metric.
#' @param ref_harmony A character string specifying the name of the Harmony reduction in the reference object to use for calculating the distance metric.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(panc8_sub)
#' srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
#' srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
#' srt_ref <- integration_scop(
#'   srt_ref,
#'   batch = "tech",
#'   integration_method = "Harmony"
#' )
#' CellDimPlot(srt_ref, group.by = c("celltype", "tech"))
#'
#' # Projection
#' srt_query <- RunSymphonyMap(
#'   srt_query = srt_query,
#'   srt_ref = srt_ref,
#'   ref_pca = "Harmonypca",
#'   ref_harmony = "Harmony",
#'   ref_umap = "HarmonyUMAP2D"
#' )
#' ProjectionPlot(
#'   srt_query = srt_query,
#'   srt_ref = srt_ref,
#'   query_group = "celltype",
#'   ref_group = "celltype"
#' )
#' }
RunSymphonyMap <- function(
    srt_query,
    srt_ref,
    query_assay = NULL,
    ref_assay = srt_ref[[ref_pca]]@assay.used,
    ref_pca = NULL,
    ref_harmony = NULL,
    ref_umap = NULL,
    ref_group = NULL,
    projection_method = c("model", "knn"),
    nn_method = NULL,
    k = 30,
    distance_metric = "cosine",
    vote_fun = "mean") {
  check_r("immunogenomics/symphony")
  query_assay <- query_assay %||% SeuratObject::DefaultAssay(srt_query)
  ref_assay <- ref_assay %||% SeuratObject::DefaultAssay(srt_ref)
  if (!is.null(ref_group)) {
    if (length(ref_group) == ncol(srt_ref)) {
      srt_ref[["ref_group"]] <- ref_group
    } else if (length(ref_group) == 1) {
      if (!ref_group %in% colnames(srt_ref@meta.data)) {
        log_message(
          "{.arg ref_group} must be one of the column names in the meta.data",
          message_type = "error"
        )
      } else {
        srt_ref[["ref_group"]] <- srt_ref[[ref_group]]
      }
    } else {
      log_message(
        "Length of {.arg ref_group} must be one or length of {.arg srt_ref}.",
        message_type = "error"
      )
    }
    ref_group <- "ref_group"
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
        "Cannot find PCA reduction in the {.arg srt_ref}",
        message_type = "error"
      )
    } else {
      log_message("Set ref_pca to {.val {ref_pca}}")
    }
  }
  if (is.null(ref_harmony)) {
    ref_harmony <- sort(
      SeuratObject::Reductions(srt_ref)[grep(
        "harmony",
        SeuratObject::Reductions(srt_ref),
        ignore.case = TRUE
      )]
    )[1]
    if (length(ref_harmony) == 0) {
      log_message(
        "Cannot find Harmony reduction in the {.arg srt_ref}",
        message_type = "error"
      )
    } else {
      log_message("Set ref_harmony to {.val {ref_harmony}}")
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
        "Cannot find UMAP reduction in the {.arg srt_ref}",
        message_type = "error"
      )
    } else {
      log_message("Set ref_umap to {.val {ref_umap}}")
    }
  }
  ref_pca_dims <- srt_ref[[ref_harmony]]@misc$reduction_dims

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

  status_query <- check_data_type(
    data = GetAssayData5(
      srt_query,
      layer = "data",
      assay = query_assay
    )
  )
  log_message("Detected srt_query data type: {.val {status_query}}")
  status_ref <- check_data_type(
    data = GetAssayData5(
      srt_ref,
      layer = "data",
      assay = ref_assay
    )
  )
  log_message("Detected srt_ref data type: {.val {status_ref}}")
  if (
    status_ref != status_query ||
      any(status_query == "unknown", status_ref == "unknown")
  ) {
    log_message(
      "Data type is unknown or different between {.arg srt_query} and {.arg srt_ref}.",
      message_type = "warning"
    )
  }

  log_message("Build reference")
  ref <- buildReferenceFromSeurat(
    obj = srt_ref,
    assay = ref_assay,
    pca = ref_pca,
    pca_dims = ref_pca_dims,
    harmony = ref_harmony,
    umap = ref_umap
  )
  log_message("Run mapQuery")
  res <- mapQuery(
    exp_query = GetAssayData5(
      srt_query,
      layer = "data",
      assay = query_assay
    ),
    metadata_query = srt_query@meta.data,
    ref_obj = ref,
    vars = NULL,
    sigma = 0.1,
    verbose = TRUE
  )
  Z_pca_query <- res$Z_pca_query
  Zq_corr <- res$Zq_corr
  R_query <- res$R_query

  srt_query[["ref.pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = Matrix::t(Z_pca_query),
    loadings = ref$loadings,
    stdev = as.numeric(apply(Z_pca_query, 1, stats::sd)),
    assay = query_assay,
    key = "refpca_"
  )
  srt_query[["ref.harmony"]] <- SeuratObject::CreateDimReducObject(
    embeddings = Matrix::t(Zq_corr),
    stdev = as.numeric(apply(Zq_corr, 1, stats::sd)),
    assay = query_assay,
    key = "refharmony_",
    misc = list(R = R_query)
  )

  log_message("Run UMAP projection")
  ref_dims <- seq_len(dim(srt_ref[[ref_harmony]])[2])
  srt_query <- RunKNNMap(
    srt_query = srt_query,
    query_assay = query_assay,
    srt_ref = srt_ref,
    ref_assay = ref_assay,
    ref_group = ref_group,
    ref_umap = ref_umap,
    query_reduction = "ref.harmony",
    ref_reduction = ref_harmony,
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
