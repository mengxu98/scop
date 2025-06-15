#' Single-cell reference mapping with CSS method
#' @inheritParams RunKNNMap
#' @param ref_css A character string specifying the name of the CSS reduction in the reference object to use for calculating the distance metric.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("panc8_sub")
#' srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
#' srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
#' srt_ref <- integration_scop(
#'   srt_ref,
#'   batch = "tech",
#'   integration_method = "CSS"
#' )
#' CellDimPlot(srt_ref, group.by = c("celltype", "tech"))
#'
#' # Projection
#' srt_query <- RunCSSMap(
#'   srt_query = srt_query,
#'   srt_ref = srt_ref,
#'   ref_css = "CSS",
#'   ref_umap = "CSSUMAP2D"
#' )
#' ProjectionPlot(
#'   srt_query = srt_query,
#'   srt_ref = srt_ref,
#'   query_group = "celltype",
#'   ref_group = "celltype"
#' )
#' }
RunCSSMap <- function(
    srt_query,
    srt_ref,
    query_assay = NULL,
    ref_assay = srt_ref[[ref_css]]@assay.used,
    ref_css = NULL,
    ref_umap = NULL,
    ref_group = NULL,
    projection_method = c("model", "knn"),
    nn_method = NULL,
    k = 30,
    distance_metric = "cosine",
    vote_fun = "mean") {
  check_r("quadbiolab/simspec")
  query_assay <- query_assay %||% SeuratObject::DefaultAssay(srt_query)
  ref_assay <- ref_assay %||% SeuratObject::DefaultAssay(srt_ref)
  if (!is.null(ref_group)) {
    if (length(ref_group) == ncol(srt_ref)) {
      srt_ref[["ref_group"]] <- ref_group
    } else if (length(ref_group) == 1) {
      if (!ref_group %in% colnames(srt_ref@meta.data)) {
        stop("ref_group must be one of the column names in the meta.data")
      } else {
        srt_ref[["ref_group"]] <- srt_ref[[ref_group]]
      }
    } else {
      stop("Length of ref_group must be one or length of srt_ref.")
    }
    ref_group <- "ref_group"
  }
  if (is.null(ref_css)) {
    ref_css <- sort(
      SeuratObject::Reductions(srt_ref)[grep(
        "css",
        SeuratObject::Reductions(srt_ref),
        ignore.case = TRUE
      )]
    )[1]
    if (length(ref_css) == 0) {
      stop("Cannot find CSS reduction in the srt_ref")
    } else {
      message("Set ref_css to ", ref_css)
      if (
        !"model" %in% names(srt_ref[[ref_css]]@misc) ||
          !"sim2profiles" %in% names(srt_ref[[ref_css]]@misc$model)
      ) {
        stop("CSS model is not in the reduction: ", ref_css)
      }
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
      stop("Cannot find UMAP reduction in the srt_ref")
    } else {
      message("Set ref_umap to ", ref_umap)
    }
  }
  projection_method <- match.arg(projection_method)
  if (
    projection_method == "model" &&
      !"model" %in% names(srt_ref[[ref_umap]]@misc)
  ) {
    message("No UMAP model detected. Set the projection_method to 'knn'")
    projection_method <- "knn"
  }
  if (
    projection_method == "model" &&
      !distance_metric %in% c("euclidean", "cosine", "manhattan", "hamming")
  ) {
    stop(
      "distance_metric must be one of euclidean, cosine, manhattan, and hamming when projection_method='model'"
    )
  }

  ref_assay <- srt_ref[[ref_css]]@assay.used
  status_query <- check_data_type(
    data = GetAssayData5(srt_query, layer = "data", assay = query_assay)
  )
  message("Detected srt_query data type: ", status_query)
  status_ref <- check_data_type(
    data = GetAssayData5(srt_ref, layer = "data", assay = ref_assay)
  )
  message("Detected srt_ref data type: ", status_ref)
  if (
    status_ref != status_query ||
      any(status_query == "unknown", status_ref == "unknown")
  ) {
    warning(
      "Data type is unknown or different between srt_query and srt_ref.",
      immediate. = TRUE
    )
  }

  message("Run CSS projection")
  CSSmodel <- srt_ref[[ref_css]]@misc$model
  raw_assay <- SeuratObject::DefaultAssay(srt_query)
  SeuratObject::DefaultAssay(srt_query) <- query_assay
  srt_query <- invoke_fun(
    .fn = get("css_project", envir = getNamespace("simspec")),
    .args = list(object = srt_query, model = CSSmodel)
  )
  SeuratObject::DefaultAssay(srt_query) <- raw_assay

  message("Run UMAP projection")
  ref_dims <- seq_len(dim(srt_ref[[ref_css]])[2])
  srt_query <- RunKNNMap(
    srt_query = srt_query,
    query_assay = query_assay,
    srt_ref = srt_ref,
    ref_assay = ref_assay,
    ref_group = ref_group,
    ref_umap = ref_umap,
    query_reduction = "css_proj",
    ref_reduction = ref_css,
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
