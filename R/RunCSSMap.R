#' @title Single-cell reference mapping with CSS method
#'
#' @inheritParams RunKNNMap
#' @param ref_css The name of the CSS reduction in the reference object to use for calculating the distance metric.
#'
#' @seealso
#' [RunKNNMap]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(panc8_sub)
#' panc8_sub <- standard_scop(panc8_sub)
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
  check_r("quadbiolab/simspec", verbose = FALSE)
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
        "Length of {.arg ref_group} must be one or length of {.arg srt_ref}",
        message_type = "error"
      )
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
      log_message(
        "Cannot find CSS reduction in the {.arg srt_ref}",
        message_type = "error"
      )
    } else {
      log_message("Set {.arg ref_css} to {.val {ref_css}}")
      if (
        !"model" %in% names(srt_ref[[ref_css]]@misc) ||
          !"sim2profiles" %in% names(srt_ref[[ref_css]]@misc$model)
      ) {
        log_message(
          "CSS model is not in the reduction: {.val {ref_css}}",
          message_type = "error"
        )
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
      log_message(
        "Cannot find UMAP reduction in the {.arg srt_ref}",
        message_type = "error"
      )
    } else {
      log_message("Set {.arg ref_umap} to {.val {ref_umap}}")
    }
  }
  projection_method <- match.arg(projection_method)
  if (
    projection_method == "model" &&
      !"model" %in% names(srt_ref[[ref_umap]]@misc)
  ) {
    log_message(
      "No UMAP model detected. Set the {.arg projection_method} to {.val knn}",
      message_type = "warning"
    )
    projection_method <- "knn"
  }
  if (
    projection_method == "model" &&
      !distance_metric %in% c("euclidean", "cosine", "manhattan", "hamming")
  ) {
    log_message(
      "{.arg distance_metric} must be one of {.val euclidean}, {.val cosine}, {.val manhattan}, and {.val hamming} when {.arg projection_method='model'}",
      message_type = "warning"
    )
  }

  ref_assay <- srt_ref[[ref_css]]@assay.used
  status_query <- CheckDataType(
    GetAssayData5(srt_query, layer = "data", assay = query_assay)
  )
  log_message("Detected {.arg srt_query} data type: {.val {status_query}}")
  status_ref <- CheckDataType(
    GetAssayData5(srt_ref, layer = "data", assay = ref_assay)
  )
  log_message("Detected {.arg srt_ref} data type: {.val {status_ref}}")
  if (
    status_ref != status_query ||
      any(status_query == "unknown", status_ref == "unknown")
  ) {
    log_message(
      "Data type is unknown or different between {.arg srt_query} and {.arg srt_ref}",
      message_type = "warning"
    )
  }

  log_message("Run {.pkg CSS} projection")
  CSSmodel <- srt_ref[[ref_css]]@misc$model
  raw_assay <- SeuratObject::DefaultAssay(srt_query)
  SeuratObject::DefaultAssay(srt_query) <- query_assay
  srt_query <- invoke_fun(
    .fn = get("css_project", envir = getNamespace("simspec")),
    .args = list(object = srt_query, model = CSSmodel)
  )
  SeuratObject::DefaultAssay(srt_query) <- raw_assay

  log_message("Run {.pkg UMAP} projection")
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
