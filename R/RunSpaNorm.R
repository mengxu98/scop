#' @title Run SpaNorm spatial normalization
#'
#' @description
#' Reserved wrapper entry for issue #257. The implementation will be added in a
#' follow-up PR after the package API is finalized.
#'
#' @md
#' @inheritParams RunSpatialVariableFeatures
#' @param new_assay Name planned for storing SpaNorm-normalized expression.
#' @param ... Reserved for arguments passed to the SpaNorm backend.
#'
#' @return This placeholder currently stops with an informative message.
#' @export
#'
#' @examples
#' \dontrun{
#' data(visium_human_pancreas_sub)
#' spatial <- RunSpaNorm(
#'   visium_human_pancreas_sub,
#'   assay = "Spatial",
#'   layer = "counts",
#'   new_assay = "SpaNorm"
#' )
#' }
RunSpaNorm <- function(
  srt,
  assay = NULL,
  layer = "counts",
  new_assay = "SpaNorm",
  ...
) {
  stop(
    "RunSpaNorm() is reserved for issue #257 and will be implemented in a follow-up PR.",
    call. = FALSE
  )
}
