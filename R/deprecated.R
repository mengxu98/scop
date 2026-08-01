#' Deprecated workflow entry points
#'
#' `standard_scop()` and `integration_scop()` were renamed to
#' [RunStandardWorkflow()] and [RunIntegration()]. The compatibility entry
#' points remain available with a warning in releases before 1.0.0 and will be
#' removed in version 1.0.0.
#'
#' @param ... Arguments forwarded unchanged to the replacement function.
#'
#' @return The value returned by the replacement function.
#' @name deprecated-workflows
#' @keywords internal
NULL

#' @rdname deprecated-workflows
#' @export
standard_scop <- function(...) {
  .Deprecated(
    new = "RunStandardWorkflow",
    package = "scop",
    msg = paste0(
      "`standard_scop()` is deprecated; use `RunStandardWorkflow()` instead. ",
      "It will be removed in scop 1.0.0."
    )
  )
  RunStandardWorkflow(...)
}

#' @rdname deprecated-workflows
#' @export
integration_scop <- function(...) {
  .Deprecated(
    new = "RunIntegration",
    package = "scop",
    msg = paste0(
      "`integration_scop()` is deprecated; use `RunIntegration()` instead. ",
      "It will be removed in scop 1.0.0."
    )
  )
  RunIntegration(...)
}

#' Deprecated doublet-calling entry points
#'
#' `db_scDblFinder()`, `db_scds()`, `db_Scrublet()`, and
#' `db_DoubletDetection()` were renamed to [RunscDblFinder()], [Runscds()],
#' [RunScrublet()], and [RunDoubletDetection()]. The compatibility entry points
#' remain available with a warning and will be removed in version 1.0.0.
#'
#' @param ... Arguments forwarded unchanged to the replacement function.
#'
#' @return The value returned by the replacement function.
#' @name deprecated-doublet-callers
#' @keywords internal
NULL

#' @rdname deprecated-doublet-callers
#' @export
db_scDblFinder <- function(...) {
  .Deprecated(
    new = "RunscDblFinder",
    package = "scop",
    msg = paste0(
      "`db_scDblFinder()` is deprecated; use `RunscDblFinder()` instead. ",
      "It will be removed in scop 1.0.0."
    )
  )
  RunscDblFinder(...)
}

#' @rdname deprecated-doublet-callers
#' @export
db_scds <- function(...) {
  .Deprecated(
    new = "Runscds",
    package = "scop",
    msg = paste0(
      "`db_scds()` is deprecated; use `Runscds()` instead. ",
      "It will be removed in scop 1.0.0."
    )
  )
  Runscds(...)
}

#' @rdname deprecated-doublet-callers
#' @export
db_Scrublet <- function(...) {
  .Deprecated(
    new = "RunScrublet",
    package = "scop",
    msg = paste0(
      "`db_Scrublet()` is deprecated; use `RunScrublet()` instead. ",
      "It will be removed in scop 1.0.0."
    )
  )
  RunScrublet(...)
}

#' @rdname deprecated-doublet-callers
#' @export
db_DoubletDetection <- function(...) {
  .Deprecated(
    new = "RunDoubletDetection",
    package = "scop",
    msg = paste0(
      "`db_DoubletDetection()` is deprecated; use `RunDoubletDetection()` instead. ",
      "It will be removed in scop 1.0.0."
    )
  )
  RunDoubletDetection(...)
}
