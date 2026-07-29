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
