#' @title Run a named gene regulatory network method
#' @md
#' @description
#' Each function selects its method in [RunGRN()]. Pass the remaining arguments
#' as you would to [RunGRN()].
#' @inheritParams RunGRN
#' @inherit RunGRN return
#' @details
#' PIDC and SCODE require `resources$source_dir`. For SCODE or SINCERITIES
#' with several trajectories, select one with `method_args$trajectory`.
#' SINCERITIES needs at least five distinct time points.
#'
#' Pando and multiome multiCSN need paired RNA and ATAC data and
#' `resources$genome`. LINGER needs paired count layers and a bulk reference
#' directory in `resources$bulk_reference`. Prepare the reference with
#' [PrepareDB()] before running LINGER.
#' @references
#' \url{https://github.com/CABSEL/SINCERITIES}
#' @seealso [RunGRN()], [PrepareEnv()], [GetGRN()], [GRNPlot()]
#' @rdname RunGRN-methods
#' @export
#' @examples
#' set.seed(41)
#' cells <- paste0("c", seq_len(80))
#' expression <- matrix(stats::runif(640), 8, dimnames = list(letters[1:8], cells))
#' pseudotime <- stats::setNames(rep(seq_len(8), each = 10), cells)
#' head(RunSINCERITIES(expression, pseudotime, verbose = FALSE))
RunSINCERITIES <- function(object, pseudotime, ...) {
  RunGRN(object, methods = "sincerities", pseudotime = pseudotime, ...)
}

#' @rdname RunGRN-methods
#' @export
RunPPCOR <- function(object, ...) {
  RunGRN(object, methods = "ppcor", ...)
}

#' @rdname RunGRN-methods
#' @export
RunPIDC <- function(object, ...) {
  RunGRN(object, methods = "pidc", ...)
}

#' @rdname RunGRN-methods
#' @export
RunSCORPION <- function(object, ...) {
  RunGRN(object, methods = "scorpion", ...)
}

#' @rdname RunGRN-methods
#' @export
RunLEAP <- function(object, pseudotime, ...) {
  RunGRN(object, methods = "leap", pseudotime = pseudotime, ...)
}

#' @rdname RunGRN-methods
#' @export
RunSCODE <- function(object, pseudotime, ...) {
  RunGRN(object, methods = "scode", pseudotime = pseudotime, ...)
}

#' @rdname RunGRN-methods
#' @export
RunMultiCSN <- function(object, ...) {
  RunGRN(object, methods = "multicsn", ...)
}

#' @rdname RunGRN-methods
#' @export
RunPando <- function(object, ...) {
  RunGRN(object, methods = "pando", ...)
}

#' @rdname RunGRN-methods
#' @export
RunLinger <- function(object, ...) {
  RunGRN(object, methods = "linger", ...)
}
