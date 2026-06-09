#' @title Run common cell-cell communication analyses
#'
#' @md
#' @inheritParams RunCellChat
#' @param methods Cell-cell communication methods to run. Currently supports
#' `"CellChat"`, `"CellphoneDB"`, and `"LIANA"` in the unified scheduler.
#' NicheNet and MultiNicheNet require explicit receiver/sender/contrast design
#' arguments and should be called through [RunNichenetr()] or
#' [RunMultiNichenetr()].
#' @param method_params Named list of method-specific arguments passed to the
#' corresponding wrapper. For example, use `method_params$CellphoneDB$pvalue`
#' for CellphoneDB-specific parameters.
#' @param backend Backend used for scop post-processing and unified CCC table
#' aggregation. The upstream CellChat, CellphoneDB, and LIANA inference logic is
#' unchanged.
#' @param skip_failed Whether to keep running remaining methods if one method
#' fails.
#' @param rebuild_unified Whether to rebuild `srt@tools[["CCC"]]` from the
#' completed methods after all requested methods finish.
#' @param thresh Significance threshold used when rebuilding unified CCC tables
#' and passed to `RunCellChat()` unless overridden in `method_params$CellChat`.
#'
#' @return A `Seurat` object with method-specific results and a unified
#' `srt@tools[["CCC"]]` bundle.
#' @export
RunCCC <- function(
  srt,
  group.by,
  methods = c("CellChat", "CellphoneDB", "LIANA"),
  method_params = list(),
  backend = c("cpp", "r"),
  skip_failed = FALSE,
  rebuild_unified = TRUE,
  thresh = 0.05,
  verbose = TRUE
) {
  backend <- match.arg(backend)
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (!is.character(group.by) || length(group.by) != 1L || !group.by %in% colnames(srt[[]])) {
    log_message(
      "{.arg group.by} must be a valid metadata column in {.cls Seurat}",
      message_type = "error"
    )
  }

  methods <- unique(vapply(methods, normalize_ccc_method, character(1)))
  supported <- c("CellChat", "CellphoneDB", "LIANA")
  unsupported <- setdiff(methods, supported)
  if (length(unsupported) > 0L) {
    log_message(
      paste0(
        "{.fn RunCCC} currently supports {.val CellChat}, ",
        "{.val CellphoneDB}, and {.val LIANA}. Call ",
        "{.fn RunNichenetr} or {.fn RunMultiNichenetr} directly for methods ",
        "that require experiment-specific receiver/sender/contrast design: ",
        paste(unsupported, collapse = ", ")
      ),
      message_type = "error"
    )
  }

  method_params <- ccc_normalize_run_params(method_params)
  status <- list()
  completed_methods <- character(0)
  started_at <- Sys.time()

  for (method in methods) {
    log_message(
      "Running CCC method: {.val {method}}",
      verbose = verbose
    )
    args <- ccc_run_method_args(
      method = method,
      srt = srt,
      group.by = group.by,
      backend = backend,
      verbose = verbose,
      thresh = thresh,
      params = method_params[[method]] %||% list()
    )
    fun <- switch(
      method,
      CellChat = RunCellChat,
      CellphoneDB = RunCellphoneDB,
      LIANA = RunLIANA
    )

    start <- proc.time()[["elapsed"]]
    result <- tryCatch(
      do.call(fun, args),
      error = function(e) e
    )
    elapsed <- proc.time()[["elapsed"]] - start

    if (inherits(result, "error")) {
      status[[method]] <- list(
        method = method,
        status = "failed",
        elapsed = elapsed,
        message = conditionMessage(result)
      )
      if (!isTRUE(skip_failed)) {
        stop(result)
      }
      log_message(
        "CCC method {.val {method}} failed and was skipped: {conditionMessage(result)}",
        message_type = "warning",
        verbose = verbose
      )
      next
    }

    srt <- result
    completed_methods <- c(completed_methods, method)
    status[[method]] <- list(
      method = method,
      status = "completed",
      elapsed = elapsed,
      message = NA_character_
    )
  }

  if (isTRUE(rebuild_unified) && length(completed_methods) > 0L) {
    srt@tools[["CCC"]] <- ccc_build_unified_bundle(
      srt = srt,
      methods = completed_methods,
      thresh = thresh,
      backend = backend
    )
  }

  status_df <- ccc_run_status_df(status)
  run_record <- list(
    method = "RunCCC",
    methods = methods,
    completed_methods = completed_methods,
    failed_methods = setdiff(methods, completed_methods),
    status = status_df,
    parameters = list(
      group.by = group.by,
      methods = methods,
      method_params = method_params,
      backend = backend,
      skip_failed = skip_failed,
      rebuild_unified = rebuild_unified,
      thresh = thresh
    ),
    started_at = as.character(started_at),
    updated_at = as.character(Sys.time())
  )
  srt@tools[["RunCCC"]] <- run_record
  if (!is.null(srt@tools[["CCC"]])) {
    srt@tools[["CCC"]]$metadata$runccc <- run_record
  }

  log_message(
    "CCC analysis completed for methods: {.val {completed_methods}}",
    message_type = "success",
    verbose = verbose
  )
  srt
}

ccc_normalize_run_params <- function(method_params) {
  if (is.null(method_params)) {
    return(list())
  }
  if (!is.list(method_params) || is.data.frame(method_params)) {
    log_message(
      "{.arg method_params} must be a named list",
      message_type = "error"
    )
  }
  if (length(method_params) == 0L) {
    return(list())
  }
  nms <- names(method_params)
  if (is.null(nms) || any(!nzchar(nms))) {
    log_message(
      "{.arg method_params} must be named by CCC method",
      message_type = "error"
    )
  }
  out <- list()
  for (nm in nms) {
    method <- normalize_ccc_method(nm)
    params <- method_params[[nm]] %||% list()
    if (!is.list(params) || is.data.frame(params)) {
      log_message(
        "Each {.arg method_params} entry must be a list",
        message_type = "error"
      )
    }
    out[[method]] <- utils::modifyList(out[[method]] %||% list(), params, keep.null = TRUE)
  }
  out
}

ccc_run_method_args <- function(
  method,
  srt,
  group.by,
  backend,
  verbose,
  thresh,
  params = list()
) {
  protected <- intersect(c("srt", "group.by", "backend"), names(params))
  if (length(protected) > 0L) {
    log_message(
      "Ignoring protected {.arg method_params} entries for {.val {method}}: {.val {protected}}",
      message_type = "warning",
      verbose = verbose
    )
    params[protected] <- NULL
  }

  base <- list(
    srt = srt,
    group.by = group.by,
    backend = backend,
    verbose = verbose
  )
  if (identical(method, "CellChat") && !"thresh" %in% names(params)) {
    base$thresh <- thresh
  }
  utils::modifyList(base, params, keep.null = TRUE)
}

ccc_run_status_df <- function(status) {
  if (length(status) == 0L) {
    return(data.frame())
  }
  rows <- lapply(status, function(x) {
    data.frame(
      method = as.character(x$method %||% NA_character_),
      status = as.character(x$status %||% NA_character_),
      elapsed = as.numeric(x$elapsed %||% NA_real_),
      message = as.character(x$message %||% NA_character_),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
