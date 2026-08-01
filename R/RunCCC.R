#' @title Run common cell-cell communication analyses
#'
#' @md
#' @inheritParams RunCellChat
#' @param methods Registered cell-cell communication methods to run. The
#' default core methods are `"CellChat"`, `"CellphoneDB"`, and `"LIANA"`.
#' NicheNet, MultiNicheNet, SpatialCellChat, and MDIC3 can be selected when
#' their design-specific arguments are supplied through `method_params`.
#' LIANA's default internal method set includes its CellPhoneDB scorer, so
#' standalone CellphoneDB and LIANA consensus results are not statistically
#' independent evidence.
#' @param method_params Named list of method-specific arguments passed to the
#' corresponding wrapper. For example, use `method_params$CellphoneDB$pvalue`
#' for CellphoneDB-specific parameters.
#' @param backend Backend used only for result post-processing and unified CCC table
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
  registry <- ccc_method_registry()
  supported <- names(registry)
  unsupported <- setdiff(methods, supported)
  if (length(unsupported) > 0L) {
    log_message(
      "Unsupported CCC methods: {.val {unsupported}}. Use {.fn ListCCCMethods} to inspect registered methods.",
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
    preflight <- tryCatch(
      {
        ccc_preflight_method(
          method = method,
          srt = srt,
          params = method_params[[method]] %||% list()
        )
        NULL
      },
      error = function(e) e
    )
    if (inherits(preflight, "error")) {
      status[[method]] <- list(
        method = method,
        status = "failed",
        elapsed = 0,
        message = conditionMessage(preflight)
      )
      if (!isTRUE(skip_failed)) stop(preflight)
      log_message(
        "CCC method {.val {method}} failed preflight and was skipped: {conditionMessage(preflight)}",
        message_type = "warning",
        verbose = verbose
      )
      next
    }
    args <- ccc_run_method_args(
      method = method,
      srt = srt,
      group.by = group.by,
      backend = backend,
      verbose = verbose,
      thresh = thresh,
      params = method_params[[method]] %||% list()
    )
    fun <- get(registry[[method]]$producer, mode = "function")

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
      backend_scope = "result aggregation and unified-table construction",
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
  entry <- ccc_registry_entry(method)
  protected <- intersect(c("srt", "object", "group.by"), names(params))
  if (length(protected) > 0L) {
    log_message(
      "Ignoring protected {.arg method_params} entries for {.val {method}}: {.val {protected}}",
      message_type = "warning",
      verbose = verbose
    )
    params[protected] <- NULL
  }

  base <- list(group.by = group.by, verbose = verbose)
  base[[entry$object_arg]] <- srt
  if (!identical(method, "MDIC3") && !"backend" %in% names(params)) {
    base$backend <- backend
  }
  if (identical(method, "CellChat") && !"thresh" %in% names(params)) {
    base$thresh <- thresh
  }
  utils::modifyList(base, params, keep.null = TRUE)
}

ccc_preflight_method <- function(method, srt, params = list()) {
  entry <- ccc_registry_entry(method)
  required <- if (identical(method, "MDIC3")) character(0) else entry$required
  missing <- required[
    !required %in% names(params) |
      vapply(required, function(nm) {
        value <- params[[nm]]
        is.null(value) || length(value) == 0L ||
          (is.character(value) && all(is.na(value) | !nzchar(value)))
      }, logical(1))
  ]
  if (length(missing) > 0L) {
    log_message(
      "CCC method {.val {method}} requires method_params fields: {.val {missing}}",
      message_type = "error"
    )
  }
  if (identical(method, "Nichenetr")) {
    mode <- params$mode %||% "aggregate"
    if (!identical(mode, "custom")) {
      design_required <- c(
        "condition.by", "condition_oi", "condition_reference"
      )
      design_missing <- design_required[
        !design_required %in% names(params) |
          vapply(design_required, function(nm) {
            value <- params[[nm]]
            is.null(value) || length(value) == 0L ||
              (is.character(value) && all(is.na(value) | !nzchar(value)))
          }, logical(1))
      ]
      if (length(design_missing) > 0L) {
        log_message(
          "CCC method {.val Nichenetr} in {.val {mode}} mode requires method_params fields: {.val {design_missing}}",
          message_type = "error"
        )
      }
    }
  }
  if (
    identical(method, "MDIC3") &&
      is.null(params$grn) &&
      is.null(params$grn_method)
  ) {
    log_message(
      "{.val MDIC3} requires an explicit {.arg grn} or {.arg grn_method} in {.arg method_params$MDIC3}",
      message_type = "error"
    )
  }
  if (identical(method, "SpatialCellChat")) {
    coord_cols <- params$coord.cols %||% c("col", "row")
    has_coords <- all(coord_cols %in% colnames(srt[[]]))
    has_images <- length(srt@images) > 0L
    if (!has_coords && !has_images) {
      log_message(
        "{.val SpatialCellChat} requires a spatial image or metadata coordinate columns supplied through {.arg method_params$SpatialCellChat$coord.cols}",
        message_type = "error"
      )
    }
  }
  invisible(TRUE)
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
