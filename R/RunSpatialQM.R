#' @title Run SpatialQM quality metrics
#'
#' @description
#' Run selected `SpatialQM` quality-control metrics on a spatial `Seurat`
#' object and store a compact, scop-style result bundle in `srt@tools`.
#' `SpatialQM` currently expects a `RNA` assay for object-first metrics, so this
#' wrapper maps the selected assay layer to a temporary `RNA` assay without
#' modifying the returned object. `SpatialQM` is an optional GitHub dependency
#' installable with `remotes::install_github("Center-for-Spatial-OMICs/SpatialQM")`.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt A `Seurat` object.
#' @param assay Assay used for SpatialQM input. If `NULL`, the default assay is
#' used.
#' @param layer Assay layer used as counts for SpatialQM.
#' @param metrics SpatialQM metrics to run. Supported aliases include
#' `"n_cells"`, `"tx_per_cell"`, `"tx_per_area"`, `"tx_per_nuc"`,
#' `"mean_expression"`, `"mean_signal_ratio"`, `"cell_tx_fraction"`,
#' `"max_ratio"`, `"max_detection"`, `"mecr"`, `"morans"`, `"silhouette"`,
#' `"sparsity"`, and `"entropy"`.
#' @param features Optional feature vector passed to metrics that support a
#' `features` argument.
#' @param sample_id,platform Optional values used to populate missing
#' `sample_id` and `platform` metadata columns expected by some SpatialQM
#' metrics.
#' @param tool_name Name used to store results in `srt@tools`.
#' @param store_results Whether to store results in `srt@tools`.
#' @param on_error Whether a failed metric should stop the wrapper (`"error"`)
#' or be recorded in the result bundle with a warning (`"warning"`).
#' @param ... Additional named arguments passed to SpatialQM metric functions
#' when the installed function exposes matching formal arguments.
#'
#' @return A `Seurat` object. When `store_results = TRUE`, results are stored in
#' `srt@tools[[tool_name]]`.
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#'
#' spatial <- RunSpatialQM(
#'   spatial,
#'   assay = "Spatial",
#'   layer = "counts",
#'   metrics = c("n_cells", "tx_per_cell", "sparsity", "entropy"),
#'   platform = "Visium",
#'   verbose = FALSE
#' )
#' spatial@tools$SpatialQM$summary
RunSpatialQM <- function(
  srt,
  assay = NULL,
  layer = "counts",
  metrics = c("n_cells", "tx_per_cell", "sparsity", "entropy"),
  features = NULL,
  sample_id = NULL,
  platform = NULL,
  tool_name = "SpatialQM",
  store_results = TRUE,
  on_error = c("error", "warning"),
  verbose = TRUE,
  ...
) {
  log_message(
    "Running SpatialQM quality metrics",
    message_type = "running",
    verbose = verbose
  )
  spatialqm_validate_srt(srt)
  spatialqm_assert_string(tool_name, "tool_name")
  spatialqm_assert_flag(store_results, "store_results")
  on_error <- match.arg(on_error)

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!assay %in% SeuratObject::Assays(srt)) {
    log_message(
      "{.arg assay} {.val {assay}} is not present in {.cls Seurat}",
      message_type = "error"
    )
  }
  metrics <- spatialqm_resolve_metrics(metrics)
  features <- spatialqm_resolve_features(srt, assay = assay, features = features)
  extra_args <- list(...)
  spatialqm_validate_named_list(extra_args, "...")

  check_r("SpatialQM", verbose = FALSE)
  input <- spatialqm_prepare_object(
    srt = srt,
    assay = assay,
    layer = layer,
    sample_id = sample_id,
    platform = platform
  )

  results <- list()
  errors <- list()
  for (metric in names(metrics)) {
    spec <- metrics[[metric]]
    log_message(
      "Run {.pkg SpatialQM} metric {.val {metric}}",
      verbose = verbose
    )
    metric_result <- tryCatch(
      spatialqm_run_metric(
        seu_obj = input$object,
        metric = metric,
        spec = spec,
        features = features,
        extra_args = extra_args
      ),
      error = function(e) e
    )
    if (inherits(metric_result, "error")) {
      errors[[metric]] <- conditionMessage(metric_result)
      msg <- paste0(
        "SpatialQM metric {.val {metric}} failed: ",
        errors[[metric]]
      )
      if (identical(on_error, "error")) {
        log_message(msg, message_type = "error")
      }
      log_message(msg, message_type = "warning", verbose = verbose)
      next
    }
    results[[metric]] <- metric_result
  }

  if (length(results) == 0L) {
    log_message(
      "No SpatialQM metrics completed successfully",
      message_type = "error"
    )
  }
  summary <- spatialqm_summary(results = results, errors = errors)

  if (isTRUE(store_results)) {
    srt@tools[[tool_name]] <- list(
      results = results,
      summary = summary,
      errors = errors,
      parameters = list(
        assay = assay,
        layer = layer,
        metrics = names(metrics),
        features = features,
        sample_id = input$sample_id,
        platform = input$platform,
        tool_name = tool_name,
        store_results = store_results,
        on_error = on_error,
        backend_args = extra_args
      )
    )
  }

  log_message(
    "{.pkg SpatialQM} completed {.val {length(results)}} metric{?s}",
    message_type = "success",
    verbose = verbose
  )
  srt
}

spatialqm_metric_specs <- function() {
  list(
    n_cells = list(fun = "getNcells", pass_features = FALSE),
    tx_per_cell = list(fun = "getTxPerCell", pass_features = TRUE),
    tx_per_area = list(fun = "getTxPerArea", pass_features = TRUE),
    tx_per_nuc = list(fun = "getTxPerNuc", pass_features = TRUE),
    mean_expression = list(fun = "getMeanExpression", pass_features = TRUE),
    mean_signal_ratio = list(fun = "getMeanSignalRatio", pass_features = TRUE),
    cell_tx_fraction = list(fun = "getCellTxFraction", pass_features = TRUE),
    max_ratio = list(fun = "getMaxRatio", pass_features = TRUE),
    max_detection = list(fun = "getMaxDetection", pass_features = TRUE),
    mecr = list(fun = "getMECR", pass_features = FALSE),
    morans = list(fun = "getMorans", pass_features = TRUE),
    silhouette = list(fun = "getSilhouetteWidth", pass_features = FALSE),
    sparsity = list(fun = "getSparsity", pass_features = TRUE),
    entropy = list(fun = "getEntropy", pass_features = TRUE)
  )
}

spatialqm_resolve_metrics <- function(metrics) {
  specs <- spatialqm_metric_specs()
  if (is.null(metrics) || length(metrics) == 0L) {
    log_message("{.arg metrics} must contain at least one metric", message_type = "error")
  }
  metrics <- unique(as.character(metrics))
  aliases <- setNames(names(specs), names(specs))
  aliases[paste0("get", gsub("(^|_)([a-z])", "\\U\\2", names(specs), perl = TRUE))] <- names(specs)
  function_aliases <- c(
    getNcells = "n_cells",
    getTxPerCell = "tx_per_cell",
    getTxPerArea = "tx_per_area",
    getTxPerNuc = "tx_per_nuc",
    getMeanExpression = "mean_expression",
    getMeanSignalRatio = "mean_signal_ratio",
    getCellTxFraction = "cell_tx_fraction",
    getMaxRatio = "max_ratio",
    getMaxDetection = "max_detection",
    getMECR = "mecr",
    getMorans = "morans",
    getSilhouetteWidth = "silhouette",
    getSparsity = "sparsity",
    getEntropy = "entropy"
  )
  aliases[names(function_aliases)] <- unname(function_aliases)
  resolved <- unname(aliases[metrics])
  missing <- metrics[is.na(resolved)]
  if (length(missing) > 0L) {
    log_message(
      "Unknown SpatialQM metric{?s}: {.val {missing}}",
      message_type = "error"
    )
  }
  resolved <- unique(resolved)
  specs[resolved]
}

spatialqm_prepare_object <- function(
  srt,
  assay,
  layer,
  sample_id = NULL,
  platform = NULL
) {
  counts <- GetAssayData5(srt, assay = assay, layer = layer)
  if (!inherits(counts, "Matrix")) {
    counts <- Matrix::Matrix(
      if (is.data.frame(counts)) as.matrix(counts) else counts,
      sparse = TRUE
    )
  }
  if (!inherits(counts, "dgCMatrix")) {
    counts <- methods::as(counts, "CsparseMatrix")
  }
  if (!inherits(counts, "dgCMatrix")) {
    counts <- methods::as(counts, "dgCMatrix")
  }
  counts@x[!is.finite(counts@x)] <- 0
  counts <- Matrix::drop0(counts)

  qobj <- srt
  suppressWarnings(
    qobj[["RNA"]] <- SeuratObject::CreateAssayObject(counts = counts)
  )
  SeuratObject::DefaultAssay(qobj) <- "RNA"

  resolved_sample_id <- sample_id %||% spatialqm_meta_scalar(qobj, "sample_id", "sample1")
  resolved_platform <- platform %||% spatialqm_meta_scalar(qobj, "platform", "scop")
  qobj$sample_id <- resolved_sample_id
  qobj$platform <- resolved_platform
  list(
    object = qobj,
    sample_id = resolved_sample_id,
    platform = resolved_platform
  )
}

spatialqm_run_metric <- function(seu_obj, metric, spec, features = NULL, extra_args = list()) {
  fun <- get_namespace_fun("SpatialQM", spec$fun)
  args <- list(seu_obj = seu_obj)
  if (isTRUE(spec$pass_features)) {
    args$features <- features
  }
  args <- spatialqm_merge_supported_args(fun, args, extra_args)
  do.call(fun, args)
}

spatialqm_merge_supported_args <- function(fun, args, extra_args) {
  if (length(extra_args) == 0L) {
    return(args)
  }
  formal_names <- names(formals(fun))
  if ("..." %in% formal_names) {
    return(utils::modifyList(args, extra_args))
  }
  supported <- intersect(names(extra_args), formal_names)
  utils::modifyList(args, extra_args[supported])
}

spatialqm_summary <- function(results, errors = list()) {
  rows <- lapply(names(results), function(metric) {
    result <- results[[metric]]
    cls <- class(result)[1L] %||% typeof(result)
    if (is.data.frame(result)) {
      value <- spatialqm_first_numeric_value(result)
      return(data.frame(
        metric = metric,
        status = "ok",
        class = cls,
        rows = nrow(result),
        columns = ncol(result),
        value = value,
        message = NA_character_,
        stringsAsFactors = FALSE
      ))
    }
    if (is.atomic(result)) {
      value <- if (length(result) == 1L && is.numeric(result)) as.numeric(result) else NA_real_
      return(data.frame(
        metric = metric,
        status = "ok",
        class = cls,
        rows = length(result),
        columns = NA_integer_,
        value = value,
        message = NA_character_,
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      metric = metric,
      status = "ok",
      class = cls,
      rows = NA_integer_,
      columns = NA_integer_,
      value = NA_real_,
      message = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  err_rows <- lapply(names(errors), function(metric) {
    data.frame(
      metric = metric,
      status = "error",
      class = NA_character_,
      rows = NA_integer_,
      columns = NA_integer_,
      value = NA_real_,
      message = errors[[metric]],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, c(rows, err_rows))
  rownames(out) <- NULL
  out
}

spatialqm_first_numeric_value <- function(x) {
  numeric_cols <- names(x)[vapply(x, is.numeric, logical(1))]
  if (length(numeric_cols) == 0L || nrow(x) == 0L) {
    return(NA_real_)
  }
  as.numeric(x[[numeric_cols[1L]]][1L])
}

spatialqm_resolve_features <- function(srt, assay, features = NULL) {
  if (is.null(features)) {
    return(NULL)
  }
  features <- unique(as.character(features))
  missing <- setdiff(features, rownames(srt[[assay]]))
  if (length(missing) > 0L) {
    log_message(
      "{.arg features} not found in assay {.val {assay}}: {.val {missing}}",
      message_type = "error"
    )
  }
  features
}

spatialqm_meta_scalar <- function(srt, column, fallback) {
  if (!column %in% colnames(srt@meta.data)) {
    return(fallback)
  }
  values <- unique(as.character(srt[[column, drop = TRUE]]))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0L) {
    fallback
  } else {
    values[1L]
  }
}

spatialqm_validate_srt <- function(srt) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  invisible(TRUE)
}

spatialqm_assert_string <- function(x, arg) {
  if (length(x) != 1L || !is.character(x) || is.na(x) || !nzchar(x)) {
    log_message("{.arg {arg}} must be a single non-empty string", message_type = "error")
  }
  invisible(TRUE)
}

spatialqm_assert_flag <- function(x, arg) {
  if (length(x) != 1L || !is.logical(x) || is.na(x)) {
    log_message("{.arg {arg}} must be TRUE or FALSE", message_type = "error")
  }
  invisible(TRUE)
}

spatialqm_validate_named_list <- function(x, arg) {
  if (length(x) == 0L) {
    return(invisible(TRUE))
  }
  nms <- names(x)
  if (is.null(nms) || any(is.na(nms) | !nzchar(nms))) {
    log_message("{.arg {arg}} must contain named arguments only", message_type = "error")
  }
  invisible(TRUE)
}
