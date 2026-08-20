#' @title Benchmark single-cell integration methods
#'
#' @description
#' Run one or more [RunIntegration()] methods on the same input, score each
#' embedding with scIB-style batch-correction and bio-conservation metrics
#' (scaled iLISI/cLISI, ASW, graph connectivity, ARI/NMI), and return a
#' compact result for [IntegrationBenchmarkPlot()]. Spatial-domain clustering
#' benchmarks stay in [RunBenchmark()].
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams RunIntegration
#' @param srt A `Seurat` object containing a batch column.
#' @param batch Metadata column used as the technical batch.
#' @param celltype Optional metadata column used as biological labels.
#' @param methods Integration methods passed to [RunIntegration()].
#' @param perplexity Neighborhood size for [RunLISI()].
#' @param k_graph Neighbor size for graph connectivity.
#' @param bio_weight,batch_weight Weights used for the overall score.
#'   Defaults follow scIB (`0.6` bio conservation, `0.4` batch correction).
#' @param skip_failed Whether a failed method is recorded and skipped.
#' @param keep_objects Whether to keep a per-method `Seurat` object in the
#'   result. The combined object used for plotting is always stored in `$srt`.
#' @param print_table Whether to print the metric tables with
#'   [thisplot::print_colored_table()] when the result is printed.
#'
#' @return An `integration_benchmark_result` list with `summary`, `metrics`,
#' `runs`, and `srt`. Printing uses [thisplot::print_colored_table()].
#'
#' @seealso [IntegrationBenchmarkPlot], [RunIntegration], [RunLISI], [LISIPlot]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(panc8_sub)
#' bench <- RunIntegrationBenchmark(
#'   panc8_sub,
#'   batch = "tech",
#'   celltype = "celltype",
#'   methods = c("Uncorrected", "Harmony"),
#'   nHVF = 500,
#'   linear_reduction_dims = 20,
#'   linear_reduction_dims_use = 1:10,
#'   perplexity = 10
#' )
#' bench
#' IntegrationBenchmarkPlot(bench)
#' }
RunIntegrationBenchmark <- function(
  srt,
  batch,
  celltype = NULL,
  methods = c("Uncorrected", "Harmony"),
  linear_reduction = "pca",
  nHVF = 2000,
  linear_reduction_dims = 50,
  linear_reduction_dims_use = NULL,
  nonlinear_reduction = "umap",
  perplexity = 30,
  k_graph = 15,
  bio_weight = 0.6,
  batch_weight = 0.4,
  skip_failed = TRUE,
  keep_objects = FALSE,
  print_table = TRUE,
  verbose = TRUE,
  ...
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (length(batch) != 1L || !batch %in% colnames(srt@meta.data)) {
    log_message(
      "{.arg batch} must be one metadata column in {.arg srt}",
      message_type = "error"
    )
  }
  if (!is.null(celltype) && !celltype %in% colnames(srt@meta.data)) {
    log_message(
      "{.arg celltype} must be a metadata column in {.arg srt}",
      message_type = "error"
    )
  }
  methods <- unique(as.character(methods))
  methods <- methods[nzchar(methods)]
  if (length(methods) == 0L) {
    log_message(
      "{.arg methods} must contain at least one integration method",
      message_type = "error"
    )
  }
  extra <- list(...)
  srt_work <- srt
  objects <- list()
  metrics_list <- list()
  runs <- list()
  lisi_labels <- unique(c(batch, celltype))
  lisi_labels <- lisi_labels[!is.na(lisi_labels) & nzchar(lisi_labels)]

  for (method in methods) {
    log_message(
      "Benchmark integration method {.val {method}}",
      verbose = verbose
    )
    t0 <- proc.time()[["elapsed"]]
    args <- c(
      list(
        srt_merge = srt_work,
        batch = batch,
        append = TRUE,
        integration_methods = method,
        linear_reduction = linear_reduction,
        nHVF = nHVF,
        linear_reduction_dims = linear_reduction_dims,
        linear_reduction_dims_use = linear_reduction_dims_use,
        nonlinear_reduction = nonlinear_reduction,
        compute_lisi = FALSE,
        compute_metrics = FALSE,
        verbose = verbose
      ),
      extra
    )
    args <- args[!duplicated(names(args))]
    integrated <- tryCatch(
      do.call(RunIntegration, args),
      error = function(e) e
    )
    elapsed <- proc.time()[["elapsed"]] - t0
    if (inherits(integrated, "error")) {
      msg <- conditionMessage(integrated)
      log_message(
        "Integration method {.val {method}} failed: {msg}",
        message_type = "warning",
        verbose = verbose
      )
      runs[[method]] <- data.frame(
        method = method,
        status = "failed",
        error = msg,
        runtime_s = elapsed,
        latent = NA_character_,
        umap = NA_character_,
        stringsAsFactors = FALSE
      )
      if (!isTRUE(skip_failed)) {
        log_message(msg, message_type = "error")
      }
      next
    }
    srt_work <- integrated
    latent <- resolve_integration_latent_reduction(
      srt_work,
      method,
      linear_reduction = linear_reduction
    )
    umap <- resolve_integration_umap_reduction(srt_work, method)
    cluster_col <- resolve_integration_cluster_col(
      srt_work,
      method,
      linear_reduction = linear_reduction
    )
    if (is.na(latent)) {
      msg <- paste0("No latent reduction found for ", method)
      runs[[method]] <- data.frame(
        method = method,
        status = "failed",
        error = msg,
        runtime_s = elapsed,
        latent = NA_character_,
        umap = umap,
        stringsAsFactors = FALSE
      )
      if (!isTRUE(skip_failed)) {
        log_message(msg, message_type = "error")
      }
      next
    }
    srt_work <- RunLISI(
      srt_work,
      reductions = latent,
      label_colnames = lisi_labels,
      prefix = method,
      tool_name = paste0(method, "_LISI"),
      perplexity = perplexity,
      verbose = verbose
    )
    raw <- collect_integration_metrics(
      srt = srt_work,
      reduction = latent,
      batch_col = batch,
      celltype_col = celltype,
      cluster_col = if (is.na(cluster_col)) NULL else cluster_col,
      lisi_tool_name = paste0(method, "_LISI"),
      lisi_prefix = method,
      k_graph = k_graph
    )
    if (!"method" %in% colnames(raw)) {
      raw$method <- method
    } else {
      raw$method <- method
    }
    metrics_list[[method]] <- raw
    if (isTRUE(keep_objects)) {
      objects[[method]] <- srt_work
    }
    runs[[method]] <- data.frame(
      method = method,
      status = "success",
      error = "",
      runtime_s = proc.time()[["elapsed"]] - t0,
      latent = latent,
      umap = umap,
      stringsAsFactors = FALSE
    )
  }

  if (length(metrics_list) == 0L) {
    log_message(
      "No integration method completed successfully",
      message_type = "error"
    )
  }
  metrics_df <- do.call(rbind, metrics_list)
  rownames(metrics_df) <- NULL
  overall_df <- integration_overall_scores(
    metrics_df,
    bio_weight = bio_weight,
    batch_weight = batch_weight
  )
  runs_df <- do.call(rbind, runs)
  rownames(runs_df) <- NULL
  summary_df <- integration_summary_wide(metrics_df, overall_df, runs_df)
  result <- structure(
    list(
      summary = summary_df,
      metrics = metrics_df,
      runs = runs_df,
      srt = srt_work,
      objects = objects,
      batch = batch,
      celltype = celltype,
      parameters = list(
        methods = methods,
        linear_reduction = linear_reduction,
        nHVF = nHVF,
        perplexity = perplexity,
        k_graph = k_graph,
        bio_weight = bio_weight,
        batch_weight = batch_weight
      )
    ),
    class = c("integration_benchmark_result", "list")
  )
  attr(result, "print_table") <- isTRUE(print_table)
  result
}

#' @export
#' @method print integration_benchmark_result
print.integration_benchmark_result <- function(x, ...) {
  n_ok <- sum(x$runs$status == "success")
  n_all <- nrow(x$runs)
  cat(
    "Integration benchmark: ",
    n_ok,
    "/",
    n_all,
    " methods succeeded\n",
    sep = ""
  )
  if (!isFALSE(attr(x, "print_table"))) {
    cat("\nSummary (scaled 0-1; higher is better)\n")
    print_df <- x$summary
    keep <- intersect(
      c(
        "method",
        "status",
        "overall",
        "bio",
        "batch",
        "iLISI",
        "cLISI",
        "celltype_ASW",
        "batch_ASW_mixing",
        "celltype_graph_connectivity",
        "celltype_ARI",
        "celltype_NMI",
        "runtime_s"
      ),
      colnames(print_df)
    )
    print_df <- print_df[, keep, drop = FALSE]
    num_cols <- vapply(print_df, is.numeric, logical(1))
    print_df[num_cols] <- lapply(print_df[num_cols], function(v) {
      round(v, 3)
    })
    thisplot::print_colored_table(print_df, by = "row", palette = "Chinese")
    cat("\nAll metrics\n")
    metrics_print <- x$metrics
    metrics_print <- metrics_print[
      !grepl("_LISI_mean$", metrics_print$metric),
      ,
      drop = FALSE
    ]
    if ("value" %in% colnames(metrics_print)) {
      metrics_print$value <- round(metrics_print$value, 4)
    }
    if ("scaled" %in% colnames(metrics_print)) {
      metrics_print$scaled <- round(metrics_print$scaled, 3)
    }
    thisplot::print_colored_table(metrics_print, by = "col", palette = "Paired")
  }
  invisible(x)
}
