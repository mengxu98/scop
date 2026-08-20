#' @title Benchmark single-cell integration methods
#'
#' @description
#' Run one or more [RunIntegration()] methods on the same input, score each
#' latent embedding with scIB-style batch-correction and bio-conservation
#' metrics (scaled iLISI/cLISI, ASW, graph connectivity, ARI/NMI), and store
#' the tables in `srt@tools$IntegrationBenchmark`. The object itself is
#' returned so later plotting can reuse the embeddings and per-cell LISI
#' scores.
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
#' @param tool_name Name of the `srt@tools` entry used to store the tables.
#'
#' @return The input `Seurat` object with integration reductions, per-cell
#' LISI columns, and `srt@tools[[tool_name]]` containing `summary`, `metrics`,
#' and `runs`. Print the tables with [thisplot::print_colored_table()].
#'
#' @seealso [IntegrationBenchmarkPlot], [RunIntegration], [RunLISI]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(panc8_sub)
#' panc8_sub <- RunIntegrationBenchmark(
#'   panc8_sub,
#'   batch = "tech",
#'   celltype = "celltype",
#'   methods = c("Uncorrected", "Harmony"),
#'   nHVF = 500,
#'   linear_reduction_dims = 20,
#'   linear_reduction_dims_use = 1:10,
#'   perplexity = 10
#' )
#' thisplot::print_colored_table(
#'   panc8_sub@tools$IntegrationBenchmark$summary,
#'   by = "row",
#'   palette = "Chinese"
#' )
#' thisplot::print_colored_table(
#'   panc8_sub@tools$IntegrationBenchmark$metrics,
#'   by = "col",
#'   palette = "Paired"
#' )
#' IntegrationBenchmarkPlot(panc8_sub)
#' IntegrationBenchmarkPlot(panc8_sub, plot_type = "box")
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
  tool_name = "IntegrationBenchmark",
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
    raw$method <- method
    metrics_list[[method]] <- raw
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
  srt_work@tools[[tool_name]] <- list(
    summary = summary_df,
    metrics = metrics_df,
    runs = runs_df,
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
  )
  srt_work
}
