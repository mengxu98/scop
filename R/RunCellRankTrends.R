#' @title Fit and cluster CellRank lineage trends
#'
#' @description
#' Refit CellRank GAM trends for a selected lineage and store the normalized
#' trend matrix and gene modules in `srt@tools$CellRank$trends`.
#'
#' @param srt A Seurat object returned by [RunCellRank].
#' @param lineage A lineage name in the stored fate-probability matrix.
#' @param features Optional genes. If `NULL`, positive lineage drivers are used.
#' @param top_n Maximum number of driver genes used for clustering.
#' @param heatmap_n Number of genes retained for the compact heatmap view.
#' @param time_key Metadata column containing pseudotime.
#' @param assay Assay containing the expression layer.
#' @param layer Expression layer used for GAM fitting.
#' @param n_points Number of points used to predict each trend.
#' @param norm Whether to z-normalize each gene trend before clustering.
#' @param resolution Leiden resolution for trend modules.
#' @param max_iter Maximum GAM iterations.
#' @param spline_order Spline order passed to CellRank GAM.
#' @param n_knots Number of GAM knots.
#' @param cores Number of Python workers.
#' @param random_state Reproducibility seed.
#' @param output_dir Optional directory for CSV exports.
#' @param envname Optional Python environment name.
#' @param conda Conda-compatible executable used by [PrepareEnv].
#' @param min_expressed_cells Minimum cells with positive finite expression
#' required for a GAM feature.
#' @param distribution Primary CellRank GAM distribution.
#' @param link Link function for the primary GAM distribution.
#' @param fallback_distribution Optional distribution used when the primary
#' model has a fatal fit failure; the selected distribution is recorded.
#' @param verbose Whether to print progress messages.
#'
#' @return The Seurat object with a `CellRank$trends[[lineage]]` result bundle.
#' @export
RunCellRankTrends <- function(
  srt,
  lineage,
  features = NULL,
  top_n = 500L,
  heatmap_n = 50L,
  time_key = "palantir_pseudotime",
  assay = "RNA",
  layer = "data",
  n_points = 200L,
  norm = TRUE,
  resolution = 0.7,
  max_iter = 1000L,
  spline_order = 3L,
  n_knots = 8L,
  cores = 1L,
  random_state = 0L,
  output_dir = NULL,
  envname = NULL,
  conda = "auto",
  min_expressed_cells = 20L,
  distribution = "gamma",
  link = "log",
  fallback_distribution = "normal",
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a Seurat object", message_type = "error")
  }
  if (is.null(srt@tools$CellRank)) {
    log_message("Run {.fn RunCellRank} before {.fn RunCellRankTrends}", message_type = "error")
  }
  fate <- cellrank_fate_matrix(srt)
  if (!lineage %in% colnames(fate)) {
    log_message(
      "{.arg lineage} must match a stored CellRank fate-probability column",
      message_type = "error"
    )
  }
  if (!time_key %in% colnames(srt@meta.data)) {
    log_message("Pseudotime column {.val {time_key}} is missing", message_type = "error")
  }
  if (!assay %in% names(srt)) {
    log_message("Assay {.val {assay}} is missing", message_type = "error")
  }
  if (!layer %in% SeuratObject::Layers(srt[[assay]])) {
    log_message("Layer {.val {layer}} is missing from assay {.val {assay}}", message_type = "error")
  }

  drivers <- srt@tools$CellRank$lineage_drivers
  if (is.null(features)) {
    corr_col <- paste0(lineage, "_corr")
    if (is.null(drivers) || !corr_col %in% colnames(drivers)) {
      log_message("No stored driver correlation column {.val {corr_col}}", message_type = "error")
    }
    corr_values <- drivers[, corr_col]
    features <- rownames(drivers)[is.finite(corr_values) & corr_values > 0]
    features <- features[order(drivers[features, corr_col], decreasing = TRUE)]
  }
  features <- unique(intersect(as.character(features), rownames(srt[[assay]])))
  if (!length(features)) {
    log_message("No trend genes overlap the selected assay", message_type = "error")
  }
  features <- head(features, as.integer(top_n))
  expression <- SeuratObject::LayerData(srt[[assay]], layer = layer)[features, , drop = FALSE]
  positive_cells <- Matrix::rowSums(expression > 0)
  row_mean <- Matrix::rowMeans(expression)
  row_var <- Matrix::rowSums(expression^2) / ncol(expression) - row_mean^2
  valid_features <- features[
    positive_cells >= as.integer(min_expressed_cells) &
      is.finite(row_var) & row_var > 0
  ]
  if (!length(valid_features)) {
    log_message(
      "No trend genes pass the finite, variable, and expressed-cell filters",
      message_type = "error"
    )
  }
  features <- valid_features

  PrepareEnv(
    envname = envname,
    conda = conda,
    modules = c("scanpy", "cellrank"),
    verbose = verbose
  )
  check_python("cellrank", envname = envname, conda = conda, verbose = verbose)
  adata <- srt_to_adata(
    srt = srt,
    assay_x = assay,
    layer_x = layer,
    prepare_env = FALSE
  )
  functions <- scop_python_import("functions", convert = TRUE)
  result <- do.call(functions$CellRankTrends, list(
    adata = adata,
    fate_probabilities = reticulate::r_to_py(fate),
    fate_columns = colnames(fate),
    genes = features,
    lineage = lineage,
    time_key = time_key,
    n_points = as.integer(n_points),
    norm = isTRUE(norm),
    resolution = as.numeric(resolution),
    max_iter = as.integer(max_iter),
    spline_order = as.integer(spline_order),
    n_knots = as.integer(n_knots),
    n_jobs = as.integer(cores),
    random_state = as.integer(random_state),
    distribution = distribution,
    link = link,
    fallback_distribution = fallback_distribution,
    show_plot = FALSE,
    save = NULL
  ))

  trend_matrix <- as.matrix(py_to_r2(result$trend_matrix))
  genes_out <- as.character(py_to_r2(result$genes))
  time_out <- as.numeric(py_to_r2(result$time))
  clusters <- as.character(py_to_r2(result$clusters))
  rownames(trend_matrix) <- genes_out
  colnames(trend_matrix) <- as.character(time_out)
  cluster_table <- data.frame(
    gene = genes_out,
    cluster = clusters,
    stringsAsFactors = FALSE
  )
  heatmap_genes <- head(genes_out, as.integer(heatmap_n))
  bundle <- list(
    lineage = lineage,
    genes = genes_out,
    heatmap_genes = heatmap_genes,
    time = time_out,
    trend_matrix = trend_matrix,
    cluster_table = cluster_table,
    parameters = py_to_r2(result$parameters),
    source = list(assay = assay, layer = layer, time_key = time_key)
  )
  srt@tools$CellRank$trends[[lineage]] <- bundle

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(
      cluster_table,
      file.path(output_dir, paste0(make.names(lineage), "_trend_clusters.csv")),
      row.names = FALSE
    )
    utils::write.csv(
      data.frame(gene = genes_out, trend_matrix, check.names = FALSE),
      file.path(output_dir, paste0(make.names(lineage), "_trend_matrix.csv")),
      row.names = FALSE
    )
  }
  srt
}

cellrank_fate_matrix <- function(srt) {
  fate <- srt@tools$CellRank$fate_probabilities %||%
    srt@tools$CellRank$absorption_probabilities
  if (is.null(fate)) {
    log_message("CellRank fate probabilities are missing", message_type = "error")
  }
  fate <- as.matrix(fate)
  cells <- colnames(srt)
  if (is.null(rownames(fate)) || anyDuplicated(rownames(fate))) {
    log_message(
      "CellRank fate probabilities must have unique cell row names",
      message_type = "error"
    )
  }
  if (is.null(colnames(fate)) || anyDuplicated(colnames(fate)) ||
    any(!nzchar(colnames(fate)))) {
    log_message(
      "CellRank fate probabilities must have unique lineage names",
      message_type = "error"
    )
  }
  missing_cells <- setdiff(cells, rownames(fate))
  if (length(missing_cells)) {
    log_message(
      "CellRank fate probabilities are missing cells: {.val {paste(head(missing_cells, 5), collapse = ', ')}}",
      message_type = "error"
    )
  }
  fate <- fate[cells, , drop = FALSE]
  if (any(!is.finite(fate)) || any(rowSums(fate) <= 0)) {
    log_message("CellRank fate probabilities contain invalid rows", message_type = "error")
  }
  fate
}
