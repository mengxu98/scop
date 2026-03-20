#' @title Plots for metabolism pathway scoring
#'
#' @md
#' @inheritParams GSVAPlot
#' @param srt A Seurat object containing the results of RunMetabolism.
#' @param group.by A character vector specifying the grouping variable used in RunMetabolism.
#' @param assay_name The name of the assay or tools slot containing metabolism results.
#' Default is `"METABOLISM"`.
#' @param ... Additional arguments passed to [GSVAPlot].
#'
#' @seealso
#' [RunMetabolism], [EnrichmentPlot]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunMetabolism(
#'   pancreas_sub,
#'   db = c("KEGG", "REACTOME"),
#'   group.by = "CellType",
#'   species = "Mus_musculus",
#'   method = "AUCell"
#' )
#'
#' ht <- MetabolismPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "heatmap",
#'   topTerm = 10,
#'   width = 1,
#'   height = 2
#' )
#'
#' ht2 <- MetabolismPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "heatmap",
#'   n_split = 3,
#'   topTerm = 100,
#'   use_raster = TRUE,
#'   width = 1,
#'   height = 2
#' )
#'
#' MetabolismPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   db = "GO_BP",
#'   plot_type = "comparison",
#'   topTerm = 5
#' )
#'
#' MetabolismPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   db = "GO_BP",
#'   group_use = "Ductal",
#'   plot_type = "bar",
#'   topTerm = 5
#' )
#'
#' MetabolismPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   group_use = "Ductal",
#'   db = "GO_BP",
#'   plot_type = "network",
#'   topTerm = 3
#' )
#'
#' MetabolismPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   group_use = "Ductal",
#'   db = "GO_BP",
#'   plot_type = "enrichmap"
#' )
#'
#' MetabolismPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   group_use = "Ductal",
#'   plot_type = "wordcloud",
#'   word_type = "feature"
#' )
#'
#' pancreas_sub <- RunMetabolism(
#'   pancreas_sub,
#'   assay_name = "METABOLISM",
#'   db = c("KEGG", "REACTOME"),
#'   species = "Mus_musculus"
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   assay = "METABOLISM",
#'   features = rownames(pancreas_sub[["METABOLISM"]])[1:2],
#'   reduction = "umap"
#' )
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = rownames(pancreas_sub[["METABOLISM"]])[1:2],
#'   group.by = "CellType",
#'   assay = "METABOLISM"
#' )
#'
#' ht <- GroupHeatmap(
#'   pancreas_sub,
#'   exp_legend_title = "Z-score",
#'   features = rownames(pancreas_sub[["METABOLISM"]])[1:10],
#'   group.by = "CellType",
#'   assay = "METABOLISM",
#'   width = 1,
#'   height = 2
#' )
MetabolismPlot <- function(
  srt = NULL,
  res = NULL,
  group.by = NULL,
  assay_name = "METABOLISM",
  ...
) {
  if (!is.null(res)) {
    return(GSVAPlot(res = res, assay_name = assay_name, ...))
  }
  if (is.null(srt)) {
    log_message(
      "Either {.arg srt} or {.arg res} must be provided",
      message_type = "error"
    )
  }
  tool_name <- NULL
  if (!is.null(group.by)) {
    for (m in c("AUCell", "GSVA", "ssGSEA", "VISION")) {
      tn <- paste0("Metabolism_", group.by, "_", m)
      if (tn %in% names(srt@tools)) {
        tool_name <- tn
        break
      }
    }
  }
  if (is.null(tool_name)) {
    for (m in c("AUCell", "GSVA", "ssGSEA", "VISION")) {
      tn <- paste0("Metabolism_", m)
      if (tn %in% names(srt@tools)) {
        tool_name <- tn
        break
      }
    }
  }
  if (!is.null(tool_name) && tool_name %in% names(srt@tools)) {
    obj <- srt@tools[[tool_name]]
    res <- list(
      scores = obj[["scores"]],
      enrichment = obj[["enrichment"]],
      group.by = obj[["group.by"]],
      db = obj[["db"]] %||% "Metabolism"
    )
  } else if (assay_name %in% SeuratObject::Assays(srt)) {
    scores <- GetAssayData5(srt, assay = assay_name, layer = "counts")
    res <- list(
      scores = as.matrix(scores),
      group.by = NULL,
      db = "Metabolism"
    )
  } else {
    log_message(
      "Metabolism results not found. Please run RunMetabolism first",
      message_type = "error"
    )
  }
  if (!is.null(res[["enrichment"]]) && is.data.frame(res[["enrichment"]]) &&
    "Database" %in% colnames(res[["enrichment"]])) {
    db_in_res <- unique(as.character(res[["enrichment"]][["Database"]]))
    db_in_res <- db_in_res[!is.na(db_in_res) & nzchar(trimws(db_in_res))]
    if (length(db_in_res) > 0) {
      res[["db"]] <- db_in_res
    }
  } else if (is.null(res[["db"]])) {
    res[["db"]] <- "Metabolism"
  }
  GSVAPlot(res = res, assay_name = assay_name, ...)
}
