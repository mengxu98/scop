#' @title Enrich CellRank trend modules
#'
#' @description
#' Run over-representation analysis for each stored CellRank trend module and
#' keep complete tables and an execution manifest in the Seurat object.
#'
#' @param srt A Seurat object returned by [RunCellRankTrends].
#' @param lineage CellRank lineage whose trend modules should be enriched.
#' @param db Annotation databases accepted by [PrepareDB].
#' @param species Species passed to [PrepareDB].
#' @param universe Background genes. `NULL` uses genes tested for CellRank
#' lineage drivers.
#' @param minGSSize Minimum gene-set size.
#' @param maxGSSize Maximum gene-set size.
#' @param pvalue_cutoff Nominal enrichment cutoff.
#' @param qvalue_cutoff Adjusted enrichment cutoff.
#' @param p_adjust_method Multiple-testing method.
#' @param show_category Number of terms shown in each dot plot.
#' @param output_dir Optional directory for CSV/PDF exports.
#' @param continue_on_error Whether to keep other modules when one enrichment
#' call fails. Failures are recorded in the manifest.
#' @param verbose Whether to print progress messages.
#'
#' @return The Seurat object with results under
#' `srt@tools$CellRank$enrichment[[lineage]]`.
#' @export
RunCellRankEnrichment <- function(
  srt,
  lineage,
  db = c("MSigDB_C2", "GO_BP", "GO_CC", "GO_MF", "MSigDB_H"),
  species = "Mus_musculus",
  universe = NULL,
  minGSSize = 10L,
  maxGSSize = 500L,
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.2,
  p_adjust_method = "BH",
  show_category = 8L,
  output_dir = NULL,
  continue_on_error = FALSE,
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a Seurat object", message_type = "error")
  }
  trends <- srt@tools$CellRank$trends[[lineage]]
  if (is.null(trends) || is.null(trends$cluster_table)) {
    log_message("Run {.fn RunCellRankTrends} before enrichment", message_type = "error")
  }
  check_r("clusterProfiler", verbose = FALSE)
  enricher_fun <- get_namespace_fun("clusterProfiler", "enricher")
  db <- unique(as.character(db))
  if (is.null(universe)) {
    drivers <- srt@tools$CellRank$lineage_drivers
    universe <- if (!is.null(drivers)) rownames(drivers) else unique(trends$cluster_table$gene)
  }
  universe <- unique(as.character(universe))
  clusters <- split(trends$cluster_table$gene, trends$cluster_table$cluster)
  results <- list()
  manifest <- list()
  if (!is.null(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  for (db_i in db) {
    db_list <- PrepareDB(
      species = species,
      db = db_i,
      db_IDtypes = "symbol",
      convert_species = TRUE,
      verbose = verbose
    )
    db_entry <- db_list[[species]][[db_i]]
    if (is.null(db_entry)) {
      manifest[[db_i]] <- list(status = "unavailable", n_clusters = length(clusters))
      next
    }
    term2gene <- db_entry$TERM2GENE
    term2name <- db_entry$TERM2NAME
    db_results <- list()
    db_errors <- list()
    for (cluster_i in names(clusters)) {
      genes_i <- intersect(unique(as.character(clusters[[cluster_i]])), universe)
      if (!length(genes_i)) {
        db_results[[cluster_i]] <- data.frame()
        next
      }
      enrich <- tryCatch(
        enricher_fun(
          gene = genes_i,
          universe = universe,
          TERM2GENE = term2gene,
          TERM2NAME = term2name,
          minGSSize = as.integer(minGSSize),
          maxGSSize = as.integer(maxGSSize),
          pvalueCutoff = pvalue_cutoff,
          qvalueCutoff = qvalue_cutoff,
          pAdjustMethod = p_adjust_method
        ),
        error = function(e) e
      )
      tab <- if (inherits(enrich, "error") || is.null(enrich)) {
        if (inherits(enrich, "error")) {
          db_errors[[cluster_i]] <- conditionMessage(enrich)
          if (!isTRUE(continue_on_error)) {
            log_message(
              "Enrichment failed for {.val {db_i}} / cluster {.val {cluster_i}}: {.val {conditionMessage(enrich)}}",
              message_type = "error"
            )
          }
        }
        data.frame()
      } else {
        as.data.frame(enrich)
      }
      if (nrow(tab)) {
        tab$cluster <- cluster_i
        tab$database <- db_i
      }
      db_results[[cluster_i]] <- tab
    }
    combined <- do.call(rbind, db_results)
    if (is.null(combined)) combined <- data.frame()
    results[[db_i]] <- combined
    manifest[[db_i]] <- list(
      status = if (length(db_errors)) "failed" else if (nrow(combined)) "ok" else "no_significant_terms",
      n_clusters = length(clusters),
      n_terms = nrow(combined),
      errors = db_errors
    )
    if (!is.null(output_dir)) {
      utils::write.csv(
        combined,
        file.path(output_dir, paste0(make.names(lineage), "_", make.names(db_i), "_enrichment.csv")),
        row.names = FALSE
      )
      if (nrow(combined)) {
        p <- ggplot2::ggplot(
          head(combined[order(combined$p.adjust), , drop = FALSE], as.integer(show_category)),
          ggplot2::aes(x = cluster, y = reorder(Description, -log10(p.adjust)), size = Count, color = p.adjust)
        ) + ggplot2::geom_point() + ggplot2::scale_color_viridis_c(trans = "log10", direction = -1) +
          ggplot2::labs(x = "Trend module", y = NULL, title = paste(lineage, db_i)) +
          ggplot2::theme_bw()
        ggplot2::ggsave(
          file.path(output_dir, paste0(make.names(lineage), "_", make.names(db_i), "_enrichment.pdf")),
          p, width = 8, height = 5
        )
      }
    }
  }
  bundle <- list(
    lineage = lineage,
    databases = results,
    manifest = manifest,
    universe = universe,
    parameters = list(
      species = species,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      pvalue_cutoff = pvalue_cutoff,
      qvalue_cutoff = qvalue_cutoff,
      p_adjust_method = p_adjust_method
    )
  )
  srt@tools$CellRank$enrichment[[lineage]] <- bundle
  srt
}
