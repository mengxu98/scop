run_ora_result <- function(
  gene,
  TERM2GENE,
  TERM2NAME,
  IDtype,
  min_gs_size = 10,
  max_gs_size = 500
) {
  max_size <- if (is.finite(max_gs_size)) {
    as.integer(max_gs_size)
  } else {
    .Machine$integer.max
  }
  term_names <- TERM2NAME
  if (!"Name" %in% colnames(term_names)) {
    term_names[["Name"]] <- term_names[["Term"]]
  }

  res <- ora_hypergeom(
    genes = as.character(gene),
    term_ids = as.character(TERM2GENE[["Term"]]),
    term_genes = as.character(TERM2GENE[[IDtype]]),
    term_name_ids = as.character(term_names[["Term"]]),
    term_names = as.character(term_names[["Name"]]),
    min_size = as.integer(min_gs_size),
    max_size = max_size
  )
  if (nrow(res) == 0L) {
    return(NULL)
  }

  res[["GeneRatio"]] <- paste0(res[["Count"]], "/", res[["query_size"]])
  res[["BgRatio"]] <- paste0(res[["term_size"]], "/", res[["universe_size"]])
  res[["RichFactor"]] <- res[["Count"]] / res[["term_size"]]
  res[["FoldEnrichment"]] <- (res[["Count"]] / res[["query_size"]]) /
    (res[["term_size"]] / res[["universe_size"]])

  p <- res[["term_size"]] / res[["universe_size"]]
  expected <- res[["query_size"]] * p
  variance <- res[["query_size"]] * p * (1 - p) *
    (res[["universe_size"]] - res[["query_size"]]) /
    pmax(res[["universe_size"]] - 1, 1)
  res[["zScore"]] <- (res[["Count"]] - expected) / sqrt(variance)
  res[["zScore"]][!is.finite(res[["zScore"]])] <- NA_real_
  res[["p.adjust"]] <- stats::p.adjust(res[["pvalue"]], method = "BH")
  res[["qvalue"]] <- res[["p.adjust"]]

  res <- res[order(res[["pvalue"]], res[["p.adjust"]], res[["ID"]]), , drop = FALSE]
  rownames(res) <- NULL
  res[, c(
    "ID", "Description", "GeneRatio", "BgRatio", "RichFactor",
    "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue",
    "geneID", "Count"
  ), drop = FALSE]
}
