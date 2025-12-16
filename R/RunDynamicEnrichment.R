#' @title RunDynamicEnrichment
#'
#' @description
#' This function calculates gene-set scores from the specified database (`db`) for each lineage using the specified scoring method (`score_method`).
#' It then treats these scores as expression values and uses them as input to the RunDynamicFeatures function to identify dynamically enriched terms along the lineage.
#'
#' @md
#' @inheritParams RunEnrichment
#' @inheritParams DynamicHeatmap
#' @inheritParams CellScoring
#' @param score_method The method to use for scoring.
#' Can be `"Seurat"`, `"AUCell"`, or `"UCell"`.
#' Default is `"Seurat"`.
#'
#' @seealso
#' [RunDynamicFeatures], [DynamicHeatmap]
#'
#' @export
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunSlingshot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP"
#' )
#' pancreas_sub <- RunDynamicFeatures(
#'   pancreas_sub,
#'   lineages = "Lineage1",
#'   n_candidates = 200
#' )
#' ht1 <- DynamicHeatmap(
#'   pancreas_sub,
#'   lineages = "Lineage1",
#'   cell_annotation = "SubCellType",
#'   n_split = 4
#' )
#' ht1$plot
#'
#' pancreas_sub <- RunDynamicEnrichment(
#'   pancreas_sub,
#'   lineages = "Lineage1",
#'   score_method = "UCell",
#'   db = "GO_BP",
#'   species = "Mus_musculus"
#' )
#' ht2 <- DynamicHeatmap(
#'   pancreas_sub,
#'   assay = "GO_BP",
#'   lineages = "Lineage1_GO_BP",
#'   cell_annotation = "SubCellType",
#'   n_split = 4,
#'   split_method = "kmeans-peaktime"
#' )
#' ht2$plot
RunDynamicEnrichment <- function(
    srt,
    lineages,
    score_method = "AUCell",
    layer = "data",
    assay = NULL,
    min_expcells = 20,
    r.sq = 0.2,
    dev.expl = 0.2,
    padjust = 0.05,
    IDtype = "symbol",
    species = "Homo_sapiens",
    db = "GO_BP",
    db_update = FALSE,
    db_version = "latest",
    convert_species = TRUE,
    Ensembl_version = NULL,
    mirror = NULL,
    TERM2GENE = NULL,
    TERM2NAME = NULL,
    minGSSize = 10,
    maxGSSize = 500,
    cores = 1,
    verbose = TRUE,
    seed = 11) {
  set.seed(seed)
  assay <- assay %||% DefaultAssay(srt)

  feature_union <- c()
  cell_union <- c()
  dynamic <- list()
  for (l in lineages) {
    if (!paste0("DynamicFeatures_", l) %in% names(srt@tools)) {
      log_message(
        "{.val {l}} info not found in the srt object. Should perform {.fn RunDynamicFeatures} first",
        message_type = "error",
        verbose = verbose
      )
    }
    DynamicFeatures <- srt@tools[[paste0("DynamicFeatures_", l)]][[
      "DynamicFeatures"
    ]]
    DynamicFeatures <- DynamicFeatures[
      DynamicFeatures$exp_ncells > min_expcells &
        DynamicFeatures$r.sq > r.sq &
        DynamicFeatures$dev.expl > dev.expl &
        DynamicFeatures$padjust < padjust, ,
      drop = FALSE
    ]
    dynamic[[l]] <- DynamicFeatures
    feature_union <- c(feature_union, DynamicFeatures[, "features"])
    cell_union <- c(
      cell_union,
      rownames(srt@tools[[paste0("DynamicFeatures_", l)]][["raw_matrix"]])
    )
  }
  feature_union <- unique(feature_union)

  if (is.null(TERM2GENE)) {
    db_list <- PrepareDB(
      species = species,
      db = db,
      db_update = db_update,
      db_version = db_version,
      db_IDtypes = IDtype,
      convert_species = convert_species,
      Ensembl_version = Ensembl_version,
      mirror = mirror
    )
  } else {
    colnames(TERM2GENE) <- c("Term", IDtype)
    db <- "custom"
    db_list <- list()
    db_list[[species]][[db]][["TERM2GENE"]] <- unique(TERM2GENE)
    if (is.null(TERM2NAME)) {
      TERM2NAME <- unique(TERM2GENE)[, c(1, 1)]
      colnames(TERM2NAME) <- c("Term", "Name")
    }
    db_list[[species]][[db]][["TERM2NAME"]] <- unique(TERM2NAME)
  }

  for (i in seq_along(db)) {
    term <- db[i]
    TERM2GENE_tmp <- db_list[[species]][[term]][["TERM2GENE"]][, c(
      "Term",
      IDtype
    )]
    TERM2NAME_tmp <- db_list[[species]][[term]][["TERM2NAME"]]
    dup <- duplicated(TERM2GENE_tmp)
    na <- Matrix::rowSums(is.na(TERM2GENE_tmp)) > 0
    TERM2GENE_tmp <- TERM2GENE_tmp[!(dup | na), , drop = FALSE]
    TERM2NAME_tmp <- TERM2NAME_tmp[
      TERM2NAME_tmp[, "Term"] %in% TERM2GENE_tmp[, "Term"], ,
      drop = FALSE
    ]

    term_use <- unique(TERM2GENE_tmp[
      TERM2GENE_tmp[, IDtype] %in% feature_union,
      "Term"
    ])
    TERM2GENE_tmp <- TERM2GENE_tmp[
      TERM2GENE_tmp[, "Term"] %in% term_use, ,
      drop = FALSE
    ]
    TERM2NAME_tmp <- TERM2NAME_tmp[
      TERM2NAME_tmp[, "Term"] %in% term_use, ,
      drop = FALSE
    ]
    rownames(TERM2NAME_tmp) <- TERM2NAME_tmp[, "Term"]
    TERM2GENE_tmp <- TERM2GENE_tmp[
      TERM2GENE_tmp[, IDtype] %in% rownames(srt[[assay]]), ,
      drop = FALSE
    ]
    feature_list <- split(
      TERM2GENE_tmp[, IDtype],
      TERM2NAME_tmp[TERM2GENE_tmp[, "Term"], "Name"]
    )
    GSSize <- sapply(feature_list, length)
    feature_list <- feature_list[GSSize >= minGSSize & GSSize <= maxGSSize]

    srt <- CellScoring(
      srt = srt,
      features = feature_list,
      method = score_method,
      classification = FALSE,
      layer = layer,
      assay = assay,
      name = term,
      new_assay = TRUE,
      cores = cores,
      verbose = verbose
    )
    srt <- RunDynamicFeatures(
      srt = srt,
      lineages = lineages,
      features = rownames(
        GetAssayData5(
          srt,
          assay = term,
          layer = "counts"
        )
      ),
      suffix = paste(lineages, term, sep = "_"),
      assay = term
    )
  }

  log_message(
    "Dynamic enrichment analysis completed",
    message_type = "success",
    verbose = verbose
  )

  return(srt)
}
