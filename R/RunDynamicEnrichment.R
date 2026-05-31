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
#' Can be `"Seurat"`, `"AUCell"`, `"UCell"`, `"GSVA"`, `"ssGSEA"`,
#' `"zscore"`, `"PLAGE"`, or `"VISION"`.
#' Multiple methods can be supplied at once; each method will be written to a
#' method-suffixed assay before dynamic-feature fitting. Default is `"AUCell"`.
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
#'   group.by = "CellType",
#'   reduction = "UMAP"
#' )
#' pancreas_sub <- RunDynamicFeatures(
#'   pancreas_sub,
#'   lineages = "Lineage1",
#'   fit_method = "pretsa",
#'   n_candidates = 200
#' )
#' ht1 <- DynamicHeatmap(
#'   pancreas_sub,
#'   lineages = "Lineage1",
#'   cell_annotation = "CellType",
#'   n_split = 3
#' )
#'
#' pancreas_sub <- RunDynamicEnrichment(
#'   pancreas_sub,
#'   lineages = "Lineage1",
#'   score_method = "AUCell",
#'   db = "GO_BP",
#'   species = "Mus_musculus"
#' )
#' ht2 <- DynamicHeatmap(
#'   pancreas_sub,
#'   assay = "GO_BP",
#'   lineages = "Lineage1_GO_BP",
#'   cell_annotation = "CellType",
#'   n_split = 3,
#'   split_method = "kmeans-peaktime"
#' )
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
  features = NULL,
  TERM2GENE = NULL,
  TERM2NAME = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  backend = c("cpp", "r"),
  cpp_strategy = c("sparse", "topk", "full"),
  cores = 1,
  verbose = TRUE,
  seed = 11,
  ...
) {
  set.seed(seed)
  backend_missing <- missing(backend)
  score_method <- unique(vapply(
    as.character(score_method),
    normalize_gene_set_scoring_method,
    FUN.VALUE = character(1),
    arg_name = "score_method"
  ))
  backend <- match.arg(backend)
  cpp_strategy <- match.arg(cpp_strategy)
  cpp_supported_methods <- c("AUCell", "Seurat", "GSVA", "ssGSEA", "zscore", "PLAGE")
  score_method_cpp_unsupported <- setdiff(score_method, cpp_supported_methods)
  if (!identical(backend, "r") && length(score_method_cpp_unsupported) > 0L && !isTRUE(backend_missing)) {
    log_message(
      "{.arg backend = 'cpp'} currently supports {.arg score_method} values {.val {cpp_supported_methods}} only",
      message_type = "error"
    )
  }
  assay <- assay %||% DefaultAssay(srt)

  feature_union <- c()
  cell_union <- c()
  dynamic <- list()
  for (l in lineages) {
    if (!paste0("DynamicFeatures_", l) %in% names(srt@tools)) {
      log_message(
        "{.val {l}} info not found in the srt object. Should perform {.fn RunDynamicFeatures} first",
        message_type = "error"
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

  if (!is.null(features)) {
    feature_list <- normalize_feature_list_input(features)
    feature_list <- lapply(feature_list, function(gs) {
      intersect(gs, rownames(srt[[assay]]))
    })
    feature_overlap <- vapply(feature_list, function(gs) {
      length(intersect(gs, feature_union)) > 0L
    }, logical(1))
    feature_list <- feature_list[feature_overlap]
    gs_size <- lengths(feature_list)
    feature_list <- feature_list[gs_size >= minGSSize & gs_size <= maxGSSize]
  } else if (is.null(TERM2GENE)) {
    db_list <- PrepareDB(
      species = species,
      db = db,
      db_update = db_update,
      db_version = db_version,
      db_IDtypes = IDtype,
      convert_species = convert_species,
      Ensembl_version = Ensembl_version,
      mirror = mirror,
      ...
    )
  } else {
    db <- "custom"
    custom_db <- create_custom_db_list(
      species = species,
      db = db,
      TERM2GENE = TERM2GENE,
      TERM2NAME = TERM2NAME,
      IDtype = IDtype
    )
    db_list <- custom_db[["db_list"]]
    TERM2GENE <- custom_db[["TERM2GENE"]]
    TERM2NAME <- custom_db[["TERM2NAME"]]
  }

  db_to_score <- if (!is.null(features)) {
    "custom"
  } else {
    db
  }

  for (term in db_to_score) {
    if (is.null(features)) {
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
      gs_size <- sapply(feature_list, length)
      feature_list <- feature_list[gs_size >= minGSSize & gs_size <= maxGSSize]
    }

    for (method_i in score_method) {
      backend_i <- backend
      if (!identical(backend, "r") && !method_i %in% cpp_supported_methods) {
        log_message(
          "{.arg score_method = {.val {method_i}}} does not have a C++ backend yet; using {.arg backend = 'r'} for this run.",
          message_type = "warning",
          verbose = verbose
        )
        backend_i <- "r"
      }
      assay_name_i <- if (length(score_method) > 1L) {
        paste(term, method_i, sep = "_")
      } else {
        term
      }
      srt <- CellScoring(
        srt = srt,
        features = feature_list,
        method = method_i,
        classification = FALSE,
        layer = layer,
        assay = assay,
        name = assay_name_i,
        new_assay = TRUE,
        backend = backend_i,
        cpp_strategy = cpp_strategy,
        cores = cores,
        verbose = verbose,
        ...
      )
      srt <- RunDynamicFeatures(
        srt = srt,
        lineages = lineages,
        features = rownames(
          GetAssayData5(
            srt,
            assay = assay_name_i,
            layer = "counts"
          )
        ),
        suffix = paste(lineages, assay_name_i, sep = "_"),
        assay = assay_name_i,
        cores = cores,
        verbose = verbose
      )
    }
  }

  log_message(
    "Dynamic enrichment analysis completed",
    message_type = "success",
    verbose = verbose
  )

  return(srt)
}
