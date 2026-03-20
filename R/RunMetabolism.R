#' @title Run metabolism pathway scoring
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams PrepareDB
#' @inheritParams RunEnrichment
#' @param assay Assay to use as expression matrix. Default is `DefaultAssay(srt)`.
#' @param layer Data layer to use, usually `"counts"` for count matrix.
#' @param db Databases to use for metabolism pathways. One or both of `"KEGG"`, `"REACTOME"`.
#' @param method Scoring method, one of `"AUCell"`, `"GSVA"`, `"ssGSEA"`, `"VISION"`.
#' @param group.by Name of metadata column to group cells by. If `NULL`, single-cell scoring.
#' If provided, expression is averaged by group before scoring (cell-type level).
#' @param assay_name Name of the assay to store metabolism scores when `new_assay = TRUE`.
#' Default is `"METABOLISM"`.
#' @param new_assay Whether to create a new assay for metabolism scores when `group.by = NULL`. Default is `TRUE`.
#'
#' @return
#' Returns a `Seurat` object. When `group.by = NULL`, stores scores in assay `assay_name` and tools.
#' When `group.by` is provided, stores in tools slot `Metabolism_<group.by>_<method>` for [MetabolismPlot].
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunMetabolism(
#'   pancreas_sub,
#'   assay = "RNA",
#'   layer = "counts",
#'   db = c("KEGG", "REACTOME"),
#'   group.by = "CellType",
#'   species = "Mus_musculus",
#'   method = "AUCell"
#' )
#' ht <- MetabolismPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "heatmap",
#'   topTerm = 10,
#'   width = 1,
#'   height = 2
#' )
RunMetabolism <- function(
  srt,
  assay = NULL,
  group.by = NULL,
  layer = "counts",
  db = c("KEGG", "REACTOME"),
  species = "Homo_sapiens",
  IDtype = "symbol",
  db_update = FALSE,
  db_version = "latest",
  convert_species = TRUE,
  Ensembl_version = NULL,
  mirror = NULL,
  method = c("AUCell", "GSVA", "ssGSEA", "VISION"),
  minGSSize = 10,
  maxGSSize = 500,
  assay_name = "METABOLISM",
  new_assay = TRUE,
  seed = 11,
  verbose = TRUE
) {
  log_message(
    "Start {.pkg metabolism pathway} scoring",
    verbose = verbose
  )
  set.seed(seed)

  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (layer == "counts") {
    status <- CheckDataType(srt, layer = "counts", assay = assay)
    if (status != "raw_counts") {
      log_message(
        "Layer {.val counts} is not raw counts in assay {.val {assay}}",
        message_type = "warning",
        verbose = verbose
      )
    }
  }

  method <- match.arg(method)
  db <- intersect(toupper(db), c("KEGG", "REACTOME"))
  if (length(db) == 0) {
    log_message(
      "{.arg db} must contain at least one of {.val KEGG} or {.val REACTOME}",
      message_type = "error"
    )
  }
  db_prepare <- unname(c("KEGG" = "KEGG", "REACTOME" = "Reactome")[db])
  db_labels <- stats::setNames(db, db_prepare)

  if (isTRUE(new_assay) && is.null(assay_name) && is.null(group.by)) {
    log_message(
      "{.arg assay_name} must be specified when {.arg new_assay = TRUE} and {.arg group.by = NULL}",
      message_type = "error"
    )
  }

  single_cell_mode <- is.null(group.by)
  if (single_cell_mode) {
    expr_counts <- GetAssayData5(srt, layer = layer, assay = assay)
  } else {
    if (length(group.by) != 1 || !group.by %in% colnames(srt@meta.data)) {
      log_message(
        "{.arg group.by} must be a single metadata column present in meta.data",
        message_type = "error"
      )
    }
    log_message(
      "Averaging expression by {.val {group.by}} ...",
      verbose = verbose
    )
    features <- rownames(srt[[assay]])
    expr_avg <- Seurat::AggregateExpression(
      object = srt,
      features = features,
      assays = assay,
      group.by = group.by,
      verbose = FALSE
    )
    expr_counts <- if (assay %in% names(expr_avg)) expr_avg[[assay]] else expr_avg[[1]]
    log_message(
      "Aggregated expression: {.val {nrow(expr_counts)}} genes x {.val {ncol(expr_counts)}} groups",
      verbose = verbose
    )
  }
  if (!inherits(expr_counts, c("dgCMatrix", "matrix", "Matrix"))) {
    expr_counts <- as_matrix(expr_counts)
  }

  db_list <- PrepareDB(
    species = species,
    db = db_prepare,
    db_update = db_update,
    db_version = db_version,
    db_IDtypes = IDtype,
    convert_species = convert_species,
    Ensembl_version = Ensembl_version,
    mirror = mirror
  )

  check_r("R.cache", verbose = FALSE)
  gmt_urls <- c(
    KEGG = "https://raw.githubusercontent.com/mengxu98/datasets/main/scMetabolism/KEGG_metabolism_nc.gmt",
    Reactome = "https://raw.githubusercontent.com/mengxu98/datasets/main/scMetabolism/REACTOME_metabolism.gmt"
  )
  metabolism_term_defs <- R.cache::loadCache(key = list("metabolism_term_defs", gmt_urls))
  if (is.null(metabolism_term_defs)) {
    metabolism_term_defs <- list()
    for (term_db in names(gmt_urls)) {
      tmp <- tempfile(fileext = ".gmt")
      on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
      thisutils::download(url = gmt_urls[[term_db]], destfile = tmp, quiet = !verbose)
      lines <- readLines(tmp, warn = FALSE)
      unlink(tmp)
      split_lines <- strsplit(lines, "\t", fixed = TRUE)
      Name <- vapply(split_lines, function(x) x[1], FUN.VALUE = character(1))
      Ref <- vapply(split_lines, function(x) if (length(x) >= 2) x[2] else NA_character_, FUN.VALUE = character(1))
      metabolism_term_defs[[term_db]] <- data.frame(Name = Name, Ref = Ref, stringsAsFactors = FALSE)
    }
    R.cache::saveCache(metabolism_term_defs, key = list("metabolism_term_defs", gmt_urls))
  }

  gene_sets_all <- list()
  term_names_all <- list()

  for (term_db in db_prepare) {
    term_db_label <- db_labels[[term_db]]
    TERM2GENE_raw <- db_list[[species]][[term_db]][["TERM2GENE"]]
    TERM2NAME_raw <- db_list[[species]][[term_db]][["TERM2NAME"]]

    if (is.null(TERM2GENE_raw) || is.null(TERM2NAME_raw)) {
      log_message(
        "Database {.val {term_db_label}} does not contain valid {.field TERM2GENE} / {.field TERM2NAME}, skip this database for metabolism scoring.",
        message_type = "warning",
        verbose = verbose
      )
      next
    }

    TERM2GENE_tmp <- as.data.frame(TERM2GENE_raw)[, c("Term", IDtype), drop = FALSE]
    TERM2NAME_tmp <- as.data.frame(TERM2NAME_raw)

    dup <- duplicated(TERM2GENE_tmp)
    na <- Matrix::rowSums(is.na(as.matrix(TERM2GENE_tmp))) > 0
    TERM2GENE_tmp <- TERM2GENE_tmp[!(dup | na), , drop = FALSE]
    TERM2NAME_tmp <- TERM2NAME_tmp[
      TERM2NAME_tmp[["Term"]] %in% TERM2GENE_tmp[["Term"]], ,
      drop = FALSE
    ]

    rownames(TERM2NAME_tmp) <- TERM2NAME_tmp[["Term"]]

    term_ids <- unique(TERM2GENE_tmp[["Term"]])
    if (term_db %in% names(metabolism_term_defs)) {
      metabolism_df <- metabolism_term_defs[[term_db]]
      term_name_norm <- thisutils::capitalize(trimws(as.character(TERM2NAME_tmp[["Name"]])), force_tolower = TRUE)
      metabolism_name_norm <- thisutils::capitalize(trimws(as.character(metabolism_df[["Name"]])), force_tolower = TRUE)
      matched_terms <- TERM2NAME_tmp[["Term"]][term_name_norm %in% metabolism_name_norm]

      if (identical(term_db, "KEGG")) {
        kegg_term_suffix <- gsub("^[[:alpha:]]+", "", TERM2NAME_tmp[["Term"]])
        matched_terms <- unique(c(
          matched_terms,
          TERM2NAME_tmp[["Term"]][kegg_term_suffix %in% metabolism_df[["Ref"]]]
        ))
      }

      term_ids <- intersect(term_ids, unique(stats::na.omit(matched_terms)))
      if (length(term_ids) == 0) {
        log_message(
          "No metabolism terms found in {.val {term_db_label}} after filtering by metabolism pathway definitions",
          message_type = "warning",
          verbose = verbose
        )
        next
      }
      TERM2GENE_tmp <- TERM2GENE_tmp[
        TERM2GENE_tmp[["Term"]] %in% term_ids, ,
        drop = FALSE
      ]
      TERM2NAME_tmp <- TERM2NAME_tmp[
        TERM2NAME_tmp[["Term"]] %in% term_ids, ,
        drop = FALSE
      ]
    }

    gene_sets_db <- split(
      TERM2GENE_tmp[[IDtype]],
      TERM2GENE_tmp[["Term"]]
    )
    gs_size <- lengths(gene_sets_db)
    gene_sets_db <- gene_sets_db[
      gs_size >= minGSSize & gs_size <= maxGSSize
    ]
    if (length(gene_sets_db) == 0) {
      log_message(
        "No metabolism gene sets remain in {.val {term_db_label}} after size filtering",
        message_type = "warning",
        verbose = verbose
      )
      next
    }

    term_ids_db <- names(gene_sets_db)
    term_names_db <- TERM2NAME_tmp[term_ids_db, "Name"]
    term_names_db[is.na(term_names_db)] <- term_ids_db[is.na(term_names_db)]
    term_names_db <- gsub(
      pattern = "^(KEGG|REACTOME|Reactome)\\.",
      replacement = "",
      x = term_names_db
    )

    genes_in_sets <- unique(unlist(gene_sets_db))
    overlap <- intersect(rownames(expr_counts), genes_in_sets)
    if (length(overlap) == 0) {
      log_message(
        "No overlapping genes found between {.val {term_db_label}} metabolism sets and assay {.val {assay}}",
        message_type = "warning",
        verbose = verbose
      )
      next
    }
    gene_sets_db <- lapply(
      gene_sets_db,
      function(gs) intersect(gs, overlap)
    )
    keep <- lengths(gene_sets_db) >= minGSSize
    gene_sets_db <- gene_sets_db[keep]
    term_names_db <- term_names_db[keep]

    if (length(gene_sets_db) == 0) {
      next
    }

    gene_sets_all[[term_db]] <- stats::setNames(
      gene_sets_db,
      term_names_db
    )
    term_names_all[[term_db]] <- term_names_db
  }

  if (length(gene_sets_all) == 0) {
    log_message(
      "No metabolism gene sets were constructed from the specified databases",
      message_type = "error"
    )
  }

  gene_sets <- do.call(c, gene_sets_all)
  names(gene_sets) <- gsub(
    pattern = "^(KEGG|REACTOME|Reactome)\\.",
    replacement = "",
    x = names(gene_sets)
  )
  names(gene_sets) <- make.unique(names(gene_sets))
  term_names_final <- unlist(term_names_all, use.names = FALSE)

  log_message(
    "Total metabolism gene sets to score: {.val {length(gene_sets)}}",
    verbose = verbose
  )

  scores_mat <- NULL

  if (method == "AUCell") {
    check_r("AUCell", verbose = FALSE)
    expr_rank <- AUCell::AUCell_buildRankings(
      as_matrix(expr_counts),
      plotStats = FALSE
    )
    cells_auc <- AUCell::AUCell_calcAUC(
      geneSets = gene_sets,
      rankings = expr_rank
    )
    auc_mat <- AUCell::getAUC(cells_auc)
    scores_mat <- Matrix::t(auc_mat)
    scores_mat <- as_matrix(scores_mat)
    colnames(scores_mat) <- rownames(auc_mat)
  } else if (method %in% c("GSVA", "ssGSEA")) {
    check_r("GSVA", verbose = FALSE)
    expr_mat <- as_matrix(expr_counts)
    expr_mat <- expr_mat[rowSums(expr_mat) > 0, , drop = FALSE]
    gene_sets_filt <- lapply(gene_sets, function(gs) intersect(gs, rownames(expr_mat)))
    gene_sets_filt <- gene_sets_filt[lengths(gene_sets_filt) >= minGSSize]
    if (identical(method, "ssGSEA")) {
      param <- GSVA::ssgseaParam(
        exprData = expr_mat,
        geneSets = gene_sets_filt,
        minSize = minGSSize,
        maxSize = maxGSSize
      )
    } else {
      param <- GSVA::gsvaParam(
        exprData = expr_mat,
        geneSets = gene_sets_filt,
        minSize = minGSSize,
        maxSize = maxGSSize,
        kcdf = "Poisson"
      )
    }
    gsva_es <- GSVA::gsva(param = param, verbose = verbose)
    if (inherits(gsva_es, "SummarizedExperiment")) {
      gsva_es <- SummarizedExperiment::assay(gsva_es)
    }
    gsva_es <- as.matrix(gsva_es)
    scores_mat <- Matrix::t(gsva_es)
  } else if (method == "VISION") {
    check_r("YosefLab/VISION", verbose = FALSE)
    Vision_fun <- thisutils::get_namespace_fun("VISION", "Vision")
    analyze_fun <- thisutils::get_namespace_fun("VISION", "analyze")

    n_umi <- Matrix::colSums(expr_counts)
    scaled_counts <- Matrix::t(
      Matrix::t(as_matrix(expr_counts)) / n_umi
    ) * stats::median(n_umi)

    signatures <- lapply(gene_sets, function(gs) {
      intersect(gs, rownames(expr_counts))
    })
    signatures <- signatures[lengths(signatures) > 0]
    if (length(signatures) == 0) {
      log_message(
        "No metabolism gene sets retain genes after intersecting with counts for {.pkg VISION}",
        message_type = "error"
      )
    }

    vis <- Vision_fun(
      scaled_counts,
      signatures = signatures
    )
    vis <- analyze_fun(vis)
    sig_scores <- vis@SigScores
    scores_mat <- Matrix::t(as_matrix(sig_scores))
  }

  if (is.null(scores_mat)) {
    log_message(
      "Failed to compute metabolism scores for method {.val {method}}",
      message_type = "error"
    )
  }

  if (single_cell_mode) {
    scores_mat <- scores_mat[intersect(rownames(scores_mat), colnames(srt)), , drop = FALSE]
    scores_mat <- scores_mat[colnames(srt), , drop = FALSE]
  }

  scores_for_plot <- Matrix::t(as_matrix(scores_mat))
  term_ids <- rownames(scores_for_plot)
  term_names <- thisutils::capitalize(
    trimws(as.character(term_ids)),
    force_tolower = TRUE
  )
  rownames(scores_for_plot) <- term_names

  if (single_cell_mode) {
    if (isTRUE(new_assay)) {
      srt[[assay_name]] <- Seurat::CreateAssayObject(
        counts = as_matrix(scores_for_plot)
      )
      srt[[assay_name]] <- Seurat::AddMetaData(
        object = srt[[assay_name]],
        metadata = data.frame(
          termnames = rownames(scores_for_plot),
          row.names = rownames(scores_for_plot)
        )
      )
    } else {
      srt <- Seurat::AddMetaData(
        object = srt,
        metadata = as.data.frame(Matrix::t(as_matrix(scores_for_plot)))
      )
    }
  }

  if (single_cell_mode) {
    enrichment <- NULL
  } else {
    term_gene_ids <- vapply(
      term_ids,
      function(x) {
        gs <- gene_sets[[x]]
        if (is.null(gs)) NA_character_ else paste0(gs, collapse = "/")
      },
      FUN.VALUE = character(1)
    )
    enrichment <- do.call(rbind, lapply(colnames(scores_for_plot), function(grp) {
      data.frame(
        ID = term_ids,
        Description = rownames(scores_for_plot),
        geneID = term_gene_ids,
        Groups = grp,
        Database = "Metabolism",
        GSVA_Score = as.numeric(scores_for_plot[, grp]),
        stringsAsFactors = FALSE
      )
    }))
    rownames(enrichment) <- NULL
  }

  tool_name <- if (single_cell_mode) {
    paste0("Metabolism_", method)
  } else {
    paste0("Metabolism_", group.by, "_", method)
  }
  srt@tools[[tool_name]] <- list(
    scores = scores_for_plot,
    enrichment = enrichment,
    gene_sets = gene_sets,
    method = method,
    db = db,
    group.by = group.by,
    species = species,
    assay = assay,
    layer = layer,
    db_version = db_version,
    db_update = db_update,
    convert_species = convert_species,
    Ensembl_version = Ensembl_version,
    mirror = mirror
  )

  log_message(
    "Metabolism scores stored in tools slot {.val {tool_name}}",
    message_type = "success",
    verbose = verbose
  )
  if (single_cell_mode && isTRUE(new_assay)) {
    log_message(
      "Metabolism scores also stored in assay {.val {assay_name}}",
      verbose = verbose
    )
  }

  return(srt)
}
