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
#' @param backend Scoring backend. `"cpp"` is the default for supported methods.
#' `"r"` uses the original R package implementation. `"cpp"` currently supports
#' `method = "AUCell"`, `method = "GSVA"`, and `method = "ssGSEA"`.
#' `method = "VISION"` falls back to `"r"` when `backend` is not explicitly set.
#' AUCell C++ scores may differ from the R backend when tied expression values
#' are randomly ranked.
#' @param cpp_strategy C++ AUCell ranking strategy. `"sparse"` ranks non-zero
#' genes and approximates zero ties, `"topk"` ranks only genes that can contribute
#' to AUCell AUC, and `"full"` ranks all genes.
#' @param cpp_chunk_size Optional cell chunk size for C++ GSVA kernels. `NULL`
#' or `"auto"` automatically chunks large matrices to reduce peak dense
#' intermediate memory; positive values set the chunk size manually.
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
  backend = c("cpp", "r"),
  cpp_strategy = c("sparse", "topk", "full"),
  cpp_chunk_size = NULL,
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

  backend_missing <- missing(backend)
  method <- match.arg(method)
  backend <- match.arg(backend)
  cpp_strategy <- match.arg(cpp_strategy)
  if (!identical(backend, "r") && !method %in% c("AUCell", "GSVA", "ssGSEA")) {
    if (isTRUE(backend_missing)) {
      log_message(
        "{.arg method = 'VISION'} does not have a C++ backend yet; using {.arg backend = 'r'} for this run.",
        message_type = "warning",
        verbose = verbose
      )
      backend <- "r"
    } else {
      log_message(
        "{.arg backend = 'cpp'} currently supports {.arg method = 'AUCell'}, {.arg method = 'GSVA'}, and {.arg method = 'ssGSEA'} only",
        message_type = "error"
      )
    }
  }
  db <- intersect(toupper(db), c("KEGG", "REACTOME"))
  if (length(db) == 0) {
    log_message(
      "{.arg db} must contain at least one of {.val KEGG} or {.val REACTOME}",
      message_type = "error"
    )
  }
  db_prepare <- unname(c("KEGG" = "KEGG", "REACTOME" = "Reactome")[db])
  db_labels <- stats::setNames(db, db_prepare)

  if (!identical(tolower(IDtype), "symbol")) {
    log_message(
      paste0(
        "{.fn RunMetabolism} now uses raw {.pkg scMetabolism} gene sets ",
        "directly, so {.arg IDtype} must currently be {.val symbol}"
      ),
      message_type = "error"
    )
  }

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

  log_message(
    paste0(
      "Using raw {.pkg scMetabolism} gene sets directly; ",
      "{.fn PrepareDB} / BioMart-based ID rebuilding is skipped"
    ),
    verbose = verbose
  )

  gmt_urls <- c(
    KEGG = "https://raw.githubusercontent.com/mengxu98/datasets/main/scMetabolism/KEGG_metabolism_nc.gmt",
    Reactome = "https://raw.githubusercontent.com/mengxu98/datasets/main/scMetabolism/REACTOME_metabolism.gmt"
  )

  gene_sets_all <- list()
  term_names_all <- list()
  expr_gene_lookup <- stats::setNames(
    rownames(expr_counts),
    toupper(rownames(expr_counts))
  )

  for (term_db in db_prepare) {
    term_db_label <- db_labels[[term_db]]
    metabolism_db <- load_scmetabolism_gmt(
      url = gmt_urls[[term_db]],
      db_name = term_db,
      verbose = verbose
    )
    if (is.null(metabolism_db)) {
      log_message(
        "Failed to load raw {.pkg scMetabolism} genesets for {.val {term_db_label}}, skip this database.",
        message_type = "warning",
        verbose = verbose
      )
      next
    }

    gene_sets_db <- lapply(
      metabolism_db[["gene_sets"]],
      function(gs) {
        mapped <- unname(expr_gene_lookup[toupper(gs)])
        unique(stats::na.omit(mapped[nzchar(mapped)]))
      }
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
    term_names_db <- metabolism_db[["term_info"]][term_ids_db, "Name"]
    term_names_db[is.na(term_names_db)] <- term_ids_db[is.na(term_names_db)]

    if (length(gene_sets_db) == 0) {
      next
    }

    gene_sets_all[[term_db]] <- gene_sets_db
    term_names_all[[term_db]] <- stats::setNames(
      term_names_db,
      term_ids_db
    )
  }

  if (length(gene_sets_all) == 0) {
    log_message(
      "No metabolism gene sets were constructed from the specified databases",
      message_type = "error"
    )
  }

  gene_sets <- do.call(c, gene_sets_all)
  term_names_final <- unlist(term_names_all, use.names = TRUE)

  log_message(
    "Total metabolism gene sets to score: {.val {length(gene_sets)}}",
    verbose = verbose
  )

  scores_mat <- NULL

  if (method == "AUCell") {
    if (identical(backend, "cpp")) {
      scores_mat <- run_aucell_cpp_scores(
        expr_counts = expr_counts,
        gene_sets = gene_sets,
        strategy = cpp_strategy
      )
    } else {
      gene_set_scoring_require_namespace("AUCell")
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
    }
  } else if (method %in% c("GSVA", "ssGSEA")) {
    if (identical(backend, "cpp")) {
      if (identical(method, "GSVA")) {
        scores_mat <- run_gsva_cpp_scores(
          expr_counts = expr_counts,
          gene_sets = gene_sets,
          min_gs_size = minGSSize,
          max_gs_size = maxGSSize,
          chunk_size = cpp_chunk_size
        )
      } else {
        scores_mat <- run_ssgsea_cpp_scores(
          expr_counts = expr_counts,
          gene_sets = gene_sets,
          min_gs_size = minGSSize,
          max_gs_size = maxGSSize
        )
      }
    } else {
      gene_set_scoring_require_namespace("GSVA")
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
    }
  } else if (method == "VISION") {
    gene_set_scoring_require_namespace("VISION", install_hint = "YosefLab/VISION")
    Vision_fun <- get_namespace_fun("VISION", "Vision")
    calc_signature_scores_fun <- get_namespace_fun("VISION", "calcSignatureScores")

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
    signatures <- Map(
      function(sig_name, sig_genes) {
        VISION::createGeneSignature(
          sig_name,
          stats::setNames(rep(1, length(sig_genes)), sig_genes)
        )
      },
      names(signatures),
      signatures
    )

    vis <- Vision_fun(
      scaled_counts,
      signatures = signatures,
      projection_methods = character()
    )
    vis <- calc_signature_scores_fun(
      vis,
      sig_gene_importance = FALSE
    )
    sig_scores <- vis@SigScores
    scores_mat <- as_matrix(sig_scores)
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
  term_labels <- term_names_final[term_ids]
  term_labels[is.na(term_labels) | !nzchar(trimws(term_labels))] <- term_ids[is.na(term_labels) | !nzchar(trimws(term_labels))]
  term_names <- capitalize(
    trimws(as.character(term_labels)),
    force_tolower = TRUE
  )
  rownames(scores_for_plot) <- make.unique(term_names)

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
    backend = backend,
    cpp_strategy = if (identical(backend, "cpp") && identical(method, "AUCell")) cpp_strategy else NA_character_,
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

load_scmetabolism_gmt <- function(url, db_name, verbose = TRUE) {
  check_r("R.cache", verbose = FALSE)
  cache_key <- list("scmetabolism_raw_gmt", db_name, url)
  cached <- R.cache::loadCache(key = cache_key)
  if (!is.null(cached)) {
    return(cached)
  }

  tmp <- tempfile(fileext = ".gmt")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  download(url = url, destfile = tmp, quiet = !verbose)
  lines <- readLines(tmp, warn = FALSE)
  split_lines <- strsplit(lines, "\t", fixed = TRUE)

  term_name_raw <- vapply(
    split_lines,
    function(x) trimws(x[1]),
    FUN.VALUE = character(1)
  )
  term_ref <- vapply(
    split_lines,
    function(x) {
      if (length(x) >= 2) trimws(x[2]) else NA_character_
    },
    FUN.VALUE = character(1)
  )
  term_ids <- paste0(db_name, ".", make.unique(term_name_raw))
  gene_sets <- stats::setNames(
    lapply(split_lines, function(x) {
      genes <- trimws(x[-c(1, 2)])
      genes <- genes[nzchar(genes)]
      unique(stats::na.omit(genes))
    }),
    term_ids
  )

  out <- list(
    gene_sets = gene_sets,
    term_info = data.frame(
      Term = term_ids,
      Name = term_name_raw,
      Ref = term_ref,
      stringsAsFactors = FALSE,
      row.names = term_ids
    )
  )
  R.cache::saveCache(out, key = cache_key)
  out
}
