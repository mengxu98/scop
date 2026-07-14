#' @title Run metabolism pathway scoring
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams PrepareDB
#' @inheritParams RunEnrichment
#' @param assay Assay to use as expression matrix. Default is `DefaultAssay(srt)`.
#' @param layer Data layer to use, usually `"counts"` for count matrix.
#' @param db Databases to use for metabolism pathways. One or both of `"KEGG"`, `"REACTOME"`.
#' `"Reactome"` is also accepted and treated identically to `"REACTOME"`.
#' When `use_preparedb = TRUE`, gene sets are built via [PrepareDB].
#' @param method Scoring method, one of `"AUCell"`, `"GSVA"`, `"ssGSEA"`, `"VISION"`.
#' @param backend Scoring backend. `"cpp"` is the default for supported methods.
#' `"r"` uses the original R package implementation. `"cpp"` currently supports
#' `method = "AUCell"`, `method = "GSVA"`, and `method = "ssGSEA"`.
#' `method = "VISION"` falls back to `"r"` when `backend` is not explicitly set.
#' @param use_preparedb When `TRUE`, gene sets are built via [PrepareDB] which
#' provides species-aware gene mapping via BioMart and KEGG/Reactome databases.
#' This automatically handles gene symbol conversion for non-human species
#' (e.g., `species = "Mus_musculus"` → mouse gene symbols in metabolism pathways).
#' When `FALSE`, raw scMetabolism GMT files are
#' downloaded and genes are matched case-insensitively with optional [GeneConvert]
#' supplementation when `convert_species = TRUE`.
#' @param cpp_chunk_size Optional cell chunk size for C++ GSVA kernels. `NULL`
#' or `"auto"` automatically chunks large matrices to reduce peak dense
#' intermediate memory; positive values set the chunk size manually.
#' @param group.by Name of metadata column to group cells by. If `NULL`, single-cell scoring.
#' If provided, expression is averaged by group before scoring (cell-type level).
#' @param assay_name Name of the assay to store metabolism scores when `new_assay = TRUE`.
#' Default is `"METABOLISM"`.
#' @param new_assay Whether to create a new assay for metabolism scores when `group.by = NULL`. Default is `TRUE`.
#' @param species Species of the input data. The scMetabolism gene sets contain human
#' gene symbols. When `species` is not `"Homo_sapiens"` and `convert_species` is `TRUE`,
#' [GeneConvert] is used to map human genes to the target species via biomaRt homolog
#' tables. Default is `"Homo_sapiens"`.
#' @param convert_species Whether to convert human gene symbols from the scMetabolism
#' gene sets to the target species using [GeneConvert]. When `TRUE` (default), genes
#' are mapped via cross-species orthologs from Ensembl BioMart. When `FALSE`, only
#' case-insensitive direct symbol matching is used.
#' @param biomart BioMart database name passed to [GeneConvert]. Default `NULL` uses
#' `"ensembl"`. Other options: `"protists_mart"`, `"fungi_mart"`, `"plants_mart"`.
#' @param max_tries Maximum retry attempts for biomaRt connections in [GeneConvert].
#' Default is `5`.
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
  biomart = NULL, # deprecated, kept for compat
  max_tries = 5, # deprecated, kept for compat
  use_preparedb = TRUE,
  method = c("AUCell", "GSVA", "ssGSEA", "VISION"),
  backend = c("cpp", "r"),
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

  if (!isTRUE(use_preparedb) && !identical(tolower(IDtype), "symbol")) {
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

  # ---- Build gene sets ----

  if (isTRUE(use_preparedb)) {
    log_message(
      "Using {.fn PrepareDB} for species-aware gene set construction",
      verbose = verbose
    )
    curated <- scmetabolism_pathway_refs(db_prepare, verbose = verbose)
    log_message(
      "  KEGG pathway refs: {.val {length(curated[['kegg_refs']])}}, ",
      "Reactome pathway names: {.val {length(curated[['reactome_names']])}}",
      verbose = verbose
    )
    result <- build_metabolism_gene_sets_from_preparedb(
      species = species,
      db_prepare = db_prepare,
      IDtype = IDtype,
      curated = curated,
      expr_gene_names = rownames(expr_counts),
      db_update = db_update,
      db_version = db_version,
      convert_species = convert_species,
      Ensembl_version = Ensembl_version,
      mirror = mirror,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      verbose = verbose
    )
    gene_sets <- result[["gene_sets"]]
    term_names_final <- result[["term_names"]]

    # Skip the GMT-based gene set construction below
    skip_gmt <- TRUE
  } else {
    need_species_conv <- isTRUE(convert_species) && !identical(species, "Homo_sapiens")
    log_message(
      if (need_species_conv) paste0("Using raw {.pkg scMetabolism} gene sets with species conversion to {.val {species}}") else "Using raw {.pkg scMetabolism} gene sets directly (human gene symbols)",
      verbose = verbose
    )
    skip_gmt <- FALSE
  }

  gmt_sources <- scmetabolism_gmt_sources()

  gene_sets_all <- list()
  term_names_all <- list()
  metabolism_db_all <- list()
  expr_gene_lookup <- stats::setNames(
    rownames(expr_counts),
    toupper(rownames(expr_counts))
  )

  if (!isTRUE(skip_gmt)) {
    for (term_db in db_prepare) {
      term_db_label <- db_labels[[term_db]]
      metabolism_db <- load_scmetabolism_gmt(
        source = gmt_sources[[term_db]],
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
      metabolism_db_all[[term_db]] <- metabolism_db

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
  } # end if (!skip_gmt)

  if (!isTRUE(skip_gmt) && isTRUE(convert_species) && !identical(species, "Homo_sapiens") &&
    length(gene_sets_all) > 0) {
    all_human_genes <- unique(unlist(
      lapply(metabolism_db_all, function(mdb) {
        unique(unlist(mdb[["gene_sets"]], use.names = FALSE))
      }),
      use.names = FALSE
    ))
    if (length(all_human_genes) > 0) {
      log_message(
        "Converting {.val {length(all_human_genes)}} human gene symbols to {.val {species}} via {.pkg biomaRt} ...",
        verbose = verbose
      )
      conv <- tryCatch(
        GeneConvert(
          geneID = all_human_genes,
          geneID_from_IDtype = "symbol",
          geneID_to_IDtype = "symbol",
          species_from = "Homo_sapiens",
          species_to = species,
          Ensembl_version = Ensembl_version,
          biomart = biomart,
          mirror = mirror,
          max_tries = max_tries,
          verbose = verbose
        ),
        error = function(e) {
          log_message(
            "{.fn GeneConvert} failed: {.val {conditionMessage(e)}}. Falling back to direct symbol matching.",
            message_type = "warning",
            verbose = verbose
          )
          NULL
        }
      )
      if (!is.null(conv) && !is.null(conv[["geneID_expand"]]) &&
        nrow(conv[["geneID_expand"]]) > 0) {
        expand <- conv[["geneID_expand"]]
        from_col <- "from_geneID"
        to_col <- "symbol"
        hs_to_target <- split(
          expand[[to_col]],
          toupper(expand[[from_col]])
        )

        n_added_total <- 0L
        for (db_name in names(gene_sets_all)) {
          mdb <- metabolism_db_all[[db_name]]
          gmt_gene_sets <- mdb[["gene_sets"]]
          gs_current <- gene_sets_all[[db_name]]
          for (pw in names(gmt_gene_sets)) {
            if (!pw %in% names(gs_current)) next
            human_genes <- toupper(gmt_gene_sets[[pw]])
            # For each human gene, find target species homologs
            target_genes <- unique(unlist(
              hs_to_target[human_genes],
              use.names = FALSE
            ))
            target_genes <- target_genes[!is.na(target_genes) & nzchar(target_genes)]
            if (length(target_genes) == 0) next
            # Match target genes to expression rownames
            matched <- intersect(target_genes, rownames(expr_counts))
            if (length(matched) > 0) {
              existing <- gs_current[[pw]]
              gs_current[[pw]] <- unique(c(existing, matched))
              n_added_total <- n_added_total +
                (length(gs_current[[pw]]) - length(existing))
            }
          }
          # Re-filter by size
          gs_size <- lengths(gs_current)
          gs_current <- gs_current[gs_size >= minGSSize & gs_size <= maxGSSize]
          gene_sets_all[[db_name]] <- gs_current
        }
        log_message(
          "Species conversion added {.val {n_added_total}} gene matches across pathways",
          message_type = "success",
          verbose = verbose
        )
      }
    }
  }

  if (!isTRUE(skip_gmt)) {
    if (length(gene_sets_all) == 0) {
      log_message(
        "No metabolism gene sets were constructed from the specified databases",
        message_type = "error"
      )
    }

    gene_sets <- do.call(c, gene_sets_all)
    term_names_final <- unlist(term_names_all, use.names = TRUE)
  } # end if (!skip_gmt)

  log_message(
    "Total metabolism gene sets to score: {.val {length(gene_sets)}}",
    verbose = verbose
  )

  scores_mat <- NULL

  if (method == "AUCell") {
    if (identical(backend, "cpp")) {
      scores_mat <- run_aucell_scores(
        expr_counts = expr_counts,
        gene_sets = gene_sets,
        strategy = "topk",
        tie_method = "first"
      )
    } else {
      scores_mat <- run_aucell_official_scores(
        expr_counts = expr_counts,
        gene_sets = gene_sets,
        tie_method = "first"
      )
    }
  } else if (method %in% c("GSVA", "ssGSEA")) {
    if (identical(backend, "cpp")) {
      if (identical(method, "GSVA")) {
        scores_mat <- run_gsva_scores(
          expr_counts = expr_counts,
          gene_sets = gene_sets,
          min_gs_size = minGSSize,
          max_gs_size = maxGSSize,
          chunk_size = cpp_chunk_size
        )
      } else {
        scores_mat <- run_ssgsea_scores(
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
    create_gene_signature_fun <- get_namespace_fun("VISION", "createGeneSignature")
    signatures <- Map(
      function(sig_name, sig_genes) {
        create_gene_signature_fun(
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

  scores_for_plot <- as_matrix(Matrix::t(scores_mat))
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
        counts = scores_for_plot
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
        metadata = as.data.frame(Matrix::t(scores_for_plot))
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

scmetabolism_gmt_sources <- function() {
  urls <- c(
    KEGG = "https://raw.githubusercontent.com/mengxu98/datasets/main/scMetabolism/KEGG_metabolism_nc.gmt",
    Reactome = "https://raw.githubusercontent.com/mengxu98/datasets/main/scMetabolism/REACTOME_metabolism.gmt"
  )
  local <- c(
    KEGG = system.file(
      "data",
      "KEGG_metabolism_nc.gmt",
      package = "scMetabolism"
    ),
    Reactome = system.file(
      "data",
      "REACTOME_metabolism.gmt",
      package = "scMetabolism"
    )
  )
  out <- urls
  use_local <- file.exists(local) & nzchar(local)
  out[use_local] <- local[use_local]
  out
}

load_scmetabolism_gmt <- function(source, db_name, verbose = TRUE) {
  check_r("R.cache", verbose = FALSE)
  cache_key <- list("scmetabolism_raw_gmt", db_name, source)
  cached <- R.cache::loadCache(key = cache_key)
  if (!is.null(cached)) {
    return(cached)
  }

  tmp <- tempfile(fileext = ".gmt")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  if (file.exists(source)) {
    file.copy(source, tmp, overwrite = TRUE)
  } else {
    download(url = source, destfile = tmp, quiet = !verbose)
  }
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

scmetabolism_pathway_refs <- function(db_prepare, verbose = TRUE) {
  gmt_sources <- scmetabolism_gmt_sources()
  kegg_refs <- character(0)
  reactome_names <- character(0)

  for (term_db in intersect(db_prepare, names(gmt_sources))) {
    db <- load_scmetabolism_gmt(
      source = gmt_sources[[term_db]],
      db_name = term_db,
      verbose = verbose
    )
    if (is.null(db)) next
    refs <- as.character(db[["term_info"]][["Ref"]])
    names_db <- as.character(db[["term_info"]][["Name"]])
    if (identical(term_db, "KEGG")) {
      kegg_refs <- unique(stats::na.omit(refs[nzchar(trimws(refs))]))
    } else if (identical(term_db, "Reactome")) {
      reactome_names <- unique(stats::na.omit(names_db[nzchar(trimws(names_db))]))
    }
  }
  list(kegg_refs = kegg_refs, reactome_names = reactome_names)
}

build_metabolism_gene_sets_from_preparedb <- function(
  species,
  db_prepare,
  IDtype,
  curated,
  expr_gene_names,
  db_update,
  db_version,
  convert_species,
  Ensembl_version,
  mirror,
  minGSSize,
  maxGSSize,
  verbose
) {
  db_terms <- intersect(db_prepare, c("KEGG", "Reactome"))
  if (length(db_terms) == 0) {
    return(list(gene_sets = list(), term_names = character(0)))
  }

  db_list <- PrepareDB(
    species = species,
    db = db_terms,
    db_IDtypes = IDtype,
    db_version = db_version,
    db_update = db_update,
    convert_species = convert_species,
    Ensembl_version = Ensembl_version,
    mirror = mirror,
    verbose = verbose
  )

  gene_sets <- list()
  term_names <- character(0)
  expr_set <- stats::setNames(rep(TRUE, length(expr_gene_names)), expr_gene_names)

  for (term_db in db_terms) {
    db_entry <- db_list[[species]][[term_db]]
    if (is.null(db_entry)) next

    tg <- db_entry[["TERM2GENE"]]
    tn <- db_entry[["TERM2NAME"]]
    if (is.null(tg) || nrow(tg) == 0L) next

    # Filter to scMetabolism-curated pathways
    if (identical(term_db, "KEGG")) {
      # Match KEGG pathway numbers (e.g., "00010" matches "hsa00010")
      kegg_nums <- curated[["kegg_refs"]]
      if (length(kegg_nums) == 0) next
      term_ids <- as.character(tg[["Term"]])
      keep <- vapply(term_ids, function(tid) {
        any(vapply(kegg_nums, function(kn) grepl(paste0(kn, "$"), tid), logical(1)))
      }, logical(1))
      tg <- tg[keep, , drop = FALSE]
    } else if (identical(term_db, "Reactome")) {
      reactome_names_curated <- curated[["reactome_names"]]
      if (length(reactome_names_curated) == 0) next
      tn_sub <- tn[tn[["Term"]] %in% tg[["Term"]], , drop = FALSE]
      term_name_to_id <- stats::setNames(
        as.character(tn_sub[["Term"]]),
        tolower(as.character(tn_sub[["Name"]]))
      )
      matched_ids <- unique(unlist(
        lapply(tolower(reactome_names_curated), function(rn) {
          term_name_to_id[rn]
        }),
        use.names = FALSE
      ))
      matched_ids <- matched_ids[!is.na(matched_ids)]
      tg <- tg[tg[["Term"]] %in% matched_ids, , drop = FALSE]
    }

    if (nrow(tg) == 0) next

    gene_col <- IDtype
    if (!gene_col %in% colnames(tg)) {
      gene_col <- setdiff(colnames(tg), "Term")[1]
      if (is.na(gene_col)) next
    }

    # Build gene sets: Term → genes present in expression matrix
    term_vec <- as.character(tg[["Term"]])
    gene_vec <- as.character(tg[[gene_col]])
    valid_gene <- !is.na(gene_vec) & nzchar(trimws(gene_vec))
    term_vec <- term_vec[valid_gene]
    gene_vec <- gene_vec[valid_gene]

    gene_in_expr <- gene_vec %in% expr_gene_names
    term_vec <- term_vec[gene_in_expr]
    gene_vec <- gene_vec[gene_in_expr]

    if (length(term_vec) == 0) next

    gs <- split(gene_vec, term_vec)
    gs <- lapply(gs, unique)
    gs_size <- lengths(gs)
    gs <- gs[gs_size >= minGSSize & gs_size <= maxGSSize]

    if (length(gs) == 0) next

    # Build term names
    tn_sub <- tn[tn[["Term"]] %in% names(gs), , drop = FALSE]
    term_name_vec <- stats::setNames(
      as.character(tn_sub[["Name"]]),
      as.character(tn_sub[["Term"]])
    )

    gene_sets <- c(gene_sets, gs)
    term_names <- c(term_names, term_name_vec[names(gs)])

    log_message(
      "  {.val {term_db}}: {.val {length(gs)}} metabolism pathways, ",
      "{.val {length(unique(unlist(gs, use.names = FALSE)))}} genes mapped",
      verbose = verbose
    )
  }

  list(gene_sets = gene_sets, term_names = term_names)
}
