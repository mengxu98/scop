#' @title Perform the enrichment analysis (GSEA) on the genes
#'
#' @md
#' @inheritParams RunEnrichment
#' @param geneScore A numeric vector that specifies the gene scores,
#' for example, the log2(fold change) values of gene expression.
#' @param scoreType This parameter defines the GSEA score type.
#' Possible options are "std", "pos", "neg".
#' By default ("std") the enrichment score is computed as in the original GSEA.
#' The "pos" and "neg" score types are intended to be used for one-tailed tests
#' (i.e. when one is interested only in positive ("pos") or negateive ("neg") enrichment).
#'
#' @return
#' If input is a Seurat object, returns the modified Seurat object with the enrichment result stored in the tools slot.
#' If input is a geneID vector with or without geneID_groups, return the enrichment result directly.
#' Enrichment result is a list with the following component:
#' \itemize{
#'  \item `enrichment`: A data.frame containing all enrichment results.
#'  \item `results`: A list of `gseaResult` objects from the DOSE package.
#'  \item `geneMap`: A data.frame containing the ID mapping table for input gene IDs.
#'  \item `input`: A data.frame containing the input gene IDs and gene ID groups.
#'  \item `DE_threshold`: A specific threshold for differential expression analysis (only returned if input is a Seurat object).
#' }
#'
#' @seealso
#' [PrepareDB], [ListDB], [GSEAPlot], [RunEnrichment], [EnrichmentPlot]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunDEtest(
#'   pancreas_sub,
#'   group.by = "CellType"
#' )
#' pancreas_sub <- RunGSEA(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   DE_threshold = "p_val_adj < 0.05",
#'   scoreType = "std",
#'   db = "GO_BP",
#'   species = "Mus_musculus"
#' )
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group.by = "CellType",
#'   plot_type = "comparison"
#' )
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group.by = "CellType",
#'   group_use = "Ductal",
#'   id_use = "GO:0006412"
#' )
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group.by = "CellType",
#'   group_use = "Ductal",
#'   id_use = c(
#'     "GO:0046903", "GO:0015031", "GO:0007600"
#'   )
#' )
#'
#' # Remove redundant GO terms
#' pancreas_sub <- RunGSEA(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   db = "GO_BP",
#'   GO_simplify = TRUE,
#'   species = "Mus_musculus"
#' )
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "GO_BP_sim",
#'   group.by = "CellType",
#'   plot_type = "comparison"
#' )
#'
#' # Or use "geneID", "geneScore" and
#' # "geneID_groups" as input to run GSEA
#' de_df <- dplyr::filter(
#'   pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox,
#'   p_val_adj < 0.05
#' )
#' gsea_out <- RunGSEA(
#'   geneID = de_df[["gene"]],
#'   geneScore = de_df[["avg_log2FC"]],
#'   geneID_groups = de_df[["group1"]],
#'   db = "GO_BP",
#'   species = "Mus_musculus"
#' )
#' GSEAPlot(
#'   res = gsea_out,
#'   db = "GO_BP",
#'   plot_type = "comparison"
#' )
#'
#' # Use a combined database
#' pancreas_sub <- RunGSEA(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   db = c(
#'     "KEGG", "WikiPathway", "Reactome", "PFAM", "MP"
#'   ),
#'   db_combine = TRUE,
#'   species = "Mus_musculus"
#' )
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "Combined",
#'   group.by = "CellType",
#'   plot_type = "comparison"
#' )
RunGSEA <- function(
  srt = NULL,
  group.by = NULL,
  test.use = "wilcox",
  DE_threshold = "p_val_adj < 0.05",
  scoreType = "std",
  geneID = NULL,
  geneScore = NULL,
  geneID_groups = NULL,
  geneID_exclude = NULL,
  IDtype = "symbol",
  result_IDtype = "symbol",
  species = "Homo_sapiens",
  db = "GO_BP",
  db_update = FALSE,
  db_version = "latest",
  db_combine = FALSE,
  convert_species = TRUE,
  Ensembl_version = NULL,
  mirror = NULL,
  features = NULL,
  TERM2GENE = NULL,
  TERM2NAME = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  unlimited_db = c("Chromosome", "GeneType", "TF", "Enzyme", "CSPA"),
  GO_simplify = FALSE,
  GO_simplify_cutoff = "p.adjust < 0.05",
  simplify_method = "Wang",
  simplify_similarityCutoff = 0.7,
  cores = 1,
  verbose = TRUE,
  ...
) {
  log_message("Start {.pkg GSEA} analysis", verbose = verbose)

  use_object <- !is.null(srt)
  if (is.null(geneID) && !is.null(srt)) {
    de_df <- resolve_detest_result(
      object = srt,
      group.by = group.by,
      test.use = test.use
    )
    de_df <- filter_de_results(
      de_results = de_df,
      DE_threshold = DE_threshold
    )
    de_use <- prepare_de_for_pathway(
      de_results = de_df,
      require_score = TRUE
    )
    if (is.null(de_use) || nrow(de_use) == 0) {
      log_message(
        "Cannot find filtered DEtest results. You may perform {.fn RunDEtest} first or relax {.arg DE_threshold}.",
        message_type = "error"
      )
    }
    geneID <- de_use$gene
    geneScore <- de_use$avg_log2FC
    geneID_groups <- de_use$comparison
  }

  if (is.null(geneID_groups)) {
    geneID_groups <- rep(" ", length(geneID))
  }
  if (!is.factor(geneID_groups)) {
    geneID_groups <- factor(geneID_groups, levels = unique(geneID_groups))
  }
  geneID_groups <- factor(
    geneID_groups,
    levels = levels(geneID_groups)[levels(geneID_groups) %in% geneID_groups]
  )
  if (length(geneID_groups) != length(geneID)) {
    log_message(
      "{.arg geneID_groups} must be the same length with {.arg geneID}",
      message_type = "error"
    )
  }
  if (length(geneScore) != length(geneID)) {
    log_message(
      "{.arg geneScore} must be the same length with {.arg geneID}",
      message_type = "error"
    )
  }
  if (all(geneScore > 0) && scoreType != "pos") {
    scoreType <- "pos"
    log_message(
      "All values in the {.arg geneScore} are greater than zero. Set scoreType = 'pos'",
      message_type = "warning"
    )
  }
  if (all(geneScore < 0) && scoreType != "neg") {
    scoreType <- "neg"
    log_message(
      "All values in the {.arg geneScore} are less than zero. Set scoreType = 'neg'",
      message_type = "warning"
    )
  }

  input <- data.frame(
    geneID = geneID,
    geneScore = geneScore,
    geneID_groups = geneID_groups
  )
  input <- input[!geneID %in% geneID_exclude, , drop = FALSE]

  na_index <- which(is.na(geneScore))
  if (length(na_index) > 0) {
    log_message("Ignore {.val {length(na_index)}} NA {.arg geneScore}")
    input <- input[-na_index, , drop = FALSE]
  }
  input[
    is.infinite(input$geneScore) & input$geneScore < 0,
    "geneScore"
  ] <- min(input[!is.infinite(input$geneScore), "geneScore"])
  input[
    is.infinite(input$geneScore) & input$geneScore > 0,
    "geneScore"
  ] <- max(input[!is.infinite(input$geneScore), "geneScore"])

  geneID <- input$geneID
  geneScore <- input$geneScore
  geneID_groups <- input$geneID_groups
  names(geneID_groups) <- geneID
  names(geneScore) <- paste(geneID, geneID_groups, sep = ".")

  if (!is.null(features)) {
    db <- "custom"
    custom_db <- create_custom_db_list_from_features(
      species = species,
      db = db,
      features = features,
      IDtype = IDtype,
      version = "custom"
    )
    db_list <- custom_db[["db_list"]]
    TERM2GENE <- custom_db[["TERM2GENE"]]
    TERM2NAME <- custom_db[["TERM2NAME"]]
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
      IDtype = IDtype,
      version = "custom"
    )
    db_list <- custom_db[["db_list"]]
    TERM2GENE <- custom_db[["TERM2GENE"]]
    TERM2NAME <- custom_db[["TERM2NAME"]]
  }
  if (isTRUE(db_combine)) {
    log_message(
      "Create {.val Combined} database ...",
      verbose = verbose
    )
    TERM2GENE <- do.call(
      rbind,
      lapply(
        db_list[[species]],
        function(x) x[["TERM2GENE"]][, c("Term", IDtype)]
      )
    )
    TERM2NAME <- do.call(
      rbind,
      lapply(names(db_list[[species]]), function(x) {
        db_list[[species]][[x]][["TERM2NAME"]][["Name"]] <- paste0(
          db_list[[species]][[x]][["TERM2NAME"]][["Name"]],
          " [",
          x,
          "]"
        )
        db_list[[species]][[x]][["TERM2NAME"]][, c("Term", "Name")]
      })
    )
    version <- unlist(lapply(
      db_list[[species]],
      function(x) as.character(x[["version"]])
    ))
    version <- paste0(names(version), ":", version, collapse = ";")
    db <- "Combined"
    db_list[[species]][[db]][["TERM2GENE"]] <- unique(TERM2GENE)
    db_list[[species]][[db]][["TERM2NAME"]] <- unique(TERM2NAME)
    db_list[[species]][[db]][["version"]] <- unique(version)
  }

  if (length(unique(c(IDtype, result_IDtype))) != 1) {
    res <- GeneConvert(
      geneID = unique(geneID),
      geneID_from_IDtype = IDtype,
      geneID_to_IDtype = result_IDtype,
      species_from = species,
      species_to = species,
      Ensembl_version = Ensembl_version,
      mirror = mirror
    )
    geneMap <- res[["geneID_collapse"]]
    colnames(geneMap)[colnames(geneMap) == "from_geneID"] <- IDtype
  } else {
    geneMap <- data.frame(
      from_geneID = unique(geneID),
      row.names = unique(geneID)
    )
    colnames(geneMap)[1] <- IDtype
  }

  input[[IDtype]] <- geneMap[as.character(input$geneID), IDtype]
  input[[result_IDtype]] <- geneMap[as.character(input$geneID), result_IDtype]
  input <- unnest_fun(input, cols = c(IDtype, result_IDtype))
  input <- input[!is.na(input[[IDtype]]), , drop = FALSE]

  comb <- expand.grid(
    group = levels(geneID_groups),
    term = db,
    stringsAsFactors = FALSE
  )

  check_r("clusterProfiler", verbose = FALSE)
  res_list <- parallelize_fun(
    seq_len(nrow(comb)),
    function(i) {
      group <- comb[i, "group"]
      term <- comb[i, "term"]
      geneList <- input[input$geneID_groups == group, "geneScore"]
      names(geneList) <- input[input$geneID_groups == group, IDtype]
      gene_mapid <- input[input$geneID_groups == group, result_IDtype]
      ord <- order(geneList, decreasing = TRUE)
      geneList <- geneList[ord]
      gene_mapid <- gene_mapid[ord]
      TERM2GENE_tmp <- db_list[[species]][[term]][["TERM2GENE"]][, c(
        "Term",
        IDtype
      )]
      TERM2NAME_tmp <- db_list[[species]][[term]][["TERM2NAME"]]
      dup <- duplicated(TERM2GENE_tmp)
      na <- Matrix::rowSums(is.na(TERM2GENE_tmp)) > 0
      TERM2GENE_tmp <- TERM2GENE_tmp[!(dup | na), , drop = FALSE]
      TERM2NAME_tmp <- TERM2NAME_tmp[
        TERM2NAME_tmp[["Term"]] %in% TERM2GENE_tmp[["Term"]], ,
        drop = FALSE
      ]
      enrich_res <- clusterProfiler::GSEA(
        geneList = geneList,
        minGSSize = ifelse(term %in% unlimited_db, 1, minGSSize),
        maxGSSize = ifelse(term %in% unlimited_db, Inf, maxGSSize),
        nPermSimple = 1e5,
        eps = 0,
        scoreType = scoreType,
        pAdjustMethod = "BH",
        pvalueCutoff = Inf,
        TERM2GENE = TERM2GENE_tmp,
        TERM2NAME = TERM2NAME_tmp,
        by = "fgsea",
        verbose = FALSE
      )

      if (!is.null(enrich_res) && nrow(enrich_res@result) > 0) {
        result <- enrich_res@result
        result[["Groups"]] <- group
        result[["Database"]] <- term
        result[["Version"]] <- as.character(db_list[[species]][[term]][[
          "version"
        ]])
        IDlist <- strsplit(result$core_enrichment, "/")
        result$core_enrichment <- unlist(lapply(IDlist, function(x) {
          x_result <- NULL
          for (i in x) {
            if (i %in% input[[IDtype]]) {
              x_result <- c(
                x_result,
                unique(geneMap[geneMap[[IDtype]] == i, result_IDtype])
              )
            } else {
              x_result <- c(x_result, i)
            }
          }
          return(paste0(x_result, collapse = "/"))
        }))
        enrich_res@result <- result
        enrich_res@gene2Symbol <- as.character(gene_mapid)

        if (isTRUE(GO_simplify) && term %in% c("GO", "GO_BP", "GO_CC", "GO_MF")) {
          sim_res <- enrich_res
          if (term == "GO") {
            sim_res@result[["ONTOLOGY"]] <- stats::setNames(
              TERM2NAME_tmp[["ONTOLOGY"]],
              TERM2NAME_tmp[["Term"]]
            )[enrich_res@result[["ID"]]]
            sim_res@setType <- "GOALL"
          } else {
            sim_res@setType <- gsub(pattern = "GO_", replacement = "", x = term)
          }
          nterm_simplify <- sum(with(
            sim_res@result,
            eval(rlang::parse_expr(GO_simplify_cutoff))
          ))
          if (nterm_simplify > 0) {
            sim_res@result <- sim_res@result[
              with(
                sim_res@result, eval(rlang::parse_expr(GO_simplify_cutoff))
              ), ,
              drop = FALSE
            ]
          }
          if (nterm_simplify <= 1) {
            log_message(
              "{.pkg {term}} | {.val {group}} has no term to simplify",
              message_type = "warning",
              verbose = verbose
            )
          } else {
            sem_data <- db_list[[species]][[term]][["semData"]]
            sim_res <- clusterProfiler::simplify(
              sim_res,
              measure = simplify_method,
              cutoff = simplify_similarityCutoff,
              semData = sem_data
            )
          }
          if (nrow(sim_res@result) > 0) {
            result_sim <- sim_res@result
            result_sim[["Groups"]] <- group
            result_sim[["Database"]] <- paste0(term, "_sim")
            result_sim[["Version"]] <- as.character(db_list[[species]][[term]][[
              "version"
            ]])
            result_sim[["ONTOLOGY"]] <- NULL
            sim_res@result <- result_sim
            enrich_res <- list(enrich_res, sim_res)
            names(enrich_res) <- paste(
              group,
              c(term, paste0(term, "_sim")),
              sep = "-"
            )
          }
        }
      } else {
        enrich_res <- NULL
      }
      enrich_res
    },
    cores = cores,
    verbose = verbose
  )

  nm <- paste(comb$group, comb$term, sep = "-")
  sim_index <- sapply(res_list, function(x) length(x) == 2)
  sim_list <- unlist(res_list[sim_index], recursive = FALSE)
  raw_list <- res_list[!sim_index]
  names(raw_list) <- nm[!sim_index]
  results <- c(raw_list, sim_list)
  results <- results[!sapply(results, is.null)]
  results <- results[intersect(c(nm, paste0(nm, "_sim")), names(results))]
  enrichment <- do.call(rbind, lapply(results, function(x) x@result))
  if (is.null(enrichment)) {
    enrichment <- data.frame()
  }
  rownames(enrichment) <- NULL

  log_message(
    "{.pkg GSEA} analysis done",
    message_type = "success",
    verbose = verbose
  )

  res <- list(
    enrichment = enrichment,
    results = results,
    geneMap = geneMap,
    input = input
  )
  if (isTRUE(use_object) && inherits(srt, "Seurat")) {
    group.by <- group.by %||% "custom"
    result_key <- paste("GSEA", group.by, test.use, sep = "_")
    if (!is.null(srt@tools[[result_key]][["enrichment"]])) {
      old_enrichment <- srt@tools[[result_key]][["enrichment"]]
      old_results <- srt@tools[[result_key]][["results"]] %||% list()
      if (
        isTRUE(GO_simplify) &&
          all(c("Database", "Groups") %in% colnames(old_enrichment))
      ) {
        go_terms <- intersect(db, c("GO", "GO_BP", "GO_CC", "GO_MF"))
        sim_missing <- paste0(go_terms, "_sim")
        sim_missing <- sim_missing[
          !sim_missing %in% unique(c(enrichment[["Database"]], old_enrichment[["Database"]]))
        ]
        if (length(sim_missing) > 0) {
          sim_rows <- old_enrichment[
            old_enrichment[["Database"]] %in% sub("_sim$", "", sim_missing), ,
            drop = FALSE
          ]
          if (nrow(sim_rows) > 0) {
            sim_rows[["Database"]] <- paste0(sim_rows[["Database"]], "_sim")
            enrichment <- dplyr::bind_rows(enrichment, sim_rows)
            for (sim_db in unique(sim_rows[["Database"]])) {
              raw_db <- sub("_sim$", "", sim_db)
              sim_groups <- unique(sim_rows[sim_rows[["Database"]] == sim_db, "Groups"])
              for (sim_group in sim_groups) {
                raw_name <- paste(sim_group, raw_db, sep = "-")
                sim_name <- paste(sim_group, sim_db, sep = "-")
                raw_res <- res[["results"]][[raw_name]]
                if (is.null(raw_res)) {
                  raw_res <- old_results[[raw_name]]
                }
                if (!is.null(raw_res)) {
                  sim_res <- raw_res
                  sim_result <- sim_res@result
                  sim_result[["Database"]] <- sim_db
                  sim_res@result <- sim_result
                  res[["results"]][[sim_name]] <- sim_res
                }
              }
            }
          }
        }
      }
      enrichment <- unique(dplyr::bind_rows(old_enrichment, enrichment))
      res[["enrichment"]] <- enrichment
      if (!is.null(srt@tools[[result_key]][["results"]])) {
        merged_results <- old_results
        merged_results[names(res[["results"]])] <- res[["results"]]
        res[["results"]] <- merged_results
      }
    }
    srt@tools[[result_key]] <- utils::modifyList(res, list(DE_threshold = DE_threshold))
    return(srt)
  }
  if (isTRUE(use_object) && inherits(srt, "SummarizedExperiment")) {
    return(
      store_meta(
        srt,
        "GSEA",
        utils::modifyList(res, list(DE_threshold = DE_threshold))
      )
    )
  } else {
    return(res)
  }
}
