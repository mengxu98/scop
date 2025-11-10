#' @title Perform the enrichment analysis (over-representation) on the genes
#'
#' @md
#' @inheritParams GeneConvert
#' @param srt A Seurat object containing the results of differential expression analysis (RunDEtest).
#' If specified, the genes and groups will be extracted from the Seurat object automatically.
#' If not specified, the \code{geneID} and \code{geneID_groups} arguments must be provided.
#' @param group_by A character vector specifying the grouping variable in the Seurat object.
#' This argument is only used if \code{srt} is specified.
#' @param test.use A character vector specifying the test to be used in differential expression analysis.
#' This argument is only used if \code{srt} is specified.
#' @param DE_threshold A character vector specifying the filter condition for differential expression analysis.
#' This argument is only used if \code{srt} is specified.
#' @param geneID A character vector specifying the gene IDs.
#' @param geneID_groups A factor vector specifying the group labels for each gene.
#' @param geneID_exclude A character vector specifying the gene IDs to be excluded from the analysis.
#' @param IDtype A character vector specifying the type of gene IDs in the \code{srt} object or \code{geneID} argument.
#' This argument is used to convert the gene IDs to a different type if \code{IDtype} is different from \code{result_IDtype}.
#' @param result_IDtype A character vector specifying the desired type of gene ID to be used in the output.
#' This argument is used to convert the gene IDs from \code{IDtype} to \code{result_IDtype}.
#' @param species A character vector specifying the species for which the analysis is performed.
#' @param db A character vector specifying the name of the database to be used for enrichment analysis.
#' @param db_update Whether the gene annotation databases should be forcefully updated.
#' If set to FALSE, the function will attempt to load the cached databases instead.
#' Default is FALSE.
#' @param db_version A character vector specifying the version of the database to be used.
#' This argument is ignored if \code{db_update} is \code{TRUE}.
#' Default is "latest".
#' @param db_combine Whether to combine multiple databases into one.
#' If TRUE, all database specified by \code{db} will be combined as one named "Combined".
#' @param convert_species Whether to use a species-converted database when the annotation is missing for the specified species.
#' The default value is TRUE.
#' @param TERM2GENE A data frame specifying the gene-term mapping for a custom database.
#' The first column should contain the term IDs, and the second column should contain the gene IDs.
#' @param TERM2NAME A data frame specifying the term-name mapping for a custom database.
#' The first column should contain the term IDs, and the second column should contain the corresponding term names.
#' @param minGSSize The minimum size of a gene set to be considered in the enrichment analysis.
#' @param maxGSSize The maximum size of a gene set to be considered in the enrichment analysis.
#' @param unlimited_db A character vector specifying the names of databases that do not have size restrictions.
#' @param GO_simplify Whether to simplify the GO terms.
#' If \code{TRUE}, additional results with simplified GO terms will be returned.
#' @param GO_simplify_cutoff A character vector specifying the filter condition for simplification of GO terms.
#' This argument is only used if \code{GO_simplify} is \code{TRUE}.
#' @param simplify_method A character vector specifying the method to be used for simplification of GO terms.
#' This argument is only used if \code{GO_simplify} is \code{TRUE}.
#' @param simplify_similarityCutoff The similarity cutoff for simplification of GO terms.
#' This argument is only used if \code{GO_simplify} is \code{TRUE}.
#' @inheritParams thisutils::parallelize_fun
#'
#' @returns
#' If input is a Seurat object, returns the modified Seurat object with the enrichment result stored in the tools slot.
#'
#' If input is a geneID vector with or without geneID_groups, return the enrichment result directly.
#'
#' Enrichment result is a list with the following component:
#' \itemize{
#'  \item \code{enrichment}: A data.frame containing all enrichment results.
#'  \item \code{results}: A list of \code{enrichResult} objects from the DOSE package.
#'  \item \code{geneMap}: A data.frame containing the ID mapping table for input gene IDs.
#'  \item \code{input}: A data.frame containing the input gene IDs and gene ID groups.
#'  \item \code{DE_threshold}: A specific threshold for differential expression analysis (only returned if input is a Seurat object).
#' }
#'
#' @seealso
#' [PrepareDB], [ListDB], [EnrichmentPlot], [RunGSEA], [GSEAPlot]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunDEtest(
#'   pancreas_sub,
#'   group_by = "CellType"
#' )
#' pancreas_sub <- RunEnrichment(
#'   pancreas_sub,
#'   group_by = "CellType",
#'   DE_threshold = "p_val_adj < 0.05",
#'   db = "GO_BP",
#'   species = "Mus_musculus"
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   plot_type = "comparison"
#' )
#'
#' \dontrun{
#' pancreas_sub <- RunEnrichment(
#'   pancreas_sub,
#'   group_by = "CellType",
#'   DE_threshold = "p_val_adj < 0.05",
#'   db = c("MSigDB", "MSigDB_MH"),
#'   species = "Mus_musculus"
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "MSigDB",
#'   group_by = "CellType",
#'   plot_type = "comparison"
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "MSigDB_MH",
#'   group_by = "CellType",
#'   plot_type = "comparison"
#' )
#'
#' # Remove redundant GO terms
#' pancreas_sub <- RunEnrichment(
#'   pancreas_sub,
#'   group_by = "CellType",
#'   db = "GO_BP",
#'   GO_simplify = TRUE,
#'   species = "Mus_musculus"
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "GO_BP_sim",
#'   group_by = "CellType",
#'   plot_type = "comparison"
#' )
#'
#' # Or use "geneID" and "geneID_groups" as input to run enrichment
#' de_df <- dplyr::filter(
#'   pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox,
#'   p_val_adj < 0.05
#' )
#' enrich_out <- RunEnrichment(
#'   geneID = de_df[["gene"]],
#'   geneID_groups = de_df[["group1"]],
#'   db = "GO_BP",
#'   species = "Mus_musculus"
#' )
#' EnrichmentPlot(
#'   res = enrich_out,
#'   db = "GO_BP",
#'   plot_type = "comparison"
#' )
#'
#' # Use a combined database
#' pancreas_sub <- RunEnrichment(
#'   pancreas_sub,
#'   group_by = "CellType",
#'   db = c(
#'     "KEGG", "WikiPathway", "Reactome", "PFAM", "MP"
#'   ),
#'   db_combine = TRUE,
#'   species = "Mus_musculus"
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "Combined",
#'   group_by = "CellType",
#'   plot_type = "comparison"
#' )
#' }
RunEnrichment <- function(
    srt = NULL,
    group_by = NULL,
    test.use = "wilcox",
    DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
    geneID = NULL,
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
    verbose = TRUE) {
  log_message("Start {.pkg Enrichment} analysis", verbose = verbose)
  check_r("clusterProfiler")
  use_srt <- FALSE
  if (is.null(geneID)) {
    if (is.null(group_by)) {
      group_by <- "custom"
    }
    layer <- paste0("DEtest_", group_by)
    if (
      !layer %in% names(srt@tools) ||
        length(grep(pattern = "AllMarkers", names(srt@tools[[layer]]))) == 0
    ) {
      log_message(
        "Cannot find the DEtest result for the group '",
        group_by,
        "'. You may perform RunDEtest first.",
        message_type = "error"
      )
    }
    index <- grep(
      pattern = paste0("AllMarkers_", test.use),
      names(srt@tools[[layer]])
    )[1]
    if (is.na(index)) {
      log_message(
        "Cannot find the 'AllMarkers_", test.use, "' in the DEtest result.",
        message_type = "error"
      )
    }
    de <- names(srt@tools[[layer]])[index]
    de_df <- srt@tools[[layer]][[de]]
    de_df <- de_df[
      with(de_df, eval(rlang::parse_expr(DE_threshold))), ,
      drop = FALSE
    ]
    rownames(de_df) <- seq_len(nrow(de_df))

    geneID <- de_df[["gene"]]
    geneID_groups <- de_df[["group1"]]
    use_srt <- TRUE
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
      "length(geneID_groups)!=length(geneID)",
      message_type = "error"
    )
  }
  names(geneID_groups) <- geneID
  input <- data.frame(geneID = geneID, geneID_groups = geneID_groups)
  input <- input[!geneID %in% geneID_exclude, , drop = FALSE]

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
    db_list[[species]][[db]][["version"]] <- "custom"
  }
  if (isTRUE(db_combine)) {
    log_message("Create 'Combined' database ...")
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
    geneMap <- res$geneID_collapse
    colnames(geneMap)[colnames(geneMap) == "from_geneID"] <- IDtype
  } else {
    geneMap <- data.frame(IDtype = unique(geneID), row.names = unique(geneID))
    colnames(geneMap)[1] <- IDtype
  }

  input[[IDtype]] <- geneMap[as.character(input$geneID), IDtype]
  input[[result_IDtype]] <- geneMap[as.character(input$geneID), result_IDtype]
  input <- unnest_fun(input, cols = c(IDtype, result_IDtype))
  input <- input[!is.na(input[[IDtype]]), , drop = FALSE]

  log_message("Permform enrichment...")
  comb <- expand.grid(
    group = levels(geneID_groups),
    term = db,
    stringsAsFactors = FALSE
  )

  res_list <- parallelize_fun(
    seq_len(nrow(comb)),
    function(i) {
      group <- comb[i, "group"]
      term <- comb[i, "term"]
      gene <- input[input$geneID_groups == group, IDtype]
      gene_mapid <- input[input$geneID_groups == group, result_IDtype]
      TERM2GENE_tmp <- db_list[[species]][[term]][["TERM2GENE"]][, c(
        "Term",
        IDtype
      )]
      TERM2NAME_tmp <- db_list[[species]][[term]][["TERM2NAME"]]
      dup <- duplicated(TERM2GENE_tmp)
      na <- rowSums(is.na(TERM2GENE_tmp)) > 0
      TERM2GENE_tmp <- TERM2GENE_tmp[!(dup | na), , drop = FALSE]
      TERM2NAME_tmp <- TERM2NAME_tmp[
        TERM2NAME_tmp[["Term"]] %in% TERM2GENE_tmp[["Term"]], ,
        drop = FALSE
      ]
      enrich_res <- clusterProfiler::enricher(
        gene = gene,
        minGSSize = ifelse(term %in% unlimited_db, 1, minGSSize),
        maxGSSize = ifelse(term %in% unlimited_db, Inf, maxGSSize),
        pAdjustMethod = "BH",
        pvalueCutoff = Inf,
        qvalueCutoff = Inf,
        universe = NULL,
        TERM2GENE = TERM2GENE_tmp,
        TERM2NAME = TERM2NAME_tmp
      )

      if (!is.null(enrich_res) && nrow(enrich_res@result) > 0) {
        result <- enrich_res@result
        result[["Groups"]] <- group
        result[["Database"]] <- term
        result[["Version"]] <- as.character(db_list[[species]][[term]][[
          "version"
        ]])
        IDlist <- strsplit(result$geneID, split = "/")
        result$geneID <- unlist(lapply(IDlist, function(x) {
          x_result <- NULL
          for (i in x) {
            if (i %in% geneMap[[IDtype]]) {
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

        if (
          isTRUE(GO_simplify) && term %in% c("GO", "GO_BP", "GO_CC", "GO_MF")
        ) {
          sim_res <- enrich_res
          if (term == "GO") {
            sim_res@result[["ONTOLOGY"]] <- stats::setNames(
              TERM2NAME_tmp[["ONTOLOGY"]],
              TERM2NAME_tmp[["Term"]]
            )[sim_res@result[["ID"]]]
            sim_res@ontology <- "GOALL"
          } else {
            sim_res@ontology <- gsub(
              pattern = "GO_",
              replacement = "",
              x = term
            )
          }
          nterm_simplify <- sum(with(
            sim_res@result,
            eval(rlang::parse_expr(GO_simplify_cutoff))
          ))
          if (nterm_simplify <= 1) {
            log_message(
              group,
              "|",
              term,
              " has no term to simplify.",
              message_type = "warning"
            )
          } else {
            sim_res@result <- sim_res@result[
              with(sim_res@result, eval(rlang::parse_expr(GO_simplify_cutoff))), ,
              drop = FALSE
            ]
            semData <- db_list[[species]][[term]][["semData"]]
            sim_res <- clusterProfiler::simplify(
              sim_res,
              measure = simplify_method,
              cutoff = simplify_similarityCutoff,
              semData = semData
            )
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
  rownames(enrichment) <- NULL

  log_message(
    "{.pkg Enrichment} analysis done",
    message_type = "success", verbose = verbose
  )

  res <- list(
    enrichment = enrichment,
    results = results,
    geneMap = geneMap,
    input = input
  )
  if (isTRUE(use_srt)) {
    res[["DE_threshold"]] <- DE_threshold
    srt@tools[[paste("Enrichment", group_by, test.use, sep = "_")]] <- res
    return(srt)
  } else {
    return(res)
  }
}
