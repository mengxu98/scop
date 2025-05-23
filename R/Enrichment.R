#' Perform the enrichment analysis (over-representation) on the genes
#'
#' @inheritParams GeneConvert
#' @param srt A Seurat object containing the results of differential expression analysis (RunDEtest).
#' If specified, the genes and groups will be extracted from the Seurat object automatically.
#' If not specified, the \code{geneID} and \code{geneID_groups} arguments must be provided.
#' @param group_by A character vector specifying the grouping variable in the Seurat object. This argument is only used if \code{srt} is specified.
#' @param test.use A character vector specifying the test to be used in differential expression analysis. This argument is only used if \code{srt} is specified.
#' @param DE_threshold A character vector specifying the filter condition for differential expression analysis. This argument is only used if \code{srt} is specified.
#' @param geneID A character vector specifying the gene IDs.
#' @param geneID_groups A factor vector specifying the group labels for each gene.
#' @param geneID_exclude A character vector specifying the gene IDs to be excluded from the analysis.
#' @param IDtype A character vector specifying the type of gene IDs in the \code{srt} object or \code{geneID} argument. This argument is used to convert the gene IDs to a different type if \code{IDtype} is different from \code{result_IDtype}.
#' @param result_IDtype A character vector specifying the desired type of gene ID to be used in the output. This argument is used to convert the gene IDs from \code{IDtype} to \code{result_IDtype}.
#' @param species A character vector specifying the species for which the analysis is performed.
#' @param db A character vector specifying the name of the database to be used for enrichment analysis.
#' @param db_update A logical value indicating whether the gene annotation databases should be forcefully updated. If set to FALSE, the function will attempt to load the cached databases instead. Default is FALSE.
#' @param db_version A character vector specifying the version of the database to be used. This argument is ignored if \code{db_update} is \code{TRUE}. Default is "latest".
#' @param db_combine A logical value indicating whether to combine multiple databases into one. If TRUE, all database specified by \code{db} will be combined as one named "Combined".
#' @param convert_species A logical value indicating whether to use a species-converted database when the annotation is missing for the specified species. The default value is TRUE.
#' @param TERM2GENE A data frame specifying the gene-term mapping for a custom database. The first column should contain the term IDs, and the second column should contain the gene IDs.
#' @param TERM2NAME A data frame specifying the term-name mapping for a custom database. The first column should contain the term IDs, and the second column should contain the corresponding term names.
#' @param minGSSize A numeric value specifying the minimum size of a gene set to be considered in the enrichment analysis.
#' @param maxGSSize A numeric value specifying the maximum size of a gene set to be considered in the enrichment analysis.
#' @param unlimited_db A character vector specifying the names of databases that do not have size restrictions.
#' @param GO_simplify A logical value indicating whether to simplify the GO terms. If \code{TRUE}, additional results with simplified GO terms will be returned.
#' @param GO_simplify_cutoff A character vector specifying the filter condition for simplification of GO terms. This argument is only used if \code{GO_simplify} is \code{TRUE}.
#' @param simplify_method A character vector specifying the method to be used for simplification of GO terms. This argument is only used if \code{GO_simplify} is \code{TRUE}.
#' @param simplify_similarityCutoff A numeric value specifying the similarity cutoff for simplification of GO terms. This argument is only used if \code{GO_simplify} is \code{TRUE}.
#' @param BPPARAM A BiocParallelParam object specifying the parallel back-end to be used for parallel computation. Defaults to BiocParallel::bpparam().
#' @param seed The random seed for reproducibility. Defaults to 11.
#'
#' @returns
#' If input is a Seurat object, returns the modified Seurat object with the enrichment result stored in the tools slot.
#'
#' If input is a geneID vector with or without geneID_groups, return the enrichment result directly.
#'
#' Enrichment result is a list with the following component:
#' \itemize{
#'  \item{\code{enrichment}:}{ A data.frame containing all enrichment results.}
#'  \item{\code{results}:}{ A list of \code{enrichResult} objects from the DOSE package.}
#'  \item{\code{geneMap}:}{ A data.frame containing the ID mapping table for input gene IDs.}
#'  \item{\code{input}:}{ A data.frame containing the input gene IDs and gene ID groups.}
#'  \item{\code{DE_threshold}:}{ A specific threshold for differential expression analysis (only returned if input is a Seurat object).}
#' }
#'
#' @seealso \code{\link{PrepareDB}} \code{\link{ListDB}} \code{\link{EnrichmentPlot}} \code{\link{RunGSEA}} \code{\link{GSEAPlot}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunDEtest(
#'   pancreas_sub,
#'   group_by = "CellType"
#' )
#' pancreas_sub <- RunEnrichment(
#'   srt = pancreas_sub,
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
#' pancreas_sub <- RunEnrichment(
#'   srt = pancreas_sub,
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
#'   srt = pancreas_sub,
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
#' # Use a combined database
#' pancreas_sub <- RunEnrichment(
#'   srt = pancreas_sub,
#'   group_by = "CellType",
#'   db = c("KEGG", "WikiPathway", "Reactome", "PFAM", "MP"),
#'   db_combine = TRUE,
#'   species = "Mus_musculus"
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "Combined",
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
#' @importFrom BiocParallel bplapply bpprogressbar<- bpRNGseed<- bpworkers ipcid ipclock ipcunlock
#' @importFrom clusterProfiler enricher simplify
#' @export
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
    Ensembl_version = 103,
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
    BPPARAM = BiocParallel::bpparam(),
    seed = 11) {
  bpprogressbar(BPPARAM) <- TRUE
  bpRNGseed(BPPARAM) <- seed
  time_start <- Sys.time()
  message(paste0("[", time_start, "] ", "Start Enrichment"))
  message("Workers: ", bpworkers(BPPARAM))

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
      stop(
        "Cannot find the DEtest result for the group '",
        group_by,
        "'. You may perform RunDEtest first."
      )
    }
    index <- grep(
      pattern = paste0("AllMarkers_", test.use),
      names(srt@tools[[layer]])
    )[1]
    if (is.na(index)) {
      stop("Cannot find the 'AllMarkers_", test.use, "' in the DEtest result.")
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
    stop("length(geneID_groups)!=length(geneID)")
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
    message("Create 'Combined' database ...")
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
  input <- unnest(input, cols = c(IDtype, result_IDtype))
  input <- input[!is.na(input[[IDtype]]), , drop = FALSE]

  message("Permform enrichment...")
  comb <- expand.grid(
    group = levels(geneID_groups),
    term = db,
    stringsAsFactors = FALSE
  )

  res_list <- bplapply(
    seq_len(nrow(comb)),
    function(i, id) {
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
      enrich_res <- enricher(
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
            warning(
              group,
              "|",
              term,
              " has no term to simplify.",
              immediate. = TRUE
            )
          } else {
            sim_res@result <- sim_res@result[
              with(sim_res@result, eval(rlang::parse_expr(GO_simplify_cutoff))), ,
              drop = FALSE
            ]
            semData <- db_list[[species]][[term]][["semData"]]
            ipclock(id)
            sim_res <- simplify(
              sim_res,
              measure = simplify_method,
              cutoff = simplify_similarityCutoff,
              semData = semData
            )
            ipcunlock(id)
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
      return(enrich_res)
    },
    BPPARAM = BPPARAM,
    id = ipcid()
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

  time_end <- Sys.time()
  message(paste0("[", time_end, "] ", "Enrichment done"))
  message(
    "Elapsed time:",
    format(
      round(difftime(time_end, time_start), 2),
      format = "%Y-%m-%d %H:%M:%S"
    )
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

#' EnrichmentPlot
#'
#' This function generates various types of plots for enrichment (over-representation) analysis.
#'
#' @param srt A Seurat object containing the results of RunDEtest and RunEnrichment.
#' If specified, enrichment results will be extracted from the Seurat object automatically.
#' If not specified, the \code{res} arguments must be provided.
#' @param group_by A character vector specifying the grouping variable in the Seurat object. This argument is only used if \code{srt} is specified.
#' @param test.use A character vector specifying the test to be used in differential expression analysis. This argument is only used if \code{srt} is specified.
#' @param res Enrichment results generated by RunEnrichment function. If provided, 'srt', 'test.use' and 'group_by' are ignored.
#' @param db The database to use for enrichment plot. Default is "GO_BP".
#' @param plot_type The type of plot to generate. Options are: "bar", "dot", "lollipop", "network", "enrichmap", "wordcloud", "comparison". Default is "bar".
#' @param split_by The splitting variable(s) for the plot. Can be "Database", "Groups", or both. Default is c("Database", "Groups") for plots.
#' @param color_by The variable used for coloring. Default is "Database".
#' @param group_use The group(s) to be used for enrichment plot. Default is NULL.
#' @param id_use List of IDs to be used to display specific terms in the enrichment plot. Default value is NULL.
#' @param pvalueCutoff The p-value cutoff. Default is NULL. Only work when \code{padjustCutoff} is NULL.
#' @param padjustCutoff The p-adjusted cutoff. Default is 0.05.
#' @param topTerm The number of top terms to display. Default is 6, or 100 if 'plot_type' is "enrichmap".
#' @param compare_only_sig Whether to compare only significant terms. Default is FALSE.
#' @param topWord The number of top words to display for wordcloud. Default is 100.
#' @param word_type The type of words to display in wordcloud. Options are "term" and "feature". Default is "term".
#' @param word_size The size range for words in wordcloud. Default is c(2, 8).
#' @param words_excluded Words to be excluded from the wordcloud. The default value is NULL, which means that the built-in words (scop::words_excluded) will be used.
#' @param network_layout The layout algorithm to use for network plot. Options are "fr", "kk","random", "circle", "tree", "grid", or other algorithm from 'igraph' package. Default is "fr".
#' @param network_labelsize The label size for network plot. Default is 5.
#' @param network_blendmode The blend mode for network plot. Default is "blend".
#' @param network_layoutadjust Whether to adjust the layout of the network plot to avoid overlapping words. Default is TRUE.
#' @param network_adjscale The scale for adjusting network plot layout. Default is 60.
#' @param network_adjiter The number of iterations for adjusting network plot layout. Default is 100.
#' @param enrichmap_layout The layout algorithm to use for enrichmap plot. Options are "fr", "kk","random", "circle", "tree", "grid", or other algorithm from 'igraph' package. Default is "fr".
#' @param enrichmap_cluster The clustering algorithm to use for enrichmap plot. Options are "walktrap", "fast_greedy", or other algorithm from 'igraph' package. Default is "fast_greedy".
#' @param enrichmap_label  The label type for enrichmap plot. Options are "term" and "feature". Default is "term".
#' @param enrichmap_labelsize The label size for enrichmap plot. Default is 5.
#' @param enrlichmap_nlabel The number of labels to display for each cluster in enrichmap plot. Default is 4.
#' @param enrichmap_show_keyword Whether to show the keyword of terms or features in enrichmap plot. Default is FALSE.
#' @param enrichmap_mark The mark shape for enrichmap plot. Options are "ellipse" and "hull". Default is "ellipse".
#' @param enrichmap_expand The expansion factor for enrichmap plot. Default is c(0.5, 0.5).
#' @param character_width  The maximum width of character of descriptions. Default is 50.
#' @param lineheight The line height for y-axis labels. Default is 0.5.
#' @param palette The color palette to use. Default is "Spectral".
#' @param palcolor Custom colors for palette. Default is NULL.
#' @param aspect.ratio The aspect ratio of the plot. Default is 1.
#' @param legend.position The position of the legend. Default is "right".
#' @param legend.direction The direction of the legend. Default is "vertical".
#' @param theme_use The theme to use for the plot. Default is "theme_scop".
#' @param theme_args The arguments to pass to the theme. Default is an empty list.
#' @param combine Whether to combine multiple plots into a single plot. Default is TRUE.
#' @param nrow The number of rows in the combined plot. Default is NULL, calculated based on the number of plots.
#' @param ncol The number of columns in the combined plot. Default is NULL, calculated based on the number of plots.
#' @param byrow  Whether to fill the combined plot by row. Default is TRUE.
#' @param seed The random seed to use. Default is 11.
#'
#' @seealso \code{\link{RunEnrichment}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunDEtest(
#'   pancreas_sub,
#'   group_by = "CellType"
#' )
#' pancreas_sub <- RunEnrichment(
#'   srt = pancreas_sub,
#'   db = c("GO_BP", "GO_CC"),
#'   group_by = "CellType",
#'   species = "Mus_musculus"
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Ductal",
#'   plot_type = "bar"
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Endocrine",
#'   plot_type = "bar",
#'   character_width = 30,
#'   theme_use = ggplot2::theme_classic,
#'   theme_args = list(base_size = 10)
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   plot_type = "bar",
#'   color_by = "Groups",
#'   ncol = 2
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   plot_type = "bar",
#'   id_use = list(
#'     "Ductal" = c("GO:0002181", "GO:0045787", "GO:0006260", "GO:0050679"),
#'     "Ngn3 low EP" = c("GO:0050678", "GO:0051101", "GO:0072091", "GO:0006631"),
#'     "Ngn3 high EP" = c("GO:0035270", "GO:0030325", "GO:0008637", "GO:0030856"),
#'     "Pre-endocrine" = c("GO:0090276", "GO:0031018", "GO:0030073", "GO:1903532"),
#'     "Endocrine" = c("GO:0009914", "GO:0030073", "GO:0009743", "GO:0042593")
#'   )
#' )
#'
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   topTerm = 3,
#'   plot_type = "comparison"
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   topTerm = 3,
#'   plot_type = "comparison",
#'   compare_only_sig = TRUE
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = c("Ductal", "Endocrine"),
#'   plot_type = "comparison"
#' )
#'
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = c("GO_BP", "GO_CC"),
#'   group_by = "CellType",
#'   group_use = c("Ductal", "Endocrine"),
#'   plot_type = "bar",
#'   split_by = "Groups"
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = c("GO_BP", "GO_CC"),
#'   group_by = "CellType",
#'   group_use = c("Ductal", "Endocrine"),
#'   plot_type = "bar",
#'   split_by = "Database",
#'   color_by = "Groups"
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = c("GO_BP", "GO_CC"),
#'   group_by = "CellType",
#'   group_use = c("Ductal", "Endocrine"),
#'   plot_type = "bar",
#'   split_by = c("Database", "Groups")
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = c("GO_BP", "GO_CC"),
#'   group_by = "CellType",
#'   group_use = c("Ductal", "Endocrine"),
#'   plot_type = "bar",
#'   split_by = c("Groups", "Database")
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = c("GO_BP", "GO_CC"),
#'   group_by = "CellType",
#'   plot_type = "bar",
#'   split_by = "Database",
#'   color_by = "Groups",
#'   palette = "Set1"
#' )
#'
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Ductal",
#'   plot_type = "dot",
#'   palette = "GdRd"
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Ductal",
#'   plot_type = "lollipop",
#'   palette = "GdRd"
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Ductal",
#'   plot_type = "wordcloud"
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Ductal",
#'   plot_type = "wordcloud",
#'   word_type = "feature"
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Ductal",
#'   plot_type = "network"
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Ductal",
#'   plot_type = "network",
#'   id_use = c(
#'     "GO:0050678",
#'     "GO:0035270",
#'     "GO:0090276",
#'     "GO:0030073"
#'   )
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Ductal",
#'   plot_type = "network",
#'   network_layoutadjust = FALSE
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Ductal",
#'   plot_type = "network",
#'   topTerm = 4,
#'   network_blendmode = "average",
#'   theme_use = "theme_blank",
#'   theme_args = list(add_coord = FALSE)
#' ) %>% panel_fix(height = 5)
#'
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Ductal",
#'   plot_type = "enrichmap"
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Ductal",
#'   plot_type = "enrichmap",
#'   enrichmap_expand = c(2, 1)
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Ductal",
#'   plot_type = "enrichmap",
#'   enrichmap_show_keyword = TRUE,
#'   character_width = 10
#' )
#' EnrichmentPlot(pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Ductal",
#'   plot_type = "enrichmap",
#'   topTerm = 200,
#'   enrichmap_mark = "hull",
#'   enrichmap_label = "feature",
#'   enrlichmap_nlabel = 3,
#'   character_width = 10,
#'   theme_use = "theme_blank",
#'   theme_args = list(add_coord = FALSE)
#' ) %>% panel_fix(height = 4)
#'
#' pancreas_sub <- RunEnrichment(
#'   srt = pancreas_sub,
#'   db = c("MP", "DO"),
#'   group_by = "CellType",
#'   convert_species = TRUE,
#'   species = "Mus_musculus"
#' )
#' EnrichmentPlot(
#'   pancreas_sub,
#'   db = c("MP", "DO"),
#'   group_by = "CellType",
#'   group_use = "Endocrine",
#'   ncol = 1
#' )
#'
#' @importFrom ggplot2 ggplot geom_bar geom_text geom_label labs scale_fill_manual scale_y_continuous scale_linewidth facet_grid coord_flip scale_color_gradientn scale_fill_gradientn scale_size guides geom_segment expansion guide_colorbar scale_color_manual guide_none draw_key_point  scale_color_identity scale_fill_identity .pt
#' @importFrom dplyr %>% group_by filter arrange desc across mutate slice_head reframe distinct n .data
#' @importFrom stats formula
#' @importFrom patchwork wrap_plots
#' @importFrom ggforce geom_mark_ellipse geom_mark_hull
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid grobTree convertUnit unit
#' @importFrom igraph as_data_frame graph_from_data_frame V layout_with_dh layout_with_drl layout_with_fr layout_with_gem layout_with_graphopt layout_with_kk layout_with_lgl layout_with_mds layout_in_circle layout_as_tree layout_on_grid cluster_fast_greedy cluster_infomap cluster_leiden cluster_louvain cluster_spinglass cluster_walktrap cluster_fluid_communities
#' @export
EnrichmentPlot <- function(
    srt,
    db = "GO_BP",
    group_by = NULL,
    test.use = "wilcox",
    res = NULL,
    plot_type = c(
      "bar",
      "dot",
      "lollipop",
      "network",
      "enrichmap",
      "wordcloud",
      "comparison"
    ),
    split_by = c("Database", "Groups"),
    color_by = "Database",
    group_use = NULL,
    id_use = NULL,
    pvalueCutoff = NULL,
    padjustCutoff = 0.05,
    topTerm = ifelse(plot_type == "enrichmap", 100, 6),
    compare_only_sig = FALSE,
    topWord = 100,
    word_type = c("term", "feature"),
    word_size = c(2, 8),
    words_excluded = NULL,
    network_layout = "fr",
    network_labelsize = 5,
    network_blendmode = "blend",
    network_layoutadjust = TRUE,
    network_adjscale = 60,
    network_adjiter = 100,
    enrichmap_layout = "fr",
    enrichmap_cluster = "fast_greedy",
    enrichmap_label = c("term", "feature"),
    enrichmap_labelsize = 5,
    enrlichmap_nlabel = 4,
    enrichmap_show_keyword = FALSE,
    enrichmap_mark = c("ellipse", "hull"),
    enrichmap_expand = c(0.5, 0.5),
    character_width = 50,
    lineheight = 0.5,
    palette = "Spectral",
    palcolor = NULL,
    aspect.ratio = 1,
    legend.position = "right",
    legend.direction = "vertical",
    theme_use = "theme_scop",
    theme_args = list(),
    combine = TRUE,
    nrow = NULL,
    ncol = NULL,
    byrow = TRUE,
    seed = 11) {
  set.seed(seed)
  plot_type <- match.arg(plot_type)
  word_type <- match.arg(word_type)
  enrichmap_label <- match.arg(enrichmap_label)
  enrichmap_mark <- match.arg(enrichmap_mark)
  words_excluded <- words_excluded %||% scop::words_excluded

  if (any(!split_by %in% c("Database", "Groups"))) {
    stop("'split_by' must be either 'Database', 'Groups', or both of them")
  }
  if (plot_type %in% c("network", "enrichmap") & length(split_by) == 1) {
    warning(
      "When 'plot_type' is 'network' or 'enrichmap', the 'split_by' parameter does not take effect.",
      immediate. = TRUE
    )
    split_by <- c("Database", "Groups")
  }

  if (is.null(res)) {
    if (is.null(group_by)) {
      stop("'group_by' must be provided.")
    }
    layer <- paste("Enrichment", group_by, test.use, sep = "_")
    if (!layer %in% names(srt@tools)) {
      stop("No enrichment result found. You may perform RunEnrichment first.")
    }
    enrichment <- srt@tools[[layer]][["enrichment"]]
  } else {
    enrichment <- res[["enrichment"]]
  }

  if (is.null(pvalueCutoff) && is.null(padjustCutoff)) {
    stop("One of 'pvalueCutoff' or 'padjustCutoff' must be specified")
  }
  if (!is.factor(enrichment["Groups"])) {
    enrichment[["Groups"]] <- factor(
      enrichment[["Groups"]],
      levels = unique(enrichment[["Groups"]])
    )
  }
  if (length(db[!db %in% enrichment[["Database"]]]) > 0) {
    stop(paste0(
      db[!db %in% enrichment[["Database"]]],
      " is not in the enrichment result."
    ))
  }
  if (!is.factor(enrichment[["Database"]])) {
    enrichment[["Database"]] <- factor(
      enrichment[["Database"]],
      levels = unique(enrichment[["Database"]])
    )
  }
  if (!is.null(group_use)) {
    enrichment <- enrichment[
      enrichment[["Groups"]] %in% group_use, ,
      drop = FALSE
    ]
  }
  if (length(id_use) > 0) {
    topTerm <- Inf
    if (is.list(id_use)) {
      if (is.null(names(id_use))) {
        stop("'id_use' must be named when it is a list.")
      }
      if (!all(names(id_use) %in% enrichment[["Groups"]])) {
        stop(paste0(
          "Names in 'id_use' is invalid: ",
          paste0(
            names(id_use)[!names(id_use) %in% enrichment[["Groups"]]],
            collapse = ","
          )
        ))
      }
      enrichment_list <- list()
      for (i in seq_along(id_use)) {
        enrichment_list[[i]] <- enrichment[
          enrichment[["ID"]] %in%
            id_use[[i]] &
            enrichment[["Groups"]] %in% names(id_use)[i], ,
          drop = FALSE
        ]
      }
      enrichment <- do.call(rbind, enrichment_list)
    } else {
      enrichment <- enrichment[
        enrichment[["ID"]] %in% unlist(id_use), ,
        drop = FALSE
      ]
    }
  }

  metric <- ifelse(is.null(padjustCutoff), "pvalue", "p.adjust")
  metric_value <- ifelse(is.null(padjustCutoff), pvalueCutoff, padjustCutoff)

  pvalueCutoff <- ifelse(is.null(pvalueCutoff), Inf, pvalueCutoff)
  padjustCutoff <- ifelse(is.null(padjustCutoff), Inf, padjustCutoff)

  if (any(db %in% c("GO_sim", "GO_BP_sim", "GO_CC_sim", "GO_MF_sim"))) {
    enrichment_sim <- enrichment[
      enrichment[["Database"]] %in% gsub("_sim", "", db), ,
      drop = FALSE
    ]
  }
  enrichment <- enrichment[enrichment[["Database"]] %in% db, , drop = FALSE]

  enrichment_sig <- enrichment[
    enrichment[[metric]] < metric_value |
      enrichment[["ID"]] %in% unlist(id_use), ,
    drop = FALSE
  ]
  enrichment_sig <- enrichment_sig[
    order(enrichment_sig[[metric]]), ,
    drop = FALSE
  ]
  if (nrow(enrichment_sig) == 0) {
    stop(
      "No term enriched using the threshold: ",
      paste0("pvalueCutoff = ", pvalueCutoff),
      "; ",
      paste0("padjustCutoff = ", padjustCutoff)
    )
  }
  df_list <- split(
    enrichment_sig,
    formula(paste0("~", split_by, collapse = "+"))
  )
  df_list <- df_list[lapply(df_list, nrow) > 0]

  facet <- switch(paste0(split_by, collapse = "~"),
    "Groups" = formula(paste0("Database ~ Groups")),
    "Database" = formula(paste0("Groups ~ Database")),
    formula(paste0(split_by, collapse = "~"))
  )

  if (plot_type == "comparison") {
    # comparison -------------------------------------------------------------------------------------------------
    ids <- NULL
    for (i in seq_along(df_list)) {
      df <- df_list[[i]]
      df_groups <- split(df, list(df$Database, df$Groups))
      df_groups <- lapply(df_groups, function(group) {
        filtered_group <- group[
          head(seq_len(nrow(group)), topTerm), ,
          drop = FALSE
        ]
        return(filtered_group)
      })
      df <- do.call(rbind, df_groups)
      ids <- unique(c(ids, df[, "ID"]))
    }
    if (any(db %in% c("GO_sim", "GO_BP_sim", "GO_CC_sim", "GO_MF_sim"))) {
      enrichment_sub <- subset(enrichment_sim, ID %in% ids)
      enrichment_sub[["Database"]][
        enrichment_sub[["Database"]] %in% c("GO", "GO_BP", "GO_CC", "GO_MF")
      ] <- paste0(
        enrichment_sub[["Database"]][
          enrichment_sub[["Database"]] %in% c("GO", "GO_BP", "GO_CC", "GO_MF")
        ],
        "_sim"
      )
    } else {
      enrichment_sub <- subset(enrichment, ID %in% ids)
    }
    enrichment_sub[["Database"]] <- factor(
      enrichment_sub[["Database"]],
      levels = db
    )
    enrichment_sub[["GeneRatio"]] <- sapply(
      enrichment_sub[["GeneRatio"]],
      function(x) {
        sp <- strsplit(x, "/")[[1]]
        GeneRatio <- as.numeric(sp[1]) / as.numeric(sp[2])
      }
    )
    enrichment_sub[["BgRatio"]] <- sapply(
      enrichment_sub[["BgRatio"]],
      function(x) {
        sp <- strsplit(x, "/")[[1]]
        BgRatio <- as.numeric(sp[1]) / as.numeric(sp[2])
        return(BgRatio)
      }
    )
    enrichment_sub[["EnrichmentScore"]] <- enrichment_sub[["GeneRatio"]] /
      enrichment_sub[["BgRatio"]]
    enrichment_sub[["Description"]] <- capitalize(enrichment_sub[[
      "Description"
    ]])
    enrichment_sub[["Description"]] <- str_wrap(
      enrichment_sub[["Description"]],
      width = character_width
    )
    terms <- stats::setNames(enrichment_sub[["Description"]], enrichment_sub[["ID"]])
    enrichment_sub[["Description"]] <- factor(
      enrichment_sub[["Description"]],
      levels = unique(rev(terms[ids]))
    )
    if (isTRUE(compare_only_sig)) {
      enrichment_sub <- enrichment_sub[
        enrichment_sub[[metric]] < metric_value, ,
        drop = FALSE
      ]
    }
    p <- ggplot(enrichment_sub, aes(x = Groups, y = Description)) +
      geom_point(
        aes(size = GeneRatio, fill = .data[[metric]], color = ""),
        shape = 21
      ) +
      scale_size_area(name = "GeneRatio", max_size = 6, n.breaks = 4) +
      guides(
        size = guide_legend(
          override.aes = list(fill = "grey30", shape = 21),
          order = 1
        )
      ) +
      scale_fill_gradientn(
        name = paste0(metric),
        limits = c(0, min(metric_value, 1)),
        n.breaks = 3,
        colors = palette_scop(
          palette = palette,
          palcolor = palcolor,
          reverse = TRUE
        ),
        na.value = "grey80",
        guide = guide_colorbar(
          frame.colour = "black",
          ticks.colour = "black",
          title.hjust = 0,
          order = 2
        )
      ) +
      scale_color_manual(values = NA, na.value = "black") +
      guides(
        colour = if (isTRUE(compare_only_sig)) {
          guide_none()
        } else {
          guide_legend(
            "Non-sig",
            override.aes = list(colour = "black", fill = "grey80", size = 3)
          )
        }
      ) +
      facet_grid(Database ~ ., scales = "free") +
      do.call(theme_use, theme_args) +
      theme(
        aspect.ratio = aspect.ratio,
        legend.position = legend.position,
        legend.direction = legend.direction,
        panel.grid.major = element_line(colour = "grey80", linetype = 2),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(
          lineheight = lineheight,
          hjust = 1,
          face = ifelse(
            grepl("\n", levels(enrichment_sub[["Description"]])),
            "italic",
            "plain"
          )
        )
      )
    plist <- list(p)
  } else if (plot_type == "bar") {
    # bar -------------------------------------------------------------------------------------------------
    plist <- suppressWarnings(lapply(df_list, function(df) {
      df_groups <- split(df, list(df$Database, df$Groups))
      df_groups <- lapply(df_groups, function(group) {
        filtered_group <- group[
          head(seq_len(nrow(group)), topTerm), ,
          drop = FALSE
        ]
        return(filtered_group)
      })
      df <- do.call(rbind, df_groups)

      df[["metric"]] <- -log10(df[[metric]])
      df[["Description"]] <- capitalize(df[["Description"]])
      df[["Description"]] <- str_wrap(
        df[["Description"]],
        width = character_width
      )
      df[["Description"]] <- factor(
        df[["Description"]],
        levels = unique(rev(df[["Description"]]))
      )

      p <- ggplot(
        df,
        aes(
          x = .data[["Description"]],
          y = .data[["metric"]],
          fill = .data[[color_by]],
          label = .data[["Count"]]
        )
      ) +
        geom_bar(width = 0.9, stat = "identity", color = "black") +
        geom_text(
          hjust = -0.5,
          size = 3.5,
          color = "white",
          fontface = "bold"
        ) +
        geom_text(hjust = -0.5, size = 3.5) +
        labs(x = "", y = paste0("-log10(", metric, ")")) +
        scale_fill_manual(
          values = palette_scop(
            levels(df[[color_by]]),
            palette = palette,
            palcolor = palcolor
          ),
          na.value = "grey80",
          guide = "none"
        ) +
        scale_y_continuous(
          limits = c(0, 1.3 * max(df[["metric"]], na.rm = TRUE)),
          expand = expansion(0, 0)
        ) +
        facet_grid(facet, scales = "free") +
        coord_flip() +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction,
          panel.grid.major = element_line(colour = "grey80", linetype = 2),
          axis.text.y = element_text(
            lineheight = lineheight,
            hjust = 1,
            face = ifelse(
              grepl("\n", levels(df[["Description"]])),
              "italic",
              "plain"
            )
          )
        )
      return(p)
    }))
  } else if (plot_type == "dot") {
    # dot -------------------------------------------------------------------------------------------------
    plist <- suppressWarnings(lapply(df_list, function(df) {
      df_groups <- split(df, list(df$Database, df$Groups))
      df_groups <- lapply(df_groups, function(group) {
        filtered_group <- group[
          head(seq_len(nrow(group)), topTerm), ,
          drop = FALSE
        ]
        return(filtered_group)
      })
      df <- do.call(rbind, df_groups)

      df[["GeneRatio"]] <- sapply(df[["GeneRatio"]], function(x) {
        sp <- strsplit(x, "/")[[1]]
        GeneRatio <- as.numeric(sp[1]) / as.numeric(sp[2])
      })
      df <- df[order(df[["GeneRatio"]], decreasing = TRUE), ]
      df[["metric"]] <- -log10(df[[metric]])
      df[["Description"]] <- capitalize(df[["Description"]])
      df[["Description"]] <- str_wrap(
        df[["Description"]],
        width = character_width
      )
      df[["Description"]] <- factor(
        df[["Description"]],
        levels = unique(rev(df[["Description"]]))
      )

      p <- ggplot(
        df,
        aes(
          x = .data[["Description"]],
          y = .data[["GeneRatio"]]
        )
      ) +
        geom_point(
          aes(fill = .data[["metric"]], size = .data[["Count"]]),
          color = "black",
          shape = 21
        ) +
        labs(x = "", y = "GeneRatio") +
        scale_size(
          name = "Count",
          range = c(3, 6),
          scales::breaks_extended(n = 4)
        ) +
        guides(
          size = guide_legend(
            override.aes = list(fill = "grey30", shape = 21),
            order = 1
          )
        ) +
        scale_fill_gradientn(
          name = paste0("-log10(", metric, ")"),
          n.breaks = 3,
          colors = palette_scop(palette = palette, palcolor = palcolor),
          na.value = "grey80",
          guide = guide_colorbar(
            frame.colour = "black",
            ticks.colour = "black",
            title.hjust = 0
          )
        ) +
        scale_y_continuous(
          limits = c(0, 1.3 * max(df[["GeneRatio"]], na.rm = TRUE)),
          expand = expansion(0, 0)
        ) +
        facet_grid(facet, scales = "free") +
        coord_flip() +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction,
          panel.grid.major = element_line(colour = "grey80", linetype = 2),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(
            lineheight = lineheight,
            hjust = 1,
            face = ifelse(
              grepl("\n", levels(df[["Description"]])),
              "italic",
              "plain"
            )
          )
        )
      return(p)
    }))
  } else if (plot_type == "lollipop") {
    # lollipop -------------------------------------------------------------------------------------------------
    plist <- suppressWarnings(lapply(df_list, function(df) {
      df_groups <- split(df, list(df$Database, df$Groups))
      df_groups <- lapply(df_groups, function(group) {
        filtered_group <- group[
          head(seq_len(nrow(group)), topTerm), ,
          drop = FALSE
        ]
        return(filtered_group)
      })
      df <- do.call(rbind, df_groups)

      df[["GeneRatio"]] <- sapply(df[["GeneRatio"]], function(x) {
        sp <- strsplit(x, "/")[[1]]
        GeneRatio <- as.numeric(sp[1]) / as.numeric(sp[2])
      })
      df[["BgRatio"]] <- sapply(df[["BgRatio"]], function(x) {
        sp <- strsplit(x, "/")[[1]]
        BgRatio <- as.numeric(sp[1]) / as.numeric(sp[2])
        return(BgRatio)
      })
      df[["FoldEnrichment"]] <- df[["GeneRatio"]] / df[["BgRatio"]]
      df[["metric"]] <- -log10(df[[metric]])
      df[["Description"]] <- capitalize(df[["Description"]])
      df[["Description"]] <- str_wrap(
        df[["Description"]],
        width = character_width
      )
      df[["Description"]] <- factor(
        df[["Description"]],
        levels = unique(df[order(df[["FoldEnrichment"]]), "Description"])
      )

      p <- ggplot(
        df,
        aes(
          x = .data[["Description"]],
          y = .data[["FoldEnrichment"]],
          fill = .data[["metric"]]
        )
      ) +
        geom_blank() +
        geom_segment(
          aes(
            y = 0,
            xend = .data[["Description"]],
            yend = .data[["FoldEnrichment"]]
          ),
          color = "black",
          linewidth = 2
        ) +
        geom_segment(
          aes(
            y = 0,
            xend = .data[["Description"]],
            yend = .data[["FoldEnrichment"]],
            color = .data[["metric"]]
          ),
          linewidth = 1
        ) +
        geom_point(
          aes(size = .data[["GeneRatio"]]),
          shape = 21,
          color = "black"
        ) +
        scale_size(
          name = "GeneRatio",
          range = c(3, 6),
          scales::breaks_extended(n = 4)
        ) +
        guides(
          size = guide_legend(
            override.aes = list(fill = "grey30", shape = 21),
            order = 1
          )
        ) +
        scale_y_continuous(
          limits = c(0, 1.2 * max(df[["FoldEnrichment"]], na.rm = TRUE)),
          expand = expansion(0, 0)
        ) +
        labs(x = "", y = "Fold Enrichment") +
        scale_fill_gradientn(
          name = paste0("-log10(", metric, ")"),
          n.breaks = 3,
          colors = palette_scop(palette = palette, palcolor = palcolor),
          na.value = "grey80",
          guide = guide_colorbar(
            frame.colour = "black",
            ticks.colour = "black",
            title.hjust = 0
          ),
          aesthetics = c("color", "fill")
        ) +
        facet_grid(facet, scales = "free") +
        coord_flip() +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction,
          panel.grid.major = element_line(colour = "grey80", linetype = 2),
          axis.text.y = element_text(
            lineheight = lineheight,
            hjust = 1,
            face = ifelse(
              grepl("\n", levels(df[["Description"]])),
              "italic",
              "plain"
            )
          )
        )
      return(p)
    }))
  } else if (plot_type == "network") {
    # network -------------------------------------------------------------------------------------------------
    plist <- suppressWarnings(lapply(df_list, function(df) {
      df_groups <- split(df, list(df$Database, df$Groups))
      df_groups <- lapply(df_groups, function(group) {
        filtered_group <- group[
          head(seq_len(nrow(group)), topTerm), ,
          drop = FALSE
        ]
        return(filtered_group)
      })
      df <- do.call(rbind, df_groups)

      df[["metric"]] <- -log10(df[[metric]])
      df[["Description"]] <- capitalize(df[["Description"]])
      df[["Description"]] <- str_wrap(
        df[["Description"]],
        width = character_width
      )
      df[["Description"]] <- factor(
        df[["Description"]],
        levels = unique(df[["Description"]])
      )
      df$geneID <- strsplit(df$geneID, "/")
      df_unnest <- unnest(df, cols = "geneID")

      nodes <- rbind(
        data.frame(
          "ID" = df[["Description"]],
          class = "term",
          metric = df[["metric"]]
        ),
        data.frame("ID" = unique(df_unnest$geneID), class = "gene", metric = 0)
      )
      nodes$Database <- df$Database[1]
      nodes$Groups <- df$Groups[1]
      edges <- as.data.frame(df_unnest[, c("Description", "geneID")])
      colnames(edges) <- c("from", "to")
      edges[["weight"]] <- 1
      graph <- graph_from_data_frame(
        d = edges,
        vertices = nodes,
        directed = FALSE
      )
      if (network_layout %in% c("circle", "tree", "grid")) {
        layout <- switch(network_layout,
          "circle" = layout_in_circle(graph),
          "tree" = layout_as_tree(graph),
          "grid" = layout_on_grid(graph)
        )
      } else {
        layout <- do.call(paste0("layout_with_", network_layout), list(graph))
      }
      df_graph <- as_data_frame(graph, what = "both")

      df_nodes <- df_graph$vertices
      if (isTRUE(network_layoutadjust)) {
        width <- nchar(df_nodes$name)
        width[df_nodes$class == "term"] <- 8
        layout <- adjustlayout(
          graph = graph,
          layout = layout,
          width = width,
          height = 2,
          scale = network_adjscale,
          iter = network_adjiter
        )
      }
      df_nodes[["dim1"]] <- layout[, 1]
      df_nodes[["dim2"]] <- layout[, 2]

      df_edges <- df_graph$edges
      df_edges[["from_dim1"]] <- df_nodes[df_edges[["from"]], "dim1"]
      df_edges[["from_dim2"]] <- df_nodes[df_edges[["from"]], "dim2"]
      df_edges[["to_dim1"]] <- df_nodes[df_edges[["to"]], "dim1"]
      df_edges[["to_dim2"]] <- df_nodes[df_edges[["to"]], "dim2"]

      colors <- palette_scop(
        levels(df[["Description"]]),
        palette = palette,
        palcolor = palcolor
      )
      df_edges[["color"]] <- colors[df_edges$from]
      node_colors <- aggregate(
        df_unnest$Description,
        by = list(df_unnest$geneID),
        FUN = function(x) {
          blendcolors(colors = colors[x], mode = network_blendmode)
        }
      )
      colors <- c(colors, stats::setNames(node_colors[, 2], node_colors[, 1]))
      label_colors <- ifelse(
        colSums(col2rgb(colors)) > 255 * 2,
        "black",
        "white"
      )
      df_nodes[["color"]] <- colors[df_nodes$name]
      df_nodes[["label_color"]] <- label_colors[df_nodes$name]
      df_nodes[["label"]] <- NA
      df_nodes[levels(df[["Description"]]), "label"] <- seq_len(nlevels(df[[
        "Description"
      ]]))

      draw_key_cust <- function(data, params, size) {
        data_text <- data
        data_text$label <- which(
          levels(df[["Description"]]) %in%
            names(colors)[colors == data_text$fill]
        )
        data_text$colour <- "black"
        data_text$alpha <- 1
        data_text$size <- 11 / .pt
        grobTree(
          draw_key_point(data, list(color = "white", shape = 21)),
          ggrepel:::shadowtextGrob(
            label = data_text$label,
            bg.colour = "black",
            bg.r = 0.1,
            gp = gpar(col = "white", fontface = "bold")
          )
        )
      }

      p <- ggplot() +
        geom_segment(
          data = df_edges,
          aes(
            x = from_dim1,
            y = from_dim2,
            xend = to_dim1,
            yend = to_dim2,
            color = color
          ),
          alpha = 1,
          lineend = "round",
          show.legend = FALSE
        ) +
        geom_label(
          data = df_nodes[df_nodes$class == "gene", ],
          aes(
            x = dim1,
            y = dim2,
            label = name,
            fill = color,
            color = label_color
          ),
          size = 3,
          show.legend = FALSE
        ) +
        geom_point(
          data = df_nodes[df_nodes$class == "term", ],
          aes(x = dim1, y = dim2),
          size = 8,
          color = "black",
          fill = "black",
          stroke = 1,
          shape = 21,
          show.legend = FALSE
        ) +
        geom_point(
          data = df_nodes[df_nodes$class == "term", ],
          aes(x = dim1, y = dim2, fill = color),
          size = 7,
          color = "white",
          stroke = 1,
          shape = 21,
          key_glyph = draw_key_cust
        ) +
        geom_text_repel(
          data = df_nodes[df_nodes$class == "term", ],
          aes(x = dim1, y = dim2, label = label),
          fontface = "bold",
          min.segment.length = 0,
          segment.color = "black",
          point.size = NA,
          max.overlaps = 100,
          force = 0,
          color = "white",
          bg.color = "black",
          bg.r = 0.1,
          size = network_labelsize
        ) +
        scale_color_identity(guide = "none") +
        scale_fill_identity(
          name = "Term:",
          guide = "legend",
          labels = levels(df[["Description"]]),
          breaks = colors[levels(df[["Description"]])]
        ) +
        guides(
          color = guide_legend(override.aes = list(color = "transparent"))
        ) +
        labs(x = "", y = "") +
        facet_grid(facet, scales = "free") +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction
        )
      return(p)
    }))
  } else if (plot_type == "enrichmap") {
    # enrichmap -------------------------------------------------------------------------------------------------
    plist <- suppressWarnings(lapply(df_list, function(df) {
      df_groups <- split(df, list(df$Database, df$Groups))
      df_groups <- lapply(df_groups, function(group) {
        filtered_group <- group[
          head(seq_len(nrow(group)), topTerm), ,
          drop = FALSE
        ]
        return(filtered_group)
      })
      df <- do.call(rbind, df_groups)

      df[["metric"]] <- -log10(df[[metric]])
      df[["Description"]] <- capitalize(df[["Description"]])
      df[["Description"]] <- factor(
        df[["Description"]],
        levels = unique(df[["Description"]])
      )
      df$geneID <- strsplit(df$geneID, "/")
      rownames(df) <- df[["ID"]]

      nodes <- df
      edges <- as.data.frame(t(combn(nodes$ID, 2)))
      colnames(edges) <- c("from", "to")
      edges[["weight"]] <- mapply(
        function(x, y) length(intersect(df[[x, "geneID"]], df[[y, "geneID"]])),
        edges$from,
        edges$to
      )
      edges <- edges[edges[["weight"]] > 0, , drop = FALSE]
      graph <- graph_from_data_frame(
        d = edges,
        vertices = nodes,
        directed = FALSE
      )
      if (enrichmap_layout %in% c("circle", "tree", "grid")) {
        layout <- switch(enrichmap_layout,
          "circle" = layout_in_circle(graph),
          "tree" = layout_as_tree(graph),
          "grid" = layout_on_grid(graph)
        )
      } else {
        layout <- do.call(paste0("layout_with_", enrichmap_layout), list(graph))
      }
      clusters <- do.call(paste0("cluster_", enrichmap_cluster), list(graph))
      df_graph <- as_data_frame(graph, what = "both")

      df_nodes <- df_graph$vertices
      df_nodes[["dim1"]] <- layout[, 1]
      df_nodes[["dim2"]] <- layout[, 2]
      df_nodes[["clusters"]] <- factor(
        paste0("C", clusters$membership),
        paste0("C", unique(sort(clusters$membership)))
      )

      if (isTRUE(enrichmap_show_keyword)) {
        df_keyword1 <- df_nodes %>%
          mutate(
            keyword = strsplit(
              tolower(as.character(.data[["Description"]])),
              "\\s|\\n",
              perl = TRUE
            )
          ) %>%
          unnest(cols = "keyword") %>%
          group_by(.data[["keyword"]], Database, Groups, clusters) %>%
          reframe(
            keyword = capitalize(.data[["keyword"]]),
            score = sum(-(log10(.data[[metric]]))),
            count = n(),
            Database = .data[["Database"]],
            Groups = .data[["Groups"]],
            .groups = "keep"
          ) %>%
          filter(!grepl(pattern = "\\[.*\\]", x = .data[["keyword"]])) %>%
          filter(nchar(.data[["keyword"]]) >= 1) %>%
          filter(!tolower(.data[["keyword"]]) %in% tolower(words_excluded)) %>%
          distinct() %>%
          group_by(Database, Groups, clusters) %>%
          arrange(desc(score)) %>%
          slice_head(n = enrlichmap_nlabel) %>%
          reframe(keyword = paste0(.data[["keyword"]], collapse = " ")) %>%
          as.data.frame()
        rownames(df_keyword1) <- as.character(df_keyword1[["clusters"]])
        df_keyword1[["keyword"]] <- str_wrap(
          df_keyword1[["keyword"]],
          width = character_width
        )
        df_keyword1[["label"]] <- paste0(
          df_keyword1[["clusters"]],
          ":\n",
          df_keyword1[["keyword"]]
        )
      } else {
        if (enrichmap_label == "term") {
          df_nodes[["Description"]] <- str_wrap(
            df_nodes[["Description"]],
            width = character_width
          )
        }
        df_keyword1 <- df_nodes %>%
          group_by(Database, Groups, clusters) %>%
          arrange(desc(metric)) %>%
          reframe(keyword = Description) %>%
          distinct() %>%
          group_by(Database, Groups, clusters) %>%
          slice_head(n = enrlichmap_nlabel) %>%
          reframe(keyword = paste0(.data[["keyword"]], collapse = "\n")) %>%
          as.data.frame()
        rownames(df_keyword1) <- as.character(df_keyword1[["clusters"]])
        df_keyword1[["label"]] <- paste0(
          df_keyword1[["clusters"]],
          ":\n",
          df_keyword1[["keyword"]]
        )
      }

      df_keyword2 <- df_nodes %>%
        mutate(keyword = .data[["geneID"]]) %>%
        unnest(cols = "keyword") %>%
        group_by(.data[["keyword"]], Database, Groups, clusters) %>%
        reframe(
          keyword = .data[["keyword"]],
          score = sum(-(log10(.data[[metric]]))),
          count = n(),
          Database = .data[["Database"]],
          Groups = .data[["Groups"]],
          .groups = "keep"
        ) %>%
        distinct() %>%
        group_by(Database, Groups, clusters) %>%
        arrange(desc(score)) %>%
        slice_head(n = enrlichmap_nlabel) %>%
        reframe(keyword = paste0(.data[["keyword"]], collapse = " ")) %>%
        as.data.frame()
      rownames(df_keyword2) <- as.character(df_keyword2[["clusters"]])
      df_keyword2[["keyword"]] <- str_wrap(
        df_keyword2[["keyword"]],
        width = character_width
      )
      df_keyword2[["label"]] <- paste0(
        df_keyword2[["clusters"]],
        ":\n",
        df_keyword2[["keyword"]]
      )

      df_nodes[["keyword1"]] <- df_keyword1[
        as.character(df_nodes$clusters),
        "keyword"
      ]
      df_nodes[["keyword2"]] <- df_keyword2[
        as.character(df_nodes$clusters),
        "keyword"
      ]

      df_edges <- df_graph$edges
      df_edges[["from_dim1"]] <- df_nodes[df_edges[["from"]], "dim1"]
      df_edges[["from_dim2"]] <- df_nodes[df_edges[["from"]], "dim2"]
      df_edges[["to_dim1"]] <- df_nodes[df_edges[["to"]], "dim1"]
      df_edges[["to_dim2"]] <- df_nodes[df_edges[["to"]], "dim2"]

      if (enrichmap_mark == "hull") {
        check_r("concaveman")
      }
      mark_layer <- do.call(
        switch(enrichmap_mark,
          "ellipse" = "geom_mark_ellipse",
          "hull" = "geom_mark_hull"
        ),
        list(
          data = df_nodes,
          aes(
            x = dim1,
            y = dim2,
            color = clusters,
            fill = clusters,
            label = clusters,
            description = if (enrichmap_label == "term") keyword1 else keyword2
          ),
          expand = unit(3, "mm"),
          alpha = 0.1,
          label.margin = margin(1, 1, 1, 1, "mm"),
          label.fontsize = enrichmap_labelsize * 2,
          label.fill = "grey95",
          label.minwidth = unit(character_width, "in"),
          label.buffer = unit(0, "mm"),
          con.size = 1,
          con.cap = 0
        )
      )

      p <- ggplot() +
        mark_layer +
        geom_segment(
          data = df_edges,
          aes(
            x = from_dim1,
            y = from_dim2,
            xend = to_dim1,
            yend = to_dim2,
            linewidth = weight
          ),
          alpha = 0.1,
          lineend = "round"
        ) +
        geom_point(
          data = df_nodes,
          aes(x = dim1, y = dim2, size = Count, fill = clusters),
          color = "black",
          shape = 21
        ) +
        labs(x = "", y = "") +
        scale_size(
          name = "Count",
          range = c(2, 6),
          scales::breaks_extended(n = 4)
        ) +
        guides(
          size = guide_legend(
            override.aes = list(fill = "grey30", shape = 21),
            order = 1
          )
        ) +
        scale_linewidth(
          name = "Intersection",
          range = c(0.3, 3),
          scales::breaks_extended(n = 4)
        ) +
        guides(
          linewidth = guide_legend(
            override.aes = list(alpha = 1, color = "grey"),
            order = 2
          )
        ) +
        scale_fill_manual(
          name = switch(enrichmap_label,
            "term" = "Feature:",
            "feature" = "Term:"
          ),
          values = palette_scop(
            levels(df_nodes[["clusters"]]),
            palette = palette,
            palcolor = palcolor
          ),
          labels = if (enrichmap_label == "term") {
            df_keyword2[levels(df_nodes[["clusters"]]), "label"]
          } else {
            df_keyword1[levels(df_nodes[["clusters"]]), "label"]
          },
          na.value = "grey80",
          aesthetics = c("colour", "fill")
        ) +
        guides(
          fill = guide_legend(
            override.aes = list(alpha = 1, color = "black", shape = NA),
            byrow = TRUE,
            order = 3
          )
        ) +
        guides(color = guide_none()) +
        scale_x_continuous(
          expand = expansion(c(enrichmap_expand[1], enrichmap_expand[1]), 0)
        ) +
        scale_y_continuous(
          expand = expansion(c(enrichmap_expand[2], enrichmap_expand[2]), 0)
        ) +
        facet_grid(facet, scales = "free") +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction
        )
      return(p)
    }))
  } else if (plot_type == "wordcloud") {
    # wordcloud -------------------------------------------------------------------------------------------------
    check_r("ggwordcloud")
    check_r("jokergoo/simplifyEnrichment")
    plist <- lapply(df_list, function(df) {
      if (word_type == "term") {
        df_groups <- split(df, list(df$Database, df$Groups))
        df_groups <- df_groups[sapply(df_groups, nrow) > 0]
        for (i in seq_along(df_groups)) {
          df_sub <- df_groups[[i]]
          if (all(df_sub$Database %in% c("GO", "GO_BP", "GO_CC", "GO_MF"))) {
            df0 <- simplifyEnrichment::keyword_enrichment_from_GO(df_sub[[
              "ID"
            ]])
            if (nrow(df0 > 0)) {
              df_sub <- df0 %>%
                reframe(
                  keyword = .data[["keyword"]],
                  score = -(log10(.data[["padj"]])),
                  count = .data[["n_term"]],
                  Database = df_sub[["Database"]][1],
                  Groups = df_sub[["Groups"]][1]
                ) %>%
                filter(!grepl(pattern = "\\[.*\\]", x = .data[["keyword"]])) %>%
                filter(nchar(.data[["keyword"]]) >= 1) %>%
                filter(
                  !tolower(.data[["keyword"]]) %in% tolower(words_excluded)
                ) %>%
                distinct() %>%
                mutate(
                  angle = 90 *
                    sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))
                ) %>%
                as.data.frame()
              df_sub <- df_sub[
                head(order(df_sub[["score"]], decreasing = TRUE), topWord), ,
                drop = FALSE
              ]
            } else {
              df_sub <- NULL
            }
          } else {
            df_sub <- df_sub %>%
              mutate(
                keyword = strsplit(
                  tolower(as.character(.data[["Description"]])),
                  " "
                )
              ) %>%
              unnest(cols = "keyword") %>%
              group_by(.data[["keyword"]], Database, Groups) %>%
              reframe(
                keyword = .data[["keyword"]],
                score = sum(-(log10(.data[[metric]]))),
                count = n(),
                Database = .data[["Database"]],
                Groups = .data[["Groups"]],
                .groups = "keep"
              ) %>%
              filter(!grepl(pattern = "\\[.*\\]", x = .data[["keyword"]])) %>%
              filter(nchar(.data[["keyword"]]) >= 1) %>%
              filter(
                !tolower(.data[["keyword"]]) %in% tolower(words_excluded)
              ) %>%
              distinct() %>%
              mutate(
                angle = 90 *
                  sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))
              ) %>%
              as.data.frame()
            df_sub <- df_sub[
              head(order(df_sub[["score"]], decreasing = TRUE), topWord), ,
              drop = FALSE
            ]
          }
          df_groups[[i]] <- df_sub
        }
        df <- do.call(rbind, df_groups)
      } else {
        df <- df %>%
          mutate(keyword = strsplit(as.character(.data[["geneID"]]), "/")) %>%
          unnest(cols = "keyword") %>%
          group_by(.data[["keyword"]], Database, Groups) %>%
          reframe(
            keyword = .data[["keyword"]],
            score = sum(-(log10(.data[[metric]]))),
            count = n(),
            Database = .data[["Database"]],
            Groups = .data[["Groups"]],
            .groups = "keep"
          ) %>%
          distinct() %>%
          mutate(
            angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))
          ) %>%
          as.data.frame()
        df <- df[
          head(order(df[["score"]], decreasing = TRUE), topWord), ,
          drop = FALSE
        ]
      }
      colors <- palette_scop(
        df[["score"]],
        type = "continuous",
        palette = palette,
        palcolor = palcolor,
        matched = FALSE
      )
      colors_value <- seq(
        min(df[["score"]], na.rm = TRUE),
        quantile(df[["score"]], 0.99, na.rm = TRUE) + 0.001,
        length.out = 100
      )
      p <- ggplot(
        df,
        aes(
          label = .data[["keyword"]],
          size = .data[["count"]],
          color = .data[["score"]],
          angle = .data[["angle"]]
        )
      ) +
        ggwordcloud::geom_text_wordcloud(
          rm_outside = TRUE,
          eccentricity = 1,
          shape = "square",
          show.legend = TRUE,
          grid_margin = 3
        ) +
        scale_color_gradientn(
          name = "Score:",
          colours = colors,
          values = rescale(colors_value),
          guide = guide_colorbar(
            frame.colour = "black",
            ticks.colour = "black",
            title.hjust = 0
          )
        ) +
        scale_size(
          name = "Count",
          range = word_size,
          breaks = ceiling(seq(
            min(df[["count"]], na.rm = TRUE),
            max(df[["count"]], na.rm = TRUE),
            length.out = 3
          ))
        ) +
        guides(
          size = guide_legend(
            override.aes = list(colour = "black", label = "G"),
            order = 1
          )
        ) +
        facet_grid(facet, scales = "free") +
        coord_flip() +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction
        )
      return(p)
    })
  }

  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- wrap_plots(
        plotlist = plist,
        nrow = nrow,
        ncol = ncol,
        byrow = byrow
      )
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
}
