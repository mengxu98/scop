#' Perform the enrichment analysis (GSEA) on the genes
#'
#' @inheritParams RunEnrichment
#' @param geneScore A numeric vector that specifies the gene scores, for example, the log2(fold change) values of gene expression.
#' @param scoreType This parameter defines the GSEA score type. Possible options are "std", "pos", "neg". By default ("std") the enrichment score is computed as in the original GSEA. The "pos" and "neg" score types are intended to be used for one-tailed tests (i.e. when one is interested only in positive ("pos") or negateive ("neg") enrichment).
#' @returns
#' If input is a Seurat object, returns the modified Seurat object with the enrichment result stored in the tools slot.
#'
#' If input is a geneID vector with or without geneID_groups, return the enrichment result directly.
#'
#' Enrichment result is a list with the following component:
#' \itemize{
#'  \item{\code{enrichment}:}{ A data.frame containing all enrichment results.}
#'  \item{\code{results}:}{ A list of \code{gseaResult} objects from the DOSE package.}
#'  \item{\code{geneMap}:}{ A data.frame containing the ID mapping table for input gene IDs.}
#'  \item{\code{input}:}{ A data.frame containing the input gene IDs and gene ID groups.}
#'  \item{\code{DE_threshold}:}{ A specific threshold for differential expression analysis (only returned if input is a Seurat object).}
#' }
#'
#' @seealso \code{\link{PrepareDB}} \code{\link{ListDB}} \code{\link{GSEAPlot}} \code{\link{RunEnrichment}} \code{\link{EnrichmentPlot}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunDEtest(
#'   pancreas_sub,
#'   group_by = "CellType",
#'   only.pos = FALSE,
#'   fc.threshold = 1
#' )
#' pancreas_sub <- RunGSEA(
#'   pancreas_sub,
#'   group_by = "CellType",
#'   DE_threshold = "p_val_adj < 0.05",
#'   scoreType = "std",
#'   db = "GO_BP",
#'   species = "Mus_musculus"
#' )
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   plot_type = "comparison"
#' )
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Ductal",
#'   id_use = "GO:0006412"
#' )
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Endocrine",
#'   id_use = c("GO:0046903", "GO:0015031", "GO:0007600")
#' )
#'
#' # Remove redundant GO terms
#' pancreas_sub <- RunGSEA(
#'   srt = pancreas_sub,
#'   group_by = "CellType",
#'   db = "GO_BP",
#'   GO_simplify = TRUE,
#'   species = "Mus_musculus"
#' )
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "GO_BP_sim",
#'   group_by = "CellType",
#'   plot_type = "comparison"
#' )
#'
#' # Use a combined database
#' pancreas_sub <- RunGSEA(
#'   srt = pancreas_sub, group_by = "CellType",
#'   db = c("KEGG", "WikiPathway", "Reactome", "PFAM", "MP"),
#'   db_combine = TRUE,
#'   species = "Mus_musculus"
#' )
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "Combined",
#'   group_by = "CellType",
#'   plot_type = "comparison"
#' )
#'
#' # Or use "geneID", "geneScore" and "geneID_groups" as input to run GSEA
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
#' @importFrom BiocParallel bplapply bpprogressbar<- bpRNGseed<- bpworkers ipcid ipclock ipcunlock
#' @importFrom clusterProfiler GSEA simplify
#' @export
RunGSEA <- function(
    srt = NULL,
    group_by = NULL,
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
  message(paste0("[", time_start, "] ", "Start GSEA"))
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
    geneScore <- de_df[["avg_log2FC"]]
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
  if (length(geneScore) != length(geneID)) {
    stop("geneScore must be the same length with geneID")
  }
  if (all(geneScore > 0) && scoreType != "pos") {
    scoreType <- "pos"
    warning(
      "All values in the geneScore are greater than zero. Set scoreType = 'pos'.",
      immediate. = TRUE
    )
  }
  if (all(geneScore < 0) && scoreType != "neg") {
    scoreType <- "neg"
    warning(
      "All values in the geneScore are less than zero. Set scoreType = 'neg'.",
      immediate. = TRUE
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
    message("Ignore ", length(na_index), " NA geneScore")
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
    geneMap <- data.frame(
      from_geneID = unique(geneID),
      row.names = unique(geneID)
    )
    colnames(geneMap)[1] <- IDtype
  }

  input[[IDtype]] <- geneMap[as.character(input$geneID), IDtype]
  input[[result_IDtype]] <- geneMap[as.character(input$geneID), result_IDtype]
  input <- unnest(input, cols = c(IDtype, result_IDtype))
  input <- input[!is.na(input[[IDtype]]), , drop = FALSE]

  message("Permform GSEA...")
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
      na <- rowSums(is.na(TERM2GENE_tmp)) > 0
      TERM2GENE_tmp <- TERM2GENE_tmp[!(dup | na), , drop = FALSE]
      TERM2NAME_tmp <- TERM2NAME_tmp[
        TERM2NAME_tmp[["Term"]] %in% TERM2GENE_tmp[["Term"]], ,
        drop = FALSE
      ]
      enrich_res <- GSEA(
        geneList = geneList,
        minGSSize = ifelse(term %in% unlimited_db, 1, minGSSize),
        maxGSSize = ifelse(term %in% unlimited_db, Inf, maxGSSize),
        nPermSimple = 1e5, # nPermSimple:fgseaMultilevel; nperm:fgseaSimple
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

        if (
          isTRUE(GO_simplify) && term %in% c("GO", "GO_BP", "GO_CC", "GO_MF")
        ) {
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
  message(paste0("[", time_end, "] ", "GSEA done"))
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
    srt@tools[[paste("GSEA", group_by, test.use, sep = "_")]] <- res
    return(srt)
  } else {
    return(res)
  }
}

#' GSEA Plot
#'
#' This function generates various types of plots for Gene Set Enrichment Analysis (GSEA) results.
#'
#' @inheritParams EnrichmentPlot
#' @param srt A Seurat object containing the results of RunDEtest and RunGSEA.
#' If specified, GSEA results will be extracted from the Seurat object automatically.
#' If not specified, the \code{res} arguments must be provided.
#' @param res Enrichment results generated by RunGSEA function. If provided, 'srt', 'test.use' and 'group_by' are ignored.
#' @param plot_type The type of plot to generate. Options are: "line", "comparison", "bar", "network", "enrichmap", "wordcloud". Default is "line".
#' @param direction The direction of enrichment to include in the plot. Must be one of "pos", "neg", or "both". The default value is "both".
#' @param line_width The linewidth for the line plot.
#' @param line_alpha The alpha value for the line plot.
#' @param line_color The color for the line plot.
#' @param n_coregene The number of core genes to label in the line plot.
#' @param sample_coregene Whether to randomly sample core genes for labeling in the line plot.
#' @param features_label A character vector of feature names to include as labels in the line plot.
#' @param label.fg The color of the labels.
#' @param label.bg The background color of the labels.
#' @param label.bg.r The radius of the rounding of the label's background.
#' @param label.size The size of the labels.
#'
#' @seealso \code{\link{RunGSEA}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunDEtest(
#'   pancreas_sub,
#'   group_by = "CellType",
#'   only.pos = FALSE,
#'   fc.threshold = 1
#' )
#' pancreas_sub <- RunGSEA(
#'   pancreas_sub,
#'   group_by = "CellType",
#'   db = "GO_BP",
#'   species = "Mus_musculus"
#' )
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Ductal"
#' )
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Ductal",
#'   id_use = "GO:0006412"
#' )
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Endocrine",
#'   id_use = c("GO:0046903", "GO:0015031", "GO:0007600")
#' ) %>%
#'   panel_fix_overall(height = 6)
#' # As the plot is created by combining,
#' # we can adjust the overall height and width directly.
#'
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   topTerm = 3,
#'   plot_type = "comparison"
#' )
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   topTerm = 3,
#'   plot_type = "comparison",
#'   direction = "neg"
#' )
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   topTerm = 3,
#'   plot_type = "comparison",
#'   compare_only_sig = TRUE
#' )
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   plot_type = "bar"
#' )
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   plot_type = "bar",
#'   direction = "both"
#' )
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Ductal",
#'   plot_type = "bar",
#'   topTerm = 20,
#'   direction = "both",
#'   palcolor = c("red3", "steelblue")
#' )
#'
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Endocrine",
#'   plot_type = "network"
#' )
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Endocrine",
#'   plot_type = "enrichmap"
#' )
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Endocrine",
#'   plot_type = "wordcloud"
#' )
#' GSEAPlot(
#'   pancreas_sub,
#'   db = "GO_BP",
#'   group_by = "CellType",
#'   group_use = "Endocrine",
#'   plot_type = "wordcloud",
#'   word_type = "feature"
#' )
#'
#' @importFrom ggplot2 ggplot aes theme theme_classic alpha element_blank element_rect margin geom_line geom_point geom_rect geom_linerange geom_hline geom_vline geom_segment annotate ggtitle labs xlab ylab scale_x_continuous scale_y_continuous scale_color_manual scale_alpha_manual guides guide_legend guide_none
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices colorRamp
#' @importFrom patchwork wrap_plots
#' @importFrom dplyr case_when filter pull %>%
#' @importFrom stats quantile
#' @importFrom gtable gtable_add_rows gtable_add_grob
#' @importFrom grid textGrob
#' @export
GSEAPlot <- function(
    srt,
    db = "GO_BP",
    group_by = NULL,
    test.use = "wilcox",
    res = NULL,
    plot_type = c(
      "line",
      "bar",
      "network",
      "enrichmap",
      "wordcloud",
      "comparison"
    ),
    group_use = NULL,
    id_use = NULL,
    pvalueCutoff = NULL,
    padjustCutoff = 0.05,
    topTerm = ifelse(plot_type == "enrichmap", 100, 6),
    direction = c("pos", "neg", "both"),
    compare_only_sig = FALSE,
    topWord = 100,
    word_type = c("term", "feature"),
    word_size = c(2, 8),
    words_excluded = NULL,
    line_width = 1.5,
    line_alpha = 1,
    line_color = "#6BB82D",
    n_coregene = 10,
    sample_coregene = FALSE,
    features_label = NULL,
    label.fg = "black",
    label.bg = "white",
    label.bg.r = 0.1,
    label.size = 4,
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
    aspect.ratio = NULL,
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
  direction <- match.arg(direction)
  enrichmap_label <- match.arg(enrichmap_label)
  enrichmap_mark <- match.arg(enrichmap_mark)
  words_excluded <- words_excluded %||% scop::words_excluded

  subplots <- 1:3
  rel_heights <- c(1.5, 0.5, 1)
  rel_width <- 3

  if (is.null(res)) {
    if (is.null(group_by)) {
      stop("'group_by' must be provided.")
    }
    layer <- paste("GSEA", group_by, test.use, sep = "_")
    if (!layer %in% names(srt@tools)) {
      stop("No enrichment result found. You may perform RunGSEA first.")
    }
    enrichment <- srt@tools[[layer]][["enrichment"]]
    res <- srt@tools[[layer]][["results"]]
  } else {
    enrichment <- res[["enrichment"]]
    res <- res[["results"]]
  }
  group_use <- group_use %||% unique(enrichment[["Groups"]])
  comb <- expand.grid(group_use, db)
  use <- names(res)[names(res) %in% paste(comb$Var1, comb$Var2, sep = "-")]
  if (length(use) == 0) {
    stop(paste0(db, " is not in the enrichment result."))
  }
  res <- res[use]
  enrichment <- enrichment[
    enrichment[["Groups"]] %in% group_use, ,
    drop = FALSE
  ]

  if (is.null(pvalueCutoff) && is.null(padjustCutoff)) {
    stop("One of 'pvalueCutoff' or 'padjustCutoff' must be specified")
  }
  if (!is.factor(enrichment[["Database"]])) {
    enrichment[["Database"]] <- factor(
      enrichment[["Database"]],
      levels = unique(enrichment[["Database"]])
    )
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

  pvalueCutoff <- ifelse(is.null(pvalueCutoff), 1, pvalueCutoff)
  padjustCutoff <- ifelse(is.null(padjustCutoff), 1, padjustCutoff)

  if (any(db %in% c("GO_sim", "GO_BP_sim", "GO_CC_sim", "GO_MF_sim"))) {
    enrichment_sim <- enrichment[
      enrichment[["Database"]] %in% gsub("_sim", "", db), ,
      drop = FALSE
    ]
  }
  enrichment <- enrichment[enrichment[["Database"]] %in% db, , drop = FALSE]

  plist <- NULL
  if (plot_type == "comparison") {
    # comparison -------------------------------------------------------------------------------------------------
    if (length(id_use) > 0) {
      ids <- unlist(id_use)
    } else {
      ids <- NULL
      for (i in group_use) {
        df <- enrichment[enrichment[["Groups"]] == i, , drop = FALSE]
        df <- df[df[[metric]] < metric_value, , drop = FALSE]
        df <- df[order(df[[metric]]), , drop = FALSE]
        df_up <- df[df[["NES"]] > 0, , drop = FALSE]
        ID_up <- df_up[head(order(df_up[[metric]]), topTerm), "ID"]
        df_down <- df[df[["NES"]] < 0, , drop = FALSE]
        ID_down <- df_down[head(order(df_down[[metric]]), topTerm), "ID"]
        ids <- switch(direction,
          "pos" = unique(c(ids, head(ID_up, topTerm))),
          "neg" = unique(c(ids, head(ID_down, topTerm))),
          "both" = unique(c(
            ids,
            head(
              c(
                head(ID_up, ceiling(topTerm / 2)),
                head(ID_down, ceiling(topTerm / 2))
              ),
              topTerm
            )
          ))
        )
      }
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
    enrichment_sub[["Significant"]] <- enrichment_sub[[metric]] < metric_value
    enrichment_sub[["Significant"]] <- factor(
      enrichment_sub[["Significant"]],
      levels = c("TRUE", "FALSE")
    )
    if (isTRUE(compare_only_sig)) {
      enrichment_sub <- enrichment_sub[
        enrichment_sub[["Significant"]] == "TRUE", ,
        drop = FALSE
      ]
    }
    enrichment_sub <- switch(direction,
      "pos" = enrichment_sub[enrichment_sub[["NES"]] > 0, , drop = FALSE],
      "neg" = enrichment_sub[enrichment_sub[["NES"]] < 0, , drop = FALSE],
      "both" = enrichment_sub
    )

    p <- ggplot(enrichment_sub, aes(x = Groups, y = Description)) +
      geom_point(
        aes(size = setSize, fill = NES, color = Significant),
        shape = 21,
        stroke = 0.8
      ) +
      scale_size_area(name = "setSize", max_size = 6, n.breaks = 4) +
      guides(
        size = guide_legend(
          override.aes = list(fill = "grey30", shape = 21),
          order = 2
        )
      ) +
      scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5)) +
      scale_fill_gradientn(
        name = "NES",
        n.breaks = 4,
        limits = c(
          -max(abs(enrichment_sub[["NES"]])),
          max(abs(enrichment_sub[["NES"]]))
        ),
        colors = palette_scop(palette = palette, palcolor = palcolor),
        guide = guide_colorbar(
          frame.colour = "black",
          ticks.colour = "black",
          title.hjust = 0,
          order = 1
        )
      ) +
      scale_color_manual(
        name = paste0(
          "Significant\n(",
          metric,
          "<",
          metric_value,
          ")",
          collapse = ""
        ),
        values = c("TRUE" = "black", "FALSE" = "grey90"),
        guide = if (isTRUE(compare_only_sig)) guide_none() else guide_legend()
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
  } else if (plot_type == "line") {
    # line -------------------------------------------------------------------------------------------------
    for (nm in names(res)) {
      res_enrich <- res[[nm]]
      if (is.null(id_use)) {
        geneSetID_filter <- res_enrich@result[
          res_enrich@result[[metric]] < metric_value, ,
          drop = FALSE
        ]
        geneSetID_filter <- geneSetID_filter[
          order(geneSetID_filter[[metric]]), ,
          drop = FALSE
        ]
        geneSetID_up <- geneSetID_filter[
          geneSetID_filter[["NES"]] > 0, ,
          drop = FALSE
        ]
        geneSetID_up <- geneSetID_up[
          head(order(geneSetID_up[[metric]]), topTerm),
          "ID"
        ]
        geneSetID_down <- geneSetID_filter[
          geneSetID_filter[["NES"]] < 0, ,
          drop = FALSE
        ]
        geneSetID_down <- geneSetID_down[
          head(order(geneSetID_down[[metric]]), topTerm),
          "ID"
        ]
        geneSetID_use <- switch(direction,
          "pos" = unique(head(geneSetID_up, topTerm)),
          "neg" = unique(head(geneSetID_down, topTerm)),
          "both" = unique(head(
            c(
              head(geneSetID_up, ceiling(topTerm / 2)),
              head(geneSetID_down, ceiling(topTerm / 2))
            ),
            topTerm
          ))
        )
      } else {
        if (is.list(id_use)) {
          geneSetID_use <- intersect(
            res_enrich@result[["ID"]],
            id_use[[unique(res_enrich@result$Groups)]]
          )
        } else {
          geneSetID_use <- id_use
        }
      }
      if (length(geneSetID_use) == 1) {
        gsdata <- gsInfo(object = res_enrich, id_use = geneSetID_use)
      } else {
        gsdata <- do.call(
          rbind,
          lapply(geneSetID_use, gsInfo, object = res_enrich)
        )
      }
      if (length(geneSetID_use) == 0) {
        plist[[nm]] <- NULL
        next
      }
      stat <- res_enrich[geneSetID_use, c("Description", "NES", metric)]
      rownames(stat) <- stat[, "Description"]
      stat$p.sig <- case_when(
        stat[[metric]] > 0.05 ~ "ns  ",
        stat[[metric]] <= 0.05 & stat[[metric]] > 0.01 ~ "*   ",
        stat[[metric]] <= 0.01 & stat[[metric]] > 0.001 ~ "**  ",
        stat[[metric]] <= 0.001 & stat[[metric]] > 0.0001 ~ "*** ",
        stat[[metric]] <= 0.0001 ~ "****"
      )
      gsdata[["NES"]] <- stat[gsdata$Description, "NES"]
      gsdata[[metric]] <- stat[gsdata$Description, metric]
      gsdata[["p.sig"]] <- stat[gsdata$Description, "p.sig"]
      gsdata[["DescriptionP"]] <- capitalize(gsdata[["Description"]])
      gsdata[["DescriptionP"]] <- str_wrap(
        gsdata[["DescriptionP"]],
        width = character_width
      )
      gsdata[["DescriptionP"]] <- paste0(
        gsdata[["DescriptionP"]],
        "\n(NES=",
        round(gsdata[["NES"]], 3),
        ", ",
        metric,
        "=",
        format(gsdata[[metric]], digits = 3, scientific = TRUE),
        ", ",
        gsdata[["p.sig"]],
        ")"
      )
      gsdata[["DescriptionP"]] <- factor(
        gsdata[["DescriptionP"]],
        levels = unique(gsdata[["DescriptionP"]])
      )
      p <- ggplot(gsdata, aes(x = x)) +
        xlab(NULL) +
        theme_classic(base_size = 12) +
        theme(
          panel.grid.major = element_line(colour = "grey90", linetype = 2),
          panel.grid.minor = element_line(colour = "grey90", linetype = 2)
        ) +
        scale_x_continuous(expand = c(0.01, 0))
      es_layer <- geom_line(
        aes(y = runningScore, color = DescriptionP),
        linewidth = line_width,
        alpha = line_alpha
      )
      bg_dat <- data.frame(
        xmin = -Inf,
        xmax = Inf,
        ymin = c(0, -Inf),
        ymax = c(Inf, 0),
        fill = c(alpha("#C40003", 0.2), alpha("#1D008F", 0.2))
      )
      p1 <- p +
        geom_rect(
          data = bg_dat,
          mapping = aes(
            xmin = xmin,
            xmax = xmax,
            ymin = ymin,
            ymax = ymax,
            fill = I(fill)
          ),
          inherit.aes = FALSE
        ) +
        geom_hline(yintercept = 0, linetype = 1, color = "grey40") +
        es_layer +
        ylab("Enrichment Score") +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line = element_blank(),
          panel.border = element_rect(
            color = "black",
            fill = "transparent",
            linewidth = 1
          ),
          plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, unit = "cm"),
          legend.position = "right",
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent")
        )

      i <- 0
      for (term in rev(levels(gsdata$DescriptionP))) {
        idx <- which(gsdata$ymin != 0 & gsdata$DescriptionP == term)
        gsdata[idx, "ymin"] <- i
        gsdata[idx, "ymax"] <- i + 1
        i <- i + 1
      }
      p2 <- ggplot(gsdata, aes(x = x)) +
        geom_linerange(
          aes(ymin = ymin, ymax = ymax, color = DescriptionP),
          alpha = line_alpha
        ) +
        xlab(NULL) +
        ylab(NULL) +
        theme_classic(base_size = 12) +
        theme(
          legend.position = "none",
          plot.margin = margin(t = -0.1, b = 0, r = 0.2, l = 0.2, unit = "cm"),
          panel.border = element_rect(
            color = "black",
            fill = "transparent",
            linewidth = 1
          ),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank()
        ) +
        scale_x_continuous(expand = c(0.01, 0)) +
        scale_y_continuous(expand = c(0, 0))
      if (length(geneSetID_use) == 1) {
        subtitle_use <- paste0(
          "(NES=",
          round(stat[["NES"]], 3),
          ", ",
          metric,
          "=",
          format(stat[[metric]], digits = 3, scientific = TRUE),
          ", ",
          stat[["p.sig"]],
          ")"
        )
        p1 <- p1 +
          annotate(
            geom = "segment",
            x = 0,
            xend = p$data$x[which.max(abs(p$data$runningScore))],
            y = p$data$runningScore[which.max(abs(p$data$runningScore))],
            yend = p$data$runningScore[which.max(abs(p$data$runningScore))],
            linetype = 2
          ) +
          annotate(
            geom = "segment",
            x = p$data$x[which.max(abs(p$data$runningScore))],
            xend = p$data$x[which.max(abs(p$data$runningScore))],
            y = 0,
            yend = p$data$runningScore[which.max(abs(p$data$runningScore))],
            linetype = 2
          ) +
          annotate(
            geom = "point",
            x = p$data$x[which.max(abs(p$data$runningScore))],
            y = p$data$runningScore[which.max(abs(p$data$runningScore))],
            fill = ifelse(stat[["NES"]] < 0, "#5E34F5", "#F52323"),
            color = "black",
            size = 2.5,
            shape = ifelse(stat[["NES"]] < 0, 25, 24)
          ) +
          labs(subtitle = subtitle_use) +
          theme(plot.subtitle = element_text(face = "italic"))

        if (
          (is.numeric(n_coregene) && n_coregene > 1) ||
            length(features_label) > 0
        ) {
          if (length(features_label) == 0) {
            features_label_tmp <- unlist(strsplit(gsdata$CoreGene[1], "/"))
            n_coregene <- min(n_coregene, length(features_label_tmp))
            if (isTRUE(sample_coregene)) {
              features_label_tmp <- sample(
                features_label_tmp,
                n_coregene,
                replace = FALSE
              )
            } else {
              features_label_tmp <- gsdata$GeneName[
                gsdata$GeneName %in% features_label_tmp
              ][1:n_coregene]
            }
          } else {
            features_label_tmp <- features_label
          }
          df_gene <- gsdata[
            gsdata$position == 1 & gsdata$GeneName %in% features_label_tmp, ,
            drop = FALSE
          ]
          gene_drop <- features_label_tmp[
            !features_label_tmp %in% df_gene$GeneName
          ]
          if (length(gene_drop) > 0) {
            warning(
              "Gene ",
              paste(gene_drop, collapse = ","),
              " is not in the geneset of the ",
              gsdata$Description[1],
              immediate. = TRUE
            )
          }
          x_nudge <- diff(range(gsdata$x)) * 0.05
          y_nudge <- diff(range(gsdata$runningScore)) * 0.05
          p1 <- p1 +
            geom_point(
              data = df_gene,
              mapping = aes(y = runningScore),
              color = "black"
            ) +
            geom_text_repel(
              data = df_gene,
              mapping = aes(y = runningScore, label = GeneName),
              min.segment.length = 0,
              max.overlaps = 100,
              segment.colour = "grey40",
              color = label.fg,
              bg.color = label.bg,
              bg.r = label.bg.r,
              size = label.size,
              nudge_x = ifelse(df_gene$runningScore >= 0, x_nudge, -x_nudge),
              nudge_y = ifelse(df_gene$runningScore > 0, -y_nudge, y_nudge)
            )
        }

        x <- p$data$x
        y <- y_raw <- p$data$geneList
        y[y > quantile(y_raw, 0.98)] <- quantile(y_raw, 0.98)
        y[y < quantile(y_raw, 0.02)] <- quantile(y_raw, 0.02)
        col <- rep("white", length(y))
        y_pos <- which(y > 0)
        if (length(y_pos) > 0) {
          y_pos_i <- cut(
            y[y_pos],
            breaks = seq(
              min(y[y_pos], na.rm = TRUE),
              max(y[y_pos], na.rm = TRUE),
              len = 100
            ),
            include.lowest = TRUE
          )
          col[y_pos] <- colorRampPalette(c("#F5DCDC", "#C40003"))(100)[y_pos_i]
        }

        y_neg <- which(y < 0)
        if (length(y_neg) > 0) {
          y_neg_i <- cut(
            y[y_neg],
            breaks = seq(
              min(y[y_neg], na.rm = TRUE),
              max(y[y_neg], na.rm = TRUE),
              len = 100
            ),
            include.lowest = TRUE
          )
          col[y_neg] <- colorRampPalette(c("#1D008F", "#DDDCF5"))(100)[y_neg_i]
        }

        ymin <- min(p2$data$ymin, na.rm = TRUE)
        ymax <- max(p2$data$ymax - p2$data$ymin, na.rm = TRUE) * 0.3
        xmin <- which(!duplicated(col))
        xmax <- xmin + as.numeric(table(col)[as.character(unique(col))])
        d <- data.frame(
          ymin = ymin,
          ymax = ymax,
          xmin = xmin,
          xmax = xmax,
          col = unique(col)
        )
        p2 <- p2 +
          geom_rect(
            aes(
              xmin = xmin,
              xmax = xmax,
              ymin = ymin,
              ymax = ymax,
              fill = I(col)
            ),
            data = d,
            alpha = 0.95,
            inherit.aes = FALSE
          )
      }
      df2 <- p$data
      df2$y <- p$data$geneList[df2$x]
      min_y <- df2$y[which.min(abs(df2$y))]
      corss_x <- median(df2$x[df2$y == min_y])
      p3 <- p +
        geom_segment(
          data = df2,
          aes(
            x = x,
            xend = x,
            y = y,
            yend = 0
          ),
          color = "grey30"
        )

      if (max(df2$y) > 0) {
        p3 <- p3 +
          annotate(
            geom = "text",
            x = 0,
            y = Inf,
            vjust = 1.3,
            hjust = 0,
            color = "#C81A1F",
            size = 4,
            label = " Positively correlated"
          )
      }
      if (min(df2$y) < 0) {
        p3 <- p3 +
          annotate(
            geom = "text",
            x = Inf,
            y = -Inf,
            vjust = -0.3,
            hjust = 1,
            color = "#3C298C",
            size = 4,
            label = "Negtively correlated "
          )
      }
      if (max(df2$y) > 0 && min(df2$y) < 0) {
        p3 <- p3 +
          geom_vline(xintercept = corss_x, linetype = 2, color = "black") +
          annotate(
            geom = "text",
            y = 0,
            x = corss_x,
            vjust = ifelse(diff(abs(range(df2$y))) > 0, -0.3, 1.3),
            size = 4,
            label = paste0("Zero cross at ", corss_x)
          )
      }
      p3 <- p3 +
        ylab("Ranked List Metric") +
        xlab("Rank in Ordered Dataset") +
        theme(
          plot.margin = margin(
            t = -0.1,
            r = 0.2,
            b = 0.2,
            l = 0.2,
            unit = "cm"
          ),
          axis.line = element_blank(),
          axis.line.x = element_blank(),
          panel.border = element_rect(
            color = "black",
            fill = "transparent",
            linewidth = 1
          )
        )
      if (length(geneSetID_use) == 1) {
        p1 <- p1 + ggtitle(gsdata$Description[1], subtitle = subtitle_use)
      }
      if (length(line_color) != length(geneSetID_use)) {
        color_use <- palette_scop(
          levels(gsdata$DescriptionP),
          palette = palette,
          palcolor = palcolor
        )
      } else {
        color_use <- line_color
      }
      p1 <- p1 + scale_color_manual(values = color_use)
      if (length(color_use) == 1) {
        p1 <- p1 + theme(legend.position = "none")
        p2 <- p2 + scale_color_manual(values = "black")
      } else {
        p2 <- p2 + scale_color_manual(values = color_use)
      }
      legend <- get_legend(
        p1 +
          guides(color = guide_legend(title = "Term:", byrow = TRUE)) +
          do.call(theme_use, theme_args) +
          theme(
            legend.position = legend.position,
            legend.direction = legend.direction
          )
      )
      plotlist <- list(p1 + theme(legend.position = "none"), p2, p3)[subplots]
      if (length(subplots) == 1) {
        plist[[nm]] <- plotlist[[1]] +
          theme(
            aspect.ratio = rel_heights[subplots] / rel_width,
            plot.margin = margin(
              t = 0.2,
              r = 0.2,
              b = 0.2,
              l = 0.2,
              unit = "cm"
            )
          )
      } else {
        plotlist <- lapply(plotlist[subplots], as_grob)
        rel_heights <- rel_heights[subplots]
        for (i in seq_along(plotlist)) {
          plotlist[[i]] <- panel_fix_overall(
            plotlist[[i]],
            height = rel_heights[i],
            units = "null",
            margin = 0,
            respect = TRUE,
            return_grob = TRUE
          )
          plotlist[[i]] <- panel_fix_overall(
            plotlist[[i]],
            width = rel_width,
            units = "null",
            margin = 0,
            respect = TRUE,
            return_grob = TRUE
          )
        }
        p_out <- do.call(rbind, c(plotlist, size = "first"))

        if (length(geneSetID_use) > 1) {
          p_out <- add_grob(p_out, legend, legend.position)
        }
        lab <- textGrob(label = nm, rot = -90, hjust = 0.5)
        p_out <- add_grob(p_out, lab, "right", clip = "off")
        p_out <- wrap_plots(p_out)
        plist[[nm]] <- p_out
      }
    }
  } else if (plot_type == "bar") {
    # bar -------------------------------------------------------------------------------------------------
    for (nm in names(res)) {
      res_enrich <- res[[nm]]
      if (is.null(id_use)) {
        geneSetID_filter <- res_enrich@result[
          res_enrich@result[[metric]] < metric_value, ,
          drop = FALSE
        ]
        geneSetID_filter <- geneSetID_filter[
          order(geneSetID_filter[[metric]]), ,
          drop = FALSE
        ]
        geneSetID_up <- geneSetID_filter[
          geneSetID_filter[["NES"]] > 0, ,
          drop = FALSE
        ]
        geneSetID_up <- geneSetID_up[
          head(order(geneSetID_up[[metric]]), topTerm),
          "ID"
        ]
        geneSetID_down <- geneSetID_filter[
          geneSetID_filter[["NES"]] < 0, ,
          drop = FALSE
        ]
        geneSetID_down <- geneSetID_down[
          head(order(geneSetID_down[[metric]]), topTerm),
          "ID"
        ]
        geneSetID_use <- switch(direction,
          "pos" = unique(head(geneSetID_up, topTerm)),
          "neg" = unique(head(geneSetID_down, topTerm)),
          "both" = unique(head(
            c(
              head(geneSetID_up, ceiling(topTerm / 2)),
              head(geneSetID_down, ceiling(topTerm / 2))
            ),
            topTerm
          ))
        )
      } else {
        if (is.list(id_use)) {
          geneSetID_use <- intersect(
            res_enrich@result[["ID"]],
            id_use[[unique(res_enrich@result$Groups)]]
          )
        } else {
          geneSetID_use <- id_use
        }
      }
      if (length(geneSetID_use) == 0) {
        plist[[nm]] <- NULL
        next
      }
      stat <- res_enrich[geneSetID_use, , drop = FALSE]
      stat <- stat[order(stat[["NES"]]), , drop = FALSE]
      rownames(stat) <- stat[, "Description"]
      stat[["Description"]] <- capitalize(stat[["Description"]])
      stat[["Description"]] <- str_wrap(
        stat[["Description"]],
        width = character_width
      )
      stat[["Description"]] <- factor(
        stat[["Description"]],
        levels = unique(stat[["Description"]])
      )
      stat[["Direction"]] <- ifelse(stat[["NES"]] > 0, "Pos", "Neg")
      stat[["Direction"]] <- factor(
        stat[["Direction"]],
        levels = c("Pos", "Neg")
      )

      p <- ggplot(
        stat,
        aes(
          x = .data[["NES"]],
          y = .data[["Description"]]
        )
      ) +
        geom_vline(xintercept = 0) +
        geom_col(
          aes(fill = .data[["Direction"]], alpha = -log10(.data[[metric]])),
          color = "black"
        ) +
        geom_text(
          aes(
            x = 0,
            y = .data[["Description"]],
            label = .data[["Description"]],
            hjust = ifelse(.data[["NES"]] > 0, 1, 0),
          ),
          nudge_x = ifelse(stat[["NES"]] > 0, -0.05, 0.05),
          lineheight = lineheight,
          fontface = ifelse(
            grepl("\n", levels(stat[["Description"]])),
            "italic",
            "plain"
          )
        ) +
        scale_fill_manual(
          values = palette_scop(
            x = rev(levels(stat[["Direction"]])),
            palette = palette,
            palcolor = rev(palcolor)
          ),
          guide = if (direction == "both") {
            guide_legend(order = 1)
          } else {
            guide_none()
          }
        ) +
        facet_grid(Database ~ Groups, scales = "free") +
        coord_cartesian(
          xlim = c(-max(abs(stat[["NES"]])), max(abs(stat[["NES"]])))
        ) +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction,
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        )
      plist[[nm]] <- p
    }
  } else if (plot_type == "network") {
    # network -------------------------------------------------------------------------------------------------
    for (nm in names(res)) {
      res_enrich <- res[[nm]]
      if (is.null(id_use)) {
        geneSetID_filter <- res_enrich@result[
          res_enrich@result[[metric]] < metric_value, ,
          drop = FALSE
        ]
        geneSetID_filter <- geneSetID_filter[
          order(geneSetID_filter[[metric]]), ,
          drop = FALSE
        ]
        geneSetID_up <- geneSetID_filter[
          geneSetID_filter[["NES"]] > 0, ,
          drop = FALSE
        ]
        geneSetID_up <- geneSetID_up[
          head(order(geneSetID_up[[metric]]), topTerm),
          "ID"
        ]
        geneSetID_down <- geneSetID_filter[
          geneSetID_filter[["NES"]] < 0, ,
          drop = FALSE
        ]
        geneSetID_down <- geneSetID_down[
          head(order(geneSetID_down[[metric]]), topTerm),
          "ID"
        ]
        geneSetID_use <- switch(direction,
          "pos" = unique(head(geneSetID_up, topTerm)),
          "neg" = unique(head(geneSetID_down, topTerm)),
          "both" = unique(head(
            c(
              head(geneSetID_up, ceiling(topTerm / 2)),
              head(geneSetID_down, ceiling(topTerm / 2))
            ),
            topTerm
          ))
        )
      } else {
        if (is.list(id_use)) {
          geneSetID_use <- intersect(
            res_enrich@result[["ID"]],
            id_use[[unique(res_enrich@result$Groups)]]
          )
        } else {
          geneSetID_use <- id_use
        }
      }
      if (length(geneSetID_use) == 0) {
        plist[[nm]] <- NULL
        next
      }
      df <- res_enrich[geneSetID_use, , drop = FALSE]
      df$p.sig <- case_when(
        df[[metric]] > 0.05 ~ "ns  ",
        df[[metric]] <= 0.05 & df[[metric]] > 0.01 ~ "*   ",
        df[[metric]] <= 0.01 & df[[metric]] > 0.001 ~ "**  ",
        df[[metric]] <= 0.001 & df[[metric]] > 0.0001 ~ "*** ",
        df[[metric]] <= 0.0001 ~ "****"
      )
      df[["metric"]] <- -log10(df[[metric]])
      df[["Description"]] <- capitalize(df[["Description"]])
      df[["Description"]] <- str_wrap(
        df[["Description"]],
        width = character_width
      )
      df[["Description"]] <- paste0(
        df[["Description"]],
        "\n(NES=",
        round(df[["NES"]], 3),
        ", ",
        metric,
        "=",
        format(df[[metric]], digits = 3, scientific = TRUE),
        ", ",
        df[["p.sig"]],
        ")"
      )
      df[["Description"]] <- factor(
        df[["Description"]],
        levels = unique(df[["Description"]])
      )
      df[["geneID"]] <- strsplit(df[["core_enrichment"]], "/")
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
        guides(fill = guide_legend(title = "Term:", byrow = TRUE)) +
        labs(x = "", y = "") +
        facet_grid(Database ~ Groups, scales = "free") +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction
        )
      plist[[nm]] <- p
    }
  } else if (plot_type == "enrichmap") {
    # enrichmap -------------------------------------------------------------------------------------------------
    for (nm in names(res)) {
      res_enrich <- res[[nm]]
      if (is.null(id_use)) {
        geneSetID_filter <- res_enrich@result[
          res_enrich@result[[metric]] < metric_value, ,
          drop = FALSE
        ]
        geneSetID_filter <- geneSetID_filter[
          order(geneSetID_filter[[metric]]), ,
          drop = FALSE
        ]
        geneSetID_up <- geneSetID_filter[
          geneSetID_filter[["NES"]] > 0, ,
          drop = FALSE
        ]
        geneSetID_up <- geneSetID_up[
          head(order(geneSetID_up[[metric]]), topTerm),
          "ID"
        ]
        geneSetID_down <- geneSetID_filter[
          geneSetID_filter[["NES"]] < 0, ,
          drop = FALSE
        ]
        geneSetID_down <- geneSetID_down[
          head(order(geneSetID_down[[metric]]), topTerm),
          "ID"
        ]
        geneSetID_use <- switch(direction,
          "pos" = unique(head(geneSetID_up, topTerm)),
          "neg" = unique(head(geneSetID_down, topTerm)),
          "both" = unique(head(
            c(
              head(geneSetID_up, ceiling(topTerm / 2)),
              head(geneSetID_down, ceiling(topTerm / 2))
            ),
            topTerm
          ))
        )
      } else {
        if (is.list(id_use)) {
          geneSetID_use <- intersect(
            res_enrich@result[["ID"]],
            id_use[[unique(res_enrich@result$Groups)]]
          )
        } else {
          geneSetID_use <- id_use
        }
      }
      if (length(geneSetID_use) == 0) {
        plist[[nm]] <- NULL
        next
      }
      df <- res_enrich[geneSetID_use, , drop = FALSE]
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
      df[["Direction"]] <- ifelse(df[["NES"]] > 0, "Pos", "Neg")
      df[["Direction"]] <- factor(df[["Direction"]], levels = c("Pos", "Neg"))
      df[["geneID"]] <- strsplit(df[["core_enrichment"]], "/")
      df[["Count"]] <- sapply(df[["geneID"]], length)
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
        facet_grid(Database ~ Groups, scales = "free") +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction
        )
      plist[[nm]] <- p
    }
  } else if (plot_type == "wordcloud") {
    # wordcloud -------------------------------------------------------------------------------------------------
    check_r("ggwordcloud")
    check_r("jokergoo/simplifyEnrichment")
    for (nm in names(res)) {
      res_enrich <- res[[nm]]
      if (is.null(id_use)) {
        geneSetID_filter <- res_enrich@result[
          res_enrich@result[[metric]] < metric_value, ,
          drop = FALSE
        ]
        geneSetID_filter <- geneSetID_filter[
          order(geneSetID_filter[[metric]]), ,
          drop = FALSE
        ]
        geneSetID_up <- geneSetID_filter[geneSetID_filter[["NES"]] > 0, "ID"]
        geneSetID_down <- geneSetID_filter[geneSetID_filter[["NES"]] < 0, "ID"]
        geneSetID_use <- switch(direction,
          "pos" = unique(geneSetID_up),
          "neg" = unique(geneSetID_down),
          "both" = unique(c(geneSetID_up, geneSetID_down))
        )
      } else {
        if (is.list(id_use)) {
          geneSetID_use <- intersect(
            res_enrich@result[["ID"]],
            id_use[[unique(res_enrich@result$Groups)]]
          )
        } else {
          geneSetID_use <- id_use
        }
      }
      if (length(geneSetID_use) == 0) {
        plist[[nm]] <- NULL
        next
      }
      df <- res_enrich[geneSetID_use, , drop = FALSE]

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
          mutate(
            keyword = strsplit(as.character(.data[["core_enrichment"]]), "/")
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
        facet_grid(Database ~ Groups, scales = "free") +
        coord_flip() +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction
        )
      plist[[nm]] <- p
    }
  }

  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- wrap_plots(plotlist = plist, nrow = nrow, ncol = ncol)
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
}

gsInfo <- function(object, id_use) {
  geneList <- object@geneList
  if (is.numeric(id_use)) {
    id_use <- object@result[id_use, "ID"]
  }
  geneSet <- object@geneSets[[id_use]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore)) / 20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  df$Description <- object@result[id_use, "Description"]
  df$CoreGene <- object@result[id_use, "core_enrichment"]
  if (length(object@gene2Symbol) == length(object@geneList)) {
    df$GeneName <- object@gene2Symbol
  } else {
    df$GeneName <- df$gene
  }
  return(df)
}

gseaScores <- function(geneList, geneSet, exponent = 1) {
  geneSet <- intersect(geneSet, names(geneList))
  N <- length(geneList)
  Nh <- length(geneSet)
  Phit <- Pmiss <- numeric(N)
  hits <- names(geneList) %in% geneSet
  Phit[hits] <- abs(geneList[hits])^exponent
  NR <- sum(Phit)
  Phit <- cumsum(Phit / NR)
  Pmiss[!hits] <- 1 / (N - Nh)
  Pmiss <- cumsum(Pmiss)
  runningES <- Phit - Pmiss
  max.ES <- max(runningES, na.rm = TRUE)
  min.ES <- min(runningES, na.rm = TRUE)
  if (abs(max.ES) > abs(min.ES)) {
    ES <- max.ES
  } else {
    ES <- min.ES
  }
  df <- data.frame(
    x = seq_along(runningES),
    runningScore = runningES,
    position = as.integer(hits),
    gene = names(geneList)
  )
  return(df)
}
