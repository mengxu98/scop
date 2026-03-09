heatmap_enrichment <- function(
    geneID,
    geneID_groups,
    feature_split_palette = "simspec",
    feature_split_palcolor = NULL,
    ha_right = NULL,
    flip = FALSE,
    anno_terms = FALSE,
    anno_keys = FALSE,
    anno_features = FALSE,
    terms_width = grid::unit(4, "in"),
    terms_fontsize = 8,
    keys_width = grid::unit(2, "in"),
    keys_fontsize = c(6, 10),
    features_width = grid::unit(2, "in"),
    features_fontsize = c(6, 10),
    IDtype = "symbol",
    species = "Homo_sapiens",
    db_update = FALSE,
    db_combine = FALSE,
    db_version = "latest",
    convert_species = FALSE,
    Ensembl_version = NULL,
    mirror = NULL,
    db = "GO_BP",
    TERM2GENE = NULL,
    TERM2NAME = NULL,
    minGSSize = 10,
    maxGSSize = 500,
    GO_simplify = FALSE,
    GO_simplify_cutoff = "p.adjust < 0.05",
    simplify_method = "Wang",
    simplify_similarityCutoff = 0.7,
    pvalueCutoff = NULL,
    padjustCutoff = 0.05,
    topTerm = 5,
    show_termid = FALSE,
    topWord = 20,
    words_excluded = NULL,
    cores = 1) {
  res <- NULL
  words_excluded <- words_excluded %||% scop::words_excluded

  if (isTRUE(anno_keys) || isTRUE(anno_features) || isTRUE(anno_terms)) {
    if (isTRUE(flip)) {
      log_message(
        "{.arg anno_keys}, {.arg anno_features} and {.arg anno_terms} can only be used when {.arg flip = FALSE}",
        message_type = "error"
      )
    }
    if (all(is.na(geneID_groups))) {
      geneID_groups <- rep(1, length(geneID))
    }
    if (!is.factor(geneID_groups)) {
      geneID_groups <- factor(geneID_groups, levels = unique(geneID_groups))
    }
    fill_split <- palette_colors(
      levels(geneID_groups),
      type = "discrete",
      palette = feature_split_palette,
      palcolor = feature_split_palcolor
    )[levels(geneID_groups) %in% geneID_groups]
    res <- RunEnrichment(
      geneID = geneID,
      geneID_groups = geneID_groups,
      IDtype = IDtype,
      species = species,
      db_update = db_update,
      db_version = db_version,
      db_combine = db_combine,
      convert_species = convert_species,
      Ensembl_version = Ensembl_version,
      mirror = mirror,
      db = db,
      TERM2GENE = TERM2GENE,
      TERM2NAME = TERM2NAME,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      GO_simplify = GO_simplify,
      GO_simplify_cutoff = GO_simplify_cutoff,
      simplify_method = simplify_method,
      simplify_similarityCutoff = simplify_similarityCutoff,
      cores = cores
    )
    if (!is.null(TERM2GENE)) {
      db <- "custom"
    }
    if (isTRUE(db_combine)) {
      db <- "Combined"
    }
    go_dbs <- c("GO_BP", "GO_CC", "GO_MF")
    if (isTRUE(GO_simplify) && any(db %in% go_dbs)) {
      db[db %in% go_dbs] <- paste0(
        db[db %in% go_dbs],
        "_sim"
      )
    }
    if (nrow(res$enrichment) == 0) {
      log_message(
        "No enrichment result found",
        message_type = "warning"
      )
    } else {
      metric <- ifelse(is.null(padjustCutoff), "pvalue", "p.adjust")
      metric_value <- ifelse(
        is.null(padjustCutoff),
        pvalueCutoff,
        padjustCutoff
      )
      pvalueCutoff <- ifelse(is.null(pvalueCutoff), 1, pvalueCutoff)
      padjustCutoff <- ifelse(is.null(padjustCutoff), 1, padjustCutoff)

      df <- res$enrichment
      df <- df[df[["Database"]] %in% db, , drop = FALSE]
      df <- df[df[[metric]] < metric_value, , drop = FALSE]
      df <- df[order(df[[metric]]), , drop = FALSE]
      if (nrow(df) == 0) {
        log_message(
          "No term enriched using the threshold:\n",
          "pvalueCutoff = {.pkg {pvalueCutoff}}\n",
          "padjustCutoff = {.pkg {padjustCutoff}}",
          message_type = "warning"
        )
      } else {
        df_list <- split.data.frame(df, ~ Database + Groups)
        df_list <- df_list[lapply(df_list, nrow) > 0]

        for (enrich in db) {
          nm <- strsplit(names(df_list), "\\.")
          subdf_list <- df_list[
            unlist(lapply(nm, function(x) x[[1]])) %in% enrich
          ]
          if (length(subdf_list) == 0) {
            log_message(
              "No {.pkg {enrich}} term enriched using the threshold:\n",
              "pvalueCutoff = {.pkg {pvalueCutoff}}\n",
              "padjustCutoff = {.pkg {padjustCutoff}}",
              message_type = "warning"
            )
            next
          }
          nm <- strsplit(names(subdf_list), "\\.")

          ha_terms <- NULL
          if (isTRUE(anno_terms)) {
            terms_list <- lapply(subdf_list, function(df) {
              if (isTRUE(show_termid)) {
                terms <- paste(
                  utils::head(df$ID, topTerm),
                  utils::head(df$Description, topTerm)
                )
              } else {
                terms <- utils::head(df$Description, topTerm)
                terms <- capitalize(terms)
              }
              df_out <- data.frame(keyword = terms)
              df_out[["col"]] <- palette_colors(
                -log10(utils::head(df[, metric], topTerm)),
                type = "continuous",
                palette = "Spectral",
                matched = TRUE
              )
              df_out[["col"]] <- sapply(
                df_out[["col"]],
                function(x) blendcolors(c(x, "black"))
              )
              df_out[["fontsize"]] <- rep(terms_fontsize, nrow(df_out))
              return(df_out)
            })
            names(terms_list) <- unlist(lapply(nm, function(x) x[[2]]))
            if (length(intersect(geneID_groups, names(terms_list))) > 0) {
              ha_terms <- ComplexHeatmap::HeatmapAnnotation(
                "terms_empty" = ComplexHeatmap::anno_empty(
                  width = grid::unit(0.05, "in"),
                  border = FALSE,
                  which = "row"
                ),
                "terms_split" = ComplexHeatmap::anno_block(
                  gp = grid::gpar(fill = fill_split),
                  width = grid::unit(0.1, "in"),
                  which = "row"
                ),
                "terms" = ComplexHeatmap::anno_textbox(
                  align_to = geneID_groups,
                  text = terms_list,
                  max_width = terms_width,
                  word_wrap = TRUE,
                  add_new_line = TRUE,
                  background_gp = grid::gpar(fill = "grey98", col = "black"),
                  round_corners = TRUE,
                  which = "row"
                ),
                which = "row",
                gap = grid::unit(0, "points")
              )
              names(ha_terms) <- paste0(names(ha_terms), "_", enrich)
            }
          }

          ha_keys <- NULL
          if (isTRUE(anno_keys)) {
            check_r("simplifyEnrichment", verbose = FALSE)
            keys_list <- lapply(subdf_list, function(df) {
              if (all(df$Database %in% c("GO", "GO_BP", "GO_CC", "GO_MF"))) {
                df0 <- simplifyEnrichment::keyword_enrichment_from_GO(df[[
                  "ID"
                ]])
                if (nrow(df0) > 0) {
                  df <- df0 %>%
                    dplyr::reframe(
                      keyword = .data[["keyword"]],
                      score = -(log10(.data[["padj"]])),
                      count = .data[["n_term"]],
                      Database = df[["Database"]][1],
                      Groups = df[["Groups"]][1]
                    ) %>%
                    dplyr::filter(
                      !grepl(pattern = "\\[.*\\]", x = .data[["keyword"]])
                    ) %>%
                    dplyr::filter(nchar(.data[["keyword"]]) >= 1) %>%
                    dplyr::filter(
                      !tolower(.data[["keyword"]]) %in% tolower(words_excluded)
                    ) %>%
                    dplyr::distinct() %>%
                    dplyr::mutate(
                      angle = 90 *
                        sample(c(0, 1), dplyr::n(), replace = TRUE, prob = c(60, 40))
                    ) %>%
                    as.data.frame()
                  df <- df[
                    utils::head(order(df[["score"]], decreasing = TRUE), topWord), ,
                    drop = FALSE
                  ]
                } else {
                  df <- NULL
                }
              } else {
                df <- df %>%
                  dplyr::mutate(
                    keyword = strsplit(
                      tolower(as.character(.data[["Description"]])),
                      " "
                    )
                  ) %>%
                  unnest_fun(cols = "keyword") %>%
                  dplyr::group_by(.data[["keyword"]], Database, Groups) %>%
                  dplyr::reframe(
                    keyword = .data[["keyword"]],
                    score = sum(-(log10(.data[[metric]]))),
                    count = dplyr::n(),
                    Database = .data[["Database"]],
                    Groups = .data[["Groups"]],
                    .groups = "keep"
                  ) %>%
                  dplyr::filter(
                    !grepl(pattern = "\\[.*\\]", x = .data[["keyword"]])
                  ) %>%
                  dplyr::filter(nchar(.data[["keyword"]]) >= 1) %>%
                  dplyr::filter(
                    !tolower(.data[["keyword"]]) %in% tolower(words_excluded)
                  ) %>%
                  dplyr::distinct() %>%
                  dplyr::mutate(
                    angle = 90 *
                      sample(c(0, 1), dplyr::n(), replace = TRUE, prob = c(60, 40))
                  ) %>%
                  as.data.frame()
                df <- df[
                  utils::head(order(df[["score"]], decreasing = TRUE), topWord), ,
                  drop = FALSE
                ]
              }
              if (isTRUE(nrow(df) > 0)) {
                df[["col"]] <- palette_colors(
                  df[, "score"],
                  type = "continuous",
                  palette = "Spectral",
                  matched = TRUE
                )
                df[["col"]] <- sapply(
                  df[["col"]],
                  function(x) blendcolors(c(x, "black"))
                )
                df[["fontsize"]] <- scales::rescale(
                  df[, "count"],
                  to = keys_fontsize
                )
                return(df)
              } else {
                return(NULL)
              }
            })
            names(keys_list) <- unlist(lapply(nm, function(x) x[[2]]))
            keys_list <- keys_list[lapply(keys_list, length) > 0]
            if (length(intersect(geneID_groups, names(keys_list))) > 0) {
              ha_keys <- ComplexHeatmap::HeatmapAnnotation(
                "keys_empty" = ComplexHeatmap::anno_empty(
                  width = grid::unit(0.05, "in"),
                  border = FALSE,
                  which = "row"
                ),
                "keys_split" = ComplexHeatmap::anno_block(
                  gp = grid::gpar(fill = fill_split),
                  width = grid::unit(0.1, "in"),
                  which = "row"
                ),
                "keys" = ComplexHeatmap::anno_textbox(
                  align_to = geneID_groups,
                  text = keys_list,
                  max_width = keys_width,
                  background_gp = grid::gpar(fill = "grey98", col = "black"),
                  round_corners = TRUE,
                  which = "row"
                ),
                which = "row",
                gap = grid::unit(0, "points")
              )
              names(ha_keys) <- paste0(names(ha_keys), "_", enrich)
            }
          }

          ha_features <- NULL
          if (isTRUE(anno_features)) {
            features_list <- lapply(subdf_list, function(df) {
              df <- df %>%
                dplyr::mutate(
                  keyword = strsplit(as.character(.data[["geneID"]]), "/")
                ) %>%
                unnest_fun(cols = "keyword") %>%
                dplyr::group_by(.data[["keyword"]], Database, Groups) %>%
                dplyr::reframe(
                  keyword = .data[["keyword"]],
                  score = sum(-(log10(.data[[metric]]))),
                  count = dplyr::n(),
                  Database = .data[["Database"]],
                  Groups = .data[["Groups"]],
                  .groups = "keep"
                ) %>%
                dplyr::distinct() %>%
                dplyr::mutate(
                  angle = 90 *
                    sample(
                      c(0, 1),
                      dplyr::n(),
                      replace = TRUE, prob = c(60, 40)
                    )
                ) %>%
                as.data.frame()
              df <- df[
                utils::head(order(df[["score"]], decreasing = TRUE), topWord), ,
                drop = FALSE
              ]
              df[["col"]] <- palette_colors(
                df[, "score"],
                type = "continuous",
                palette = "Spectral",
                matched = TRUE
              )
              df[["col"]] <- sapply(
                df[["col"]],
                function(x) blendcolors(c(x, "black"))
              )
              df[["fontsize"]] <- scales::rescale(
                df[, "count"],
                to = features_fontsize
              )
              return(df)
            })
            names(features_list) <- unlist(lapply(nm, function(x) x[[2]]))
            if (length(intersect(geneID_groups, names(features_list))) > 0) {
              ha_features <- ComplexHeatmap::HeatmapAnnotation(
                "features_empty" = ComplexHeatmap::anno_empty(
                  width = grid::unit(0.05, "in"),
                  border = FALSE,
                  which = "row"
                ),
                "features_split" = ComplexHeatmap::anno_block(
                  gp = grid::gpar(fill = fill_split),
                  width = grid::unit(0.1, "in"),
                  which = "row"
                ),
                "features" = ComplexHeatmap::anno_textbox(
                  align_to = geneID_groups,
                  text = features_list,
                  max_width = features_width,
                  background_gp = grid::gpar(fill = "grey98", col = "black"),
                  round_corners = TRUE,
                  which = "row"
                ),
                which = "row",
                gap = grid::unit(0, "points")
              )
              names(ha_features) <- paste0(names(ha_features), "_", enrich)
            }
          }

          ha_enrichment <- list(ha_terms, ha_keys, ha_features)
          ha_enrichment <- ha_enrichment[sapply(ha_enrichment, length) > 0]
          ha_enrichment <- do.call(c, ha_enrichment)

          if (is.null(ha_right)) {
            ha_right <- ha_enrichment
          } else {
            ha_right <- c(ha_right, ha_enrichment)
          }
        }
      }
    }
  }
  list(ha_right = ha_right, res = res)
}
