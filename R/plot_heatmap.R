compact_heatmap_feature_annotation <- function(values, annotation_name) {
  if (is.numeric(values) || is.logical(values)) {
    return(values)
  }
  value_chr <- as.character(values)
  has_value <- !is.na(value_chr) & nzchar(value_chr)
  if (!any(has_value)) {
    return(values)
  }

  needs_compact <- any(grepl(";", value_chr[has_value], fixed = TRUE)) ||
    any(nchar(value_chr[has_value], type = "width") > 80)
  if (!isTRUE(needs_compact)) {
    return(values)
  }

  out <- rep(NA_character_, length(value_chr))
  names(out) <- names(values)
  out[has_value] <- annotation_name
  out
}

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
  terms_stat = "none",
  terms_stat_digits = 2,
  terms_stat_test = TRUE,
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
  cores = 1,
  ...
) {
  res <- NULL
  lgd <- list()
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
      cores = cores,
      ...
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
            use_graphic_terms <- !heatmap_enrichment_terms_stat_none(terms_stat)
            terms_list <- lapply(subdf_list, function(df) {
              heatmap_enrichment_terms_data(
                df = df,
                metric = metric,
                topTerm = topTerm,
                show_termid = show_termid,
                terms_fontsize = terms_fontsize,
                terms_stat = if (isTRUE(use_graphic_terms)) terms_stat else "none",
                terms_stat_digits = terms_stat_digits
              )
            })
            names(terms_list) <- unlist(lapply(nm, function(x) x[[2]]))
            lgd[[paste0(enrich, "_terms_score")]] <- heatmap_enrichment_score_legend_from_terms(
              terms = terms_list,
              title = paste0(enrich, " terms\n-log10(", metric, ")")
            )
            terms_annotation <- if (isTRUE(use_graphic_terms)) {
              terms_list <- heatmap_enrichment_rescale_terms_stats(terms_list)
              heatmap_enrichment_terms_annotation(
                align_to = geneID_groups,
                terms = terms_list,
                width = terms_width,
                show_stat_text = terms_stat_test,
                which = "row"
              )
            } else {
              textbox_terms <- lapply(terms_list, heatmap_enrichment_terms_textbox_df)
              ComplexHeatmap::anno_textbox(
                align_to = geneID_groups,
                text = textbox_terms,
                max_width = terms_width,
                word_wrap = TRUE,
                add_new_line = TRUE,
                background_gp = grid::gpar(fill = "grey98", col = "black"),
                round_corners = TRUE,
                which = "row"
              )
            }
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
                "terms" = terms_annotation,
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
  list(ha_right = ha_right, res = res, lgd = lgd)
}

heatmap_enrichment_terms_data <- function(
  df,
  metric,
  topTerm,
  show_termid,
  terms_fontsize,
  terms_stat = "none",
  terms_stat_digits = 2
) {
  df <- utils::head(df, topTerm)
  terms <- if (isTRUE(show_termid)) {
    paste(df[["ID"]], df[["Description"]])
  } else {
    capitalize(df[["Description"]])
  }
  score <- -log10(df[[metric]])
  term_col <- heatmap_enrichment_score_colors(score)
  stat <- heatmap_enrichment_term_stat(
    df = df,
    metric = metric,
    score = score,
    terms_stat = terms_stat,
    digits = terms_stat_digits
  )
  if (is.null(stat)) {
    return(data.frame(
      term = terms,
      col = term_col,
      term_score = score,
      fontsize = rep(terms_fontsize, length(terms)),
      stat_label = NA_character_,
      stat_display = NA_character_,
      stat_strength = NA_real_,
      stat_scaled = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  strength <- stat[["strength"]]
  finite_strength <- strength[is.finite(strength)]
  if (length(finite_strength) == 0 || diff(range(finite_strength)) == 0) {
    scaled <- rep(1, length(strength))
  } else {
    scaled <- (strength - min(finite_strength)) / diff(range(finite_strength))
    scaled[!is.finite(scaled)] <- 0
  }
  data.frame(
    term = terms,
    col = term_col,
    term_score = score,
    fontsize = rep(terms_fontsize, length(terms)),
    stat_label = stat[["label"]],
    stat_display = stat[["display"]],
    stat_strength = strength,
    stat_scaled = scaled,
    stringsAsFactors = FALSE
  )
}

heatmap_enrichment_terms_df <- function(
  df,
  metric,
  topTerm,
  show_termid,
  terms_fontsize,
  terms_stat = "none",
  terms_stat_digits = 2
) {
  data <- heatmap_enrichment_terms_data(
    df = df,
    metric = metric,
    topTerm = topTerm,
    show_termid = show_termid,
    terms_fontsize = terms_fontsize,
    terms_stat = terms_stat,
    terms_stat_digits = terms_stat_digits
  )
  keyword <- data[["term"]]
  has_stat <- !is.na(data[["stat_label"]])
  keyword[has_stat] <- paste0(
    keyword[has_stat],
    " (",
    data[["stat_label"]][has_stat],
    "=",
    data[["stat_display"]][has_stat],
    ")"
  )
  df_out <- data.frame(keyword = keyword, stringsAsFactors = FALSE)
  df_out[["col"]] <- data[["col"]]
  df_out[["fontsize"]] <- data[["fontsize"]]
  df_out
}

heatmap_enrichment_terms_textbox_df <- function(data) {
  df_out <- data.frame(keyword = data[["term"]], stringsAsFactors = FALSE)
  df_out[["col"]] <- data[["col"]]
  df_out[["fontsize"]] <- data[["fontsize"]]
  df_out
}

heatmap_enrichment_rescale_terms_stats <- function(terms) {
  strength <- unlist(lapply(terms, function(data) {
    data[["stat_strength"]]
  }), use.names = FALSE)
  finite_strength <- strength[is.finite(strength)]
  if (length(finite_strength) == 0 || max(finite_strength) <= 0) {
    return(lapply(terms, function(data) {
      data[["stat_scaled"]] <- ifelse(
        is.finite(data[["stat_strength"]]),
        1,
        0
      )
      data
    }))
  }
  max_strength <- max(finite_strength)
  lapply(terms, function(data) {
    scaled <- data[["stat_strength"]] / max_strength
    scaled[!is.finite(scaled)] <- 0
    data[["stat_scaled"]] <- scaled
    data
  })
}

heatmap_enrichment_terms_annotation <- function(
  align_to,
  terms,
  width,
  show_stat_text = TRUE,
  which = "row"
) {
  force(width)
  force(show_stat_text)
  force(which)
  align_to <- heatmap_enrichment_terms_align_to(
    align_to = align_to,
    terms = terms
  )
  terms <- terms[names(align_to)]
  if (length(terms) == 0) {
    return(ComplexHeatmap::anno_empty(width = width, border = FALSE, which = which))
  }
  stat_width_npc <- 0.995
  size <- heatmap_enrichment_terms_size(
    terms = terms,
    max_width = width
  )
  width <- heatmap_enrichment_terms_width(
    terms = terms,
    max_width = width,
    show_stat_text = show_stat_text
  )
  ComplexHeatmap::anno_link(
    align_to = align_to,
    which = which,
    side = "right",
    size = size,
    gap = grid::unit(2, "mm"),
    width = width,
    link_gp = grid::gpar(fill = "grey98", col = "black", lwd = 0.8),
    internal_line = FALSE,
    panel_fun = function(index, nm = NULL) {
      if (is.null(nm) || !nm %in% names(terms)) {
        return(invisible(NULL))
      }
      data <- terms[[nm]]
      if (is.null(data) || nrow(data) == 0) {
        return(invisible(NULL))
      }
      heatmap_enrichment_draw_terms_box(
        data = data,
        stat_width_npc = stat_width_npc,
        max_width = width,
        show_stat_text = show_stat_text
      )
      invisible(NULL)
    }
  )
}

heatmap_enrichment_terms_align_to <- function(align_to, terms) {
  term_names <- names(terms)
  if (is.null(term_names)) {
    term_names <- character(0)
  }
  if (is.list(align_to)) {
    if (is.null(names(align_to))) {
      names(align_to) <- term_names[seq_along(align_to)]
    }
    return(align_to[intersect(names(align_to), term_names)])
  }
  align_to_names <- as.character(align_to)
  split(seq_along(align_to_names), align_to_names)[term_names]
}

heatmap_enrichment_terms_size <- function(terms, max_width) {
  sizes <- lapply(terms, function(data) {
    fontsize <- data[["fontsize"]]
    n_lines <- heatmap_enrichment_terms_n_lines(data, max_width = max_width)
    grid::unit(
      sum(n_lines * fontsize * 1.55) + 10,
      "pt"
    )
  })
  do.call(grid::unit.c, sizes)
}

heatmap_enrichment_terms_width <- function(terms, max_width, show_stat_text = TRUE) {
  widths <- lapply(terms, function(data) {
    label_widths <- lapply(seq_len(nrow(data)), function(i) {
      grid::grobWidth(grid::textGrob(
        label = data[["term"]][[i]],
        gp = grid::gpar(fontsize = data[["fontsize"]][[i]])
      ))
    })
    label_widths <- do.call(grid::unit.c, label_widths)
    text_width <- min(max(label_widths) + grid::unit(8, "mm"), max_width)
    stat_width <- if (isTRUE(show_stat_text)) {
      heatmap_enrichment_terms_stat_width(data)
    } else {
      grid::unit(0, "mm")
    }
    text_width + stat_width
  })
  max(do.call(grid::unit.c, widths))
}

heatmap_enrichment_terms_stat_width <- function(data) {
  has_stat <- !is.na(data[["stat_label"]])
  if (!any(has_stat)) {
    return(grid::unit(0, "mm"))
  }
  stat_widths <- lapply(which(has_stat), function(i) {
    grid::grobWidth(grid::textGrob(
      label = data[["stat_display"]][[i]],
      gp = grid::gpar(fontsize = data[["fontsize"]][[i]])
    ))
  })
  stat_widths <- do.call(grid::unit.c, stat_widths)
  max(stat_widths) + grid::unit(7, "mm")
}

heatmap_enrichment_terms_n_lines <- function(data, max_width) {
  max_width_in <- grid::convertWidth(max_width, "in", valueOnly = TRUE)
  vapply(seq_len(nrow(data)), function(i) {
    term <- data[["term"]][[i]]
    wrap_width <- heatmap_enrichment_terms_wrap_chars(
      width_in = max_width_in,
      fontsize = data[["fontsize"]][[i]]
    )
    max(length(strwrap(term, width = wrap_width)), 1L)
  }, integer(1))
}

heatmap_enrichment_terms_wrap_chars <- function(width_in, fontsize) {
  pmax(20L, floor(width_in * 72 * 0.95 / (fontsize * 0.35)))
}

heatmap_enrichment_draw_terms_box <- function(
  data,
  stat_width_npc = 0.3,
  max_width = NULL,
  show_stat_text = TRUE
) {
  grid::grid.rect(
    x = grid::unit(0.5, "npc"),
    y = grid::unit(0.5, "npc"),
    width = grid::unit(1, "npc"),
    height = grid::unit(1, "npc"),
    gp = grid::gpar(fill = "grey98", col = NA)
  )
  grid::grid.lines(
    x = c(0, 1, 1, 0),
    y = c(0, 0, 1, 1),
    default.units = "npc",
    gp = grid::gpar(col = "black", lwd = 0.8)
  )
  box_width_in <- grid::convertWidth(grid::unit(1, "npc"), "in", valueOnly = TRUE)
  fontsize <- data[["fontsize"]]
  stat_width_in <- if (isTRUE(show_stat_text)) {
    grid::convertWidth(
      heatmap_enrichment_terms_stat_width(data),
      "in",
      valueOnly = TRUE
    )
  } else {
    0
  }
  value_width_npc <- if (box_width_in > 0) {
    min(0.35, stat_width_in / box_width_in)
  } else {
    0
  }
  term_width <- max(0.4, 0.98 - value_width_npc)
  bar_width_npc <- max(0.2, 0.97 - value_width_npc)
  wrap_width <- pmax(16L, floor(box_width_in * 72 * term_width / (fontsize * 0.35)))
  term_lines <- lapply(seq_len(nrow(data)), function(i) {
    term <- data[["term"]][[i]]
    out <- strwrap(term, width = wrap_width[[i]])
    if (length(out) == 0) {
      out <- ""
    }
    out
  })
  has_graphic <- !is.na(data[["stat_label"]])
  block_units <- vapply(term_lines, length, integer(1))
  gap_units <- rep(0, length(block_units))
  total_units <- sum(block_units) + sum(gap_units[-length(gap_units)])
  if (!is.finite(total_units) || total_units <= 0) {
    total_units <- 1
  }
  bottom_reserved <- 0.04
  line_height <- (0.95 - bottom_reserved) / total_units
  y <- 0.95
  x_text <- 0.005

  for (i in seq_len(nrow(data))) {
    block_top <- y
    block_lines <- length(term_lines[[i]])
    block_bottom <- y - line_height * block_lines
    if (isTRUE(has_graphic[[i]])) {
      heatmap_enrichment_draw_term_stat(
        data = data[i, , drop = FALSE],
        y = block_top - line_height / 2,
        height_npc = max(line_height * 0.82, 0.016),
        stat_width_npc = bar_width_npc
      )
    }
    for (line in term_lines[[i]]) {
      y <- y - line_height / 2
      grid::grid.text(
        label = line,
        x = grid::unit(x_text, "npc"),
        y = grid::unit(y, "npc"),
        just = c("left", "center"),
        gp = grid::gpar(col = data[["col"]][[i]], fontsize = fontsize[[i]])
      )
      y <- y - line_height / 2
    }
    if (isTRUE(show_stat_text) && isTRUE(has_graphic[[i]])) {
      grid::grid.text(
        label = data[["stat_display"]][[i]],
        x = grid::unit(0.97, "npc"),
        y = grid::unit(block_top - line_height / 2, "npc"),
        just = c("right", "center"),
        gp = grid::gpar(col = "grey20", fontsize = fontsize[[i]])
      )
    }
    y <- y
  }
}

heatmap_enrichment_draw_term_stat <- function(
  data,
  y,
  height_npc = NULL,
  stat_width_npc = 0.3
) {
  scaled <- data[["stat_scaled"]][[1]]
  if (!is.finite(scaled)) {
    scaled <- 0
  }
  scaled <- max(0, min(1, scaled))
  col <- data[["col"]][[1]]
  x0 <- 0
  x1 <- min(0.96, x0 + stat_width_npc)
  grid::grid.rect(
    x = grid::unit(x0, "npc"),
    y = grid::unit(y, "npc"),
    width = grid::unit((x1 - x0) * scaled, "npc"),
    height = grid::unit(height_npc %||% 0.03, "npc"),
    just = c("left", "center"),
    gp = grid::gpar(fill = col, col = NA, alpha = 0.18)
  )
}

heatmap_enrichment_term_stat <- function(
  df,
  metric,
  score,
  terms_stat = "none",
  digits = 2
) {
  if (identical(terms_stat, FALSE) || identical(terms_stat, "none")) {
    return(NULL)
  }
  terms_stat <- as.character(terms_stat %||% "none")[[1]]
  if (identical(terms_stat, "none")) {
    return(NULL)
  }
  aliases <- c(
    score = "score",
    `-log10` = "score",
    `-log10(metric)` = "score",
    pvalue = "pvalue",
    p_value = "pvalue",
    padj = "p.adjust",
    p_adjust = "p.adjust"
  )
  stat_key <- if (terms_stat %in% names(aliases)) {
    aliases[[terms_stat]]
  } else {
    terms_stat
  }
  if (identical(stat_key, "score")) {
    value <- score
    strength <- score
    label <- paste0("-log10(", metric, ")")
    display_value <- strength
  } else {
    if (!stat_key %in% colnames(df)) {
      log_message(
        "{.arg terms_stat} must be {.val none}, {.val score}, or one enrichment column. ",
        "Available columns include: {.val {colnames(df)}}",
        message_type = "error"
      )
    }
    value <- df[[stat_key]]
    strength <- if (stat_key %in% c("pvalue", "p.adjust", "qvalue")) {
      -log10(suppressWarnings(as.numeric(value)))
    } else {
      suppressWarnings(as.numeric(value))
    }
    label <- stat_key
    display_value <- if (stat_key %in% c("pvalue", "p.adjust", "qvalue")) {
      strength
    } else {
      value
    }
  }
  display <- heatmap_enrichment_format_stat(display_value, digits = digits)
  list(label = label, display = display, strength = strength)
}

heatmap_enrichment_terms_stat_none <- function(terms_stat) {
  if (is.null(terms_stat) || identical(terms_stat, FALSE)) {
    return(TRUE)
  }
  identical(as.character(terms_stat)[[1]], "none")
}

heatmap_enrichment_format_stat <- function(value, digits = 2) {
  if (is.numeric(value)) {
    out <- vapply(value, function(x) {
      if (!is.finite(x)) {
        return(as.character(x))
      }
      if (x != 0 && (abs(x) < 0.01 || abs(x) >= 1e4)) {
        return(format(signif(x, digits = digits), scientific = TRUE, trim = TRUE))
      }
      formatted <- format(signif(x, digits = digits), scientific = FALSE, trim = TRUE)
      formatted <- sub("\\.?0+$", "", formatted)
      if (identical(formatted, "")) {
        formatted <- "0"
      }
      formatted
    }, character(1))
    return(out)
  }
  as.character(value)
}

heatmap_enrichment_score_colors <- function(score, palette = "Spectral", limits = NULL) {
  finite <- is.finite(score)
  out <- rep("black", length(score))
  if (!any(finite)) {
    return(out)
  }
  if (is.null(limits)) {
    limits <- range(score[finite])
  }
  limits <- as.numeric(limits)
  limits <- limits[is.finite(limits)]
  if (length(limits) < 2) {
    limits <- range(score[finite])
  } else {
    limits <- range(limits)
  }
  if (!is.finite(diff(limits)) || diff(limits) == 0) {
    limits <- limits[[1]] + c(-0.5, 0.5)
  }
  cols <- palette_colors(
    c(limits, score[finite]),
    type = "continuous",
    palette = palette,
    matched = TRUE
  )
  cols <- utils::tail(cols, sum(finite))
  out[finite] <- sapply(cols, function(x) blendcolors(c(x, x, "black")))
  out
}

heatmap_enrichment_score_legend_from_terms <- function(terms, title) {
  score <- unlist(lapply(terms, function(data) {
    data[["term_score"]]
  }), use.names = FALSE)
  keep <- is.finite(score)
  score <- score[keep]
  if (length(score) == 0) {
    return(NULL)
  }
  at <- sort(unique(round(score, digits = 1)))
  if (length(at) > 6) {
    at <- at[unique(round(seq(1, length(at), length.out = 6)))]
  }
  cols <- heatmap_enrichment_score_colors(at, limits = range(score))
  ComplexHeatmap::Legend(
    title = title,
    labels = heatmap_enrichment_format_stat(at, digits = 2),
    legend_gp = grid::gpar(fill = cols, col = cols),
    border = TRUE
  )
}
