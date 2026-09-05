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

heatmap_wrap_row_labels <- function(labels, width = NULL) {
  if (is.null(width)) {
    return(labels)
  }
  if (
    length(width) != 1L ||
      !is.numeric(width) ||
      is.na(width) ||
      !is.finite(width) ||
      width < 1
  ) {
    cli::cli_abort(
      "{.arg row_names_wrap} must be a single positive number or {.code NULL}."
    )
  }
  width <- as.integer(floor(width))
  wrapped <- vapply(as.character(labels), function(label) {
    if (is.na(label)) {
      return(NA_character_)
    }
    display_label <- gsub("_", " ", label, fixed = TRUE)
    paragraphs <- strsplit(display_label, "\n", fixed = TRUE)[[1]]
    lines <- unlist(lapply(paragraphs, function(paragraph) {
      out <- strwrap(paragraph, width = width)
      if (length(out) == 0L) "" else out
    }), use.names = FALSE)
    paste(lines, collapse = "\n")
  }, character(1))
  names(wrapped) <- names(labels)
  wrapped
}

heatmap_row_labels_max_width <- function(labels, gp = NULL) {
  gp <- gp %||% grid::gpar(fontsize = 10)
  lines <- unlist(strsplit(as.character(labels), "\n", fixed = TRUE), use.names = FALSE)
  lines <- lines[!is.na(lines)]
  if (length(lines) == 0L) {
    return(grid::unit(0, "mm"))
  }
  ComplexHeatmap::max_text_width(lines, gp = gp) + grid::unit(2, "mm")
}

heatmap_discrete_legend_gp <- function(
  fill,
  border,
  border_color,
  border_size
) {
  do.call(
    grid::gpar,
    c(
      list(fill = fill),
      heatmap_border_gp(border, border_color, border_size)
    )
  )
}

build_heatmap_feature_annotation <- function(
  annotation_name,
  values,
  palette,
  palcolor = NULL,
  flip = FALSE,
  border = TRUE,
  annotation_border = border,
  border_color = "black",
  border_size = 1,
  discrete_legend_border = heatmap_legend_border(border, border_color),
  discrete_legend_gp_border = TRUE,
  params = list()
) {
  values <- compact_heatmap_feature_annotation(values, annotation_name)
  which <- ifelse(flip, "column", "row")
  annotation <- list()

  if (!is.numeric(values)) {
    if (is.logical(values)) {
      values <- factor(values, levels = c(TRUE, FALSE))
    } else if (!is.factor(values)) {
      values <- factor(values, levels = unique(values))
    }
    annotation[[annotation_name]] <- ComplexHeatmap::anno_simple(
      x = as.character(values),
      col = palette_colors(values, palette = palette, palcolor = palcolor),
      which = which,
      na_col = "transparent",
      border = annotation_border
    )
    levels <- levels(values)
    levels <- levels[!is.na(levels) & nzchar(levels)]
    legend <- if (length(levels) > 0) {
      legend_gp <- list(
        fill = palette_colors(levels, palette = palette, palcolor = palcolor)
      )
      if (isTRUE(discrete_legend_gp_border)) {
        legend_gp <- c(
          legend_gp,
          col = heatmap_border_color(border, border_color),
          lwd = border_size
        )
      }
      ComplexHeatmap::Legend(
        title = annotation_name,
        labels = levels,
        legend_gp = do.call(grid::gpar, legend_gp),
        border = discrete_legend_border
      )
    } else {
      NULL
    }
  } else {
    col_fun <- circlize::colorRamp2(
      breaks = seq(
        min(values, na.rm = TRUE),
        max(values, na.rm = TRUE),
        length = 100
      ),
      colors = palette_colors(palette = palette, palcolor = palcolor)
    )
    annotation[[annotation_name]] <- ComplexHeatmap::anno_simple(
      x = values,
      col = col_fun,
      which = which,
      na_col = "transparent",
      border = annotation_border
    )
    legend <- ComplexHeatmap::Legend(
      title = annotation_name,
      col_fun = col_fun,
      legend_gp = heatmap_border_gp(border, border_color, border_size),
      border = heatmap_legend_border(border, border_color)
    )
  }

  annotation <- build_heatmap_annotation(
    annotations = annotation,
    which = which,
    show_annotation_name = TRUE,
    annotation_name_side = ifelse(flip, "left", "top"),
    border = annotation_border,
    params = params
  )
  list(annotation = annotation, legend = legend)
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
  terms_stat_width = grid::unit(1.35, "in"),
  terms_fontsize = 8,
  terms_stat = "none",
  terms_stat_digits = 2,
  terms_stat_label = "value",
  terms_stat_axis = FALSE,
  terms_stat_background_palcolor = NULL,
  terms_stat_border = NULL,
  terms_stat_border_palcolor = NULL,
  terms_stat_border_size = NULL,
  terms_stat_label_palcolor = NULL,
  terms_group_background = FALSE,
  terms_background_palcolor = "grey98",
  terms_background_alpha = 1,
  terms_border = TRUE,
  terms_border_palcolor = "black",
  terms_border_size = 0.8,
  terms_text_palcolor = NULL,
  terms_bar_palcolor = NULL,
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
                terms_stat = if (isTRUE(use_graphic_terms)) terms_stat else "none"
              )
            })
            names(terms_list) <- unlist(lapply(nm, function(x) x[[2]]))
            uses_term_score_colors <-
              !isTRUE(use_graphic_terms) ||
                is.null(terms_text_palcolor)
            if (isTRUE(uses_term_score_colors)) {
              lgd[[paste0(enrich, "_terms_score")]] <-
                heatmap_enrichment_score_legend_from_terms(
                  terms = terms_list,
                  title = paste0(enrich, " terms\n-log10(", metric, ")")
                )
            }
            terms_annotation <- if (isTRUE(use_graphic_terms)) {
              terms_list <- heatmap_enrichment_rescale_terms_stats(terms_list)
              heatmap_enrichment_terms_annotation(
                align_to = geneID_groups,
                terms = terms_list,
                width = terms_width,
                stat_width = terms_stat_width,
                axis_digits = terms_stat_digits,
                stat_label = terms_stat_label,
                group_colors = fill_split,
                show_axis = terms_stat_axis,
                stat_background = terms_stat_background_palcolor,
                stat_border = terms_stat_border,
                stat_border_color = terms_stat_border_palcolor,
                stat_border_size = terms_stat_border_size,
                stat_label_color = terms_stat_label_palcolor,
                group_background = terms_group_background,
                background_color = terms_background_palcolor,
                background_alpha = terms_background_alpha,
                border = terms_border,
                border_color = terms_border_palcolor,
                border_size = terms_border_size,
                text_color = terms_text_palcolor,
                bar_color = terms_bar_palcolor,
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
              ha_terms_args <- list(
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
                "terms" = if (is.list(terms_annotation)) {
                  terms_annotation[["terms"]]
                } else {
                  terms_annotation
                },
                which = "row",
                gap = grid::unit(0, "points")
              )
              if (is.list(terms_annotation)) {
                ha_terms_args[["terms_stat"]] <- terms_annotation[["stat"]]
                ha_terms_args[["gap"]] <- grid::unit(c(0, 0, 1), "mm")
              }
              ha_terms <- do.call(
                ComplexHeatmap::HeatmapAnnotation,
                ha_terms_args
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
  terms_stat = "none"
) {
  df <- utils::head(df, topTerm)
  terms <- if (isTRUE(show_termid)) {
    paste(df[["ID"]], df[["Description"]])
  } else {
    capitalize(df[["Description"]])
  }
  metric_value <- suppressWarnings(as.numeric(df[[metric]]))
  score <- -log10(pmax(metric_value, .Machine$double.xmin))
  term_col <- heatmap_enrichment_score_colors(score)
  stat <- heatmap_enrichment_term_stat(
    df = df,
    metric = metric,
    score = score,
    terms_stat = terms_stat
  )
  if (is.null(stat)) {
    return(data.frame(
      term = terms,
      col = term_col,
      term_score = score,
      fontsize = rep(terms_fontsize, length(terms)),
      stat_label = NA_character_,
      stat_signif = NA_character_,
      stat_axis_label = NA_character_,
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
    stat_signif = heatmap_enrichment_p_significance(stat[["pvalue"]]),
    stat_axis_label = stat[["axis_label"]],
    stat_strength = strength,
    stat_scaled = scaled,
    stringsAsFactors = FALSE
  )
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
      data[["stat_scaled"]] <- rep(0, nrow(data))
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
  stat_width = grid::unit(1.35, "in"),
  axis_digits = 2,
  stat_label = "value",
  group_colors = NULL,
  show_axis = FALSE,
  stat_background = NULL,
  stat_border = NULL,
  stat_border_color = NULL,
  stat_border_size = NULL,
  stat_label_color = NULL,
  group_background = FALSE,
  background_color = "grey98",
  background_alpha = 1,
  border = TRUE,
  border_color = "black",
  border_size = 0.8,
  text_color = NULL,
  bar_color = NULL,
  which = "row"
) {
  force(width)
  force(stat_width)
  force(axis_digits)
  force(which)
  stat_label <- match.arg(
    as.character(stat_label %||% "none")[[1]],
    c("none", "value", "significance", "both")
  )
  align_to <- heatmap_enrichment_terms_align_to(
    align_to = align_to,
    terms = terms
  )
  terms <- terms[names(align_to)]
  if (length(terms) == 0) {
    return(ComplexHeatmap::anno_empty(width = width, border = FALSE, which = which))
  }
  term_groups <- names(terms)
  last_group <- names(which.max(vapply(
    align_to,
    function(index) max(index, na.rm = TRUE),
    numeric(1)
  )))
  panel_colors <- if (isTRUE(group_background)) {
    heatmap_enrichment_group_colors(
      groups = term_groups,
      colors = group_colors,
      alpha = background_alpha
    )
  } else {
    colors <- rep(background_color, length(term_groups))
    colors <- grDevices::adjustcolor(colors, alpha.f = background_alpha)
    stats::setNames(as.list(colors), term_groups)
  }
  stat_strength <- unlist(lapply(terms, function(data) {
    data[["stat_strength"]]
  }), use.names = FALSE)
  stat_strength <- stat_strength[is.finite(stat_strength)]
  stat_max <- if (length(stat_strength) == 0) 1 else max(stat_strength)
  stat_axis_label <- unique(unlist(lapply(terms, function(data) {
    data[["stat_axis_label"]]
  }), use.names = FALSE))
  stat_axis_label <- stat_axis_label[
    !is.na(stat_axis_label) & nzchar(stat_axis_label)
  ]
  stat_axis_label <- if (length(stat_axis_label) == 0) "" else stat_axis_label[[1]]
  size <- heatmap_enrichment_terms_size(
    terms = terms,
    max_width = width,
    axis_group = if (isTRUE(show_axis)) last_group else NULL
  )
  terms_width <- heatmap_enrichment_terms_width(
    terms = terms,
    max_width = width
  )
  terms_annotation <- ComplexHeatmap::anno_link(
    align_to = align_to,
    which = which,
    side = "right",
    size = size,
    gap = grid::unit(2, "mm"),
    width = terms_width,
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
        max_width = terms_width,
        background = panel_colors[[nm]],
        border_color = if (isTRUE(border)) border_color else NA,
        border_size = border_size,
        term_color = text_color,
        reserve_axis = isTRUE(show_axis) && identical(nm, last_group)
      )
      invisible(NULL)
    }
  )
  stat_annotation <- ComplexHeatmap::anno_link(
    align_to = align_to,
    which = which,
    side = "right",
    size = size,
    gap = grid::unit(2, "mm"),
    link_width = grid::unit(0, "mm"),
    width = stat_width,
    link_gp = grid::gpar(fill = NA, col = NA),
    internal_line = FALSE,
    panel_fun = function(index, nm = NULL) {
      if (is.null(nm) || !nm %in% names(terms)) {
        return(invisible(NULL))
      }
      data <- terms[[nm]]
      if (is.null(data) || nrow(data) == 0) {
        return(invisible(NULL))
      }
      heatmap_enrichment_draw_stat_box(
        data = data,
        max_width = width,
        background = stat_background %||% panel_colors[[nm]],
        border_color = if (isTRUE(stat_border %||% border)) {
          stat_border_color %||% border_color
        } else {
          NA
        },
        border_size = stat_border_size %||% border_size,
        bar_color = bar_color %||% text_color,
        label_mode = stat_label,
        label_color = stat_label_color,
        show_axis = isTRUE(show_axis) && identical(nm, last_group),
        axis_max = stat_max,
        axis_label = stat_axis_label,
        axis_digits = axis_digits
      )
      invisible(NULL)
    }
  )
  list(terms = terms_annotation, stat = stat_annotation)
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

heatmap_enrichment_group_colors <- function(groups, colors = NULL, alpha = 0.16) {
  groups <- as.character(groups)
  if (is.null(colors) || length(colors) == 0) {
    colors <- palette_colors(groups, type = "discrete", palette = "simspec")
  }
  if (!is.null(names(colors))) {
    colors <- colors[groups]
  } else {
    colors <- rep(colors, length.out = length(groups))
  }
  colors[is.na(colors) | !nzchar(colors)] <- "grey85"
  colors <- grDevices::adjustcolor(colors, alpha.f = alpha)
  stats::setNames(as.list(colors), groups)
}

heatmap_enrichment_terms_size <- function(terms, max_width, axis_group = NULL) {
  sizes <- lapply(names(terms), function(group) {
    data <- terms[[group]]
    fontsize <- data[["fontsize"]]
    n_lines <- heatmap_enrichment_terms_n_lines(data, max_width = max_width)
    grid::unit(
      sum(n_lines * fontsize * 1.55) +
        10 +
        ifelse(group %in% axis_group, 24, 0),
      "pt"
    )
  })
  do.call(grid::unit.c, sizes)
}

heatmap_enrichment_terms_width <- function(terms, max_width) {
  widths <- lapply(terms, function(data) {
    label_widths <- lapply(seq_len(nrow(data)), function(i) {
      grid::grobWidth(grid::textGrob(
        label = data[["term"]][[i]],
        gp = grid::gpar(fontsize = data[["fontsize"]][[i]])
      ))
    })
    label_widths <- do.call(grid::unit.c, label_widths)
    min(max(label_widths) + grid::unit(8, "mm"), max_width)
  })
  max(do.call(grid::unit.c, widths))
}

heatmap_enrichment_terms_n_lines <- function(data, max_width) {
  max_width_in <- grid::convertWidth(max_width, "in", valueOnly = TRUE)
  vapply(seq_len(nrow(data)), function(i) {
    term <- data[["term"]][[i]]
    wrap_width <- {
      .inline0 <- max_width_in
      .inline1 <- data[["fontsize"]][[i]]
      pmax(20L, floor(.inline0 * 72 * 0.95 / (.inline1 * 0.35)))
    }
    max(length(strwrap(term, width = wrap_width)), 1L)
  }, integer(1))
}


heatmap_enrichment_draw_terms_box <- function(
  data,
  max_width = NULL,
  background = "grey98",
  border_color = "black",
  border_size = 0.8,
  term_color = NULL,
  reserve_axis = FALSE
) {
  grid::grid.rect(
    x = grid::unit(0.5, "npc"),
    y = grid::unit(0.5, "npc"),
    width = grid::unit(1, "npc"),
    height = grid::unit(1, "npc"),
    gp = grid::gpar(fill = background, col = border_color, lwd = border_size)
  )
  layout <- heatmap_enrichment_terms_layout(
    data = data,
    max_width = max_width,
    reserve_axis = reserve_axis
  )
  fontsize <- data[["fontsize"]]
  x_text <- 0.005

  for (i in seq_len(nrow(data))) {
    y <- layout[["block_top"]][[i]]
    for (line in layout[["term_lines"]][[i]]) {
      y <- y - layout[["line_height"]] / 2
      grid::grid.text(
        label = line,
        x = grid::unit(x_text, "npc"),
        y = grid::unit(y, "npc"),
        just = c("left", "center"),
        gp = grid::gpar(
          col = term_color %||% data[["col"]][[i]],
          fontsize = fontsize[[i]]
        )
      )
      y <- y - layout[["line_height"]] / 2
    }
  }
}

heatmap_enrichment_terms_layout <- function(
  data,
  max_width,
  reserve_axis = FALSE
) {
  width_in <- grid::convertWidth(max_width, "in", valueOnly = TRUE)
  fontsize <- data[["fontsize"]]
  wrap_width <- pmax(16L, floor(width_in * 72 * 0.95 / (fontsize * 0.35)))
  term_lines <- lapply(seq_len(nrow(data)), function(i) {
    out <- strwrap(data[["term"]][[i]], width = wrap_width[[i]])
    if (length(out) == 0) "" else out
  })
  block_units <- vapply(term_lines, length, integer(1))
  total_units <- sum(block_units)
  if (!is.finite(total_units) || total_units <= 0) {
    total_units <- 1
  }
  bottom_reserved <- if (isTRUE(reserve_axis)) 0.42 else 0.04
  line_height <- (0.95 - bottom_reserved) / total_units
  block_top <- numeric(length(block_units))
  block_center <- numeric(length(block_units))
  y <- 0.95
  for (i in seq_along(block_units)) {
    block_top[[i]] <- y
    block_height <- line_height * block_units[[i]]
    block_center[[i]] <- y - block_height / 2
    y <- y - block_height
  }
  list(
    term_lines = term_lines,
    line_height = line_height,
    block_top = block_top,
    block_center = block_center
  )
}

heatmap_enrichment_draw_stat_box <- function(
  data,
  max_width,
  background = "white",
  border_color = "grey60",
  border_size = 0.8,
  bar_color = NULL,
  label_mode = "value",
  label_color = NULL,
  show_axis = FALSE,
  axis_max = NULL,
  axis_label = "",
  axis_digits = 2
) {
  # The statistic chart is an independent annotation panel. Its shared scale
  # occupies an unfilled strip outside the bordered plotting region.
  panel_bottom <- if (isTRUE(show_axis)) 0.34 else 0
  grid::grid.rect(
    x = grid::unit(0.5, "npc"),
    y = grid::unit(panel_bottom + (1 - panel_bottom) / 2, "npc"),
    width = grid::unit(1, "npc"),
    height = grid::unit(1 - panel_bottom, "npc"),
    gp = grid::gpar(fill = background, col = border_color, lwd = border_size)
  )
  layout <- heatmap_enrichment_terms_layout(
    data = data,
    max_width = max_width,
    reserve_axis = show_axis
  )
  labels <- heatmap_enrichment_stat_labels(
    data = data,
    mode = label_mode,
    digits = axis_digits
  )
  has_graphic <- !is.na(data[["stat_label"]])
  for (i in seq_len(nrow(data))) {
    if (!isTRUE(has_graphic[[i]])) {
      next
    }
    heatmap_enrichment_draw_term_stat(
      data = data[i, , drop = FALSE],
      y = layout[["block_center"]][[i]],
      height_npc = max(layout[["line_height"]] * 0.72, 0.016),
      stat_x_npc = 0.04,
      stat_width_npc = 0.91,
      bar_color = bar_color,
      label = labels[[i]],
      label_color = label_color,
      fontsize = data[["fontsize"]][[i]]
    )
  }
  if (isTRUE(show_axis)) {
    heatmap_enrichment_draw_stat_axis(
      x_npc = 0.04,
      width_npc = 0.91,
      max_value = axis_max,
      label = axis_label,
      fontsize = min(data[["fontsize"]]),
      digits = axis_digits
    )
  }
}

heatmap_enrichment_draw_term_stat <- function(
  data,
  y,
  height_npc = NULL,
  stat_x_npc = 0,
  stat_width_npc = 0.3,
  bar_color = NULL,
  label = "",
  label_color = NULL,
  fontsize = 8
) {
  scaled <- data[["stat_scaled"]][[1]]
  if (!is.finite(scaled)) {
    scaled <- 0
  }
  scaled <- max(0, min(1, scaled))
  col <- data[["col"]][[1]]
  fill_color <- bar_color %||% col
  x0 <- stat_x_npc
  x1 <- min(0.96, x0 + stat_width_npc)
  grid::grid.rect(
    x = grid::unit(x0, "npc"),
    y = grid::unit(y, "npc"),
    width = grid::unit((x1 - x0) * scaled, "npc"),
    height = grid::unit(height_npc %||% 0.03, "npc"),
    just = c("left", "center"),
    gp = grid::gpar(
      fill = fill_color,
      col = NA,
      alpha = 1
    )
  )
  if (!is.null(label) && !is.na(label) && nzchar(label)) {
    bar_length <- (x1 - x0) * scaled
    label_inside <- bar_length >= 0.24
    resolved_label_color <- label_color
    if (is.null(resolved_label_color) || is.na(resolved_label_color)) {
      resolved_label_color <- if (isTRUE(label_inside)) {
        heatmap_enrichment_contrast_color(fill_color)
      } else {
        "black"
      }
    }
    grid::grid.text(
      label = label,
      x = grid::unit(
        if (isTRUE(label_inside)) {
          x0 + bar_length / 2
        } else {
          min(x1, x0 + bar_length + 0.025)
        },
        "npc"
      ),
      y = grid::unit(y, "npc"),
      just = if (isTRUE(label_inside)) "center" else c("left", "center"),
      gp = grid::gpar(
        col = resolved_label_color,
        fontsize = max(fontsize - 1, 5)
      )
    )
  }
}

heatmap_enrichment_contrast_color <- function(
  color,
  dark = "black",
  light = "white",
  threshold = 0.45
) {
  rgb <- grDevices::col2rgb(color)[, 1] / 255
  rgb <- ifelse(
    rgb <= 0.04045,
    rgb / 12.92,
    ((rgb + 0.055) / 1.055)^2.4
  )
  luminance <- sum(c(0.2126, 0.7152, 0.0722) * rgb)
  if (luminance < threshold) light else dark
}

heatmap_enrichment_draw_stat_axis <- function(
  x_npc = 0,
  width_npc,
  max_value,
  label,
  fontsize = 8,
  digits = 2
) {
  if (!is.finite(max_value) || max_value <= 0) {
    max_value <- 1
  }
  at <- c(0, 0.5, 1)
  values <- at * max_value
  y <- 0.24
  grid::grid.lines(
    x = grid::unit(x_npc + c(0, width_npc), "npc"),
    y = grid::unit(c(y, y), "npc"),
    gp = grid::gpar(col = "grey25", lwd = 0.8)
  )
  for (i in seq_along(at)) {
    x <- x_npc + width_npc * at[[i]]
    grid::grid.lines(
      x = grid::unit(c(x, x), "npc"),
      y = grid::unit(c(y, y + 0.035), "npc"),
      gp = grid::gpar(col = "grey25", lwd = 0.8)
    )
    grid::grid.text(
      label = heatmap_enrichment_format_stat(values[[i]], digits = digits),
      x = grid::unit(x, "npc"),
      y = grid::unit(y - 0.06, "npc"),
      just = c("center", "top"),
      gp = grid::gpar(col = "grey20", fontsize = max(fontsize - 1, 5))
    )
  }
  if (!is.null(label) && nzchar(label)) {
    grid::grid.text(
      label = label,
      x = grid::unit(x_npc + width_npc / 2, "npc"),
      y = grid::unit(0.055, "npc"),
      just = c("center", "center"),
      gp = grid::gpar(col = "grey20", fontsize = max(fontsize - 2, 5))
    )
  }
}

heatmap_enrichment_term_stat <- function(
  df,
  metric,
  score,
  terms_stat = "none"
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
    pvalue <- if (metric %in% colnames(df)) {
      suppressWarnings(as.numeric(df[[metric]]))
    } else {
      rep(NA_real_, length(score))
    }
    axis_label <- heatmap_enrichment_p_axis_label(metric)
  } else {
    if (!stat_key %in% colnames(df)) {
      log_message(
        "{.arg terms_stat} must be {.val none}, {.val score}, or one enrichment column. ",
        "Available columns include: {.val {colnames(df)}}",
        message_type = "error"
      )
    }
    value <- df[[stat_key]]
    is_pvalue <- stat_key %in% c("pvalue", "p.adjust", "qvalue")
    pvalue <- if (is_pvalue) {
      suppressWarnings(as.numeric(value))
    } else {
      rep(NA_real_, length(value))
    }
    strength <- if (is_pvalue) {
      -log10(pmax(pvalue, .Machine$double.xmin))
    } else {
      suppressWarnings(as.numeric(value))
    }
    label <- stat_key
    axis_label <- if (is_pvalue) {
      heatmap_enrichment_p_axis_label(stat_key)
    } else {
      stat_key
    }
  }
  list(
    label = label,
    strength = strength,
    pvalue = pvalue,
    axis_label = axis_label
  )
}

heatmap_enrichment_p_axis_label <- function(metric) {
  if (metric %in% c("p.adjust", "qvalue")) {
    return("-log10(Q)")
  }
  if (metric %in% "pvalue") {
    return("-log10(P)")
  }
  paste0("-log10(", metric, ")")
}

heatmap_enrichment_p_significance <- function(
  pvalue,
  cutoffs = c(`*` = 5e-2, `**` = 1e-2, `***` = 1e-3, `****` = 1e-4)
) {
  pvalue <- suppressWarnings(as.numeric(pvalue))
  out <- rep("", length(pvalue))
  for (label in names(cutoffs)) {
    out[!is.na(pvalue) & pvalue < cutoffs[[label]]] <- label
  }
  out
}

heatmap_enrichment_stat_labels <- function(data, mode = "value", digits = 2) {
  mode <- match.arg(
    as.character(mode %||% "none")[[1]],
    c("none", "value", "significance", "both")
  )
  value <- heatmap_enrichment_format_stat(data[["stat_strength"]], digits = digits)
  significance <- data[["stat_signif"]]
  significance[is.na(significance)] <- ""
  switch(mode,
    none = rep("", nrow(data)),
    value = value,
    significance = significance,
    both = trimws(paste(value, significance))
  )
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
