enrichment_wordcloud_plot <- function(
  df, metric, word_type, words_excluded, topWord, word_size, palette, palcolor,
  score_mode = FALSE, score_label = "Score", gene_col = "geneID"
) {
  term_words <- word_type == "term"
  if (term_words && !score_mode && any(df$Database %in% c("GO", "GO_BP", "GO_CC", "GO_MF"))) {
    check_r("simplifyEnrichment", verbose = FALSE)
  }
  df_groups <- if (term_words) {
    groups <- split(df, list(df$Database, df$Groups))
    groups[sapply(groups, nrow) > 0]
  } else {
    list(df)
  }
  df_groups <- lapply(df_groups, function(df_sub) {
    if (term_words && !score_mode && all(df_sub$Database %in% c("GO", "GO_BP", "GO_CC", "GO_MF"))) {
      keywords <- enrichment_go_keywords(df_sub[["ID"]])
      if (nrow(keywords) == 0L) {
        return(NULL)
      }
      words <- keywords |>
        dplyr::reframe(
          keyword = .data[["keyword"]],
          score = enrichment_metric_weight(.data[["padj"]]),
          count = .data[["n_term"]],
          Database = df_sub[["Database"]][1],
          Groups = df_sub[["Groups"]][1]
        )
    } else {
      text <- if (term_words) {
        tolower(as.character(df_sub[["Description"]]))
      } else {
        as.character(df_sub[[gene_col]])
      }
      df_sub[["keyword"]] <- strsplit(text, if (term_words) " " else "/")
      words <- df_sub |>
        thisutils::unnest_fun(cols = "keyword") |>
        dplyr::group_by(.data[["keyword"]], Database, Groups) |>
        dplyr::reframe(
          keyword = .data[["keyword"]],
          score = sum(enrichment_metric_weight(.data[[metric]], score_mode)),
          count = dplyr::n(),
          Database = .data[["Database"]],
          Groups = .data[["Groups"]],
          .groups = "keep"
        )
    }
    if (term_words) {
      words <- words |>
        dplyr::filter(!grepl(pattern = "\\[.*\\]", x = .data[["keyword"]])) |>
        dplyr::filter(nchar(.data[["keyword"]]) >= 1) |>
        dplyr::filter(!tolower(.data[["keyword"]]) %in% tolower(words_excluded))
    }
    words <- words |>
      dplyr::distinct() |>
      dplyr::mutate(
        angle = 90 * sample(c(0, 1), dplyr::n(), replace = TRUE, prob = c(60, 40))
      ) |>
      as.data.frame()
    words[
      utils::head(order(words[["score"]], decreasing = TRUE), topWord), ,
      drop = FALSE
    ]
  })
  df <- do.call(rbind, df_groups)
  colors <- palette_colors(
    df[["score"]],
    type = "continuous",
    palette = palette,
    palcolor = palcolor,
    matched = FALSE
  )
  colors_value <- seq(
    min(df[["score"]], na.rm = TRUE),
    stats::quantile(df[["score"]], 0.99, na.rm = TRUE) + 0.001,
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
      name = paste0(if (score_mode) score_label else "Score", ":"),
      colours = colors,
      values = scales::rescale(colors_value),
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
    )
  p
}
