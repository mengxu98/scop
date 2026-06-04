#' @title Ferris Wheel Plot
#'
#' @description
#' Draw a Ferris wheel-style donut plot for up/down counts across pathways or
#' other categories. The center donut shows the total counts, and outer donuts
#' show category-level counts with radius scaled by category size.
#'
#' @md
#' @param data A `data.frame` containing category labels and two non-negative
#' count columns. If `NULL`, the table is built from `res`.
#' @param res A result list returned by [RunEnrichment()] or an enrichment
#' `data.frame`. When provided, `FerrisWheelPlot()` summarizes enriched terms
#' into up/down counts internally.
#' @param label_col Column name for category labels.
#' Default is `"pathway"`.
#' @param up_col Column name for up-regulated counts.
#' Default is `"up"`.
#' @param down_col Column name for down-regulated counts.
#' Default is `"down"`.
#' @param total_up,total_down Total counts shown in the center donut. When
#' `NULL`, the corresponding column sum is used. For enrichment input, these
#' can be inferred from `de_results`.
#' @param de_results Optional differential expression result table used to infer
#' `total_up` and `total_down` from `avg_log2FC`.
#' @param term_col,group_col,count_col,padj_col Column names used when `res` is
#' an enrichment result. Defaults match [RunEnrichment()] output.
#' @param up_group,down_group Group labels in enrichment results. Default
#' matches examples that use `"Up"` and `"Down"` as `geneID_groups`.
#' @param padj_cutoff Adjusted p-value cutoff used to select enriched terms
#' from `res`.
#' @param top_n Number of enriched terms to show when `res` is provided.
#' Default is `10`.
#' @param up_label,down_label Legend labels and center text labels.
#' Defaults are `"Up-regulated Genes"` and `"Down-regulated Genes"`.
#' @param up_color,down_color Optional colors for up/down counts. When `NULL`,
#' colors are selected from the `"Chinese"` palette used by `scop`.
#' @param palette,palcolor Palette passed to [thisplot::palette_colors()] when
#' `up_color` or `down_color` is not provided.
#' @param center_label Text shown inside the center donut.
#' Default is `"Significant Genes"`.
#' @param center_radius Radius of the center donut.
#' Default is `1.25`.
#' @param center_width Width of the center donut band.
#' Default is `0.28`.
#' @param outer_distance Distance from plot center to outer donut centers.
#' Default is `2.55`.
#' @param outer_width Width of outer donut bands.
#' Default is `0.12`.
#' @param min_outer_radius,max_outer_radius Minimum and maximum outer donut
#' radii after scaling by category count.
#' @param label_wrap Maximum characters per category label line. Set `NULL` to
#' disable wrapping.
#' @param label_case Label case for category names. `"title"` capitalizes words
#' by default; `"none"` keeps labels unchanged.
#' @param label.stroke White outline width around text labels. Default is
#' `0.1`. Set `0` to disable the outline.
#' @param label.stroke.color Outline color for text labels.
#' @param label.size,number.size,center.size,legend.size Text sizes.
#' @param line.color,line.linewidth Connector line styling.
#' @param legend.position Position of the compact legend. Use `"bottom"` or
#' `"none"`.
#' @param direction Drawing direction for the outer donuts. `1` is clockwise and
#' `-1` is counterclockwise.
#' @param start Start angle in radians for the first category.
#' Default is `pi / 2`.
#' @param theme_use Theme function name.
#' Default is `"theme_scop"`.
#' @param theme_args Additional arguments passed to `theme_use`.
#'
#' @return A `ggplot` object.
#'
#' @seealso [EnrichmentPlot], [DEtestPlot]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunDEtest(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   only.pos = FALSE
#' )
#'
#' de_df <- pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox
#' de_df <- de_df[
#'   de_df$p_val_adj < 0.05 & abs(de_df$avg_log2FC) > 0.25,
#'   ,
#'   drop = FALSE
#' ]
#' de_df$direction <- ifelse(de_df$avg_log2FC > 0, "Up", "Down")
#'
#' enrich_out <- RunEnrichment(
#'   geneID = de_df$gene,
#'   geneID_groups = de_df$direction,
#'   db = "GO_BP",
#'   species = "Mus_musculus"
#' )
#' FerrisWheelPlot(
#'   res = enrich_out,
#'   de_results = de_df
#' )
FerrisWheelPlot <- function(
  data = NULL,
  res = NULL,
  label_col = "pathway",
  up_col = "up",
  down_col = "down",
  total_up = NULL,
  total_down = NULL,
  de_results = NULL,
  term_col = "Description",
  group_col = "Groups",
  count_col = "Count",
  padj_col = "p.adjust",
  up_group = "Up",
  down_group = "Down",
  padj_cutoff = 0.05,
  top_n = 10,
  up_label = "Up-regulated Genes",
  down_label = "Down-regulated Genes",
  up_color = NULL,
  down_color = NULL,
  palette = "Chinese",
  palcolor = NULL,
  center_label = "Significant Genes",
  center_radius = 1.25,
  center_width = 0.28,
  outer_distance = 2.55,
  outer_width = 0.12,
  min_outer_radius = 0.28,
  max_outer_radius = 0.48,
  label_wrap = 26,
  label_case = c("title", "none"),
  label.stroke = 0.2,
  label.stroke.color = "white",
  label.size = 3,
  number.size = 3,
  center.size = 4.5,
  legend.size = 3.5,
  line.color = "grey35",
  line.linewidth = 0.35,
  legend.position = c("bottom", "none"),
  direction = 1,
  start = pi / 2,
  theme_use = "theme_scop",
  theme_args = list()
) {
  legend.position <- match.arg(legend.position)
  label_case <- match.arg(label_case)
  if (is.null(data) && is.null(res)) {
    log_message("Either {.arg data} or {.arg res} must be provided", message_type = "error")
  }
  if (!is.null(res)) {
    data <- ferris_wheel_enrichment_data(
      res = res,
      term_col = term_col,
      group_col = group_col,
      count_col = count_col,
      padj_col = padj_col,
      up_group = up_group,
      down_group = down_group,
      padj_cutoff = padj_cutoff,
      top_n = top_n
    )
    label_col <- "pathway"
    up_col <- "up"
    down_col <- "down"
    if (!is.null(de_results)) {
      totals <- ferris_wheel_de_totals(de_results)
      total_up <- total_up %||% totals[["up"]]
      total_down <- total_down %||% totals[["down"]]
    }
  }
  if (!is.data.frame(data)) {
    log_message("{.arg data} must be a data.frame", message_type = "error")
  }
  required_cols <- c(label_col, up_col, down_col)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0L) {
    log_message(
      "{.arg data} is missing required column(s): {.val {missing_cols}}",
      message_type = "error"
    )
  }
  if (!direction %in% c(-1, 1)) {
    log_message(
      "{.arg direction} must be {.val 1} or {.val -1}",
      message_type = "error"
    )
  }
  radius_args <- c(
    center_radius,
    center_width,
    outer_distance,
    outer_width,
    min_outer_radius,
    max_outer_radius
  )
  if (any(!is.finite(radius_args)) || any(radius_args <= 0)) {
    log_message(
      "Radius and width arguments must be finite positive numbers",
      message_type = "error"
    )
  }
  if (center_width >= center_radius) {
    log_message(
      "{.arg center_width} must be smaller than {.arg center_radius}",
      message_type = "error"
    )
  }
  if (outer_width >= min_outer_radius) {
    log_message(
      "{.arg outer_width} must be smaller than {.arg min_outer_radius}",
      message_type = "error"
    )
  }
  if (!is.finite(label.stroke) || label.stroke < 0) {
    log_message("{.arg label.stroke} must be a finite non-negative number", message_type = "error")
  }

  df <- data.frame(
    label_raw = as.character(data[[label_col]]),
    up = suppressWarnings(as.numeric(data[[up_col]])),
    down = suppressWarnings(as.numeric(data[[down_col]])),
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(df$up)) || any(!is.finite(df$down))) {
    log_message(
      "{.arg up_col} and {.arg down_col} must contain finite numeric values",
      message_type = "error"
    )
  }
  if (any(df$up < 0 | df$down < 0)) {
    log_message(
      "{.arg up_col} and {.arg down_col} must contain non-negative values",
      message_type = "error"
    )
  }
  df$total <- df$up + df$down
  df <- df[df$total > 0, , drop = FALSE]
  if (nrow(df) == 0L) {
    log_message("No category has positive total counts", message_type = "error")
  }

  total_up <- total_up %||% sum(df$up, na.rm = TRUE)
  total_down <- total_down %||% sum(df$down, na.rm = TRUE)
  total_up <- suppressWarnings(as.numeric(total_up))
  total_down <- suppressWarnings(as.numeric(total_down))
  if (
    length(total_up) != 1L ||
      length(total_down) != 1L ||
      any(!is.finite(c(total_up, total_down))) ||
      any(c(total_up, total_down) < 0)
  ) {
    log_message(
      "{.arg total_up} and {.arg total_down} must be finite non-negative numbers",
      message_type = "error"
    )
  }

  cols <- palette_colors(palette = palette, palcolor = palcolor, n = 8)
  up_color <- up_color %||% cols[min(6L, length(cols))]
  down_color <- down_color %||% cols[min(2L, length(cols))]

  n <- nrow(df)
  angle <- start - direction * seq(0, 2 * pi - 2 * pi / n, length.out = n)
  df$x <- outer_distance * cos(angle)
  df$y <- outer_distance * sin(angle)
  df$radius <- if (n == 1L || diff(range(sqrt(df$total))) == 0) {
    rep(mean(c(min_outer_radius, max_outer_radius)), n)
  } else {
    scales::rescale(sqrt(df$total), to = c(min_outer_radius, max_outer_radius))
  }
  df$label_x <- (outer_distance + 0.5) * cos(angle)
  df$label_y <- (outer_distance + 0.5) * sin(angle)
  df$label_display <- if (identical(label_case, "title")) {
    ferris_wheel_title_case(df$label_raw)
  } else {
    df$label_raw
  }
  df$label <- if (is.null(label_wrap)) {
    df$label_display
  } else {
    vapply(
      df$label_display,
      function(x) paste(strwrap(x, width = label_wrap), collapse = "\n"),
      character(1)
    )
  }
  df$hjust <- ifelse(cos(angle) > 0.15, 0, ifelse(cos(angle) < -0.15, 1, 0.5))
  df$line_x <- (center_radius + 0.05) * cos(angle)
  df$line_y <- (center_radius + 0.05) * sin(angle)
  df$line_xend <- (outer_distance - df$radius - 0.04) * cos(angle)
  df$line_yend <- (outer_distance - df$radius - 0.04) * sin(angle)

  center_total <- total_up + total_down
  center_split <- if (center_total > 0) total_up / center_total else 0.5
  arc_df <- rbind(
    ferris_wheel_arc(
      0,
      0,
      center_radius - center_width,
      center_radius,
      start,
      start - direction * 2 * pi * center_split,
      up_color,
      "center_up"
    ),
    ferris_wheel_arc(
      0,
      0,
      center_radius - center_width,
      center_radius,
      start - direction * 2 * pi * center_split,
      start - direction * 2 * pi,
      down_color,
      "center_down"
    ),
    do.call(
      rbind,
      lapply(seq_len(n), function(i) {
        split <- df$up[i] / df$total[i]
        rbind(
          ferris_wheel_arc(
            df$x[i],
            df$y[i],
            df$radius[i] - outer_width,
            df$radius[i],
            start,
            start - direction * 2 * pi * split,
            up_color,
            paste0("outer_up_", i)
          ),
          ferris_wheel_arc(
            df$x[i],
            df$y[i],
            df$radius[i] - outer_width,
            df$radius[i],
            start - direction * 2 * pi * split,
            start - direction * 2 * pi,
            down_color,
            paste0("outer_down_", i)
          )
        )
      })
    )
  )

  xlim <- c(
    min(df$label_x, -outer_distance - max_outer_radius - 1.7, na.rm = TRUE),
    max(df$label_x, outer_distance + max_outer_radius + 1.7, na.rm = TRUE)
  )
  ylim <- c(
    min(
      df$label_y,
      -outer_distance -
        max_outer_radius -
        ifelse(legend.position == "bottom", 0.95, 0.45),
      na.rm = TRUE
    ),
    max(df$label_y, outer_distance + max_outer_radius + 0.45, na.rm = TRUE)
  )
  stroke_unit <- min(diff(xlim), diff(ylim)) / 350

  label_stroke_layers <- ferris_wheel_text_stroke_layer(
    data = df,
    mapping = ggplot2::aes(
      x = .data$label_x,
      y = .data$label_y,
      label = .data$label,
      hjust = .data$hjust
    ),
    label.size = label.size,
    label.stroke = label.stroke,
    label.stroke.color = label.stroke.color,
    stroke.unit = stroke_unit,
    lineheight = 0.95
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = df,
      ggplot2::aes(
        x = .data$line_x,
        y = .data$line_y,
        xend = .data$line_xend,
        yend = .data$line_yend
      ),
      linewidth = line.linewidth,
      color = line.color
    ) +
    ggforce::geom_arc_bar(
      data = arc_df,
      ggplot2::aes(
        x0 = .data$x0,
        y0 = .data$y0,
        r0 = .data$r0,
        r = .data$r,
        start = .data$start,
        end = .data$end,
        fill = .data$fill,
        group = .data$group
      ),
      color = "white",
      linewidth = 0.55
    ) +
    ferris_wheel_text_stroke_layer(
      data = df,
      ggplot2::aes(x = .data$x, y = .data$y + 0.08, label = .data$up),
      label.size = number.size,
      label.stroke = label.stroke,
      label.stroke.color = label.stroke.color,
      stroke.unit = stroke_unit
    ) +
    ggplot2::geom_text(
      data = df,
      ggplot2::aes(x = .data$x, y = .data$y + 0.08, label = .data$up),
      color = up_color,
      size = number.size
    ) +
    ferris_wheel_text_stroke_layer(
      data = df,
      ggplot2::aes(x = .data$x, y = .data$y - 0.12, label = .data$down),
      label.size = number.size,
      label.stroke = label.stroke,
      label.stroke.color = label.stroke.color,
      stroke.unit = stroke_unit
    ) +
    ggplot2::geom_text(
      data = df,
      ggplot2::aes(x = .data$x, y = .data$y - 0.12, label = .data$down),
      color = down_color,
      size = number.size
    ) +
    label_stroke_layers +
    ggplot2::geom_text(
      data = df,
      ggplot2::aes(
        x = .data$label_x,
        y = .data$label_y,
        label = .data$label,
        hjust = .data$hjust
      ),
      size = label.size,
      lineheight = 0.95
    ) +
    ferris_wheel_text_stroke_layer(
      data = data.frame(x = 0, y = 0.18, label = center_label),
      mapping = ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      label.size = center.size,
      label.stroke = label.stroke,
      label.stroke.color = label.stroke.color,
      stroke.unit = stroke_unit
    ) +
    ggplot2::geom_text(
      data = data.frame(x = 0, y = 0.18, label = center_label),
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = center.size
    ) +
    ferris_wheel_text_stroke_layer(
      data = data.frame(x = 0, y = -0.07, label = paste0("Up = ", total_up)),
      mapping = ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      label.size = center.size * 0.9,
      label.stroke = label.stroke,
      label.stroke.color = label.stroke.color,
      stroke.unit = stroke_unit
    ) +
    ggplot2::geom_text(
      data = data.frame(x = 0, y = -0.07, label = paste0("Up = ", total_up)),
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      color = up_color,
      size = center.size * 0.9
    ) +
    ferris_wheel_text_stroke_layer(
      data = data.frame(x = 0, y = -0.32, label = paste0("Down = ", total_down)),
      mapping = ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      label.size = center.size * 0.9,
      label.stroke = label.stroke,
      label.stroke.color = label.stroke.color,
      stroke.unit = stroke_unit
    ) +
    ggplot2::geom_text(
      data = data.frame(x = 0, y = -0.32, label = paste0("Down = ", total_down)),
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      color = down_color,
      size = center.size * 0.9
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::xlim(xlim) +
    ggplot2::ylim(ylim) +
    ggplot2::labs(x = NULL, y = NULL)

  if (legend.position == "bottom") {
    legend_df <- data.frame(
      x = min(xlim) + diff(xlim) * 0.42,
      y = c(min(ylim) + 0.42, min(ylim) + 0.18),
      label = c(up_label, down_label),
      fill = c(up_color, down_color)
    )
    p <- p +
      ggplot2::geom_tile(
        data = legend_df,
        ggplot2::aes(x = .data$x, y = .data$y, fill = .data$fill),
        width = 0.28,
        height = 0.12
      ) +
      ferris_wheel_text_stroke_layer(
        data = legend_df,
        ggplot2::aes(x = .data$x + 0.22, y = .data$y, label = .data$label),
        label.size = legend.size,
        label.stroke = label.stroke,
        label.stroke.color = label.stroke.color,
        stroke.unit = stroke_unit,
        hjust = 0
      ) +
      ggplot2::geom_text(
        data = legend_df,
        ggplot2::aes(x = .data$x + 0.22, y = .data$y, label = .data$label),
        hjust = 0,
        size = legend.size
      )
  }

  p +
    do.call(theme_use, theme_args) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      legend.position = "none",
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

ferris_wheel_arc <- function(x, y, r0, r, start, end, fill, group) {
  data.frame(
    x0 = x,
    y0 = y,
    r0 = r0,
    r = r,
    start = start,
    end = end,
    fill = fill,
    group = group,
    stringsAsFactors = FALSE
  )
}

ferris_wheel_enrichment_data <- function(
  res,
  term_col = "Description",
  group_col = "Groups",
  count_col = "Count",
  padj_col = "p.adjust",
  up_group = "Up",
  down_group = "Down",
  padj_cutoff = 0.05,
  top_n = 10
) {
  enrichment <- if (is.list(res) && !is.null(res[["enrichment"]])) {
    res[["enrichment"]]
  } else {
    res
  }
  if (!is.data.frame(enrichment)) {
    log_message("{.arg res} must be a RunEnrichment result list or an enrichment data.frame", message_type = "error")
  }
  required_cols <- c(term_col, group_col, count_col)
  missing_cols <- setdiff(required_cols, colnames(enrichment))
  if (length(missing_cols) > 0L) {
    log_message(
      "{.arg res} is missing required column(s): {.val {missing_cols}}",
      message_type = "error"
    )
  }
  if (padj_col %in% colnames(enrichment)) {
    padj <- suppressWarnings(as.numeric(enrichment[[padj_col]]))
    enrichment <- enrichment[is.finite(padj) & padj <= padj_cutoff, , drop = FALSE]
  }
  enrichment <- enrichment[
    enrichment[[group_col]] %in% c(up_group, down_group),
    ,
    drop = FALSE
  ]
  if (nrow(enrichment) == 0L) {
    log_message("No enrichment rows remain after filtering", message_type = "error")
  }
  term_order <- if (padj_col %in% colnames(enrichment)) {
    unique(enrichment[[term_col]][order(enrichment[[padj_col]])])
  } else {
    unique(enrichment[[term_col]])
  }
  term_use <- head(term_order, top_n)
  data.frame(
    pathway = term_use,
    up = vapply(term_use, function(term) {
      sum(as.numeric(enrichment[[count_col]][enrichment[[term_col]] == term & enrichment[[group_col]] == up_group]), na.rm = TRUE)
    }, numeric(1)),
    down = vapply(term_use, function(term) {
      sum(as.numeric(enrichment[[count_col]][enrichment[[term_col]] == term & enrichment[[group_col]] == down_group]), na.rm = TRUE)
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
}

ferris_wheel_de_totals <- function(de_results) {
  if (!is.data.frame(de_results) || !"avg_log2FC" %in% colnames(de_results)) {
    log_message("{.arg de_results} must be a data.frame containing {.field avg_log2FC}", message_type = "error")
  }
  fc <- suppressWarnings(as.numeric(de_results[["avg_log2FC"]]))
  c(up = sum(fc > 0, na.rm = TRUE), down = sum(fc < 0, na.rm = TRUE))
}

ferris_wheel_title_case <- function(x) {
  x <- as.character(x)
  vapply(x, function(value) {
    pieces <- strsplit(value, "\\s+")[[1]]
    pieces <- vapply(pieces, function(piece) {
      if (!nzchar(piece)) {
        return(piece)
      }
      if (grepl("[[:upper:]0-9]", piece)) {
        return(piece)
      }
      sub("^([[:alpha:]])", "\\U\\1", tolower(piece), perl = TRUE)
    }, character(1))
    paste(pieces, collapse = " ")
  }, character(1))
}

ferris_wheel_text_stroke_layer <- function(
  data,
  mapping,
  label.size,
  label.stroke = 0.1,
  label.stroke.color = "white",
  stroke.unit = 0.01,
  lineheight = 1,
  ...
) {
  if (label.stroke <= 0) {
    return(NULL)
  }
  stroke <- label.stroke * stroke.unit
  offsets <- data.frame(
    dx = c(-stroke, stroke, 0, 0, -stroke, -stroke, stroke, stroke),
    dy = c(0, 0, -stroke, stroke, -stroke, stroke, -stroke, stroke)
  )
  lapply(seq_len(nrow(offsets)), function(i) {
    ggplot2::geom_text(
      data = data,
      mapping = mapping,
      nudge_x = offsets$dx[i],
      nudge_y = offsets$dy[i],
      size = label.size,
      lineheight = lineheight,
      color = label.stroke.color,
      ...
    )
  })
}
