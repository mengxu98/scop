# Render already-filtered ORA results without recomputing enrichment.
enrichment_rank_plot <- function(
  res, db, group_use, score_col, topTerm, label_n, label_terms,
  palcolor, palette, point_size, point_alpha, label_size, label_width,
  shared_y, theme_use, theme_args, combine, nrow, ncol, byrow, seed,
  legend.position, legend.direction
) {
  d <- res$enrichment
  required <- c("Groups", "Database", "ID", "Description", "p.adjust", score_col)
  if (!all(required %in% names(d)) || length(score_col) != 1L ||
    !is.numeric(d[[score_col]]) || !is.numeric(d$p.adjust)) {
    log_message("Rank plots require Groups, Database, ID, Description, numeric p.adjust and rank_by columns", message_type = "error")
  }
  if (length(db) != 1L) {
    log_message("Rank plots require exactly one database; select one with db", message_type = "error")
  }
  if (length(topTerm) != 1L || is.na(topTerm) || topTerm < 1 ||
    (is.finite(topTerm) && topTerm != floor(topTerm))) {
    log_message("topTerm must be a positive integer or Inf", message_type = "error")
  }
  if (length(label_n) != 1L || !is.finite(label_n) || label_n < 0 || label_n != floor(label_n)) {
    log_message("rank_nlabel must be a nonnegative integer", message_type = "error")
  }
  if (!is.numeric(point_size) || length(point_size) != 2L ||
    any(!is.finite(point_size)) || any(point_size <= 0) || point_size[1] > point_size[2]) {
    log_message("rank_pt_size must contain two positive increasing sizes", message_type = "error")
  }
  if (length(point_alpha) != 1L || !is.finite(point_alpha) || point_alpha < 0 || point_alpha > 1 ||
    length(label_size) != 1L || !is.finite(label_size) || label_size <= 0 ||
    length(label_width) != 1L || !is.finite(label_width) || label_width < 1 ||
    length(shared_y) != 1L || is.na(shared_y) || !is.logical(shared_y)) {
    log_message("Invalid rank opacity, label size, character width, or shared-y setting", message_type = "error")
  }
  if (is.null(group_use)) group_use <- unique(as.character(d$Groups))
  group_use <- unique(as.character(group_use))
  if (is.list(label_terms) && (is.null(names(label_terms)) ||
    any(!names(label_terms) %in% group_use))) {
    log_message("rank_label_terms lists must be named by selected groups", message_type = "error")
  }
  valid <- is.finite(d[[score_col]]) & is.finite(d$p.adjust) &
    d$p.adjust >= 0 & d$p.adjust <= 1 & !is.na(d$ID) & !is.na(d$Groups)
  if (any(!valid)) log_message("Omitting rank rows with missing identifiers or invalid scores/p.adjust", message_type = "warning")
  d <- d[valid, , drop = FALSE]
  if (!nrow(d)) log_message("No valid terms remain for the rank plot", message_type = "error")
  d$.score <- d[[score_col]]
  # Preserve the original adjusted P; cap only its plotted log transform at zero.
  d$.significance <- -log10(pmax(d$p.adjust, .Machine$double.xmin))
  if (is.null(palcolor) || is.null(names(palcolor))) {
    palcolor <- thisplot::palette_colors(factor(group_use, levels = group_use),
      palette = palette, palcolor = palcolor, type = "discrete"
    )
  }
  if (!all(group_use %in% names(palcolor))) {
    log_message("Named palcolor must contain every selected group", message_type = "error")
  }
  title <- "Enrichment rank"
  label.fg <- "black"
  label.bg <- "white"
  label.bg.r <- 0.1
  parts <- lapply(group_use, function(g) {
    z <- d[d$Groups == g, , drop = FALSE]
    z <- z[order(-z$.score, z$p.adjust, z$ID), , drop = FALSE]
    z <- z[!duplicated(z$ID), , drop = FALSE]
    z <- utils::head(z, topTerm)
    z$.rank <- seq_len(nrow(z))
    terms <- if (is.list(label_terms)) label_terms[[g]] else label_terms
    z$.highlight <- if (is.null(label_terms)) {
      z$.rank <= label_n
    } else {
      z$ID %in% terms | z$Description %in% terms
    }
    z$.label <- vapply(as.character(z$Description), function(s) {
      paste(strwrap(s, width = label_width), collapse = "\n")
    }, character(1))
    z
  })
  plotted <- do.call(rbind, parts)
  ylim <- range(plotted$.score)
  span <- diff(ylim)
  if (span == 0) span <- max(abs(ylim), 1) * 0.1
  ylim <- ylim + c(-0.06, 0.65) * span
  size_limits <- range(plotted$.significance)
  if (diff(size_limits) == 0) size_limits <- size_limits + c(-0.01, 0.01)
  panels <- lapply(seq_along(group_use), function(i) {
    z <- parts[[i]]
    g <- group_use[i]
    color <- palcolor[[g]]
    if (!nrow(z)) {
      return(ggplot2::ggplot() +
        ggplot2::annotate("text",
          x = 0, y = 0,
          label = "No terms passing cutoff"
        ) +
        ggplot2::theme_void() +
        ggplot2::labs(title = g))
    }
    z$.panel <- g
    hi <- z[z$.highlight, , drop = FALSE]
    panel_ylim <- ylim
    if (!shared_y) {
      yr <- range(z$.score)
      ys <- max(diff(yr), max(abs(yr), 1) * 0.1)
      panel_ylim <- yr + c(-0.06, 0.65) * ys
    }
    xmax <- max(z$.rank)
    # Ordered callouts, inspired by heatmap anno_mark. Lines anchor to the
    # original point coordinates; labels never change rank or enrichment.
    if (nrow(hi)) {
      line_count <- lengths(strsplit(hi$.label, "\n", fixed = TRUE))
      weights <- line_count + 1.2
      centers <- cumsum(weights) - weights / 2
      hi$.label_y <- panel_ylim[2] - diff(panel_ylim) *
        (0.06 + 0.66 * centers / sum(weights))
      hi$.label_x <- max(xmax * 0.48, max(hi$.rank) + xmax * 0.12)
    }
    p <- ggplot2::ggplot(z, ggplot2::aes(.data[[".rank"]], .data[[".score"]])) +
      ggplot2::geom_point(ggplot2::aes(size = .data[[".significance"]]),
        shape = 16, colour = color, alpha = point_alpha
      ) +
      ggplot2::geom_point(
        data = hi, ggplot2::aes(size = .data[[".significance"]]),
        shape = 16, colour = color
      ) +
      ggplot2::facet_wrap(~.panel) +
      ggplot2::scale_size_continuous(
        name = expression(-log[10](adjusted ~ P)),
        range = point_size, limits = size_limits
      ) +
      ggplot2::guides(size = ggplot2::guide_legend(
        override.aes = list(colour = "grey45", alpha = 1, shape = 16)
      )) +
      ggplot2::scale_x_continuous(breaks = unique(round(seq(1, xmax, length.out = 4)))) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.35))) +
      ggplot2::labs(x = "Pathway rank", y = if (score_col == "FoldEnrichment") {
        "Fold enrichment"
      } else {
        score_col
      }) +
      apply_plot_theme(theme_use, theme_args) +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(fill = NA, colour = "#30343B", linewidth = 0.45),
        panel.grid = ggplot2::element_blank(),
        axis.text = ggplot2::element_text(colour = "#4F5661", size = 9),
        axis.title = ggplot2::element_text(colour = "#202328", size = 11),
        axis.ticks = ggplot2::element_line(colour = "#30343B", linewidth = 0.4),
        strip.background = ggplot2::element_rect(fill = color, colour = "#30343B", linewidth = 0.45),
        strip.text = ggplot2::element_text(size = 12, face = "bold"),
        plot.margin = ggplot2::margin(8, 12, 8, 8)
      )
    if (nrow(hi)) {
      p <- p +
        ggrepel::geom_text_repel(
          data = hi,
          ggplot2::aes(x = .data[[".rank"]], y = .data[[".score"]], label = .data[[".label"]]),
          inherit.aes = FALSE, hjust = 0, size = label_size,
          lineheight = 1.05, colour = label.fg,
          bg.color = label.bg, bg.r = label.bg.r,
          family = "sans", fontface = "plain",
          nudge_x = hi$.label_x - hi$.rank,
          nudge_y = hi$.label_y - hi$.score,
          direction = "y", force = 1, force_pull = 0,
          box.padding = 0.3, point.padding = 0.15,
          min.segment.length = 0, segment.size = 0.25,
          max.overlaps = Inf, seed = seed, segment.colour = "grey40"
        )
    }
    xright <- if (nrow(hi)) max(xmax * 1.42, hi$.label_x[1] + xmax * 0.94) else xmax * 1.1
    p <- p + ggplot2::coord_cartesian(
      xlim = c(1 - xmax * 0.045, xright), ylim = panel_ylim, expand = FALSE
    )
    p
  })
  names(panels) <- group_use
  if (!combine) {
    return(panels)
  }
  p <- patchwork::wrap_plots(panels, nrow = nrow, ncol = ncol, byrow = byrow, guides = "collect") +
    patchwork::plot_annotation(title = paste0(title, " | ", db))
  p <- p & ggplot2::theme(legend.position = legend.position, legend.direction = legend.direction)
  attr(p, "enrichment_data") <- plotted
  p
}
