#' Mirror stored enrichment results for two groups
#'
#' Draw independent ranked GO or other database results on opposite sides of zero.
#' No enrichment or differential comparison is computed. Group colors and shared
#' significance opacity encode the saved results; full black term labels may
#' extend past bars. Hit genes appear below each bar and database bands are rounded.
#' @param object Seurat object containing saved RunEnrichment results; NULL when res is supplied.
#' @param res Result list with an enrichment data.frame, as returned by RunEnrichment. Overrides object.
#' @param group.by Saved enrichment grouping key; NULL selects custom results.
#' @param group_use Character vector of exactly two groups, left then right.
#' @param test.use Saved differential-test key; defaults to wilcox.
#' @param db Database names, in display order. GO aliases are expanded using stored results.
#' @param topTerm Positive integer number of terms per database and group, default 5.
#' @param padjustCutoff Adjusted P cutoff in (0, 1]; terms strictly below it are retained.
#' @param palcolor Two group colors; unnamed colors follow group_use, named colors match group names.
#' @param row_step,db_gap Positive row spacing and non-negative extra database gap, in y-axis units.
#' @param bar_height Positive bar height in y-axis units.
#' @param term_size,gene_size Positive term and gene text sizes in mm.
#' @param gene_width Maximum gene-label characters, including ellipsis; at least 4 or Inf. Term names are never truncated.
#' @param alpha_range Increasing two-element opacity range within [0, 1].
#' @param corner_mm Positive database-band corner radius in mm.
#' @param title Optional plot title; NULL creates a title from the selected groups.
#' @param symmetric_axis Whether both sides use the same extent. FALSE adapts
#' each extent while preserving identical numeric unit lengths on both sides.
#' @param band_gap Non-negative gap between the database band and estimated left content extent, in x-axis units.
#' @param base_size Positive theme font size in points.
#' @param font_family Font family used by the theme, term and gene text.
#' @return A ggplot. Its data retains selected original rows plus plotting columns.
#' Terms are ranked independently by adjusted P within each database and group.
#' Missing adjusted P values are omitted with a warning; zero values are capped
#' at .Machine$double.xmin for plotting. Both groups must retain significant terms.
#' @details Negative x positions denote the left group, not negative enrichment.
#' Opposite rows need not represent the same pathway. Asymmetric axis extents
#' include a character-count estimate of label width; large fonts or narrow exports
#' may need additional space. No automatic font-size fitting is performed.
#' @seealso [RunEnrichment()], [EnrichmentPlot()]
#' @export
#' @examples
#' tab <- data.frame(ID = c("a", "b"), Description = c("Process A", "Process B"),
#'   Database = "GO_BP", Groups = c("A", "B"), p.adjust = c(0.01, 0.001),
#'   geneID = c("Gene1/Gene2", "Gene3/Gene4"))
#' EnrichmentMirrorPlot(res = list(enrichment = tab), group_use = c("A", "B"), db = "GO_BP")
EnrichmentMirrorPlot <- function(
  object = NULL, group.by = NULL, group_use = NULL,
  test.use = "wilcox", db = c("GO_BP", "GO_CC", "GO_MF"),
  topTerm = 5, padjustCutoff = 0.05,
  palcolor = c("#2468A2", "#C83E4D"),
  row_step = 0.84, db_gap = 0.40, bar_height = 0.46,
  term_size = 3, gene_size = 2.2, gene_width = 82,
  alpha_range = c(0.55, 1), corner_mm = 2, title = NULL,
  symmetric_axis = FALSE, band_gap = 0.6, base_size = 11, font_family = "sans", res = NULL
) {
  if (length(group_use) != 2L || anyNA(group_use) || anyDuplicated(group_use)) {
    log_message("group_use must specify exactly two distinct groups, left then right", message_type = "error")
  }
  group_use <- as.character(group_use)
  values <- list(topTerm = topTerm, row_step = row_step, bar_height = bar_height,
    term_size = term_size, gene_size = gene_size, corner_mm = corner_mm, base_size = base_size)
  if (any(!vapply(values, function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0, logical(1))) ||
      topTerm != floor(topTerm)) {
    log_message("topTerm must be a positive integer and size/spacing parameters must be positive finite numbers", message_type = "error")
  }
  for (x in list(db_gap, band_gap)) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) {
      log_message("db_gap and band_gap must be finite non-negative numbers", message_type = "error")
    }
  }
  if (!is.logical(symmetric_axis) || length(symmetric_axis) != 1L || is.na(symmetric_axis) ||
      !is.numeric(alpha_range) || length(alpha_range) != 2L || anyNA(alpha_range) ||
      any(alpha_range < 0 | alpha_range > 1) || alpha_range[1] > alpha_range[2]) {
    log_message("symmetric_axis must be logical and alpha_range must be increasing within [0, 1]", message_type = "error")
  }
  if (!is.numeric(gene_width) || length(gene_width) != 1L || is.na(gene_width) || gene_width < 4 ||
      !is.numeric(padjustCutoff) || length(padjustCutoff) != 1L || !is.finite(padjustCutoff) ||
      padjustCutoff <= 0 || padjustCutoff > 1) {
    log_message("gene_width must be at least 4 or Inf; padjustCutoff must be in (0, 1]", message_type = "error")
  }
  if (length(palcolor) != 2L || anyNA(palcolor)) {
    log_message("palcolor must provide two colors", message_type = "error")
  }
  if (is.null(names(palcolor))) names(palcolor) <- group_use
  if (!all(group_use %in% names(palcolor))) log_message("palcolor names must match group_use", message_type = "error")
  if (is.null(res)) {
    if (is.null(object)) log_message("Provide object or res", message_type = "error")
    key <- paste("Enrichment", group.by %||% "custom", test.use, sep = "_")
    if (!key %in% names(object@tools)) log_message("Missing saved enrichment: ", key, message_type = "error")
    res <- object@tools[[key]]
  }
  tab <- res[["enrichment"]]
  required <- c("ID", "Description", "Database", "Groups", "p.adjust", "geneID")
  if (!is.data.frame(tab) || !all(required %in% names(tab)) || !is.numeric(tab$p.adjust)) {
    log_message("Enrichment table requires ID, Description, Database, Groups, numeric p.adjust and geneID", message_type = "error")
  }
  db <- resolve_enrichment_plot_db(db, tab)
  tab <- tab[tab$Groups %in% group_use & tab$Database %in% db, , drop = FALSE]
  if (any(!is.na(tab$p.adjust) & (!is.finite(tab$p.adjust) | tab$p.adjust < 0 | tab$p.adjust > 1))) {
    log_message("Adjusted P values must be between zero and one", message_type = "error")
  }
  if (anyNA(tab$p.adjust)) log_message("Omitting terms with missing adjusted P values", message_type = "warning")
  tab <- tab[!is.na(tab$p.adjust) & tab$p.adjust < padjustCutoff, , drop = FALSE]
  if (!all(group_use %in% tab$Groups)) log_message("Both selected groups must have significant terms", message_type = "error")
  d <- do.call(rbind, lapply(db, function(database) {
    do.call(rbind, lapply(group_use, function(group) {
      z <- tab[tab$Database == database & tab$Groups == group, , drop = FALSE]
      z <- utils::head(z[order(z$p.adjust, z$ID), , drop = FALSE], topTerm)
      z$.rank <- seq_len(nrow(z))
      z
    }))
  }))
  rownames(d) <- NULL
  present <- db[db %in% d$Database]
  n_rows <- vapply(present, function(database) {
    max(table(factor(d$Groups[d$Database == database], levels = group_use)))
  }, numeric(1))
  total_height <- sum(n_rows) * row_step + (length(present) - 1) * db_gap
  offset <- 0
  d$.y <- NA_real_
  bands <- list()
  for (i in seq_along(present)) {
    keep <- d$Database == present[i]
    d$.y[keep] <- total_height - offset - (d$.rank[keep] - 1) * row_step
    yy <- d$.y[keep]
    bands[[i]] <- data.frame(label = sub("^GO_", "", present[i]),
      ymin = min(yy) - 0.42, ymax = max(yy) + 0.38, ymid = mean(range(yy)))
    offset <- offset + n_rows[i] * row_step + db_gap
  }
  bands <- do.call(rbind, bands)
  d$.score <- -log10(pmax(d$p.adjust, .Machine$double.xmin))
  d$.side <- ifelse(d$Groups == group_use[1], -1, 1)
  d$.xmin <- ifelse(d$.side == -1, -d$.score, 0)
  d$.xmax <- ifelse(d$.side == -1, 0, d$.score)
  xmax <- max(1, d$.score)
  d$.text_x <- d$.side * 0.018 * xmax
  d$.hjust <- ifelse(d$.side == -1, 1, 0)
  left_span <- max(d$.score[d$.side == -1], 0.12 * max(nchar(d$Description[d$.side == -1])))
  right_span <- max(d$.score[d$.side == 1], 0.12 * max(nchar(d$Description[d$.side == 1])))
  if (symmetric_axis) left_span <- right_span <- max(left_span, right_span)
  band_right <- -left_span - band_gap
  band_left <- band_right - 0.045 * xmax
  # Pathway names remain complete and black, even when longer than their bars.
  d$.term <- gsub("[\r\n]+", " ", as.character(d$Description))
  d$.genes <- gsub("/", ", ", as.character(d$geneID), fixed = TRUE)
  d$.genes[is.na(d$.genes)] <- ""
  if (is.finite(gene_width)) {
    long <- nchar(d$.genes) > gene_width
    d$.genes[long] <- paste0(substr(d$.genes[long], 1, gene_width - 3), "...")
  }
  RoundBand <- ggplot2::ggproto("RoundBand", ggplot2::GeomRect,
    draw_panel = function(data, panel_params, coord) {
      z <- coord$transform(data, panel_params)
      do.call(grid::grobTree, lapply(seq_len(nrow(z)), function(i) {
        grid::roundrectGrob(
          x = (z$xmin[i] + z$xmax[i]) / 2,
          y = (z$ymin[i] + z$ymax[i]) / 2,
          width = abs(z$xmax[i] - z$xmin[i]),
          height = abs(z$ymax[i] - z$ymin[i]),
          default.units = "native", r = grid::unit(corner_mm, "mm"),
          gp = grid::gpar(fill = "#EEF1F5", col = NA))
      }))
    })
  if (is.null(title)) title <- paste(paste(group_use, collapse = " vs "), "| GO enrichment")
  p <- ggplot2::ggplot(d) +
    ggplot2::layer(data = bands, mapping = ggplot2::aes(ymin = ymin, ymax = ymax),
      geom = RoundBand, stat = "identity", position = "identity", inherit.aes = FALSE,
      params = list(xmin = band_left, xmax = band_right)) +
    ggplot2::geom_text(data = bands, ggplot2::aes(y = ymid, label = label),
      x = (band_left + band_right) / 2, angle = 90, fontface = "bold", colour = "#485466", size = 3.4) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, colour = "#6E7781") +
    ggplot2::geom_rect(ggplot2::aes(xmin = .xmin, xmax = .xmax,
      ymin = .y + 0.10 - bar_height / 2, ymax = .y + 0.10 + bar_height / 2,
      fill = Groups, alpha = .score), colour = NA) +
    ggplot2::geom_text(ggplot2::aes(x = .text_x, y = .y + 0.10,
      label = .term, hjust = .hjust), size = term_size, colour = "black", family = font_family) +
    ggplot2::geom_text(ggplot2::aes(x = .text_x,
      y = .y + 0.10 - bar_height / 2 - 0.15,
      label = .genes, hjust = .hjust, colour = Groups, alpha = .score),
      size = gene_size, family = font_family, fontface = "italic", show.legend = FALSE) +
    ggplot2::annotate("text", x = -0.5 * left_span, y = max(d$.y) + 0.9,
      label = group_use[1], colour = palcolor[group_use[1]], fontface = "bold", size = 4.7) +
    ggplot2::annotate("text", x = 0.5 * right_span, y = max(d$.y) + 0.9,
      label = group_use[2], colour = palcolor[group_use[2]], fontface = "bold", size = 4.7) +
    ggplot2::scale_fill_manual(values = palcolor, guide = "none") +
    ggplot2::scale_colour_manual(values = palcolor, guide = "none") +
    ggplot2::scale_alpha_continuous(name = "-log10(adjusted P)", range = alpha_range,
      limits = range(d$.score), n.breaks = 3,
      guide = ggplot2::guide_legend(override.aes = list(fill = "grey30", colour = "grey30"))) +
    ggplot2::scale_x_continuous(limits = c(band_left - 0.15, right_span * 1.06),
      breaks = pretty(c(-left_span, right_span), n = 7), labels = function(x) abs(x),
      expand = ggplot2::expansion(mult = 0)) +
    ggplot2::scale_y_continuous(breaks = NULL,
      expand = ggplot2::expansion(add = c(0.35, 0.15))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = title, x = "-log10(adjusted P)", y = NULL,
      subtitle = paste0("Top ", topTerm, " terms per database and group; adjusted P < ", padjustCutoff),
      caption = "Independent group-wise rankings; opposite rows need not represent the same pathway.\nLeftward bars do not indicate negative enrichment or downregulation.") +
    ggplot2::theme_classic(base_size = base_size, base_family = font_family) +
    ggplot2::theme(axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(linewidth = 0.35),
      plot.title = ggplot2::element_text(size = 16),
      plot.subtitle = ggplot2::element_text(size = 10, colour = "grey35"),
      plot.caption = ggplot2::element_text(size = 8, hjust = 0, colour = "grey40"),
      legend.position = "bottom", legend.direction = "horizontal",
      plot.margin = ggplot2::margin(10, 20, 8, 20))
  p
}
