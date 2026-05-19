#' @title Differential Expression Test Plot
#'
#' @md
#' @inheritParams CellDimPlot
#' @param srt A `Seurat` object or `SummarizedExperiment` object containing the
#' results of differential expression analysis.
#' @param res A `data.frame` or `data.table` with differential expression results.
#' When `res` is provided, `srt` will be ignored.
#' The data.frame must contain columns: `gene`, `group1` (factor or character),
#' `avg_log2FC`, `p_val_adj`, and optionally `pct.1` and `pct.2` for calculating `diff_pct`.
#' @param test.use A character string specifying the type of statistical test to use.
#' Default is `"wilcox"`.
#' @param plot_type Type of plot to create. Options: `"volcano"`, `"manhattan"`, or `"ring"`.
#' Default is `"volcano"`.
#' @param group_use Groups to plot.
#' Default is `NULL` (all groups).
#' @param DE_threshold A character string specifying the threshold for differential expression (used to highlight significant genes in all plot types).
#' Default is `"p_val < 0.05"` for sample-level methods (`"edgeR"` and `"limma"`)
#' and `"avg_log2FC > 0 & p_val_adj < 0.05"` otherwise.
#' @param x_metric A character string specifying the metric to use for the x-axis (only for volcano plot).
#' Default is `NULL`, which uses `"avg_log2FC"` for sample-level methods (`"edgeR"` and `"limma"`)
#' and `"diff_pct"` otherwise.
#' @param threshold_method Volcano significance threshold method.
#' Options are `"rectangular"` (legacy DE_threshold) or `"hyperbolic"` (`|log2FC * -log10(padj)| > c`).
#' Default is `"rectangular"`.
#' @param hyperbola_c Numeric cutoff `c` for hyperbolic volcano threshold.
#' Default is `6`.
#' @param annotate_enrichment Whether to annotate enrichment-hit genes on volcano plots.
#' Enrichment results are read from existing results in `srt@tools` only.
#' Default is `FALSE`.
#' @param enrich_from Character vector specifying enrichment result source(s) to annotate.
#' Options are `"Enrichment"`, `"GSEA"`, `"GSVA"`.
#' Default is `c("Enrichment", "GSEA", "GSVA")`.
#' @param enrich_db Optional database filter for enrichment annotation, e.g. `"GO_BP"` or `"KEGG"`.
#' Default is `NULL`.
#' @param enrich_terms Optional whitelist of enrichment term IDs or names for annotation.
#' Default is `NULL`.
#' @param enrich_top_terms Number of top enriched terms selected per source/group/database.
#' Default is `3`.
#' @param enrich_padj_cutoff Adjusted p-value cutoff for `"Enrichment"` and `"GSEA"` annotation.
#' Default is `0.05`.
#' @param enrich_gsva_score_cutoff Optional absolute GSVA score cutoff for `"GSVA"` annotation.
#' Default is `NULL`.
#' @param gsva_method Optional GSVA method filter (e.g. `"gsva"` or `"ssgsea"`) when multiple GSVA tool slots exist.
#' Default is `NULL`.
#' @param enrich_nlabel Maximum number of enrichment-derived labels added per group.
#' Labels from `features_label` are always retained.
#' Default is `15`.
#' @param y_metric A character string specifying the metric to use for the y-axis.
#' Options: `"p_val"` or `"p_val_adj"`.
#' Default is `"p_val"` for sample-level methods (`"edgeR"` and `"limma"`) and `"p_val_adj"` otherwise.
#' @param x_order A character string specifying how to order genes on x-axis (only for Manhattan plot, not used currently).
#' Options: `"gene"` (alphabetical by gene name) or `"index"` (by data order). Default is `"gene"`.
#' @param palette Color palette name.
#' Available palettes can be found in [thisplot::show_palettes].
#' Default is `"RdBu"`.
#' @param group_palette Palette for cell types (groups) in Manhattan plot.
#' Default is `"Chinese"`.
#' @param group_palcolor Custom colors for cell types (groups) in Manhattan plot.
#' Default is `NULL`.
#' @param pt.size The size of the points.
#' Default is `1`.
#' @param cols.background A character string specifying the color for non-DE background points in volcano plots.
#' Default is `"grey80"`.
#' @param cols.highlight A character string specifying the color for highlighted points.
#' Default is `"black"`.
#' @param sizes.highlight The size of the highlighted points.
#' Default is `1`.
#' @param alpha.highlight The transparency of the highlighted points.
#' Default is `1`.
#' @param stroke.highlight The stroke width for the highlighted points.
#' Default is `0.5`.
#' @param nlabel An integer value specifying the number of labeled points per group.
#' Default is `5`.
#' @param features_label A character vector specifying the feature labels to plot.
#' Default is `NULL`.
#' @param label.fg A character string specifying the color for the labels' foreground.
#' Default is `"black"`.
#' @param label.bg A character string specifying the color for the labels' background.
#' Default is `"white"`.
#' @param label.bg.r The radius of the rounding of the labels' background.
#' Default is `0.1`.
#' @param label.size The size of the labels.
#' Default is `4`.
#' @param aspect.ratio Aspect ratio of the panel.
#' Default is `NULL`.
#' @param xlab A character string specifying the x-axis label.
#' @param ylab A character string specifying the y-axis label.
#' @param theme_use Theme to use for the plot.
#' Default is `"theme_scop"`.
#' @param theme_args A list of additional arguments to pass to the theme function.
#' Default is `list()`.
#' @param combine Whether to combine multiple plots into one.
#' Default is `TRUE`.
#' @param nrow Number of rows for combined plots.
#' Default is `NULL`.
#' @param ncol Number of columns for combined plots.
#' Default is `NULL`.
#' @param byrow Whether to fill plots by row.
#' Default is `TRUE`.
#' @param manhattan.bg Background color for Manhattan plot.
#' Default is `"white"`.
#' @param jitter_width Horizontal jitter range for points in Manhattan plot.
#' Default is `0.5`.
#' @param jitter_height Vertical jitter range for points in Manhattan plot.
#' Default is `0.4`.
#' @param tile_height Height of the cell-type track in ring plot.
#' Default is `0.3`.
#' @param tile_gap Gap between the track and nudged points in ring plot.
#' Default is `0.1`.
#' @param ring_segments Whether to draw segment lines between cell types in ring plot.
#' Default is `TRUE`.
#' @param seed Random seed for jitter in ring plot.
#' Default is `11`.
#'
#' @seealso [RunDEtest], [VolcanoPlot], [DEtestManhattanPlot], [DEtestRingPlot]
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
#' DEtestPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "volcano",
#'   ncol = 2
#' )
#'
#' DEtestPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "volcano",
#'   group_use = c("Ductal", "Endocrine"),
#'   ncol = 2
#' )
#' 
#' DEtestPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   test.use = "wilcox",
#'   plot_type = "volcano",
#'   x_metric = "avg_log2FC",
#'   y_metric = "p_val_adj",
#'   group_use = c("Ductal", "Endocrine"),
#'   DE_threshold = "abs(avg_log2FC) > 0.25 & p_val_adj < 0.05"
#' )
#'
#' DEtestPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "volcano",
#'   threshold_method = "hyperbolic",
#'   hyperbola_c = 6,
#'   ncol = 2
#' )
#'
#' pancreas_sub <- RunEnrichment(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   db = "GO_BP",
#'   species = "Mus_musculus"
#' )
#' DEtestPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "volcano",
#'   threshold_method = "hyperbolic",
#'   hyperbola_c = 6,
#'   annotate_enrichment = TRUE,
#'   enrich_from = "Enrichment",
#'   enrich_db = "GO_BP",
#'   enrich_top_terms = 3,
#'   enrich_nlabel = 15,
#'   ncol = 2
#' )
#'
#' DEtestPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "manhattan"
#' )
#'
#' DEtestPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "ring"
#' )
#'
#' de_results1 <- pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox
#' DEtestPlot(
#'   res = de_results1,
#'   plot_type = "volcano",
#'   ncol = 2
#' )
#'
#' de_results2 <- Seurat::FindMarkers(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   ident.1 = "Ductal",
#'   ident.2 = "Endocrine"
#' )
#' DEtestPlot(
#'   res = de_results2,
#'   plot_type = "volcano"
#' )
#'
#' de_results3 <- Seurat::FindAllMarkers(
#'   pancreas_sub,
#'   group.by = "CellType"
#' )
#' DEtestPlot(
#'   res = de_results3,
#'   plot_type = "volcano",
#'   ncol = 2
#' )
DEtestPlot <- function(
    srt,
    group.by = NULL,
    test.use = "wilcox",
    res = NULL,
    plot_type = c("volcano", "manhattan", "ring"),
    group_use = NULL,
    DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
    x_metric = NULL,
    y_metric = c("p_val_adj", "p_val"),
    x_order = c("gene", "index"),
    palette = "RdBu",
    palcolor = NULL,
    group_palette = "Chinese",
    group_palcolor = NULL,
    pt.size = 1,
    pt.alpha = 1,
    cols.background = "grey80",
    cols.highlight = "black",
    sizes.highlight = 1,
    alpha.highlight = 1,
    stroke.highlight = 0.5,
    nlabel = 5,
    features_label = NULL,
    label.fg = "black",
    label.bg = "white",
    label.bg.r = 0.1,
    label.size = 4,
    aspect.ratio = NULL,
    xlab = NULL,
    ylab = NULL,
    theme_use = "theme_scop",
    theme_args = list(),
    combine = TRUE,
    nrow = NULL,
    ncol = NULL,
    byrow = TRUE,
    manhattan.bg = "white",
    jitter_width = 0.5,
    jitter_height = 0.4,
    tile_height = 0.3,
    tile_gap = 0.1,
    ring_segments = TRUE,
    seed = 11,
    threshold_method = c("rectangular", "hyperbolic"),
    hyperbola_c = 6,
    annotate_enrichment = FALSE,
    enrich_from = c("Enrichment", "GSEA", "GSVA"),
    enrich_db = NULL,
    enrich_terms = NULL,
    enrich_top_terms = 3,
    enrich_padj_cutoff = 0.05,
    enrich_gsva_score_cutoff = NULL,
    gsva_method = NULL,
    enrich_nlabel = 15) {
  DE_threshold_missing <- missing(DE_threshold)
  y_metric_missing <- missing(y_metric)
  plot_type <- match.arg(plot_type)
  threshold_method <- match.arg(threshold_method)
  y_metric <- match.arg(y_metric)
  x_order <- match.arg(x_order)

  if (plot_type == "volcano") {
    DE_threshold_use <- if (
      DE_threshold_missing &&
        test.use %in% c("edgeR", "limma")
    ) {
      "p_val < 0.05"
    } else {
      DE_threshold
    }
    x_metric_use <- x_metric %||% ifelse(
      test.use %in% c("edgeR", "limma") || threshold_method == "hyperbolic",
      "avg_log2FC",
      "diff_pct"
    )
    y_metric_use <- if (y_metric_missing) {
      ifelse(test.use %in% c("edgeR", "limma"), "p_val", "p_val_adj")
    } else {
      y_metric
    }
    return(
      VolcanoPlot(
        srt = srt,
        group.by = group.by,
        test.use = test.use,
        res = res,
        group_use = group_use,
        DE_threshold = DE_threshold_use,
        x_metric = x_metric_use,
        y_metric = y_metric_use,
        threshold_method = threshold_method,
        hyperbola_c = hyperbola_c,
        annotate_enrichment = annotate_enrichment,
        enrich_from = enrich_from,
        enrich_db = enrich_db,
        enrich_terms = enrich_terms,
        enrich_top_terms = enrich_top_terms,
        enrich_padj_cutoff = enrich_padj_cutoff,
        enrich_gsva_score_cutoff = enrich_gsva_score_cutoff,
        gsva_method = gsva_method,
        enrich_nlabel = enrich_nlabel,
        palette = palette,
        palcolor = palcolor,
        pt.size = pt.size,
        pt.alpha = pt.alpha,
        cols.background = cols.background,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        alpha.highlight = alpha.highlight,
        stroke.highlight = stroke.highlight,
        nlabel = nlabel,
        features_label = features_label,
        label.fg = label.fg,
        label.bg = label.bg,
        label.bg.r = label.bg.r,
        label.size = label.size,
        aspect.ratio = aspect.ratio,
        xlab = xlab %||% ifelse(
          threshold_method == "hyperbolic",
          "avg_log2FC",
          x_metric_use
        ),
        ylab = ylab %||% ifelse(
          y_metric_use == "p_val",
          "-log10(p-value)",
          "-log10(p-adjust)"
        ),
        theme_use = theme_use,
        theme_args = theme_args,
        combine = combine,
        nrow = nrow,
        ncol = ncol,
        byrow = byrow
      )
    )
  }

  if (plot_type == "manhattan") {
    return(
      DEtestManhattanPlot(
        srt = srt,
        group.by = group.by,
        test.use = test.use,
        res = res,
        group_use = group_use,
        DE_threshold = DE_threshold,
        group_palette = group_palette,
        group_palcolor = group_palcolor,
        pt.size = pt.size,
        pt.alpha = pt.alpha,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        alpha.highlight = alpha.highlight,
        stroke.highlight = stroke.highlight,
        nlabel = nlabel,
        features_label = features_label,
        label.fg = label.fg,
        label.bg = label.bg,
        label.bg.r = label.bg.r,
        label.size = label.size,
        palette = palette,
        palcolor = palcolor,
        theme_use = theme_use,
        theme_args = theme_args,
        manhattan.bg = manhattan.bg,
        jitter_width = jitter_width,
        jitter_height = jitter_height,
        aspect.ratio = aspect.ratio,
        xlab = xlab,
        ylab = ylab
      )
    )
  }
  if (plot_type == "ring") {
    return(
      DEtestRingPlot(
        srt = srt,
        group.by = group.by,
        test.use = test.use,
        res = res,
        group_use = group_use,
        DE_threshold = DE_threshold,
        group_palette = group_palette,
        group_palcolor = group_palcolor,
        pt.size = pt.size,
        pt.alpha = pt.alpha,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        alpha.highlight = alpha.highlight,
        stroke.highlight = stroke.highlight,
        nlabel = nlabel,
        features_label = features_label,
        label.fg = label.fg,
        label.bg = label.bg,
        label.bg.r = label.bg.r,
        label.size = label.size,
        palette = palette,
        palcolor = palcolor,
        theme_use = theme_use,
        theme_args = theme_args,
        tile_height = tile_height,
        tile_gap = tile_gap,
        jitter_width = jitter_width,
        ring_segments = ring_segments,
        seed = seed
      )
    )
  }
  invisible(NULL)
}

get_de_data <- function(
  srt,
  group.by,
  test.use,
  DE_threshold,
  res = NULL,
  group_use = NULL
) {
  de_df <- resolve_detest_result(
    object = srt,
    group.by = group.by,
    test.use = test.use,
    res = res
  )
  if (!"group1" %in% colnames(de_df)) {
    if ("cluster" %in% colnames(de_df)) {
      de_df[, "group1"] <- de_df[, "cluster"]
    } else {
      de_df[, "group1"] <- "All"
    }
  }
  if ("comparison" %in% colnames(de_df)) {
    comparison <- as.character(de_df[["comparison"]])
    has_comparison <- !is.na(comparison) & nzchar(comparison)
    if (any(has_comparison)) {
      de_df[has_comparison, "group1"] <- comparison[has_comparison]
    }
  } else if ("group2" %in% colnames(de_df)) {
    g1_levels <- unique(as.character(de_df[["group1"]]))
    g1_levels <- g1_levels[!is.na(g1_levels) & nzchar(g1_levels)]
    g2_levels <- unique(as.character(de_df[["group2"]]))
    g2_levels <- g2_levels[!is.na(g2_levels) & nzchar(g2_levels)]
    if (length(g1_levels) <= 1 && length(g2_levels) >= 1) {
      de_df[, "group1"] <- as.character(de_df[["group2"]])
    } else if (length(g2_levels) > 1) {
      g1 <- as.character(de_df[["group1"]])
      g2 <- as.character(de_df[["group2"]])
      has_g2 <- !is.na(g2) & nzchar(g2)
      has_g1 <- !is.na(g1) & nzchar(g1)
      de_df[has_g2, "group1"] <- ifelse(
        has_g1[has_g2],
        paste0(g1[has_g2], "::", g2[has_g2]),
        g2[has_g2]
      )
    }
  }
  if (!is.factor(de_df[["group1"]])) {
    de_df[["group1"]] <- factor(
      de_df[["group1"]],
      levels = unique(de_df[["group1"]])
    )
  }
  de_df <- filter_de_data_group_use(de_df, group_use = group_use)
  if ("pct.1" %in% colnames(de_df) && "pct.2" %in% colnames(de_df)) {
    de_df[, "diff_pct"] <- de_df[, "pct.1"] - de_df[, "pct.2"]
  } else {
    de_df[, "diff_pct"] <- 0
  }
  de_df[, "-log10padj"] <- -log10(de_df[, "p_val_adj"])
  de_df[, "DE"] <- FALSE
  de_df[with(de_df, eval(rlang::parse_expr(DE_threshold))), "DE"] <- TRUE
  list(de_df = de_df)
}

filter_de_data_group_use <- function(de_df, group_use = NULL) {
  if (is.null(group_use)) {
    return(de_df)
  }

  available_groups <- unique(as.character(de_df[["group1"]]))
  group_levels <- if (is.factor(de_df[["group1"]])) {
    levels(de_df[["group1"]])
  } else {
    available_groups
  }
  group_levels <- group_levels[group_levels %in% available_groups]
  group_use <- unique(as.character(group_use))
  group_use <- group_use[nzchar(group_use)]
  group_use <- intersect(group_use, group_levels)
  if (length(group_use) == 0) {
    log_message(
      "No matching groups found in DEtest results",
      message_type = "error"
    )
  }

  de_df <- de_df[
    as.character(de_df[["group1"]]) %in% group_use, ,
    drop = FALSE
  ]
  de_df[["group1"]] <- factor(as.character(de_df[["group1"]]), levels = group_use)
  de_df
}

get_bulk_de_plot_data <- function(srt) {
  resolve_detest_result(object = srt)
}

clip_log2fc_symmetric <- function(df, fc_col = "avg_log2FC") {
  res <- thisplot::clip_symmetric_range(data = df, value_col = fc_col)
  list(df = res$data, fc_lim = res$limits)
}

filter_de_markers <- function(de_df, log2FC_cutoff, pvalue_cutoff) {
  de_df_marker <- de_df[
    abs(de_df[, "avg_log2FC"]) >= log2FC_cutoff &
      de_df[, "p_val"] < pvalue_cutoff, ,
    drop = FALSE
  ]
  if (nrow(de_df_marker) == 0) {
    log_message(
      "No genes pass the threshold. Please adjust the threshold.",
      message_type = "warning"
    )
    return(NULL)
  }
  de_df_marker
}

get_top_markers_for_label <- function(de_df_marker, cluster_levels, nlabel, features_label) {
  top_marker_list <- lapply(cluster_levels, function(x) {
    tmp <- de_df_marker[de_df_marker[["group1"]] == x, , drop = FALSE]
    if (nrow(tmp) == 0) {
      return(NULL)
    }
    if (is.null(features_label)) {
      top_max <- if (nrow(tmp) > 0) {
        tmp[utils::head(order(tmp[, "avg_log2FC"], decreasing = TRUE), nlabel), , drop = FALSE]
      } else {
        tmp[0, , drop = FALSE]
      }
      top_min <- if (nrow(tmp) > 0) {
        tmp[utils::head(order(tmp[, "avg_log2FC"], decreasing = FALSE), nlabel), , drop = FALSE]
      } else {
        tmp[0, , drop = FALSE]
      }
      rbind(top_max, top_min)
    } else {
      tmp[tmp[, "gene"] %in% features_label, , drop = FALSE]
    }
  })
  do.call(rbind, top_marker_list[!sapply(top_marker_list, is.null)])
}

add_volcano_plot_coords <- function(df, jitter_width = 0.2, jitter_height = 0.2, seed = 11) {
  df[, "x_plot"] <- df[, "x"]
  df[, "y_plot"] <- df[, "y"]

  border_idx <- which(df[, "border"] & is.finite(df[, "x"]) & is.finite(df[, "y"]))
  if (length(border_idx) == 0) {
    return(df)
  }

  idx <- seq_along(border_idx)
  x_offset <- ((((idx * 0.61803398875) + (seed * 0.01)) %% 1) - 0.5) * 2 * jitter_width
  y_offset <- ((((idx * 0.41421356237) + (seed * 0.01)) %% 1) - 0.5) * 2 * jitter_height
  df[border_idx, "x_plot"] <- df[border_idx, "x"] + x_offset
  df[border_idx, "y_plot"] <- df[border_idx, "y"] + y_offset

  df
}

get_safe_neglog10 <- function(padj, p_floor = .Machine$double.xmin) {
  p_floor <- max(as.numeric(p_floor), .Machine$double.xmin)
  padj_num <- suppressWarnings(as.numeric(padj))
  padj_num[is.na(padj_num)] <- NA_real_
  padj_num[padj_num <= 0 & !is.na(padj_num)] <- p_floor
  padj_num[padj_num > 1] <- 1
  -log10(padj_num)
}

compute_hyperbolic_de_flags <- function(
    avg_log2FC,
    p_val_adj,
    hyperbola_c = 6,
    p_floor = .Machine$double.xmin) {
  neglog10 <- get_safe_neglog10(padj = p_val_adj, p_floor = p_floor)
  prod_value <- abs(as.numeric(avg_log2FC)) * neglog10
  is.finite(prod_value) & (prod_value > as.numeric(hyperbola_c))
}

build_hyperbola_curve_df <- function(
    x_range,
    hyperbola_c = 6,
    y_max = NULL,
    n = 400,
    x_eps = NULL) {
  if (length(x_range) != 2 || any(!is.finite(x_range))) {
    return(data.frame(x = numeric(0), y = numeric(0)))
  }
  xmin <- min(x_range)
  xmax <- max(x_range)
  if (xmin >= xmax) {
    return(data.frame(x = numeric(0), y = numeric(0)))
  }
  x_eps <- x_eps %||% max((xmax - xmin) / n, 1e-4)
  x_eps <- max(abs(x_eps), 1e-8)
  x_left <- if (xmin < -x_eps) {
    seq(xmin, -x_eps, length.out = max(10L, floor(n / 2)))
  } else {
    numeric(0)
  }
  x_right <- if (xmax > x_eps) {
    seq(x_eps, xmax, length.out = max(10L, floor(n / 2)))
  } else {
    numeric(0)
  }
  x <- c(x_left, x_right)
  if (length(x) == 0) {
    return(data.frame(x = numeric(0), y = numeric(0)))
  }
  y <- as.numeric(hyperbola_c) / abs(x)
  if (!is.null(y_max) && is.finite(y_max)) {
    y[y > y_max] <- NA_real_
  }
  data.frame(x = x, y = y)
}

split_enrichment_genes <- function(gene_text) {
  if (length(gene_text) == 0 || is.na(gene_text) || !nzchar(trimws(gene_text))) {
    return(character(0))
  }
  genes <- unlist(strsplit(gene_text, "[/;,]"))
  genes <- trimws(genes)
  unique(genes[nzchar(genes)])
}

get_enrichment_overlay_colors <- function(keys) {
  keys <- unique(as.character(keys))
  if (length(keys) == 0) {
    return(character(0))
  }
  # Dedicated pastel categorical colors for enrichment overlays.
  # Keep away from the strong RdBu-like contrast used by volcano points.
  base_cols <- c(
    "#CFE8CC", # mint
    "#E9D8A6", # sand
    "#D9C2F0", # lavender
    "#FAD4C0", # peach
    "#C7EAE4", # aqua
    "#E6D5B8", # beige
    "#F4C6D7", # blush
    "#D5E6F2", # ice blue
    "#E3F0CC", # light lime
    "#F2E2CE", # apricot cream
    "#DCCFE6", # soft mauve
    "#CCE3D9"  # sage
  )
  cols <- if (length(keys) <= length(base_cols)) {
    base_cols[seq_along(keys)]
  } else {
    grDevices::colorRampPalette(base_cols)(length(keys))
  }
  names(cols) <- keys
  cols
}

get_gsva_tool_names <- function(srt, group.by = NULL, gsva_method = NULL) {
  tool_names <- grep("^GSVA_", names(srt@tools), value = TRUE)
  if (length(tool_names) == 0) {
    return(character(0))
  }
  if (!is.null(gsva_method)) {
    tool_names <- tool_names[grepl(
      pattern = paste0("_", gsva_method, "$"),
      x = tool_names
    )]
  }
  if (!is.null(group.by)) {
    pattern1 <- paste0("^GSVA_", group.by, "_")
    pattern2 <- paste0("^GSVA_.*_", group.by, "$")
    tool_names <- tool_names[grepl(pattern1, tool_names) | grepl(pattern2, tool_names)]
  }
  unique(tool_names)
}

extract_volcano_enrichment_source <- function(
    srt,
    source = c("Enrichment", "GSEA", "GSVA"),
    group.by = NULL,
    test.use = "wilcox",
    gsva_method = NULL) {
  source <- match.arg(source)
  group.by <- group.by %||% "custom"
  output <- list()

  if (source %in% c("Enrichment", "GSEA")) {
    layer <- paste(source, group.by, test.use, sep = "_")
    if (!layer %in% names(srt@tools)) {
      return(data.frame())
    }
    enrichment <- srt@tools[[layer]][["enrichment"]]
    if (is.null(enrichment) || nrow(enrichment) == 0) {
      return(data.frame())
    }
    enrichment <- as.data.frame(enrichment)
    gene_col <- if (source == "GSEA" && "core_enrichment" %in% colnames(enrichment)) {
      "core_enrichment"
    } else if ("geneID" %in% colnames(enrichment)) {
      "geneID"
    } else {
      NULL
    }
    if (is.null(gene_col)) {
      return(data.frame())
    }
    metric_col <- if ("p.adjust" %in% colnames(enrichment)) {
      "p.adjust"
    } else if ("pvalue" %in% colnames(enrichment)) {
      "pvalue"
    } else {
      NULL
    }
    if (is.null(metric_col)) {
      return(data.frame())
    }
    df <- data.frame(
      source = source,
      Groups = as.character(enrichment[["Groups"]]),
      Database = as.character(enrichment[["Database"]]),
      term_id = as.character(enrichment[["ID"]] %||% enrichment[["Description"]]),
      term_name = as.character(enrichment[["Description"]] %||% enrichment[["ID"]]),
      gene_text = as.character(enrichment[[gene_col]]),
      metric_value = suppressWarnings(as.numeric(enrichment[[metric_col]])),
      stringsAsFactors = FALSE
    )
    return(df)
  }

  tool_names <- get_gsva_tool_names(
    srt = srt,
    group.by = group.by,
    gsva_method = gsva_method
  )
  if (length(tool_names) == 0) {
    return(data.frame())
  }
  for (tool_name in tool_names) {
    enrichment <- srt@tools[[tool_name]][["enrichment"]]
    if (is.null(enrichment) || nrow(enrichment) == 0) {
      next
    }
    enrichment <- as.data.frame(enrichment)
    if (!"geneID" %in% colnames(enrichment)) {
      next
    }
    if (!"GSVA_Score" %in% colnames(enrichment)) {
      next
    }
    default_db <- srt@tools[[tool_name]][["db"]]
    default_db <- default_db %||% "GSVA"
    if (length(default_db) > 1) {
      default_db <- default_db[1]
    }
    db_col <- enrichment[["Database"]] %||% rep(default_db, nrow(enrichment))
    df <- data.frame(
      source = "GSVA",
      Groups = as.character(enrichment[["Groups"]]),
      Database = as.character(db_col),
      term_id = as.character(enrichment[["ID"]] %||% enrichment[["Description"]]),
      term_name = as.character(enrichment[["Description"]] %||% enrichment[["ID"]]),
      gene_text = as.character(enrichment[["geneID"]]),
      metric_value = suppressWarnings(as.numeric(enrichment[["GSVA_Score"]])),
      stringsAsFactors = FALSE
    )
    output[[tool_name]] <- df
  }
  if (length(output) == 0) {
    return(data.frame())
  }
  do.call(rbind, output)
}

collect_volcano_enrichment_annotations <- function(
    srt,
    group.by = NULL,
    test.use = "wilcox",
    enrich_from = c("Enrichment", "GSEA", "GSVA"),
    enrich_db = NULL,
    enrich_terms = NULL,
    enrich_top_terms = 3,
    enrich_padj_cutoff = 0.05,
    enrich_gsva_score_cutoff = NULL,
    gsva_method = NULL) {
  if (is.null(srt) || !inherits(srt, "Seurat")) {
    return(data.frame())
  }
  enrich_from <- unique(intersect(
    enrich_from,
    c("Enrichment", "GSEA", "GSVA")
  ))
  if (length(enrich_from) == 0) {
    return(data.frame())
  }

  res_list <- lapply(enrich_from, function(src) {
    extract_volcano_enrichment_source(
      srt = srt,
      source = src,
      group.by = group.by,
      test.use = test.use,
      gsva_method = gsva_method
    )
  })
  res_list <- res_list[sapply(res_list, nrow) > 0]
  if (length(res_list) == 0) {
    return(data.frame())
  }
  enrichment <- do.call(rbind, res_list)
  rownames(enrichment) <- NULL

  if (!is.null(enrich_db)) {
    enrichment <- enrichment[
      enrichment[["Database"]] %in% enrich_db, ,
      drop = FALSE
    ]
  }
  if (!is.null(enrich_terms) && nrow(enrichment) > 0) {
    terms_lower <- tolower(enrich_terms)
    id_match <- tolower(as.character(enrichment[["term_id"]])) %in% terms_lower
    name_match <- tolower(as.character(enrichment[["term_name"]])) %in% terms_lower
    enrichment <- enrichment[id_match | name_match, , drop = FALSE]
  }
  if (nrow(enrichment) == 0) {
    return(data.frame())
  }

  enrich_top_terms <- as.integer(enrich_top_terms)
  if (is.na(enrich_top_terms) || enrich_top_terms <= 0) {
    enrich_top_terms <- 1L
  }

  selected_terms <- list()
  split_key <- paste(
    enrichment[["source"]],
    enrichment[["Groups"]],
    enrichment[["Database"]],
    sep = "|||"
  )
  split_df <- split(enrichment, split_key)
  for (key in names(split_df)) {
    df <- split_df[[key]]
    if (nrow(df) == 0) {
      next
    }
    if (df[["source"]][1] == "GSVA") {
      metric <- abs(df[["metric_value"]])
      if (!is.null(enrich_gsva_score_cutoff)) {
        df <- df[metric >= abs(enrich_gsva_score_cutoff), , drop = FALSE]
        metric <- abs(df[["metric_value"]])
      }
      if (nrow(df) == 0) {
        next
      }
      ord <- order(metric, decreasing = TRUE, na.last = NA)
    } else {
      metric <- df[["metric_value"]]
      df <- df[is.finite(metric), , drop = FALSE]
      metric <- df[["metric_value"]]
      df <- df[metric <= enrich_padj_cutoff, , drop = FALSE]
      if (nrow(df) == 0) {
        next
      }
      ord <- order(metric, decreasing = FALSE, na.last = NA)
    }
    df <- df[ord, , drop = FALSE]
    df <- df[!duplicated(df[["term_id"]]), , drop = FALSE]
    selected_terms[[key]] <- utils::head(df, enrich_top_terms)
  }
  if (length(selected_terms) == 0) {
    return(data.frame())
  }
  selected_terms <- do.call(rbind, selected_terms)
  rownames(selected_terms) <- NULL

  gene_map <- lapply(seq_len(nrow(selected_terms)), function(i) {
    term_row <- selected_terms[i, , drop = FALSE]
    genes <- split_enrichment_genes(term_row[["gene_text"]])
    if (length(genes) == 0) {
      return(NULL)
    }
    data.frame(
      source = term_row[["source"]],
      group1 = term_row[["Groups"]],
      Database = term_row[["Database"]],
      term_id = term_row[["term_id"]],
      term_name = term_row[["term_name"]],
      metric_value = term_row[["metric_value"]],
      gene = genes,
      stringsAsFactors = FALSE
    )
  })
  gene_map <- gene_map[!sapply(gene_map, is.null)]
  if (length(gene_map) == 0) {
    return(data.frame())
  }
  gene_map <- do.call(rbind, gene_map)
  rownames(gene_map) <- NULL

  source_priority <- c("Enrichment" = 1L, "GSEA" = 2L, "GSVA" = 3L)
  gene_map[["source_priority"]] <- source_priority[gene_map[["source"]]]
  gene_map[["source_priority"]][is.na(gene_map[["source_priority"]])] <- 99L
  gene_map[["metric_order"]] <- ifelse(
    gene_map[["source"]] == "GSVA",
    -abs(gene_map[["metric_value"]]),
    gene_map[["metric_value"]]
  )
  gene_map <- gene_map[order(
    gene_map[["group1"]],
    gene_map[["gene"]],
    gene_map[["source_priority"]],
    gene_map[["metric_order"]],
    na.last = TRUE
  ), , drop = FALSE]

  gene_map <- gene_map[
    !duplicated(paste(gene_map[["group1"]], gene_map[["gene"]], sep = "|||")),
    ,
    drop = FALSE
  ]
  rownames(gene_map) <- NULL
  gene_map
}

#' @title DEtest Manhattan Plot
#'
#' @description
#' Draw a Manhattan-style plot of differential expression results by cell type.
#'
#' @md
#' @inheritParams DEtestPlot
#'
#' @seealso [DEtestPlot], [RunDEtest], [VolcanoPlot], [DEtestRingPlot]
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
#' DEtestManhattanPlot(
#'   pancreas_sub,
#'   group.by = "CellType"
#' )
DEtestManhattanPlot <- function(
  srt,
  group.by = NULL,
  test.use = "wilcox",
  res = NULL,
  group_use = NULL,
  DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
  group_palette = "Chinese",
  group_palcolor = NULL,
  pt.size = 1,
  pt.alpha = 1,
  cols.highlight = "black",
  sizes.highlight = 1,
  alpha.highlight = 1,
  stroke.highlight = 0.5,
  nlabel = 5,
  features_label = NULL,
  label.fg = "black",
  label.bg = "white",
  label.bg.r = 0.1,
  label.size = 4,
  palette = "RdBu",
  palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  manhattan.bg = "white",
  jitter_width = 0.5,
  jitter_height = 0.4,
  aspect.ratio = NULL,
  xlab = NULL,
  ylab = NULL
) {
  if (is.null(group.by)) {
    group.by <- "custom"
  }
  data_res <- get_de_data(srt, group.by, test.use, DE_threshold, res, group_use)
  de_df <- data_res$de_df
  clip_res <- clip_log2fc_symmetric(de_df)
  de_df_marker <- clip_res$df
  fc_lim <- clip_res$fc_lim
  de_df_marker[, "type"] <- ifelse(
    de_df_marker[, "avg_log2FC"] >= 0,
    "sigUp",
    "sigDown"
  )
  cluster_levels <- levels(de_df_marker[["group1"]])
  show_group_track <- length(cluster_levels) > 1
  n_m <- nrow(de_df_marker)
  de_df_marker[, "x_num"] <- as.numeric(de_df_marker[, "group1"])
  de_df_marker[, "x_plot"] <- de_df_marker[, "x_num"] + (stats::runif(n_m) - 0.5) * jitter_width
  de_df_marker[, "y_plot"] <- de_df_marker[, "avg_log2FC"] + (stats::runif(n_m) - 0.5) * jitter_height
  top_marker <- get_top_markers_for_label(de_df_marker, cluster_levels, nlabel, features_label)
  back_data_list <- lapply(cluster_levels, function(x) {
    tmp <- de_df_marker[de_df_marker[["group1"]] == x, , drop = FALSE]
    if (nrow(tmp) == 0) {
      return(NULL)
    }
    data.frame(
      cluster = x,
      min = min(tmp[, "avg_log2FC"], na.rm = TRUE) - 0.2,
      max = max(tmp[, "avg_log2FC"], na.rm = TRUE) + 0.2
    )
  })
  back_data <- do.call(rbind, back_data_list[!sapply(back_data_list, is.null)])
  back_data[, "x_num"] <- match(back_data[, "cluster"], cluster_levels)
  tile_colors <- palette_colors(
    x = cluster_levels,
    palette = group_palette,
    palcolor = group_palcolor
  )
  tile_data <- data.frame(
    group1 = cluster_levels,
    x = seq_along(cluster_levels),
    y = 0
  )
  p1 <- ggplot(de_df_marker, aes(x = x_plot, y = y_plot)) +
    geom_col(
      data = back_data,
      aes(x = x_num, y = min),
      fill = manhattan.bg,
      inherit.aes = FALSE
    ) +
    geom_col(
      data = back_data,
      aes(x = x_num, y = max),
      fill = manhattan.bg,
      inherit.aes = FALSE
    )
  p2 <- p1 +
    geom_point(
      aes(x = x_plot, y = y_plot),
      color = cols.highlight,
      size = sizes.highlight + stroke.highlight,
      alpha = alpha.highlight,
      inherit.aes = FALSE,
      data = de_df_marker
    ) +
    geom_point(aes(color = avg_log2FC), size = pt.size, alpha = pt.alpha) +
    scale_color_gradientn(
      name = "log2FC",
      colors = palette_colors(palette = palette, palcolor = palcolor),
      values = scales::rescale(unique(c(fc_lim[1], 0, fc_lim[2]))),
      limits = fc_lim,
      guide = guide_colorbar(
        frame.colour = "black",
        ticks.colour = "black",
        title.hjust = 0,
        order = 1
      )
    )
  xlab_use <- if (!is.null(res) && group.by == "custom") NULL else (xlab %||% group.by)
  p3 <- p2 +
    scale_x_continuous(
      breaks = seq_along(cluster_levels), labels = cluster_levels
    ) +
    scale_y_continuous(n.breaks = 6) +
    labs(
      x = xlab_use,
      y = ylab %||% "Average log2FoldChange"
    ) +
    do.call(theme_use, theme_args) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_line(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(),
      legend.position = "right",
      legend.background = element_blank(),
      aspect.ratio = aspect.ratio
    )
  p4 <- p3 +
    ggrepel::geom_text_repel(
      data = top_marker,
      aes(x = x_plot, y = y_plot, label = gene),
      inherit.aes = FALSE,
      min.segment.length = 0,
      max.overlaps = 100,
      segment.colour = "grey40",
      color = label.fg,
      bg.color = label.bg,
      bg.r = label.bg.r,
      size = label.size,
      force = 20
    )
  if (isTRUE(show_group_track)) {
    p4 <- p4 +
      geom_tile(
        aes(x = x, y = y, fill = group1),
        color = "black",
        height = 0.5,
        show.legend = FALSE,
        inherit.aes = FALSE,
        data = tile_data
      ) +
      scale_fill_manual(values = tile_colors) +
      geom_text(
        aes(x = x, y = y, label = group1),
        inherit.aes = FALSE,
        data = tile_data,
        size = 3,
        color = "black"
      )
  }
  p4 +
    theme(
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

#' @title DEtest Ring Plot
#'
#' @description
#' Draw a circular (ring) plot of differential expression results by cell type.
#'
#' @md
#' @inheritParams DEtestPlot
#'
#' @seealso [DEtestPlot], [RunDEtest], [VolcanoPlot], [DEtestManhattanPlot]
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
#' DEtestRingPlot(
#'   pancreas_sub,
#'   group.by = "CellType"
#' )
DEtestRingPlot <- function(
  srt,
  group.by = NULL,
  test.use = "wilcox",
  res = NULL,
  group_use = NULL,
  DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
  group_palette = "Chinese",
  group_palcolor = NULL,
  pt.size = 1,
  pt.alpha = 1,
  cols.highlight = "black",
  sizes.highlight = 1,
  alpha.highlight = 1,
  stroke.highlight = 0.5,
  nlabel = 5,
  features_label = NULL,
  label.fg = "black",
  label.bg = "white",
  label.bg.r = 0.1,
  label.size = 4,
  palette = "RdBu",
  palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  tile_height = 0.3,
  tile_gap = 0.1,
  jitter_width = 0.5,
  ring_segments = TRUE,
  seed = 11
) {
  check_r("geomtextpath", verbose = FALSE)
  if (is.null(group.by)) {
    group.by <- "custom"
  }
  data_res <- get_de_data(srt, group.by, test.use, DE_threshold, res, group_use)
  de_df <- data_res$de_df
  clip_res <- clip_log2fc_symmetric(de_df)
  de_df_marker <- clip_res$df
  fc_lim <- clip_res$fc_lim
  cluster_levels <- levels(de_df_marker[["group1"]])
  top_marker <- get_top_markers_for_label(
    de_df_marker, cluster_levels, nlabel, features_label
  )
  n_grp <- length(cluster_levels)
  show_group_track <- n_grp > 1
  ring_r0 <- 2
  max_abs_fc <- max(abs(de_df_marker[, "avg_log2FC"]), na.rm = TRUE)
  ring_k <- if (max_abs_fc > 0) 1.2 / max_abs_fc else 0.2
  set.seed(seed)
  n_m <- nrow(de_df_marker)
  de_df_marker[, "x_num"] <- as.numeric(de_df_marker[, "group1"])
  de_df_marker[, "x_angle"] <- de_df_marker[, "x_num"] + (stats::runif(n_m) - 0.5) * jitter_width
  de_df_marker[, "y_radius"] <- ring_r0 + ring_k * de_df_marker[, "avg_log2FC"]
  band_lo <- ring_r0 - tile_height / 2
  band_hi <- ring_r0 + tile_height / 2
  in_band <- de_df_marker[, "y_radius"] >= band_lo & de_df_marker[, "y_radius"] <= band_hi
  de_df_marker[in_band & de_df_marker[, "y_radius"] < ring_r0, "y_radius"] <- band_lo - tile_gap
  de_df_marker[in_band & de_df_marker[, "y_radius"] >= ring_r0, "y_radius"] <- band_hi + tile_gap
  top_marker <- merge(
    top_marker,
    de_df_marker[, c("gene", "group1", "x_angle", "y_radius")],
    by = c("gene", "group1"),
    all.x = TRUE
  )
  tile_colors <- palette_colors(
    x = cluster_levels,
    palette = group_palette,
    palcolor = group_palcolor
  )
  tile_data <- data.frame(
    x = seq_len(n_grp),
    y = ring_r0,
    group1 = cluster_levels,
    label = cluster_levels
  )
  npt <- 40
  path_margin <- 0.04
  path_df <- do.call(rbind, lapply(seq_len(n_grp), function(i) {
    x_start <- i - 0.5 + path_margin
    x_end <- i + 0.5 - path_margin
    data.frame(
      x = x_start + (0:(npt - 1)) / (npt - 1) * (x_end - x_start),
      y = ring_r0,
      label = cluster_levels[i],
      group = i
    )
  }))
  p_ring <- ggplot(de_df_marker, aes(x = x_angle, y = y_radius))
  if (isTRUE(ring_segments) && isTRUE(show_group_track)) {
    p_ring <- p_ring +
      geom_vline(
        xintercept = seq(0.5, n_grp - 0.5, by = 1),
        color = "grey85",
        linewidth = 0.5
      )
  }
  p_ring <- p_ring +
    geom_hline(
      yintercept = ring_r0,
      linetype = 2,
      color = "grey40",
      linewidth = 0.5
    ) +
    geom_point(
      aes(x = x_angle, y = y_radius),
      color = cols.highlight,
      size = sizes.highlight + stroke.highlight,
      alpha = alpha.highlight,
      inherit.aes = FALSE,
      data = de_df_marker
    ) +
    geom_point(
      aes(color = avg_log2FC),
      size = pt.size,
      alpha = pt.alpha
    ) +
    scale_color_gradientn(
      name = "log2FC",
      colors = palette_colors(palette = palette, palcolor = palcolor),
      values = scales::rescale(unique(c(fc_lim[1], 0, fc_lim[2]))),
      limits = fc_lim,
      guide = guide_colorbar(
        frame.colour = "black",
        ticks.colour = "black",
        title.hjust = 0,
        order = 1
      )
    ) +
    scale_x_continuous(
      limits = c(0.5, n_grp + 0.5),
      breaks = seq_len(n_grp), labels = NULL
    ) +
    scale_y_continuous(limits = c(0, NA), n.breaks = 5) +
    coord_polar(theta = "x", start = -pi / 2, direction = 1) +
    ggrepel::geom_text_repel(
      data = top_marker,
      aes(x = x_angle, y = y_radius, label = gene),
      inherit.aes = FALSE,
      min.segment.length = 0,
      max.overlaps = 100,
      segment.colour = "grey40",
      color = label.fg,
      bg.color = label.bg,
      bg.r = label.bg.r,
      size = label.size,
      force = 5
    ) +
    labs(x = NULL, y = NULL, title = NULL, subtitle = NULL, caption = NULL) +
    do.call(theme_use, theme_args) +
    theme(
      aspect.ratio = 1,
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      legend.position = "right",
      legend.background = element_blank()
    )
  if (isTRUE(show_group_track)) {
    p_ring <- p_ring +
      geom_tile(
        data = tile_data,
        aes(x = x, y = y, fill = group1),
        color = "black",
        height = tile_height,
        show.legend = FALSE,
        inherit.aes = FALSE
      ) +
      scale_fill_manual(values = tile_colors) +
      geomtextpath::geom_textpath(
        aes(x = x, y = y, label = label, group = group),
        data = path_df,
        inherit.aes = FALSE,
        size = 3,
        color = "black",
        upright = TRUE
      )
  }
  p_ring
}

#' @title Volcano Plot
#'
#' @description
#' Generate a volcano plot based on differential expression analysis results.
#'
#' @md
#' @inheritParams DEtestPlot
#'
#' @seealso [DEtestPlot], [RunDEtest], [DEtestManhattanPlot], [DEtestRingPlot]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunDEtest(
#'   pancreas_sub,
#'   group.by = "CellType"
#' )
#' VolcanoPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   ncol = 2
#' )
#'
#' VolcanoPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   group_use = c("Ductal", "Endocrine"),
#'   ncol = 2
#' )
#'
#' VolcanoPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   DE_threshold = "abs(diff_pct) > 0.3 & p_val_adj < 0.05",
#'   ncol = 2
#' )
#'
#' VolcanoPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   x_metric = "avg_log2FC",
#'   y_metric = "p_val",
#'   DE_threshold = "abs(avg_log2FC) > log2(1.5) & p_val < 0.05",
#'   ncol = 2
#' )
#'
#' VolcanoPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   threshold_method = "hyperbolic",
#'   hyperbola_c = 6,
#'   ncol = 2
#' )
#'
#' pancreas_sub <- RunEnrichment(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   db = "GO_BP",
#'   species = "Mus_musculus"
#' )
#' VolcanoPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   threshold_method = "hyperbolic",
#'   hyperbola_c = 6,
#'   annotate_enrichment = TRUE,
#'   enrich_from = "Enrichment",
#'   enrich_db = "GO_BP",
#'   ncol = 2
#' )
VolcanoPlot <- function(
    srt,
    group.by = NULL,
    test.use = "wilcox",
    res = NULL,
    group_use = NULL,
    DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
    x_metric = NULL,
    y_metric = NULL,
    palette = "RdBu",
    palcolor = NULL,
    pt.size = 1,
    pt.alpha = 1,
    cols.background = "grey80",
    cols.highlight = "black",
    sizes.highlight = 1,
    alpha.highlight = 1,
    stroke.highlight = 0.5,
    nlabel = 5,
    features_label = NULL,
    label.fg = "black",
    label.bg = "white",
    label.bg.r = 0.1,
    label.size = 4,
    aspect.ratio = NULL,
    xlab = NULL,
    ylab = NULL,
    theme_use = "theme_scop",
    theme_args = list(),
    combine = TRUE,
    nrow = NULL,
    ncol = NULL,
    byrow = TRUE,
    threshold_method = c("rectangular", "hyperbolic"),
    hyperbola_c = 6,
    annotate_enrichment = FALSE,
    enrich_from = c("Enrichment", "GSEA", "GSVA"),
    enrich_db = NULL,
    enrich_terms = NULL,
    enrich_top_terms = 3,
    enrich_padj_cutoff = 0.05,
    enrich_gsva_score_cutoff = NULL,
    gsva_method = NULL,
    enrich_nlabel = 15) {
  DE_threshold_missing <- missing(DE_threshold)
  threshold_method <- match.arg(threshold_method)
  if (
    DE_threshold_missing &&
      test.use %in% c("edgeR", "limma")
  ) {
    DE_threshold <- "p_val < 0.05"
  }
  x_metric <- x_metric %||% ifelse(
    test.use %in% c("edgeR", "limma") || threshold_method == "hyperbolic",
    "avg_log2FC",
    "diff_pct"
  )
  y_metric <- y_metric %||% ifelse(
    test.use %in% c("edgeR", "limma"),
    "p_val",
    "p_val_adj"
  )
  if (!x_metric %in% c("diff_pct", "avg_log2FC")) {
    log_message(
      "'x_metric' must be either 'diff_pct' or 'avg_log2FC'",
      message_type = "error"
    )
  }
  if (!y_metric %in% c("p_val", "p_val_adj")) {
    log_message(
      "'y_metric' must be either 'p_val' or 'p_val_adj'",
      message_type = "error"
    )
  }
  if (threshold_method == "hyperbolic") {
    x_metric <- "avg_log2FC"
  }
  traditional_volcano <- x_metric == "avg_log2FC"
  xlab <- xlab %||% x_metric
  ylab <- ylab %||% ifelse(
    y_metric == "p_val",
    "-log10(p-value)",
    "-log10(p-adjust)"
  )
  enrich_nlabel <- as.integer(enrich_nlabel)[1]
  if (!is.finite(enrich_nlabel) || is.na(enrich_nlabel) || enrich_nlabel < 0) {
    enrich_nlabel <- 0L
  }
  enrich_from <- unique(intersect(
    enrich_from,
    c("Enrichment", "GSEA", "GSVA")
  ))

  data_res <- get_de_data(srt, group.by, test.use, DE_threshold, res, group_use)
  de_df <- data_res$de_df
  de_df[, "avg_log2FC_raw"] <- de_df[, "avg_log2FC"]
  de_df[, "p_val_adj_raw"] <- de_df[, "p_val_adj"]
  de_df[, "neglog10_pvalue"] <- get_safe_neglog10(de_df[, y_metric])

  clip_res <- clip_log2fc_symmetric(de_df)
  de_df <- clip_res$df
  x_upper <- clip_res$fc_lim[2]
  x_lower <- clip_res$fc_lim[1]
  de_df[, "border"] <- (de_df[["avg_log2FC"]] >= x_upper) | (de_df[["avg_log2FC"]] <= x_lower)
  volcano_groups <- unique(as.character(de_df[["group1"]]))
  volcano_groups <- volcano_groups[!is.na(volcano_groups) & nzchar(volcano_groups)]
  single_group_volcano <- length(volcano_groups) <= 1
  signed_diff_pct_volcano <- FALSE

  if (threshold_method == "hyperbolic") {
    de_df[, "DE"] <- compute_hyperbolic_de_flags(
      avg_log2FC = de_df[, "avg_log2FC_raw"],
      p_val_adj = de_df[, "p_val_adj_raw"],
      hyperbola_c = hyperbola_c
    )
    de_df[, "x"] <- de_df[, "avg_log2FC"]
    de_df[, "y"] <- de_df[, "neglog10_pvalue"]
    de_df <- de_df[
      order(abs(de_df[, "avg_log2FC_raw"]), decreasing = FALSE, na.last = FALSE), ,
      drop = FALSE
    ]
  } else if (x_metric == "diff_pct") {
    de_df[, "y"] <- de_df[, "neglog10_pvalue"]
    de_df[, "x"] <- de_df[, "diff_pct"]
    signed_diff_pct_volcano <- !isTRUE(single_group_volcano)
    if (isTRUE(signed_diff_pct_volcano)) {
      de_df[de_df[, "avg_log2FC"] < 0, "y"] <- -de_df[
        de_df[, "avg_log2FC"] < 0,
        "y"
      ]
    }
    de_df <- de_df[
      order(abs(de_df[, "avg_log2FC"]), decreasing = FALSE, na.last = FALSE), ,
      drop = FALSE
    ]
  } else if (x_metric == "avg_log2FC") {
    de_df[, "y"] <- de_df[, "neglog10_pvalue"]
    de_df[, "x"] <- de_df[, "avg_log2FC"]
    de_df <- de_df[
      order(abs(de_df[, "avg_log2FC"]), decreasing = FALSE, na.last = FALSE), ,
      drop = FALSE
    ]
  }
  de_df[, "distance"] <- de_df[, "x"]^2 + de_df[, "y"]^2

  enrichment_map <- data.frame()
  if (isTRUE(annotate_enrichment)) {
    enrichment_map <- collect_volcano_enrichment_annotations(
      srt = srt,
      group.by = group.by,
      test.use = test.use,
      enrich_from = enrich_from,
      enrich_db = enrich_db,
      enrich_terms = enrich_terms,
      enrich_top_terms = enrich_top_terms,
      enrich_padj_cutoff = enrich_padj_cutoff,
      enrich_gsva_score_cutoff = enrich_gsva_score_cutoff,
      gsva_method = gsva_method
    )
  }
  enrich_key_levels <- character(0)
  enrich_colors <- NULL
  if (nrow(enrichment_map) > 0) {
    enrichment_map[["enrich_key"]] <- paste0(
      enrichment_map[["source"]],
      " | ",
      enrichment_map[["Database"]]
    )
    enrich_key_levels <- unique(enrichment_map[["enrich_key"]])
    enrich_colors <- get_enrichment_overlay_colors(enrich_key_levels)
  }

  plist <- list()
  for (group in levels(de_df[["group1"]])) {
    df <- de_df[de_df[["group1"]] == group, , drop = FALSE]
    if (nrow(df) == 0) {
      next
    }
    df <- add_volcano_plot_coords(df = df, jitter_width = 0.2, jitter_height = 0.2, seed = 11)
    df_enrich <- enrichment_map[
      enrichment_map[["group1"]] == group & enrichment_map[["gene"]] %in% df[["gene"]], ,
      drop = FALSE
    ]
    if (nrow(df_enrich) > 0) {
      match_idx <- match(df_enrich[["gene"]], df[["gene"]])
      match_ok <- !is.na(match_idx)
      df_enrich <- df_enrich[match_ok, , drop = FALSE]
      match_idx <- match_idx[match_ok]
      df_enrich[["x_plot"]] <- df[match_idx, "x_plot"]
      df_enrich[["y_plot"]] <- df[match_idx, "y_plot"]
      df_enrich[["enrich_key"]] <- factor(
        as.character(df_enrich[["enrich_key"]]),
        levels = enrich_key_levels
      )
    }

    x_nudge <- diff(range(df[, "x_plot"], na.rm = TRUE)) * 0.05
    if (!is.finite(x_nudge) || x_nudge <= 0) {
      x_nudge <- 0.1
    }
    df[, "label"] <- FALSE
    if (is.null(features_label)) {
      if (threshold_method == "hyperbolic") {
        df[
          utils::head(order(df[, "distance"], decreasing = TRUE), nlabel * 2),
          "label"
        ] <- TRUE
      } else if (traditional_volcano) {
        right_idx <- which(df[["x"]] >= 0)
        left_idx <- which(df[["x"]] < 0)
        if (length(right_idx) > 0) {
          df[right_idx[
            utils::head(order(df[right_idx, "distance"], decreasing = TRUE), nlabel)
          ], "label"] <- TRUE
        }
        if (length(left_idx) > 0) {
          df[left_idx[
            utils::head(order(df[left_idx, "distance"], decreasing = TRUE), nlabel)
          ], "label"] <- TRUE
        }
      } else if (isTRUE(signed_diff_pct_volcano)) {
        up_idx <- which(df[["y"]] >= 0)
        down_idx <- which(df[["y"]] < 0)
        if (length(up_idx) > 0) {
          df[up_idx[
            utils::head(order(df[up_idx, "distance"], decreasing = TRUE), nlabel)
          ], "label"] <- TRUE
        }
        if (length(down_idx) > 0) {
          df[down_idx[
            utils::head(order(df[down_idx, "distance"], decreasing = TRUE), nlabel)
          ], "label"] <- TRUE
        }
      } else {
        up_idx <- which(df[["avg_log2FC_raw"]] >= 0)
        down_idx <- which(df[["avg_log2FC_raw"]] < 0)
        if (length(up_idx) > 0) {
          df[up_idx[
            utils::head(order(df[up_idx, "distance"], decreasing = TRUE), nlabel)
          ], "label"] <- TRUE
        }
        if (length(down_idx) > 0) {
          df[down_idx[
            utils::head(order(df[down_idx, "distance"], decreasing = TRUE), nlabel)
          ], "label"] <- TRUE
        }
      }
    } else {
      df[df[["gene"]] %in% features_label, "label"] <- TRUE
    }
    if (nrow(df_enrich) > 0 && is.finite(enrich_nlabel) && enrich_nlabel > 0) {
      enrich_genes <- df_enrich[order(
        df_enrich[["source_priority"]],
        df_enrich[["metric_order"]],
        na.last = TRUE
      ), "gene"]
      enrich_genes <- unique(enrich_genes)
      enrich_genes <- utils::head(enrich_genes, enrich_nlabel)
      df[df[["gene"]] %in% enrich_genes, "label"] <- TRUE
    }

    color_by <- "avg_log2FC"
    color_metric <- suppressWarnings(as.numeric(df[, color_by]))
    if (!any(is.finite(color_metric))) {
      color_metric <- c(-1, 1)
    }
    color_values <- unique(c(
      min(c(color_metric, 0), na.rm = TRUE),
      0,
      max(color_metric, na.rm = TRUE)
    ))
    if (length(color_values) == 1) {
      color_values <- c(color_values - 1, color_values, color_values + 1)
    }

    label_df <- df[df[["label"]], , drop = FALSE]
    label_nudge <- if (nrow(label_df) > 0) {
      if (threshold_method == "hyperbolic" || traditional_volcano || !isTRUE(signed_diff_pct_volcano)) {
        ifelse(label_df[["x_plot"]] >= 0, x_nudge, -x_nudge)
      } else {
        ifelse(label_df[["y_plot"]] >= 0, -x_nudge, x_nudge)
      }
    } else {
      numeric(0)
    }

    p <- ggplot() +
      geom_point(
        data = df[!df[["DE"]] & !df[["border"]], , drop = FALSE],
        aes(x = x_plot, y = y_plot),
        color = cols.background,
        size = pt.size,
        alpha = pt.alpha
      ) +
      geom_point(
        data = df[!df[["DE"]] & df[["border"]], , drop = FALSE],
        aes(x = x_plot, y = y_plot),
        color = cols.background,
        size = pt.size,
        alpha = pt.alpha
      ) +
      geom_point(
        data = df[df[["DE"]] & !df[["border"]], , drop = FALSE],
        aes(x = x_plot, y = y_plot),
        color = cols.highlight,
        size = sizes.highlight + stroke.highlight,
        alpha = alpha.highlight
      ) +
      geom_point(
        data = df[df[["DE"]] & df[["border"]], , drop = FALSE],
        aes(x = x_plot, y = y_plot),
        color = cols.highlight,
        size = sizes.highlight + stroke.highlight,
        alpha = alpha.highlight
      ) +
      geom_point(
        data = df[df[["DE"]] & !df[["border"]], , drop = FALSE],
        aes(x = x_plot, y = y_plot, color = .data[[color_by]]),
        size = pt.size,
        alpha = pt.alpha
      ) +
      geom_point(
        data = df[df[["DE"]] & df[["border"]], , drop = FALSE],
        aes(x = x_plot, y = y_plot, color = .data[[color_by]]),
        size = pt.size,
        alpha = pt.alpha
      ) +
      geom_hline(yintercept = 0, color = "black", linetype = 1) +
      geom_vline(xintercept = 0, color = "grey", linetype = 2) +
      ggrepel::geom_text_repel(
        data = label_df,
        aes(x = x_plot, y = y_plot, label = gene),
        min.segment.length = 0,
        max.overlaps = 100,
        segment.colour = "grey40",
        color = label.fg,
        bg.color = label.bg,
        bg.r = label.bg.r,
        size = label.size,
        force = 20,
        nudge_x = label_nudge
      ) +
      labs(x = xlab, y = ylab) +
      scale_color_gradientn(
        name = "log2FC",
        colors = palette_colors(palette = palette, palcolor = palcolor),
        values = scales::rescale(color_values),
        guide = guide_colorbar(
          frame.colour = "black",
          ticks.colour = "black",
          title.hjust = 0,
          order = 1
        )
      )

    if (threshold_method == "hyperbolic") {
      curve_df <- build_hyperbola_curve_df(
        x_range = range(df[, "x_plot"], na.rm = TRUE),
        hyperbola_c = hyperbola_c,
        y_max = max(df[, "y_plot"], na.rm = TRUE) * 1.05
      )
      if (nrow(curve_df) > 0) {
        p <- p + geom_line(
          data = curve_df,
          aes(x = x, y = y),
          linetype = "22",
          color = "black",
          linewidth = 0.4
        )
      }
      p <- p + scale_y_continuous()
    } else if (isTRUE(signed_diff_pct_volcano)) {
      p <- p + scale_y_continuous(labels = abs)
    } else {
      p <- p + scale_y_continuous()
    }

    if (nrow(df_enrich) > 0 && !is.null(enrich_colors)) {
      p <- p +
        geom_point(
          data = df_enrich,
          aes(x = x_plot, y = y_plot, fill = enrich_key),
          shape = 21,
          size = pt.size + 0.8,
          stroke = 0.5,
          color = "black",
          alpha = 1
        ) +
        scale_fill_manual(
          name = "Enrichment",
          values = enrich_colors
        )
    }

    p <- p +
      do.call(theme_use, theme_args) +
      theme(aspect.ratio = aspect.ratio)

    if (length(levels(de_df[["group1"]])) > 1) {
      p <- p + facet_wrap(~group1)
    }
    plist[[group]] <- p
  }
  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- patchwork::wrap_plots(
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
