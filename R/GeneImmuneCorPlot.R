#' @title Gene-immune correlation butterfly plot
#'
#' @description
#' Draw a linkET-style correlation butterfly plot. The heatmap shows
#' immune-cell or immune-signature correlations, and curved links show
#' correlations between target genes and immune abundances.
#'
#' @md
#' @param object Optional `Seurat`, `SummarizedExperiment`, deconvolution bundle,
#' or matrix-like object used to resolve expression and immune abundance data.
#' @param gene.data Optional gene expression matrix with genes in rows and
#' samples in columns.
#' @param immune.data Optional immune abundance matrix with samples in rows and
#' immune cell types or signatures in columns.
#' @param features Target genes to correlate with immune abundance.
#' @param immune.cols Metadata columns to extract from a `Seurat` object.
#' @param assay Assay used for `Seurat` or `SummarizedExperiment` expression.
#' @param layer Assay layer used for `Seurat` expression.
#' @param cor_method Correlation method. Default is `"spearman"`.
#' @param p_cutoff P-value cutoff used for positive/negative link categories.
#' @param abs_cor_breaks Two numeric breakpoints for link width categories.
#' Default `c(0.2, 0.4)` creates `< 0.2`, `0.2 - 0.4`, and `>= 0.4`.
#' @param heatmap_colors Continuous colors for immune-cell correlations.
#' @param link_colors Colors for `Positive`, `Negative`, and `Not` links.
#' @param link_sizes Link sizes for strong, medium, and weak absolute
#' correlations.
#' @param title,subtitle Plot title and subtitle.
#' @param theme_use Theme function name.
#' @param theme_args Additional theme arguments.
#' @param verbose Whether to print progress messages.
#'
#' @return A `ggplot` object.
#'
#' @export
GeneImmuneCorPlot <- function(
  object = NULL,
  gene.data = NULL,
  immune.data = NULL,
  features,
  immune.cols = NULL,
  assay = NULL,
  layer = "data",
  cor_method = c("spearman", "pearson", "kendall"),
  p_cutoff = 0.05,
  abs_cor_breaks = c(0.2, 0.4),
  heatmap_colors = c(
    "#2166AC",
    "#67A9CF",
    "#F7F7F7",
    "#F4A582",
    "#B2182B"
  ),
  link_colors = c(
    Positive = "#E71D36",
    Negative = "#0073C2",
    Not = "grey82"
  ),
  link_sizes = c(">= 0.4" = 1.5, "0.2 - 0.4" = 0.8, "< 0.2" = 0.3),
  title = NULL,
  subtitle = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  verbose = TRUE
) {
  if (!requireNamespace("linkET", quietly = TRUE)) {
    log_message(
      paste(
        "{.pkg linkET} is required for {.fn GeneImmuneCorPlot}.",
        "Install it before drawing the butterfly plot."
      ),
      message_type = "error"
    )
  }
  cor_method <- match.arg(cor_method)
  if (missing(features) || length(features) == 0L) {
    log_message(
      "{.arg features} must contain at least one target gene.",
      message_type = "error"
    )
  }
  features <- unique(as.character(features))
  if (length(abs_cor_breaks) != 2L || any(!is.finite(abs_cor_breaks))) {
    log_message(
      "{.arg abs_cor_breaks} must contain two finite numeric values.",
      message_type = "error"
    )
  }
  abs_cor_breaks <- sort(abs_cor_breaks)

  immune_mat <- resolve_immune_abundance(
    object = object,
    immune.data = immune.data,
    immune.cols = immune.cols
  )
  gene_mat <- resolve_gene_expression_matrix(
    object = object,
    gene.data = gene.data,
    features = features,
    assay = assay,
    layer = layer
  )

  common_samples <- intersect(rownames(gene_mat), rownames(immune_mat))
  if (length(common_samples) < 3L) {
    log_message(
      "Need at least three matched samples/spots/cells between gene and immune matrices.",
      message_type = "error"
    )
  }
  gene_mat <- gene_mat[common_samples, , drop = FALSE]
  immune_mat <- immune_mat[common_samples, , drop = FALSE]
  keep_immune <- apply(immune_mat, 2, function(x) stats::var(x, na.rm = TRUE) > 0)
  keep_gene <- apply(gene_mat, 2, function(x) stats::var(x, na.rm = TRUE) > 0)
  if (!all(keep_immune)) {
    log_message(
      "Drop immune columns with zero variance: {.val {colnames(immune_mat)[!keep_immune]}}",
      verbose = verbose
    )
    immune_mat <- immune_mat[, keep_immune, drop = FALSE]
  }
  if (!all(keep_gene)) {
    log_message(
      "Drop target genes with zero variance: {.val {colnames(gene_mat)[!keep_gene]}}",
      verbose = verbose
    )
    gene_mat <- gene_mat[, keep_gene, drop = FALSE]
  }
  if (ncol(immune_mat) < 2L || ncol(gene_mat) == 0L) {
    log_message(
      "Need at least two immune columns and one variable target gene.",
      message_type = "error"
    )
  }

  edge_df <- gene_immune_cor_edges(
    gene_mat = gene_mat,
    immune_mat = immune_mat,
    method = cor_method,
    p_cutoff = p_cutoff,
    abs_cor_breaks = abs_cor_breaks
  )
  abs_levels <- levels(edge_df$abs_Cor)
  if (length(link_sizes) == length(abs_levels)) {
    link_sizes <- stats::setNames(as.numeric(link_sizes), abs_levels)
  }
  theme_obj <- immune_plot_theme(theme_use = theme_use, theme_args = theme_args)

  linkET::qcorrplot(
    linkET::correlate(immune_mat, method = cor_method),
    type = "lower",
    diag = FALSE
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    linkET::geom_couple(
      ggplot2::aes(colour = pvalue, size = abs_Cor),
      data = edge_df,
      curvature = linkET::nice_curvature()
    ) +
    ggplot2::scale_fill_gradientn(
      colours = heatmap_colors,
      limits = c(-1, 1),
      name = "Cell-cell cor"
    ) +
    ggplot2::scale_size_manual(values = link_sizes, drop = FALSE) +
    ggplot2::scale_colour_manual(values = link_colors, drop = FALSE) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        title = "pvalue",
        override.aes = list(size = 3),
        order = 1
      ),
      size = ggplot2::guide_legend(
        title = "abs(Cor)",
        override.aes = list(colour = "grey35"),
        order = 2
      ),
      fill = ggplot2::guide_colorbar(title = "Cell-cell cor", order = 3)
    ) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    theme_obj
}

resolve_gene_expression_matrix <- function(
  object = NULL,
  gene.data = NULL,
  features,
  assay = NULL,
  layer = "data"
) {
  if (!is.null(gene.data)) {
    mat <- as.matrix(gene.data)
  } else if (inherits(object, "Seurat")) {
    assay_use <- assay %||% SeuratObject::DefaultAssay(object)
    missing_features <- setdiff(features, rownames(object@assays[[assay_use]]))
    if (length(missing_features) > 0L) {
      log_message(
        "Some target genes are missing and will be ignored: {.val {missing_features}}",
        message_type = "warning"
      )
    }
    features <- intersect(features, rownames(object@assays[[assay_use]]))
    if (length(features) == 0L) {
      log_message(
        "No target genes are available in the selected assay.",
        message_type = "error"
      )
    }
    mat <- GetAssayData5(
      object,
      assay = assay_use,
      layer = layer
    )[features, , drop = FALSE]
  } else if (methods::is(object, "SummarizedExperiment")) {
    assay_use <- assay %||% SummarizedExperiment::assayNames(object)[1]
    mat <- SummarizedExperiment::assay(object, assay_use)
  } else if (inherits(object, c("matrix", "data.frame", "Matrix"))) {
    mat <- as.matrix(object)
  } else {
    log_message(
      "Provide {.arg gene.data}, a {.cls Seurat}, a {.cls SummarizedExperiment}, or an expression matrix.",
      message_type = "error"
    )
  }

  mat <- as.matrix(mat)
  missing_features <- setdiff(features, rownames(mat))
  if (length(missing_features) > 0L) {
    log_message(
      "Some target genes are missing and will be ignored: {.val {missing_features}}",
      message_type = "warning"
    )
  }
  features <- intersect(features, rownames(mat))
  if (length(features) == 0L) {
    log_message(
      "No target genes are available in {.arg gene.data}.",
      message_type = "error"
    )
  }
  mat <- mat[features, , drop = FALSE]
  dim_names <- dimnames(mat)
  mat <- suppressWarnings(matrix(
    as.numeric(mat),
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dim_names
  ))
  out <- t(mat)
  out[!is.finite(out)] <- NA_real_
  out
}

gene_immune_cor_edges <- function(
  gene_mat,
  immune_mat,
  method = "spearman",
  p_cutoff = 0.05,
  abs_cor_breaks = c(0.2, 0.4)
) {
  out <- list()
  idx <- 1L
  for (gene in colnames(gene_mat)) {
    for (cell_type in colnames(immune_mat)) {
      cr <- tryCatch(
        {
          cor_args <- list(
            x = as.numeric(gene_mat[, gene]),
            y = as.numeric(immune_mat[, cell_type]),
            method = method
          )
          if (method %in% c("spearman", "kendall")) {
            cor_args$exact <- FALSE
          }
          do.call(stats::cor.test, cor_args)
        },
        error = function(e) NULL
      )
      r <- if (is.null(cr)) NA_real_ else as.numeric(cr$estimate)
      p <- if (is.null(cr)) NA_real_ else as.numeric(cr$p.value)
      out[[idx]] <- data.frame(
        gene = gene,
        path = cell_type,
        r = r,
        p = p,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  df <- do.call(rbind, out)
  df$p[is.na(df$p)] <- 1
  df$r[is.na(df$r)] <- 0
  df$pvalue <- ifelse(
    df$p >= p_cutoff,
    "Not",
    ifelse(df$r > 0, "Positive", "Negative")
  )
  df$pvalue <- factor(df$pvalue, levels = c("Positive", "Negative", "Not"))

  mid_label <- paste0(abs_cor_breaks[1], " - ", abs_cor_breaks[2])
  high_label <- paste0(">= ", abs_cor_breaks[2])
  low_label <- paste0("< ", abs_cor_breaks[1])
  df$abs_Cor <- ifelse(
    abs(df$r) >= abs_cor_breaks[2],
    high_label,
    ifelse(abs(df$r) >= abs_cor_breaks[1], mid_label, low_label)
  )
  df$abs_Cor <- factor(df$abs_Cor, levels = c(high_label, mid_label, low_label))
  df
}
