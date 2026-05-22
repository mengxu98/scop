#' @title Compute LISI scores on a Seurat object
#'
#' @description
#' Compute per-cell Local Inverse Simpson's Index (LISI) scores from a
#' dimensional reduction and store them in the `meta.data` and `tools` slots
#' of a `Seurat` object.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams CellDimPlot
#' @param srt A `Seurat` object.
#' @param reductions Character vector of dimensional reductions used to compute LISI.
#' If `NULL`, [DefaultReduction()] is used.
#' @param reduction Deprecated alias of `reductions`.
#' @param dims Dimensions to use from the reduction. Default is `NULL`,
#' which uses all available dimensions.
#' @param label_colnames Character vector of metadata columns used for LISI.
#' If `NULL`, `RunLISI()` will try to use `srt@misc[["integration_batch"]]`.
#' @param prefix Prefix used for the stored LISI metadata columns.
#' If `NULL`, the reduction names are used.
#' @param tool_name Name used to store detailed results in `srt@tools`.
#' Default is `"LISI"` when multiple reductions are provided, otherwise
#' `paste0(prefix, "_LISI")`.
#' @param perplexity Effective neighborhood size. Default is `30`.
#' @param nn_method Nearest-neighbor backend. One of `"auto"` or `"exact"`.
#' Default is `"auto"`, which lets `thisutils` choose the fastest exact
#' backend available.
#' Requires the accelerated `thisutils::compute_lisi()` interface that exposes
#' `nn_method = c("auto", "exact")`.
#' @param tol Tolerance used in the binary search for the target perplexity.
#' Default is `1e-5`.
#' @param max_iter Maximum number of binary-search iterations. Default is `50`.
#' @param overwrite Whether to overwrite existing metadata columns. Default is `TRUE`.
#'
#' @return A modified `Seurat` object.
#' @export
#'
#' @seealso
#' [thisutils::compute_lisi], [LISIPlot]
#'
#' @examples
#' data(panc8_sub)
#' panc8_sub <- integration_scop(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_method = "Harmony5"
#' )
#' names(panc8_sub@reductions)
#'
#' panc8_sub <- RunLISI(
#'   panc8_sub,
#'   reductions = c("pcaUMAP2D", "Harmony5UMAP2D")
#' )
#' LISIPlot(
#'   panc8_sub,
#'   combine = TRUE
#' )
RunLISI <- function(
  srt,
  reductions = NULL,
  reduction = NULL,
  dims = NULL,
  label_colnames = NULL,
  prefix = NULL,
  tool_name = NULL,
  perplexity = 30,
  nn_method = c("auto", "exact"),
  tol = 1e-5,
  max_iter = 50,
  overwrite = TRUE,
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat}",
      message_type = "error"
    )
  }

  reductions <- reductions %||% reduction %||% DefaultReduction(srt)
  reductions <- unique(as.character(reductions))
  missing_reductions <- setdiff(reductions, SeuratObject::Reductions(srt))
  if (length(missing_reductions) > 0) {
    log_message(
      "Reductions not found in {.cls Seurat}: {.val {missing_reductions}}",
      message_type = "error"
    )
  }

  if (is.null(label_colnames)) {
    label_colnames <- srt@misc[["integration_batch"]] %||% NULL
  }
  if (is.null(label_colnames) || length(label_colnames) == 0) {
    log_message(
      "{.arg label_colnames} must contain at least one metadata column, or {.val integration_batch} must be stored in {.arg srt@misc}. Objects returned by {.fn integration_scop} store this automatically.",
      message_type = "error"
    )
  }

  if (!all(label_colnames %in% colnames(srt@meta.data))) {
    missing_cols <- setdiff(label_colnames, colnames(srt@meta.data))
    log_message(
      "The following metadata columns are missing: {.val {missing_cols}}",
      message_type = "error"
    )
  }
  nn_method <- match.arg(nn_method)
  compute_lisi_args <- names(formals(thisutils::compute_lisi))
  if (!("nn_method" %in% compute_lisi_args)) {
    log_message(
      "{.fn RunLISI} requires the accelerated {.fn thisutils::compute_lisi} interface with {.arg nn_method = c('auto', 'exact')}. Please install the updated {.pkg thisutils}.",
      message_type = "error"
    )
  }

  if (is.null(prefix)) {
    prefix <- reductions
  }
  if (length(prefix) == 1 && length(reductions) > 1) {
    prefix <- rep(prefix, length(reductions))
  }
  if (length(prefix) != length(reductions)) {
    log_message(
      "{.arg prefix} must have length 1 or the same length as {.arg reductions}",
      message_type = "error"
    )
  }
  tool_name <- tool_name %||% if (length(reductions) > 1) {
    "LISI"
  } else {
    paste0(prefix[[1]], "_LISI")
  }

  lisi_df_all <- list()
  lisi_cols_all <- character(0)
  dims_all <- list()
  for (i in seq_along(reductions)) {
    reduction_i <- reductions[[i]]
    prefix_i <- prefix[[i]] %||% reduction_i
    if (!nzchar(prefix_i)) {
      prefix_i <- reduction_i
    }

    emb <- Seurat::Embeddings(srt, reduction = reduction_i)
    dims_i <- dims
    if (is.null(dims_i)) {
      dims_i <- seq_len(ncol(emb))
    }
    dims_i <- unique(as.integer(dims_i))
    if (anyNA(dims_i) || any(dims_i < 1) || any(dims_i > ncol(emb))) {
      log_message(
        "{.arg dims} must be within the available dimensions of {.val {reduction_i}}",
        message_type = "error"
      )
    }

    lisi_cols <- make.names(
      paste0(prefix_i, "_", label_colnames, "_LISI"),
      unique = TRUE
    )
    if (!isTRUE(overwrite) && any(lisi_cols %in% colnames(srt@meta.data))) {
      existing_cols <- intersect(lisi_cols, colnames(srt@meta.data))
      log_message(
        "Metadata columns already exist: {.val {existing_cols}}. Set {.arg overwrite = TRUE} to replace them.",
        message_type = "error"
      )
    }

    log_message(
      "Compute {.pkg LISI} scores from reduction {.val {reduction_i}}",
      verbose = verbose
    )
    lisi_df <- thisutils::compute_lisi(
      X = emb[, dims_i, drop = FALSE],
      meta_data = srt@meta.data,
      label_colnames = label_colnames,
      perplexity = perplexity,
      nn_method = nn_method,
      tol = tol,
      max_iter = max_iter
    )
    colnames(lisi_df) <- lisi_cols
    lisi_df <- lisi_df[colnames(srt), , drop = FALSE]

    srt@meta.data[, lisi_cols] <- lisi_df
    lisi_df_all[[reduction_i]] <- lisi_df
    lisi_cols_all <- c(lisi_cols_all, lisi_cols)
    dims_all[[reduction_i]] <- dims_i
  }

  lisi_scores <- do.call(cbind, lisi_df_all)
  srt@tools[[tool_name]] <- list(
    scores = lisi_scores,
    reductions = reductions,
    reduction = if (length(reductions) == 1) reductions[[1]] else reductions,
    dims = dims_all,
    label_colnames = label_colnames,
    colnames = lisi_cols_all,
    perplexity = perplexity,
    nn_method = nn_method,
    tol = tol,
    max_iter = max_iter
  )

  log_message(
    "Stored {.pkg LISI} scores in metadata: {.val {lisi_cols_all}}",
    message_type = "success",
    text_color = "green",
    verbose = verbose
  )
  srt
}

#' @title Plot LISI scores
#'
#' @description
#' Backward-compatible wrapper around [BenchmarkPlot()] for LISI scores.
#' Visualize LISI scores on a dimensional reduction and compare methods with a
#' summary boxplot.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams CellDimPlot
#' @param srt A `Seurat` object.
#' @param features Metadata columns containing LISI scores.
#' Default is `NULL`, which will use columns stored in `tool_name`, or all
#' metadata columns ending with `"_LISI"` when `tool_name` is `NULL`.
#' @param tool_name Tool entry created by [RunLISI()]. Default is `NULL`.
#' @param reduction Dimensional reduction used for feature plots.
#' If `NULL`, the reduction recorded in `tool_name` is used when available;
#' otherwise [DefaultReduction()] is used.
#' @param plot_boxplot Whether to add boxplots. Default is `TRUE`.
#' @param boxplot_jitter Whether to overlay jittered points on boxplots.
#' Default is `FALSE`.
#'
#' @return
#' If `combine = TRUE`, returns a combined `patchwork` plot.
#' If `combine = FALSE`, returns a named list of ggplot objects.
#'
#' @export
#'
#' @seealso
#' [RunLISI], [FeatureDimPlot]
LISIPlot <- function(
  srt,
  features = NULL,
  tool_name = NULL,
  reduction = NULL,
  plot_boxplot = TRUE,
  boxplot_jitter = FALSE,
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  pt.size = NULL,
  pt.alpha = 1,
  palette = "Chinese",
  palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  verbose = TRUE,
  ...
) {
  BenchmarkPlot(
    srt = srt,
    features = features,
    tool_name = tool_name,
    reduction = reduction,
    plot_type = "auto",
    plot_boxplot = plot_boxplot,
    boxplot_jitter = boxplot_jitter,
    combine = combine,
    nrow = nrow,
    ncol = ncol,
    byrow = byrow,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    palette = palette,
    palcolor = palcolor,
    theme_use = theme_use,
    theme_args = theme_args,
    verbose = verbose,
    ...
  )
}

lisi_feature_boxplot <- function(
  srt,
  features,
  palette = "Chinese",
  palcolor = NULL,
  boxplot_jitter = FALSE,
  theme_use = "theme_scop",
  theme_args = list(),
  verbose = TRUE
) {
  clean_method_label <- function(x) {
    x <- gsub("_+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    raw_index <- grepl(
      "pca(UMAP2D|TSNE2D|DM2D|PHATE2D|PACMAP2D|TRIMAP2D|LARGEVIS2D|FR2D)?$",
      x,
      ignore.case = FALSE
    )
    x <- sub("UMAP2D$", "", x, ignore.case = TRUE)
    x <- sub("TSNE2D$", "", x, ignore.case = TRUE)
    x <- sub("DM2D$", "", x, ignore.case = TRUE)
    x <- sub("PHATE2D$", "", x, ignore.case = TRUE)
    x <- sub("PACMAP2D$", "", x, ignore.case = TRUE)
    x <- sub("TRIMAP2D$", "", x, ignore.case = TRUE)
    x <- sub("LARGEVIS2D$", "", x, ignore.case = TRUE)
    x <- sub("FR2D$", "", x, ignore.case = TRUE)
    x <- sub("(?<![A-Za-z])pca$", "", x, perl = TRUE)
    x <- sub("(?<![A-Za-z])PCA$", "", x, perl = TRUE)
    x <- sub("lsi$", "", x, ignore.case = FALSE)
    x <- sub("LSI$", "", x, ignore.case = FALSE)
    x <- gsub("^_+|_+$", "", x)
    x[raw_index | grepl("^pca$", x, ignore.case = TRUE)] <- "Raw"
    x[nchar(x) == 0] <- NA_character_
    x
  }

  common_suffix <- function(x) {
    if (length(x) == 0) {
      return("")
    }
    rev_split <- lapply(x, function(val) strsplit(paste(rev(strsplit(val, "")[[1]]), collapse = ""), "")[[1]])
    min_len <- min(vapply(rev_split, length, integer(1)))
    chars <- character(0)
    for (i in seq_len(min_len)) {
      current <- vapply(rev_split, `[`, character(1), i)
      if (length(unique(current)) != 1) {
        break
      }
      chars <- c(chars, current[1])
    }
    if (length(chars) == 0) {
      return("")
    }
    paste(rev(chars), collapse = "")
  }

  plot_list <- lapply(features, function(feature) {
    data.frame(
      feature = feature,
      score = srt@meta.data[[feature]],
      stringsAsFactors = FALSE
    )
  })
  plot_df <- do.call(rbind, plot_list)
  plot_df <- plot_df[!is.na(plot_df$score), , drop = FALSE]

  if (nrow(plot_df) == 0) {
    log_message(
      "No valid observations available for LISI boxplot",
      message_type = "warning",
      verbose = verbose
    )
    return(ggplot2::ggplot() +
      ggplot2::theme_void())
  }

  feature_suffix <- common_suffix(features)
  feature_labels <- features
  if (nzchar(feature_suffix)) {
    feature_labels <- sub(
      paste0(gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", feature_suffix), "$"),
      "",
      features
    )
  }
  feature_labels <- gsub("_+$", "", feature_labels)
  feature_labels <- clean_method_label(feature_labels)
  feature_labels[is.na(feature_labels)] <- features[is.na(feature_labels)]
  feature_labels[!nzchar(feature_labels)] <- features[!nzchar(feature_labels)]
  label_map <- stats::setNames(feature_labels, features)
  plot_df$feature_label <- label_map[plot_df$feature]

  mean_df <- stats::aggregate(
    score ~ feature_label,
    data = plot_df,
    FUN = mean
  )
  mean_df <- mean_df[order(mean_df$score, decreasing = FALSE), , drop = FALSE]
  plot_df$feature_label <- factor(
    plot_df$feature_label,
    levels = mean_df$feature_label
  )

  fill_cols <- palette_colors(
    levels(plot_df$feature_label),
    palette = palette,
    palcolor = palcolor
  )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data[["feature_label"]],
      y = .data[["score"]],
      fill = .data[["feature_label"]]
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.5,
      outlier.shape = NA,
      linewidth = 0.4,
      alpha = 0.9
    )

  if (isTRUE(boxplot_jitter)) {
    p <- p +
      ggplot2::geom_jitter(
        width = 0.12,
        size = 0.25,
        alpha = 0.12,
        color = "grey50"
      )
  }

  if (
    length(levels(plot_df$feature_label)) >= 2 &&
      "Raw" %in% levels(plot_df$feature_label)
  ) {
    check_r("ggpubr", verbose = FALSE)
    comparisons <- lapply(
      setdiff(levels(plot_df$feature_label), "Raw"),
      function(method) c("Raw", method)
    )
    y_top <- max(plot_df$score, na.rm = TRUE) * 1.2
    p <- p +
      ggpubr::stat_compare_means(
        comparisons = comparisons,
        method = "wilcox.test",
        label = "p.signif",
        step.increase = 0.08,
        tip.length = 0.03,
        size = 3
      ) +
      ggplot2::coord_cartesian(ylim = c(0, y_top))
  }

  p +
    ggplot2::scale_fill_manual(values = fill_cols, drop = FALSE) +
    ggplot2::labs(
      x = "",
      y = "LISI"
    ) +
    do.call(theme_use, theme_args) +
    ggplot2::theme(
      aspect.ratio = 1,
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "none",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}
