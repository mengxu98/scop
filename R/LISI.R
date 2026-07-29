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
#' panc8_sub <- RunIntegration(
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
      "{.arg label_colnames} must contain at least one metadata column, or {.val integration_batch} must be stored in {.arg srt@misc}. Objects returned by {.fn RunIntegration} store this automatically.",
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
  tool_name <- tool_name %||%
    if (length(reductions) > 1) {
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
  feature_boxplot(
    srt = srt,
    features = features,
    palette = palette,
    palcolor = palcolor,
    boxplot_jitter = boxplot_jitter,
    theme_use = theme_use,
    theme_args = theme_args,
    verbose = verbose,
    y_label = "LISI",
    empty_message = "No valid observations available for LISI boxplot"
  )
}
