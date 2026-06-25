marker_all_supported <- function(
  extra_args,
  features,
  group.by,
  node,
  latent.vars,
  mean.fxn,
  fc.name,
  densify,
  test.use,
  slot,
  max.cells.per.ident,
  min.diff.pct,
  base
) {
  identical(test.use, "wilcox") &&
    identical(slot, "data") &&
    is.null(features) &&
    is.null(node) &&
    is.null(latent.vars) &&
    is.null(mean.fxn) &&
    is.null(fc.name) &&
    !isTRUE(densify) &&
    is.infinite(max.cells.per.ident) &&
    identical(min.diff.pct, -Inf) &&
    isTRUE(base == 2) &&
    length(extra_args) == 0L &&
    (is.null(group.by) || identical(group.by, "ident"))
}

marker_context <- function(object, assay, slot, base, fc.name) {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  assay_obj <- object[[assay]]
  if (length(SeuratObject::Layers(assay_obj, search = slot)) > 1L) {
    stop(slot, " layers are not joined. Please run JoinLayers")
  }
  data_use <- assay_obj@layers[[slot]]
  if (is.null(data_use) || !inherits(data_use, "dgCMatrix")) {
    return(NULL)
  }
  dimnames(data_use) <- list(rownames(assay_obj), colnames(assay_obj))

  cell_names <- colnames(data_use)
  ident <- SeuratObject::Idents(object)
  if (is.null(names(ident)) || !all(cell_names %in% names(ident))) {
    return(NULL)
  }

  labels <- as.character(sort(unique(ident[cell_names])))
  group <- factor(as.character(ident[cell_names]), levels = labels)
  sizes <- tabulate(as.integer(group), nbins = length(labels))
  list(
    data = data_use,
    features = rownames(data_use),
    labels = labels,
    group = group,
    sizes = sizes,
    cells = length(cell_names),
    feature_count = nrow(data_use),
    fc_name = fc.name %||% paste0("avg_log", base, "FC")
  )
}

marker_bind <- function(pieces) {
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  out <- if (length(pieces) == 0L) data.frame() else do.call(rbind, pieces)
  if (nrow(out) == 0L) {
    warning("No DE genes identified", call. = FALSE, immediate. = TRUE)
  }
  out
}

marker_one_cluster <- function(
  label,
  label_index,
  ctx,
  totals,
  stats,
  min.cells.group,
  min.pct,
  logfc.threshold,
  only.pos,
  return.thresh,
  base
) {
  n1 <- ctx$sizes[[label_index]]
  n2 <- ctx$cells - n1
  if (n1 < min.cells.group || n2 < min.cells.group) {
    return(NULL)
  }

  sums_1 <- stats$sum[, label_index]
  counts_1 <- stats$detected[, label_index]
  counts_2 <- totals$detected - counts_1
  fc <- log((sums_1 + 1) / n1, base = base) -
    log((totals$sum - sums_1 + 1) / n2, base = base)
  pct_1 <- round(counts_1 / n1, digits = 3)
  pct_2 <- round(counts_2 / n2, digits = 3)
  pass <- pmax(pct_1, pct_2) >= min.pct
  pass <- pass &
    if (only.pos) fc >= logfc.threshold else abs(fc) >= logfc.threshold
  selected <- which(pass)
  if (length(selected) == 0L) {
    return(NULL)
  }

  out <- data.frame(
    p_val = stats$p[selected, label_index],
    fc[selected],
    pct_1[selected],
    pct_2[selected],
    row.names = ctx$features[selected],
    check.names = FALSE
  )
  colnames(out) <- c("p_val", ctx$fc_name, "pct.1", "pct.2")
  if (only.pos) {
    out <- out[out[[ctx$fc_name]] > 0, , drop = FALSE]
    if (nrow(out) == 0L) {
      return(NULL)
    }
  }
  out <- out[order(out$p_val, -abs(out$pct.1 - out$pct.2)), , drop = FALSE]
  out$p_val_adj <- p.adjust(
    out$p_val,
    method = "bonferroni",
    n = ctx$feature_count
  )
  out <- out[out$p_val < return.thresh, , drop = FALSE]
  if (nrow(out) == 0L) {
    return(NULL)
  }
  out$cluster <- factor(rep(label, nrow(out)), levels = ctx$labels)
  out$gene <- rownames(out)
  out
}

#' @export
FindAllMarkers.Seurat <- function(
  object,
  assay = NULL,
  features = NULL,
  group.by = NULL,
  logfc.threshold = 0.1,
  test.use = "wilcox",
  slot = "data",
  min.pct = 0.01,
  min.diff.pct = -Inf,
  node = NULL,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  mean.fxn = NULL,
  fc.name = NULL,
  base = 2,
  return.thresh = 1e-2,
  densify = FALSE,
  ...
) {
  if (
    !marker_all_supported(
      extra_args = list(...),
      features = features,
      group.by = group.by,
      node = node,
      latent.vars = latent.vars,
      mean.fxn = mean.fxn,
      fc.name = fc.name,
      densify = densify,
      test.use = test.use,
      slot = slot,
      max.cells.per.ident = max.cells.per.ident,
      min.diff.pct = min.diff.pct,
      base = base
    )
  ) {
    stop("FindAllMarkers.Seurat received unsupported arguments for the scop implementation.", call. = FALSE)
  }

  ctx <- marker_context(object, assay, slot, base, fc.name)
  if (is.null(ctx)) {
    stop("FindAllMarkers.Seurat could not prepare marker context.", call. = FALSE)
  }
  marker_result <- parallel_all_in_one_dgc(
    x_sexp = ctx$data,
    groups = as.integer(ctx$group),
    group_sizes = as.integer(ctx$sizes)
  )
  stats <- list(
    sum = marker_result$sum_by_group,
    detected = marker_result$detected_by_group,
    p = marker_result$pval_by_group
  )
  totals <- list(sum = rowSums(stats$sum), detected = rowSums(stats$detected))
  marker_bind(lapply(seq_along(ctx$labels), function(i) {
    marker_one_cluster(
      label = ctx$labels[[i]],
      label_index = i,
      ctx = ctx,
      totals = totals,
      stats = stats,
      min.cells.group = min.cells.group,
      min.pct = min.pct,
      logfc.threshold = logfc.threshold,
      only.pos = only.pos,
      return.thresh = return.thresh,
      base = base
    )
  }))
}

#' Find markers for all groups
#'
#' @param object Object containing expression data.
#' @param ... Passed to methods.
#'
#' @return A data frame of marker statistics.
#' @export
FindAllMarkers <- function(object, ...) {
  UseMethod("FindAllMarkers")
}

#' @export
FindAllMarkers.default <- function(object, ...) {
  stop("FindAllMarkers supports Seurat objects.", call. = FALSE)
}
