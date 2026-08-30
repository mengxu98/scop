marker_get_data <- function(object, assay_obj, layer) {
  if (inherits(assay_obj, "Assay") && !inherits(assay_obj, "StdAssay")) {
    if (!layer %in% methods::slotNames(assay_obj)) {
      return(NULL)
    }
    data_use <- methods::slot(assay_obj, layer)
  } else {
    layer_names <- SeuratObject::Layers(assay_obj, search = layer)
    if (length(layer_names) != 1L) {
      return(NULL)
    }
    data_use <- methods::slot(assay_obj, "layers")[[layer_names]]
  }
  if (is.null(data_use)) {
    return(NULL)
  }
  if (!inherits(data_use, "dgCMatrix")) {
    data_use <- tryCatch(
      methods::as(data_use, "dgCMatrix"),
      error = function(e) NULL
    )
  }
  if (!is.null(data_use)) {
    dimnames(data_use) <- list(rownames(assay_obj), colnames(assay_obj))
  }
  data_use
}

marker_assay_is_chromatin <- function(object, assay = NULL) {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  inherits(object[[assay]], "ChromatinAssay")
}

marker_all_supported <- function(
  extra_args,
  features,
  group.by,
  node,
  latent.vars,
  mean.fxn,
  densify,
  test.use,
  slot,
  max.cells.per.ident,
  random.seed,
  base
) {
  identical(test.use, "wilcox") &&
    identical(slot, "data") &&
    (is.null(features) || is.character(features)) &&
    is.null(node) &&
    is.null(latent.vars) &&
    is.null(mean.fxn) &&
    !isTRUE(densify) &&
    is.numeric(max.cells.per.ident) &&
    length(max.cells.per.ident) == 1L &&
    !is.na(max.cells.per.ident) &&
    max.cells.per.ident > 0 &&
    is.numeric(random.seed) &&
    length(random.seed) == 1L &&
    is.finite(random.seed) &&
    is.numeric(base) &&
    length(base) == 1L &&
    is.finite(base) &&
    base > 0 &&
    base != 1 &&
    length(extra_args) == 0L &&
    (
      is.null(group.by) ||
        identical(group.by, "ident") ||
        (is.character(group.by) && length(group.by) == 1L)
    )
}

marker_context <- function(
  object,
  assay,
  slot,
  base,
  fc.name,
  features = NULL,
  groups = NULL
) {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  assay_obj <- object[[assay]]
  data_use <- marker_get_data(object, assay_obj, slot)
  if (is.null(data_use) || !inherits(data_use, "dgCMatrix")) {
    return(NULL)
  }
  dimnames(data_use) <- list(rownames(assay_obj), colnames(assay_obj))
  if (!is.null(features)) {
    if (!all(features %in% rownames(data_use))) {
      return(NULL)
    }
    data_use <- data_use[features, , drop = FALSE]
  }
  if (any(!is.finite(data_use@x)) || any(data_use@x < 0)) {
    return(NULL)
  }

  cell_names <- colnames(data_use)
  if (is.null(groups)) {
    groups <- SeuratObject::Idents(object)
    if (is.null(names(groups)) || !all(cell_names %in% names(groups))) {
      return(NULL)
    }
    labels <- as.character(sort(unique(groups[cell_names])))
  } else {
    if (is.null(names(groups))) {
      return(NULL)
    }
    cell_names <- cell_names[cell_names %in% names(groups)]
    if (length(cell_names) == 0L) {
      return(NULL)
    }
    data_use <- data_use[, cell_names, drop = FALSE]
    labels <- if (is.factor(groups)) {
      levels(groups)
    } else {
      unique(as.character(groups))
    }
    labels <- labels[labels %in% as.character(groups[cell_names])]
  }

  group <- factor(as.character(groups[cell_names]), levels = labels)
  keep <- !is.na(group)
  if (!all(keep)) {
    data_use <- data_use[, keep, drop = FALSE]
    cell_names <- cell_names[keep]
    group <- droplevels(group[keep])
    labels <- levels(group)
  }
  sizes <- tabulate(as.integer(group), nbins = length(labels))
  list(
    data = data_use,
    features = rownames(data_use),
    labels = labels,
    group = group,
    sizes = sizes,
    cells = length(cell_names),
    feature_count = nrow(data_use),
    fc_name = fc.name %||% paste0(
      "avg_log",
      if (isTRUE(base == exp(1))) "" else base,
      "FC"
    )
  )
}

marker_bind <- function(pieces) {
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  out <- if (length(pieces) == 0L) data.frame() else do.call(rbind, pieces)
  if (nrow(out) > 0L && "gene" %in% colnames(out)) {
    rownames(out) <- make.unique(as.character(out$gene))
  }
  if (nrow(out) == 0L) {
    log_message("No DE genes identified", message_type = "warning")
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
  min.diff.pct,
  logfc.threshold,
  only.pos,
  return.thresh,
  base,
  p.adjust.method = "bonferroni",
  p.adjust.n = ctx$feature_count
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
  pass <- pass & abs(pct_1 - pct_2) >= min.diff.pct
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
    method = p.adjust.method,
    n = p.adjust.n
  )
  out <- out[out$p_val < return.thresh, , drop = FALSE]
  if (nrow(out) == 0L) {
    return(NULL)
  }
  out$cluster <- factor(rep(label, nrow(out)), levels = ctx$labels)
  out$gene <- rownames(out)
  out
}

marker_all_from_context <- function(
  ctx,
  min.cells.group,
  min.pct,
  min.diff.pct,
  logfc.threshold,
  only.pos,
  return.thresh,
  base,
  max.cells.per.ident = Inf,
  random.seed = 1,
  p.adjust.method = "bonferroni",
  p.adjust.n = ctx$feature_count
) {
  if (is.finite(max.cells.per.ident)) {
    p_by_group <- vapply(seq_along(ctx$labels), function(i) {
      cells_1 <- which(ctx$group == ctx$labels[[i]])
      cells_2 <- which(ctx$group != ctx$labels[[i]])
      set.seed(random.seed)
      if (length(cells_1) > max.cells.per.ident) {
        cells_1 <- sample(cells_1, size = max.cells.per.ident)
      }
      if (length(cells_2) > max.cells.per.ident) {
        cells_2 <- sample(cells_2, size = max.cells.per.ident)
      }
      cells <- c(cells_1, cells_2)
      sampled <- ctx$data[, cells, drop = FALSE]
      sampled_group <- factor(
        c(rep("Group1", length(cells_1)), rep("Group2", length(cells_2))),
        levels = c("Group1", "Group2")
      )
      presto_result <- presto::wilcoxauc(
        X = sampled,
        y = sampled_group,
        verbose = FALSE
      )
      block <- presto_result[
        as.character(presto_result$group) == "Group1",
        ,
        drop = FALSE
      ]
      block$pval[match(ctx$features, block$feature)]
    }, numeric(ctx$feature_count))
  } else {
    presto_result <- presto::wilcoxauc(
      X = ctx$data,
      y = ctx$group,
      verbose = FALSE
    )
    p_by_group <- vapply(ctx$labels, function(label) {
      block <- presto_result[as.character(presto_result$group) == label, , drop = FALSE]
      block$pval[match(ctx$features, block$feature)]
    }, numeric(ctx$feature_count))
  }
  if (!is.matrix(p_by_group)) {
    p_by_group <- matrix(p_by_group, ncol = length(ctx$labels))
  }
  p_by_group[is.nan(p_by_group)] <- 1

  membership <- Matrix::sparseMatrix(
    i = seq_len(ctx$cells),
    j = as.integer(ctx$group),
    x = 1,
    dims = c(ctx$cells, length(ctx$labels))
  )
  linear_data <- ctx$data
  linear_data@x <- expm1(linear_data@x)
  detected_data <- ctx$data
  detected_data@x <- as.numeric(detected_data@x > 0)
  detected_data <- Matrix::drop0(detected_data)
  stats <- list(
    sum = as.matrix(linear_data %*% membership),
    detected = as.matrix(detected_data %*% membership),
    p = p_by_group
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
      min.diff.pct = min.diff.pct,
      logfc.threshold = logfc.threshold,
      only.pos = only.pos,
      return.thresh = return.thresh,
      base = base,
      p.adjust.method = p.adjust.method,
      p.adjust.n = p.adjust.n
    )
  }))
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
  dots <- list(...)
  run_seurat <- function() {
    do.call(
      get_namespace_fun("Seurat", "FindAllMarkers"),
      c(
        list(
          object = object,
          assay = assay,
          features = features,
          group.by = group.by,
          logfc.threshold = logfc.threshold,
          test.use = test.use,
          slot = slot,
          min.pct = min.pct,
          min.diff.pct = min.diff.pct,
          node = node,
          verbose = verbose,
          only.pos = only.pos,
          max.cells.per.ident = max.cells.per.ident,
          random.seed = random.seed,
          latent.vars = latent.vars,
          min.cells.feature = min.cells.feature,
          min.cells.group = min.cells.group,
          mean.fxn = mean.fxn,
          fc.name = fc.name,
          base = base,
          return.thresh = return.thresh,
          densify = densify
        ),
        dots
      )
    )
  }
  if (
    marker_assay_is_chromatin(object, assay) ||
      !requireNamespace("presto", quietly = TRUE) ||
      !marker_all_supported(
        extra_args = list(...),
        features = features,
        group.by = group.by,
        node = node,
        latent.vars = latent.vars,
        mean.fxn = mean.fxn,
        densify = densify,
        test.use = test.use,
        slot = slot,
        max.cells.per.ident = max.cells.per.ident,
        random.seed = random.seed,
        base = base
      )
  ) {
    return(run_seurat())
  }

  groups <- NULL
  if (!is.null(group.by) && !identical(group.by, "ident")) {
    metadata <- methods::slot(object, "meta.data")
    if (!group.by %in% colnames(metadata)) {
      return(run_seurat())
    }
    groups <- metadata[[group.by]]
    names(groups) <- rownames(metadata)
  }
  ctx <- marker_context(
    object,
    assay,
    slot,
    base,
    fc.name,
    features = features,
    groups = groups
  )
  if (is.null(ctx) || length(ctx$labels) < 2L) {
    return(run_seurat())
  }
  marker_all_from_context(
    ctx = ctx,
    min.cells.group = min.cells.group,
    min.pct = min.pct,
    min.diff.pct = min.diff.pct,
    logfc.threshold = logfc.threshold,
    only.pos = only.pos,
    return.thresh = return.thresh,
    base = base,
    max.cells.per.ident = max.cells.per.ident,
    random.seed = random.seed,
    p.adjust.n = nrow(object[[assay %||% SeuratObject::DefaultAssay(object)]])
  )
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
  log_message("FindAllMarkers supports Seurat objects.", message_type = "error")
}
