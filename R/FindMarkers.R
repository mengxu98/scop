marker_pair_options <- function(dots) {
  list(
    test = dots[["test.use"]] %||% "wilcox",
    layer = dots[["slot"]] %||% "data",
    logfc = dots[["logfc.threshold"]] %||% 0.1,
    min_pct = dots[["min.pct"]] %||% 0.01,
    base = dots[["base"]] %||% 2,
    positive = dots[["only.pos"]] %||% FALSE,
    min_diff = dots[["min.diff.pct"]] %||% -Inf,
    cell_cap = dots[["max.cells.per.ident"]] %||% Inf,
    min_group = dots[["min.cells.group"]] %||% 3,
    min_feature = dots[["min.cells.feature"]] %||% 3,
    pseudocount = dots[["pseudocount.use"]] %||% 1
  )
}

marker_pair_supported <- function(
  opts,
  ident.1,
  cells.1,
  latent.vars,
  group.by,
  subset.ident,
  reduction,
  features,
  mean.fxn,
  fc.name,
  norm.method,
  extra_args
) {
  is.null(latent.vars) &&
    is.null(group.by) &&
    is.null(subset.ident) &&
    is.null(reduction) &&
    is.null(mean.fxn) &&
    is.null(fc.name) &&
    identical(opts$test, "wilcox") &&
    opts$layer %in% c("data", "counts") &&
    (is.null(ident.1) || length(ident.1) == 1L) &&
    (!is.null(ident.1) || !is.null(cells.1)) &&
    (is.null(features) || is.character(features)) &&
    identical(norm.method, "LogNormalize") &&
    is.numeric(opts$min_group) &&
    length(opts$min_group) == 1L &&
    is.finite(opts$min_group) &&
    is.numeric(opts$min_feature) &&
    length(opts$min_feature) == 1L &&
    is.finite(opts$min_feature) &&
    is.numeric(opts$pseudocount) &&
    length(opts$pseudocount) == 1L &&
    is.finite(opts$pseudocount) &&
    opts$pseudocount > 0 &&
    length(extra_args) == 0L
}

marker_pair_cells <- function(object, assay, layer, ident.1, ident.2, cells.1, cells.2) {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  assay.obj <- object[[assay]]
  if (length(SeuratObject::Layers(assay.obj, search = layer)) > 1L) {
    return(NULL)
  }
  data.use <- SeuratObject::LayerData(object, layer = layer, assay = assay)
  if (is.null(data.use) || !inherits(data.use, "dgCMatrix")) {
    return(NULL)
  }

  cellnames <- colnames(data.use)
  if (!is.null(cells.1)) {
    if (!all(cells.1 %in% cellnames)) {
      return(NULL)
    }
    if (is.null(cells.2)) {
      cells.2 <- setdiff(cellnames, cells.1)
    }
    if (!all(cells.2 %in% cellnames) || length(intersect(cells.1, cells.2)) > 0L) {
      return(NULL)
    }
    return(list(cells.1 = cells.1, cells.2 = cells.2))
  }

  idents <- SeuratObject::Idents(object)
  if (is.null(names(idents)) || !all(cellnames %in% names(idents))) {
    return(NULL)
  }
  idents <- as.character(idents[cellnames])
  label_a <- as.character(ident.1)
  if (!(label_a %in% idents)) {
    return(NULL)
  }

  in_a <- idents == label_a
  if (is.null(ident.2)) {
    cells.1 <- cellnames[in_a]
    cells.2 <- cellnames[!in_a]
  } else {
    label_b <- as.character(ident.2)
    if (!all(label_b %in% idents)) {
      return(NULL)
    }
    in_b <- idents %in% label_b
    if (any(in_a & in_b)) {
      return(NULL)
    }
    cells.1 <- cellnames[in_a]
    cells.2 <- cellnames[in_b]
  }

  list(cells.1 = cells.1, cells.2 = cells.2)
}

marker_pair_parallel <- function(
  object,
  assay,
  layer,
  cells.1,
  cells.2,
  logfc.threshold,
  base,
  min.pct,
  min.diff.pct,
  min.cells.group,
  min.cells.feature,
  only.pos,
  pseudocount.use
) {
  if (
    !identical(layer, "data") ||
      is.finite(min.diff.pct) ||
      !identical(as.numeric(min.cells.group), 3) ||
      !identical(as.numeric(min.cells.feature), 3)
  ) {
    return(NULL)
  }
  data.use <- SeuratObject::LayerData(object, layer = layer, assay = assay)
  if (!inherits(data.use, "dgCMatrix")) {
    return(NULL)
  }
  cell.use <- c(cells.1, cells.2)
  data.use <- data.use[, cell.use, drop = FALSE]
  group <- factor(
    c(rep("group1", length(cells.1)), rep("group2", length(cells.2))),
    levels = c("group1", "group2")
  )
  group_sizes <- tabulate(as.integer(group), nbins = 2L)
  if (any(group_sizes < min.cells.group)) {
    return(NULL)
  }
  marker_result <- parallel_all_in_one_dgc(
    x_sexp = data.use,
    groups = as.integer(group),
    group_sizes = as.integer(group_sizes)
  )
  n1 <- group_sizes[[1L]]
  n2 <- group_sizes[[2L]]
  sums1 <- marker_result$sum_by_group[, 1L]
  sums2 <- marker_result$sum_by_group[, 2L]
  counts1 <- marker_result$detected_by_group[, 1L]
  counts2 <- marker_result$detected_by_group[, 2L]
  fc <- log((sums1 + pseudocount.use) / n1, base = base) -
    log((sums2 + pseudocount.use) / n2, base = base)
  pct1 <- round(counts1 / n1, digits = 3)
  pct2 <- round(counts2 / n2, digits = 3)
  pass <- pmax(pct1, pct2) >= min.pct
  pass <- pass & if (isTRUE(only.pos)) fc >= logfc.threshold else abs(fc) >= logfc.threshold
  selected <- which(pass)
  if (!length(selected)) {
    return(data.frame(
      p_val = numeric(),
      avg_log2FC = numeric(),
      pct.1 = numeric(),
      pct.2 = numeric(),
      p_val_adj = numeric()
    ))
  }
  out <- data.frame(
    p_val = marker_result$pval_by_group[selected, 1L],
    avg_log2FC = fc[selected],
    pct.1 = pct1[selected],
    pct.2 = pct2[selected],
    row.names = rownames(data.use)[selected],
    check.names = FALSE
  )
  if (isTRUE(only.pos)) {
    out <- out[out$avg_log2FC > 0, , drop = FALSE]
  }
  out <- out[order(out$p_val, -abs(out$pct.1 - out$pct.2)), , drop = FALSE]
  out$p_val_adj <- stats::p.adjust(out$p_val, method = "bonferroni", n = nrow(data.use))
  out
}

#' @export
FindMarkers.Seurat <- function(
  object,
  ident.1 = NULL,
  ident.2 = NULL,
  latent.vars = NULL,
  group.by = NULL,
  subset.ident = NULL,
  assay = NULL,
  reduction = NULL,
  ...
) {
  dots <- list(...)
  run_seurat <- function() {
    do.call(
      get_namespace_fun("Seurat", "FindMarkers.Seurat"),
      c(
        list(
          object = object,
          ident.1 = ident.1,
          ident.2 = ident.2,
          latent.vars = latent.vars,
          group.by = group.by,
          subset.ident = subset.ident,
          assay = assay,
          reduction = reduction
        ),
        dots
      )
    )
  }

  cells.1 <- dots[["cells.1"]]
  cells.2 <- dots[["cells.2"]]
  features <- dots[["features"]]
  test.use <- dots[["test.use"]] %||% "wilcox"
  layer <- dots[["layer"]] %||% dots[["slot"]] %||% "data"
  logfc.threshold <- dots[["logfc.threshold"]] %||% 0.1
  base <- dots[["base"]] %||% 2
  min.pct <- dots[["min.pct"]] %||% 0.01
  min.diff.pct <- dots[["min.diff.pct"]] %||% -Inf
  max.cells.per.ident <- dots[["max.cells.per.ident"]] %||% Inf
  min.cells.feature <- dots[["min.cells.feature"]] %||% 3
  min.cells.group <- dots[["min.cells.group"]] %||% 3
  only.pos <- dots[["only.pos"]] %||% FALSE
  norm.method <- dots[["norm.method"]] %||% "LogNormalize"
  pseudocount.use <- dots[["pseudocount.use"]] %||% 1
  mean.fxn <- dots[["mean.fxn"]]
  fc.name <- dots[["fc.name"]]
  verbose <- dots[["verbose"]] %||% TRUE
  random.seed <- dots[["random.seed"]] %||% 1
  supported_args <- c(
    "cells.1", "cells.2", "features", "test.use", "slot", "layer",
    "logfc.threshold", "base", "min.pct", "min.diff.pct",
    "max.cells.per.ident", "min.cells.feature", "min.cells.group",
    "only.pos", "norm.method", "pseudocount.use", "mean.fxn", "fc.name",
    "verbose", "random.seed"
  )
  extra_args <- dots[setdiff(names(dots), supported_args)]
  opts <- list(
    test = test.use,
    layer = layer,
    logfc = logfc.threshold,
    min_pct = min.pct,
    base = base,
    positive = only.pos,
    min_diff = min.diff.pct,
    cell_cap = max.cells.per.ident,
    min_group = min.cells.group,
    min_feature = min.cells.feature,
    pseudocount = pseudocount.use
  )
  if (
    !marker_pair_supported(
      opts = opts,
      ident.1 = ident.1,
      cells.1 = cells.1,
      latent.vars = latent.vars,
      group.by = group.by,
      subset.ident = subset.ident,
      reduction = reduction,
      features = features,
      mean.fxn = mean.fxn,
      fc.name = fc.name,
      norm.method = norm.method,
      extra_args = extra_args
    )
  ) {
    return(run_seurat())
  }
  cells <- marker_pair_cells(object, assay, opts$layer, ident.1, ident.2, cells.1, cells.2)
  if (is.null(cells)) {
    return(run_seurat())
  }
  fast <- tryCatch(
    marker_pair_parallel(
      object = object,
      assay = assay %||% SeuratObject::DefaultAssay(object),
      layer = opts$layer,
      cells.1 = cells$cells.1,
      cells.2 = cells$cells.2,
      logfc.threshold = opts$logfc,
      base = opts$base,
      min.pct = opts$min_pct,
      min.diff.pct = opts$min_diff,
      min.cells.group = opts$min_group,
      min.cells.feature = opts$min_feature,
      only.pos = opts$positive,
      pseudocount.use = opts$pseudocount
    ),
    error = function(e) NULL
  )
  if (!is.null(fast)) {
    return(fast)
  }
  RunDEtestSparseWilcoxMarkers(
    srt = object,
    assay = assay %||% SeuratObject::DefaultAssay(object),
    layer = opts$layer,
    cells.1 = cells$cells.1,
    cells.2 = cells$cells.2,
    features = features,
    logfc.threshold = opts$logfc,
    base = opts$base,
    min.pct = opts$min_pct,
    min.diff.pct = opts$min_diff,
    max.cells.per.ident = opts$cell_cap,
    min.cells.group = opts$min_group,
    only.pos = opts$positive,
    pseudocount.use = opts$pseudocount,
    random.seed = random.seed,
    verbose = verbose
  )
}

#' Find markers between groups
#'
#' @param object Object containing expression data.
#' @param ... Passed to methods.
#'
#' @return A data frame of marker statistics.
#' @export
FindMarkers <- function(object, ...) {
  UseMethod("FindMarkers")
}

#' @export
FindMarkers.default <- function(object, ...) {
  Seurat::FindMarkers(object = object, ...)
}
