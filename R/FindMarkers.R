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
    identical(opts$layer, "data") &&
    is.logical(opts$positive) &&
    length(opts$positive) == 1L &&
    !is.na(opts$positive) &&
    is.numeric(opts$logfc) && length(opts$logfc) == 1L &&
    is.finite(opts$logfc) &&
    is.numeric(opts$min_pct) && length(opts$min_pct) == 1L &&
    is.finite(opts$min_pct) &&
    is.numeric(opts$min_diff) && length(opts$min_diff) == 1L &&
    !is.na(opts$min_diff) &&
    is.numeric(opts$base) && length(opts$base) == 1L &&
    is.finite(opts$base) && opts$base > 0 && opts$base != 1 &&
    is.numeric(opts$cell_cap) && length(opts$cell_cap) == 1L &&
    !is.na(opts$cell_cap) && opts$cell_cap > 0 &&
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
  cellnames <- colnames(assay.obj)
  if (is.null(cellnames)) {
    return(NULL)
  }
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
  unnamed_dots <- length(dots) > 0L &&
    (is.null(names(dots)) || any(!nzchar(names(dots))))
  cells.1 <- dots[["cells.1"]]
  cells.2 <- dots[["cells.2"]]
  run_seurat <- function() {
    object.use <- object
    ident.1.use <- ident.1
    ident.2.use <- ident.2
    dots.use <- dots
    if (!is.null(dots.use[["layer"]])) {
      dots.use[["slot"]] <- dots.use[["layer"]]
      dots.use[["layer"]] <- NULL
    }
    dots.use[["norm.method"]] <- NULL
    if (
      is.null(ident.1.use) &&
        !is.null(cells.1) &&
        (is.null(group.by) || identical(group.by, "ident")) &&
        is.null(subset.ident)
    ) {
      cellnames.use <- if (is.null(reduction)) {
        assay.use <- assay %||% SeuratObject::DefaultAssay(object.use)
        colnames(object.use[[assay.use]])
      } else {
        rownames(object.use[[reduction]])
      }
      cells.1.use <- as.character(cells.1)
      cells.2.use <- if (is.null(cells.2)) {
        setdiff(cellnames.use, cells.1.use)
      } else {
        as.character(cells.2)
      }
      if (
        !all(cells.1.use %in% cellnames.use) ||
          !all(cells.2.use %in% cellnames.use)
      ) {
        stop("Cells in one or both requested groups are not present in the selected data.", call. = FALSE)
      }
      if (length(intersect(cells.1.use, cells.2.use)) > 0L) {
        stop("The requested cell groups must not overlap.", call. = FALSE)
      }
      temporary_idents <- stats::setNames(
        rep(".scop_other", ncol(object.use)),
        colnames(object.use)
      )
      temporary_idents[cells.1.use] <- ".scop_group1"
      temporary_idents[cells.2.use] <- ".scop_group2"
      SeuratObject::Idents(object.use) <- temporary_idents
      ident.1.use <- ".scop_group1"
      ident.2.use <- ".scop_group2"
      dots.use[["cells.1"]] <- NULL
      dots.use[["cells.2"]] <- NULL
    }
    do.call(
      get_namespace_fun("Seurat", "FindMarkers.Seurat"),
      c(
        list(
          object = object.use,
          ident.1 = ident.1.use,
          ident.2 = ident.2.use,
          latent.vars = latent.vars,
          group.by = group.by,
          subset.ident = subset.ident,
          assay = assay,
          reduction = reduction
        ),
        dots.use
      )
    )
  }

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
  extra_args <- if (unnamed_dots) dots else dots[setdiff(names(dots), supported_args)]
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
    marker_assay_is_chromatin(object, assay) ||
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
  presto_fun <- presto_get_fun(
    install = FALSE,
    error_on_missing = FALSE
  )
  if (is.null(presto_fun)) {
    return(run_seurat())
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
    verbose = verbose,
    presto_fun = presto_fun
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
  get("FindMarkers.default", asNamespace("Seurat"))(
    object = object,
    ...
  )
}
