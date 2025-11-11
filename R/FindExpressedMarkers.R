#' @title FindExpressedMarkers
#'
#' @inheritParams Seurat::FindMarkers
#' @param layer The layer used.
#' @param min.expression The min.expression used.
#' @param seed The seed used.
#' @export
#'
#' @examples
#' markers <- FindExpressedMarkers(
#'   pancreas_sub,
#'   cells.1 = SeuratObject::WhichCells(
#'     pancreas_sub,
#'     expression = Phase == "G2M"
#'   )
#' )
#' head(markers)
#' FeatureStatPlot(
#'   pancreas_sub,
#'   rownames(markers)[1],
#'   group.by = "Phase",
#'   add_point = TRUE
#' )
FindExpressedMarkers <- function(
    object,
    ident.1 = NULL,
    ident.2 = NULL,
    cells.1 = NULL,
    cells.2 = NULL,
    features = NULL,
    assay = NULL,
    layer = "data",
    min.expression = 0,
    test.use = "wilcox",
    logfc.threshold = 0.25,
    base = 2,
    pseudocount.use = 1,
    mean.fxn = NULL,
    fc.name = NULL,
    min.pct = 0.1,
    min.diff.pct = -Inf,
    max.cells.per.ident = Inf,
    latent.vars = NULL,
    only.pos = FALSE,
    min.cells.group = 3,
    min.cells.feature = 3,
    norm.method = "LogNormalize",
    seed = 11,
    verbose = TRUE,
    ...) {
  assay <- assay %||% DefaultAssay(object)

  if (!is.null(cells.1)) {
    if (is.null(cells.2)) {
      cells.2 <- setdiff(colnames(object), cells.1)
    }
  } else {
    cells.1 <- SeuratObject::WhichCells(
      object = object,
      idents = ident.1
    )
    cells.2 <- SeuratObject::WhichCells(
      object = object,
      idents = ident.2
    )
  }

  if (!is.null(x = latent.vars)) {
    latent.vars <- SeuratObject::FetchData(
      object = object,
      vars = latent.vars,
      cells = c(cells.1, cells.2)
    )
  }

  data.use <- Seurat::GetAssayData(
    object = object,
    assay = assay,
    layer = layer
  )
  object <- data.use

  pseudocount.use <- pseudocount.use %||% 1
  data_layer <- ifelse(
    test = test.use %in% c("negbinom", "poisson", "DESeq2"),
    yes = "counts",
    no = layer
  )

  data.use <- data.use[Matrix::rowSums(data.use) > 0, ]
  data.use <- as_matrix(data.use)
  data.use[data.use <= min.expression] <- NA
  counts <- switch(
    EXPR = data_layer,
    "scale.data" = GetAssayData5(
      object = object,
      layer = "counts"
    ),
    numeric()
  )

  ########## FoldChange.Assay ##########
  features <- features %||% rownames(x = data.use)
  layer <- data_layer

  # By default run as if LogNormalize is done
  log1pdata.mean.fxn <- function(x) {
    return(log(
      x = rowMeans(x = expm1(x), na.rm = TRUE) + pseudocount.use,
      base = base
    ))
  }
  scaledata.mean.fxn <- function(x) {
    return(rowMeans(x, na.rm = TRUE))
  }
  counts.mean.fxn <- function(x) {
    return(
      log(
        x = rowMeans(x = x, na.rm = TRUE) + pseudocount.use,
        base = base
      )
    )
  }
  if (!is.null(x = norm.method)) {
    # For anything apart from log normalization set to rowMeans
    if (norm.method != "LogNormalize") {
      new.mean.fxn <- counts.mean.fxn
    } else {
      new.mean.fxn <- switch(
        EXPR = layer,
        "data" = log1pdata.mean.fxn,
        "scale.data" = scaledata.mean.fxn,
        "counts" = counts.mean.fxn,
        counts.mean.fxn
      )
    }
  } else {
    # If no normalization method is passed use slots to decide mean function
    new.mean.fxn <- switch(
      EXPR = layer,
      "data" = log1pdata.mean.fxn,
      "scale.data" = scaledata.mean.fxn,
      "counts" = counts.mean.fxn,
      log1pdata.mean.fxn
    )
  }
  mean.fxn <- mean.fxn %||% new.mean.fxn
  # Omit the decimal value of e from the column name if base == exp(1)
  base.text <- ifelse(
    test = base == exp(1),
    yes = "",
    no = base
  )
  fc.name <- fc.name %||%
    ifelse(
      test = layer == "scale.data",
      yes = "avg_diff",
      no = paste0("avg_log", base.text, "FC")
    )
  fc.results <- FoldChange.default(
    object = data.use,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    mean.fxn = mean.fxn,
    fc.name = fc.name
  )

  ########## FindMarkers.default ##########

  object <- data.use
  layer <- data_layer

  ValidateCellGroups <- get_namespace_fun("Seurat", "ValidateCellGroups")
  ValidateCellGroups(
    object = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    min.cells.group = min.cells.group
  )

  # reset parameters so no feature filtering is performed
  if (test.use %in% "DESeq2") {
    features <- rownames(x = object)
    min.diff.pct <- -Inf
    logfc.threshold <- 0
  }
  data <- switch(
    EXPR = layer,
    "scale.data" = counts,
    object
  )
  # feature selection (based on percentages)
  alpha.min <- pmax(fc.results$pct.1, fc.results$pct.2)
  names(x = alpha.min) <- rownames(x = fc.results)
  features <- names(x = which(x = alpha.min >= min.pct))
  if (length(x = features) == 0) {
    log_message(
      "No features pass min.pct threshold; returning empty data.frame",
      message_type = "warning"
    )
    return(fc.results[features, ])
  }
  alpha.diff <- alpha.min - pmin(fc.results$pct.1, fc.results$pct.2)
  features <- names(
    x = which(x = alpha.min >= min.pct & alpha.diff >= min.diff.pct)
  )
  if (length(x = features) == 0) {
    log_message(
      "No features pass min.diff.pct threshold; returning empty data.frame",
      message_type = "warning"
    )
    return(fc.results[features, ])
  }
  # feature selection (based on logFC)
  if (layer != "scale.data") {
    total.diff <- fc.results[, 1] # first column is logFC
    names(total.diff) <- rownames(fc.results)
    features.diff <- if (only.pos) {
      names(x = which(x = total.diff >= logfc.threshold))
    } else {
      names(x = which(x = abs(x = total.diff) >= logfc.threshold))
    }
    features <- intersect(x = features, y = features.diff)
    if (length(x = features) == 0) {
      log_message(
        "No features pass logfc.threshold threshold; returning empty data.frame",
        message_type = "warning"
      )
      return(fc.results[features, ])
    }
  }
  # subsample cell groups if they are too large
  if (max.cells.per.ident < Inf) {
    set.seed(seed)
    if (length(x = cells.1) > max.cells.per.ident) {
      cells.1 <- sample(x = cells.1, size = max.cells.per.ident)
    }
    if (length(x = cells.2) > max.cells.per.ident) {
      cells.2 <- sample(x = cells.2, size = max.cells.per.ident)
    }
    if (!is.null(x = latent.vars)) {
      latent.vars <- latent.vars[c(cells.1, cells.2), , drop = FALSE]
    }
  }

  de.results <- PerformDE(
    object = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    test.use = test.use,
    verbose = verbose,
    min.cells.feature = min.cells.feature,
    latent.vars = latent.vars,
    ...
  )
  de.results <- cbind(
    de.results,
    fc.results[rownames(x = de.results), , drop = FALSE]
  )
  if (only.pos) {
    de.results <- de.results[de.results[, 2] > 0, , drop = FALSE]
  }
  if (test.use %in% "roc") {
    de.results <- de.results[order(-de.results$power, -de.results[, 1]), ]
  } else {
    de.results <- de.results[order(de.results$p_val, -de.results[, 1]), ]
    de.results$p_val_adj <- stats::p.adjust(
      p = de.results$p_val,
      method = "bonferroni",
      n = nrow(x = object)
    )
  }

  return(de.results)
}

FindConservedMarkers2 <- function(
    object,
    grouping.var,
    ident.1,
    ident.2 = NULL,
    cells.1 = NULL,
    cells.2 = NULL,
    features = NULL,
    test.use = "wilcox",
    logfc.threshold = 0.25,
    base = 2,
    pseudocount.use = 1,
    mean.fxn = NULL,
    min.pct = 0.1,
    min.diff.pct = -Inf,
    max.cells.per.ident = Inf,
    latent.vars = NULL,
    only.pos = FALSE,
    assay = NULL,
    layer = "data",
    min.cells.group = 3,
    min.cells.feature = 3,
    meta.method = c(
      "maximump",
      "minimump",
      "wilkinsonp",
      "meanp",
      "sump",
      "votep"
    ),
    norm.method = "LogNormalize",
    verbose = TRUE,
    ...) {
  meta.method <- match.arg(meta.method)
  object.var <- SeuratObject::FetchData(object = object, vars = grouping.var)
  levels.split <- names(x = sort(x = table(object.var[, 1])))
  num.groups <- length(levels.split)
  assay <- assay %||% DefaultAssay(object)

  cells <- list()
  for (i in 1:num.groups) {
    cells[[i]] <- rownames(
      x = object.var[object.var[, 1] == levels.split[i], , drop = FALSE]
    )
  }
  marker.test <- list()
  if (is.null(cells.1)) {
    cells.1 <- SeuratObject::WhichCells(object = object, idents = ident.1)
    if (!is.null(ident.2)) {
      cells.2 <- cells.2 %||% SeuratObject::WhichCells(object = object, idents = ident.2)
    } else {
      cells.2 <- setdiff(colnames(object), cells.1)
    }
    object <- SeuratObject::SetIdent(
      object = object,
      cells = colnames(x = object),
      value = paste(SeuratObject::Idents(object = object), object.var[, 1], sep = "_")
    )
    ident.2.save <- ident.2
    for (i in 1:num.groups) {
      level.use <- levels.split[i]
      ident.use.1 <- paste(ident.1, level.use, sep = "_")
      ident.use.1.exists <- ident.use.1 %in% SeuratObject::Idents(object = object)
      if (!all(ident.use.1.exists)) {
        bad.ids <- ident.1[!ident.use.1.exists]
        log_message(
          "Identity: {.val {bad.ids}} not present in group {.val {level.use}}. Skipping {.val {level.use}}",
          message_type = "warning"
        )
        next
      }
      ident.2 <- ident.2.save
      cells.1.use <- SeuratObject::WhichCells(object = object, idents = ident.use.1)
      if (length(cells.1.use) < min.cells.group) {
        log_message(
          "{.val {level.use}} has fewer than {.val {min.cells.group}} cells in Identity: {.val {ident.1}}. Skipping {.val {level.use}}",
          message_type = "warning"
        )
        next
      }
      if (is.null(x = ident.2)) {
        cells.2.use <- setdiff(
          x = cells[[i]],
          y = cells.1.use
        )
        ident.use.2 <- names(
          x = which(
            x = table(SeuratObject::Idents(object = object)[cells.2.use]) > 0
          )
        )
        ident.2 <- gsub(
          pattern = paste0("_", level.use),
          replacement = "",
          x = ident.use.2
        )
        if (length(x = ident.use.2) == 0) {
          log_message(
            "Only one identity class present: {.val {ident.1}}",
            message_type = "error"
          )
        }
      } else {
        ident.use.2 <- paste(ident.2, level.use, sep = "_")
        cells.2.use <- SeuratObject::WhichCells(
          object = object,
          idents = ident.use.2
        )
      }
      if (length(cells.2.use) < min.cells.group) {
        log_message(
          "{.val {level.use}} has fewer than {.val {min.cells.group}} cells. Skipping {.val {level.use}}",
          message_type = "warning"
        )
        next
      }

      log_message(
        "Testing group {.val {level.use}}: ({.val {ident.1}}) vs ({.val {ident.2}})",
        verbose = verbose
      )

      ident.use.2.exists <- ident.use.2 %in% SeuratObject::Idents(object = object)
      if (!all(ident.use.2.exists)) {
        bad.ids <- ident.2[!ident.use.2.exists]
        log_message(
          "Identity: {.val {bad.ids}} not present in group {.val {level.use}}. Skipping {.val {level.use}}",
          message_type = "warning"
        )
        next
      }
      marker.test[[level.use]] <- Seurat::FindMarkers(
        object = SeuratObject::Assays(object, assay),
        layer = layer,
        cells.1 = cells.1.use,
        cells.2 = cells.2.use,
        features = features,
        test.use = test.use,
        logfc.threshold = logfc.threshold,
        min.pct = min.pct,
        min.diff.pct = min.diff.pct,
        max.cells.per.ident = max.cells.per.ident,
        min.cells.group = min.cells.group,
        min.cells.feature = min.cells.feature,
        norm.method = norm.method,
        base = base,
        pseudocount.use = pseudocount.use,
        mean.fxn = mean.fxn,
        latent.vars = latent.vars,
        only.pos = only.pos,
        verbose = verbose,
        ...
      )
    }
  } else {
    for (i in 1:num.groups) {
      level.use <- levels.split[i]
      cells.1.use <- intersect(cells[[i]], cells.1)
      if (length(cells.1.use) < min.cells.group) {
        log_message(
          "{.val {level.use}} has fewer than {.val {min.cells.group}} cells. Skipping {.val {level.use}}",
          message_type = "warning"
        )
        next
      }
      if (is.null(x = cells.2)) {
        cells.2.use <- setdiff(x = cells[[i]], y = cells.1.use)
      } else {
        cells.2.use <- intersect(cells[[i]], cells.2)
      }
      if (length(cells.2.use) < min.cells.group) {
        log_message(
          "{.val {level.use}} has fewer than {.val {min.cells.group}} cells. Skipping {.val {level.use}}",
          message_type = "warning"
        )
        next
      }
      log_message(
        "Testing group {.val {level.use}}: ({.val {cells.1}}) vs ({.val {cells.2}})",
        verbose = verbose
      )
      marker.test[[level.use]] <- Seurat::FindMarkers(
        object = SeuratObject::Assays(object, assay),
        layer = layer,
        cells.1 = cells.1.use,
        cells.2 = cells.2.use,
        features = features,
        test.use = test.use,
        logfc.threshold = logfc.threshold,
        min.pct = min.pct,
        min.diff.pct = min.diff.pct,
        max.cells.per.ident = max.cells.per.ident,
        min.cells.group = min.cells.group,
        min.cells.feature = min.cells.feature,
        norm.method = norm.method,
        base = base,
        pseudocount.use = pseudocount.use,
        mean.fxn = mean.fxn,
        latent.vars = latent.vars,
        only.pos = only.pos,
        verbose = verbose,
        ...
      )
    }
  }
  marker.test <- marker.test[!sapply(marker.test, is.null)]
  if (length(marker.test) == 0) {
    log_message(
      "No group was tested",
      message_type = "warning"
    )
    return(NULL)
  }
  genes.conserved <- Reduce(
    f = intersect,
    x = lapply(
      X = marker.test,
      FUN = function(x) {
        return(rownames(x = x))
      }
    )
  )
  markers.conserved <- list()
  for (i in 1:length(x = marker.test)) {
    markers.conserved[[i]] <- marker.test[[i]][genes.conserved, , drop = FALSE]
    colnames(
      x = markers.conserved[[i]]
    ) <- paste(
      names(x = marker.test)[i],
      colnames(x = markers.conserved[[i]]),
      sep = "_"
    )
  }
  markers.combined <- Reduce(cbind, markers.conserved)
  fc <- Seurat::FoldChange(
    SeuratObject::Assays(object, assay),
    layer = layer,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = genes.conserved,
    norm.method = norm.method,
    base = base,
    pseudocount.use = pseudocount.use,
    mean.fxn = mean.fxn
  )
  markers.combined <- cbind(
    markers.combined,
    fc[genes.conserved, , drop = FALSE]
  )
  logFC.codes <- colnames(
    x = markers.combined
  )[grepl(pattern = "*avg_log.*FC$", x = colnames(x = markers.combined))]
  if (isTRUE(only.pos)) {
    markers.combined <- markers.combined[
      apply(markers.combined[, logFC.codes] > 0, 1, all), ,
      drop = FALSE
    ]
  } else {
    markers.combined <- markers.combined[
      apply(markers.combined[, logFC.codes] < 0, 1, all) |
        apply(markers.combined[, logFC.codes] > 0, 1, all), ,
      drop = FALSE
    ]
  }
  pval.codes <- colnames(
    x = markers.combined
  )[grepl(pattern = "*_p_val$", x = colnames(x = markers.combined))]
  if (length(x = pval.codes) > 1) {
    markers.combined[["max_pval"]] <- apply(
      X = markers.combined[, pval.codes, drop = FALSE],
      MARGIN = 1,
      FUN = max
    )
    combined.pval <- data.frame(
      cp = apply(
        X = markers.combined[, pval.codes, drop = FALSE],
        MARGIN = 1,
        FUN = function(x) {
          return(metap(x, method = meta.method)$p)
        }
      )
    )
    meta.method.name <- meta.method
    colnames(x = combined.pval) <- paste0(meta.method.name, "_p_val")
    markers.combined <- cbind(markers.combined, combined.pval)
    markers.combined[, "p_val"] <- markers.combined[, paste0(
      meta.method.name,
      "_p_val"
    )]
    markers.combined <- markers.combined[
      order(markers.combined[, paste0(meta.method.name, "_p_val")]), ,
      drop = FALSE
    ]
  } else {
    markers.combined[, "max_pval"] <- markers.combined[
      ,
      "p_val"
    ] <- markers.combined[, pval.codes]
    log_message(
      "Only a single group was tested",
      message_type = "warning"
    )
  }
  return(markers.combined)
}

metap <- function(
    p,
    method = c(
      "maximump",
      "minimump",
      "wilkinsonp",
      "meanp",
      "sump",
      "votep"
    ),
    ...) {
  method <- match.arg(method)
  res <- do.call(method, args = list(p = p, ...))
  return(res)
}
