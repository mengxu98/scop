#' @export
FoldChange.default <- function(
  object,
  cells.1,
  cells.2,
  mean.fxn,
  fc.name,
  features = NULL,
  ...
) {
  features <- features %||% rownames(x = object)
  thresh.min <- 0
  pct.1 <- round(
    x = Matrix::rowSums(
      x = object[features, cells.1, drop = FALSE] > thresh.min,
      na.rm = TRUE
    ) /
      length(cells.1),
    digits = 3
  )
  pct.2 <- round(
    x = Matrix::rowSums(
      x = object[features, cells.2, drop = FALSE] > thresh.min,
      na.rm = TRUE
    ) /
      length(cells.2),
    digits = 3
  )
  data.1 <- mean.fxn(object[features, cells.1, drop = FALSE])
  data.2 <- mean.fxn(object[features, cells.2, drop = FALSE])
  fc <- (data.1 - data.2)
  fc.results <- as.data.frame(x = cbind(fc, pct.1, pct.2))
  colnames(fc.results) <- c(fc.name, "pct.1", "pct.2")
  return(fc.results)
}


PerformDE <- function(
  object,
  cells.1,
  cells.2,
  features,
  test.use,
  verbose,
  min.cells.feature,
  latent.vars,
  ...
) {
  if (
    !(test.use %in% c("negbinom", "poisson", "MAST", "LR")) &&
      !is.null(latent.vars)
  ) {
    log_message(
      "'latent.vars' is only used for the following tests: ",
      paste(c("negbinom", "poisson", "MAST", "LR"), collapse = ", "),
      message_type = "warning",
      verbose = verbose
    )
  }
  data.use <- object[
    features,
    c(cells.1, cells.2),
    drop = FALSE
  ]
  data.use <- as_matrix(data.use)

  de.results <- switch(
    EXPR = test.use,
    "wilcox" = WilcoxDETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose,
      ...
    ),
    "bimod" = get_namespace_fun(
      "Seurat", "DiffExpTest"
    )(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    "roc" = get_namespace_fun(
      "Seurat", "MarkerTest"
    )(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    "t" = get_namespace_fun(
      "Seurat", "DiffTTest"
    )(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    "negbinom" = get_namespace_fun(
      "Seurat", "GLMDETest"
    )(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      min.cells = min.cells.feature,
      latent.vars = latent.vars,
      test.use = test.use,
      verbose = verbose
    ),
    "poisson" = get_namespace_fun(
      "Seurat", "GLMDETest"
    )(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      min.cells = min.cells.feature,
      latent.vars = latent.vars,
      test.use = test.use,
      verbose = verbose
    ),
    "MAST" = get_namespace_fun(
      "Seurat", "MASTDETest"
    )(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      latent.vars = latent.vars,
      verbose = verbose,
      ...
    ),
    "DESeq2" = get_namespace_fun(
      "Seurat", "DESeq2DETest"
    )(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose,
      ...
    ),
    "LR" = get_namespace_fun(
      "Seurat", "LRDETest"
    )(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      latent.vars = latent.vars,
      verbose = verbose
    ),
    log_message(
      "Unknown test: {.pkg {test.use}}",
      message_type = "error"
    )
  )
  de.results
}

WilcoxDETest <- function(
  data.use,
  cells.1,
  cells.2,
  verbose = TRUE,
  ...
) {
  data.use <- data.use[, c(cells.1, cells.2), drop = FALSE]
  check_r("limma", verbose = FALSE)
  p_val <- parallelize_fun(
    seq_len(nrow(data.use)),
    fun = function(x) {
      keep <- colnames(data.use)[!is.na(data.use[x, ])]
      j <- seq_len(length.out = length(intersect(cells.1, keep)))
      statistics <- data.use[x, keep]
      min(
        2 * min(
          limma::rankSumTestWithCorrelation(
            index = j,
            statistics = statistics
          )
        ),
        1
      )
    }
  ) |> purrr::list_c()
  data.frame(
    p_val,
    row.names = rownames(data.use)
  )
}

aggregate_counts_by_group <- function(counts, groups) {
  groups <- as.character(groups)
  groups <- groups[!is.na(groups)]
  group_levels <- unique(groups)
  agg <- lapply(group_levels, function(group) {
    Matrix::rowSums(counts[, groups == group, drop = FALSE])
  })
  agg <- do.call(cbind, agg)
  if (is.null(dim(agg))) {
    agg <- matrix(agg, ncol = 1)
  }
  rownames(agg) <- rownames(counts)
  colnames(agg) <- group_levels
  agg
}

RunDEtest_limma_voom <- function(
  count_matrix,
  condition,
  condition1 = NULL,
  condition2 = NULL,
  only.pos = TRUE,
  logfc.threshold = 0,
  p.adjust.method = "bonferroni"
) {
  check_r("limma", verbose = FALSE)
  check_r("edgeR", verbose = FALSE)

  condition_all <- as.character(condition)
  if (is.null(condition1) || is.null(condition2)) {
    condition_levels <- unique(condition_all)
    if (length(condition_levels) < 2) {
      return(NULL)
    }
    condition1 <- condition1 %||% condition_levels[[1]]
    condition2 <- condition2 %||% condition_levels[[2]]
  }

  keep <- condition_all %in% c(condition1, condition2)
  count_matrix <- count_matrix[, keep, drop = FALSE]
  condition <- factor(condition_all[keep], levels = c(condition1, condition2))
  if (ncol(count_matrix) < 2 || any(table(condition) < 2)) {
    return(NULL)
  }

  dge <- edgeR::DGEList(counts = count_matrix)
  keep_features <- edgeR::filterByExpr(dge, group = condition)
  if (!any(keep_features)) {
    return(data.frame())
  }
  dge <- dge[keep_features, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  design <- stats::model.matrix(~condition)
  v <- limma::voom(dge, design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = ncol(design), number = Inf, sort.by = "none")
  if (is.null(tt) || nrow(tt) == 0) {
    return(data.frame())
  }

  detect <- dge$counts[rownames(tt), , drop = FALSE] > 0
  pct.1 <- round(rowMeans(detect[, condition == condition1, drop = FALSE]), 3)
  pct.2 <- round(rowMeans(detect[, condition == condition2, drop = FALSE]), 3)
  out <- data.frame(
    p_val = tt$P.Value,
    avg_log2FC = tt$logFC,
    ave_expr = tt$AveExpr,
    pct.1 = pct.1,
    pct.2 = pct.2,
    row.names = rownames(tt)
  )
  out$p_val_adj <- stats::p.adjust(out$p_val, method = p.adjust.method)
  if (isTRUE(only.pos)) {
    out <- out[out$avg_log2FC >= logfc.threshold, , drop = FALSE]
  } else {
    out <- out[abs(out$avg_log2FC) >= logfc.threshold, , drop = FALSE]
  }
  out
}

RunDEtest_edgeR <- function(
  count_matrix,
  condition,
  condition1 = NULL,
  condition2 = NULL,
  only.pos = TRUE,
  logfc.threshold = 0,
  p.adjust.method = "bonferroni"
) {
  check_r("edgeR", verbose = FALSE)

  condition_all <- as.character(condition)
  if (is.null(condition1) || is.null(condition2)) {
    condition_levels <- unique(condition_all)
    if (length(condition_levels) < 2) {
      return(NULL)
    }
    condition1 <- condition1 %||% condition_levels[[1]]
    condition2 <- condition2 %||% condition_levels[[2]]
  }

  keep <- condition_all %in% c(condition1, condition2)
  count_matrix <- count_matrix[, keep, drop = FALSE]
  condition <- factor(condition_all[keep], levels = c(condition1, condition2))
  if (ncol(count_matrix) < 2 || any(table(condition) < 2)) {
    return(NULL)
  }

  dge <- edgeR::DGEList(counts = count_matrix)
  keep_features <- edgeR::filterByExpr(dge, group = condition)
  if (!any(keep_features)) {
    return(data.frame())
  }
  dge <- dge[keep_features, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  design <- stats::model.matrix(~condition)
  dge <- edgeR::estimateDisp(dge, design = design)
  fit <- edgeR::glmQLFit(dge, design = design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = ncol(design))
  tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  if (is.null(tt) || nrow(tt) == 0) {
    return(data.frame())
  }

  detect <- dge$counts[rownames(tt), , drop = FALSE] > 0
  pct.1 <- round(rowMeans(detect[, condition == condition1, drop = FALSE]), 3)
  pct.2 <- round(rowMeans(detect[, condition == condition2, drop = FALSE]), 3)
  out <- data.frame(
    p_val = tt$PValue,
    avg_log2FC = tt$logFC,
    ave_expr = tt$logCPM,
    pct.1 = pct.1,
    pct.2 = pct.2,
    row.names = rownames(tt)
  )
  out$p_val_adj <- stats::p.adjust(out$p_val, method = p.adjust.method)
  if (isTRUE(only.pos)) {
    out <- out[out$avg_log2FC >= logfc.threshold, , drop = FALSE]
  } else {
    out <- out[abs(out$avg_log2FC) >= logfc.threshold, , drop = FALSE]
  }
  out
}

RunDEtest_pseudobulk <- function(
  srt,
  group.by = NULL,
  group1 = NULL,
  group2 = NULL,
  features = NULL,
  feature_type = "gene",
  markers_type = "all",
  test.use = "limma_voom",
  only.pos = TRUE,
  fc.threshold = 1.5,
  base = 2,
  sample_col = NULL,
  condition_col = NULL,
  p.adjust.method = "bonferroni",
  layer = "counts",
  assay = NULL,
  verbose = TRUE,
  cores = 1,
  ...
) {
  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(srt)
  }
  if (!markers_type %in% c("all")) {
    log_message(
      "Pseudobulk differential testing currently supports only {.val all} markers.",
      message_type = "error"
    )
  }
  if (!test.use %in% c("limma_voom", "edgeR")) {
    log_message(
      "Pseudobulk differential testing currently supports only {.val limma_voom} and {.val edgeR}.",
      message_type = "error"
    )
  }
  if (is.null(sample_col) || !sample_col %in% colnames(srt@meta.data)) {
    log_message(
      "{.arg sample_col} must be a metadata column in {.cls Seurat} for pseudobulk analysis",
      message_type = "error"
    )
  }
  if (is.null(condition_col) || !condition_col %in% colnames(srt@meta.data)) {
    log_message(
      "{.arg condition_col} must be a metadata column in {.cls Seurat} for pseudobulk analysis",
      message_type = "error"
    )
  }
  if (!is.null(group.by) && !group.by %in% colnames(srt@meta.data)) {
    log_message(
      "{.arg group.by} must be a metadata column in {.cls Seurat} for pseudobulk analysis",
      message_type = "error"
    )
  }

  counts <- GetAssayData5(srt, layer = layer, assay = assay)
  features <- features %||% rownames(counts)
  features <- intersect(features, rownames(counts))
  if (length(features) == 0) {
    log_message(
      "No valid features available for pseudobulk differential testing",
      message_type = "error"
    )
  }
  counts <- counts[features, , drop = FALSE]

  condition_levels_all <- unique(as.character(stats::na.omit(srt[[condition_col, drop = TRUE]])))
  condition1 <- group1 %||% NULL
  condition2 <- group2 %||% NULL
  if (is.null(condition1) || is.null(condition2)) {
    if (length(condition_levels_all) < 2) {
      log_message(
        "Pseudobulk differential testing requires at least two condition levels",
        message_type = "error"
      )
    }
    condition1 <- condition1 %||% condition_levels_all[[1]]
    condition2 <- condition2 %||% condition_levels_all[[2]]
  }
  if (identical(condition1, condition2)) {
    log_message(
      "{.arg group1} and {.arg group2} must refer to two different condition labels in pseudobulk analysis",
      message_type = "error"
    )
  }

  target_groups <- if (is.null(group.by)) {
    "All"
  } else {
    unique(as.character(stats::na.omit(srt[[group.by, drop = TRUE]])))
  }

  log_message(
    "Start pseudobulk differential testing",
    verbose = verbose
  )

  run_group_fun <- if (isTRUE(cores == 1)) {
    function(X, FUN) lapply(X, FUN)
  } else {
    function(X, FUN) parallelize_fun(X, FUN, cores = cores, verbose = verbose)
  }

  all_markers <- run_group_fun(
    target_groups,
    function(current_group) {
      cells_use <- colnames(srt)
      if (!is.null(group.by)) {
        cells_use <- colnames(srt)[srt[[group.by, drop = TRUE]] %in% current_group]
      }
      meta_use <- srt@meta.data[cells_use, c(sample_col, condition_col), drop = FALSE]
      meta_use <- meta_use[stats::complete.cases(meta_use), , drop = FALSE]
      cells_use <- rownames(meta_use)
      if (length(cells_use) == 0) {
        return(NULL)
      }

      sample_condition <- unique(meta_use[, c(sample_col, condition_col), drop = FALSE])
      sample_tab <- table(meta_use[[sample_col]], meta_use[[condition_col]])
      if (any(rowSums(sample_tab > 0) > 1)) {
        log_message(
          "Each sample must map to a single condition in pseudobulk analysis",
          message_type = "error"
        )
      }

      counts_use <- counts[, cells_use, drop = FALSE]
      count_matrix <- aggregate_counts_by_group(
        counts = counts_use,
        groups = meta_use[[sample_col]]
      )
      condition_map <- stats::setNames(
        as.character(sample_condition[[condition_col]]),
        sample_condition[[sample_col]]
      )
      condition_use <- condition_map[colnames(count_matrix)]
      markers <- switch(
        test.use,
        limma_voom = RunDEtest_limma_voom(
          count_matrix = count_matrix,
          condition = condition_use,
          condition1 = condition1,
          condition2 = condition2,
          only.pos = only.pos,
          logfc.threshold = log(fc.threshold, base = base),
          p.adjust.method = p.adjust.method
        ),
        edgeR = RunDEtest_edgeR(
          count_matrix = count_matrix,
          condition = condition_use,
          condition1 = condition1,
          condition2 = condition2,
          only.pos = only.pos,
          logfc.threshold = log(fc.threshold, base = base),
          p.adjust.method = p.adjust.method
        )
      )
      if (is.null(markers) || nrow(markers) == 0) {
        return(NULL)
      }
      markers[, "gene"] <- rownames(markers)
      markers[, "feature"] <- rownames(markers)
      markers[, "feature_type"] <- feature_type
      markers[, "group1"] <- current_group
      markers[, "group2"] <- paste0(condition2, "_vs_", condition1)
      markers[, "condition1"] <- condition1
      markers[, "condition2"] <- condition2
      markers[, "sample_number1"] <- sum(condition_use == condition1, na.rm = TRUE)
      markers[, "sample_number2"] <- sum(condition_use == condition2, na.rm = TRUE)
      markers
    }
  )
  all_markers <- do.call(rbind.data.frame, all_markers)
  tool_name <- if (is.null(group.by)) "DEtest_pseudobulk" else paste0("DEtest_", group.by)
  if (is.null(srt@tools[[tool_name]])) {
    srt@tools[[tool_name]] <- list()
  }
  srt@tools[[tool_name]][["analysis_level"]] <- "pseudobulk"
  srt@tools[[tool_name]][["feature_type"]] <- feature_type
  srt@tools[[tool_name]][["sample_col"]] <- sample_col
  srt@tools[[tool_name]][["condition_col"]] <- condition_col
  srt@tools[[tool_name]][["condition1"]] <- condition1
  srt@tools[[tool_name]][["condition2"]] <- condition2
  srt@tools[[tool_name]][["test.use"]] <- test.use
  srt@tools[[tool_name]][["assay"]] <- assay
  srt@tools[[tool_name]][["layer"]] <- layer
  srt@tools[[tool_name]][["group.by"]] <- group.by
  srt@tools[[tool_name]][["groups_tested"]] <- target_groups

  if (is.null(all_markers) || nrow(all_markers) == 0) {
    srt@tools[[tool_name]][[paste0("AllMarkers_", test.use)]] <- data.frame()
  } else {
    rownames(all_markers) <- NULL
    all_markers[, "group1"] <- factor(all_markers[, "group1"], levels = target_groups)
    all_markers[, "test_group_number"] <- as.integer(
      table(all_markers[["gene"]])[all_markers[, "gene"]]
    )
    all_markers_matrix <- as.data.frame.matrix(
      table(all_markers[, c("gene", "group1")])
    )
    all_markers[, "test_group"] <- apply(all_markers_matrix, 1, function(x) {
      paste0(colnames(all_markers_matrix)[x > 0], collapse = ";")
    })[all_markers[, "gene"]]
    srt@tools[[tool_name]][[paste0("AllMarkers_", test.use)]] <- all_markers
  }

  log_message(
    "Pseudobulk differential testing completed",
    message_type = "success",
    verbose = verbose
  )
  srt
}

#' @title Differential gene test
#'
#' @description
#' This function utilizes the Seurat package to perform a differential expression (DE) test on gene expression data.
#' Users have the flexibility to specify custom cell groups, marker types, and various options for DE analysis.
#'
#' @md
#' @inheritParams thisutils::parallelize_fun
#' @inheritParams Seurat::FindMarkers
#' @inheritParams standard_scop
#' @inheritParams FeatureDimPlot
#' @param group.by A grouping variable in the dataset to define the groups or conditions for the differential test.
#' If not provided, the function uses the "active.ident" variable in the Seurat object.
#' @param group1 A vector of cell IDs or a character vector specifying the cells that belong to the first group.
#' If both group.by and group1 are provided, group1 takes precedence.
#' For pseudobulk analysis, this parameter is interpreted as the first condition label.
#' @param group2 A vector of cell IDs or a character vector specifying the cells that belong to the second group.
#' This parameter is only used when group.by or group1 is provided.
#' For pseudobulk analysis, this parameter is interpreted as the second condition label.
#' @param cells1 A vector of cell IDs specifying the cells that belong to group1. If provided, group1 is ignored.
#' @param cells2 A vector of cell IDs specifying the cells that belong to group2.
#' This parameter is only used when cells1 is provided.
#' @param features A vector of feature names specifying the features to consider for the differential test.
#' If not provided, all features in the dataset are considered.
#' @param feature_type Feature type used for differential testing.
#' Default is `"gene"`.
#' @param analysis_level Analysis level used for differential testing.
#' Default is `"cell"`.
#' @param markers_type A character value specifying the type of markers to find.
#' Possible values are "all", "paired", "conserved", and "disturbed".
#' Pseudobulk analysis currently supports only `"all"`.
#' @param grouping.var A character value specifying the grouping variable for finding conserved or disturbed markers.
#' This parameter is only used when markers_type is "conserved" or "disturbed".
#' @param fc.threshold A numeric value used to filter genes for testing based on their average fold change between/among the two groups.
#' Default is `1.5`.
#' @param meta.method A character value specifying the method to use for combining p-values in the conserved markers test.
#' Possible values are "maximump", "minimump", "wilkinsonp", "meanp", "sump", and "votep".
#' @param norm.method Normalization method for fold change calculation when layer is 'data'.
#' Default is `"LogNormalize"`.
#' @param sample_col Metadata column storing biological sample IDs for pseudobulk analysis.
#' Required when `analysis_level = "pseudobulk"`.
#' @param condition_col Metadata column storing condition labels for pseudobulk analysis.
#' Required when `analysis_level = "pseudobulk"`.
#' @param p.adjust.method A character value specifying the method to use for adjusting p-values.
#' Default is `"bonferroni"`.
#' @param test.use Differential testing method.
#' For pseudobulk analysis, only `"limma_voom"` and `"edgeR"` are currently supported.
#' @param ... Additional arguments to pass to the [Seurat::FindMarkers] function.
#'
#' @export
#'
#' @seealso
#' [VolcanoPlot], [RunEnrichment], [RunGSEA], [GroupHeatmap]
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunDEtest(
#'   pancreas_sub,
#'   group.by = "SubCellType"
#' )
#' AllMarkers <- dplyr::filter(
#'   pancreas_sub@tools$DEtest_SubCellType$AllMarkers_wilcox,
#'   p_val_adj < 0.05 & avg_log2FC > 1
#' )
#' ht1 <- GroupHeatmap(
#'   pancreas_sub,
#'   features = AllMarkers$gene,
#'   feature_split = AllMarkers$group1,
#'   group.by = "SubCellType"
#' )
#' ht1$plot
#'
#' TopMarkers <- AllMarkers |>
#'   dplyr::group_by(gene) |>
#'   dplyr::top_n(1, avg_log2FC) |>
#'   dplyr::group_by(group1) |>
#'   dplyr::top_n(3, avg_log2FC)
#' ht2 <- GroupHeatmap(
#'   pancreas_sub,
#'   features = TopMarkers$gene,
#'   feature_split = TopMarkers$group1,
#'   group.by = "SubCellType",
#'   show_row_names = TRUE
#' )
#' ht2$plot
#'
#' pancreas_sub <- RunDEtest(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   markers_type = "paired",
#'   cores = 2
#' )
#' PairedMarkers <- dplyr::filter(
#'   pancreas_sub@tools$DEtest_SubCellType$PairedMarkers_wilcox,
#'   p_val_adj < 0.05 & avg_log2FC > 1
#' )
#' ht3 <- GroupHeatmap(
#'   pancreas_sub,
#'   features = PairedMarkers$gene,
#'   feature_split = PairedMarkers$group1,
#'   group.by = "SubCellType"
#' )
#' ht3$plot
#'
#' data(panc8_sub)
#' panc8_sub <- integration_scop(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_method = "Uncorrected"
#' )
#' CellDimPlot(
#'   panc8_sub,
#'   group.by = c("celltype", "tech")
#' )
#'
#' panc8_sub <- RunDEtest(
#'   srt = panc8_sub,
#'   group.by = "celltype",
#'   grouping.var = "tech",
#'   markers_type = "conserved",
#'   cores = 2
#' )
#' ConservedMarkers1 <- dplyr::filter(
#'   panc8_sub@tools$DEtest_celltype$ConservedMarkers_wilcox,
#'   p_val_adj < 0.05 & avg_log2FC > 1
#' )
#' ht4 <- GroupHeatmap(
#'   panc8_sub,
#'   layer = "data",
#'   features = ConservedMarkers1$gene,
#'   feature_split = ConservedMarkers1$group1,
#'   group.by = "tech",
#'   split.by = "celltype",
#'   within_groups = TRUE
#' )
#' ht4$plot
#'
#' panc8_sub <- RunDEtest(
#'   srt = panc8_sub,
#'   group.by = "tech",
#'   grouping.var = "celltype",
#'   markers_type = "conserved",
#'   cores = 2
#' )
#' ConservedMarkers2 <- dplyr::filter(
#'   panc8_sub@tools$DEtest_tech$ConservedMarkers_wilcox,
#'   p_val_adj < 0.05 & avg_log2FC > 1
#' )
#' ht4 <- GroupHeatmap(
#'   srt = panc8_sub,
#'   layer = "data",
#'   features = ConservedMarkers2$gene,
#'   feature_split = ConservedMarkers2$group1,
#'   group.by = "tech",
#'   split.by = "celltype"
#' )
#' ht4$plot
#'
#' panc8_sub <- RunDEtest(
#'   srt = panc8_sub,
#'   group.by = "celltype",
#'   grouping.var = "tech",
#'   markers_type = "disturbed",
#'   cores = 2
#' )
#' DisturbedMarkers <- dplyr::filter(
#'   panc8_sub@tools$DEtest_celltype$DisturbedMarkers_wilcox,
#'   p_val_adj < 0.05 & avg_log2FC > 1 & var1 == "smartseq2"
#' )
#' ht5 <- GroupHeatmap(
#'   srt = panc8_sub,
#'   layer = "data",
#'   features = DisturbedMarkers$gene,
#'   feature_split = DisturbedMarkers$group1,
#'   group.by = "celltype",
#'   split.by = "tech"
#' )
#' ht5$plot
#'
#' gene_specific <- names(which(table(DisturbedMarkers$gene) == 1))
#' DisturbedMarkers_specific <- DisturbedMarkers[
#'   DisturbedMarkers$gene %in% gene_specific,
#' ]
#' ht6 <- GroupHeatmap(
#'   srt = panc8_sub,
#'   layer = "data",
#'   features = DisturbedMarkers_specific$gene,
#'   feature_split = DisturbedMarkers_specific$group1,
#'   group.by = "celltype",
#'   split.by = "tech"
#' )
#' ht6$plot
#'
#' ht7 <- GroupHeatmap(
#'   srt = panc8_sub,
#'   layer = "data",
#'   aggregate_fun = function(x) mean(expm1(x)) + 1,
#'   features = DisturbedMarkers_specific$gene,
#'   feature_split = DisturbedMarkers_specific$group1,
#'   group.by = "celltype",
#'   grouping.var = "tech",
#'   numerator = "smartseq2"
#' )
#' ht7$plot
#'
#' pbmc_small <- UpdateSeuratObject(pbmc_small)
#' pbmc_small[["sample"]] <- rep(c("S1", "S2", "S3", "S4"), length.out = ncol(pbmc_small))
#' pbmc_small[["condition"]] <- rep(c("ctrl", "ctrl", "case", "case"), length.out = ncol(pbmc_small))
#' pbmc_small <- RunDEtest(
#'   pbmc_small,
#'   analysis_level = "pseudobulk",
#'   sample_col = "sample",
#'   condition_col = "condition",
#'   test.use = "limma_voom",
#'   layer = "counts",
#'   cores = 1
#' )
#' pbmc_small <- RunDEtest(
#'   pbmc_small,
#'   analysis_level = "pseudobulk",
#'   sample_col = "sample",
#'   condition_col = "condition",
#'   test.use = "edgeR",
#'   layer = "counts",
#'   cores = 1
#' )
#' edgeR_markers <- pbmc_small@tools$DEtest_pseudobulk$AllMarkers_edgeR
#'
#' \donttest{
#' data(pbmcmultiome_sub)
#' pbmcmultiome_sub[["sample"]] <- rep(
#'   c("S1", "S2", "S3", "S4"),
#'   length.out = ncol(pbmcmultiome_sub)
#' )
#' pbmcmultiome_sub[["condition"]] <- rep(
#'   c("ctrl", "ctrl", "case", "case"),
#'   length.out = ncol(pbmcmultiome_sub)
#' )
#' pbmcmultiome_sub <- RunDEtest(
#'   pbmcmultiome_sub,
#'   assay = "peaks",
#'   layer = "counts",
#'   feature_type = "peak",
#'   analysis_level = "pseudobulk",
#'   sample_col = "sample",
#'   condition_col = "condition",
#'   test.use = "edgeR",
#'   cores = 1
#' )
#' peak_markers <- pbmcmultiome_sub@tools$DEtest_pseudobulk$AllMarkers_edgeR
#' }
RunDEtest <- function(
  srt,
  group.by = NULL,
  group1 = NULL,
  group2 = NULL,
  cells1 = NULL,
  cells2 = NULL,
  features = NULL,
  feature_type = c("gene", "peak", "cCRE"),
  analysis_level = c("cell", "pseudobulk"),
  markers_type = c(
    "all",
    "paired",
    "conserved",
    "disturbed"
  ),
  grouping.var = NULL,
  meta.method = c(
    "maximump",
    "minimump",
    "wilkinsonp",
    "meanp",
    "sump",
    "votep"
  ),
  test.use = "wilcox",
  only.pos = TRUE,
  fc.threshold = 1.5,
  base = 2,
  pseudocount.use = 1,
  mean.fxn = NULL,
  min.pct = 0.1,
  min.diff.pct = -Inf,
  max.cells.per.ident = Inf,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  norm.method = "LogNormalize",
  sample_col = NULL,
  condition_col = NULL,
  p.adjust.method = "bonferroni",
  layer = "data",
  assay = NULL,
  seed = 11,
  verbose = TRUE,
  cores = 1,
  ...
) {
  set.seed(seed)
  feature_type <- match.arg(feature_type)
  analysis_level <- match.arg(analysis_level)
  markers_type <- match.arg(markers_type)
  meta.method <- match.arg(meta.method)
  if (markers_type %in% c("conserved", "disturbed")) {
    if (is.null(grouping.var)) {
      log_message(
        "'grouping.var' must be provided when finding conserved or disturbed markers",
        message_type = "error"
      )
    }
  }
  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(srt)
  }

  if (analysis_level == "pseudobulk") {
    if (!is.null(cells1) || !is.null(cells2)) {
      log_message(
        "{.arg cells1} and {.arg cells2} are not supported when {.arg analysis_level = 'pseudobulk'}. Use {.arg group1} and {.arg group2} to specify condition labels.",
        message_type = "error"
      )
    }
    if (!is.null(grouping.var)) {
      log_message(
        "{.arg grouping.var} is not supported when {.arg analysis_level = 'pseudobulk'}.",
        message_type = "error"
      )
    }
    if (!identical(layer, "counts")) {
      log_message(
        "Pseudobulk differential testing uses the {.arg counts} layer. Reset {.arg layer = 'counts'}.",
        message_type = "warning",
        verbose = verbose
      )
      layer <- "counts"
    }
    return(RunDEtest_pseudobulk(
      srt = srt,
      group.by = group.by,
      group1 = group1,
      group2 = group2,
      features = features,
      feature_type = feature_type,
      markers_type = markers_type,
      test.use = test.use,
      only.pos = only.pos,
      fc.threshold = fc.threshold,
      base = base,
      sample_col = sample_col,
      condition_col = condition_col,
      p.adjust.method = p.adjust.method,
      layer = layer,
      assay = assay,
      verbose = verbose,
      cores = cores,
      ...
    ))
  }

  check_r("immunogenomics/presto", verbose = FALSE)

  status <- CheckDataType(srt, layer = layer, assay = assay)
  if (layer == "counts" && status != "raw_counts") {
    log_message(
      "Data in the {.arg counts} layer is not raw counts",
      message_type = "error"
    )
  }
  if (layer == "data" && status != "log_normalized_counts") {
    if (status == "raw_counts") {
      log_message(
        "Data in the {.arg data} layer is raw counts. Perform {.fun NormalizeData}({.val LogNormalize})",
        message_type = "warning",
        verbose = verbose
      )
      srt <- NormalizeData(
        object = srt,
        assay = assay,
        normalization.method = "LogNormalize",
        verbose = FALSE
      )
    }
    if (status == "raw_normalized_counts") {
      log_message(
        "Data in the {.arg data} layer is raw_normalized_counts. Perform {.fun NormalizeData}({.val LogNormalize})",
        message_type = "warning",
        verbose = verbose
      )
      srt <- NormalizeData(
        object = srt,
        assay = assay,
        normalization.method = "LogNormalize",
        verbose = FALSE
      )
    }
    if (status == "unknown") {
      log_message(
        "Data in the 'data' layer is unknown. Please check the data type",
        message_type = "warning",
        verbose = verbose
      )
    }
  }

  log_message(
    "Start differential expression test",
    verbose = verbose
  )
  if (fc.threshold < 1) {
    log_message(
      "{.arg fc.threshold} must be greater than or equal to 1",
      message_type = "error"
    )
  }

  if (!is.null(cells1) || !is.null(group1)) {
    if (is.null(cells1)) {
      if (is.null(group.by)) {
        log_message(
          "{.arg group.by} must be provided when {.arg group1} specified",
          message_type = "error"
        )
      }
      cells1 <- colnames(srt)[srt[[group.by, drop = TRUE]] %in% group1]
    }
    if (is.null(cells2) && !is.null(group2)) {
      cells2 <- colnames(srt)[srt[[group.by, drop = TRUE]] %in% group2]
    }
    if (!all(cells1 %in% colnames(srt))) {
      log_message(
        "{.arg cells1} has some cells not in {.cls Seurat}",
        message_type = "error"
      )
    }
    if (is.null(cells2)) {
      cells2 <- colnames(srt)[!colnames(srt) %in% cells1]
      group2 <- "others"
    }
    if (!all(cells2 %in% colnames(srt))) {
      log_message(
        "{.arg cells2} has some cells not in {.cls Seurat}",
        message_type = "error"
      )
    }
    if (length(cells1) < 3 || length(cells2) < 3) {
      log_message(
        "Cell groups must have more than 3 cells",
        message_type = "error"
      )
    }

    log_message(
      "Find ",
      markers_type,
      " markers(",
      test.use,
      ") for custom cell groups...",
      verbose = verbose
    )

    if (markers_type == "all") {
      markers <- Seurat::FindMarkers(
        object = Seurat::GetAssay(srt, assay),
        layer = layer,
        cells.1 = cells1,
        cells.2 = cells2,
        features = features,
        test.use = test.use,
        logfc.threshold = log(fc.threshold, base = base),
        base = base,
        min.pct = min.pct,
        min.diff.pct = min.diff.pct,
        max.cells.per.ident = max.cells.per.ident,
        min.cells.feature = min.cells.feature,
        min.cells.group = min.cells.group,
        latent.vars = latent.vars,
        only.pos = only.pos,
        norm.method = norm.method,
        pseudocount.use = pseudocount.use,
        mean.fxn = mean.fxn,
        verbose = FALSE,
        ...
      )

      if (!is.null(markers) && nrow(markers) > 0) {
        markers[, "gene"] <- rownames(markers)
        group1_str <- if (is.null(group1)) {
          "group1"
        } else if (length(group1) > 1) {
          paste(group1, collapse = ";")
        } else {
          as.character(group1)
        }
        group2_str <- if (is.null(group2)) {
          "group2"
        } else if (length(group2) > 1) {
          paste(group2, collapse = ";")
        } else {
          as.character(group2)
        }
        markers[, "group1"] <- group1_str
        markers[, "group2"] <- group2_str
        rownames(markers) <- NULL
        markers[, "group1"] <- factor(
          markers[, "group1"],
          levels = unique(markers[, "group1"])
        )
        if ("p_val" %in% colnames(markers)) {
          markers[, "p_val_adj"] <- stats::p.adjust(
            markers[, "p_val"],
            method = p.adjust.method
          )
        }
        markers[, "test_group_number"] <- as.integer(
          table(markers[["gene"]])[markers[, "gene"]]
        )
        markers_matrix <- as.data.frame.matrix(
          table(markers[, c("gene", "group1")])
        )
        markers[, "test_group"] <- apply(markers_matrix, 1, function(x) {
          paste0(colnames(markers_matrix)[x > 0], collapse = ";")
        })[markers[, "gene"]]
        srt@tools[["DEtest_custom"]][[paste0(
          "AllMarkers_",
          test.use
        )]] <- markers
        srt@tools[["DEtest_custom"]][["cells1"]] <- cells1
        srt@tools[["DEtest_custom"]][["cells2"]] <- cells2
      } else {
        log_message(
          "No markers found.",
          message_type = "warning",
          verbose = verbose
        )
      }
    }

    if (markers_type == "conserved") {
      markers <- FindConservedMarkers2(
        object = srt,
        assay = assay,
        layer = layer,
        cells.1 = cells1,
        cells.2 = cells2,
        features = features,
        grouping.var = grouping.var,
        test.use = test.use,
        logfc.threshold = log(fc.threshold, base = base),
        base = base,
        min.pct = min.pct,
        min.diff.pct = min.diff.pct,
        max.cells.per.ident = max.cells.per.ident,
        min.cells.feature = min.cells.feature,
        min.cells.group = min.cells.group,
        latent.vars = latent.vars,
        only.pos = only.pos,
        norm.method = norm.method,
        meta.method = meta.method,
        pseudocount.use = pseudocount.use,
        mean.fxn = mean.fxn,
        verbose = FALSE,
        ...
      )
      if (!is.null(markers) && nrow(markers) > 0) {
        markers[, "gene"] <- rownames(markers)
        group1_str <- if (is.null(group1)) {
          "group1"
        } else if (length(group1) > 1) {
          paste(group1, collapse = ";")
        } else {
          as.character(group1)
        }
        group2_str <- if (is.null(group2)) {
          "group2"
        } else if (length(group2) > 1) {
          paste(group2, collapse = ";")
        } else {
          as.character(group2)
        }
        markers[, "group1"] <- group1_str
        markers[, "group2"] <- group2_str
        rownames(markers) <- NULL
        markers[, "group1"] <- factor(
          markers[, "group1"],
          levels = unique(markers[, "group1"])
        )
        if ("p_val" %in% colnames(markers)) {
          markers[, "p_val_adj"] <- stats::p.adjust(
            markers[, "p_val"],
            method = p.adjust.method
          )
        }
        markers[, "test_group_number"] <- as.integer(
          table(markers[["gene"]])[markers[, "gene"]]
        )
        markers_matrix <- as.data.frame.matrix(
          table(markers[, c("gene", "group1")])
        )
        markers[, "test_group"] <- apply(markers_matrix, 1, function(x) {
          paste0(colnames(markers_matrix)[x > 0], collapse = ";")
        })[markers[, "gene"]]
        srt@tools[["DEtest_custom"]][[paste0(
          "ConservedMarkers_",
          test.use
        )]] <- markers
        srt@tools[["DEtest_custom"]][["cells1"]] <- cells1
        srt@tools[["DEtest_custom"]][["cells2"]] <- cells2
      } else {
        log_message(
          "No markers found.",
          message_type = "warning",
          verbose = verbose
        )
      }
    }
    if (markers_type == "disturbed") {
      srt_tmp <- srt
      srt_tmp[[grouping.var, drop = TRUE]][setdiff(
        colnames(srt_tmp),
        cells1
      )] <- NA
      srt_tmp <- RunDEtest(
        srt = srt_tmp,
        assay = assay,
        layer = layer,
        group.by = grouping.var,
        markers_type = "all",
        features = features,
        test.use = test.use,
        fc.threshold = fc.threshold,
        base = base,
        min.pct = min.pct,
        min.diff.pct = min.diff.pct,
        max.cells.per.ident = max.cells.per.ident,
        min.cells.feature = min.cells.feature,
        min.cells.group = min.cells.group,
        latent.vars = latent.vars,
        only.pos = only.pos,
        norm.method = norm.method,
        p.adjust.method = p.adjust.method,
        pseudocount.use = pseudocount.use,
        mean.fxn = mean.fxn,
        cores = cores,
        seed = seed,
        verbose = FALSE,
        ...
      )
      markers <- srt_tmp@tools[[paste0("DEtest_", grouping.var)]][[paste0(
        "AllMarkers_",
        test.use
      )]]
      if (!is.null(markers) && nrow(markers) > 0) {
        colnames(markers) <- gsub("group", "var", colnames(markers))
        markers[["group1"]] <- group1 %||% "group1"
        srt@tools[["DEtest_custom"]][[paste0(
          "DisturbedMarkers_",
          test.use
        )]] <- markers
        srt@tools[["DEtest_custom"]][["cells1"]] <- cells1
      } else {
        log_message(
          "No markers found",
          message_type = "warning",
          verbose = verbose
        )
      }
    }
  } else {
    if (is.null(group.by)) {
      cell_group <- Seurat::Idents(srt)
      group.by <- "active.ident"
    } else {
      cell_group <- srt[[group.by, drop = TRUE]]
    }
    if (!is.factor(cell_group)) {
      cell_group <- factor(cell_group, levels = unique(cell_group))
    }
    cell_group <- lapply(levels(cell_group), function(x) {
      cell <- cell_group[cell_group == x]
      out <- sample(
        cell,
        size = min(max.cells.per.ident, length(cell)),
        replace = FALSE
      )
      return(out)
    })
    cell_group <- stats::setNames(
      unlist(
        lapply(cell_group, function(x) x),
        use.names = FALSE
      ),
      unlist(lapply(cell_group, names))
    )

    args1 <- list(
      object = Seurat::GetAssay(srt, assay),
      layer = layer,
      features = features,
      test.use = test.use,
      logfc.threshold = log(fc.threshold, base = base),
      base = base,
      min.pct = min.pct,
      min.diff.pct = min.diff.pct,
      max.cells.per.ident = max.cells.per.ident,
      min.cells.feature = min.cells.feature,
      min.cells.group = min.cells.group,
      latent.vars = latent.vars,
      only.pos = only.pos,
      norm.method = norm.method,
      pseudocount.use = pseudocount.use,
      mean.fxn = mean.fxn,
      verbose = FALSE,
      ...
    )


    log_message(
      "Find ",
      markers_type,
      " markers(",
      test.use,
      ") among ",
      nlevels(cell_group),
      " groups..."
    )

    if (markers_type == "all") {
      AllMarkers <- parallelize_fun(
        levels(cell_group),
        function(group) {
          cells.1 <- names(cell_group)[which(cell_group == group)]
          cells.2 <- names(cell_group)[which(cell_group != group)]
          if (length(cells.1) < 3 || length(cells.2) < 3) {
            return(NULL)
          } else {
            args1[["cells.1"]] <- cells.1
            args1[["cells.2"]] <- cells.2
            markers <- do.call(FindMarkers, args1)
            if (!is.null(markers) && nrow(markers) > 0) {
              markers[, "gene"] <- rownames(markers)
              markers[, "group1"] <- as.character(group)
              markers[, "group2"] <- "others"
              if ("p_val" %in% colnames(markers)) {
                markers[, "p_val_adj"] <- stats::p.adjust(
                  markers[, "p_val"],
                  method = p.adjust.method
                )
              }
              return(markers)
            } else {
              return(NULL)
            }
          }
        },
        cores = cores,
        verbose = verbose
      )
      AllMarkers <- do.call(rbind.data.frame, AllMarkers)
      if (!is.null(AllMarkers) && nrow(AllMarkers) > 0) {
        rownames(AllMarkers) <- NULL
        AllMarkers[, "group1"] <- factor(
          AllMarkers[, "group1"],
          levels = levels(cell_group)
        )
        AllMarkers[, "test_group_number"] <- as.integer(
          table(AllMarkers[["gene"]])[AllMarkers[, "gene"]]
        )
        AllMarkersMatrix <- as.data.frame.matrix(
          table(AllMarkers[, c("gene", "group1")])
        )
        AllMarkers[, "test_group"] <- apply(
          AllMarkersMatrix,
          1,
          function(x) {
            paste0(colnames(AllMarkersMatrix)[x > 0], collapse = ";")
          }
        )[AllMarkers[, "gene"]]
        srt@tools[[paste0("DEtest_", group.by)]][[paste0(
          "AllMarkers_",
          test.use
        )]] <- AllMarkers
      } else {
        srt@tools[[paste0("DEtest_", group.by)]][[paste0(
          "AllMarkers_",
          test.use
        )]] <- data.frame()
      }
    }

    if (markers_type == "paired") {
      pair <- expand.grid(x = levels(cell_group), y = levels(cell_group))
      pair <- pair[pair[, 1] != pair[, 2], , drop = FALSE]
      PairedMarkers <- parallelize_fun(
        seq_len(nrow(pair)),
        function(i) {
          cells.1 <- names(cell_group)[which(cell_group == pair[i, 1])]
          cells.2 <- names(cell_group)[which(cell_group == pair[i, 2])]
          if (length(cells.1) < 3 || length(cells.2) < 3) {
            return(NULL)
          } else {
            args1[["cells.1"]] <- cells.1
            args1[["cells.2"]] <- cells.2
            markers <- do.call(FindMarkers, args1)
            if (!is.null(markers) && nrow(markers) > 0) {
              markers[, "gene"] <- rownames(markers)
              markers[, "group1"] <- as.character(pair[i, 1])
              markers[, "group2"] <- as.character(pair[i, 2])
              if ("p_val" %in% colnames(markers)) {
                markers[, "p_val_adj"] <- stats::p.adjust(
                  markers[, "p_val"],
                  method = p.adjust.method
                )
              }
              return(markers)
            } else {
              return(NULL)
            }
          }
        },
        cores = cores,
        verbose = verbose
      )
      PairedMarkers <- do.call(rbind.data.frame, PairedMarkers)
      if (!is.null(PairedMarkers) && nrow(PairedMarkers) > 0) {
        rownames(PairedMarkers) <- NULL
        PairedMarkers[, "group1"] <- factor(
          PairedMarkers[, "group1"],
          levels = levels(cell_group)
        )
        PairedMarkers[, "test_group_number"] <- as.integer(
          table(PairedMarkers[["gene"]])[PairedMarkers[, "gene"]]
        )
        PairedMarkersMatrix <- as.data.frame.matrix(
          table(
            PairedMarkers[, c("gene", "group1")]
          )
        )
        PairedMarkers[, "test_group"] <- apply(
          PairedMarkersMatrix,
          1,
          function(x) {
            paste0(colnames(PairedMarkersMatrix)[x > 0], collapse = ";")
          }
        )[PairedMarkers[, "gene"]]
        srt@tools[[paste0("DEtest_", group.by)]][[paste0(
          "PairedMarkers_",
          test.use
        )]] <- PairedMarkers
        srt@tools[[paste0("DEtest_", group.by)]][[paste0(
          "PairedMarkersMatrix_",
          test.use
        )]] <- PairedMarkersMatrix
      } else {
        log_message(
          "No markers found",
          message_type = "warning",
          verbose = verbose
        )
        srt@tools[[paste0("DEtest_", group.by)]][[paste0(
          "PairedMarkers_",
          test.use
        )]] <- data.frame()
        srt@tools[[paste0("DEtest_", group.by)]][[paste0(
          "PairedMarkersMatrix_",
          test.use
        )]] <- NULL
      }
    }

    if (markers_type == "conserved") {
      ConservedMarkers <- parallelize_fun(
        levels(cell_group),
        function(group) {
          cells.1 <- names(cell_group)[which(cell_group == group)]
          cells.2 <- names(cell_group)[which(cell_group != group)]
          if (length(cells.1) < 3 || length(cells.2) < 3) {
            return(NULL)
          } else {
            args1[["cells.1"]] <- cells.1
            args1[["cells.2"]] <- cells.2
            args1[["object"]] <- srt
            args1[["assay"]] <- assay
            args1[["grouping.var"]] <- grouping.var
            args1[["meta.method"]] <- meta.method
            markers <- do.call(FindConservedMarkers2, args1)
            if (!is.null(markers) && nrow(markers) > 0) {
              markers[, "gene"] <- rownames(markers)
              markers[, "group1"] <- as.character(group)
              markers[, "group2"] <- "others"
              if ("p_val" %in% colnames(markers)) {
                markers[, "p_val_adj"] <- stats::p.adjust(
                  markers[, "p_val"],
                  method = p.adjust.method
                )
              }
              return(markers)
            } else {
              return(NULL)
            }
          }
        },
        cores = cores,
        verbose = verbose
      )
      ConservedMarkers <- do.call(
        rbind.data.frame,
        lapply(
          ConservedMarkers,
          function(x) {
            x[, c(
              "avg_log2FC",
              "pct.1",
              "pct.2",
              "max_pval",
              "p_val",
              "p_val_adj",
              "gene",
              "group1",
              "group2"
            )]
          }
        )
      )
      if (!is.null(ConservedMarkers) && nrow(ConservedMarkers) > 0) {
        rownames(ConservedMarkers) <- NULL
        ConservedMarkers[, "group1"] <- factor(
          ConservedMarkers[, "group1"],
          levels = levels(cell_group)
        )
        ConservedMarkers[
          ,
          "test_group_number"
        ] <- as.integer(table(ConservedMarkers[["gene"]])[ConservedMarkers[
          ,
          "gene"
        ]])
        ConservedMarkersMatrix <- as.data.frame.matrix(table(ConservedMarkers[, c(
          "gene",
          "group1"
        )]))
        ConservedMarkers[, "test_group"] <- apply(
          ConservedMarkersMatrix,
          1,
          function(x) {
            paste0(colnames(ConservedMarkersMatrix)[x > 0], collapse = ";")
          }
        )[ConservedMarkers[, "gene"]]
        ConservedMarkers <- ConservedMarkers[, c(
          "avg_log2FC",
          "pct.1",
          "pct.2",
          "max_pval",
          "p_val",
          "p_val_adj",
          "gene",
          "group1",
          "group2",
          "test_group_number",
          "test_group"
        )]
        srt@tools[[paste0("DEtest_", group.by)]][[paste0(
          "ConservedMarkers_",
          test.use
        )]] <- ConservedMarkers
      } else {
        log_message(
          "No markers found",
          message_type = "warning",
          verbose = verbose
        )
        srt@tools[[paste0("DEtest_", group.by)]][[paste0(
          "ConservedMarkers_",
          test.use
        )]] <- data.frame()
      }
    }

    if (markers_type == "disturbed") {
      DisturbedMarkers <- parallelize_fun(
        levels(cell_group),
        function(group) {
          cells.1 <- names(cell_group)[which(cell_group == group)]
          srt_tmp <- srt
          srt_tmp[[grouping.var, drop = TRUE]][setdiff(
            colnames(srt_tmp),
            cells.1
          )] <- NA
          if (
            length(stats::na.omit(unique(srt_tmp[[grouping.var, drop = TRUE]]))) < 2
          ) {
            return(NULL)
          } else {
            srt_tmp <- RunDEtest(
              srt = srt_tmp,
              assay = assay,
              layer = layer,
              group.by = grouping.var,
              markers_type = "all",
              features = features,
              test.use = test.use,
              fc.threshold = fc.threshold,
              base = base,
              min.pct = min.pct,
              min.diff.pct = min.diff.pct,
              max.cells.per.ident = max.cells.per.ident,
              min.cells.feature = min.cells.feature,
              min.cells.group = min.cells.group,
              latent.vars = latent.vars,
              only.pos = only.pos,
              norm.method = norm.method,
              p.adjust.method = p.adjust.method,
              pseudocount.use = pseudocount.use,
              mean.fxn = mean.fxn,
              cores = cores,
              seed = seed,
              verbose = FALSE,
              ...
            )
            markers <- srt_tmp@tools[[paste0("DEtest_", grouping.var)]][[paste0(
              "AllMarkers_",
              test.use
            )]]
            if (!is.null(markers) && nrow(markers) > 0) {
              colnames(markers) <- gsub("group", "var", colnames(markers))
              markers[["group1"]] <- as.character(group)
              return(markers)
            } else {
              return(NULL)
            }
          }
        },
        cores = cores,
        verbose = verbose
      )

      DisturbedMarkers <- do.call(rbind.data.frame, DisturbedMarkers)
      if (!is.null(DisturbedMarkers) && nrow(DisturbedMarkers) > 0) {
        rownames(DisturbedMarkers) <- NULL
        DisturbedMarkers[, "group1"] <- factor(
          DisturbedMarkers[, "group1"],
          levels = levels(cell_group)
        )
        DisturbedMarkers[
          ,
          "test_group_number"
        ] <- as.integer(table(unique(DisturbedMarkers[, c("gene", "group1")])[[
          "gene"
        ]])[DisturbedMarkers[, "gene"]])
        DisturbedMarkersMatrix <- as.data.frame.matrix(table(DisturbedMarkers[, c(
          "gene",
          "group1"
        )]))
        DisturbedMarkers[, "test_group"] <- apply(
          DisturbedMarkersMatrix,
          1,
          function(x) {
            paste0(colnames(DisturbedMarkersMatrix)[x > 0], collapse = ";")
          }
        )[DisturbedMarkers[, "gene"]]
        srt@tools[[paste0("DEtest_", group.by)]][[paste0(
          "DisturbedMarkers_",
          test.use
        )]] <- DisturbedMarkers
      } else {
        log_message(
          "No markers found.",
          message_type = "warning",
          verbose = verbose
        )
        srt@tools[[paste0("DEtest_", group.by)]][[paste0(
          "DisturbedMarkers_",
          test.use
        )]] <- data.frame()
      }
    }
  }

  log_message(
    "Differential expression test completed",
    message_type = "success",
    verbose = verbose
  )

  return(srt)
}
