#' @export
FoldChange.default <- function(
    object,
    cells.1,
    cells.2,
    mean.fxn,
    fc.name,
    features = NULL,
    ...) {
  features <- features %||% rownames(x = object)
  thresh.min <- 0
  pct.1 <- round(
    x = Matrix::rowSums(
      x = object[features, cells.1, drop = FALSE] > thresh.min,
      na.rm = TRUE
    ) /
      length(x = cells.1),
    digits = 3
  )
  pct.2 <- round(
    x = Matrix::rowSums(
      x = object[features, cells.2, drop = FALSE] > thresh.min,
      na.rm = TRUE
    ) /
      length(x = cells.2),
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
    ...) {
  if (
    !(test.use %in% c("negbinom", "poisson", "MAST", "LR")) &&
      !is.null(x = latent.vars)
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
  data.use <- Matrix::as.matrix(data.use)

  de.results <- switch(
    EXPR = test.use,
    "wilcox" = WilcoxDETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose,
      ...
    ),
    "bimod" = Seurat:::DiffExpTest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    "roc" = Seurat:::MarkerTest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    "t" = Seurat:::DiffTTest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    "negbinom" = Seurat:::GLMDETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      min.cells = min.cells.feature,
      latent.vars = latent.vars,
      test.use = test.use,
      verbose = verbose
    ),
    "poisson" = Seurat:::GLMDETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      min.cells = min.cells.feature,
      latent.vars = latent.vars,
      test.use = test.use,
      verbose = verbose
    ),
    "MAST" = Seurat:::MASTDETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      latent.vars = latent.vars,
      verbose = verbose,
      ...
    ),
    "DESeq2" = Seurat:::DESeq2DETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose,
      ...
    ),
    "LR" = Seurat:::LRDETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      latent.vars = latent.vars,
      verbose = verbose
    ),
    log_message(
      "Unknown test: ", test.use,
      message_type = "error"
    )
  )
  return(de.results)
}

WilcoxDETest <- function(
    data.use,
    cells.1,
    cells.2,
    verbose = TRUE,
    ...) {
  data.use <- data.use[, c(cells.1, cells.2), drop = FALSE]
  j <- seq_len(length.out = length(x = cells.1))
  my.sapply <- ifelse(
    test = verbose && future::nbrOfWorkers() == 1,
    yes = pbapply::pbsapply,
    no = future.apply::future_sapply
  )
  check_r("limma")
  p_val <- my.sapply(
    X = 1:nrow(x = data.use),
    FUN = function(x) {
      keep <- colnames(data.use)[!is.na(data.use[x, ])]
      j <- seq_len(length.out = length(x = intersect(cells.1, keep)))
      statistics <- data.use[x, keep]
      return(
        min(
          2 * min(
            limma::rankSumTestWithCorrelation(
              index = j,
              statistics = statistics
            )
          ),
          1
        )
      )
    }
  )

  return(data.frame(p_val, row.names = rownames(x = data.use)))
}


#' Differential gene test
#'
#' This function utilizes the Seurat package to perform a differential expression (DE) test on gene expression data.
#' Users have the flexibility to specify custom cell groups, marker types, and various options for DE analysis.
#'
#' @md
#' @inheritParams Seurat::FindMarkers
#' @param srt A Seurat object.
#' @param layer The layer used.
#' @param group_by A grouping variable in the dataset to define the groups or conditions for the differential test.
#' If not provided, the function uses the "active.ident" variable in the Seurat object.
#' @param group1 A vector of cell IDs or a character vector specifying the cells that belong to the first group.
#' If both group_by and group1 are provided, group1 takes precedence.
#' @param group2 A vector of cell IDs or a character vector specifying the cells that belong to the second group.
#' This parameter is only used when group_by or group1 is provided.
#' @param cells1 A vector of cell IDs specifying the cells that belong to group1. If provided, group1 is ignored.
#' @param cells2 A vector of cell IDs specifying the cells that belong to group2.
#' This parameter is only used when cells1 is provided.
#' @param features A vector of gene names specifying the features to consider for the differential test.
#' If not provided, all features in the dataset are considered.
#' @param markers_type A character value specifying the type of markers to find.
#' Possible values are "all", "paired", "conserved", and "disturbed".
#' @param grouping.var A character value specifying the grouping variable for finding conserved or disturbed markers.
#' This parameter is only used when markers_type is "conserved" or "disturbed".
#' @param fc.threshold A numeric value used to filter genes for testing based on their average fold change between/among the two groups.
#' Default is 1.5.
#' @param meta.method A character value specifying the method to use for combining p-values in the conserved markers test.
#' Possible values are "maximump", "minimump", "wilkinsonp", "meanp", "sump", and "votep".
#' @param norm.method Normalization method for fold change calculation when layer is 'data'.
#' Default is "LogNormalize".
#' @param p.adjust.method A character value specifying the method to use for adjusting p-values.
#' Default is "bonferroni".
#' @param BPPARAM A BiocParallelParam object specifying the parallelization parameters for the differential test.
#' Default is [BiocParallel::bpparam()].
#' @param seed An integer value specifying the seed.
#' Default is 11.
#' @param verbose A logical value indicating whether to display progress messages during the differential test.
#' Default is TRUE.
#' @param ... Additional arguments to pass to the [Seurat::FindMarkers] function.
#'
#' @export
#'
#' @seealso \code{\link{RunEnrichment}} \code{\link{RunGSEA}} \code{\link{GroupHeatmap}}
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- RunDEtest(
#'   srt = pancreas_sub,
#'   group_by = "SubCellType"
#' )
#' AllMarkers <- dplyr::filter(
#'   pancreas_sub@tools$DEtest_SubCellType$AllMarkers_wilcox,
#'   p_val_adj < 0.05 & avg_log2FC > 1
#' )
#' table(AllMarkers$group1)
#' ht1 <- GroupHeatmap(
#'   srt = pancreas_sub,
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
#'   srt = pancreas_sub,
#'   features = TopMarkers$gene,
#'   feature_split = TopMarkers$group1,
#'   group.by = "SubCellType",
#'   show_row_names = TRUE
#' )
#' ht2$plot
#'
#' pancreas_sub <- RunDEtest(
#'   srt = pancreas_sub,
#'   group_by = "SubCellType",
#'   markers_type = "paired"
#' )
#' PairedMarkers <- dplyr::filter(
#'   pancreas_sub@tools$DEtest_SubCellType$PairedMarkers_wilcox,
#'   p_val_adj < 0.05 & avg_log2FC > 1
#' )
#' ht3 <- GroupHeatmap(
#'   srt = pancreas_sub,
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
#'   integration_method = "Seurat"
#' )
#' CellDimPlot(
#'   panc8_sub,
#'   group.by = c("celltype", "tech")
#' )
#'
#' panc8_sub <- RunDEtest(
#'   srt = panc8_sub,
#'   group_by = "celltype",
#'   grouping.var = "tech",
#'   markers_type = "conserved"
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
#'   group_by = "tech",
#'   grouping.var = "celltype",
#'   markers_type = "conserved"
#' )
#' ConservedMarkers2 <- dplyr::filter(
#'   panc8_sub@tools$DEtest_tech$ConservedMarkers_wilcox,
#'   p_val_adj < 0.05 & avg_log2FC > 1
#' )
#' table(ConservedMarkers2$group1)
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
#'   group_by = "celltype",
#'   grouping.var = "tech",
#'   markers_type = "disturbed"
#' )
#' DisturbedMarkers <- dplyr::filter(
#'   panc8_sub@tools$DEtest_celltype$DisturbedMarkers_wilcox,
#'   p_val_adj < 0.05 & avg_log2FC > 1 & var1 == "smartseq2"
#' )
#' table(DisturbedMarkers$group1)
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
#' table(DisturbedMarkers_specific$group1)
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
#' }
RunDEtest <- function(
    srt,
    group_by = NULL,
    group1 = NULL,
    group2 = NULL,
    cells1 = NULL,
    cells2 = NULL,
    features = NULL,
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
    p.adjust.method = "bonferroni",
    layer = "data",
    assay = NULL,
    BPPARAM = BiocParallel::bpparam(),
    seed = 11,
    verbose = TRUE,
    ...) {
  set.seed(seed)
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
  assay <- assay %||% DefaultAssay(srt)

  status <- check_data_type(srt, layer = layer, assay = assay)
  if (layer == "counts" && status != "raw_counts") {
    log_message(
      "Data in the 'counts' layer is not raw counts.",
      message_type = "error"
    )
  }
  if (layer == "data" && status != "log_normalized_counts") {
    if (status == "raw_counts") {
      log_message(
        "Data in the 'data' layer is raw counts. Perform NormalizeData(LogNormalize) on the data.",
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
        "Data in the 'data' layer is raw_normalized_counts. Perform NormalizeData(LogNormalize) on the data.",
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
        "Data in the 'data' layer is unknown. Please check the data type.",
        message_type = "warning",
        verbose = verbose
      )
    }
  }
  BiocParallel::bpprogressbar(BPPARAM) <- FALSE
  BiocParallel::bpRNGseed(BPPARAM) <- seed

  time_start <- Sys.time()

  log_message("Start DEtest", verbose = verbose)
  log_message("Workers: ", bpworkers(BPPARAM), verbose = verbose)


  if (fc.threshold < 1) {
    log_message(
      "fc.threshold must be greater than or equal to 1",
      message_type = "error"
    )
  }

  if (!is.null(cells1) || !is.null(group1)) {
    if (is.null(cells1)) {
      if (is.null(group_by)) {
        log_message(
          "'group_by' must be provided when 'group1' specified",
          message_type = "error"
        )
      }
      cells1 <- colnames(srt)[srt[[group_by, drop = TRUE]] %in% group1]
    }
    if (is.null(cells2) && !is.null(group2)) {
      cells2 <- colnames(srt)[srt[[group_by, drop = TRUE]] %in% group2]
    }
    if (!all(cells1 %in% colnames(srt))) {
      log_message(
        "cells1 has some cells not in the Seurat object.",
        message_type = "error"
      )
    }
    if (is.null(cells2)) {
      cells2 <- colnames(srt)[!colnames(srt) %in% cells1]
      group2 <- "others"
    }
    if (!all(cells2 %in% colnames(srt))) {
      log_message(
        "cells2 has some cells not in the Seurat object.",
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
        object = SeuratObject::Assays(srt, assay),
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
        markers[, "group1"] <- group1 %||% "group1"
        markers[, "group2"] <- group2 %||% "group2"
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
        MarkersMatrix <- as.data.frame.matrix(
          table(markers[, c("gene", "group1")])
        )
        markers[, "test_group"] <- apply(MarkersMatrix, 1, function(x) {
          paste0(colnames(MarkersMatrix)[x > 0], collapse = ";")
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
        markers[, "group1"] <- group1 %||% "group1"
        markers[, "group2"] <- group2 %||% "group2"
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
        MarkersMatrix <- as.data.frame.matrix(
          table(markers[, c("gene", "group1")])
        )
        markers[, "test_group"] <- apply(MarkersMatrix, 1, function(x) {
          paste0(colnames(MarkersMatrix)[x > 0], collapse = ";")
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
      bpprogressbar(BPPARAM) <- FALSE
      srt_tmp <- RunDEtest(
        srt = srt_tmp,
        assay = assay,
        layer = layer,
        group_by = grouping.var,
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
        BPPARAM = BPPARAM,
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
          "No markers found.",
          message_type = "warning",
          verbose = verbose
        )
      }
    }
  } else {
    if (is.null(group_by)) {
      cell_group <- Seurat::Idents(srt)
      group_by <- "active.ident"
    } else {
      cell_group <- srt[[group_by, drop = TRUE]]
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
      object = SeuratObject::Assays(srt, assay),
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
      AllMarkers <- BiocParallel::bplapply(
        levels(cell_group),
        FUN = function(group) {
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
        BPPARAM = BPPARAM
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
        srt@tools[[paste0("DEtest_", group_by)]][[paste0(
          "AllMarkers_",
          test.use
        )]] <- AllMarkers
      } else {
        srt@tools[[paste0("DEtest_", group_by)]][[paste0(
          "AllMarkers_",
          test.use
        )]] <- data.frame()
      }
    }

    if (markers_type == "paired") {
      pair <- expand.grid(x = levels(cell_group), y = levels(cell_group))
      pair <- pair[pair[, 1] != pair[, 2], , drop = FALSE]
      PairedMarkers <- BiocParallel::bplapply(
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
        BPPARAM = BPPARAM
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
        srt@tools[[paste0("DEtest_", group_by)]][[paste0(
          "PairedMarkers_",
          test.use
        )]] <- PairedMarkers
        srt@tools[[paste0("DEtest_", group_by)]][[paste0(
          "PairedMarkersMatrix_",
          test.use
        )]] <- PairedMarkersMatrix
      } else {
        log_message(
          "No markers found.",
          message_type = "warning",
          verbose = verbose
        )
        srt@tools[[paste0("DEtest_", group_by)]][[paste0(
          "PairedMarkers_",
          test.use
        )]] <- data.frame()
        srt@tools[[paste0("DEtest_", group_by)]][[paste0(
          "PairedMarkersMatrix_",
          test.use
        )]] <- NULL
      }
    }

    if (markers_type == "conserved") {
      ConservedMarkers <- BiocParallel::bplapply(
        levels(cell_group),
        FUN = function(group) {
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
        BPPARAM = BPPARAM
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
        srt@tools[[paste0("DEtest_", group_by)]][[paste0(
          "ConservedMarkers_",
          test.use
        )]] <- ConservedMarkers
      } else {
        log_message(
          "No markers found.",
          message_type = "warning",
          verbose = verbose
        )
        srt@tools[[paste0("DEtest_", group_by)]][[paste0(
          "ConservedMarkers_",
          test.use
        )]] <- data.frame()
      }
    }

    if (markers_type == "disturbed") {
      sub_BPPARAM <- BiocParallel::SerialParam()
      BiocParallel::bpprogressbar(sub_BPPARAM) <- FALSE
      DisturbedMarkers <- BiocParallel::bplapply(
        levels(cell_group),
        FUN = function(group) {
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
              group_by = grouping.var,
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
              BPPARAM = sub_BPPARAM,
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
        BPPARAM = BPPARAM
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
        srt@tools[[paste0("DEtest_", group_by)]][[paste0(
          "DisturbedMarkers_",
          test.use
        )]] <- DisturbedMarkers
      } else {
        log_message(
          "No markers found.",
          message_type = "warning",
          verbose = verbose
        )
        srt@tools[[paste0("DEtest_", group_by)]][[paste0(
          "DisturbedMarkers_",
          test.use
        )]] <- data.frame()
      }
    }
  }
  time_end <- Sys.time()

  log_message("DEtest done", verbose = verbose)
  log_message(
    "Elapsed time:",
    format(
      round(difftime(time_end, time_start), 2),
      format = "%Y-%m-%d %H:%M:%S"
    ),
    verbose = verbose
  )

  return(srt)
}
