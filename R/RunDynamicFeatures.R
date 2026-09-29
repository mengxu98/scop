#' @title Calculates dynamic features for lineages
#'
#' @md
#' @inheritParams thisutils::parallelize_fun
#' @inheritParams RunStandardWorkflow
#' @inheritParams GroupHeatmap
#' @param lineages Lineage names for which dynamic features should be calculated.
#' @param features Features to use.
#' If `NULL`, n_candidates must be provided.
#' @param suffix Suffix to append to the output layer names for each lineage.
#' Default is the lineage names.
#' @param n_candidates A number of candidate features to select when features is `NULL`.
#' @param minfreq Minimum frequency threshold for candidate features.
#' Features with a frequency less than minfreq will be excluded.
#' @param libsize A numeric or numeric vector specifying the library size correction factors for each cell.
#' If NULL, the library size correction factors will be calculated based on the expression matrix.
#' If length(libsize) is 1, the same value will be used for all cells.
#' Otherwise, libsize must have the same length as the number of cells in srt.
#' @param fit_method The method used for fitting features.
#' Either `"gam"` (generalized additive models) or `"pretsa"` (Pattern recognition in Temporal and Spatial Analyses).
#' @param family A character or character vector specifying the family of distributions to use for the GAM.
#' If family is set to NULL, the appropriate family will be automatically determined based on the data.
#' If length(family) is 1, the same family will be used for all features.
#' Otherwise, family must have the same length as features.
#' @param knot For `fit_method = "pretsa"`: B-spline knots. `0` or `"auto"`.
#' @param max_knot_allowed For `fit_method = "pretsa"` when `knot = "auto"`: max knots.
#' @param padjust_method The method used for p-value adjustment.
#'
#' @return
#' Returns the modified Seurat object with the calculated dynamic features stored in the tools slot.
#'
#' @seealso
#' [DynamicHeatmap], [DynamicPlot], [RunDynamicEnrichment]
#'
#' @export
#'
#' @references
#' Zhuang, H., Ji, Z. PreTSA: computationally efficient modeling of temporal and spatial gene expression patterns.
#' Genome Biol (2026). https://doi.org/10.1186/s13059-026-03994-3
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#' pancreas_sub <- RunSlingshot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP"
#' )
#'
#' pancreas_sub <- RunDynamicFeatures(
#'   pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   n_candidates = 200,
#'   fit_method = "gam"
#' )
#'
#' names(
#'   pancreas_sub@tools$DynamicFeatures_Lineage1
#' )
#' head(
#'   pancreas_sub@tools$DynamicFeatures_Lineage1$DynamicFeatures
#' )
#' ht <- DynamicHeatmap(
#'   pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   cell_annotation = "SubCellType",
#'   n_split = 3,
#'   reverse_ht = "Lineage1"
#' )
#' ht$plot
#'
#' DynamicPlot(
#'   pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   features = c("Arxes1", "Ncoa2"),
#'   group.by = "SubCellType",
#'   compare_lineages = TRUE,
#'   compare_features = FALSE
#' )
#'
#' pancreas_sub <- RunDynamicFeatures(
#'   pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   n_candidates = 200,
#'   fit_method = "pretsa"
#' )
#' head(
#'   pancreas_sub@tools$DynamicFeatures_Lineage1$DynamicFeatures
#' )
#' ht <- DynamicHeatmap(
#'   pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   cell_annotation = "SubCellType",
#'   n_split = 3,
#'   reverse_ht = "Lineage1"
#' )
#' ht$plot
#'
#' DynamicPlot(
#'   pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   features = c("Arxes1", "Ncoa2"),
#'   group.by = "SubCellType",
#'   compare_lineages = TRUE,
#'   compare_features = FALSE
#' )
RunDynamicFeatures <- function(
  object,
  lineages,
  features = NULL,
  suffix = lineages,
  n_candidates = 1000,
  minfreq = 5,
  family = NULL,
  layer = "counts",
  assay = NULL,
  libsize = NULL,
  fit_method = c("gam", "pretsa"),
  knot = 0,
  max_knot_allowed = 10,
  padjust_method = "fdr",
  cores = 1,
  verbose = TRUE,
  seed = 11,
  srt = NULL
) {
  srt <- resolve_deprecated_srt(object, srt, missing(object))
  set.seed(seed)
  assay <- assay %||% DefaultAssay(srt)
  fit_method <- match.arg(fit_method)

  log_message(
    "Start find dynamic features",
    verbose = verbose
  )

  if (fit_method == "gam") {
    check_r("mgcv", verbose = FALSE)
  }
  meta <- c()
  gene <- c()
  if (!is.null(features)) {
    gene <- features[features %in% rownames(srt[[assay]])]
    meta <- features[features %in% colnames(srt@meta.data)]
    isnum <- sapply(
      srt@meta.data[, meta, drop = FALSE], is.numeric
    )
    if (!all(isnum)) {
      log_message(
        "{.val {meta[!isnum]}} is not numeric and will be dropped",
        message_type = "warning",
        verbose = verbose
      )
      meta <- meta[isnum]
    }
    features <- c(gene, meta)
    if (length(features) == 0) {
      log_message(
        "No feature found in the srt object",
        message_type = "error"
      )
    }
  }

  y_mat <- GetAssayData5(
    srt,
    layer = layer,
    assay = assay
  )
  if (is.null(libsize)) {
    status <- CheckDataType(
      srt,
      assay = assay,
      layer = "counts",
      verbose = verbose
    )
    if (status != "raw_counts") {
      y_libsize <- stats::setNames(
        rep(1, ncol(srt)),
        colnames(srt)
      )
    } else {
      y_libsize <- Matrix::colSums(
        GetAssayData5(
          srt,
          assay = assay,
          layer = "counts"
        )
      )
    }
  } else {
    if (length(libsize) == 1) {
      y_libsize <- stats::setNames(
        rep(libsize, ncol(srt)),
        colnames(srt)
      )
    } else if (length(libsize) == ncol(srt)) {
      y_libsize <- stats::setNames(libsize, colnames(srt))
    } else {
      log_message(
        "{.arg libsize} must be length of 1 or the number of cells",
        message_type = "error"
      )
    }
  }

  if (length(meta) > 0) {
    y_mat <- rbind(
      y_mat, Matrix::t(srt@meta.data[, meta, drop = FALSE])
    )
  }

  features_list <- c()
  srt_sub_list <- list()
  for (l in lineages) {
    srt_sub <- subset(
      srt,
      cell = rownames(srt@meta.data)[is.finite(srt@meta.data[[l]])]
    )
    if (is.null(features)) {
      if (is.null(n_candidates)) {
        log_message(
          "{.arg features} or {.arg n_candidates} must be provided at least one",
          message_type = "error"
        )
      }
      HVF <- SeuratObject::VariableFeatures(
        FindVariableFeatures(
          srt_sub,
          nfeatures = n_candidates,
          assay = assay,
          verbose = FALSE
        ),
        assay = assay
      )
      HVF_counts <- GetAssayData5(
        srt_sub,
        assay = assay,
        layer = "counts"
      )[HVF, , drop = FALSE]
      HVF <- HVF[dynamic_row_unique_counts(HVF_counts) >= minfreq]
      features_list[[l]] <- HVF
    } else {
      features_list[[l]] <- features
    }
    srt_sub_list[[l]] <- srt_sub
  }
  features <- unique(unlist(features_list))
  gene <- features[features %in% rownames(srt[[assay]])]
  meta <- features[features %in% colnames(srt@meta.data)]
  log_message(
    "Number of candidate features (union): {.val {length(features)}}",
    verbose = verbose
  )

  gene_status <- CheckDataType(srt, assay = assay, layer = layer)
  meta_status <- sapply(meta, function(x) {
    CheckDataType(srt[[x]])
  })
  if (is.null(family)) {
    family <- rep("gaussian", length(features))
    names(family) <- features
    family[names(meta_status)[meta_status == "raw_counts"]] <- "nb"
    if (gene_status == "raw_counts") {
      family[gene] <- "nb"
    }
  } else {
    if (length(family) == 1) {
      family <- rep(family, length(features))
      names(family) <- features
    }
    if (length(family) != length(features)) {
      log_message(
        "{.arg family} must be one character or a vector of the same length as features",
        message_type = "error"
      )
    }
  }

  for (i in seq_along(lineages)) {
    l <- lineages[i]
    srt_sub <- srt_sub_list[[l]]
    t <- srt_sub[[l, drop = TRUE]]
    t <- t[is.finite(t)]
    t_ordered <- t[order(t)]
    y_ordered <- as_matrix(
      y_mat[features, names(t_ordered), drop = FALSE]
    )
    l_libsize <- y_libsize[names(t_ordered)]
    raw_matrix <- dynamic_raw_matrix(y_ordered, t_ordered)

    log_message(
      "Calculating dynamic features for {.val {l}}...",
      verbose = verbose
    )
    if (fit_method == "gam") {
      reference <- stats::median(y_libsize[is.finite(y_libsize) & y_libsize > 0])
      fitted <- thisutils::fit_trends(
        y_ordered, t_ordered,
        method = "gam", family = family,
        exposure = l_libsize, reference_exposure = reference,
        use_exposure = layer == "counts" & !features %in% meta,
        padjust_method = padjust_method, cores = cores, verbose = verbose
      )
    } else {
      if (gene_status == "raw_counts" || layer == "counts") {
        y_ordered[gene, ] <- log1p(y_ordered[gene, , drop = FALSE])
      }
      fitted <- thisutils::fit_trends(
        y_ordered, t_ordered,
        method = "pretsa", knot = knot,
        max_knot_allowed = max_knot_allowed,
        padjust_method = padjust_method, verbose = verbose
      )
      for (group in list(gene, meta)) {
        if (length(group)) {
          fitted$statistics[group, "padjust"] <- stats::p.adjust(
            fitted$statistics[group, "pvalue"],
            method = padjust_method
          )
        }
      }
    }
    DF <- fitted$statistics
    names(DF)[names(DF) == "n_above_min"] <- "exp_ncells"
    out <- list(
      DynamicFeatures = DF,
      fitted_matrix = cbind(pseudotime = t_ordered, t(fitted$fitted)),
      upr_matrix = cbind(pseudotime = t_ordered, t(fitted$upper)),
      lwr_matrix = cbind(pseudotime = t_ordered, t(fitted$lower))
    )
    res <- list(
      DynamicFeatures = out$DynamicFeatures,
      raw_matrix = raw_matrix,
      fitted_matrix = out$fitted_matrix,
      upr_matrix = out$upr_matrix,
      lwr_matrix = out$lwr_matrix,
      libsize = l_libsize,
      lineages = l,
      family = family
    )
    srt@tools[[paste0("DynamicFeatures_", suffix[i])]] <- res
  }

  log_message(
    "Find dynamic features done",
    message_type = "success",
    verbose = verbose
  )

  return(srt)
}

dynamic_row_unique_counts <- function(x) {
  out <- if (inherits(x, "sparseMatrix")) {
    dynamic_row_unique_counts_sparse_cpp(methods::as(x, "dgCMatrix"))
  } else {
    dynamic_row_unique_counts_dense_cpp(as.matrix(x))
  }
  names(out) <- rownames(x)
  out
}

dynamic_raw_matrix <- function(y_ordered, t_ordered) {
  raw_matrix <- matrix(
    NA_real_,
    nrow = length(t_ordered),
    ncol = nrow(y_ordered) + 1L,
    dimnames = list(
      names(t_ordered),
      c("pseudotime", rownames(y_ordered))
    )
  )
  raw_matrix[, "pseudotime"] <- t_ordered
  if (nrow(y_ordered) > 0) {
    raw_matrix[, rownames(y_ordered)] <- t(y_ordered)
  }
  raw_matrix
}
