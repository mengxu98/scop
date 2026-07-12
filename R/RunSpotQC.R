#' @title Run spot-level quality control
#'
#' @description
#' Calculate common spot-level QC metrics for spatial transcriptomics data and
#' label failed spots without running single-cell-specific checks such as
#' doublet calling or ambient RNA decontamination.
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param return_filtered Logical indicating whether to return a spot-filtered
#' Seurat object. Default is `FALSE`.
#' @param qc_metrics QC metrics to apply. Available metrics are `"outlier"`,
#' `"umi"`, `"gene"`, and `"mito"`.
#' @param outlier_threshold Character vector specifying outlier thresholds as
#' `"metric:direction:nmads"`. Available default metrics are
#' `"log10_nCount"`, `"log10_nFeature"`, and `"spot_featurecount_dist"`.
#' @param outlier_n Minimum number of outlier metrics required to fail a spot.
#' @param UMI_threshold Minimum UMI count required to pass `"umi"` QC.
#' @param gene_threshold Minimum detected gene count required to pass `"gene"`
#' QC.
#' @param mito_threshold Maximum mitochondrial percentage allowed by `"mito"`
#' QC.
#' @param mito_pattern Regex patterns used to identify mitochondrial genes.
#' @param mito_gene Optional explicit mitochondrial gene vector. When provided,
#' `mito_pattern` is ignored.
#' @param seed Random seed for reproducibility.
#'
#' @return A `Seurat` object with spot QC metadata columns.
#' @concept spatial-producer
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- RunSpotQC(
#'   visium_human_pancreas_sub,
#'   assay = "Spatial"
#' )
#' SpatialSpotPlot(spatial, group.by = "SpotQC")
RunSpotQC <- function(
  srt,
  assay = NULL,
  return_filtered = FALSE,
  qc_metrics = c("outlier", "umi", "gene", "mito"),
  outlier_threshold = c(
    "log10_nCount:lower:3",
    "log10_nFeature:lower:3",
    "spot_featurecount_dist:lower:3"
  ),
  outlier_n = 1,
  UMI_threshold = 500,
  gene_threshold = 200,
  mito_threshold = 20,
  mito_pattern = c("MT-", "Mt-", "mt-"),
  mito_gene = NULL,
  verbose = TRUE,
  seed = 11
) {
  log_message(
    "Running spot-level quality control",
    message_type = "running",
    verbose = verbose
  )
  set.seed(seed)

  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!assay %in% SeuratObject::Assays(srt)) {
    log_message(
      "{.arg assay} {.val {assay}} is not present in {.cls Seurat}",
      message_type = "error"
    )
  }
  qc_metrics_available <- c("outlier", "umi", "gene", "mito")
  if (any(!qc_metrics %in% qc_metrics_available)) {
    log_message(
      "{.arg qc_metrics} must be one of {.val {qc_metrics_available}}",
      message_type = "error"
    )
  }

  counts <- GetAssayData5(srt, assay = assay, layer = "counts")
  nCount <- Matrix::colSums(counts)
  nFeature <- Matrix::colSums(counts > 0)
  srt[[paste0("nCount_", assay)]] <- nCount
  srt[[paste0("nFeature_", assay)]] <- nFeature

  mito_features <- spot_qc_mito_features(
    features = rownames(counts),
    mito_pattern = mito_pattern,
    mito_gene = mito_gene
  )
  percent_mito <- rep(0, ncol(counts))
  names(percent_mito) <- colnames(counts)
  if (length(mito_features) > 0L) {
    mito_counts <- Matrix::colSums(counts[mito_features, , drop = FALSE])
    percent_mito <- ifelse(nCount > 0, mito_counts / nCount * 100, 0)
  }
  srt[["percent.mito"]] <- percent_mito

  log10_nCount <- log10(nCount)
  log10_nFeature <- log10(nFeature)
  log10_nCount[is.infinite(log10_nCount)] <- NA_real_
  log10_nFeature[is.infinite(log10_nFeature)] <- NA_real_
  srt[[paste0("log10_nCount_", assay)]] <- log10_nCount
  srt[[paste0("log10_nFeature_", assay)]] <- log10_nFeature

  spot_featurecount_dist <- spot_qc_featurecount_dist(
    log10_nCount = log10_nCount,
    log10_nFeature = log10_nFeature
  )
  srt[["spot_featurecount_dist"]] <- spot_featurecount_dist

  spot_umi_qc <- spot_gene_qc <- spot_mito_qc <- spot_outlier_qc <- character()
  if ("umi" %in% qc_metrics) {
    spot_umi_qc <- colnames(srt)[which(nCount < UMI_threshold)]
  }
  if ("gene" %in% qc_metrics) {
    spot_gene_qc <- colnames(srt)[which(nFeature < gene_threshold)]
  }
  if ("mito" %in% qc_metrics) {
    spot_mito_qc <- colnames(srt)[which(percent_mito > mito_threshold)]
  }
  if ("outlier" %in% qc_metrics) {
    outlier <- lapply(
      strsplit(outlier_threshold, ":", fixed = TRUE),
      function(rule) {
        if (length(rule) != 3L) {
          log_message(
            "{.arg outlier_threshold} entries must use {.val metric:direction:nmads}",
            message_type = "error"
          )
        }
        metric <- spot_qc_metric(
          metric = rule[[1L]],
          srt = srt,
          assay = assay,
          log10_nCount = log10_nCount,
          log10_nFeature = log10_nFeature,
          spot_featurecount_dist = spot_featurecount_dist
        )
        colnames(srt)[spot_qc_is_outlier(
          metric,
          nmads = as.numeric(rule[[3L]]),
          type = rule[[2L]]
        )]
      }
    )
    names(outlier) <- outlier_threshold
    outlier_tb <- table(unlist(outlier))
    spot_outlier_qc <- names(outlier_tb)[outlier_tb >= outlier_n]
    for (nm in names(outlier)) {
      srt[[make.names(paste0("spot_", nm))]] <- colnames(srt) %in% outlier[[nm]]
    }
  }

  SpotQC <- unique(c(
    spot_umi_qc,
    spot_gene_qc,
    spot_mito_qc,
    spot_outlier_qc
  ))
  qc_map <- list(
    spot_umi_qc = spot_umi_qc,
    spot_gene_qc = spot_gene_qc,
    spot_mito_qc = spot_mito_qc,
    spot_outlier_qc = spot_outlier_qc,
    SpotQC = SpotQC
  )
  for (qc in names(qc_map)) {
    srt[[qc]] <- ifelse(colnames(srt) %in% qc_map[[qc]], "Fail", "Pass")
    srt[[qc]] <- factor(srt[[qc, drop = TRUE]], levels = c("Pass", "Fail"))
  }

  log_message(
    "{.val {ncol(srt) - length(SpotQC)}} spots passed QC and {.val {length(SpotQC)}} spots failed QC",
    message_type = "success",
    verbose = verbose
  )

  if (isTRUE(return_filtered)) {
    srt <- srt[, srt$SpotQC == "Pass"]
  }
  srt
}

spot_qc_mito_features <- function(features, mito_pattern, mito_gene = NULL) {
  if (!is.null(mito_gene)) {
    return(intersect(mito_gene, features))
  }
  pattern <- paste0("^(", paste(mito_pattern, collapse = "|"), ")")
  grep(pattern = pattern, x = features, value = TRUE)
}

spot_qc_featurecount_dist <- function(log10_nCount, log10_nFeature) {
  keep <- is.finite(log10_nCount) & is.finite(log10_nFeature)
  out <- rep(NA_real_, length(log10_nCount))
  names(out) <- names(log10_nCount)
  if (sum(keep) < 4L || length(unique(log10_nCount[keep])) < 2L) {
    return(out)
  }
  dat <- data.frame(
    x = log10_nCount[keep],
    y = log10_nFeature[keep]
  )
  pred <- tryCatch(
    stats::predict(
      stats::loess(y ~ x, data = dat),
      newdata = data.frame(x = dat$x)
    ),
    error = function(e) rep(NA_real_, sum(keep))
  )
  out[keep] <- log10_nFeature[keep] - pred
  out
}

spot_qc_metric <- function(
  metric,
  srt,
  assay,
  log10_nCount,
  log10_nFeature,
  spot_featurecount_dist
) {
  if (identical(metric, "log10_nCount")) {
    return(log10_nCount)
  }
  if (identical(metric, "log10_nFeature")) {
    return(log10_nFeature)
  }
  if (identical(metric, "spot_featurecount_dist")) {
    return(spot_featurecount_dist)
  }
  assay_metric <- paste0(metric, "_", assay)
  if (assay_metric %in% colnames(srt@meta.data)) {
    return(srt[[assay_metric, drop = TRUE]])
  }
  if (metric %in% colnames(srt@meta.data)) {
    return(srt[[metric, drop = TRUE]])
  }
  log_message(
    "{.arg outlier_threshold} metric {.val {metric}} was not found",
    message_type = "error"
  )
}

spot_qc_is_outlier <- function(x, nmads = 3, type = c("lower", "higher")) {
  type <- match.arg(type)
  if (!is.numeric(nmads) || length(nmads) != 1L || is.na(nmads) || nmads < 0) {
    log_message(
      "{.arg nmads} must be a non-negative number",
      message_type = "error"
    )
  }
  keep <- is.finite(x)
  out <- rep(FALSE, length(x))
  if (sum(keep) == 0L) {
    return(out)
  }
  med <- stats::median(x[keep], na.rm = TRUE)
  mad_value <- stats::mad(x[keep], na.rm = TRUE)
  if (!is.finite(mad_value) || mad_value == 0) {
    return(out)
  }
  if (identical(type, "lower")) {
    out[keep] <- x[keep] < med - nmads * mad_value
  } else {
    out[keep] <- x[keep] > med + nmads * mad_value
  }
  out
}
