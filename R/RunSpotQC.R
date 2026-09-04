#' @title Run spot-level quality control
#'
#' @description
#' Calculate common spot-level QC metrics for spatial transcriptomics data and
#' label failed spots without running single-cell-specific checks such as
#' doublet calling or ambient RNA decontamination.
#'
#' @md
#' @inheritParams RunStandardWorkflow
#' @inheritParams thisutils::log_message
#' @param srt A `Seurat` object. The same object may be supplied as
#' `object =` for consistency with spatial plotting APIs.
#' @param object Optional alias for `srt`. Supply exactly one of `srt` or
#' `object`.
#' @param return_filtered Whether to return a spot-filtered
#' Seurat object.
#' @param qc_metrics QC metrics to apply. Available metrics are `"outlier"`,
#' `"umi"`, `"gene"`, and `"mito"`.
#' @param outlier_threshold Character vector specifying outlier thresholds as
#' `"metric:direction:nmads"`. Available default metrics are
#' `"log10_nCount"`, `"log10_nFeature"`, and `"spot_featurecount_dist"`.
#' @param outlier_n Positive integer giving the minimum number of outlier rules
#' required to fail a spot. It cannot exceed the number of
#' `outlier_threshold` rules.
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
#' @return A `Seurat` object with spot QC metadata columns. Cells that are not
#' present in the selected assay receive `NA` QC labels. A selected-assay cell
#' also receives `NA` when at least one requested QC is unresolved and no
#' requested QC has a known failure. Only cells labelled `Pass` are retained when
#' `return_filtered = TRUE`.
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- RunSpotQC(
#'   object = visium_human_pancreas_sub,
#'   assay = "Spatial"
#' )
#' SpatialSpotPlot(object = spatial, group.by = "SpotQC")
RunSpotQC <- function(
  srt = NULL,
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
  seed = 11,
  object = NULL
) {
  srt <- spatial_resolve_srt(srt = srt, object = object)
  log_message(
    "Running spot-level quality control",
    message_type = "running",
    verbose = verbose
  )
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
  validate_scalar_flag(return_filtered, "return_filtered")
  if ("umi" %in% qc_metrics) {
    UMI_threshold <- spot_qc_validate_cutoff(
      UMI_threshold,
      arg = "UMI_threshold",
      minimum = 0
    )
  }
  if ("gene" %in% qc_metrics) {
    gene_threshold <- spot_qc_validate_cutoff(
      gene_threshold,
      arg = "gene_threshold",
      minimum = 0
    )
  }
  if ("mito" %in% qc_metrics) {
    mito_threshold <- spot_qc_validate_cutoff(
      mito_threshold,
      arg = "mito_threshold"
    )
  }
  if ("outlier" %in% qc_metrics) {
    valid_outlier_threshold <- is.character(outlier_threshold) &&
      length(outlier_threshold) > 0L &&
      !anyNA(outlier_threshold) &&
      all(nzchar(outlier_threshold))
    if (!isTRUE(valid_outlier_threshold)) {
      log_message(
        paste0(
          "{.arg outlier_threshold} must be a non-empty character vector ",
          "of {.val metric:direction:nmads} rules"
        ),
        message_type = "error"
      )
    }
    outlier_n_valid <- is.numeric(outlier_n) &&
      length(outlier_n) == 1L &&
      !is.na(outlier_n) &&
      is.finite(outlier_n) &&
      outlier_n >= 1 &&
      outlier_n == floor(outlier_n)
    if (!isTRUE(outlier_n_valid)) {
      log_message(
        "{.arg outlier_n} must be a finite positive integer",
        message_type = "error"
      )
    }
    n_outlier_rules <- length(outlier_threshold)
    if (outlier_n > n_outlier_rules) {
      log_message(
        paste0(
          "{.arg outlier_n} cannot exceed the number of ",
          "{.arg outlier_threshold} rules ({.val {n_outlier_rules}})"
        ),
        message_type = "error"
      )
    }
    outlier_n <- as.integer(outlier_n)
  }
  set.seed(seed)

  counts <- GetAssayData5(srt, assay = assay, layer = "counts")
  evaluated_spots <- colnames(counts)
  n_evaluated <- ncol(counts)
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
    spot_umi_qc <- evaluated_spots[which(nCount < UMI_threshold)]
  }
  if ("gene" %in% qc_metrics) {
    spot_gene_qc <- evaluated_spots[which(nFeature < gene_threshold)]
  }
  if ("mito" %in% qc_metrics) {
    spot_mito_qc <- evaluated_spots[which(percent_mito > mito_threshold)]
  }
  if ("outlier" %in% qc_metrics) {
    outlier_flags <- lapply(
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
          evaluated_spots = evaluated_spots,
          log10_nCount = log10_nCount,
          log10_nFeature = log10_nFeature,
          spot_featurecount_dist = spot_featurecount_dist
        )
        nmads <- suppressWarnings(as.numeric(rule[[3L]]))
        spot_qc_is_outlier(
          metric,
          nmads = nmads,
          type = rule[[2L]]
        )
      }
    )
    names(outlier_flags) <- outlier_threshold
    outlier_matrix <- do.call(cbind, unname(outlier_flags))
    rownames(outlier_matrix) <- evaluated_spots
    colnames(outlier_matrix) <- outlier_threshold
    known_fail <- rowSums(outlier_matrix == TRUE, na.rm = TRUE)
    missing_rules <- rowSums(is.na(outlier_matrix))
    spot_outlier_status <- rep(NA_character_, n_evaluated)
    names(spot_outlier_status) <- evaluated_spots
    spot_outlier_status[known_fail >= outlier_n] <- "Fail"
    spot_outlier_status[
      is.na(spot_outlier_status) &
        known_fail + missing_rules < outlier_n
    ] <- "Pass"
    spot_outlier_qc <- evaluated_spots[which(spot_outlier_status == "Fail")]
    for (i in seq_along(outlier_flags)) {
      nm <- names(outlier_flags)[[i]]
      outlier_flag <- rep(NA, ncol(srt))
      names(outlier_flag) <- colnames(srt)
      outlier_flag[evaluated_spots] <- outlier_flags[[i]]
      srt[[make.names(paste0("spot_", nm))]] <- outlier_flag
    }
  } else {
    spot_outlier_status <- setNames(
      rep("Pass", n_evaluated),
      evaluated_spots
    )
  }

  pass_fail_status <- function(failed_spots) {
    out <- ifelse(evaluated_spots %in% failed_spots, "Fail", "Pass")
    names(out) <- evaluated_spots
    out
  }
  qc_status_map <- list(
    spot_umi_qc = pass_fail_status(spot_umi_qc),
    spot_gene_qc = pass_fail_status(spot_gene_qc),
    spot_mito_qc = pass_fail_status(spot_mito_qc),
    spot_outlier_qc = spot_outlier_status
  )
  requested_status_names <- paste0("spot_", qc_metrics, "_qc")
  if (length(requested_status_names) == 0L) {
    overall_status <- setNames(rep("Pass", n_evaluated), evaluated_spots)
  } else {
    requested_status <- do.call(
      cbind,
      unname(qc_status_map[requested_status_names])
    )
    known_fail <- rowSums(requested_status == "Fail", na.rm = TRUE) > 0L
    not_evaluated <- rowSums(is.na(requested_status)) > 0L
    overall_status <- ifelse(
      known_fail,
      "Fail",
      ifelse(not_evaluated, NA_character_, "Pass")
    )
    names(overall_status) <- evaluated_spots
  }
  qc_status_map$SpotQC <- overall_status
  for (qc in names(qc_status_map)) {
    qc_status <- rep(NA_character_, ncol(srt))
    names(qc_status) <- colnames(srt)
    qc_status[evaluated_spots] <- qc_status_map[[qc]]
    srt[[qc]] <- qc_status
    srt[[qc]] <- factor(srt[[qc, drop = TRUE]], levels = c("Pass", "Fail"))
  }

  selected_labels <- as.character(
    srt@meta.data[evaluated_spots, "SpotQC", drop = TRUE]
  )
  n_passed <- sum(selected_labels == "Pass", na.rm = TRUE)
  n_failed <- sum(selected_labels == "Fail", na.rm = TRUE)
  n_not_evaluated <- sum(is.na(selected_labels))
  n_evaluated <- n_passed + n_failed
  if (isTRUE(return_filtered)) {
    keep <- !is.na(srt$SpotQC) & srt$SpotQC == "Pass"
    if (!any(keep)) {
      log_message(
        paste0(
          "No spots passed QC; set {.arg return_filtered = FALSE} ",
          "to inspect the unfiltered result"
        ),
        message_type = "error"
      )
    }
    srt <- srt[, keep]
  }

  n_returned <- ncol(srt)
  done <- paste0(
    "Spot QC completed: {.val {n_evaluated}} evaluated, ",
    "{.val {n_passed}} Pass, {.val {n_failed}} Fail"
  )
  if (n_not_evaluated > 0L) {
    done <- paste0(
      done,
      "; {.val {n_not_evaluated}} NotEvaluated"
    )
  }
  if (isTRUE(return_filtered)) {
    done <- paste0(done, "; {.val {n_returned}} returned")
  }
  receipt_verbose <- thisutils::get_verbose(verbose)
  plot_call <- inspect_call <- NULL
  if (isTRUE(receipt_verbose)) {
    image_names <- tryCatch(
      SeuratObject::Images(srt),
      error = function(e) character()
    )
    can_plot <- !isTRUE(return_filtered) && spot_qc_has_plottable_coordinates(
      srt = srt,
      image_names = image_names
    )
    plot_call <- if (isTRUE(can_plot) && length(image_names) > 1L) {
      paste0(
        "lapply(SeuratObject::Images(<returned_object>), function(image) ",
        "SpatialSpotPlot(<returned_object>, group.by = \"SpotQC\", image = image))"
      )
    } else if (isTRUE(can_plot)) {
      "SpatialSpotPlot(<returned_object>, group.by = \"SpotQC\")"
    } else {
      NULL
    }
    inspect_call <- if (!isTRUE(return_filtered) && !isTRUE(can_plot)) {
      "table(<returned_object>[[\"SpotQC\", drop = TRUE]], useNA = \"ifany\")"
    } else {
      NULL
    }
  }
  spatial_run_receipt(
    done = done,
    scope = "assay {.val {assay}}, layer {.val counts}",
    saved = "metadata column {.var SpotQC}",
    plot = plot_call,
    inspect = inspect_call,
    status = if (n_not_evaluated > 0L) "partial" else "completed",
    verbose = receipt_verbose,
    .envir = environment()
  )
  srt
}

spot_qc_has_plottable_coordinates <- function(srt, image_names = NULL) {
  if (is.null(image_names)) {
    image_names <- tryCatch(
      SeuratObject::Images(srt),
      error = function(e) character()
    )
  }
  targets <- if (length(image_names) > 0L) as.list(image_names) else list(NULL)
  all(vapply(
    targets,
    function(image) {
      tryCatch(
        {
          coords <- spatial_coords_raw(
            srt = srt,
            image = image,
            coord.cols = c("col", "row"),
            image.scale = "lowres",
            require_scale = TRUE,
            image_policy = "strict"
          )$data
          nrow(coords) > 0L
        },
        error = function(e) FALSE
      )
    },
    logical(1)
  ))
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
  evaluated_spots,
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
    metric_values <- srt@meta.data[
      evaluated_spots,
      assay_metric,
      drop = TRUE
    ]
    return(spot_qc_validate_metric(
      metric_values,
      metric = assay_metric,
      n_expected = length(evaluated_spots)
    ))
  }
  if (metric %in% colnames(srt@meta.data)) {
    metric_values <- srt@meta.data[evaluated_spots, metric, drop = TRUE]
    return(spot_qc_validate_metric(
      metric_values,
      metric = metric,
      n_expected = length(evaluated_spots)
    ))
  }
  log_message(
    "{.arg outlier_threshold} metric {.val {metric}} was not found",
    message_type = "error"
  )
}

spot_qc_validate_metric <- function(x, metric, n_expected) {
  if (!is.numeric(x) || length(x) != n_expected) {
    log_message(
      paste0(
        "{.arg outlier_threshold} metric {.val {metric}} must be a numeric ",
        "vector with one value per selected-assay spot"
      ),
      message_type = "error"
    )
  }
  x
}

spot_qc_validate_cutoff <- function(x, arg, minimum = -Inf) {
  if (
    !is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      x < minimum
  ) {
    qualifier <- if (is.finite(minimum)) "non-negative " else ""
    log_message(
      paste0("{.arg {arg}} must be a finite ", qualifier, "number"),
      message_type = "error"
    )
  }
  as.numeric(x)
}

spot_qc_is_outlier <- function(x, nmads = 3, type = c("lower", "higher")) {
  type <- match.arg(type)
  if (
    !is.numeric(nmads) ||
      length(nmads) != 1L ||
      is.na(nmads) ||
      !is.finite(nmads) ||
      nmads < 0
  ) {
    log_message(
      "{.arg nmads} must be a non-negative number",
      message_type = "error"
    )
  }
  keep <- is.finite(x)
  out <- rep(NA, length(x))
  names(out) <- names(x)
  out[keep] <- FALSE
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
