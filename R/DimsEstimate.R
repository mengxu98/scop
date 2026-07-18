#' @title Estimate useful dimensions from a reduction
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt A `Seurat` object.
#' @param reduction Name of the dimensional reduction to inspect.
#' Default is `NULL`, which automatically selects a PCA-like reduction via
#' [DefaultReduction()] with `pattern = "pca"`.
#' @param reduction_method Optional reduction method name.
#' When set to `"nmf"` or `"glmpca"`, all available dimensions will be retained.
#' @param k Number of neighbors used by [intrinsicDimension::maxLikGlobalDimEst].
#' Default is `30`.
#' @param method Dimension-selection method. `"scree"` uses PCA standard
#' deviations with broken-stick, elbow, cumulative-variance, and marginal-gain
#' criteria. `"intrinsic"` uses [intrinsicDimension::maxLikGlobalDimEst].
#' `"ensemble"` keeps the larger recommendation from both methods when both are
#' available. Default is `"scree"`.
#' @param min_dims Minimum number of dimensions kept when intrinsic-dimension estimation succeeds.
#' Default is `5`.
#' @param variance_threshold Cumulative variance threshold used by
#' `method = "scree"`. Default is `0.8`.
#' @param marginal_gain_threshold Stop point for marginal variance gain
#' (percentage points) used by `method = "scree"`. Default is `0.5`.
#' @param skip_first Whether to drop the first dimension from the returned result.
#' Useful for `TFIDF/LSI` workflows. Default is `FALSE`.
#' @param use_stored Whether to use `misc$dims_estimate` already stored in the reduction when available.
#' Default is `TRUE`.
#'
#' @return An integer vector of dimensions to use.
#'
#' @export
#' @seealso [DimsEstimatePlot]
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' RunDimsEstimate(pancreas_sub)
#'
#' DimsEstimatePlot(pancreas_sub)
RunDimsEstimate <- function(
  srt,
  reduction = NULL,
  reduction_method = NULL,
  k = 30L,
  method = c("scree", "intrinsic", "ensemble"),
  min_dims = 5L,
  variance_threshold = 0.8,
  marginal_gain_threshold = 0.5,
  skip_first = FALSE,
  use_stored = TRUE,
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  method <- match.arg(method)
  if (is.null(reduction)) {
    reduction <- DefaultReduction(
      srt = srt,
      pattern = "pca",
      min_dim = 2L
    )
  } else {
    if (!is.character(reduction) || length(reduction) != 1L) {
      log_message(
        "{.arg reduction} must be a single reduction name",
        message_type = "error"
      )
    }
    reduction <- DefaultReduction(
      srt = srt,
      pattern = reduction,
      min_dim = 2L
    )
  }

  reduction_obj <- srt@reductions[[reduction]]
  reduction_ncol <- ncol(reduction_obj@cell.embeddings)
  if (is.null(reduction_ncol) || reduction_ncol < 1L) {
    log_message(
      "{.val {reduction}} does not contain valid embeddings",
      message_type = "error"
    )
  }

  reduction_method <- tolower(reduction_method %||% reduction)
  dims_source <- NULL
  dims_use <- integer()
  if (isTRUE(use_stored)) {
    dims_use_stored <- dims_estimate_validate(
      dims_use = reduction_obj@misc[["dims_estimate"]],
      reduction_ncol = reduction_ncol,
      skip_first = skip_first
    )
    if (length(dims_use_stored) > 0L) {
      dims_use <- dims_estimate_upgrade(
        dims_use = dims_use_stored,
        reduction_obj = reduction_obj,
        reduction_ncol = reduction_ncol,
        reduction_method = reduction_method,
        min_dims = min_dims,
        skip_first = skip_first
      )
      dims_source <- if (
        length(dims_use) == length(dims_use_stored) &&
          identical(dims_use, dims_use_stored)
      ) {
        "stored"
      } else {
        "stored_adjusted"
      }
    }
  }

  if (is.null(dims_source)) {
    if (reduction_method %in% c("glmpca", "nmf")) {
      dims_use <- dims_estimate_validate(
        dims_use = seq_len(reduction_ncol),
        reduction_ncol = reduction_ncol,
        skip_first = skip_first
      )
      dims_source <- "all"
    } else {
      dim_est <- NA_real_
      scree_est <- NA_integer_
      if (method %in% c("scree", "ensemble")) {
        scree_est <- dims_estimate_scree_recommendation(
          stdev = reduction_obj@stdev,
          max_pcs = min(50L, reduction_ncol),
          variance_threshold = variance_threshold,
          marginal_gain_threshold = marginal_gain_threshold,
          min_dims = min_dims,
          skip_first = skip_first
        )
      }
      if (method %in% c("intrinsic", "ensemble")) {
        check_r("intrinsicDimension", verbose = FALSE)
        dim_est <- dims_estimate_intrinsic(
          srt = srt,
          reduction = reduction,
          k = k,
          reduction_ncol = reduction_ncol,
          verbose = verbose
        )
      }

      if (method == "scree" && !is.na(scree_est)) {
        dims_use <- dims_estimate_validate(
          dims_use = seq_len(scree_est),
          reduction_ncol = reduction_ncol,
          skip_first = skip_first
        )
        dims_source <- "scree"
      } else if (method == "intrinsic" && !is.na(dim_est)) {
        dims_use <- dims_estimate_validate(
          dims_use = seq_len(max(
            min(reduction_ncol, as.integer(min_dims)),
            ceiling(dim_est)
          )),
          reduction_ncol = reduction_ncol,
          skip_first = skip_first
        )
        dims_source <- "intrinsic"
      } else if (method == "ensemble") {
        estimate_candidates <- c(scree_est, ceiling(dim_est))
        estimate_candidates <- estimate_candidates[!is.na(estimate_candidates)]
        if (length(estimate_candidates) > 0L) {
          dims_use <- dims_estimate_validate(
            dims_use = seq_len(max(
              min(reduction_ncol, as.integer(min_dims)),
              estimate_candidates
            )),
            reduction_ncol = reduction_ncol,
            skip_first = skip_first
          )
          dims_source <- "ensemble"
        }
      }
    }
  }

  if (length(dims_use) == 0L) {
    log_message(
      "No valid dimensions can be estimated for {.pkg {reduction}}",
      message_type = "error"
    )
  }

  if (isTRUE(verbose)) {
    if (dims_source == "stored") {
      log_message(
        "Use stored estimated dimensions {.val {min(dims_use)}}:{.val {max(dims_use)}} for {.pkg {reduction}}",
        verbose = verbose
      )
    } else if (dims_source == "stored_adjusted") {
      log_message(
        "Use adjusted stored dimensions {.val {min(dims_use)}}:{.val {max(dims_use)}} for {.pkg {reduction}}",
        verbose = verbose
      )
    } else if (dims_source %in% c("scree", "intrinsic", "ensemble")) {
      log_message(
        "Estimated dimensions {.val {min(dims_use)}}:{.val {max(dims_use)}} for {.pkg {reduction}} with {.val {dims_source}} method",
        verbose = verbose
      )
    } else if (dims_source == "all") {
      log_message(
        "Use all available dimensions {.val {min(dims_use)}}:{.val {max(dims_use)}} for {.pkg {reduction}}",
        verbose = verbose
      )
    } else {
      log_message(
        "No valid estimated dimensions found for {.pkg {reduction}}. Use fallback dimensions {.val {min(dims_use)}}:{.val {max(dims_use)}}",
        message_type = "warning",
        verbose = verbose
      )
    }
  }

  dims_use
}

dims_estimate_validate <- function(
  dims_use,
  reduction_ncol,
  skip_first = FALSE
) {
  dims_use <- sort(unique(as.integer(dims_use)))
  dims_use <- dims_use[!is.na(dims_use)]
  lower_bound <- if (isTRUE(skip_first)) 2L else 1L
  dims_use[dims_use >= lower_bound & dims_use <= reduction_ncol]
}

dims_estimate_upgrade <- function(
  dims_use,
  reduction_obj,
  reduction_ncol,
  reduction_method = "",
  min_dims = 10L,
  skip_first = FALSE
) {
  lower_bound <- if (isTRUE(skip_first)) 2L else 1L
  if (reduction_method %in% c("glmpca", "nmf")) {
    return(seq.int(lower_bound, reduction_ncol))
  }

  target_max_dim <- max(
    max(dims_use),
    min(reduction_ncol, as.integer(min_dims)),
    dims_estimate_scree_lower_bound(
      stdev = reduction_obj@stdev,
      max_pcs = min(50L, reduction_ncol),
      skip_first = skip_first
    )
  )
  dims_estimate_validate(
    dims_use = seq_len(target_max_dim),
    reduction_ncol = reduction_ncol,
    skip_first = skip_first
  )
}

dims_estimate_scree_lower_bound <- function(
  stdev,
  max_pcs = 50L,
  skip_first = FALSE
) {
  lower_bound <- if (isTRUE(skip_first)) 2L else 1L
  if (is.null(stdev) || length(stdev) == 0L) {
    return(lower_bound)
  }

  stdev <- as.numeric(stdev)
  stdev <- stdev[is.finite(stdev)]
  if (length(stdev) == 0L) {
    return(lower_bound)
  }

  max_pcs <- min(as.integer(max_pcs), length(stdev))
  if (is.na(max_pcs) || max_pcs < lower_bound) {
    return(lower_bound)
  }

  stdev_use <- stdev[seq_len(max_pcs)]
  total_var <- sum(stdev_use^2)
  if (!is.finite(total_var) || total_var <= 0) {
    return(lower_bound)
  }

  pct_var <- stdev_use^2 / total_var * 100
  broken_stick <- rev(cumsum(1 / seq_len(max_pcs))) / max_pcs * 100
  broken_stick_hits <- which(pct_var >= broken_stick)
  broken_stick_point <- if (length(broken_stick_hits) == 0L) {
    lower_bound
  } else {
    max(broken_stick_hits)
  }

  elbow_point <- lower_bound
  if (length(stdev_use) >= 3L) {
    curvature <- diff(diff(stdev_use))
    curvature_idx <- which.max(abs(curvature))[1]
    if (!is.na(curvature_idx)) {
      elbow_point <- max(lower_bound, curvature_idx + 2L)
    }
  }

  max(lower_bound, broken_stick_point, elbow_point)
}

dims_estimate_scree_recommendation <- function(
  stdev,
  max_pcs = 50L,
  variance_threshold = 0.8,
  marginal_gain_threshold = 0.5,
  min_dims = 5L,
  skip_first = FALSE
) {
  lower_bound <- if (isTRUE(skip_first)) 2L else 1L
  if (is.null(stdev) || length(stdev) == 0L) {
    return(NA_integer_)
  }

  stdev <- as.numeric(stdev)
  stdev <- stdev[is.finite(stdev)]
  if (length(stdev) == 0L) {
    return(NA_integer_)
  }

  max_pcs <- min(as.integer(max_pcs), length(stdev))
  if (is.na(max_pcs) || max_pcs < lower_bound) {
    return(NA_integer_)
  }

  stdev_use <- stdev[seq_len(max_pcs)]
  variance <- stdev_use^2
  total_var <- sum(stdev^2)
  if (!is.finite(total_var) || total_var <= 0) {
    return(NA_integer_)
  }

  cumulative_var <- cumsum(variance) / total_var
  marginal_gain <- variance / total_var * 100
  variance_point <- which(cumulative_var >= variance_threshold)[1]
  if (is.na(variance_point)) {
    variance_point <- max_pcs
  }
  marginal_point <- which(marginal_gain < marginal_gain_threshold)[1]
  if (is.na(marginal_point)) {
    marginal_point <- max_pcs
  }

  max(
    lower_bound,
    min(max_pcs, as.integer(min_dims)),
    dims_estimate_scree_lower_bound(
      stdev = stdev,
      max_pcs = max_pcs,
      skip_first = skip_first
    ),
    min(variance_point, marginal_point)
  )
}

dims_estimate_intrinsic <- function(
  srt,
  reduction,
  k = 30L,
  reduction_ncol,
  verbose = TRUE
) {
  tryCatch(
    expr = {
      min(
        intrinsicDimension::maxLikGlobalDimEst(
          data = Seurat::Embeddings(srt, reduction = reduction),
          k = as.integer(k)
        )[["dim.est"]],
        reduction_ncol
      )
    },
    error = function(e) {
      log_message(
        "Can not estimate intrinsic dimensions with {.pkg maxLikGlobalDimEst}",
        message_type = "warning",
        verbose = verbose
      )
      NA_real_
    }
  )
}

pc_selection_stats <- function(
  srt,
  reduction = NULL,
  max_pcs = 50,
  variance_thresholds = c(0.60, 0.70, 0.80, 0.90)
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  if (is.null(reduction)) {
    reduction <- DefaultReduction(
      srt = srt,
      pattern = "pca",
      min_dim = 2L
    )
  } else {
    if (!is.character(reduction) || length(reduction) != 1L) {
      log_message(
        "{.arg reduction} must be a single reduction name",
        message_type = "error"
      )
    }
    reduction <- DefaultReduction(
      srt = srt,
      pattern = reduction,
      min_dim = 2L
    )
  }

  reduction_obj <- srt@reductions[[reduction]]
  stdev_all <- reduction_obj@stdev
  if (is.null(stdev_all) || length(stdev_all) == 0L) {
    log_message(
      "{.pkg {reduction}} does not contain standard deviations. This visualization currently supports PCA-like reductions",
      message_type = "error"
    )
  }

  max_pcs <- min(as.integer(max_pcs), length(stdev_all))
  stdev <- stdev_all[seq_len(max_pcs)]
  variance <- stdev^2
  total_var <- sum(stdev_all^2)
  pct_var <- variance / total_var * 100
  cumulative_var <- cumsum(variance) / total_var * 100
  marginal_gain <- c(pct_var[1], diff(cumulative_var))
  diff1 <- diff(stdev)
  diff2 <- diff(diff1)
  curvature <- c(0, 0, diff2)
  if (length(curvature) < max_pcs) {
    curvature <- c(curvature, rep(NA_real_, max_pcs - length(curvature)))
  }

  cumulative_pct_full <- cumsum(stdev_all^2) / sum(stdev_all^2)
  pc_at_thresholds <- sapply(variance_thresholds, function(thresh) {
    which(cumulative_pct_full > thresh)[1]
  })

  annotation_data <- data.frame(
    pc = pc_at_thresholds,
    threshold = variance_thresholds * 100,
    label = paste0(pc_at_thresholds, " PCs\n(", variance_thresholds * 100, "%)")
  )

  plot_data <- data.frame(
    PC = seq_len(max_pcs),
    stdev = stdev,
    individual_var = pct_var,
    cumulative_var = cumulative_var,
    marginal_gain = marginal_gain,
    curvature = curvature[seq_len(max_pcs)]
  )

  list(
    plot_data = plot_data,
    annotation_data = annotation_data,
    threshold_point = which(marginal_gain < 0.5)[1],
    elbow_point = if (max_pcs >= 3L) {
      which.max(abs(curvature[3:min(max_pcs, 30L)])) + 2L
    } else {
      NA_integer_
    }
  )
}

#' @title Dimension estimate diagnostic plot
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt A `Seurat` object with a PCA-like reduction computed.
#' @param max_pcs Maximum number of PCs to visualize. Default is `50`.
#' @param variance_thresholds Numeric vector of variance thresholds to mark.
#' Default is `c(0.60, 0.70, 0.80, 0.90)`.
#' @param reduction Reduction name to inspect. Default is `NULL`, which
#' automatically selects a PCA-like reduction via [DefaultReduction()] with
#' `pattern = "pca"`.
#' @param palcolor Colors for the selected-PC line, curves, and bars,
#' respectively. Default is `c("#D70440", "#0AA344", "#1772B4")`.
#' @param aspect.ratio Aspect ratio of the plot. Default is `NULL`.
#' @param title Plot title. When `NULL` (default),
#' reports the selected number of PCs.
#' @param subtitle Plot subtitle. Default is
#' `NULL`.
#' @param xlab X-axis label. Default is `"Principal component"`.
#' @param theme_use Theme function used to style the plot.
#' Default is `"theme_scop"`.
#' @param theme_args Other arguments passed to the `theme_use`.
#' @param seed Random seed. Default is `11`.
#'
#' @return A `ggplot` object showing per-PC explained variance (bars, left
#' axis) and cumulative explained variance (line, right axis).
#'
#' @export
#'
#' @seealso [RunDimsEstimate]
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' DimsEstimatePlot(pancreas_sub)
DimsEstimatePlot <- function(
  srt,
  max_pcs = 50,
  variance_thresholds = c(0.60, 0.70, 0.80, 0.90),
  reduction = NULL,
  palcolor = c("#D70440", "#0AA344", "#1772B4"),
  aspect.ratio = NULL,
  title = NULL,
  subtitle = NULL,
  xlab = "Principal component",
  theme_use = "theme_scop",
  theme_args = list(),
  seed = 11,
  verbose = TRUE
) {
  set.seed(seed)
  if (is.null(reduction)) {
    reduction <- DefaultReduction(
      srt = srt,
      pattern = "pca",
      min_dim = 2L
    )
  } else {
    if (!is.character(reduction) || length(reduction) != 1L) {
      log_message(
        "{.arg reduction} must be a single reduction name",
        message_type = "error"
      )
    }
    reduction <- DefaultReduction(
      srt = srt,
      pattern = reduction,
      min_dim = 2L
    )
  }
  stats_use <- pc_selection_stats(
    srt = srt,
    reduction = reduction,
    max_pcs = max_pcs,
    variance_thresholds = variance_thresholds
  )
  plot_data <- stats_use[["plot_data"]]

  recommended_dims <- RunDimsEstimate(
    srt = srt,
    reduction = reduction,
    reduction_method = reduction,
    use_stored = TRUE,
    verbose = FALSE
  )
  recommended_pcs <- max(recommended_dims)

  if (!is.character(palcolor) || length(palcolor) < 3L) {
    log_message(
      "{.arg palcolor} must provide three colors: selected-PC line, curves, and bars",
      message_type = "error"
    )
  }
  selection_color <- palcolor[1]
  cumulative_color <- palcolor[2]
  variance_color <- palcolor[3]
  if (identical(theme_use, "theme_scop")) {
    theme_use <- "theme_this"
  }

  apply_theme <- function(p) {
    p +
      do.call(theme_use, theme_args) +
      ggplot2::theme(
        aspect.ratio = aspect.ratio,
        legend.position = "none"
      )
  }

  max_individual_var <- max(plot_data$individual_var, na.rm = TRUE)
  max_cumulative_var <- max(plot_data$cumulative_var, na.rm = TRUE)
  scale_factor <- max_individual_var / max_cumulative_var
  y_upper <- max_individual_var * 1.08
  plot_data$variance_series <- "Variance explained"
  plot_data$cumulative_series <- "Cumulative"
  threshold_data <- data.frame(
    scaled_threshold = variance_thresholds * 100 * scale_factor
  )

  p1 <- apply_theme(
    ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$PC)) +
      ggplot2::geom_col(
        ggplot2::aes(y = .data$individual_var, fill = .data$variance_series),
        alpha = 0.75,
        width = 0.82
      ) +
      ggplot2::geom_line(
        ggplot2::aes(
          y = .data$cumulative_var * scale_factor,
          color = .data$cumulative_series,
          group = 1
        ),
        linewidth = 0.9
      ) +
      ggplot2::geom_point(
        ggplot2::aes(
          y = .data$cumulative_var * scale_factor,
          color = .data$cumulative_series
        ),
        size = 1.8
      ) +
      ggplot2::geom_hline(
        data = threshold_data,
        ggplot2::aes(yintercept = .data$scaled_threshold),
        color = cumulative_color,
        linetype = "dashed",
        alpha = 0.4,
        show.legend = FALSE
      ) +
      ggplot2::geom_vline(
        xintercept = recommended_pcs,
        color = selection_color,
        linetype = "dashed",
        linewidth = 0.8,
        show.legend = FALSE
      ) +
      ggplot2::scale_y_continuous(
        name = "Variance explained (%)",
        limits = c(0, y_upper),
        sec.axis = ggplot2::sec_axis(
          transform = ~ . / scale_factor,
          name = "Cumulative variance (%)"
        )
      ) +
      ggplot2::scale_fill_manual(
        values = c("Variance explained" = variance_color)
      ) +
      ggplot2::scale_color_manual(
        values = c("Cumulative" = cumulative_color)
      ) +
      ggplot2::scale_linetype_manual(values = c("Selected PCs" = "dashed")) +
      ggplot2::labs(
        title = title %||% paste0("PCA elbow plot | Selected: ", recommended_pcs, " PCs"),
        subtitle = subtitle,
        x = xlab,
        fill = NULL,
        color = NULL,
        linetype = NULL
      ) +
      ggplot2::theme(
        legend.position = "none",
        legend.box = "horizontal",
        axis.title.y.left = ggplot2::element_text(color = variance_color),
        axis.text.y.left = ggplot2::element_text(color = variance_color),
        axis.title.y.right = ggplot2::element_text(color = cumulative_color),
        axis.text.y.right = ggplot2::element_text(color = cumulative_color)
      )
  )

  attr(p1, "data") <- plot_data
  attr(p1, "recommended_dims") <- recommended_dims
  attr(p1, "recommended_pcs") <- recommended_pcs
  attr(p1, "scale_factor") <- scale_factor
  return(p1)
}
