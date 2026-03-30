#' @title Estimate useful dimensions from a reduction
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt A `Seurat` object.
#' @param reduction Name of the dimensional reduction to inspect.
#' @param reduction_method Optional reduction method name.
#' When set to `"nmf"` or `"glmpca"`, all available dimensions will be retained.
#' @param k Number of neighbors used by [intrinsicDimension::maxLikGlobalDimEst].
#' Default is `20`.
#' @param min_dims Minimum number of dimensions kept when intrinsic-dimension estimation succeeds.
#' Default is `10`.
#' @param fallback_max_dims Maximum number of dimensions kept when no valid estimate is available.
#' Default is `50`.
#' @param skip_first Whether to drop the first dimension from the returned result.
#' Useful for `TFIDF/LSI` workflows. Default is `FALSE`.
#' @param use_stored Whether to use `misc$dims_estimate` already stored in the reduction when available.
#' Default is `TRUE`.
#'
#' @return An integer vector of dimensions to use.
#'
#' @export
RunDimsEstimate <- function(
  srt,
  reduction,
  reduction_method = NULL,
  k = 20L,
  min_dims = 10L,
  fallback_max_dims = 50L,
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
  check_r("intrinsicDimension", verbose = FALSE)
  if (!is.character(reduction) || length(reduction) != 1L) {
    log_message(
      "{.arg reduction} must be a single reduction name",
      message_type = "error"
    )
  }
  if (!reduction %in% names(srt@reductions)) {
    log_message(
      "{.val {reduction}} is not in the reduction names",
      message_type = "error"
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

  dims_source <- NULL
  dims_use <- integer()
  if (isTRUE(use_stored)) {
    dims_use <- dims_estimate_validate(
      dims_use = reduction_obj@misc[["dims_estimate"]],
      reduction_ncol = reduction_ncol,
      skip_first = skip_first
    )
    if (length(dims_use) > 0L) {
      dims_source <- "stored"
    }
  }

  reduction_method <- tolower(reduction_method %||% "")
  if (is.null(dims_source)) {
    if (reduction_method %in% c("glmpca", "nmf")) {
      dims_use <- dims_estimate_validate(
        dims_use = seq_len(reduction_ncol),
        reduction_ncol = reduction_ncol,
        skip_first = skip_first
      )
      dims_source <- "all"
    } else {
      dim_est <- tryCatch(
        expr = {
          min(
            intrinsicDimension::maxLikGlobalDimEst(
              data = Embeddings(srt, reduction = reduction),
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

      if (!is.na(dim_est)) {
        dims_use <- dims_estimate_validate(
          dims_use = seq_len(max(
            min(reduction_ncol, as.integer(min_dims)),
            ceiling(dim_est)
          )),
          reduction_ncol = reduction_ncol,
          skip_first = skip_first
        )
        dims_source <- "estimated"
      }
    }
  }

  if (length(dims_use) == 0L) {
    dims_use <- dims_estimate_fallback(
      reduction_ncol = reduction_ncol,
      fallback_max_dims = as.integer(fallback_max_dims),
      skip_first = skip_first
    )
    dims_source <- "fallback"
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
    } else if (dims_source == "estimated") {
      log_message(
        "Estimated dimensions {.val {min(dims_use)}}:{.val {max(dims_use)}} for {.pkg {reduction}}",
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

dims_estimate_fallback <- function(
  reduction_ncol,
  fallback_max_dims = 50L,
  skip_first = FALSE
) {
  lower_bound <- if (isTRUE(skip_first)) 2L else 1L
  if (reduction_ncol < lower_bound) {
    return(integer())
  }
  seq.int(lower_bound, min(reduction_ncol, fallback_max_dims))
}

pc_selection_stats <- function(
  srt,
  reduction = "pca",
  max_pcs = 50,
  variance_thresholds = c(0.60, 0.70, 0.80, 0.90)
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  if (!reduction %in% names(srt@reductions)) {
    log_message(
      "{.val {reduction}} is not in the reduction names",
      message_type = "error"
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
#' @param srt A `Seurat` object with PCA reduction computed.
#' @param max_pcs Maximum number of PCs to visualize. Default is `50`.
#' @param variance_thresholds Numeric vector of variance thresholds to mark.
#' Default is `c(0.60, 0.70, 0.80, 0.90)`.
#' @param reduction Reduction name to inspect. Default is `"pca"`.
#' @param palette Palette used for the main curves. Default is `"Chinese"`.
#' @param palcolor Optional palette colors.
#' @param aspect.ratio Aspect ratio of each panel. Default is `NULL`.
#' @param title Title for the combined plot. When `NULL` (default), an
#' auto-generated summary line is used.
#' @param subtitle Subtitle for the combined plot. Default is `NULL`.
#' @param xlab X-axis label shared by all panels. Default is
#' `"Principal component"`.
#' @param theme_use Theme function used to style the plot. 
#' Default is `"theme_scop"`.
#' @param theme_args Other arguments passed to the `theme_use`.
#' @param combine Whether to combine the four panels into one plot. Default is
#' `TRUE`. When `FALSE`, returns a named list of ggplot objects.
#' @param nrow Number of rows in the combined layout. Default is `NULL`.
#' @param ncol Number of columns in the combined layout. Default is `NULL`
#' @param seed Random seed. Default is `11`.
#'
#' @return A patchwork plot object when `combine = TRUE`,
#' or a named list of ggplot objects when `combine = FALSE`.
#'
#' @export
DimsEstimatePlot <- function(
  srt,
  max_pcs = 50,
  variance_thresholds = c(0.60, 0.70, 0.80, 0.90),
  reduction = "pca",
  palette = "Chinese",
  palcolor = NULL,
  aspect.ratio = NULL,
  title = NULL,
  subtitle = NULL,
  xlab = "Principal component",
  theme_use = "theme_scop",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  seed = 11
) {
  set.seed(seed)
  stats_use <- pc_selection_stats(
    srt = srt,
    reduction = reduction,
    max_pcs = max_pcs,
    variance_thresholds = variance_thresholds
  )
  plot_data <- stats_use[["plot_data"]]
  annotation_data <- stats_use[["annotation_data"]]
  threshold_point <- stats_use[["threshold_point"]]
  elbow_point <- stats_use[["elbow_point"]]

  recommended_dims <- RunDimsEstimate(
    srt = srt,
    reduction = reduction,
    reduction_method = reduction,
    use_stored = TRUE,
    verbose = FALSE
  )
  recommended_pcs <- max(recommended_dims)

  main_colors <- palette_colors(
    type = "discrete",
    palette = palette,
    n = 9,
    palcolor = palcolor
  )
  threshold_colors <- palette_colors(
    seq_len(max(1L, nrow(annotation_data))),
    palette = "Chinese"
  )
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

  p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = PC, y = stdev)) +
    ggplot2::geom_line(color = main_colors[1], linewidth = 0.9) +
    ggplot2::geom_point(color = main_colors[2], size = 1.6) +
    ggplot2::geom_vline(
      xintercept = recommended_pcs,
      color = "#1F1F1F",
      linetype = "longdash",
      linewidth = 0.7
    )

  for (i in seq_len(nrow(annotation_data))) {
    if (!is.na(annotation_data$pc[i]) && annotation_data$pc[i] <= max_pcs) {
      p1 <- p1 +
        ggplot2::geom_vline(
          xintercept = annotation_data$pc[i],
          linetype = "dashed",
          color = threshold_colors[i],
          alpha = 0.55
        ) +
        ggplot2::annotate(
          "text",
          x = annotation_data$pc[i],
          y = max(plot_data$stdev) * (1 - i * 0.11),
          label = paste0(annotation_data$threshold[i], "%"),
          color = threshold_colors[i],
          size = 3.1,
          hjust = -0.15
        )
    }
  }

  p1 <- apply_theme(
    p1 +
      ggplot2::labs(
        title = "Elbow",
        x = xlab,
        y = "Standard deviation"
      )
  )

  p2 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = PC, y = cumulative_var)) +
    ggplot2::geom_line(color = main_colors[3], linewidth = 0.9) +
    ggplot2::geom_point(color = main_colors[4], size = 1.4) +
    ggplot2::geom_vline(
      xintercept = recommended_pcs,
      color = "#1F1F1F",
      linetype = "longdash",
      linewidth = 0.7
    )
  for (i in seq_along(variance_thresholds)) {
    p2 <- p2 +
      ggplot2::geom_hline(
        yintercept = variance_thresholds[i] * 100,
        linetype = "dashed",
        color = threshold_colors[i],
        alpha = 0.55
      )
  }
  p2 <- apply_theme(
    p2 +
      ggplot2::scale_y_continuous(limits = c(0, 100)) +
      ggplot2::labs(
        title = "Cumulative variance",
        x = xlab,
        y = "Cumulative variance (%)"
      )
  )

  p3 <- apply_theme(
    ggplot2::ggplot(plot_data, ggplot2::aes(x = PC, y = marginal_gain)) +
      ggplot2::geom_col(fill = main_colors[5], alpha = 0.85, width = 0.85) +
      ggplot2::geom_hline(
        yintercept = 0.5,
        linetype = "dashed",
        color = "#B22222"
      ) +
      ggplot2::geom_vline(
        xintercept = recommended_pcs,
        color = "#1F1F1F",
        linetype = "longdash",
        linewidth = 0.7
      ) +
      ggplot2::labs(
        title = "Marginal gain",
        x = xlab,
        y = "Variance gained (%)"
      )
  )

  p4 <- apply_theme(
    ggplot2::ggplot(
      plot_data[
        seq.int(min(3L, nrow(plot_data)), nrow(plot_data)),
        ,
        drop = FALSE
      ],
      ggplot2::aes(x = PC, y = abs(curvature))
    ) +
      ggplot2::geom_line(color = main_colors[6], linewidth = 0.9) +
      ggplot2::geom_point(color = main_colors[7], size = 1.4) +
      ggplot2::geom_vline(
        xintercept = recommended_pcs,
        color = "#1F1F1F",
        linetype = "longdash",
        linewidth = 0.7
      ) +
      ggplot2::labs(
        title = "Curvature",
        x = xlab,
        y = "Absolute curvature"
      )
  )

  plist <- list(elbow = p1, cumvar = p2, marginal = p3, curvature = p4)

  if (isTRUE(combine)) {
    plot <- patchwork::wrap_plots(
      plotlist = plist,
      nrow = nrow,
      ncol = ncol %||% 2L
    )
    auto_title <- title %||%
      paste0(
        "Recommended PCs: ",
        recommended_pcs,
        if (!is.na(threshold_point)) {
          paste0(" | marginal <0.5% at ", threshold_point)
        } else {
          ""
        }
      )
    plot <- plot +
      patchwork::plot_annotation(
        title = auto_title,
        subtitle = subtitle,
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold"),
          plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
        )
      )
    attr(plot, "data") <- plot_data
    attr(plot, "recommended_dims") <- recommended_dims
    attr(plot, "recommended_pcs") <- recommended_pcs
    return(plot)
  }

  return(plist)
}
