#' @title Proportion Test Plot
#'
#' @description
#' Generate differential-abundance plots based on results from [RunProportionTest].
#' Supports both legacy storage and method-layer storage from the multi-method
#' proportion-test workflow.
#'
#' @md
#' @inheritParams CellDimPlot
#' @inheritParams scop-params
#' @param srt Deprecated alias for `object`; supply exactly one of the two. It
#' will be removed in scop 1.0.0.
#' @param object A Seurat object containing proportion-test results.
#' @param comparison Which comparison to plot.
#' If `NULL`, plots all comparisons.
#' @param proportion_method Optional method to select from
#' `srt@tools[['ProportionTest']][['methods']]`.
#' If `NULL`, uses the active/most recent method.
#' @param result_level Result level to draw. Use `"group"` for group-level
#' results or `"neighborhood"` for stored Milo neighborhood-level results.
#' A missing requested level is an error. Milo group results are sample-level
#' cell-type proportion tests, distinct from the Milo neighborhood test.
#' @param plot_type Plot type. One of `"effect"` or `"umap"`.
#' @param umap_mode UMAP projection mode for `plot_type = "umap"`.
#' `"discrete"` maps cells to DA direction categories;
#' `"continuous"` maps group-level `obs_log2FD` or the mean effect of
#' overlapping Milo neighborhoods to cells. Neighborhood directions are based
#' on significant memberships: opposite directions are `"Mixed"`, while cells
#' in no neighborhood are `"Uncovered"`. This is a member-cell projection,
#' not a separate Milo test per cell.
#' @param reduction Reduction name used by UMAP projection.
#' @param projection_args Additional arguments passed to [CellDimPlot]
#' (`umap_mode = "discrete"`) or [FeatureDimPlot]
#' (`umap_mode = "continuous"`).
#' @param FDR_threshold FDR value cutoff for frequentist methods. Milo
#' neighborhood results use `SpatialFDR`. scCODA uses its stored `credible`
#' decision instead.
#' @param log2FD_threshold Absolute value of log2FD cutoff for frequentist
#' methods; not applied to scCODA credibility.
#' @param order_by Method to order clusters.
#' Options: `"name"` (alphabetical), `"value"` (by log2FD value).
#' @param palette Color palette name for continuous effect coloring.
#' @param palcolor Custom colors for `palette`.
#' @param group_palette Palette for cluster/group coloring.
#' @param group_palcolor Custom colors for `group_palette`.
#' @param pt.size The size of the points.
#' @param pt.alpha Point transparency.
#' @param cols.sig Color for significant/credible points and intervals.
#' @param cols.ns Color for non-significant points and intervals.
#' @param cols.increase Default color for increased DA groups.
#' @param cols.decrease Default color for decreased DA groups.
#' @param effect_color_mode Coloring mode for `plot_type = "effect"`.
#' Use `"directional"` (default) for increased/decreased/NS colors,
#' or `"classic"` for legacy significant/non-significant coloring.
#' @param nlabel Number of labels added when `label = TRUE` and
#' `features_label = NULL`.
#' @param features_label Character vector specifying points to label.
#' @param label Whether to add labels.
#' @param label.fg Label foreground color.
#' @param label.bg Label background color.
#' @param label.bg.r Label background radius.
#' @param label.size Label text size.
#' @param aspect.ratio Aspect ratio of the panel.
#' @param xlab X-axis label. For
#' `plot_type = "umap"`, this is forwarded to the projection plot when set.
#' @param ylab Y-axis label. For
#' `plot_type = "umap"`, this is forwarded to the projection plot when set.
#' @param legend.position The position of legends,
#' one of `"none"`, `"left"`, `"right"`, `"bottom"`, `"top"`.
#' @param legend.title Title of the legend.
#'
#' @seealso [RunProportionTest]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#' pancreas_sub <- RunProportionTest(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   split.by = "Phase",
#'   proportion_method = "permutation"
#' )
#'
#' ProportionTestPlot(pancreas_sub)
#'
#' ProportionTestPlot(
#'   pancreas_sub,
#'   reduction = "UMAP",
#'   plot_type = "umap",
#'   xlab = "UMAP_1",
#'   ylab = "UMAP_2"
#' )
ProportionTestPlot <- function(
  object,
  comparison = NULL,
  proportion_method = NULL,
  result_level = c("group", "neighborhood"),
  plot_type = c("effect", "umap"),
  umap_mode = c("discrete", "continuous"),
  reduction = "UMAP",
  projection_args = list(),
  FDR_threshold = 0.05,
  log2FD_threshold = log2(1.5),
  order_by = c("value", "name"),
  palette = "RdBu",
  palcolor = NULL,
  group_palette = "Chinese",
  group_palcolor = NULL,
  pt.size = 1,
  pt.alpha = 1,
  cols.sig = "red",
  cols.ns = "grey",
  cols.increase = "#d7301f",
  cols.decrease = "#2b8cbe",
  effect_color_mode = c("directional", "classic"),
  nlabel = 5,
  features_label = NULL,
  label = FALSE,
  label.fg = "black",
  label.bg = "white",
  label.bg.r = 0.1,
  label.size = 4,
  aspect.ratio = NULL,
  xlab = "Cell Type",
  ylab = "log2 (FD)",
  theme_use = "theme_scop",
  theme_args = list(),
  legend.position = "bottom",
  legend.direction = "vertical",
  legend.title = "Significance",
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  seed = 11,
  verbose = TRUE,
  srt = NULL
) {
  srt <- resolve_deprecated_srt(object, srt, missing(object))
  order_by <- match.arg(order_by)
  plot_type <- match.arg(plot_type)
  umap_mode <- match.arg(umap_mode)
  effect_color_mode <- match.arg(effect_color_mode)

  target_level <- match.arg(result_level)

  resolved <- get_proportion_plot_results(
    srt = srt,
    proportion_method = proportion_method,
    result_level = target_level
  )
  results <- resolved$results
  method_used <- resolved$method
  method_bundle <- resolved$method_bundle
  result_level <- resolved$result_level

  if (!is.null(comparison)) {
    comparison <- normalize_plot_comparison_input(comparison)
    missing_comparisons <- comparison[!comparison %in% names(results)]
    if (length(missing_comparisons) > 0) {
      log_message(
        "Specified comparisons not found in results: {.val {missing_comparisons}}",
        message_type = "error"
      )
    }
    results <- results[comparison]
  }

  if (is.null(results) || length(results) == 0) {
    log_message(
      "No proportion test results found for plotting",
      message_type = "error"
    )
  }

  std_results <- list()
  for (comp_name in names(results)) {
    comp_groups <- plot_comparison_groups(comp_name)
    plot_data <- standardize_proportion_result(
      results[[comp_name]],
      cluster_1 = comp_groups[1],
      cluster_2 = comp_groups[2],
      comparison_name = comp_name,
      method = method_used
    )
    plot_data <- prepare_proportion_plot_data(
      plot_data,
      FDR_threshold = FDR_threshold,
      log2FD_threshold = log2FD_threshold,
      order_by = order_by,
      nlabel = nlabel,
      features_label = features_label,
      label = label,
      credible_threshold = method_bundle[["parameters"]][["credible_effect_threshold"]] %||% 0.95,
      seed = seed
    )
    std_results[[comp_name]] <- plot_data
  }

  plot_args <- list(
    FDR_threshold = FDR_threshold,
    log2FD_threshold = log2FD_threshold,
    palette = palette,
    palcolor = palcolor,
    group_palette = group_palette,
    group_palcolor = group_palcolor,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    cols.sig = cols.sig,
    cols.ns = cols.ns,
    cols.increase = cols.increase,
    cols.decrease = cols.decrease,
    effect_color_mode = effect_color_mode,
    label = label,
    label.fg = label.fg,
    label.bg = label.bg,
    label.bg.r = label.bg.r,
    label.size = label.size,
    xlab = if (identical(result_level, "neighborhood") && identical(xlab, "Cell Type")) {
      "Milo neighborhood"
    } else {
      xlab
    },
    ylab = if (identical(result_level, "neighborhood") && identical(ylab, "log2 (FD)")) {
      "Milo log2 FC"
    } else {
      ylab
    },
    aspect.ratio = aspect.ratio,
    legend.position = legend.position,
    legend.direction = legend.direction,
    legend.title = legend.title,
    theme_use = theme_use,
    theme_args = theme_args
  )

  plist <- switch(plot_type,
    effect = lapply(std_results, function(df) do.call(plot_proportion_effect, c(list(df = df), plot_args))),
    umap = plot_proportion_umap(
      srt = srt,
      std_results = std_results,
      method_bundle = method_bundle,
      result_level = result_level,
      umap_mode = umap_mode,
      reduction = reduction,
      projection_args = projection_args,
      FDR_threshold = FDR_threshold,
      log2FD_threshold = log2FD_threshold,
      cols.increase = cols.increase,
      cols.decrease = cols.decrease,
      cols.ns = cols.ns,
      palette = palette,
      palcolor = palcolor,
      xlab = if (identical(xlab, "Cell Type")) NULL else xlab,
      ylab = if (identical(ylab, "log2 (FD)")) NULL else ylab,
      legend.title = legend.title,
      theme_use = theme_use,
      theme_args = theme_args
    )
  )

  if (!isTRUE(combine)) {
    return(plist)
  }

  combine_proportion_plots(
    plist = plist,
    nrow = nrow,
    ncol = ncol,
    byrow = byrow,
    legend.position = legend.position,
    legend.direction = legend.direction,
    pair_reverse = identical(plot_type, "effect")
  )
}

prepare_proportion_plot_data <- function(
  plot_data,
  FDR_threshold,
  log2FD_threshold,
  order_by,
  nlabel,
  features_label,
  label,
  credible_threshold = 0.95,
  seed = 11
) {
  plot_data$clusters <- as.character(plot_data$clusters)
  plot_data$FDR <- suppressWarnings(as.numeric(plot_data$FDR))
  plot_data$pval <- suppressWarnings(as.numeric(plot_data$pval))
  plot_data$obs_log2FD <- suppressWarnings(as.numeric(plot_data$obs_log2FD))
  plot_data$effect <- suppressWarnings(as.numeric(plot_data$effect))
  plot_data$boot_CI_2.5 <- suppressWarnings(as.numeric(plot_data$boot_CI_2.5))
  plot_data$boot_CI_97.5 <- suppressWarnings(as.numeric(plot_data$boot_CI_97.5))
  plot_data$hdi_2.5 <- suppressWarnings(as.numeric(plot_data$hdi_2.5))
  plot_data$hdi_97.5 <- suppressWarnings(as.numeric(plot_data$hdi_97.5))
  plot_data$inclusion_prob <- suppressWarnings(as.numeric(plot_data$inclusion_prob))

  method <- unique(as.character(plot_data$method))
  is_sccoda <- length(method) == 1L && identical(method, "sccoda")
  not_tested <- if ("inference_valid" %in% colnames(plot_data)) {
    !is.na(plot_data$inference_valid) & !as.logical(plot_data$inference_valid)
  } else {
    rep(FALSE, nrow(plot_data))
  }
  if (is_sccoda) {
    credible <- as.logical(plot_data$credible)
    missing_credible <- is.na(credible)
    credible[missing_credible] <- !is.na(plot_data$inclusion_prob[missing_credible]) &
      plot_data$inclusion_prob[missing_credible] >= credible_threshold
    significant <- credible & is.finite(plot_data$obs_log2FD) & plot_data$obs_log2FD != 0
    sig_label <- "Credible effect"
    stat_p <- pmax(1 - plot_data$inclusion_prob, .Machine$double.eps)
  } else {
    significant <- !is.na(plot_data$FDR) & is.finite(plot_data$obs_log2FD) &
      plot_data$FDR < FDR_threshold &
      abs(plot_data$obs_log2FD) > log2FD_threshold
    sig_label <- proportion_sig_label(FDR_threshold, log2FD_threshold)
    stat_p <- ifelse(is.finite(plot_data$FDR), plot_data$FDR, plot_data$pval)
  }
  plot_data$significance <- ifelse(not_tested, "Not tested",
    ifelse(significant, sig_label, "n.s."))
  plot_data$significance <- factor(plot_data$significance,
    levels = c(sig_label, "n.s.", "Not tested"))

  plot_data$direction <- ifelse(not_tested, "Not tested",
    ifelse(significant & plot_data$obs_log2FD > 0, "Increased",
      ifelse(significant & plot_data$obs_log2FD < 0, "Decreased", "NS")))
  plot_data$direction <- factor(
    plot_data$direction,
    levels = c("Increased", "Decreased", "NS", "Not tested")
  )

  stat_p[!is.finite(stat_p) | stat_p <= 0] <- 1
  plot_data$minus_log10 <- -log10(stat_p)

  if (!"neighborhood" %in% colnames(plot_data)) {
    plot_data$neighborhood <- NA_character_
  }
  plot_data$label_id <- ifelse(
    !is.na(plot_data$neighborhood) & nzchar(plot_data$neighborhood),
    as.character(plot_data$neighborhood),
    as.character(plot_data$clusters)
  )

  if (order_by == "value") {
    clusters_factor <- factor(plot_data$clusters)
    cluster_means <- tapply(plot_data$obs_log2FD, clusters_factor, mean, na.rm = TRUE)
    cluster_order <- names(sort(cluster_means, decreasing = TRUE))
  } else {
    cluster_order <- sort(unique(plot_data$clusters))
  }
  plot_data$clusters <- factor(plot_data$clusters, levels = cluster_order)

  plot_data$label <- FALSE
  if (isTRUE(label)) {
    if (!is.null(features_label)) {
      plot_data$label <- plot_data$label_id %in% features_label
    } else {
      score <- abs(plot_data$obs_log2FD) * pmax(plot_data$minus_log10, 0)
      ord <- order(score, decreasing = TRUE)
      top_idx <- utils::head(ord[is.finite(score[ord])], nlabel)
      plot_data$label[top_idx] <- TRUE
    }
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  plot_data
}

plot_proportion_effect <- function(
  df,
  FDR_threshold,
  log2FD_threshold,
  pt.size,
  pt.alpha,
  cols.sig,
  cols.ns,
  cols.increase,
  cols.decrease,
  effect_color_mode = c("directional", "classic"),
  xlab,
  ylab,
  aspect.ratio,
  legend.position,
  legend.direction,
  legend.title,
  theme_use,
  theme_args,
  ...
) {
  effect_color_mode <- match.arg(effect_color_mode)
  sig_label <- if (identical(unique(as.character(df$method)), "sccoda")) {
    "Credible effect"
  } else {
    proportion_sig_label(FDR_threshold, log2FD_threshold)
  }
  has_ci <- any(is.finite(df$boot_CI_2.5) & is.finite(df$boot_CI_97.5))

  p <- ggplot(df, aes(x = clusters, y = obs_log2FD))
  if (identical(effect_color_mode, "classic")) {
    if (has_ci) {
      p <- p +
        geom_pointrange(
          aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = significance),
          size = pt.size,
          alpha = pt.alpha,
          na.rm = TRUE
        )
    } else {
      p <- p +
        geom_point(
          aes(color = significance),
          size = pt.size * 1.8,
          alpha = pt.alpha,
          na.rm = TRUE
        )
    }

    p <- p +
      scale_color_manual(
        name = legend.title %||% "Significance",
        labels = c(sig_label, "n.s.", "Not tested"),
        values = stats::setNames(c(cols.sig, cols.ns, cols.ns),
          c(sig_label, "n.s.", "Not tested")),
        drop = FALSE
      )
  } else {
    if (has_ci) {
      p <- p +
        geom_pointrange(
          aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = direction),
          size = pt.size,
          alpha = pt.alpha,
          na.rm = TRUE
        )
    } else {
      p <- p +
        geom_point(
          aes(color = direction),
          size = pt.size * 1.8,
          alpha = pt.alpha,
          na.rm = TRUE
        )
    }

    p <- p +
      scale_color_manual(
        name = legend.title %||% "Direction",
        values = c(
          Increased = cols.increase,
          Decreased = cols.decrease,
          NS = cols.ns,
          `Not tested` = cols.ns
        ),
        drop = FALSE
      )
  }

  thresholds <- if (identical(unique(as.character(df$method)), "sccoda")) 0 else
    c(-log2FD_threshold, 0, log2FD_threshold)
  subtitle <- if ("inference_valid" %in% colnames(df) &&
      any(!is.na(df$inference_valid) & !as.logical(df$inference_valid))) {
    "Descriptive virtual samples; no inferential test"
  } else if (identical(unique(as.character(df$method)), "milo")) {
    if (any(!is.na(df$neighborhood))) {
      "Milo neighborhood test; spatial FDR"
    } else {
      "Sample-level cell-type summary"
    }
  } else if (identical(unique(as.character(df$method)), "propeller")) {
    "Internal logit-transformed sample-proportion test"
  } else if (identical(unique(as.character(df$method)), "sccoda")) {
    "scCODA credible effects"
  } else {
    NULL
  }
  p +
    geom_hline(yintercept = thresholds,
      linetype = if (length(thresholds) == 1L) 1 else c(2, 1, 2),
      color = if (length(thresholds) == 1L) "black" else c("grey", "black", "grey")) +
    labs(
      title = proportion_title(df),
      subtitle = subtitle,
      x = xlab,
      y = ylab
    ) +
    coord_flip() +
    apply_plot_theme(theme_use, theme_args) +
    theme(
      aspect.ratio = aspect.ratio,
      legend.position = legend.position,
      legend.direction = legend.direction
    )
}

plot_proportion_umap <- function(
  srt,
  std_results,
  method_bundle,
  result_level = "group",
  umap_mode = c("discrete", "continuous"),
  reduction = "UMAP",
  projection_args = list(),
  FDR_threshold = 0.05,
  log2FD_threshold = log2(1.5),
  cols.increase = "#d7301f",
  cols.decrease = "#2b8cbe",
  cols.ns = "grey80",
  palette = "RdBu",
  palcolor = NULL,
  xlab = NULL,
  ylab = NULL,
  legend.title = NULL,
  theme_use = "theme_scop",
  theme_args = list()
) {
  umap_mode <- match.arg(umap_mode)
  reduction_use <- if (is.null(reduction)) {
    DefaultReduction(srt)
  } else {
    DefaultReduction(srt, pattern = reduction)
  }

  if (identical(result_level, "group")) {
    group.by <- method_bundle[["parameters"]][["group.by"]] %||%
      srt@tools[["ProportionTest"]][["parameters"]][["group.by"]]
    if (is.null(group.by) || !nzchar(group.by) ||
        !group.by %in% colnames(srt@meta.data)) {
      log_message(
        "Cannot determine a valid {.arg group.by} from proportion test metadata for {.val plot_type = 'umap'}",
        message_type = "error"
      )
    }
    cluster_by_cell <- as.character(srt@meta.data[[group.by]])
  } else {
    members <- method_bundle[["details"]][["milo_graph_data"]][[".metadata"]][["members"]]
    if (is.null(members) || length(members) == 0L) {
      log_message(
        "Milo neighborhood membership is unavailable; rerun {.fn RunProportionTest} before projecting neighborhood results",
        message_type = "error"
      )
    }
  }

  projection_args <- projection_args %||% list()
  legend_discrete <- if (is.null(legend.title) || identical(legend.title, "Significance")) {
    "DA Direction"
  } else {
    legend.title
  }
  legend_continuous <- if (is.null(legend.title) || identical(legend.title, "Significance")) {
    if (identical(result_level, "neighborhood")) "Mean neighborhood log2 FC" else "log2 (FD)"
  } else {
    legend.title
  }
  plist <- list()

  for (comp_name in names(std_results)) {
    df_comp <- std_results[[comp_name]]
    if (identical(result_level, "neighborhood")) {
      projected <- project_milo_neighborhood_to_cells(
        df_comp, members, rownames(srt@meta.data)
      )
      cell_direction <- projected$direction
      cell_effect <- projected$effect
    } else {
      proj_df <- summarize_proportion_projection(
        df = df_comp,
        FDR_threshold = FDR_threshold,
        log2FD_threshold = log2FD_threshold
      )
      effect_map <- stats::setNames(proj_df$obs_log2FD, proj_df$clusters)
      direction_map <- stats::setNames(as.character(proj_df$direction), proj_df$clusters)
      cell_direction <- direction_map[cluster_by_cell]
      cell_direction[is.na(cell_direction)] <- "Uncovered"
      cell_effect <- effect_map[cluster_by_cell]
    }

    suffix <- proportion_projection_suffix(comp_name)
    direction_col <- paste0(".proportion_da_direction_", suffix)
    effect_col <- paste0(".proportion_da_log2fd_", suffix)

    srt@meta.data[[direction_col]] <- cell_direction
    srt@meta.data[[direction_col]] <- factor(
      srt@meta.data[[direction_col]],
      levels = c("Increased", "Decreased", "Mixed", "NS", "Not tested", "Uncovered")
    )

    srt@meta.data[[effect_col]] <- cell_effect

    title_use <- proportion_title(df_comp)
    if (identical(umap_mode, "discrete")) {
      plot_call <- utils::modifyList(
        projection_args,
        list(
          object = srt,
          group.by = direction_col,
          reduction = reduction_use,
          palette = "Set1",
          palcolor = c(cols.increase, cols.decrease, "#7b3294", cols.ns,
            "grey65", "grey95"),
          show_stat = FALSE,
          title = title_use,
          subtitle = if (identical(result_level, "neighborhood")) {
            "Overlapping Milo neighborhoods projected to member cells"
          } else NULL,
          xlab = xlab,
          ylab = ylab,
          legend.title = legend_discrete,
          theme_use = theme_use,
          theme_args = theme_args,
          combine = TRUE
        )
      )
      p <- invoke_fun(
        CellDimPlot,
        plot_call[names(plot_call) %in% names(formals(CellDimPlot))]
      )
    } else {
      plot_call <- utils::modifyList(
        projection_args,
        list(
          object = srt,
          features = effect_col,
          reduction = reduction_use,
          palette = palette,
          palcolor = palcolor,
          bg_cutoff = -Inf,
          show_stat = FALSE,
          title = title_use,
          subtitle = if (identical(result_level, "neighborhood")) {
            "Mean effect of overlapping Milo neighborhoods"
          } else NULL,
          xlab = xlab,
          ylab = ylab,
          legend.title = legend_continuous,
          theme_use = theme_use,
          theme_args = theme_args,
          combine = TRUE
        )
      )
      p <- invoke_fun(
        FeatureDimPlot,
        plot_call[names(plot_call) %in% names(formals(FeatureDimPlot))]
      )
    }

    plist[[comp_name]] <- p
  }

  plist
}

project_milo_neighborhood_to_cells <- function(df, members, cells) {
  if (is.null(names(members)) || any(!nzchar(names(members)))) {
    log_message("Milo neighborhood membership must have neighborhood IDs", message_type = "error")
  }
  neighborhood <- as.character(df$neighborhood)
  if (anyNA(neighborhood) || any(!neighborhood %in% names(members))) {
    log_message(
      "Milo neighborhood result IDs do not match stored cell membership",
      message_type = "error"
    )
  }
  effect_sum <- numeric(length(cells))
  effect_count <- integer(length(cells))
  covered <- rep(FALSE, length(cells))
  increased <- rep(FALSE, length(cells))
  decreased <- rep(FALSE, length(cells))
  not_tested <- rep(FALSE, length(cells))
  for (i in seq_len(nrow(df))) {
    cell_idx <- match(as.character(members[[neighborhood[[i]]]]), cells)
    cell_idx <- unique(cell_idx[!is.na(cell_idx)])
    if (length(cell_idx) == 0L) next
    covered[cell_idx] <- TRUE
    value <- as.numeric(df$obs_log2FD[[i]])
    if (is.finite(value)) {
      effect_sum[cell_idx] <- effect_sum[cell_idx] + value
      effect_count[cell_idx] <- effect_count[cell_idx] + 1L
    }
    direction <- as.character(df$direction[[i]])
    if (identical(direction, "Increased")) increased[cell_idx] <- TRUE
    if (identical(direction, "Decreased")) decreased[cell_idx] <- TRUE
    if (identical(direction, "Not tested")) not_tested[cell_idx] <- TRUE
  }
  if (!any(covered)) {
    log_message(
      "No Milo neighborhood members match cells in the object",
      message_type = "error"
    )
  }
  effect <- rep(NA_real_, length(cells))
  effect[effect_count > 0L] <- effect_sum[effect_count > 0L] / effect_count[effect_count > 0L]
  direction <- rep("NS", length(cells))
  direction[not_tested] <- "Not tested"
  direction[increased] <- "Increased"
  direction[decreased] <- "Decreased"
  direction[increased & decreased] <- "Mixed"
  direction[!covered] <- "Uncovered"
  names(effect) <- cells
  names(direction) <- cells
  list(effect = effect, direction = direction)
}

summarize_proportion_projection <- function(
  df,
  FDR_threshold,
  log2FD_threshold
) {
  if (!"clusters" %in% colnames(df)) {
    log_message(
      "Cannot map UMAP projection without {.field clusters} in proportion results",
      message_type = "error"
    )
  }

  tmp <- data.frame(
    clusters = as.character(df$clusters),
    obs_log2FD = suppressWarnings(as.numeric(df$obs_log2FD)),
    FDR = suppressWarnings(as.numeric(df$FDR)),
    direction = as.character(df$direction),
    stringsAsFactors = FALSE
  )
  tmp <- tmp[!is.na(tmp$clusters) & nzchar(tmp$clusters), , drop = FALSE]
  if (nrow(tmp) == 0) {
    return(data.frame(
      clusters = character(0),
      obs_log2FD = numeric(0),
      FDR = numeric(0),
      direction = factor(character(0), levels = c("Increased", "Decreased", "NS")),
      stringsAsFactors = FALSE
    ))
  }

  cluster_order <- unique(tmp$clusters)
  merged <- do.call(
    rbind,
    lapply(split(tmp, tmp$clusters), function(x) {
      effect <- if (all(is.na(x$obs_log2FD))) NA_real_ else mean(x$obs_log2FD, na.rm = TRUE)
      fdr <- if (all(is.na(x$FDR))) NA_real_ else min(x$FDR, na.rm = TRUE)
      direction <- if (all(c("Increased", "Decreased") %in% x$direction)) "Mixed" else
        if ("Increased" %in% x$direction) "Increased" else
          if ("Decreased" %in% x$direction) "Decreased" else
            if ("Not tested" %in% x$direction) "Not tested" else "NS"

      data.frame(
        clusters = x$clusters[1],
        obs_log2FD = effect,
        FDR = fdr,
        direction = direction,
        stringsAsFactors = FALSE
      )
    })
  )

  merged$clusters <- factor(merged$clusters, levels = cluster_order)
  merged <- merged[order(merged$clusters), , drop = FALSE]
  merged$clusters <- as.character(merged$clusters)
  merged$direction <- factor(merged$direction,
    levels = c("Increased", "Decreased", "Mixed", "NS", "Not tested", "Uncovered"))
  rownames(merged) <- NULL
  merged
}

proportion_projection_suffix <- function(comparison_name) {
  out <- gsub("[^A-Za-z0-9]+", "_", comparison_name)
  out <- gsub("^_+|_+$", "", out)
  if (!nzchar(out)) {
    out <- "comparison"
  }
  out
}

proportion_sig_label <- function(FDR_threshold, log2FD_threshold) {
  paste("FDR <", FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2))
}

proportion_title <- function(df) {
  if (!is.null(df$group1) && !is.null(df$group2) && !all(is.na(df$group1)) && !all(is.na(df$group2))) {
    paste0(df$group1[1], " vs ", df$group2[1])
  } else if (!is.null(df$comparison) && !all(is.na(df$comparison))) {
    as.character(df$comparison[1])
  } else {
    "Proportion Test"
  }
}

combine_proportion_plots <- function(
  plist,
  nrow,
  ncol,
  byrow,
  legend.position,
  legend.direction,
  pair_reverse = FALSE
) {
  if (length(plist) == 0) {
    return(plist)
  }

  if (length(plist) == 1) {
    return(plist[[1]])
  }

  check_r("patchwork", verbose = FALSE)

  plot_names <- names(plist)
  order_idx <- seq_along(plist)

  if (isTRUE(pair_reverse)) {
    paired_groups <- list()
    used <- rep(FALSE, length(plot_names))

    for (i in seq_along(plot_names)) {
      if (used[i]) {
        next
      }
      groups <- plot_comparison_groups(plot_names[i])
      if (all(!is.na(groups))) {
        reverse_name <- paste0(groups[2], "_vs_", groups[1])
        j <- which(plot_names == reverse_name)[1]
        if (!is.na(j)) {
          paired_groups[[length(paired_groups) + 1]] <- c(i, j)
          used[c(i, j)] <- TRUE
        }
      }
    }

    remaining <- which(!used)
    order_idx <- c(unlist(paired_groups), remaining)
  }

  reordered <- plist[order_idx]

  if (is.null(nrow) && is.null(ncol)) {
    n_plots <- length(reordered)
    ncol <- ceiling(sqrt(n_plots))
    nrow <- ceiling(n_plots / ncol)
  }

  patchwork::wrap_plots(
    plotlist = reordered,
    nrow = nrow,
    ncol = ncol,
    byrow = byrow
  ) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      theme = theme(
        plot.margin = margin(10, 10, 10, 10),
        legend.position = legend.position,
        legend.direction = legend.direction
      )
    )
}

get_proportion_plot_results <- function(
  srt,
  proportion_method = NULL,
  result_level = c("group", "neighborhood")
) {
  result_level <- match.arg(result_level)

  if (!"ProportionTest" %in% names(srt@tools)) {
    log_message(
      "Cannot find the ProportionTest result. Perform {.fn RunProportionTest} first",
      message_type = "error"
    )
  }

  pt <- srt@tools[["ProportionTest"]]
  methods_store <- pt[["methods"]]

  if (is.null(methods_store) || length(methods_store) == 0) {
    if (identical(result_level, "neighborhood")) {
      log_message("Legacy proportion results do not contain neighborhood results", message_type = "error")
    }
    results <- pt[["results"]]
    if (is.null(results) || length(results) == 0) {
      log_message(
        "No proportion test results found",
        message_type = "error"
      )
    }
    return(list(results = results, method = "legacy", result_level = "group", method_bundle = list()))
  }

  method_use <- proportion_method
  if (!is.null(method_use)) {
    method_use <- normalize_proportion_method(method_use)
  } else {
    method_use <- pt[["active_method"]] %||% pt[["parameters"]][["proportion_method"]]
    if (is.null(method_use) || !method_use %in% names(methods_store)) {
      method_use <- names(methods_store)[length(methods_store)]
    }
  }

  if (!method_use %in% names(methods_store)) {
    log_message(
      "Method {.val {method_use}} is not found in stored proportion test results",
      message_type = "error"
    )
  }

  method_bundle <- methods_store[[method_use]]
  results <- NULL

  if (identical(result_level, "neighborhood")) {
    results <- method_bundle[["neighborhood_results"]] %||%
      method_bundle[["details"]][["neighborhood_results"]]
    if (is.null(results) || length(results) == 0) {
      log_message(
        "Neighborhood-level results are not available for method {.val {method_use}}",
        message_type = "error"
      )
    }
  }

  if (is.null(results) || length(results) == 0) {
    results <- method_bundle[["results"]]
  }

  if (is.null(results) || length(results) == 0) {
    log_message(
      "No proportion test results found for method {.val {method_use}}",
      message_type = "error"
    )
  }

  list(
    results = results,
    method = method_use,
    result_level = result_level,
    method_bundle = method_bundle
  )
}

normalize_plot_comparison_input <- function(comparison) {
  if (is.list(comparison)) {
    if (all(vapply(comparison, function(x) is.character(x) && length(x) == 2, logical(1)))) {
      comparison_names <- c()
      for (pair in comparison) {
        comparison_names <- c(
          comparison_names,
          paste0(pair[1], "_vs_", pair[2]),
          paste0(pair[2], "_vs_", pair[1])
        )
      }
      return(unique(comparison_names))
    }
    return(unique(as.character(unlist(comparison))))
  }

  unique(as.character(comparison))
}

plot_comparison_groups <- function(comparison_name) {
  groups <- strsplit(comparison_name, "_vs_", fixed = TRUE)[[1]]
  if (length(groups) != 2) {
    return(c(NA_character_, NA_character_))
  }
  groups
}
