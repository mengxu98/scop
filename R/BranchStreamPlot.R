#' @title Branch Stream Plot
#'
#' @description
#' Draw branch-aware pseudotime ribbons from cell-state annotations. The plot is
#' a visualization of pseudotime density by group; it does not infer a lineage
#' tree by itself.
#'
#' @md
#' @inheritParams CellDimPlot
#' @param object A `Seurat` object or a `data.frame` containing cell metadata.
#' @param group.by Column containing cell-state or branch labels.
#' @param lineages Pseudotime columns to use. If `NULL`, lineage-like
#' pseudotime columns such as `Lineage1`, `Lineage2`, `prefix_Lineage1`, or
#' `pseudotime` are detected and merged into one global pseudotime for a
#' single panel. Use `"all"` to plot each detected lineage in separate panels.
#' @param lineage.merge How to merge multiple lineage pseudotime columns when
#' `lineages` is `NULL` or `"merge"`.
#' @param labels Optional order of groups to plot. Factor levels are used when
#' omitted.
#' @param trunk_groups Groups placed on the shared trunk. If `NULL`, groups are
#' inferred from KDE peak positions.
#' @param branch_groups Optional named numeric vector/list mapping group labels
#' to branch amplitudes. Positive and negative amplitudes bend branches in
#' opposite directions.
#' @param branch_center,branch_steepness,branch_power,branch_amplitude
#' Parameters controlling automatic branch centerlines. Use `"auto"` for
#' `branch_center` to place the split point on the current pseudotime scale.
#' @param n_branches Number of branch amplitudes used when `branch_groups` is
#' inferred.
#' @param n_grid Number of pseudotime grid points.
#' @param bw,pad SciPy-style Gaussian KDE bandwidth factor and taper padding.
#' @param count_power Exponent used to scale profiles by group abundance.
#' @param normalize_pseudotime Whether to min-max normalize pseudotime before
#' KDE. Default is `FALSE`, so the original pseudotime scale is shown. Use
#' `"auto"` to normalize only when values fall outside `[0, 1]`.
#' @param scale_to Maximum ribbon thickness after global rescaling. Use
#' `"auto"` for a sensible default or `NULL` to disable rescaling.
#' @param min_visible_width Ribbons thinner than this threshold are hidden.
#' @param ribbon_alpha Ribbon transparency.
#' @param xlim,ylim Plot limits. `xlim = NULL` uses the observed pseudotime
#' range.
#' @param xlabel X-axis label.
#' @param xticks Pseudotime tick locations. Use `"auto"` for pretty breaks
#' on the current pseudotime scale.
#' @param axis_y Y position of the arrow-style pseudotime axis. Use `"auto"` or
#' `NULL`.
#' @param label_positions Optional data.frame/list with label, x, y, and
#' optional size columns. `"auto"` places labels at ribbon maxima.
#' @param label_size,label_fontface,label_fill,label_fill_alpha Label styling.
#' If `label_fill` is `NULL`, labels are drawn as OmicVerse-style colored
#' text with a white outline.
#' @param label_outline_width,label_outline_alpha Outline styling used when
#' `label_fill` is `NULL`.
#' @param palette Color palette name.
#' @param palcolor Custom colors used to create a color palette.
#' @param axis_arrow Arrow used for the pseudotime axis.
#' @param title Plot title. `NULL` hides the title for merged/single panels.
#' When multiple lineages are plotted and `title` is `NULL`, each panel is
#' titled with its lineage column.
#' @param combine Whether to combine multiple lineage plots with `patchwork`.
#' @param nrow,ncol,byrow Layout controls for combined lineage plots.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' branch_df <- data.frame(
#'   cell_type = rep(c("Root", "Branch A", "Branch B"), each = 30),
#'   pseudotime = c(
#'     seq(0, 30, length.out = 30),
#'     seq(25, 80, length.out = 30),
#'     seq(25, 100, length.out = 30)
#'   )
#' )
#' BranchStreamPlot(
#'   branch_df,
#'   group.by = "cell_type",
#'   lineages = "pseudotime"
#' )
BranchStreamPlot <- function(
  object,
  group.by,
  lineages = NULL,
  lineage.merge = c("min", "mean", "max"),
  labels = NULL,
  trunk_groups = NULL,
  branch_groups = NULL,
  branch_center = "auto",
  branch_steepness = 11,
  branch_power = 1.1,
  branch_amplitude = 0.28,
  n_branches = 2,
  n_grid = 800,
  bw = 0.15,
  pad = 0.035,
  count_power = 0.3,
  normalize_pseudotime = FALSE,
  scale_to = "auto",
  min_visible_width = 1e-4,
  xlim = NULL,
  ylim = NULL,
  xlabel = "Pseudotime",
  xticks = "auto",
  axis_y = "auto",
  label_positions = "auto",
  label_size = 4,
  label_fontface = "bold",
  label_fill = NULL,
  label_fill_alpha = 0.85,
  label_outline_width = 0.0018,
  label_outline_alpha = 0.8,
  ribbon_alpha = 1,
  palette = "Chinese",
  palcolor = NULL,
  axis_arrow = grid::arrow(length = grid::unit(0.08, "inches")),
  title = NULL,
  legend.position = "none",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE
) {
  obs <- branch_stream_resolve_obs(object, group.by)
  lineage_info <- branch_stream_resolve_lineages(
    obs,
    lineages = lineages,
    lineage.merge = lineage.merge
  )
  obs <- lineage_info$obs
  lineages <- lineage_info$lineages

  plots <- lapply(lineages, function(lineage) {
    panel_title <- title
    if (is.null(panel_title)) {
      panel_title <- lineage_info$titles[[lineage]]
      if (is.na(panel_title)) {
        panel_title <- NULL
      }
    }
    branch_stream_single_plot(
      obs = obs,
      group.by = group.by,
      pseudotime.by = lineage,
      labels = labels,
      trunk_groups = trunk_groups,
      branch_groups = branch_groups,
      branch_center = branch_center,
      branch_steepness = branch_steepness,
      branch_power = branch_power,
      branch_amplitude = branch_amplitude,
      n_branches = n_branches,
      n_grid = n_grid,
      bw = bw,
      pad = pad,
      count_power = count_power,
      normalize_pseudotime = normalize_pseudotime,
      scale_to = scale_to,
      min_visible_width = min_visible_width,
      ribbon_alpha = ribbon_alpha,
      xlim = xlim,
      ylim = ylim,
      xlabel = xlabel,
      xticks = xticks,
      axis_y = axis_y,
      label_positions = label_positions,
      label_size = label_size,
      label_fontface = label_fontface,
      label_fill = label_fill,
      label_fill_alpha = label_fill_alpha,
      label_outline_width = label_outline_width,
      label_outline_alpha = label_outline_alpha,
      palette = palette,
      palcolor = palcolor,
      axis_arrow = axis_arrow,
      title = panel_title,
      legend.position = legend.position,
      legend.direction = legend.direction,
      theme_use = theme_use,
      theme_args = theme_args
    )
  })
  names(plots) <- lineages

  if (isFALSE(combine)) {
    return(plots)
  }
  if (length(plots) == 1L) {
    return(plots[[1L]])
  }
  patchwork::wrap_plots(plots, nrow = nrow, ncol = ncol, byrow = byrow)
}

branch_stream_single_plot <- function(
  obs,
  group.by,
  pseudotime.by,
  labels = NULL,
  trunk_groups = NULL,
  branch_groups = NULL,
  branch_center = "auto",
  branch_steepness = 11,
  branch_power = 1.1,
  branch_amplitude = 0.28,
  n_branches = 2,
  n_grid = 800,
  bw = 0.15,
  pad = 0.035,
  count_power = 0.3,
  normalize_pseudotime = FALSE,
  scale_to = "auto",
  min_visible_width = 1e-4,
  xlim = NULL,
  ylim = NULL,
  xlabel = "Pseudotime",
  xticks = "auto",
  axis_y = "auto",
  label_positions = "auto",
  label_size = 4,
  label_fontface = "bold",
  label_fill = NULL,
  label_fill_alpha = 0.85,
  label_outline_width = 0.0018,
  label_outline_alpha = 0.8,
  ribbon_alpha = 1,
  palette = "Chinese",
  palcolor = NULL,
  axis_arrow = grid::arrow(length = grid::unit(0.08, "inches")),
  title = NULL,
  legend.position = "none",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list()
) {
  obs <- branch_stream_filter_obs(obs, group.by, pseudotime.by)
  obs <- branch_stream_prepare_pseudotime(
    obs,
    pseudotime.by = pseudotime.by,
    normalize_pseudotime = normalize_pseudotime
  )
  x_scale <- branch_stream_resolve_x_scale(
    obs[[pseudotime.by]],
    xlim = xlim,
    xticks = xticks,
    n_grid = n_grid
  )
  x <- x_scale$x
  xlim <- x_scale$xlim
  xticks <- x_scale$xticks
  branch_center <- branch_stream_resolve_branch_center(branch_center, x)
  branch_steepness <- branch_stream_resolve_branch_steepness(
    branch_steepness,
    x
  )
  groups <- obs[[group.by]]
  if (is.null(labels)) {
    labels <- branch_stream_group_order(groups)
  }
  labels <- as.character(labels)

  profiles <- branch_stream_kde_profiles(
    obs = obs,
    group.by = group.by,
    pseudotime.by = pseudotime.by,
    labels = labels,
    x = x,
    bw = bw,
    pad = pad,
    count_power = count_power
  )
  inferred <- branch_stream_infer_groups(
    x = x,
    profiles = profiles,
    labels = labels,
    trunk_groups = trunk_groups,
    branch_groups = branch_groups,
    branch_center = branch_center,
    branch_amplitude = branch_amplitude,
    n_branches = n_branches
  )

  branches <- list(list(
    center = rep(0, length(x)),
    layers = lapply(inferred$trunk_groups, function(group) {
      list(label = group, width = profiles[[group]])
    })
  ))
  amplitudes <- unique(unname(inferred$branch_groups))
  for (amp in amplitudes) {
    branch_labels <- names(inferred$branch_groups)[
      unname(inferred$branch_groups) == amp
    ]
    branches[[length(branches) + 1L]] <- list(
      center = branch_stream_centerline(
        x,
        amplitude = amp,
        center = branch_center,
        steepness = branch_steepness,
        power = branch_power
      ),
      layers = lapply(branch_labels, function(group) {
        list(label = group, width = profiles[[group]])
      })
    )
  }
  branches <- branch_stream_scale_widths(branches, scale_to)

  ribbon_df <- branch_stream_ribbon_data(
    x = x,
    branches = branches,
    min_visible_width = min_visible_width
  )
  bounds <- branch_stream_resolve_bounds(
    x,
    branches,
    ylim = ylim,
    axis_y = axis_y
  )
  ylim <- bounds$ylim
  axis_y <- bounds$axis_y
  label_df <- branch_stream_label_data(
    x = x,
    branches = branches,
    label_positions = label_positions,
    label_size = label_size
  )

  colors <- palette_colors(labels, palette = palette, palcolor = palcolor)
  theme_obj <- do.call(theme_use, theme_args)

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = ribbon_df,
      ggplot2::aes(
        x = .data$x,
        ymin = .data$ymin,
        ymax = .data$ymax,
        fill = .data$label
      ),
      linewidth = 0,
      alpha = ribbon_alpha,
      na.rm = TRUE
    ) +
    ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, clip = "off") +
    ggplot2::scale_x_continuous(breaks = xticks) +
    ggplot2::labs(title = title, x = xlabel, y = NULL, fill = NULL) +
    theme_obj +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(
        size = 12,
        margin = ggplot2::margin(t = 2)
      ),
      axis.text.x = ggplot2::element_text(
        size = 11,
        color = "black",
        margin = ggplot2::margin(t = 1)
      ),
      legend.position = legend.position,
      legend.direction = legend.direction
    )

  if (!is.null(axis_y)) {
    axis_x_pad <- diff(xlim) * 0.01
    p <- p +
      ggplot2::geom_segment(
        data = data.frame(
          x = xlim[1] + axis_x_pad,
          xend = xlim[2] - axis_x_pad,
          y = axis_y,
          yend = axis_y
        ),
        ggplot2::aes(
          x = .data$x,
          xend = .data$xend,
          y = .data$y,
          yend = .data$yend
        ),
        inherit.aes = FALSE,
        linewidth = 0.8,
        arrow = axis_arrow
      )
  }
  if (!is.null(label_df) && nrow(label_df) > 0L) {
    if (is.null(label_fill)) {
      outline_df <- branch_stream_label_outline_data(
        label_df,
        xlim = xlim,
        ylim = ylim,
        width = label_outline_width
      )
      p <- p +
        ggplot2::geom_text(
          data = outline_df,
          ggplot2::aes(
            x = .data$x,
            y = .data$y,
            label = .data$label,
            size = .data$size
          ),
          color = grDevices::adjustcolor(
            "white",
            alpha.f = label_outline_alpha
          ),
          fontface = label_fontface,
          show.legend = FALSE
        ) +
        ggplot2::geom_text(
          data = label_df,
          ggplot2::aes(
            x = .data$x,
            y = .data$y,
            label = .data$label,
            color = .data$label,
            size = .data$size
          ),
          fontface = label_fontface,
          show.legend = FALSE
        ) +
        ggplot2::scale_size_identity() +
        ggplot2::scale_color_manual(values = colors, guide = "none")
    } else {
      p <- p +
        ggplot2::geom_label(
          data = label_df,
          ggplot2::aes(
            x = .data$x,
            y = .data$y,
            label = .data$label,
            color = .data$label,
            size = .data$size
          ),
          fill = grDevices::adjustcolor(label_fill, alpha.f = label_fill_alpha),
          label.size = 0,
          fontface = label_fontface,
          show.legend = FALSE
        ) +
        ggplot2::scale_size_identity() +
        ggplot2::scale_color_manual(values = colors, guide = "none")
    }
  }
  p
}

branch_stream_label_outline_data <- function(
  label_df,
  xlim,
  ylim,
  width = 0.004
) {
  if (is.null(width) || !is.finite(width) || width <= 0) {
    return(label_df[0, , drop = FALSE])
  }
  x_offset <- diff(xlim) * width
  y_offset <- diff(ylim) * width
  offsets <- expand.grid(
    dx = c(-x_offset, 0, x_offset),
    dy = c(-y_offset, 0, y_offset)
  )
  offsets <- offsets[offsets$dx != 0 | offsets$dy != 0, , drop = FALSE]
  out <- do.call(
    rbind,
    lapply(seq_len(nrow(offsets)), function(i) {
      df <- label_df
      df$x <- df$x + offsets$dx[i]
      df$y <- df$y + offsets$dy[i]
      df
    })
  )
  rownames(out) <- NULL
  out
}

branch_stream_resolve_obs <- function(object, group.by) {
  if (inherits(object, "Seurat")) {
    obs <- object[[]]
  } else {
    obs <- as.data.frame(object, check.names = FALSE)
  }
  if (!group.by %in% colnames(obs)) {
    log_message(
      "Missing column in {.arg object}: {.val {group.by}}",
      message_type = "error"
    )
  }
  obs <- obs[!is.na(obs[[group.by]]), , drop = FALSE]
  if (nrow(obs) == 0L) {
    log_message(
      "No valid cells are available for {.fn BranchStreamPlot}",
      message_type = "error"
    )
  }
  obs
}

branch_stream_resolve_lineages <- function(
  obs,
  lineages = NULL,
  lineage.merge = c("min", "mean", "max")
) {
  lineage.merge <- match.arg(lineage.merge)
  detected <- branch_stream_detect_lineage_columns(obs)

  use_merge <- is.null(lineages) || identical(lineages, "merge")
  if (use_merge) {
    if (length(detected) == 0L) {
      log_message(
        "No lineage pseudotime columns were found. Provide {.arg lineages} or add columns such as {.val Lineage1} or {.val pseudotime}.",
        message_type = "error"
      )
    }
    if (length(detected) == 1L) {
      titles <- stats::setNames(NA_character_, detected)
      return(list(
        obs = obs,
        lineages = detected,
        titles = titles,
        source_lineages = detected
      ))
    }
    merged_col <- ".branch_stream_merged_pseudotime"
    obs[[merged_col]] <- branch_stream_merge_lineage_values(
      obs[detected],
      method = lineage.merge
    )
    titles <- stats::setNames(NA_character_, merged_col)
    return(list(
      obs = obs,
      lineages = merged_col,
      titles = titles,
      source_lineages = detected
    ))
  }

  if (identical(lineages, "all")) {
    lineages <- detected
  }
  if (length(lineages) == 0L) {
    log_message(
      "No lineage pseudotime columns were found. Provide {.arg lineages} or add columns such as {.val Lineage1} or {.val pseudotime}.",
      message_type = "error"
    )
  }
  missing <- setdiff(lineages, colnames(obs))
  if (length(missing) > 0L) {
    log_message(
      "Missing lineage columns in {.arg object}: {.val {missing}}",
      message_type = "error"
    )
  }
  lineages <- as.character(lineages)
  titles <- if (length(lineages) > 1L) lineages else NA_character_
  titles <- stats::setNames(titles, lineages)
  list(
    obs = obs,
    lineages = lineages,
    titles = titles,
    source_lineages = lineages
  )
}

branch_stream_detect_lineage_columns <- function(obs) {
  cols <- colnames(obs)
  detected <- grep(
    "(^|_)Lineage[0-9]+$|(^|_)pseudotime([0-9]+|_[0-9]+)?$",
    cols,
    value = TRUE,
    ignore.case = TRUE
  )
  if (length(detected) > 0L) {
    numeric_cols <- vapply(
      obs[detected],
      function(x) {
        suppressWarnings(any(is.finite(as.numeric(x))))
      },
      logical(1)
    )
    detected <- detected[numeric_cols]
  }
  detected
}

branch_stream_merge_lineage_values <- function(
  lineage_df,
  method = c("min", "mean", "max")
) {
  method <- match.arg(method)
  mat <- as.data.frame(lineage_df, check.names = FALSE)
  mat[] <- lapply(mat, function(x) suppressWarnings(as.numeric(x)))
  mat <- as.matrix(mat)
  apply(mat, 1, function(values) {
    values <- values[is.finite(values)]
    if (length(values) == 0L) {
      return(NA_real_)
    }
    switch(
      method,
      min = min(values),
      mean = mean(values),
      max = max(values)
    )
  })
}

branch_stream_filter_obs <- function(obs, group.by, pseudotime.by) {
  values <- suppressWarnings(as.numeric(obs[[pseudotime.by]]))
  obs <- obs[is.finite(values), , drop = FALSE]
  values <- values[is.finite(values)]
  if (nrow(obs) == 0L) {
    log_message(
      "No valid cells are available for lineage {.val {pseudotime.by}}",
      message_type = "error"
    )
  }
  obs[[pseudotime.by]] <- as.numeric(obs[[pseudotime.by]])
  obs
}

branch_stream_prepare_pseudotime <- function(
  obs,
  pseudotime.by,
  normalize_pseudotime = FALSE
) {
  values <- suppressWarnings(as.numeric(obs[[pseudotime.by]]))
  obs[[pseudotime.by]] <- values
  do_normalize <- if (identical(normalize_pseudotime, "auto")) {
    min(values, na.rm = TRUE) < 0 || max(values, na.rm = TRUE) > 1
  } else if (
    is.logical(normalize_pseudotime) &&
      length(normalize_pseudotime) == 1L
  ) {
    isTRUE(normalize_pseudotime)
  } else {
    log_message(
      "{.arg normalize_pseudotime} must be TRUE, FALSE, or {.val auto}",
      message_type = "error"
    )
  }
  if (do_normalize) {
    rng <- range(values, na.rm = TRUE)
    if (diff(rng) > 0) {
      obs[[pseudotime.by]] <- (values - rng[1]) / diff(rng)
    }
  }
  obs
}

branch_stream_resolve_x_scale <- function(values, xlim, xticks, n_grid) {
  rng <- range(values, na.rm = TRUE)
  if (!all(is.finite(rng))) {
    rng <- c(0, 1)
  }
  if (rng[1] == rng[2]) {
    offset <- max(abs(rng[1]) * 0.05, 0.5)
    rng <- rng + c(-offset, offset)
  }
  if (is.null(xlim)) {
    pad <- diff(rng) * 0.02
    xlim <- rng + c(-pad, pad)
  } else {
    xlim <- as.numeric(xlim)
    if (length(xlim) != 2L || any(!is.finite(xlim)) || xlim[1] >= xlim[2]) {
      log_message(
        "{.arg xlim} must be a numeric vector of length 2 with increasing finite values",
        message_type = "error"
      )
    }
  }
  if (identical(xticks, "auto")) {
    xticks <- pretty(xlim, n = 6)
    xticks <- xticks[xticks >= xlim[1] & xticks <= xlim[2]]
  }
  list(
    x = seq(xlim[1], xlim[2], length.out = n_grid),
    xlim = xlim,
    xticks = xticks
  )
}

branch_stream_resolve_branch_center <- function(branch_center, x) {
  if (identical(branch_center, "auto")) {
    rng <- range(x, na.rm = TRUE)
    return(rng[1] + diff(rng) * 0.56)
  }
  branch_center <- as.numeric(branch_center)
  if (length(branch_center) != 1L || !is.finite(branch_center)) {
    log_message(
      "{.arg branch_center} must be a single finite number or {.val auto}",
      message_type = "error"
    )
  }
  branch_center
}

branch_stream_resolve_branch_steepness <- function(branch_steepness, x) {
  branch_steepness <- as.numeric(branch_steepness)
  if (length(branch_steepness) != 1L || !is.finite(branch_steepness)) {
    log_message(
      "{.arg branch_steepness} must be a single finite number",
      message_type = "error"
    )
  }
  x_span <- diff(range(x, na.rm = TRUE))
  if (!is.finite(x_span) || x_span <= 0) {
    return(branch_steepness)
  }
  branch_steepness / x_span
}

branch_stream_group_order <- function(groups) {
  if (is.factor(groups)) {
    return(as.character(levels(groups)))
  }
  unique(as.character(groups))
}

branch_stream_sigmoid <- function(x, center, steepness) {
  1 / (1 + exp(-(x - center) * steepness))
}

branch_stream_tapered_density <- function(values, x, bw = 0.15, pad = 0.035) {
  values <- values[is.finite(values)]
  if (length(values) < 3L || length(unique(values)) < 2L) {
    return(rep(0, length(x)))
  }
  bw_use <- bw
  if (is.numeric(bw) && length(bw) == 1L && is.finite(bw) && bw > 0) {
    # Match scipy.stats.gaussian_kde(bw_method = bw) for a scalar bw_method.
    bw_use <- stats::sd(values) * bw
    if (!is.finite(bw_use) || bw_use <= 0) {
      bw_use <- bw
    }
  }
  density <- tryCatch(
    {
      z <- outer(x, values, "-") / bw_use
      rowMeans(stats::dnorm(z)) / bw_use
    },
    error = function(e) rep(0, length(x))
  )
  x_span <- diff(range(x, na.rm = TRUE))
  pad_use <- if (
    is.numeric(pad) &&
      length(pad) == 1L &&
      is.finite(pad) &&
      abs(pad) < 1 &&
      is.finite(x_span)
  ) {
    pad * x_span
  } else {
    pad
  }
  qs <- stats::quantile(values, c(0.05, 0.95), na.rm = TRUE, names = FALSE)
  start <- max(min(x), qs[1] - pad_use)
  end <- min(max(x), qs[2] + pad_use)
  gate_steepness <- if (is.finite(x_span) && x_span > 0) {
    40 / x_span
  } else {
    40
  }
  gate <- branch_stream_sigmoid(x, start, gate_steepness) *
    (1 - branch_stream_sigmoid(x, end, gate_steepness))
  density * gate
}

branch_stream_kde_profiles <- function(
  obs,
  group.by,
  pseudotime.by,
  labels,
  x,
  bw = 0.15,
  pad = 0.035,
  count_power = 0.3
) {
  groups <- as.character(obs[[group.by]])
  counts <- table(groups)
  size_factor <- as.numeric(counts)^count_power
  names(size_factor) <- names(counts)

  profiles <- lapply(labels, function(label) {
    values <- obs[[pseudotime.by]][groups == label]
    branch_stream_tapered_density(values, x, bw = bw, pad = pad) *
      (size_factor[[label]] %||% 1)
  })
  names(profiles) <- labels
  profiles
}

branch_stream_peak_positions <- function(x, profiles, labels) {
  peaks <- rep(NA_real_, length(labels))
  names(peaks) <- labels
  for (label in labels) {
    profile <- profiles[[label]]
    if (
      length(profile) > 0L &&
        any(is.finite(profile)) &&
        max(profile, na.rm = TRUE) > 0
    ) {
      peaks[[label]] <- x[which.max(profile)]
    }
  }
  peaks
}

branch_stream_infer_groups <- function(
  x,
  profiles,
  labels,
  trunk_groups = NULL,
  branch_groups = NULL,
  branch_center = 0.56,
  branch_amplitude = 0.28,
  n_branches = 2
) {
  if (!is.null(trunk_groups)) {
    trunk <- as.character(trunk_groups)
  } else {
    peaks <- branch_stream_peak_positions(x, profiles, labels)
    trunk <- labels[is.finite(peaks[labels]) & peaks[labels] <= branch_center]
    if (length(trunk) == length(labels)) {
      trunk <- labels[-length(labels)]
    }
    if (length(trunk) == 0L) {
      trunk <- labels[seq_len(max(1, ceiling(length(labels) * 0.5)))]
    }
  }

  if (!is.null(branch_groups)) {
    branches <- unlist(branch_groups)
    branches <- stats::setNames(as.numeric(branches), names(branches))
  } else {
    branch_labels <- setdiff(labels, trunk)
    if (length(branch_labels) == 0L) {
      branch_labels <- labels[length(labels)]
      trunk <- setdiff(trunk, branch_labels)
    }
    n_branches <- max(1L, as.integer(n_branches))
    amplitudes <- if (n_branches == 1L) {
      branch_amplitude
    } else {
      seq(branch_amplitude, -branch_amplitude, length.out = n_branches)
    }
    branches <- stats::setNames(
      amplitudes[((seq_along(branch_labels) - 1L) %% length(amplitudes)) + 1L],
      branch_labels
    )
  }
  list(trunk_groups = trunk, branch_groups = branches)
}

branch_stream_centerline <- function(
  x,
  amplitude,
  center,
  steepness,
  power = 1.1,
  baseline = 0
) {
  baseline + amplitude * (branch_stream_sigmoid(x, center, steepness)^power)
}

branch_stream_scale_widths <- function(branches, scale_to = "auto") {
  if (is.null(scale_to)) {
    return(branches)
  }
  max_total <- max(
    vapply(
      branches,
      function(branch) {
        total <- Reduce(
          `+`,
          lapply(branch$layers, `[[`, "width"),
          init = rep(0, length(branch$center))
        )
        max(total, na.rm = TRUE)
      },
      numeric(1)
    ),
    na.rm = TRUE
  )
  if (!is.finite(max_total) || max_total <= 0) {
    return(branches)
  }
  if (identical(scale_to, "auto")) {
    scale_to <- 0.32
  }
  scale <- as.numeric(scale_to) / max_total
  lapply(branches, function(branch) {
    branch$layers <- lapply(branch$layers, function(layer) {
      layer$width <- layer$width * scale
      layer
    })
    branch
  })
}

branch_stream_ribbon_data <- function(x, branches, min_visible_width = 1e-4) {
  out <- list()
  for (branch_id in seq_along(branches)) {
    branch <- branches[[branch_id]]
    total <- Reduce(
      `+`,
      lapply(branch$layers, `[[`, "width"),
      init = rep(0, length(x))
    )
    cumulative <- rep(0, length(x))
    for (layer in branch$layers) {
      width <- layer$width
      ymin <- branch$center - total / 2 + cumulative
      ymax <- ymin + width
      keep <- width > min_visible_width
      ymin[!keep] <- NA_real_
      ymax[!keep] <- NA_real_
      out[[length(out) + 1L]] <- data.frame(
        x = x,
        ymin = ymin,
        ymax = ymax,
        label = layer$label,
        branch = branch_id,
        stringsAsFactors = FALSE
      )
      cumulative <- cumulative + width
    }
  }
  do.call(rbind, out)
}

branch_stream_bounds <- function(x, branches) {
  lower <- upper <- numeric()
  for (branch in branches) {
    total <- Reduce(
      `+`,
      lapply(branch$layers, `[[`, "width"),
      init = rep(0, length(x))
    )
    lower <- c(lower, branch$center - total / 2)
    upper <- c(upper, branch$center + total / 2)
  }
  if (length(lower) == 0L || length(upper) == 0L) {
    return(c(-0.5, 0.5))
  }
  bounds <- c(min(lower, na.rm = TRUE), max(upper, na.rm = TRUE))
  if (!all(is.finite(bounds)) || bounds[1] == bounds[2]) {
    bounds <- c(-0.5, 0.5)
  }
  bounds
}

branch_stream_resolve_bounds <- function(
  x,
  branches,
  ylim = NULL,
  axis_y = "auto"
) {
  bounds <- branch_stream_bounds(x, branches)
  span <- max(diff(bounds), 0.2)
  resolved_axis_y <- axis_y
  if (identical(axis_y, "auto")) {
    resolved_axis_y <- bounds[1] - span * 0.08
  }
  if (is.null(ylim)) {
    ylim <- c(bounds[1] - span * 0.10, bounds[2] + span * 0.18)
    if (!is.null(resolved_axis_y)) {
      ylim[1] <- min(ylim[1], resolved_axis_y - span * 0.035)
    }
  }
  list(ylim = ylim, axis_y = resolved_axis_y)
}

branch_stream_label_data <- function(
  x,
  branches,
  label_positions = "auto",
  label_size = 4
) {
  if (is.null(label_positions)) {
    return(NULL)
  }
  if (!identical(label_positions, "auto")) {
    if (
      is.list(label_positions) &&
        !is.data.frame(label_positions) &&
        !is.null(names(label_positions)) &&
        !all(c("label", "x", "y") %in% names(label_positions))
    ) {
      df <- do.call(
        rbind,
        lapply(names(label_positions), function(label) {
          pos <- unlist(label_positions[[label]])
          size <- as.numeric(pos[3])
          if (!is.finite(size)) {
            size <- label_size
          }
          data.frame(
            label = label,
            x = as.numeric(pos[1]),
            y = as.numeric(pos[2]),
            size = size,
            stringsAsFactors = FALSE
          )
        })
      )
    } else {
      df <- as.data.frame(label_positions, check.names = FALSE)
    }
    if (!all(c("label", "x", "y") %in% colnames(df))) {
      log_message(
        "{.arg label_positions} must contain {.val label}, {.val x}, and {.val y}",
        message_type = "error"
      )
    }
    df$size <- df$size %||% label_size
    return(df)
  }

  out <- list()
  for (branch in branches) {
    total <- Reduce(
      `+`,
      lapply(branch$layers, `[[`, "width"),
      init = rep(0, length(x))
    )
    cumulative <- rep(0, length(x))
    for (layer in branch$layers) {
      width <- layer$width
      if (
        length(width) > 0L &&
          any(is.finite(width)) &&
          max(width, na.rm = TRUE) > 0
      ) {
        idx <- which.max(width)
        out[[length(out) + 1L]] <- data.frame(
          label = layer$label,
          x = x[idx],
          y = branch$center[idx] -
            total[idx] / 2 +
            cumulative[idx] +
            width[idx] / 2,
          size = label_size,
          stringsAsFactors = FALSE
        )
      }
      cumulative <- cumulative + width
    }
  }
  do.call(rbind, out)
}
