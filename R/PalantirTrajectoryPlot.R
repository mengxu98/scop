#' @title Plot Palantir trajectories
#'
#' @description
#' Plot branch-aware Palantir trajectories on a two-dimensional embedding.
#'
#' @md
#' @inheritParams CellDimPlot
#' @param pseudotime_key Name of the metadata column containing Palantir
#' pseudotime.
#' @param branch_cols Metadata columns containing Palantir branch probabilities.
#' If `NULL`, columns ending with `"_diff_potential"` are used, excluding
#' `pseudotime_key` and `diff_potential_key`.
#' @param diff_potential_key Name of the Palantir entropy/differentiation
#' potential column to exclude from branch auto-detection.
#' @param pseudotime_interval Numeric vector of length 2 specifying the
#' pseudotime range to plot.
#' @param branch_min_prob Minimum branch probability used to select cells for a
#' branch trajectory.
#' @param n_bins Number of pseudotime bins used to summarize each trajectory.
#' @param min_cells_per_bin Minimum number of cells required in a bin.
#' @param smooth Whether to smooth the trajectory with [stats::loess].
#' @param trajectory_method Method used to fit trajectory coordinates along
#' pseudotime. `"loess"` uses a fully R-native smoother over a Palantir-style
#' pseudotime grid; `"bin"` uses binned median coordinates.
#' @param smoothness Smoothing multiplier for the R-native loess span. Higher
#' values yield smoother curves.
#' @param span Base span used for loess smoothing.
#' @param n_path_points Number of pseudotime points used to draw each smoothed
#' trajectory.
#' @param cell_color Cell coloring mode. Use `"pseudotime"` for Palantir
#' pseudotime, `"branch_selection"` for the branch with highest probability,
#' `"none"` to hide cell coloring, or any metadata column.
#' @param pt.size Point size for cells.
#' @param pt.alpha Point alpha for cells.
#' @param trajectory_palette Color palette for trajectories.
#' @param trajectory_palcolor Custom colors for trajectories.
#' @param trajectory_linewidth Line width of trajectories.
#' @param trajectory_bg Color for the trajectory background stroke.
#' @param trajectory_bg_stroke Width added to the trajectory background stroke.
#' @param trajectory_arrow Arrow used for trajectories. See [grid::arrow].
#' @param return_layer Logical. If `TRUE`, returns ggplot2 layers instead of a
#' complete plot.
#'
#' @seealso [RunPalantir], [FeatureDimPlot], [CellDimPlot]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunPalantir(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   linear_reduction = "PCA",
#'   nonlinear_reduction = "UMAP",
#'   early_group = "Ductal",
#'   terminal_groups = c("Alpha", "Beta", "Delta", "Epsilon")
#' )
#' PalantirTrajectoryPlot(
#'   pancreas_sub,
#'   reduction = "UMAP",
#'   pseudotime_interval = c(0, 0.9)
#' )
#' PalantirTrajectoryPlot(
#'   pancreas_sub,
#'   reduction = "UMAP",
#'   cell_color = "branch_selection",
#'   pseudotime_interval = c(0, 0.9)
#' )
#' }
PalantirTrajectoryPlot <- function(
  srt,
  reduction = NULL,
  dims = c(1, 2),
  cells = NULL,
  pseudotime_key = "palantir_pseudotime",
  branch_cols = NULL,
  diff_potential_key = "palantir_diff_potential",
  pseudotime_interval = c(0, 1),
  branch_min_prob = 0.05,
  n_bins = 60,
  min_cells_per_bin = 3,
  smooth = TRUE,
  trajectory_method = c("loess", "bin"),
  smoothness = 1,
  span = 0.75,
  n_path_points = 200,
  cell_color = "pseudotime",
  pt.size = 0.5,
  pt.alpha = 0.8,
  palette = "Dark2",
  palcolor = NULL,
  trajectory_palette = "Dark2",
  trajectory_palcolor = NULL,
  trajectory_linewidth = 1.2,
  trajectory_bg = "black",
  trajectory_bg_stroke = 0.7,
  trajectory_arrow = grid::arrow(
    length = grid::unit(0.12, "inches"),
    type = "closed"
  ),
  aspect.ratio = 1,
  title = "Palantir",
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  return_layer = FALSE,
  seed = 11,
  verbose = TRUE
) {
  set.seed(seed)
  trajectory_method <- match.arg(trajectory_method)

  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  if (!reduction %in% names(srt@reductions)) {
    log_message(
      "{.val {reduction}} is not in the srt reduction names",
      message_type = "error"
    )
  }
  if (!pseudotime_key %in% colnames(srt@meta.data)) {
    log_message(
      "Cannot find the Palantir pseudotime column {.val {pseudotime_key}}",
      message_type = "error"
    )
  }

  if (is.null(branch_cols)) {
    branch_cols <- grep(
      "_diff_potential$",
      colnames(srt@meta.data),
      value = TRUE
    )
    branch_cols <- setdiff(branch_cols, c(pseudotime_key, diff_potential_key))
  }
  if (!length(branch_cols)) {
    log_message(
      "No Palantir branch probability columns were found. Provide {.arg branch_cols}.",
      message_type = "error"
    )
  }
  missing_cols <- setdiff(branch_cols, colnames(srt@meta.data))
  if (length(missing_cols)) {
    log_message(
      "Cannot find branch columns: {.val {missing_cols}}",
      message_type = "error"
    )
  }
  if (!is.numeric(pseudotime_interval) || length(pseudotime_interval) != 2L) {
    log_message(
      "{.arg pseudotime_interval} must be a numeric vector of length 2",
      message_type = "error"
    )
  }
  pseudotime_interval <- sort(pseudotime_interval)

  reduction_key <- srt@reductions[[reduction]]@key
  emb <- srt@reductions[[reduction]]@cell.embeddings[, dims, drop = FALSE]
  colnames(emb) <- paste0(reduction_key, dims)
  rownames(emb) <- rownames(emb) %||% colnames(srt@assays[[1]])

  meta_cols <- unique(c(pseudotime_key, branch_cols))
  if (
    !cell_color %in% c("pseudotime", "branch_selection", "none") &&
      cell_color %in% colnames(srt@meta.data)
  ) {
    meta_cols <- unique(c(meta_cols, cell_color))
  }
  dat <- cbind.data.frame(
    emb,
    srt@meta.data[rownames(emb), meta_cols, drop = FALSE]
  )
  dat[["cell"]] <- rownames(dat)
  if (!is.null(cells)) {
    dat <- dat[intersect(rownames(dat), cells), , drop = FALSE]
  }

  dat[[pseudotime_key]] <- suppressWarnings(as.numeric(dat[[pseudotime_key]]))
  for (branch in branch_cols) {
    dat[[branch]] <- suppressWarnings(as.numeric(dat[[branch]]))
  }
  valid <- !is.na(dat[[pseudotime_key]]) &
    dat[[pseudotime_key]] >= pseudotime_interval[1] &
    dat[[pseudotime_key]] <= pseudotime_interval[2]
  dat <- dat[valid, , drop = FALSE]
  if (!nrow(dat)) {
    log_message(
      "No cells remain after applying {.arg pseudotime_interval}",
      message_type = "error"
    )
  }

  axes <- paste0(reduction_key, dims)
  xlab <- xlab %||% axes[1]
  ylab <- ylab %||% axes[2]
  if (identical(theme_use, "theme_blank")) {
    theme_args[["xlab"]] <- xlab
    theme_args[["ylab"]] <- ylab
  }

  dat[["branch_selection"]] <- palantir_branch_selection(dat, branch_cols)
  path_df <- palantir_trajectory_paths(
    dat = dat,
    axes = axes,
    pseudotime_key = pseudotime_key,
    pseudotime_interval = pseudotime_interval,
    branch_cols = branch_cols,
    branch_min_prob = branch_min_prob,
    n_bins = n_bins,
    min_cells_per_bin = min_cells_per_bin,
    smooth = smooth,
    smoothness = smoothness,
    trajectory_method = trajectory_method,
    span = span,
    n_path_points = n_path_points,
    verbose = verbose
  )
  branch_levels <- unique(path_df[["branch"]])
  branch_labels <- palantir_branch_labels(branch_levels)
  path_df[["branch"]] <- factor(path_df[["branch"]], levels = branch_levels)

  trajectory_colors <- palette_colors(
    branch_levels,
    palette = trajectory_palette,
    palcolor = trajectory_palcolor
  )
  names(trajectory_colors) <- branch_levels

  cell_layer <- palantir_cell_layer(
    dat = dat,
    axes = axes,
    pseudotime_key = pseudotime_key,
    branch_cols = branch_cols,
    cell_color = cell_color,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    palette = palette,
    palcolor = palcolor
  )

  trajectory_layer <- list(
    ggnewscale::new_scale_color(),
    geom_path(
      data = path_df,
      mapping = aes(
        x = .data[["Axis_1"]],
        y = .data[["Axis_2"]],
        group = .data[["branch"]]
      ),
      color = trajectory_bg,
      linewidth = trajectory_linewidth + trajectory_bg_stroke,
      arrow = trajectory_arrow,
      lineend = "round",
      linejoin = "round",
      show.legend = FALSE,
      inherit.aes = FALSE
    ),
    geom_path(
      data = path_df,
      mapping = aes(
        x = .data[["Axis_1"]],
        y = .data[["Axis_2"]],
        color = .data[["branch"]],
        group = .data[["branch"]]
      ),
      linewidth = trajectory_linewidth,
      lineend = "round",
      linejoin = "round",
      inherit.aes = FALSE
    ),
    geom_path(
      data = path_df,
      mapping = aes(
        x = .data[["Axis_1"]],
        y = .data[["Axis_2"]],
        color = .data[["branch"]],
        group = .data[["branch"]]
      ),
      linewidth = trajectory_linewidth,
      arrow = trajectory_arrow,
      lineend = "round",
      linejoin = "round",
      show.legend = FALSE,
      inherit.aes = FALSE
    ),
    scale_color_manual(
      name = "Palantir branch",
      values = trajectory_colors,
      labels = branch_labels,
      guide = guide_legend(
        title.hjust = 0,
        order = 2,
        override.aes = list(linewidth = trajectory_linewidth + 0.8)
      )
    )
  )

  lab_layer <- list(labs(
    title = title,
    subtitle = subtitle,
    x = xlab,
    y = ylab
  ))
  theme_layer <- list(
    do.call(theme_use, theme_args) +
      theme(
        aspect.ratio = aspect.ratio,
        legend.position = legend.position,
        legend.direction = legend.direction
      )
  )

  if (isTRUE(return_layer)) {
    return(list(
      cell_layer = cell_layer,
      trajectory_layer = trajectory_layer,
      lab_layer = lab_layer,
      theme_layer = theme_layer
    ))
  }

  ggplot() +
    cell_layer +
    trajectory_layer +
    lab_layer +
    theme_layer
}

palantir_branch_selection <- function(dat, branch_cols) {
  probs <- as.matrix(dat[, branch_cols, drop = FALSE])
  probs[is.na(probs)] <- -Inf
  selected <- branch_cols[max.col(probs, ties.method = "first")]
  has_prob <- rowSums(is.finite(probs)) > 0L
  selected[!has_prob] <- NA_character_
  factor(selected, levels = branch_cols)
}

palantir_trajectory_paths <- function(
  dat,
  axes,
  pseudotime_key,
  pseudotime_interval,
  branch_cols,
  branch_min_prob,
  n_bins,
  min_cells_per_bin,
  smooth,
  smoothness,
  trajectory_method,
  span,
  n_path_points,
  verbose = TRUE
) {
  paths <- lapply(branch_cols, function(branch) {
    branch_prob <- dat[[branch]]
    branch_dat <- dat[
      !is.na(branch_prob) & branch_prob >= branch_min_prob,
      ,
      drop = FALSE
    ]
    if (nrow(branch_dat) < min_cells_per_bin * 2L) {
      branch_dat <- dat[
        !is.na(branch_prob) & dat[["branch_selection"]] == branch,
        ,
        drop = FALSE
      ]
    }
    if (nrow(branch_dat) < min_cells_per_bin * 2L) {
      log_message(
        "Skipping Palantir branch {.val {branch}} because too few cells are available.",
        message_type = "warning",
        verbose = verbose
      )
      return(NULL)
    }
    branch_dat <- branch_dat[order(branch_dat[[pseudotime_key]]), , drop = FALSE]

    if (
      isTRUE(smooth) &&
        trajectory_method == "loess" &&
        nrow(branch_dat) >= max(10L, min_cells_per_bin * 2L)
    ) {
      path <- palantir_smooth_branch_path(
        branch_dat = branch_dat,
        axes = axes,
        pseudotime_key = pseudotime_key,
        branch = branch,
        pseudotime_interval = pseudotime_interval,
        n_path_points = n_path_points,
        span = min(1, max(0.05, span * smoothness))
      )
      if (!is.null(path) && nrow(path) >= 2L) {
        return(path)
      }
    }

    breaks <- seq(
      min(branch_dat[[pseudotime_key]], na.rm = TRUE),
      max(branch_dat[[pseudotime_key]], na.rm = TRUE),
      length.out = n_bins + 1L
    )
    branch_dat[[".bin"]] <- cut(
      branch_dat[[pseudotime_key]],
      breaks = unique(breaks),
      include.lowest = TRUE,
      labels = FALSE
    )
    binned <- stats::aggregate(
      branch_dat[, c(axes, pseudotime_key), drop = FALSE],
      by = list(.bin = branch_dat[[".bin"]]),
      FUN = stats::median,
      na.rm = TRUE
    )
    counts <- stats::aggregate(
      branch_dat[["cell"]],
      by = list(.bin = branch_dat[[".bin"]]),
      FUN = length
    )
    binned <- binned[counts[["x"]] >= min_cells_per_bin, , drop = FALSE]
    binned <- binned[order(binned[[pseudotime_key]]), , drop = FALSE]
    if (nrow(binned) < 2L) {
      return(NULL)
    }
    path <- data.frame(
      pseudotime = binned[[pseudotime_key]],
      Axis_1 = binned[[axes[1]]],
      Axis_2 = binned[[axes[2]]],
      branch = branch,
      stringsAsFactors = FALSE
    )
    unique(stats::na.omit(path))
  })
  paths <- do.call(rbind, paths)
  if (is.null(paths) || !nrow(paths)) {
    log_message(
      "No Palantir trajectories could be fitted. Try lowering {.arg branch_min_prob}.",
      message_type = "error"
    )
  }
  paths
}

palantir_smooth_branch_path <- function(
  branch_dat,
  axes,
  pseudotime_key,
  branch,
  pseudotime_interval,
  n_path_points,
  span
) {
  fit_dat <- data.frame(
    pseudotime = branch_dat[[pseudotime_key]],
    Axis_1 = branch_dat[[axes[1]]],
    Axis_2 = branch_dat[[axes[2]]],
    weight = branch_dat[[branch]],
    stringsAsFactors = FALSE
  )
  fit_dat <- fit_dat[stats::complete.cases(fit_dat), , drop = FALSE]
  if (nrow(fit_dat) < 10L || length(unique(fit_dat[["pseudotime"]])) < 5L) {
    return(NULL)
  }
  fit_dat[["weight"]] <- pmax(fit_dat[["weight"]], .Machine$double.eps)
  path_interval <- pseudotime_interval %||% c(
    min(fit_dat[["pseudotime"]], na.rm = TRUE),
    max(fit_dat[["pseudotime"]], na.rm = TRUE)
  )
  grid <- data.frame(
    pseudotime = seq(path_interval[1], path_interval[2], length.out = max(2L, as.integer(n_path_points)))
  )
  path <- tryCatch(
    {
      fit_x <- stats::loess(
        Axis_1 ~ pseudotime,
        data = fit_dat,
        weights = fit_dat[["weight"]],
        span = span,
        degree = 2,
        control = stats::loess.control(surface = "direct")
      )
      fit_y <- stats::loess(
        Axis_2 ~ pseudotime,
        data = fit_dat,
        weights = fit_dat[["weight"]],
        span = span,
        degree = 2,
        control = stats::loess.control(surface = "direct")
      )
      data.frame(
        pseudotime = grid[["pseudotime"]],
        Axis_1 = stats::predict(fit_x, newdata = grid),
        Axis_2 = stats::predict(fit_y, newdata = grid),
        branch = branch,
        stringsAsFactors = FALSE
      )
    },
    error = function(e) NULL
  )
  if (is.null(path)) {
    return(NULL)
  }
  unique(stats::na.omit(path))
}

palantir_branch_labels <- function(branches) {
  labels <- sub("_diff_potential$", "", branches)
  labels <- sub("_branch_mask$", "", labels)
  labels <- gsub("\\.", " ", labels)
  paste("Branch", labels)
}

palantir_cell_layer <- function(
  dat,
  axes,
  pseudotime_key,
  branch_cols,
  cell_color,
  pt.size,
  pt.alpha,
  palette,
  palcolor
) {
  if (identical(cell_color, "none")) {
    return(list())
  }
  if (identical(cell_color, "pseudotime")) {
    dat[[".cell_color"]] <- dat[[pseudotime_key]]
    return(list(
      geom_point(
        data = dat,
        mapping = aes(x = .data[[axes[1]]], y = .data[[axes[2]]], color = .data[[".cell_color"]]),
        size = pt.size,
        alpha = pt.alpha,
        inherit.aes = FALSE
      ),
      scale_color_viridis_c(name = pseudotime_key)
    ))
  }
  if (identical(cell_color, "branch_selection")) {
    dat[[".cell_color"]] <- dat[["branch_selection"]]
    colors <- palette_colors(branch_cols, palette = palette, palcolor = palcolor)
    names(colors) <- branch_cols
    return(list(
      geom_point(
        data = dat,
        mapping = aes(x = .data[[axes[1]]], y = .data[[axes[2]]], color = .data[[".cell_color"]]),
        size = pt.size,
        alpha = pt.alpha,
        inherit.aes = FALSE
      ),
      scale_color_manual(
        name = "branch_selection",
        values = colors,
        labels = palantir_branch_labels(branch_cols),
        na.value = "grey80"
      )
    ))
  }
  if (!cell_color %in% colnames(dat)) {
    log_message(
      "{.arg cell_color} must be {.val pseudotime}, {.val branch_selection}, {.val none}, or a metadata column in the plotted data.",
      message_type = "error"
    )
  }
  dat[[".cell_color"]] <- dat[[cell_color]]
  if (is.numeric(dat[[".cell_color"]])) {
    list(
      geom_point(
        data = dat,
        mapping = aes(x = .data[[axes[1]]], y = .data[[axes[2]]], color = .data[[".cell_color"]]),
        size = pt.size,
        alpha = pt.alpha,
        inherit.aes = FALSE
      ),
      scale_color_viridis_c(name = cell_color)
    )
  } else {
    list(
      geom_point(
        data = dat,
        mapping = aes(x = .data[[axes[1]]], y = .data[[axes[2]]], color = .data[[".cell_color"]]),
        size = pt.size,
        alpha = pt.alpha,
        inherit.aes = FALSE
      ),
      scale_color_manual(
        name = cell_color,
        values = palette_colors(dat[[".cell_color"]], palette = palette, palcolor = palcolor),
        na.value = "grey80"
      )
    )
  }
}
