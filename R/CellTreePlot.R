#' @title Plot a circular cell annotation tree
#'
#' @description
#' Show cell proportions and optional numeric rings. Requires `ggforce`.
#'
#' @md
#' @param object A Seurat object or a data frame containing cell annotations.
#' @param levels One or more unique annotation column names, broadest first.
#'   Labels must be nonmissing and nonblank. The names `.count` and
#'   `.proportion` are reserved.
#' @param count Column of nonnegative whole counts with a finite positive
#'   total. With `NULL`, each row counts as one cell.
#' @param tracks Data frame keyed by unique full annotation paths in `object`,
#'   with one numeric column per ring (inner to outer). Rings share a
#'   zero-centered scale. Missing values and omitted paths appear gray.
#' @param palette,palcolor Palette (default `"Dark2"`) and optional custom
#'   colors. Named colors override matching broad groups; other groups use
#'   the palette. Unnamed colors follow group order.
#' @param label.size Positive, finite leaf label size (default `3`).
#' @param show_proportion Logical. Label each leaf with its share of the total
#'   count (default `TRUE`).
#' @param title Plot title (default `NULL`, no title).
#' @param return_data Logical. Include the leaf and track data (default `FALSE`).
#' @param show_internal_labels Logical. Label intermediate levels (default `FALSE`).
#'
#' @return A `ggplot` object that accepts layers with `+`. With
#'   `return_data = TRUE`, a list with `plot`, `leaves` (one row per annotation
#'   path with `.count` and `.proportion`), and `tracks` (aligned to the leaves,
#'   or `NULL`).
#'
#' @references
#' \url{https://doi.org/10.1126/science.adt3130}
#'
#' \url{https://mp.weixin.qq.com/s/F7cgHFmpEAhwhrzAN0NvUA}
#'
#' @seealso [ClusterTreePlot()], [CellStatPlot()]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' hierarchy <- c("CellType", "SubCellType")
#' CellTreePlot(pancreas_sub, hierarchy)
#'
#' meta <- pancreas_sub[[]]
#' scores <- stats::aggregate(meta[c("S_score", "G2M_score")], meta[hierarchy], mean)
#' CellTreePlot(pancreas_sub, hierarchy, tracks = scores)
#'
#' data(pbmcmultiome_sub)
#' pbmc <- pbmcmultiome_sub[[]]
#' groups <- c(
#'   B = "B cells", "CD4 T" = "T cells", "CD8 T" = "T cells",
#'   NK = "NK cells", DC = "Myeloid", Mono = "Myeloid"
#' )
#' pbmc$Group <- unname(groups[as.character(pbmc$CellType)])
#' hierarchy <- c("Group", "CellType", "SubCellType")
#'
#' metrics <- c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC")
#' tracks <- stats::aggregate(pbmc[metrics], pbmc[hierarchy], stats::median)
#' subset_medians <- vapply(pbmc[metrics], stats::median, numeric(1))
#' tracks[metrics] <- Map(
#'   function(x, center) log2(x / center),
#'   tracks[metrics], subset_medians
#' )
#' names(tracks)[match(metrics, names(tracks))] <- c(
#'   "RNA counts", "RNA genes", "ATAC counts", "ATAC peaks"
#' )
#'
#' CellTreePlot(pbmc, hierarchy, tracks = tracks, label.size = 2.25)
#'
#' genes <- c("MS4A1", "CD3D", "NKG7", "LST1")
#' rna <- SeuratObject::LayerData(
#'   pbmcmultiome_sub,
#'   assay = "RNA", layer = "counts", features = genes
#' )
#' detected <- as.data.frame(t(as.matrix(rna[genes, rownames(pbmc)] > 0)))
#' gene_tracks <- stats::aggregate(detected, pbmc[hierarchy], mean)
#' subset_detection <- colMeans(detected)
#' gene_tracks[genes] <- Map(
#'   function(x, center) 100 * (x - center),
#'   gene_tracks[genes], subset_detection
#' )
#' CellTreePlot(pbmc, hierarchy, tracks = gene_tracks, label.size = 2.25)
CellTreePlot <- function(
  object,
  levels,
  count = NULL,
  tracks = NULL,
  palette = "Dark2",
  palcolor = NULL,
  label.size = 3,
  show_proportion = TRUE,
  title = NULL,
  return_data = FALSE,
  show_internal_labels = FALSE
) {
  if (!is.character(levels) || !is.null(dim(levels)) || !length(levels) ||
    anyNA(levels) || any(!nzchar(trimws(levels))) || anyDuplicated(levels)) {
    stop("`levels` must be a nonempty vector of unique annotation column names.", call. = FALSE)
  }
  levels <- unname(levels)
  if (any(levels %in% c(".count", ".proportion"))) {
    stop("`.count` and `.proportion` are reserved and cannot be used in `levels`.", call. = FALSE)
  }
  if (!is.numeric(label.size) || is.complex(label.size) || length(label.size) != 1L ||
    !is.finite(label.size) || label.size <= 0) {
    stop("`label.size` must be one positive, finite number.", call. = FALSE)
  }
  if (!is.logical(show_proportion) || length(show_proportion) != 1L || is.na(show_proportion)) {
    stop("`show_proportion` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(return_data) || length(return_data) != 1L || is.na(return_data)) {
    stop("`return_data` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(show_internal_labels) || length(show_internal_labels) != 1L ||
    !is.null(dim(show_internal_labels)) || is.na(show_internal_labels)) {
    stop("`show_internal_labels` must be TRUE or FALSE.", call. = FALSE)
  }

  metadata <- if (inherits(object, "Seurat")) object[[]] else object
  if (!is.data.frame(metadata) || !nrow(metadata)) {
    stop("`object` must be a Seurat object or a nonempty data frame.", call. = FALSE)
  }
  metadata <- as.data.frame(metadata)
  paths <- cell_tree_paths(metadata, levels, "object")
  if (is.null(count)) {
    weights <- rep(1, nrow(metadata))
  } else {
    if (!is.character(count) || length(count) != 1L || is.na(count) ||
      !nzchar(trimws(count)) || sum(names(metadata) %in% count) != 1L) {
      stop("`count` must name one column in `object`.", call. = FALSE)
    }
    weights <- metadata[[count]]
    if (!cell_tree_numeric(weights) || any(!is.finite(weights)) ||
      any(weights < 0 | weights != trunc(weights))) {
      stop("The `count` column must contain finite nonnegative whole numbers.", call. = FALSE)
    }
    weights <- as.numeric(weights)
  }
  total <- sum(weights)
  if (!is.finite(total) || total <= 0) {
    stop("The total cell count must be finite and positive.", call. = FALSE)
  }

  path_order <- do.call(order, c(unname(as.list(paths)), list(method = "radix")))
  paths <- paths[path_order, , drop = FALSE]
  prefix_starts <- cell_tree_prefix_starts(paths)
  leaf_starts <- prefix_starts[[length(levels)]]
  leaves <- paths[leaf_starts, , drop = FALSE]
  leaf_sizes <- diff(c(leaf_starts, nrow(paths) + 1L))
  leaves$.count <- as.numeric(rowsum(
    weights[path_order], rep.int(seq_along(leaf_starts), leaf_sizes),
    reorder = FALSE
  ))
  leaves$.proportion <- leaves$.count / total
  rownames(leaves) <- NULL
  aligned <- cell_tree_align_tracks(tracks, leaves[levels])
  check_r("ggforce", verbose = FALSE)

  plot <- cell_tree_draw(
    leaves, lapply(prefix_starts, match, table = leaf_starts), levels, aligned,
    palette, palcolor, label.size, show_proportion, title, show_internal_labels
  )
  if (return_data) {
    return(list(plot = plot, leaves = leaves, tracks = aligned))
  }
  plot
}

cell_tree_paths <- function(data, levels, arg) {
  paths <- lapply(levels, function(level) {
    if (sum(names(data) %in% level) != 1L) {
      stop(sprintf("`%s` must contain exactly one `%s` column.", arg, level), call. = FALSE)
    }
    value <- data[[level]]
    if (!is.atomic(value) || !is.null(dim(value))) {
      stop(sprintf("`%s` in `%s` must be a vector of annotation labels.", level, arg), call. = FALSE)
    }
    labels <- enc2utf8(as.character(value))
    if (anyNA(value) || anyNA(labels) || any(!nzchar(trimws(labels, whitespace = "[\\h\\v]")))) {
      stop(sprintf("`%s` in `%s` contains missing or blank annotations.", level, arg), call. = FALSE)
    }
    labels
  })
  names(paths) <- levels
  as.data.frame(paths, check.names = FALSE, stringsAsFactors = FALSE)
}

cell_tree_numeric <- function(x) {
  is.numeric(x) && !is.complex(x) && is.null(dim(x))
}

cell_tree_prefix_starts <- function(paths) {
  changed <- rep(FALSE, nrow(paths) - 1L)
  starts <- vector("list", ncol(paths))
  for (i in seq_along(paths)) {
    labels <- paths[[i]]
    changed <- changed | (labels[-1L] != labels[-length(labels)])
    starts[[i]] <- which(c(TRUE, changed))
  }
  starts
}

cell_tree_align_tracks <- function(tracks, paths) {
  if (is.null(tracks)) {
    return(NULL)
  }
  if (!is.data.frame(tracks)) {
    stop("`tracks` must be a data frame.", call. = FALSE)
  }
  tracks <- as.data.frame(tracks)
  if (anyNA(names(tracks)) || any(!nzchar(trimws(names(tracks)))) || anyDuplicated(names(tracks))) {
    stop("`tracks` must have unique, nonempty column names.", call. = FALSE)
  }
  track_paths <- cell_tree_paths(tracks, names(paths), "tracks")
  track_names <- setdiff(names(tracks), names(paths))
  if (!length(track_names)) {
    stop("`tracks` must contain at least one numeric track column.", call. = FALSE)
  }
  if (nrow(dplyr::distinct(track_paths)) != nrow(track_paths)) {
    stop("`tracks` contains duplicate annotation paths.", call. = FALSE)
  }
  if (nrow(dplyr::anti_join(track_paths, paths, by = names(paths)))) {
    stop("`tracks` contains unknown annotation paths.", call. = FALSE)
  }
  for (name in track_names) {
    value <- tracks[[name]]
    if (is.logical(value) && is.null(dim(value)) && all(is.na(value))) {
      value <- as.numeric(value)
    }
    if (!cell_tree_numeric(value) || any(is.infinite(value) | is.nan(value))) {
      stop(sprintf("Track `%s` must contain finite numeric values or NA.", name), call. = FALSE)
    }
    tracks[[name]] <- value
  }
  tracks[names(paths)] <- track_paths
  aligned <- dplyr::left_join(paths, tracks, by = names(paths))
  aligned[names(tracks)]
}

cell_tree_draw <- function(
  leaves, prefix_starts, levels, tracks, palette, palcolor, label.size,
  show_proportion, title, show_internal_labels
) {
  track_names <- setdiff(names(tracks), levels)
  gap <- if (length(track_names)) pi / 8 else 0
  boundaries <- seq(gap, 2 * pi - gap, length.out = nrow(leaves) + 1L)
  tip_angle <- (boundaries[-1L] + boundaries[-length(boundaries)]) / 2
  if (nrow(leaves) == 1L) {
    tip_angle <- pi / 2
  }
  depth <- length(levels)
  tip_label <- leaves[[levels[[depth]]]]
  if (show_proportion) {
    percent <- trimws(formatC(100 * leaves$.proportion, digits = 3, format = "fg"))
    tip_label <- paste0(tip_label, " (", percent, "%)")
  }
  coarse <- leaves[[levels[[1L]]]]
  groups <- unique(coarse)
  if (is.list(palcolor)) {
    palcolor <- unlist(palcolor)
  }
  color_names <- names(palcolor)
  if (is.null(color_names)) {
    colors <- thisplot::palette_colors(groups, palette = palette, palcolor = palcolor, NA_keep = TRUE)
  } else {
    if (anyNA(color_names) || any(!nzchar(color_names)) || anyDuplicated(color_names)) {
      stop("Names in `palcolor` must be unique and nonempty.", call. = FALSE)
    }
    colors <- thisplot::palette_colors(groups, palette = palette, NA_keep = TRUE)
    matched <- intersect(groups, color_names)
    colors[matched] <- palcolor[matched]
  }
  radii <- seq(0, 1, length.out = depth + 1L)
  branch_widths <- seq(1.8, 0.65, length.out = depth)
  coarse_starts <- c(prefix_starts[[1L]], length(boundaries))
  sector_angles <- boundaries[coarse_starts]
  sector_inset <- min(0.012, min(diff(sector_angles)) * 0.1)
  sectors <- data.frame(
    start = sector_angles[-length(sector_angles)] + sector_inset,
    end = sector_angles[-1L] - sector_inset, coarse = groups
  )
  nodes <- vector("list", depth)
  edges <- vector("list", depth)
  for (level in seq_len(depth)) {
    first <- prefix_starts[[level]]
    last <- c(first[-1L] - 1L, nrow(leaves))
    angle <- (tip_angle[first] + tip_angle[last]) / 2
    radius <- radii[[level + 1L]]
    parent_radius <- radii[[level]]
    parent_angle <- if (level == 1L) {
      angle
    } else {
      parent <- findInterval(first, prefix_starts[[level - 1L]])
      nodes[[level - 1L]]$angle[parent]
    }
    middle <- (parent_radius + radius) / 2
    nodes[[level]] <- data.frame(
      x = radius * sin(angle), y = radius * cos(angle), angle = angle,
      label = leaves[[levels[[level]]]][first], coarse = coarse[first]
    )
    edges[[level]] <- data.frame(
      x = as.vector(rbind(
        parent_radius * sin(parent_angle), middle * sin(parent_angle),
        middle * sin(angle), radius * sin(angle)
      )),
      y = as.vector(rbind(
        parent_radius * cos(parent_angle), middle * cos(parent_angle),
        middle * cos(angle), radius * cos(angle)
      )),
      group = rep(seq_along(first), each = 4L), coarse = rep(coarse[first], each = 4L),
      depth = level, linewidth = branch_widths[[level]]
    )
  }
  edges <- do.call(rbind, edges)
  tip <- nodes[[depth]]
  band_end <- 1.04
  leaf_inset <- diff(boundaries) * 0.03
  leaf_tiles <- data.frame(
    start = boundaries[-length(boundaries)] + leaf_inset,
    end = boundaries[-1L] - leaf_inset, coarse = coarse,
    group = seq_len(nrow(leaves))
  )
  ring_width <- 0.12
  ring_start <- 1.055
  outer_radius <- if (length(track_names)) {
    ring_start + ring_width * length(track_names)
  } else {
    band_end
  }
  label_radius <- outer_radius + 0.045
  left <- sin(tip_angle) < -1e-8
  tip$angle <- 90 - tip_angle * 180 / pi + ifelse(left, 180, 0)
  tip$hjust <- as.numeric(left)
  tip$end_x <- label_radius * sin(tip_angle)
  tip$end_y <- label_radius * cos(tip_angle)
  tip$label <- tip_label
  tick_end <- label_radius - 0.019
  tick_start <- if (length(track_names)) outer_radius - 0.012 else band_end
  ticks <- data.frame(
    x = tick_start * sin(tip_angle), y = tick_start * cos(tip_angle),
    xend = tick_end * sin(tip_angle), yend = tick_end * cos(tip_angle), coarse = coarse
  )

  p <- ggplot2::ggplot() +
    get_namespace_fun("ggforce", "geom_arc_bar")(
      data = sectors,
      ggplot2::aes(
        x0 = 0, y0 = 0, r0 = 0, r = 1, start = .data$start,
        end = .data$end, fill = .data$coarse
      ),
      color = NA, alpha = 0.13, show.legend = FALSE
    ) +
    get_namespace_fun("ggforce", "geom_arc_bar")(
      data = leaf_tiles,
      ggplot2::aes(
        x0 = 0, y0 = 0, r0 = 1, r = band_end, start = .data$start,
        end = .data$end, fill = .data$coarse, group = .data$group
      ),
      color = NA, show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = colors) +
    get_namespace_fun("ggforce", "geom_bezier")(
      data = edges,
      ggplot2::aes(
        x = .data$x, y = .data$y, color = .data$coarse,
        group = interaction(.data$depth, .data$group), linewidth = .data$linewidth
      ),
      n = 25, lineend = "round"
    ) +
    ggplot2::scale_linewidth_identity() +
    ggplot2::geom_segment(
      data = ticks,
      ggplot2::aes(
        x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
        color = .data$coarse
      ),
      linewidth = 0.38, show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(name = NULL, values = colors, breaks = groups) +
    ggplot2::guides(color = ggplot2::guide_legend(
      order = 1, nrow = ceiling(length(groups) / 5),
      override.aes = list(linewidth = 1.3)
    ))
  if (show_internal_labels && depth > 2L) {
    internal <- do.call(rbind, nodes[seq.int(2L, depth - 1L)])
    left <- sin(internal$angle) < -1e-8
    internal$angle <- 90 - internal$angle * 180 / pi + ifelse(left, 180, 0)
    internal$hjust <- ifelse(left, 1.04, -0.04)
    p <- p + ggplot2::geom_text(
      data = internal,
      ggplot2::aes(
        x = .data$x, y = .data$y, label = .data$label,
        angle = .data$angle, hjust = .data$hjust
      ),
      size = 0.75 * label.size, vjust = -0.45, color = "grey30"
    )
  }
  if (length(track_names)) {
    p <- cell_tree_rings(
      p + ggnewscale::new_scale_fill(),
      tracks, track_names, boundaries, ring_start, ring_width, label.size
    )
  }
  padding <- max(0.12, min(0.8, max(nchar(tip_label, type = "width")) * label.size * 0.009))
  limit <- c(-1, 1) * (label_radius + padding)
  p + ggplot2::geom_text(
    data = tip,
    ggplot2::aes(
      x = .data$end_x, y = .data$end_y, label = .data$label,
      angle = .data$angle, hjust = .data$hjust
    ),
    size = label.size, color = "grey15"
  ) +
    ggplot2::coord_fixed(xlim = limit, ylim = limit, expand = FALSE, clip = "off") +
    ggplot2::labs(title = title) +
    ggplot2::theme_void() +
    theme_scop(base_size = 10) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.margin = ggplot2::margin(t = 22),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.caption = ggplot2::element_text(color = "grey30", hjust = 0.5),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(10, 14, 10, 14)
    )
}

cell_tree_rings <- function(p, tracks, track_names, boundaries, ring_start, ring_width,
                            label.size) {
  tiles <- do.call(rbind, lapply(seq_along(track_names), function(i) {
    data.frame(
      x0 = 0, y0 = 0, r0 = ring_start + ring_width * (i - 1L),
      r = ring_start + ring_width * i - 0.012, start = boundaries[-length(boundaries)],
      end = boundaries[-1L], value = tracks[[track_names[[i]]]],
      group = (i - 1L) * nrow(tracks) + seq_len(nrow(tracks))
    )
  }))
  bound <- max(0, abs(tiles$value), na.rm = TRUE)
  if (bound == 0) {
    bound <- 1
  }
  breaks <- pretty(c(0, bound), n = 3)
  outer_break <- max(breaks[breaks <= bound])
  if (outer_break <= 0) {
    outer_break <- bound
  }
  breaks <- c(-outer_break, 0, outer_break)
  scale_by <- 1
  labels <- ggplot2::waiver()
  if (bound > .Machine$double.xmax / 2) {
    scale_by <- bound
    labels <- format(breaks, trim = TRUE)
    labels[2L] <- "0"
  }
  p <- p + get_namespace_fun("ggforce", "geom_arc_bar")(
    data = tiles,
    ggplot2::aes(
      x0 = .data$x0, y0 = .data$y0, r0 = .data$r0, r = .data$r,
      start = .data$start, end = .data$end, fill = .data$value / scale_by, group = .data$group
    ),
    color = "white", linewidth = 0.23
  ) +
    ggplot2::scale_fill_gradient2(
      low = "#2C6B9E", mid = "#F8F7F4", high = "#BD4E42",
      limits = c(-bound, bound) / scale_by, na.value = "#A9ADB4",
      breaks = breaks / scale_by, labels = labels,
      guide = ggplot2::guide_colorbar(
        order = 2,
        theme = ggplot2::theme(
          legend.key.width = grid::unit(96, "pt"),
          legend.key.height = grid::unit(8, "pt"),
          legend.title.position = "top"
        )
      )
    ) +
    ggplot2::labs(fill = NULL)
  caption <- if (anyNA(tiles$value)) "Gray tiles: missing values"
  if (length(track_names) <= 8L && all(nchar(track_names, type = "width") <= 24L)) {
    track_labels <- data.frame(
      x = 0, y = ring_start + ring_width * (seq_along(track_names) - 0.5), label = track_names
    )
    p <- p + ggplot2::geom_text(
      data = track_labels,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = min(2.65, label.size), color = "grey30"
    )
  } else {
    caption <- c(paste0("Tracks (inner to outer): ", paste(track_names, collapse = ", ")), caption)
  }
  if (length(caption)) {
    p <- p + ggplot2::labs(caption = paste(caption, collapse = "\n"))
  }
  p
}
