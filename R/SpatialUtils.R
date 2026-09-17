# Internal helpers shared by spatial workflow wrappers.

# Resolve the common Seurat-object input used by native spatial functions.
#
# Spatial plotting and analysis APIs take the Seurat object as `object`. `srt`
# remains an accepted deprecated alias until scop 1.0.0.
spatial_resolve_object <- function(object = NULL, srt = NULL) {
  input_arg <- if (is.null(object)) "srt" else "object"
  resolved <- resolve_deprecated_srt(object, srt, is.null(object))
  if (is.null(resolved)) {
    log_message(
      "Provide a {.cls Seurat} object through {.arg object} or {.arg srt}",
      message_type = "error"
    )
  }
  if (!inherits(resolved, "Seurat")) {
    log_message(
      "{.arg {input_arg}} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  resolved
}

spatial_resolve_coord_cols <- function(srt, coord.cols = c("col", "row")) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  meta_cols <- colnames(srt@meta.data)
  requested <- coord.cols
  if (!is.null(requested)) {
    requested <- requested[seq_len(min(2L, length(requested)))]
  }

  default_requested <- is.null(requested) ||
    length(requested) < 2L ||
    identical(requested, c("col", "row"))

  if (!default_requested) {
    missing <- setdiff(requested, meta_cols)
    if (length(missing) > 0L) {
      log_message(
        "Spatial coordinates were not found. Missing metadata column{?s}: {.val {missing}}.",
        message_type = "error"
      )
    }
    return(requested)
  }

  candidates <- list(
    c("x", "y"),
    c("col", "row"),
    c("imagecol", "imagerow"),
    c("pxl_col_in_fullres", "pxl_row_in_fullres")
  )
  for (candidate in candidates) {
    if (all(candidate %in% meta_cols)) {
      return(candidate)
    }
  }

  log_message(
    "Spatial coordinates were not found. Provide a Seurat image or metadata columns {.val x/y} or {.val col/row}.",
    message_type = "error"
  )
}

spatial_metadata_coords <- function(srt, coord.cols = c("col", "row")) {
  coord.cols <- spatial_resolve_coord_cols(srt, coord.cols = coord.cols)
  data.frame(
    x = spatial_coordinate_numeric(srt@meta.data[[coord.cols[1L]]]),
    y = spatial_coordinate_numeric(srt@meta.data[[coord.cols[2L]]]),
    row.names = rownames(srt@meta.data),
    stringsAsFactors = FALSE
  )
}

spatial_coordinate_numeric <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (!is.numeric(x) && !is.character(x)) {
    log_message("Spatial coordinates must be numeric or numeric text", message_type = "error")
  }
  suppressWarnings(as.numeric(x))
}

# Resolve one image for every sample. A scalar image must cover every sample;
# a named image map must cover exactly the requested sample names. `preferred`
# is an internal provenance hint and is used only when it names a covering
# image; it never overrides an explicit user image.
spatial_resolve_sample_images <- function(
  srt,
  sample.by,
  image = NULL,
  preferred = NULL,
  require_image = FALSE
) {
  labels <- as.character(srt[[]][[sample.by]])
  if (length(labels) != ncol(srt) || anyNA(labels) || any(!nzchar(labels))) {
    log_message("{.arg sample.by} must identify every cell or spot", message_type = "error")
  }
  samples <- unique(labels)
  if (!is.null(image)) {
    if (!is.character(image) || length(image) == 0L || anyNA(image) || any(!nzchar(image))) {
      log_message("{.arg image} must be an image name or a named image map", message_type = "error")
    }
    if (length(image) > 1L || !is.null(names(image))) {
      if (is.null(names(image)) || anyDuplicated(names(image)) || !setequal(names(image), samples)) {
        log_message("The named {.arg image} map must cover every sample exactly once", message_type = "error")
      }
    }
  }
  images <- SeuratObject::Images(srt)
  image_cells <- lapply(images, function(nm) SeuratObject::Cells(srt[[nm]]))
  names(image_cells) <- images
  if (isTRUE(require_image) && length(images) == 0L) {
    log_message("One image covering every sample is required", message_type = "error")
  }
  out <- stats::setNames(rep(NA_character_, length(samples)), samples)
  for (sample in samples) {
    cells <- colnames(srt)[labels == sample]
    selected <- if (is.null(image)) NULL else if (is.null(names(image))) image else unname(image[[sample]])
    candidates <- images[vapply(image_cells, function(ids) all(cells %in% ids), logical(1))]
    if (is.null(selected) && !is.null(preferred) && sample %in% names(preferred)) {
      preferred_image <- preferred[[sample]]
      if (length(preferred_image) == 1L && isTRUE(preferred_image %in% candidates)) {
        selected <- preferred_image
      }
    }
    if (is.null(selected) && length(images) > 0L) {
      if (length(candidates) != 1L) {
        log_message(
          "Sample {.val {sample}} requires one image covering all its cells; supply a named {.arg image} map when ambiguous",
          message_type = "error"
        )
      }
      selected <- candidates[[1L]]
    }
    if (isTRUE(require_image) && (is.null(selected) || length(selected) != 1L)) {
      log_message("Sample {.val {sample}} requires one covering image", message_type = "error")
    }
    if (!is.null(selected)) {
      if (length(selected) != 1L || !selected %in% images || !all(cells %in% image_cells[[selected]])) {
        log_message("Selected image is not a covering image and does not cover every cell in sample {.val {sample}}", message_type = "error")
      }
      out[[sample]] <- selected
    }
  }
  out
}

# Resolve every sample independently and retain the source metadata needed by
# schema-v1 spatial results.
spatial_sample_coords <- function(srt, sample.by, image = NULL,
                                  coord.cols = c("col", "row"), coordinate_space = "raw") {
  labels <- as.character(srt[[]][[sample.by]])
  samples <- unique(labels)
  image_map <- spatial_resolve_sample_images(srt, sample.by, image = image)
  data <- sources <- vector("list", length(samples))
  names(data) <- names(sources) <- samples
  for (sample in samples) {
    cells <- colnames(srt)[labels == sample]
    selected <- image_map[[sample]]
    if (is.na(selected)) selected <- NULL
    resolved <- spatial_analysis_coords(
      srt,
      image = selected,
      coord.cols = coord.cols,
      coordinate_space = coordinate_space
    )
    if (!all(cells %in% resolved$data$cell_id)) {
      log_message("Selected image does not cover every cell in sample {.val {sample}}", message_type = "error")
    }
    data[[sample]] <- resolved$data[match(cells, resolved$data$cell_id), , drop = FALSE]
    sources[[sample]] <- c(resolved$source, list(transform = resolved$transform))
  }
  coords <- do.call(rbind, unname(data))
  rownames(coords) <- coords$cell_id
  list(data = coords[colnames(srt), , drop = FALSE], sources = sources, images = image_map)
}

spatial_empty_plot <- function(
  message,
  title = NULL,
  theme_use = "theme_spatial",
  theme_args = list()
) {
  theme_obj <- apply_plot_theme(theme_use = theme_use, theme_args = theme_args)
  ggplot2::ggplot(data.frame(x = 0, y = 0, label = message), ggplot2::aes(x, y)) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 3.6, color = "grey35") +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
    theme_obj
}

#' ggplot2 theme for spatial plots
#'
#' Built on [theme_scop()] with axes, ticks, and panel grid hidden by default.
#' Spatial plot helpers default to `theme_use = "theme_spatial"`.
#'
#' @md
#' @param show_axes Whether to keep axis titles, text, ticks, and grid lines.
#' @param aspect.ratio Optional panel aspect ratio. The default `NULL` lets
#' spatial coordinates determine the geometry through `coord_equal()`.
#' @inheritParams theme_scop
#' @inherit theme_scop return
#' @seealso [theme_scop()]
#' @export
#' @examples
#' theme_spatial()
theme_spatial <- function(
  show_axes = FALSE,
  aspect.ratio = NULL,
  base_size = 12,
  ...
) {
  out <- theme_scop(aspect.ratio = aspect.ratio, base_size = base_size, ...)
  if (isFALSE(show_axes)) {
    out <- out +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )
  }
  out
}

spatial_crop_limits <- function(x, y, pad_fraction = 0.04, min_pad = 0) {
  xr <- range(x, na.rm = TRUE)
  yr <- range(y, na.rm = TRUE)
  xpad <- diff(xr) * pad_fraction
  ypad <- diff(yr) * pad_fraction
  if (!is.finite(xpad) || xpad <= 0) {
    xpad <- max(min_pad, 0.5)
  }
  if (!is.finite(ypad) || ypad <= 0) {
    ypad <- max(min_pad, 0.5)
  }
  list(xlim = xr + c(-xpad, xpad), ylim = yr + c(-ypad, ypad))
}

spatial_weight_summary <- function(weights) {
  weights <- as.matrix(weights)
  assigned <- if (ncol(weights) > 0L) {
    rowSums(is.finite(weights)) == ncol(weights) & rowSums(weights) > 0
  } else {
    rep(FALSE, nrow(weights))
  }
  assigned[is.na(assigned)] <- FALSE
  if (!any(assigned)) {
    return(list(
      n_spots = nrow(weights),
      n_types = ncol(weights),
      dominant_counts = data.frame(type = character(), count = integer()),
      max_prop = c(min = NA_real_, median = NA_real_, mean = NA_real_, max = NA_real_),
      n_assigned = 0L,
      n_unassigned = nrow(weights)
    ))
  }
  valid_weights <- weights[assigned, , drop = FALSE]
  max_idx <- max.col(valid_weights, ties.method = "first")
  dominant <- colnames(weights)[max_idx]
  dominant_counts <- as.data.frame(table(dominant), stringsAsFactors = FALSE)
  colnames(dominant_counts) <- c("type", "count")
  max_prop <- apply(valid_weights, 1, max)
  list(
    n_spots = nrow(weights),
    n_types = ncol(weights),
    dominant_counts = dominant_counts,
    max_prop = c(
      min = unname(min(max_prop, na.rm = TRUE)),
      median = unname(stats::median(max_prop, na.rm = TRUE)),
      mean = unname(mean(max_prop, na.rm = TRUE)),
      max = unname(max(max_prop, na.rm = TRUE))
    ),
    n_assigned = sum(assigned),
    n_unassigned = sum(!assigned)
  )
}

spatial_normalize_weights <- function(weights) {
  weights <- as.matrix(weights)
  if (!is.numeric(weights) || any(!is.finite(weights))) {
    log_message("Deconvolution weights must be finite numeric values", message_type = "error")
  }
  if (any(weights < -sqrt(.Machine$double.eps))) {
    log_message("Deconvolution weights must be non-negative", message_type = "error")
  }
  # Only tolerate floating-point noise around zero; do not hide invalid output.
  weights[weights < 0] <- 0
  totals <- rowSums(weights)
  if (any(!is.finite(totals))) {
    log_message("Deconvolution weight totals must be finite", message_type = "error")
  }
  keep <- is.finite(totals) & totals > 0
  weights[keep, ] <- weights[keep, , drop = FALSE] / totals[keep]
  weights[!keep, ] <- 0
  weights
}

# Keep named colors attached to their categories, independent of row order.
spatial_palette_colors <- function(
  x, palette = "Chinese", palcolor = NULL,
  type = c("auto", "discrete", "continuous"), NA_keep = FALSE, ...
) {
  type <- match.arg(type)
  continuous <- missing(x) || identical(type, "continuous") ||
    (identical(type, "auto") && is.numeric(x))
  if (!continuous && !is.null(palcolor) && !is.null(names(palcolor))) {
    if (anyNA(names(palcolor)) || any(!nzchar(names(palcolor))) || anyDuplicated(names(palcolor))) {
      log_message("Named {.arg palcolor} must have unique non-empty category names", message_type = "error")
    }
    labels <- if (is.factor(x)) levels(x) else unique(as.character(x[!is.na(x)]))
    missing_colors <- setdiff(labels, names(palcolor))
    if (length(missing_colors)) {
      log_message("Named {.arg palcolor} is missing colors for {.val {missing_colors}}", message_type = "error")
    }
    colors <- stats::setNames(as.character(unlist(palcolor[labels])), labels)
    if (isTRUE(NA_keep) && anyNA(x)) colors <- c(colors, "NA" = "grey80")
    return(colors)
  }
  if (missing(x)) {
    return(palette_colors(palette = palette, palcolor = palcolor, type = type, NA_keep = NA_keep, ...))
  }
  palette_colors(x, palette = palette, palcolor = palcolor, type = type, NA_keep = NA_keep, ...)
}

spatial_plot_factor <- function(x) {
  labels <- if (is.factor(x)) levels(x) else unique(as.character(x[!is.na(x)]))
  # Keep missingness separate from a genuine category whose name is "NA".
  factor(as.character(x), levels = labels, ordered = is.ordered(x))
}

spatial_colorbar_guide <- function(direction = c("vertical", "horizontal")) {
  direction <- match.arg(direction)
  horizontal <- identical(direction, "horizontal")
  ggplot2::guide_colorbar(direction = direction, theme = ggplot2::theme(
    legend.key.width = grid::unit(if (horizontal) 38 else 3, "mm"),
    legend.key.height = grid::unit(if (horizontal) 3 else 28, "mm")
  ))
}

spatial_finalize_weights <- function(weights, all_spots) {
  weights <- spatial_normalize_weights(weights)
  valid_names <- function(x) !is.null(x) && !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x)
  if (any(dim(weights) == 0L) || !valid_names(rownames(weights)) ||
      !valid_names(colnames(weights)) || !valid_names(all_spots)) {
    log_message("Deconvolution weights must be non-empty with unique row and column names", message_type = "error")
  }
  if (length(setdiff(rownames(weights), all_spots))) {
    log_message("Deconvolution weights contain unknown spatial spot names", message_type = "error")
  }
  full_weights <- matrix(
    NA_real_,
    nrow = length(all_spots),
    ncol = ncol(weights),
    dimnames = list(all_spots, colnames(weights))
  )
  matched <- intersect(all_spots, rownames(weights))
  full_weights[matched, ] <- weights[matched, , drop = FALSE]

  dominant <- rep(NA_character_, length(all_spots))
  max_prop <- rep(NA_real_, length(all_spots))
  names(dominant) <- all_spots
  names(max_prop) <- all_spots
  for (spot in matched) {
    values <- full_weights[spot, ]
    if (all(is.na(values))) {
      next
    }
    values_cmp <- values
    values_cmp[is.na(values_cmp)] <- 0
    idx <- which.max(values_cmp)
    max_prop[spot] <- values_cmp[idx]
    if (is.finite(values_cmp[idx]) && values_cmp[idx] > 0) {
      dominant[spot] <- colnames(full_weights)[idx]
    }
  }

  list(
    weights = weights,
    full_weights = full_weights,
    dominant = dominant,
    max_prop = max_prop
  )
}

spatial_add_deconv_metadata <- function(srt, weights, prefix, metadata = NULL) {
  if (is.null(metadata)) {
    metadata <- spatial_finalize_weights(weights, all_spots = colnames(srt))
  }
  full_weights <- metadata$full_weights
  meta <- as.data.frame(full_weights, check.names = FALSE)
  prop_cols <- paste0(
    prefix,
    "_prop_",
    make.unique(make.names(colnames(full_weights)), sep = "_")
  )
  colnames(meta) <- prop_cols
  meta[[paste0(prefix, "_dominant_type")]] <- metadata$dominant
  meta[[paste0(prefix, "_max_prop")]] <- metadata$max_prop
  Seurat::AddMetaData(srt, metadata = meta)
}

spatial_domain_summary <- function(labels) {
  labels <- as.character(labels)
  labels <- labels[!is.na(labels) & nzchar(labels)]
  tab <- as.data.frame(table(labels), stringsAsFactors = FALSE)
  colnames(tab) <- c("domain", "count")
  tab[order(tab$count, decreasing = TRUE), , drop = FALSE]
}

spatial_feature_summary <- function(features, scores = NULL, n = 20L) {
  features <- as.character(features)
  features <- features[!is.na(features) & nzchar(features)]
  features <- utils::head(features, n)
  out <- data.frame(feature = features, rank = seq_along(features), stringsAsFactors = FALSE)
  if (!is.null(scores)) {
    out$score <- as.numeric(scores[features])
  }
  out
}

spatial_svf_selection_marker_key <- function() {
  ".scop_RunSpatialVariableFeatures_selection"
}

spatial_begin_svf_selection <- function(srt, assay) {
  assay_object <- srt[[assay]]
  key <- spatial_svf_selection_marker_key()
  old_exists <- key %in% names(assay_object@misc)
  old_value <- if (old_exists) assay_object@misc[[key]] else NULL
  token <- basename(tempfile("scop-svf-selection-"))
  assay_object@misc[[key]] <- structure(
    list(pending_token = token),
    class = "scop_svf_selection_pending"
  )
  srt[[assay]] <- assay_object
  list(
    srt = srt,
    state = list(
      key = key,
      token = token,
      old_exists = old_exists,
      old_value = old_value
    )
  )
}

spatial_restore_svf_selection_marker <- function(srt, assay, state) {
  if (is.null(state)) {
    return(srt)
  }
  assay_object <- srt[[assay]]
  if (isTRUE(state$old_exists)) {
    assay_object@misc[[state$key]] <- state$old_value
  } else {
    assay_object@misc[[state$key]] <- NULL
  }
  srt[[assay]] <- assay_object
  srt
}

spatial_set_active_variable_features <- function(srt, assay, features) {
  features <- unique(as.character(features))
  features <- features[!is.na(features) & nzchar(features)]
  assay_object <- srt[[assay]]

  if (length(features) > 0L) {
    features <- intersect(features, rownames(assay_object))
    if (length(features) == 0L) {
      stop("None of the features specified are present in this assay", call. = FALSE)
    }
  }
  if (inherits(assay_object, "StdAssay")) {
    # Write full, ID-aligned metadata: some Assay5 setters recycle short
    # character selections and treat an empty assignment as a no-op.
    feature_names <- rownames(assay_object)
    empty_labels <- feature_names %in% features
    empty_ranks <- match(feature_names, features)
    names(empty_labels) <- names(empty_ranks) <- feature_names
    assay_object[["var.features"]] <- empty_labels
    assay_object[["var.features.rank"]] <- empty_ranks
  } else {
    SeuratObject::VariableFeatures(assay_object) <- features
  }

  marker_key <- spatial_svf_selection_marker_key()
  marker <- assay_object@misc[[marker_key]]
  if (inherits(marker, "scop_svf_selection_pending") &&
    is.character(marker$pending_token) &&
    length(marker$pending_token) == 1L &&
    !is.na(marker$pending_token) &&
    nzchar(marker$pending_token)) {
    assay_object@misc[[marker_key]] <- structure(
      list(completed_token = marker$pending_token),
      class = "scop_svf_selection_completed"
    )
  }
  srt[[assay]] <- assay_object
  srt
}

spatial_has_explicit_empty_variable_features <- function(
  srt,
  assay,
  expected_token = NULL
) {
  assay_object <- srt[[assay]]
  active_features <- suppressWarnings(
    SeuratObject::VariableFeatures(assay_object)
  )
  active_features <- as.character(active_features)
  active_features <- active_features[
    !is.na(active_features) & nzchar(active_features)
  ]
  if (length(active_features) > 0L) {
    return(FALSE)
  }

  if (!is.null(expected_token)) {
    marker <- assay_object@misc[[spatial_svf_selection_marker_key()]]
    if (
      !inherits(marker, "scop_svf_selection_completed") ||
        !identical(marker$completed_token, expected_token)
    ) {
      return(FALSE)
    }
  }

  if (!inherits(assay_object, "StdAssay")) {
    return(TRUE)
  }

  feature_metadata <- assay_object[[]]
  if (!all(c("var.features", "var.features.rank") %in%
    colnames(feature_metadata))) {
    return(FALSE)
  }
  labels <- feature_metadata[["var.features"]]
  ranks <- feature_metadata[["var.features.rank"]]
  is.logical(labels) &&
    length(labels) == nrow(feature_metadata) &&
    !anyNA(labels) &&
    !any(labels) &&
    is.integer(ranks) &&
    length(ranks) == nrow(feature_metadata) &&
    all(is.na(ranks))
}

spatial_neighborhood_summary <- function(pair_table, edge_table = NULL) {
  list(
    n_pairs = nrow(pair_table),
    n_edges = if (is.null(edge_table)) NA_integer_ else nrow(edge_table),
    top_pairs = utils::head(pair_table[order(abs(pair_table$estimate), decreasing = TRUE), , drop = FALSE], 10L)
  )
}
