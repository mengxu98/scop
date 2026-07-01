# Internal spatial helpers shared by SCOP spatial wrappers.

scop_spatial_resolve_coord_cols <- function(srt, coord.cols = c("col", "row")) {
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

scop_spatial_metadata_coords <- function(srt, coord.cols = c("col", "row")) {
  coord.cols <- scop_spatial_resolve_coord_cols(srt, coord.cols = coord.cols)
  data.frame(
    x = suppressWarnings(as.numeric(srt@meta.data[[coord.cols[1L]]])),
    y = suppressWarnings(as.numeric(srt@meta.data[[coord.cols[2L]]])),
    row.names = rownames(srt@meta.data),
    stringsAsFactors = FALSE
  )
}

scop_spatial_empty_plot <- function(
  message,
  title = NULL,
  theme_use = "theme_blank",
  theme_args = list()
) {
  theme_obj <- scop_spatial_theme(theme_use = theme_use, theme_args = theme_args)
  ggplot2::ggplot(data.frame(x = 0, y = 0, label = message), ggplot2::aes(x, y)) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 3.6, color = "grey35") +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
    theme_obj
}

scop_spatial_theme <- function(
  theme_use = "theme_blank",
  theme_args = list(),
  show_axes = FALSE
) {
  if (is.null(theme_use)) {
    theme_obj <- ggplot2::theme_minimal()
  } else if (is.function(theme_use)) {
    theme_obj <- do.call(theme_use, theme_args)
  } else {
    theme_obj <- do.call(theme_use, theme_args)
  }
  if (isFALSE(show_axes)) {
    theme_obj <- theme_obj +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )
  }
  theme_obj
}

scop_spatial_crop_limits <- function(x, y, pad_fraction = 0.04, min_pad = 0) {
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

scop_spatial_weight_summary <- function(weights) {
  weights <- as.data.frame(weights, check.names = FALSE)
  if (nrow(weights) == 0L || ncol(weights) == 0L) {
    return(list(
      n_spots = nrow(weights),
      n_types = ncol(weights),
      dominant_counts = data.frame(type = character(), count = integer()),
      max_prop = c(min = NA_real_, median = NA_real_, mean = NA_real_, max = NA_real_)
    ))
  }
  max_idx <- max.col(as.matrix(weights), ties.method = "first")
  dominant <- colnames(weights)[max_idx]
  dominant_counts <- as.data.frame(table(dominant), stringsAsFactors = FALSE)
  colnames(dominant_counts) <- c("type", "count")
  max_prop <- apply(weights, 1, max, na.rm = TRUE)
  list(
    n_spots = nrow(weights),
    n_types = ncol(weights),
    dominant_counts = dominant_counts,
    max_prop = c(
      min = unname(min(max_prop, na.rm = TRUE)),
      median = unname(stats::median(max_prop, na.rm = TRUE)),
      mean = unname(mean(max_prop, na.rm = TRUE)),
      max = unname(max(max_prop, na.rm = TRUE))
    )
  )
}

scop_spatial_normalize_weights <- function(weights) {
  weights <- as.matrix(weights)
  weights[!is.finite(weights) | weights < 0] <- 0
  totals <- rowSums(weights)
  keep <- is.finite(totals) & totals > 0
  weights[keep, ] <- weights[keep, , drop = FALSE] / totals[keep]
  weights[!keep, ] <- 0
  weights
}

scop_spatial_finalize_weights <- function(weights, all_spots) {
  weights <- scop_spatial_normalize_weights(weights)
  if (is.null(rownames(weights)) || is.null(colnames(weights))) {
    log_message("Deconvolution weights must have row and column names", message_type = "error")
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

scop_spatial_add_deconv_metadata <- function(srt, weights, prefix, metadata = NULL) {
  if (is.null(metadata)) {
    metadata <- scop_spatial_finalize_weights(weights, all_spots = colnames(srt))
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

scop_spatial_domain_summary <- function(labels) {
  labels <- as.character(labels)
  labels <- labels[!is.na(labels) & nzchar(labels)]
  tab <- as.data.frame(table(labels), stringsAsFactors = FALSE)
  colnames(tab) <- c("domain", "count")
  tab[order(tab$count, decreasing = TRUE), , drop = FALSE]
}

scop_spatial_feature_summary <- function(features, scores = NULL, n = 20L) {
  features <- as.character(features)
  features <- features[!is.na(features) & nzchar(features)]
  features <- utils::head(features, n)
  out <- data.frame(feature = features, rank = seq_along(features), stringsAsFactors = FALSE)
  if (!is.null(scores)) {
    out$score <- as.numeric(scores[features])
  }
  out
}

scop_spatial_neighborhood_summary <- function(pair_table, edge_table = NULL) {
  list(
    n_pairs = nrow(pair_table),
    n_edges = if (is.null(edge_table)) NA_integer_ else nrow(edge_table),
    top_pairs = utils::head(pair_table[order(abs(pair_table$estimate), decreasing = TRUE), , drop = FALSE], 10L)
  )
}
