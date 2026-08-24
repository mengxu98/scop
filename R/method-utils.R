resolve_expression_input <- function(object, assay = NULL, layer = "data") {
  if (inherits(object, "Seurat")) {
    assay <- assay %||% SeuratObject::DefaultAssay(object)
    return(list(
      matrix = GetAssayData5(object, assay = assay, layer = layer),
      assay = assay
    ))
  }
  if (inherits(object, c("matrix", "data.frame", "Matrix"))) {
    return(list(matrix = object, assay = assay))
  }
  log_message(
    "{.arg object} must be a Seurat object or expression matrix.",
    message_type = "error"
  )
}

validate_named_param_list <- function(
  x,
  arg,
  require_list = FALSE,
  type_message = "must be a list",
  names_message = "must contain named arguments only"
) {
  if (isTRUE(require_list) && !is.list(x)) {
    log_message(
      paste("{.arg {arg}}", type_message),
      message_type = "error"
    )
  }
  if (length(x) == 0L) {
    return(invisible(TRUE))
  }
  nms <- names(x)
  if (is.null(nms) || any(is.na(nms) | !nzchar(nms))) {
    log_message(
      paste("{.arg {arg}}", names_message),
      message_type = "error"
    )
  }
  invisible(TRUE)
}

validate_named_list <- function(x, arg) {
  validate_named_param_list(
    x,
    arg,
    require_list = TRUE,
    type_message = "must be a named list",
    names_message = "must be a named list"
  )
}

validate_scalar_string <- function(
  x,
  arg,
  require_character = TRUE,
  message = "must be a single non-empty string"
) {
  if (
    (isTRUE(require_character) && !is.character(x)) ||
      is.null(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !nzchar(x)
  ) {
    log_message(
      paste("{.arg {arg}}", message),
      message_type = "error"
    )
  }
  invisible(TRUE)
}

validate_scalar_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    log_message(
      "{.arg {arg}} must be TRUE or FALSE",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

validate_scalar_integer <- function(
  x,
  arg,
  minimum = 1L,
  message = "must be a positive integer"
) {
  if (
    !is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      x < minimum ||
      x != as.integer(x)
  ) {
    log_message(
      paste("{.arg {arg}}", message),
      message_type = "error"
    )
  }
  as.integer(x)
}

validate_positive_integer <- function(x, arg) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 1) {
    log_message(
      "{.arg {arg}} must be a positive integer",
      message_type = "error"
    )
  }
  as.integer(x)
}

resolve_reference_labels <- function(reference, reference_label) {
  if (
    missing(reference_label) ||
      is.null(reference_label) ||
      length(reference_label) != 1L
  ) {
    log_message(
      "{.arg reference_label} must be a single reference metadata column",
      message_type = "error"
    )
  }
  if (!reference_label %in% colnames(reference[[]])) {
    log_message(
      "{.arg reference_label} {.val {reference_label}} is not present in {.arg reference}",
      message_type = "error"
    )
  }
  reference[[reference_label, drop = TRUE]]
}

resolve_common_features <- function(
  srt,
  reference,
  assay,
  reference_assay,
  features = NULL
) {
  common <- intersect(
    rownames(srt[[assay]]),
    rownames(reference[[reference_assay]])
  )
  if (!is.null(features)) {
    common <- intersect(features, common)
  }
  common
}

resolve_spatial_spot_coords <- function(
  srt,
  spot_ids,
  image = NULL,
  coord.cols = c("x", "y"),
  coordinate_space = c("raw", "legacy_display")
) {
  coordinate_space <- match.arg(coordinate_space)
  resolved <- spatial_analysis_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    coordinate_space = coordinate_space,
    image_policy = "strict"
  )
  matched <- match(spot_ids, resolved$data$cell_id)
  if (anyNA(matched)) {
    log_message(
      "Spatial coordinates are missing for one or more requested spots",
      message_type = "error"
    )
  }
  coords <- resolved$data[matched, c("x", "y"), drop = FALSE]
  rownames(coords) <- spot_ids
  attr(coords, "spatial_source") <- resolved$source
  attr(coords, "spatial_transform") <- resolved$transform
  coords
}

collapse_parameter_value <- function(x) {
  if (is.null(x)) {
    return(NA_character_)
  }
  if (length(x) == 0L) {
    return("")
  }
  paste(as.character(x), collapse = ",")
}

pick_case_insensitive_column <- function(x, candidates) {
  nm <- colnames(x)
  hit <- candidates[tolower(candidates) %in% tolower(nm)][1]
  if (is.na(hit)) {
    return(NULL)
  }
  nm[match(tolower(hit), tolower(nm))]
}

pick_numeric_column <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0L) {
    return(rep(NA_real_, nrow(df)))
  }
  suppressWarnings(as.numeric(df[[hit[[1L]]]]))
}

resource_ref <- function(base, path) {
  path <- gsub("^/+", "", path)
  if (dir.exists(base)) {
    return(file.path(base, path))
  }
  paste0(sub("/+$", "", base), "/", path)
}

resolve_method_features <- function(mat, features = NULL, nfeatures = 2000) {
  if (!is.null(features)) {
    features <- intersect(features, rownames(mat))
    if (length(features) == 0L) {
      log_message(
        "No requested {.arg features} are present in the expression matrix.",
        message_type = "error"
      )
    }
    return(features)
  }
  vars <- fast_row_vars(mat, unbiased = FALSE)
  head(
    names(sort(vars, decreasing = TRUE)),
    min(length(vars), as.integer(nfeatures))
  )
}

scale_feature_matrix <- function(mat) {
  mat <- methods::as(Matrix::Matrix(mat, sparse = TRUE), "dgCMatrix")
  mu <- Matrix::rowMeans(mat)
  mu2 <- Matrix::rowMeans(mat^2)
  sd <- sqrt(pmax(mu2 - mu^2, 1e-8))
  x <- scale_sparse_rows_from_stats(
    mat,
    center = as.numeric(mu),
    scale = as.numeric(sd)
  )
  dimnames(x) <- dimnames(mat)
  x
}

fitdevo_spearman_weights <- function(
  scaled,
  target,
  backend = c("cpp", "r")
) {
  backend <- match.arg(backend)
  if (any(!is.finite(target))) {
    return(stats::setNames(numeric(nrow(scaled)), rownames(scaled)))
  }
  target_rank <- rank(target, ties.method = "average")
  target_centered <- target_rank - mean(target_rank)
  if (identical(backend, "cpp")) {
    out <- fitdevo_spearman_weights_cpp(
      scaled = as.matrix(scaled),
      target_centered = target_centered
    )
  } else {
    check_r("matrixStats", verbose = FALSE)
    ranked <- matrixStats::rowRanks(as.matrix(scaled), ties.method = "average")
    rank_sum <- rowSums(ranked)
    rank_sumsq <- rowSums(ranked^2)
    rank_ss <- rank_sumsq - rank_sum^2 / ncol(ranked)
    out <- as.vector(ranked %*% target_centered) /
      sqrt(rank_ss * sum(target_centered^2))
  }
  out[!is.finite(out)] <- 0
  stats::setNames(out, rownames(scaled))
}

fitdevo_score <- function(mat, target = NULL, backend = c("cpp", "r")) {
  backend <- match.arg(backend)
  scaled <- scale_feature_matrix(mat)
  if (is.null(target)) {
    detection <- Matrix::rowMeans(
      methods::as(Matrix::Matrix(mat, sparse = TRUE), "dgCMatrix") > 0
    )
    weights <- sqrt(pmax(fast_row_vars(mat, unbiased = FALSE), 0)) *
      (1 - pmin(detection, 0.99))
  } else {
    weights <- fitdevo_spearman_weights(scaled, target, backend = backend)
  }
  names(weights) <- rownames(mat)
  weights <- weights / (sqrt(sum(weights^2)) + 1e-8)
  score <- as.numeric(crossprod(weights, scaled))
  names(score) <- colnames(mat)
  score <- scale01(score)
  relative <- rank(score, ties.method = "average") / length(score)
  names(relative) <- names(score)
  list(
    scores = score,
    relative = relative,
    weights = weights,
    status = "success"
  )
}

vector_field <- function(
  emb,
  pca,
  grid.n = 30,
  arrow.p = 0.9,
  arrow.ol = 1.5,
  backend = c("cpp", "r")
) {
  backend <- match.arg(backend)
  common <- intersect(rownames(emb), rownames(pca))
  emb <- emb[common, , drop = FALSE]
  pca <- pca[common, , drop = FALSE]
  colnames(emb) <- c("Dim1", "Dim2")
  pca_rank <- apply(pca, 2, rank, ties.method = "average")
  signal <- rowMeans(scale(pca_rank))
  signal <- scale01(signal)
  names(signal) <- common
  delta <- 1e-6
  x_breaks <- seq(
    min(emb[, 1]) - delta,
    max(emb[, 1]) + delta,
    length.out = grid.n + 1L
  )
  y_breaks <- seq(
    min(emb[, 2]) - delta,
    max(emb[, 2]) + delta,
    length.out = grid.n + 1L
  )
  x_centers <- (x_breaks[-1L] + x_breaks[-length(x_breaks)]) / 2
  y_centers <- (y_breaks[-1L] + y_breaks[-length(y_breaks)]) / 2
  gx <- cut(emb[, 1], breaks = x_breaks, labels = FALSE, include.lowest = TRUE)
  gy <- cut(emb[, 2], breaks = y_breaks, labels = FALSE, include.lowest = TRUE)
  cell_grid <- stats::setNames(paste(gx, gy, sep = "_"), common)
  grid_index <- split(seq_len(nrow(emb)), paste(gx, gy, sep = "_"))
  grid_df <- do.call(
    rbind,
    lapply(names(grid_index), function(id) {
      idx <- grid_index[[id]]
      parts <- as.integer(strsplit(id, "_", fixed = TRUE)[[1L]])
      data.frame(
        grid = id,
        x = x_centers[parts[1]],
        y = y_centers[parts[2]],
        score = mean(signal[idx], na.rm = TRUE),
        n = length(idx),
        row.names = NULL
      )
    })
  )
  arrows <- vector_weighted_arrows(
    grid_df,
    emb,
    p = arrow.p,
    ol = arrow.ol,
    backend = backend
  )
  list(
    score = signal,
    embedding = emb,
    grid = grid_df,
    arrows = arrows,
    cell_grid = cell_grid,
    status = "success"
  )
}

vector_weighted_arrows <- function(
  grid_df,
  emb,
  p = 0.9,
  ol = 1.5,
  backend = c("cpp", "r")
) {
  backend <- match.arg(backend)
  if (nrow(grid_df) < 2L) {
    return(NULL)
  }
  centers <- as.matrix(grid_df[, c("x", "y"), drop = FALSE])
  scores <- grid_df$score
  if (identical(backend, "cpp")) {
    emb_range <- c(
      min(emb[, 1], na.rm = TRUE),
      max(emb[, 1], na.rm = TRUE),
      min(emb[, 2], na.rm = TRUE),
      max(emb[, 2], na.rm = TRUE)
    )
    native <- vector_weighted_arrows_cpp(
      centers = centers,
      scores = scores,
      embedding_range = emb_range,
      p = p,
      arrow_length = ol
    )
    if (nrow(native) == 0L) {
      return(NULL)
    }
    native <- as.data.frame(native)
    colnames(native) <- c(
      "index", "x", "y", "dx", "dy", "xend", "yend", "length"
    )
    native$grid <- grid_df$grid[as.integer(native$index)]
    native$index <- NULL
    return(native[, c(
      "grid", "x", "y", "dx", "dy", "xend", "yend", "length"
    ), drop = FALSE])
  }
  # Distances are computed row-wise (dx * dx + dy * dy), matching the
  # C++ backend exactly, so tied distances rank identically on every
  # platform instead of drifting at the last bit inside stats::dist().
  row_distances <- t(apply(centers, 1, function(ci) {
    dx <- centers[, 1] - ci[1]
    dy <- centers[, 2] - ci[2]
    sqrt(dx * dx + dy * dy)
  }))
  positive_dist <- row_distances[row_distances > 0]
  if (length(positive_dist) == 0L) {
    return(NULL)
  }
  one <- min(positive_dist) * ol
  emb_range <- apply(emb, 2, range, na.rm = TRUE)
  out <- lapply(seq_len(nrow(centers)), function(i) {
    vec <- sweep(centers, 2, centers[i, ], "-")
    distances <- row_distances[i, ]
    vec_norm <- t(vapply(seq_len(nrow(vec)), function(j) {
      len <- distances[j]
      if (!is.finite(len) || len == 0) {
        return(c(0, 0))
      }
      vec[j, ] / len * one
    }, numeric(2)))
    distance_weight <- p^(rank(distances, ties.method = "first") - 1)
    score_weight <- scores[i] - scores
    weight <- distance_weight * score_weight
    denom <- sum(abs(weight), na.rm = TRUE)
    if (!is.finite(denom) || denom == 0) {
      return(NULL)
    }
    final_vec <- as.numeric(t(vec_norm) %*% (weight / denom))
    if (any(!is.finite(final_vec)) || sqrt(sum(final_vec^2)) == 0) {
      return(NULL)
    }
    end <- centers[i, ] + final_vec
    end[1] <- min(max(end[1], emb_range[1, 1]), emb_range[2, 1])
    end[2] <- min(max(end[2], emb_range[1, 2]), emb_range[2, 2])
    data.frame(
      grid = grid_df$grid[i],
      x = centers[i, 1],
      y = centers[i, 2],
      dx = final_vec[1],
      dy = final_vec[2],
      xend = end[1],
      yend = end[2],
      length = sqrt(sum(final_vec^2)),
      row.names = NULL
    )
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0L) {
    return(NULL)
  }
  do.call(rbind, out)
}

fwp_score <- function(mat, y = NULL, weights = NULL) {
  scaled <- scale_feature_matrix(mat)
  if (is.null(weights)) {
    weights <- rowMeans(scaled[, y == 1, drop = FALSE]) -
      rowMeans(scaled[, y == 0, drop = FALSE])
  } else {
    weights <- weights[intersect(names(weights), rownames(scaled))]
    scaled <- scaled[names(weights), , drop = FALSE]
  }
  weights[!is.finite(weights)] <- 0
  weights <- weights / (sqrt(sum(weights^2)) + 1e-8)
  score <- as.numeric(crossprod(weights, scaled))
  names(score) <- colnames(scaled)
  score <- scale01(score)
  list(score = score, weights = weights, status = "success")
}

ordered_numeric <- function(x) {
  if (is.numeric(x)) {
    return(as.numeric(x))
  }
  if (is.factor(x)) {
    return(as.numeric(x))
  }
  as.numeric(factor(x, levels = unique(x)))
}

binary_numeric <- function(x, positive = NULL) {
  if (is.logical(x)) {
    return(as.integer(x))
  }
  if (!is.null(positive)) {
    if (!positive %in% x) {
      log_message(
        "{.arg positive} must be present in {.arg phenotype.by}.",
        message_type = "error"
      )
    }
    return(as.integer(x == positive))
  }
  vals <- unique(x[!is.na(x)])
  if (length(vals) != 2L) {
    log_message(
      "{.arg phenotype.by} must contain exactly two non-missing classes.",
      message_type = "error"
    )
  }
  as.integer(x == vals[2])
}

scale01 <- function(x) {
  rng <- range(x, finite = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(stats::setNames(rep(0.5, length(x)), names(x)))
  }
  (x - rng[1]) / diff(rng)
}

resolve_reduction_name <- function(object, reduction) {
  reductions <- names(object@reductions)
  if (reduction %in% reductions) {
    return(reduction)
  }
  hit <- reductions[tolower(reductions) == tolower(reduction)]
  if (length(hit) > 0L) {
    return(hit[1])
  }
  hit <- grep(reduction, reductions, ignore.case = TRUE, value = TRUE)
  if (length(hit) > 0L) {
    return(hit[1])
  }
  NULL
}

tool_bundle_get_result <- function(object, tool, result.name = NULL) {
  if (!inherits(object, "Seurat")) {
    log_message("{.arg object} must be a {.cls Seurat} object", message_type = "error")
  }
  bundle <- object@tools[[tool]]
  if (is.null(bundle)) log_message("{.val {tool}} results are absent", message_type = "error")
  if (is.null(result.name) && length(bundle$results) > 1L) {
    log_message("Multiple {.val {tool}} results are stored; select {.arg result.name}", message_type = "error")
  }
  result.name <- result.name %||% bundle$active_result
  result <- bundle$results[[result.name]]
  if (is.null(result)) log_message("Unknown {.val {tool}} result {.val {result.name}}", message_type = "error")
  list(bundle = bundle, result = result, result.name = result.name)
}

reorder_first_columns <- function(df, first_cols) {
  first_cols <- intersect(first_cols, colnames(df))
  df[, c(first_cols, setdiff(colnames(df), first_cols)), drop = FALSE]
}

matrix_values <- function(x) {
  if (methods::is(x, "sparseMatrix")) {
    methods::slot(x, "x")
  } else {
    as.numeric(x)
  }
}

write_tsv <- function(x, path) {
  utils::write.table(x, file = path, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = "")
  invisible(path)
}

pkg_version_safe <- function(pkg) {
  tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA_character_)
}

rename_first_column <- function(df, target, candidates) {
  hit <- intersect(candidates, colnames(df))
  hit <- setdiff(hit, target)
  if (length(hit) > 0L && !target %in% colnames(df)) {
    colnames(df)[match(hit[1L], colnames(df))] <- target
  }
  df
}

set_continuous_color_scale <- function(plot, limits, title, context = "proportion") {
  scale_index <- which(vapply(plot$scales$scales, function(scale) any(scale$aesthetics %in%
    c("colour", "color")), logical(1)))
  if (length(scale_index) != 1L) {
    log_message(
      paste0("Unable to identify the continuous ", context, " color scale"),
      message_type = "error"
    )
  }
  plot$scales$scales[[scale_index]]$limits <- limits
  plot$scales$scales[[scale_index]]$name <- title
  plot$labels$colour <- title
  plot
}

parameters_summary_df <- function(parameters, version_pkgs = character(0)) {
  parameters$scop_version <- pkg_version_safe("scop")
  for (nm in names(version_pkgs)) {
    parameters[[nm]] <- pkg_version_safe(version_pkgs[[nm]])
  }
  data.frame(
    key = names(parameters),
    value = vapply(parameters, collapse_parameter_value, character(1)),
    stringsAsFactors = FALSE
  )
}

validate_result_bundle <- function(bundle, label = "CCC", empty_message = NULL) {
  required <- c(
    "method", "results", "active_result", "long_table", "primary_table",
    "cells", "coordinates", "parameters", "summary", "provenance"
  )
  if (!is.list(bundle) || !all(required %in% names(bundle))) {
    log_message(
      paste0(label, " result bundle is incomplete"),
      message_type = "error"
    )
  }
  if (!is.data.frame(bundle$long_table) || nrow(bundle$long_table) == 0L) {
    log_message(
      empty_message %||% paste0(label, " result bundle is incomplete"),
      message_type = "error"
    )
  }
  invisible(TRUE)
}

numeric_matrix_result <- function(x, row_label = "matrix", col_prefix = "score_",
                                  check_null = FALSE) {
  if (isTRUE(check_null) && is.null(x)) {
    log_message(
      "{.arg {row_label}} is empty.",
      message_type = "error"
    )
  }
  mat <- as.matrix(x)
  dim_names <- dimnames(mat)
  mat <- suppressWarnings(matrix(
    as.numeric(mat),
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dim_names
  ))
  if (is.null(rownames(mat))) {
    log_message(
      "{.arg {row_label}} must have rownames for sample alignment.",
      message_type = "error"
    )
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- paste0(col_prefix, seq_len(ncol(mat)))
  }
  mat[!is.finite(mat)] <- NA_real_
  mat
}
