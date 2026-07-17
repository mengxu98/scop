scop_expr_input <- function(object, assay = NULL, layer = "data") {
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

resolve_scop_features <- function(mat, features = NULL, nfeatures = 2000) {
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
  vars <- sparse_row_vars(mat)
  head(
    names(sort(vars, decreasing = TRUE)),
    min(length(vars), as.integer(nfeatures))
  )
}

scop_scale_features <- function(mat) {
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

sparse_row_vars <- function(mat) {
  mat <- methods::as(Matrix::Matrix(mat, sparse = TRUE), "dgCMatrix")
  mu <- Matrix::rowMeans(mat)
  mu2 <- Matrix::rowMeans(mat^2)
  out <- pmax(mu2 - mu^2, 0)
  names(out) <- rownames(mat)
  out
}

fitdevo_spearman_weights <- function(scaled, target) {
  if (any(!is.finite(target))) {
    return(stats::setNames(numeric(nrow(scaled)), rownames(scaled)))
  }
  check_r("matrixStats", verbose = FALSE)
  target_rank <- rank(target, ties.method = "average")
  target_centered <- target_rank - mean(target_rank)
  ranked <- matrixStats::rowRanks(as.matrix(scaled), ties.method = "average")
  rank_sum <- rowSums(ranked)
  rank_sumsq <- rowSums(ranked^2)
  rank_ss <- rank_sumsq - rank_sum^2 / ncol(ranked)
  out <- as.vector(ranked %*% target_centered) /
    sqrt(rank_ss * sum(target_centered^2))
  out[!is.finite(out)] <- 0
  stats::setNames(out, rownames(scaled))
}

fitdevo_score <- function(mat, target = NULL) {
  scaled <- scop_scale_features(mat)
  if (is.null(target)) {
    detection <- Matrix::rowMeans(
      methods::as(Matrix::Matrix(mat, sparse = TRUE), "dgCMatrix") > 0
    )
    weights <- sqrt(pmax(sparse_row_vars(mat), 0)) * (1 - pmin(detection, 0.99))
  } else {
    weights <- fitdevo_spearman_weights(scaled, target)
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

vector_field <- function(emb, pca, grid.n = 30, arrow.p = 0.9, arrow.ol = 1.5) {
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
  arrows <- vector_weighted_arrows(grid_df, emb, p = arrow.p, ol = arrow.ol)
  list(
    score = signal,
    embedding = emb,
    grid = grid_df,
    arrows = arrows,
    cell_grid = cell_grid,
    status = "success"
  )
}

vector_weighted_arrows <- function(grid_df, emb, p = 0.9, ol = 1.5) {
  if (nrow(grid_df) < 2L) {
    return(NULL)
  }
  centers <- as.matrix(grid_df[, c("x", "y"), drop = FALSE])
  scores <- grid_df$score
  dist_mat <- as.matrix(stats::dist(centers))
  positive_dist <- dist_mat[dist_mat > 0]
  if (length(positive_dist) == 0L) {
    return(NULL)
  }
  one <- min(positive_dist) * ol
  emb_range <- apply(emb, 2, range, na.rm = TRUE)
  out <- lapply(seq_len(nrow(centers)), function(i) {
    vec <- sweep(centers, 2, centers[i, ], "-")
    vec_norm <- t(apply(vec, 1, function(x) {
      len <- sqrt(sum(x^2))
      if (!is.finite(len) || len == 0) {
        return(c(0, 0))
      }
      x / len * one
    }))
    distance_weight <- p^(rank(dist_mat[i, ], ties.method = "first") - 1)
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
  scaled <- scop_scale_features(mat)
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
