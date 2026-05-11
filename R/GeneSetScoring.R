gene_set_scoring_to_dgC <- function(expr) {
  if (inherits(expr, "Matrix") && !inherits(expr, "dgCMatrix")) {
    methods::as(expr, "dgCMatrix")
  } else if (is.matrix(expr)) {
    methods::as(Matrix::Matrix(expr, sparse = TRUE), "dgCMatrix")
  } else if (!inherits(expr, "dgCMatrix")) {
    methods::as(Matrix::Matrix(as_matrix(expr), sparse = TRUE), "dgCMatrix")
  } else {
    expr
  }
}

gene_set_scoring_indices <- function(gene_sets, feature_names) {
  feature_idx <- stats::setNames(seq_along(feature_names), feature_names)
  stats::setNames(
    lapply(gene_sets, function(gs) {
      idx <- unname(feature_idx[gs])
      idx <- idx[!is.na(idx)]
      unique(as.integer(idx))
    }),
    names(gene_sets)
  )
}

gene_set_scoring_require_namespace <- function(pkg, install_hint = pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    log_message(
      "{.pkg {install_hint}} is required for gene-set scoring",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

gene_set_scoring_drop_invalid_score_sets <- function(scores) {
  keep <- colSums(!is.na(scores)) > 0
  scores[, keep, drop = FALSE]
}

gene_set_scoring_normalize_chunk_size <- function(
  chunk_size,
  n_genes = NULL,
  n_cells = NULL,
  max_dense_entries = 1e8
) {
  if (is.null(chunk_size) || length(chunk_size) == 0L || is.na(chunk_size[[1]])) {
    if (
      is.null(n_genes) ||
        is.null(n_cells) ||
        !is.finite(n_genes) ||
        !is.finite(n_cells) ||
        n_genes <= 0 ||
        n_cells <= 0
    ) {
      return(0L)
    }
    chunk_size_auto <- floor(max_dense_entries / as.numeric(n_genes))
    if (!is.finite(chunk_size_auto) || chunk_size_auto <= 0 || chunk_size_auto >= n_cells) {
      return(0L)
    }
    return(max(1L, as.integer(chunk_size_auto)))
  }
  if (is.character(chunk_size[[1]]) && identical(tolower(chunk_size[[1]]), "auto")) {
    return(gene_set_scoring_normalize_chunk_size(
      chunk_size = NULL,
      n_genes = n_genes,
      n_cells = n_cells,
      max_dense_entries = max_dense_entries
    ))
  }
  if (is.na(chunk_size[[1]])) {
    return(0L)
  }
  chunk_size_num <- suppressWarnings(as.numeric(chunk_size[[1]]))
  if (!is.finite(chunk_size_num) || chunk_size_num <= 0) {
    return(0L)
  }
  as.integer(chunk_size_num)
}

run_aucell_scores <- function(expr_counts, gene_sets, strategy = c("sparse", "topk", "full")) {
  if (!run_aucell_available()) {
    log_message(
      "{.arg backend = 'cpp'} requires the compiled {.pkg scop} shared library. Reinstall the package to build native code.",
      message_type = "error"
    )
  }
  strategy <- match.arg(strategy)
  strategy_id <- c("topk" = 1L, "sparse" = 2L, "full" = 3L)[[strategy]]

  gene_set_idx <- gene_set_scoring_indices(
    gene_sets = gene_sets,
    feature_names = rownames(expr_counts)
  )
  gene_set_idx <- gene_set_idx[lengths(gene_set_idx) > 0L]
  if (length(gene_set_idx) == 0L) {
    log_message(
      "No gene sets retain genes after intersecting with the expression matrix",
      message_type = "error"
    )
  }

  expr_counts <- gene_set_scoring_to_dgC(expr_counts)
  auc_max_rank <- ceiling(0.05 * nrow(expr_counts))
  scores <- aucell_auc_sparse(
    expr = expr_counts,
    gene_sets = gene_set_idx,
    auc_max_rank = as.integer(auc_max_rank),
    norm_auc = TRUE,
    strategy = strategy_id
  )
  dimnames(scores) <- list(colnames(expr_counts), names(gene_set_idx))
  scores
}

run_seurat_module_scores <- function(
  expr_data,
  features,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  seed = 11
) {
  if (!run_seurat_module_available()) {
    log_message(
      "{.arg backend = 'cpp'} requires the compiled {.pkg scop} shared library. Reinstall the package to build native code.",
      message_type = "error"
    )
  }

  set.seed(seed = seed)
  expr_data <- gene_set_scoring_to_dgC(expr_data)
  pool <- pool %||% rownames(expr_data)
  pool <- intersect(pool, rownames(expr_data))
  if (length(pool) == 0L) {
    log_message(
      "{.arg pool} has no overlap with expression features",
      message_type = "error"
    )
  }

  features <- lapply(features, function(x) intersect(x, rownames(expr_data)))
  keep <- lengths(features) > 0L
  features <- features[keep]
  if (length(features) == 0L) {
    log_message(
      "No feature sets retain genes after intersecting with the expression matrix",
      message_type = "error"
    )
  }

  data_avg <- Matrix::rowMeans(expr_data[pool, , drop = FALSE])
  data_avg <- data_avg[order(data_avg)]
  data_cut <- ggplot2::cut_number(
    data_avg + stats::rnorm(n = length(data_avg)) / 1e+30,
    n = nbin,
    labels = FALSE,
    right = FALSE
  )
  names(data_cut) <- names(data_avg)

  control_sets <- lapply(features, function(features_use) {
    ctrl_use <- unlist(
      lapply(
        seq_len(length(features_use)),
        function(j) {
          data_cut[which(data_cut == data_cut[features_use[j]])]
        }
      )
    )
    names(
      sample(
        ctrl_use,
        size = min(ctrl * length(features_use), length(ctrl_use)),
        replace = FALSE
      )
    )
  })

  feature_idx <- lapply(features, function(x) {
    idx <- match(x, rownames(expr_data))
    as.integer(idx[!is.na(idx)])
  })
  control_idx <- lapply(control_sets, function(x) {
    idx <- match(x, rownames(expr_data))
    as.integer(idx[!is.na(idx)])
  })

  scores <- module_score_sparse(
    expr = expr_data,
    feature_sets = feature_idx,
    control_sets = control_idx
  )
  dimnames(scores) <- list(colnames(expr_data), names(features))
  scores
}

run_gsva_scores <- function(
  expr_counts,
  gene_sets,
  kcdf = c("Poisson", "Gaussian"),
  min_gs_size = 10,
  max_gs_size = 500,
  max_diff = TRUE,
  abs_ranking = FALSE,
  tau = 1,
  chunk_size = NULL
) {
  kcdf <- match.arg(kcdf)
  if (!run_gsva_available(kcdf = kcdf)) {
    log_message(
      "{.arg backend = 'cpp'} requires the compiled {.pkg scop} shared library. Reinstall the package to build native code.",
      message_type = "error"
    )
  }

  expr_counts <- gene_set_scoring_to_dgC(expr_counts)
  keep_features <- Matrix::rowSums(expr_counts) > 0
  expr_counts <- expr_counts[keep_features, , drop = FALSE]
  gene_set_idx <- gene_set_scoring_indices(
    gene_sets = gene_sets,
    feature_names = rownames(expr_counts)
  )
  gs_size <- lengths(gene_set_idx)
  gene_set_idx <- gene_set_idx[
    gs_size >= min_gs_size & gs_size <= max_gs_size
  ]
  if (length(gene_set_idx) == 0L) {
    log_message(
      "No gene sets remain after intersecting with the expression matrix",
      message_type = "error"
    )
  }
  chunk_size <- gene_set_scoring_normalize_chunk_size(
    chunk_size = chunk_size,
    n_genes = nrow(expr_counts),
    n_cells = ncol(expr_counts)
  )

  scores <- switch(
    kcdf,
    Gaussian = gsva_gaussian_dense(
      expr = expr_counts,
      gene_sets = gene_set_idx,
      max_diff = max_diff,
      abs_ranking = abs_ranking,
      tau = tau,
      chunk_size = chunk_size
    ),
    Poisson = gsva_poisson_dense(
      expr = expr_counts,
      gene_sets = gene_set_idx,
      max_diff = max_diff,
      abs_ranking = abs_ranking,
      tau = tau,
      chunk_size = chunk_size
    )
  )
  dimnames(scores) <- list(colnames(expr_counts), names(gene_set_idx))
  scores
}

run_ssgsea_scores <- function(
  expr_counts,
  gene_sets,
  min_gs_size = 10,
  max_gs_size = 500,
  alpha = 0.25,
  normalize = TRUE
) {
  if (!run_ssgsea_available()) {
    log_message(
      "{.arg backend = 'cpp'} requires the compiled {.pkg scop} shared library. Reinstall the package to build native code.",
      message_type = "error"
    )
  }

  expr_counts <- gene_set_scoring_to_dgC(expr_counts)
  keep_features <- Matrix::rowSums(expr_counts) > 0
  expr_counts <- expr_counts[keep_features, , drop = FALSE]
  gene_set_idx <- gene_set_scoring_indices(
    gene_sets = gene_sets,
    feature_names = rownames(expr_counts)
  )
  gs_size <- lengths(gene_set_idx)
  gene_set_idx <- gene_set_idx[
    gs_size >= min_gs_size & gs_size <= max_gs_size
  ]
  if (length(gene_set_idx) == 0L) {
    log_message(
      "No gene sets remain after intersecting with the expression matrix",
      message_type = "error"
    )
  }

  scores <- ssgsea_rank_dense(
    expr = expr_counts,
    gene_sets = gene_set_idx,
    alpha = alpha,
    normalize = normalize
  )
  dimnames(scores) <- list(colnames(expr_counts), names(gene_set_idx))
  scores
}

run_zscore_scores <- function(
  expr_counts,
  gene_sets,
  min_gs_size = 10,
  max_gs_size = 500
) {
  if (!run_zscore_available()) {
    log_message(
      "{.arg backend = 'cpp'} requires the compiled {.pkg scop} shared library. Reinstall the package to build native code.",
      message_type = "error"
    )
  }

  expr_counts <- gene_set_scoring_to_dgC(expr_counts)
  keep_features <- Matrix::rowSums(expr_counts) > 0
  expr_counts <- expr_counts[keep_features, , drop = FALSE]
  gene_set_idx <- gene_set_scoring_indices(
    gene_sets = gene_sets,
    feature_names = rownames(expr_counts)
  )
  gs_size <- lengths(gene_set_idx)
  gene_set_idx <- gene_set_idx[
    gs_size >= min_gs_size & gs_size <= max_gs_size
  ]
  if (length(gene_set_idx) == 0L) {
    log_message(
      "No gene sets remain after intersecting with the expression matrix",
      message_type = "error"
    )
  }

  max_size <- if (is.finite(max_gs_size)) {
    as.integer(max_gs_size)
  } else {
    .Machine$integer.max
  }
  scores <- zscore_dense(
    expr = expr_counts,
    gene_sets = gene_set_idx,
    min_size = as.integer(min_gs_size),
    max_size = max_size
  )
  dimnames(scores) <- list(colnames(expr_counts), names(gene_set_idx))
  gene_set_scoring_drop_invalid_score_sets(scores)
}

run_plage_scores <- function(
  expr_counts,
  gene_sets,
  min_gs_size = 10,
  max_gs_size = 500
) {
  if (!run_plage_available()) {
    log_message(
      "{.arg backend = 'cpp'} requires the compiled {.pkg scop} shared library. Reinstall the package to build native code.",
      message_type = "error"
    )
  }

  expr_counts <- gene_set_scoring_to_dgC(expr_counts)
  keep_features <- Matrix::rowSums(expr_counts) > 0
  expr_counts <- expr_counts[keep_features, , drop = FALSE]
  gene_set_idx <- gene_set_scoring_indices(
    gene_sets = gene_sets,
    feature_names = rownames(expr_counts)
  )
  gs_size <- lengths(gene_set_idx)
  gene_set_idx <- gene_set_idx[
    gs_size >= min_gs_size & gs_size <= max_gs_size
  ]
  if (length(gene_set_idx) == 0L) {
    log_message(
      "No gene sets remain after intersecting with the expression matrix",
      message_type = "error"
    )
  }

  max_size <- if (is.finite(max_gs_size)) {
    as.integer(max_gs_size)
  } else {
    .Machine$integer.max
  }
  scores <- plage_dense(
    expr = expr_counts,
    gene_sets = gene_set_idx,
    min_size = as.integer(min_gs_size),
    max_size = max_size
  )
  dimnames(scores) <- list(colnames(expr_counts), names(gene_set_idx))
  gene_set_scoring_drop_invalid_score_sets(scores)
}

orient_plage_scores <- function(scores, expr, gene_sets) {
  if (is.null(scores) || length(gene_sets) == 0L) {
    return(scores)
  }
  expr <- as_matrix(expr)
  scores <- as_matrix(scores)
  feature_names <- rownames(expr)

  for (set_name in intersect(rownames(scores), names(gene_sets))) {
    genes <- intersect(gene_sets[[set_name]], feature_names)
    if (length(genes) == 0L) {
      next
    }
    z <- t(scale(t(expr[genes, , drop = FALSE])))
    z <- z[stats::complete.cases(z), , drop = FALSE]
    if (nrow(z) == 0L) {
      next
    }
    ref <- colMeans(z)
    score <- as.numeric(scores[set_name, ])
    dot <- sum(score * ref, na.rm = TRUE)
    if (is.finite(dot) && dot < 0) {
      scores[set_name, ] <- -scores[set_name, ]
    }
  }
  scores
}

run_aucell_available <- function() {
  exists("aucell_auc_sparse", mode = "function") &&
    isTRUE(is.loaded("_scop_aucell_auc_sparse"))
}

run_seurat_module_available <- function() {
  exists("module_score_sparse", mode = "function") &&
    isTRUE(is.loaded("_scop_module_score_sparse"))
}

run_gsva_available <- function(kcdf = c("Poisson", "Gaussian")) {
  kcdf <- match.arg(kcdf)
  if (identical(kcdf, "Gaussian")) {
    exists("gsva_gaussian_dense", mode = "function") &&
      isTRUE(is.loaded("_scop_gsva_gaussian_dense"))
  } else {
    exists("gsva_poisson_dense", mode = "function") &&
      isTRUE(is.loaded("_scop_gsva_poisson_dense"))
  }
}

run_ssgsea_available <- function() {
  exists("ssgsea_rank_dense", mode = "function") &&
    isTRUE(is.loaded("_scop_ssgsea_rank_dense"))
}

run_zscore_available <- function() {
  exists("zscore_dense", mode = "function") &&
    isTRUE(is.loaded("_scop_zscore_dense"))
}

run_plage_available <- function() {
  exists("plage_dense", mode = "function") &&
    isTRUE(is.loaded("_scop_plage_dense"))
}

run_metabolism_auc <- function(expr_counts, gene_sets, strategy = c("sparse", "topk", "full")) {
  run_aucell_scores(
    expr_counts = expr_counts,
    gene_sets = gene_sets,
    strategy = strategy
  )
}

run_metabolism_gene_set_indices <- gene_set_scoring_indices
run_metabolism_require_namespace <- gene_set_scoring_require_namespace
