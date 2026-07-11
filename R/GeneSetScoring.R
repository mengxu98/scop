gene_set_scoring_to_dgC <- function(expr) {
  if (inherits(expr, "Matrix") && !inherits(expr, "dgCMatrix")) {
    methods::as(expr, "dgCMatrix")
  } else if (is.matrix(expr)) {
    methods::as(Matrix::Matrix(expr, sparse = TRUE), "dgCMatrix")
  } else if (!inherits(expr, "dgCMatrix")) {
    expr_mat <- as_matrix(expr)
    methods::as(Matrix::Matrix(expr_mat, sparse = TRUE), "dgCMatrix")
  } else {
    expr
  }
}

normalize_gene_set_scoring_method <- function(method, arg_name = "method") {
  method_map <- c(
    "seurat" = "Seurat",
    "aucell" = "AUCell",
    "ucell" = "UCell",
    "gsva" = "GSVA",
    "ssgsea" = "ssGSEA",
    "zscore" = "zscore",
    "plage" = "PLAGE",
    "vision" = "VISION"
  )
  method_key <- tolower(as.character(method[[1]]))
  if (!method_key %in% names(method_map)) {
    log_message(
      "{.arg {arg_name}} must be one of {.val {unname(method_map)}}",
      message_type = "error"
    )
  }
  unname(method_map[[method_key]])
}

gene_set_scoring_make_score_prefix <- function(prefix = NULL, fallback = NULL) {
  prefix_use <- prefix %||% ""
  prefix_use <- as.character(prefix_use[[1]])
  if (!is.na(prefix_use) && nzchar(prefix_use)) {
    return(prefix_use)
  }
  fallback_use <- fallback %||% ""
  fallback_use <- as.character(fallback_use[[1]])
  if (!is.na(fallback_use) && nzchar(fallback_use)) {
    return(fallback_use)
  }
  ""
}

gene_set_scoring_make_score_colnames <- function(
  feature_names,
  prefix = NULL,
  fallback = NULL
) {
  prefix_use <- gene_set_scoring_make_score_prefix(
    prefix = prefix,
    fallback = fallback
  )
  if (nzchar(prefix_use)) {
    return(make.names(paste(prefix_use, feature_names, sep = "_")))
  }
  make.names(feature_names)
}

gene_set_scoring_make_classification_colname <- function(
  prefix = NULL,
  fallback = NULL
) {
  prefix_use <- gene_set_scoring_make_score_prefix(
    prefix = prefix,
    fallback = fallback
  )
  if (nzchar(prefix_use)) {
    return(make.names(paste0(prefix_use, "_classification")))
  }
  "classification"
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

gene_set_scoring_subset_to_index_union <- function(expr, gene_set_idx) {
  union_idx <- sort(unique(unlist(gene_set_idx, use.names = FALSE)))
  if (!length(union_idx)) {
    return(list(expr = expr, gene_sets = gene_set_idx))
  }
  index_map <- integer(nrow(expr))
  index_map[union_idx] <- seq_along(union_idx)
  list(
    expr = expr[union_idx, , drop = FALSE],
    gene_sets = stats::setNames(
      lapply(gene_set_idx, function(idx) {
        as.integer(index_map[idx])
      }),
      names(gene_set_idx)
    )
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

gene_set_scoring_keep_variable_rows <- function(expr) {
  if (inherits(expr, "dgCMatrix")) {
    keep <- rep(FALSE, nrow(expr))
    if (length(expr@x) > 0L) {
      values_by_row <- split(expr@x, expr@i + 1L)
      keep[as.integer(names(values_by_row))] <- vapply(
        values_by_row,
        function(x) {
          x <- x[is.finite(x)]
          length(x) > 1L && diff(range(x)) > 0
        },
        logical(1)
      )
    }
    return(expr[keep, , drop = FALSE])
  }
  expr_mat <- as_matrix(expr)
  keep <- apply(expr_mat, 1L, function(x) {
    x <- x[is.finite(x)]
    length(x) > 1L && diff(range(x)) > 0
  })
  expr[keep, , drop = FALSE]
}

gene_set_scoring_normalize_chunk_size <- function(
  chunk_size,
  n_genes = NULL,
  n_cells = NULL,
  max_dense_entries = 1e8
) {
  if (
    is.null(chunk_size) || length(chunk_size) == 0L || is.na(chunk_size[[1]])
  ) {
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
    if (
      !is.finite(chunk_size_auto) ||
        chunk_size_auto <= 0 ||
        chunk_size_auto >= n_cells
    ) {
      return(0L)
    }
    return(max(1L, as.integer(chunk_size_auto)))
  }
  if (
    is.character(chunk_size[[1]]) && identical(tolower(chunk_size[[1]]), "auto")
  ) {
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

run_aucell_scores <- function(
  expr_counts,
  gene_sets,
  strategy = c("sparse", "topk", "full"),
  algorithm = c("aucell", "ctxcore"),
  seed = 0L,
  tie_method = c("first", "hash")
) {
  strategy <- match.arg(strategy)
  algorithm <- match.arg(algorithm)
  tie_method <- match.arg(tie_method)
  strategy_id <- c("topk" = 1L, "sparse" = 2L, "full" = 3L)[[strategy]]
  algorithm_id <- c("aucell" = 1L, "ctxcore" = 2L)[[algorithm]]

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
  rank_seed <- if (identical(tie_method, "first")) -1L else as.integer(seed %||% 0L)
  scores <- aucell_auc_sparse(
    expr = expr_counts,
    gene_sets = gene_set_idx,
    auc_max_rank = as.integer(auc_max_rank),
    norm_auc = TRUE,
    strategy = strategy_id,
    algorithm = algorithm_id,
    seed = rank_seed
  )
  dimnames(scores) <- list(colnames(expr_counts), names(gene_set_idx))
  names(dimnames(scores)) <- c("cells", "gene sets")
  scores
}

run_aucell_official_scores <- function(
  expr_counts,
  gene_sets,
  tie_method = c("first", "random"),
  ...
) {
  gene_set_scoring_require_namespace("AUCell")
  tie_method <- match.arg(tie_method)
  expr_rank <- if (identical(tie_method, "first")) {
    expr_mat <- as_matrix(expr_counts)
    rankings <- vapply(
      seq_len(ncol(expr_mat)),
      function(cell) as.integer(rank(-expr_mat[, cell], ties.method = "first")),
      integer(nrow(expr_mat))
    )
    dimnames(rankings) <- dimnames(expr_mat)
    new(
      "aucellResults",
      SummarizedExperiment::SummarizedExperiment(assays = list(ranking = rankings))
    )
  } else {
    AUCell::AUCell_buildRankings(
      as_matrix(expr_counts),
      plotStats = FALSE
    )
  }
  cells_auc <- AUCell::AUCell_calcAUC(
    geneSets = gene_sets,
    rankings = expr_rank,
    ...
  )
  scores <- Matrix::t(AUCell::getAUC(cells_auc))
  as_matrix(scores)
}

# Calculate AUCell scores from an official aucellResults ranking object while
# retaining AUCell's gene-set filtering rules. This is used for compatibility
# routes that require AUCell_buildRankings() semantics (random ties, blocks,
# or keepZeroesAsNA) but can still use the native AUC kernel.
run_aucell_scores_from_official_rankings <- function(
  rankings,
  gene_sets,
  auc_max_rank = ceiling(0.05 * nrow(rankings)),
  norm_auc = TRUE
) {
  gene_set_scoring_require_namespace("AUCell")
  if (!methods::is(rankings, "aucellResults")) {
    stop("rankings must be an aucellResults object from AUCell_buildRankings()", call. = FALSE)
  }
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
    stop("gene_sets must be a named list", call. = FALSE)
  }
  auc_max_rank <- suppressWarnings(as.integer(round(auc_max_rank[[1L]])))
  if (is.na(auc_max_rank) || auc_max_rank <= 0L) {
    stop("auc_max_rank must be a positive value", call. = FALSE)
  }

  rank_mat <- AUCell::getRanking(rankings)
  gene_sets <- lapply(gene_sets, unique)
  gene_sets <- gene_sets[lengths(gene_sets) > 0L]
  if (length(gene_sets) == 0L) {
    stop("No gene sets provided or remaining.", call. = FALSE)
  }
  gene_set_idx <- gene_set_scoring_indices(gene_sets, rownames(rank_mat))
  n_genes <- lengths(gene_sets)
  missing_percent <- (n_genes - lengths(gene_set_idx)) / n_genes
  if (all(missing_percent >= 0.8)) {
    stop(
      "Fewer than 20% of the genes in the gene sets are included in the rankings.",
      " Check whether the gene IDs in the rankings and gene sets match.",
      call. = FALSE
    )
  }
  keep <- missing_percent < 0.8
  gene_set_idx <- gene_set_idx[keep]
  if (length(gene_set_idx) == 0L) {
    stop("No gene sets remain after AUCell gene-ID filtering.", call. = FALSE)
  }

  scores <- aucell_auc_ranked_full(
    rankings = as.matrix(rank_mat),
    gene_sets = gene_set_idx,
    auc_max_rank = auc_max_rank,
    norm_auc = isTRUE(norm_auc)
  )
  dimnames(scores) <- list(colnames(rank_mat), names(gene_set_idx))
  scores
}

run_aucell_scores_from_rankings <- function(rankings, gene_sets) {
  if (!is.matrix(rankings)) {
    rankings <- as.matrix(rankings)
  }
  if (is.null(rownames(rankings)) || is.null(colnames(rankings))) {
    log_message(
      "{.arg rankings} must have cell row names and gene column names",
      message_type = "error"
    )
  }
  gene_set_idx <- gene_set_scoring_indices(
    gene_sets = gene_sets,
    feature_names = colnames(rankings)
  )
  gene_set_idx <- gene_set_idx[lengths(gene_set_idx) > 0L]
  if (length(gene_set_idx) == 0L) {
    log_message(
      "No gene sets retain genes after intersecting with the ranking matrix",
      message_type = "error"
    )
  }
  auc_max_rank <- round(0.05 * ncol(rankings))
  scores <- aucell_auc_ranked(
    rankings = rankings,
    gene_sets = gene_set_idx,
    auc_max_rank = as.integer(auc_max_rank)
  )
  dimnames(scores) <- list(rownames(rankings), names(gene_set_idx))
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

  feature_membership <- Matrix::sparseMatrix(
    i = unlist(feature_idx, use.names = FALSE),
    j = rep.int(seq_along(feature_idx), lengths(feature_idx)),
    x = 1 / rep.int(lengths(feature_idx), lengths(feature_idx)),
    dims = c(nrow(expr_data), length(feature_idx))
  )
  control_membership <- Matrix::sparseMatrix(
    i = unlist(control_idx, use.names = FALSE),
    j = rep.int(seq_along(control_idx), lengths(control_idx)),
    x = 1 / rep.int(lengths(control_idx), lengths(control_idx)),
    dims = c(nrow(expr_data), length(control_idx))
  )
  scores <- Matrix::t(feature_membership) %*%
    expr_data -
    Matrix::t(control_membership) %*% expr_data
  scores <- Matrix::t(scores)
  dimnames(scores) <- list(colnames(expr_data), names(features))
  scores
}

run_gsva_scores <- function(
  expr_counts,
  gene_sets,
  kcdf = c("Poisson", "Gaussian", "none"),
  min_gs_size = 10,
  max_gs_size = 500,
  max_diff = TRUE,
  abs_ranking = FALSE,
  tau = 1,
  chunk_size = NULL
) {
  kcdf <- match.arg(kcdf)

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

  gene_set_scoring_require_namespace("GSVA")
  gene_sets_exact <- lapply(gene_set_idx, function(idx) {
    rownames(expr_counts)[idx]
  })
  param <- GSVA::gsvaParam(
    exprData = expr_counts,
    geneSets = gene_sets_exact,
    minSize = min_gs_size,
    maxSize = max_gs_size,
    kcdf = kcdf,
    tau = tau,
    maxDiff = max_diff,
    absRanking = abs_ranking,
    sparse = identical(kcdf, "none")
  )
  ranks <- GSVA::gsvaRanks(
    param = param,
    verbose = FALSE
  )
  scores <- GSVA::gsvaScores(
    param = ranks,
    verbose = FALSE
  )
  if (inherits(scores, "SummarizedExperiment")) {
    scores <- SummarizedExperiment::assay(scores)
  }
  if (!is.matrix(scores)) {
    scores <- as.matrix(scores)
  }
  score_terms <- rownames(scores)
  if (is.null(score_terms)) {
    score_terms <- names(gene_set_idx)[seq_len(nrow(scores))]
  }
  scores <- Matrix::t(scores)
  dimnames(scores) <- list(colnames(expr_counts), score_terms)
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
  subset <- gene_set_scoring_subset_to_index_union(expr_counts, gene_set_idx)
  expr_counts <- subset$expr
  gene_set_idx <- subset$gene_sets

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
  subset <- gene_set_scoring_subset_to_index_union(expr_counts, gene_set_idx)
  expr_counts <- subset$expr
  gene_set_idx <- subset$gene_sets

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
  expr <- gene_set_scoring_to_dgC(expr)
  if (!is.matrix(scores)) {
    scores <- as_matrix(scores)
  }
  feature_names <- rownames(expr)
  n_cells <- ncol(expr)
  if (n_cells <= 1L) {
    return(scores)
  }

  row_sum <- Matrix::rowSums(expr)
  expr_sq <- expr
  expr_sq@x <- expr_sq@x * expr_sq@x
  row_sq_sum <- Matrix::rowSums(expr_sq)
  row_mean <- row_sum / n_cells
  row_var <- (row_sq_sum - n_cells * row_mean * row_mean) / (n_cells - 1L)
  row_var[row_var < 0 & row_var > -1e-12] <- 0
  row_sd <- sqrt(row_var)
  valid_rows <- is.finite(row_sd) & row_sd > 0

  for (set_name in intersect(rownames(scores), names(gene_sets))) {
    idx <- match(gene_sets[[set_name]], feature_names)
    idx <- unique(idx[!is.na(idx)])
    idx <- idx[valid_rows[idx]]
    if (length(idx) == 0L) {
      next
    }
    inv_sd <- 1 / row_sd[idx]
    ref <- rep(sum(-row_mean[idx] * inv_sd), n_cells)
    ref <- ref +
      as.numeric(Matrix::colSums(
        Matrix::Diagonal(x = inv_sd) %*% expr[idx, , drop = FALSE]
      ))
    ref <- ref / length(idx)
    score <- as.numeric(scores[set_name, ])
    dot <- sum(score * ref, na.rm = TRUE)
    if (is.finite(dot) && dot < 0) {
      scores[set_name, ] <- -scores[set_name, ]
    }
  }
  scores
}

run_metabolism_auc <- function(
  expr_counts,
  gene_sets,
  strategy = c("sparse", "topk", "full")
) {
  run_aucell_scores(
    expr_counts = expr_counts,
    gene_sets = gene_sets,
    strategy = strategy
  )
}

run_vision_scores <- function(
  expr_counts,
  gene_sets,
  scale_by_library = FALSE,
  sig_gene_importance = FALSE
) {
  gene_set_scoring_require_namespace("VISION", install_hint = "YosefLab/VISION")
  Vision_fun <- get_namespace_fun("VISION", "Vision")
  calc_signature_scores_fun <- get_namespace_fun(
    "VISION",
    "calcSignatureScores"
  )
  create_gene_signature_fun <- get_namespace_fun(
    "VISION",
    "createGeneSignature"
  )

  expr_mat <- as_matrix(expr_counts)
  if (isTRUE(scale_by_library)) {
    n_umi <- Matrix::colSums(expr_counts)
    expr_mat <- Matrix::t(Matrix::t(expr_mat) / n_umi) * stats::median(n_umi)
  }

  signatures <- lapply(gene_sets, function(gs) {
    intersect(gs, rownames(expr_mat))
  })
  signatures <- signatures[lengths(signatures) > 0L]
  if (length(signatures) == 0L) {
    log_message(
      "No gene sets retain genes after intersecting with the expression matrix",
      message_type = "error"
    )
  }

  signatures <- Map(
    function(sig_name, sig_genes) {
      create_gene_signature_fun(
        sig_name,
        stats::setNames(rep(1, length(sig_genes)), sig_genes)
      )
    },
    names(signatures),
    signatures
  )

  vis <- Vision_fun(
    expr_mat,
    signatures = signatures,
    projection_methods = character()
  )
  vis <- calc_signature_scores_fun(
    vis,
    sig_gene_importance = sig_gene_importance
  )
  as_matrix(vis@SigScores)
}

run_metabolism_gene_set_indices <- gene_set_scoring_indices
run_metabolism_require_namespace <- gene_set_scoring_require_namespace
