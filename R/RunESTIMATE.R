#' @title Run ESTIMATE tumor microenvironment scoring
#'
#' @description
#' Compute ESTIMATE stromal, immune, combined ESTIMATE, and tumor-purity scores
#' from bulk or pseudo-bulk expression data. The implementation uses the
#' original stromal and immune ESTIMATE signatures with an in-package scoring
#' routine, so it does not require the external `estimate` package.
#'
#' @md
#' @param object Optional expression matrix, `SummarizedExperiment`, or `Seurat`
#' object.
#' @param count_matrix Optional expression matrix with genes in rows and samples
#' in columns. Used when `object` is not provided as a matrix.
#' @param assay Assay used for `Seurat` input.
#' @param layer Assay layer used for `Seurat` input.
#' @param bulk_assay Assay name used for `SummarizedExperiment` input.
#' @param sample.by Metadata column used to aggregate `Seurat` cells into
#' pseudo-bulk samples.
#' @param group.by Optional metadata column used for grouping pseudo-bulk
#' samples. If `sample.by` is not supplied, `group.by` is used as the aggregation
#' variable.
#' @param aggregate_fun Pseudo-bulk aggregation function for `Seurat` input.
#' @param platform Platform label stored in the result metadata.
#' @param filter_common_genes Whether to restrict the expression matrix to the
#' ESTIMATE common-gene universe before scoring.
#' @param min_sig_genes Minimum number of observed stromal and immune signature
#' genes required after filtering.
#' @param purity Whether to calculate ESTIMATE tumor purity. The purity formula
#' was calibrated in the original ESTIMATE work primarily for Affymetrix data;
#' for RNA-seq it is best interpreted cautiously.
#' @param verbose Whether to print progress messages.
#'
#' @return
#' A result bundle for matrix input, the modified `SummarizedExperiment` for
#' `SummarizedExperiment` input, or the modified `Seurat` object for `Seurat`
#' input.
#'
#' @references
#' Yoshihara et al. (2013) <doi:10.1038/ncomms3612>.
#' Barbie et al. (2009) <doi:10.1038/nature08460>.
#'
#' @export
RunESTIMATE <- function(
  object = NULL,
  count_matrix = NULL,
  assay = NULL,
  layer = "data",
  bulk_assay = "counts",
  sample.by = NULL,
  group.by = NULL,
  aggregate_fun = c("mean", "sum"),
  platform = c("rnaseq", "affymetrix", "agilent", "illumina"),
  filter_common_genes = TRUE,
  min_sig_genes = 10,
  purity = TRUE,
  verbose = TRUE
) {
  aggregate_fun <- match.arg(aggregate_fun)
  platform <- match.arg(platform)
  input_object <- object
  sample_metadata <- NULL

  if (is.null(count_matrix)) {
    if (methods::is(object, "SummarizedExperiment")) {
      count_matrix <- SummarizedExperiment::assay(object, bulk_assay)
    } else if (inherits(object, c("matrix", "data.frame", "Matrix"))) {
      count_matrix <- object
      input_object <- NULL
    } else if (inherits(object, "Seurat")) {
      pseudo <- estimate_seurat_pseudobulk(
        object = object,
        assay = assay,
        layer = layer,
        sample.by = sample.by,
        group.by = group.by,
        aggregate_fun = aggregate_fun
      )
      count_matrix <- pseudo$matrix
      sample_metadata <- pseudo$sample_metadata
    }
  }

  if (is.null(count_matrix)) {
    log_message(
      paste(
        "{.arg object} must be a {.cls SummarizedExperiment}, {.cls Seurat},",
        "or expression matrix, or {.arg count_matrix} must be provided."
      ),
      message_type = "error"
    )
  }

  bundle <- run_estimate_bundle(
    count_matrix = count_matrix,
    platform = platform,
    filter_common_genes = filter_common_genes,
    min_sig_genes = min_sig_genes,
    purity = purity,
    sample_metadata = sample_metadata,
    input = list(
      assay = assay,
      layer = layer,
      bulk_assay = bulk_assay,
      sample.by = sample.by,
      group.by = group.by,
      aggregate_fun = aggregate_fun
    ),
    verbose = verbose
  )

  if (methods::is(input_object, "SummarizedExperiment")) {
    return(store_meta(input_object, "ESTIMATE", bundle))
  }
  if (inherits(input_object, "Seurat")) {
    input_object@tools$ESTIMATE <- bundle
    return(input_object)
  }
  bundle
}

run_estimate_bundle <- function(
  count_matrix,
  platform = "rnaseq",
  filter_common_genes = TRUE,
  min_sig_genes = 10,
  purity = TRUE,
  sample_metadata = NULL,
  input = list(),
  verbose = TRUE
) {
  expr <- estimate_prepare_matrix(count_matrix)
  log_message(
    "ESTIMATE input: {.val {nrow(expr)}} genes x {.val {ncol(expr)}} samples",
    verbose = verbose
  )

  signatures <- estimate_get_signatures()
  if (isTRUE(filter_common_genes)) {
    common <- intersect(rownames(expr), signatures$common_genes)
    if (length(common) == 0L) {
      log_message(
        "No genes overlap the ESTIMATE common-gene set.",
        message_type = "error"
      )
    }
    expr <- expr[common, , drop = FALSE]
    log_message(
      "Use {.val {length(common)}} ESTIMATE common genes for scoring",
      verbose = verbose
    )
  }

  overlaps <- list(
    StromalScore = intersect(signatures$stromal_signature, rownames(expr)),
    ImmuneScore = intersect(signatures$immune_signature, rownames(expr))
  )
  n_overlap <- lengths(overlaps)
  if (any(n_overlap < min_sig_genes)) {
    log_message(
      paste0(
        "Too few ESTIMATE signature genes after filtering: ",
        "StromalScore=", n_overlap[["StromalScore"]],
        ", ImmuneScore=", n_overlap[["ImmuneScore"]],
        ". Increase gene overlap or lower {.arg min_sig_genes}."
      ),
      message_type = "error"
    )
  }
  log_message(
    "ESTIMATE signature overlap: stromal {.val {n_overlap[['StromalScore']]}}, immune {.val {n_overlap[['ImmuneScore']]}}",
    verbose = verbose
  )

  scores <- estimate_ssgsea_scores(
    expr = expr,
    gene_sets = overlaps
  )
  scores <- as.data.frame(scores, check.names = FALSE)
  scores$ESTIMATEScore <- scores$StromalScore + scores$ImmuneScore
  scores$TumorPurity <- if (isTRUE(purity)) {
    purity_score <- cos(0.6049872018 + 0.0001467884 * scores$ESTIMATEScore)
    purity_score[purity_score < 0] <- NA_real_
    purity_score
  } else {
    rep(NA_real_, nrow(scores))
  }
  scores <- scores[, estimate_score_columns(), drop = FALSE]

  list(
    status = "success",
    scores = scores,
    results = estimate_scores_long(scores),
    details = list(
      engine = "scop::RunESTIMATE",
      package_backend = FALSE,
      filtered_expression_dim = dim(expr),
      signature_overlap = n_overlap,
      sample_metadata = sample_metadata,
      source = signatures$source
    ),
    parameters = c(
      list(
        method = "ESTIMATE",
        backend = "scop",
        platform = platform,
        filter_common_genes = filter_common_genes,
        min_sig_genes = min_sig_genes,
        purity = purity
      ),
      input
    )
  )
}

estimate_score_columns <- function() {
  c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity")
}

estimate_get_signatures <- function() {
  sig <- estimate_signatures
  list(
    stromal_signature = unique(toupper(sig$stromal_signature)),
    immune_signature = unique(toupper(sig$immune_signature)),
    common_genes = unique(toupper(sig$common_genes)),
    source = sig$source
  )
}

estimate_prepare_matrix <- function(x) {
  if (is.null(x)) {
    log_message("{.arg count_matrix} must be provided.", message_type = "error")
  }
  mat <- as.matrix(x)
  dim_names <- dimnames(mat)
  mat <- suppressWarnings(matrix(
    as.numeric(mat),
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dim_names
  ))
  if (is.null(rownames(mat)) || any(!nzchar(rownames(mat)))) {
    log_message(
      "{.arg count_matrix} must have gene names in rownames.",
      message_type = "error"
    )
  }
  if (is.null(colnames(mat)) || any(!nzchar(colnames(mat)))) {
    log_message(
      "{.arg count_matrix} must have sample names in colnames.",
      message_type = "error"
    )
  }
  mat[!is.finite(mat)] <- 0
  gene_names <- toupper(rownames(mat))
  mat <- rowsum(mat, group = gene_names, reorder = FALSE) /
    as.numeric(table(factor(gene_names, levels = unique(gene_names))))
  mat <- as.matrix(mat)
  keep <- rowSums(abs(mat), na.rm = TRUE) > 0
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) == 0L || ncol(mat) == 0L) {
    log_message(
      "No usable expression values remain after preprocessing.",
      message_type = "error"
    )
  }
  mat
}

estimate_ssgsea_scores <- function(expr, gene_sets) {
  ranked <- apply(expr, 2, function(x) {
    rank(as.numeric(x), ties.method = "average", na.last = "keep")
  })
  ranked <- as.matrix(ranked)
  ranked[!is.finite(ranked)] <- 0
  ranked <- 10000 * ranked / nrow(ranked)
  rownames(ranked) <- rownames(expr)
  colnames(ranked) <- colnames(expr)

  scores <- matrix(
    NA_real_,
    nrow = ncol(ranked),
    ncol = length(gene_sets),
    dimnames = list(colnames(ranked), names(gene_sets))
  )
  for (set_i in seq_along(gene_sets)) {
    common_genes <- intersect(gene_sets[[set_i]], rownames(ranked))
    if (length(common_genes) == 0L) {
      next
    }
    for (sample_i in seq_len(ncol(ranked))) {
      ord <- order(ranked[, sample_i], decreasing = TRUE)
      ordered <- ranked[ord, sample_i]
      names(ordered) <- rownames(ranked)[ord]
      hit_ind <- names(ordered) %in% common_genes
      no_hit_ind <- !hit_ind
      if (!any(hit_ind) || !any(no_hit_ind)) {
        next
      }
      ordered_weight <- ordered^0.25
      hit_exp <- ordered_weight[hit_ind]
      if (sum(hit_exp) <= 0) {
        next
      }
      no_hit_penalty <- cumsum(as.numeric(no_hit_ind) / sum(no_hit_ind))
      hit_reward <- cumsum((as.numeric(hit_ind) * ordered_weight) / sum(hit_exp))
      scores[sample_i, set_i] <- sum(hit_reward - no_hit_penalty)
    }
  }
  scores
}

estimate_scores_long <- function(scores) {
  df <- as.data.frame(as.table(as.matrix(scores)), stringsAsFactors = FALSE)
  colnames(df) <- c("sample", "score_type", "score")
  df$method <- "ESTIMATE"
  df
}

estimate_seurat_pseudobulk <- function(
  object,
  assay = NULL,
  layer = "data",
  sample.by = NULL,
  group.by = NULL,
  aggregate_fun = c("mean", "sum"),
  features = NULL
) {
  aggregate_fun <- match.arg(aggregate_fun)
  if (is.null(sample.by) && is.null(group.by)) {
    log_message(
      "{.arg sample.by} or {.arg group.by} must be provided for {.cls Seurat} pseudo-bulk ESTIMATE scoring.",
      message_type = "error"
    )
  }
  meta <- object@meta.data
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  if (!is.null(sample.by) && !sample.by %in% colnames(meta)) {
    log_message("{.arg sample.by} is not in Seurat metadata.", message_type = "error")
  }
  if (!is.null(group.by) && !group.by %in% colnames(meta)) {
    log_message("{.arg group.by} is not in Seurat metadata.", message_type = "error")
  }

  sample_vec <- if (!is.null(sample.by)) meta[[sample.by]] else meta[[group.by]]
  sample_vec <- as.character(sample_vec)
  if (any(is.na(sample_vec) | !nzchar(sample_vec))) {
    log_message(
      "Pseudo-bulk sample labels contain missing or empty values.",
      message_type = "error"
    )
  }
  group_vec <- if (!is.null(group.by)) as.character(meta[[group.by]]) else sample_vec
  if (any(is.na(group_vec) | !nzchar(group_vec))) {
    log_message(
      "Pseudo-bulk group labels contain missing or empty values.",
      message_type = "error"
    )
  }

  sample_group <- unique(data.frame(
    sample = sample_vec,
    group = group_vec,
    stringsAsFactors = FALSE
  ))
  if (any(duplicated(sample_group$sample))) {
    bad <- unique(sample_group$sample[duplicated(sample_group$sample)])
    log_message(
      "Each {.arg sample.by} level must map to one {.arg group.by} value. Ambiguous samples: {.val {bad}}",
      message_type = "error"
    )
  }

  expr <- GetAssayData5(object, assay = assay, layer = layer)
  if (!is.null(features)) {
    feature_key <- toupper(rownames(expr))
    features_key <- toupper(features)
    keep <- feature_key %in% features_key
    expr <- expr[keep, , drop = FALSE]
  }
  if (!inherits(expr, c("dgCMatrix", "matrix", "Matrix"))) {
    expr <- as_matrix(expr)
  }
  sample_factor <- factor(sample_vec, levels = unique(sample_vec))
  design <- Matrix::sparseMatrix(
    i = seq_along(sample_factor),
    j = as.integer(sample_factor),
    x = 1,
    dims = c(length(sample_factor), nlevels(sample_factor)),
    dimnames = list(NULL, levels(sample_factor))
  )
  pseudo <- expr %*% design
  if (identical(aggregate_fun, "mean")) {
    pseudo <- sweep(
      pseudo,
      2,
      as.numeric(table(sample_factor)[colnames(pseudo)]),
      "/"
    )
  }
  pseudo <- as.matrix(pseudo)
  colnames(pseudo) <- levels(sample_factor)
  sample_metadata <- sample_group[match(colnames(pseudo), sample_group$sample), , drop = FALSE]
  rownames(sample_metadata) <- sample_metadata$sample
  list(matrix = pseudo, sample_metadata = sample_metadata)
}
