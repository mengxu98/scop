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
#' Yoshihara K, Shahmoradgoli M, Martínez E, et al. (2013). Inferring tumour
#' purity and stromal and immune cell admixture from expression data.
#' Nature Communications, 4, 2612. doi:10.1038/ncomms3612
#'
#' Barbie DA, Tamayo P, Boehm JS, et al. (2009). Systematic RNA interference
#' reveals that oncogenic KRAS-driven cancers require TBK1.
#' Nature, 462(7269), 108-112. doi:10.1038/nature08460
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
  signatures <- NULL
  sparse_common <- NULL
  if (isTRUE(filter_common_genes) && identical(class(count_matrix)[[1L]], "dgCMatrix")) {
    sparse_common <- function() {
      signatures <<- estimate_get_signatures()
      signatures$common_genes
    }
  }
  expr <- estimate_prepare_matrix(count_matrix, sparse_common = sparse_common)
  log_message(
    "ESTIMATE input: {.val {nrow(expr)}} genes x {.val {ncol(expr)}} samples",
    verbose = verbose
  )

  if (is.null(signatures)) signatures <- estimate_get_signatures()
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
    expr = as.matrix(expr),
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
  env <- new.env(parent = emptyenv())
  utils::data("estimate_signatures", package = "scop", envir = env)
  sig <- env$estimate_signatures
  list(
    stromal_signature = unique(toupper(sig$stromal_signature)),
    immune_signature = unique(toupper(sig$immune_signature)),
    common_genes = unique(toupper(sig$common_genes)),
    source = sig$source
  )
}

estimate_prepare_matrix <- function(x, sparse_common = NULL) {
  if (is.null(x)) {
    log_message("{.arg count_matrix} must be provided.", message_type = "error")
  }
  if (!is.null(sparse_common) && identical(class(x)[[1L]], "dgCMatrix")) {
    sparse <- estimate_prepare_sparse_matrix(x, sparse_common)
    if (!is.null(sparse)) return(sparse)
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

estimate_prepare_sparse_matrix <- function(x, sparse_common) {
  genes <- rownames(x)
  samples <- colnames(x)
  if (is.null(genes) || is.null(samples) || anyNA(genes) ||
    !is.null(attributes(genes)) || !is.null(attributes(samples)) ||
    any(!nzchar(genes)) || any(!nzchar(samples)) || !length(x@x) ||
    !is.null(attributes(x@Dim)) || !is.null(attributes(x@p)) ||
    !is.null(attributes(x@i)) || !is.null(attributes(x@x)) ||
    !isTRUE(methods::validObject(x, test = TRUE))) {
    return(NULL)
  }
  n_stored <- x@p[[length(x@p)]]
  if (length(x@x) != n_stored || length(x@i) != n_stored) return(NULL)
  genes <- toupper(genes)
  if (anyDuplicated(genes)) return(NULL)
  common <- genes %in% sparse_common()
  n_common <- sum(common)
  if (n_common == 0L || n_common > 0.55 * length(genes)) return(NULL)

  values <- x@x
  nonzero <- is.finite(values) & values != 0
  keep <- logical(nrow(x))
  keep[x@i[nonzero] + 1L] <- TRUE
  n_common <- sum(common & keep)
  if (n_common == 0L || n_common > 0.55 * sum(keep)) return(NULL)

  if (!all(nonzero)) x@x[!nonzero] <- 0
  dimnames(x) <- list(genes, samples)
  if (all(keep)) x else x[keep, , drop = FALSE]
}

estimate_ssgsea_scores <- function(expr, gene_sets) {
  check_r("matrixStats", verbose = FALSE)
  ranked <- t(matrixStats::colRanks(
    expr,
    ties.method = "average"
  ))
  ranked <- as.matrix(ranked)
  ranked[!is.finite(ranked)] <- 0
  ranked <- 10000 * ranked / nrow(ranked)
  rownames(ranked) <- rownames(expr)
  colnames(ranked) <- colnames(expr)
  sample_order <- vapply(
    seq_len(ncol(ranked)),
    function(sample_i) order(ranked[, sample_i], decreasing = TRUE),
    integer(nrow(ranked))
  )
  gene_index <- seq_len(nrow(ranked))
  names(gene_index) <- rownames(ranked)
  gene_sets_index <- lapply(gene_sets, function(gene_set) {
    unname(gene_index[intersect(gene_set, names(gene_index))])
  })
  scores <- estimate_ssgsea_scores_cpp(
    ranked = ranked,
    sample_order = sample_order,
    gene_sets = gene_sets_index
  )
  dimnames(scores) <- list(colnames(ranked), names(gene_sets))
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

  expr <- GetAssayData5(object, assay = assay, layer = layer)
  meta_rows <- match(colnames(expr), rownames(meta))
  if (!ncol(expr)) {
    log_message("Selected Seurat assay layer has no cells.", message_type = "error")
  }
  if (length(meta_rows) != ncol(expr) || anyNA(meta_rows) || anyDuplicated(meta_rows)) {
    log_message(
      "Selected Seurat assay layer cells must match distinct Seurat metadata rows.",
      message_type = "error"
    )
  }

  sample_vec <- as.character(meta[[sample.by %||% group.by]][meta_rows])
  if (any(is.na(sample_vec) | !nzchar(sample_vec))) {
    log_message(
      "Pseudo-bulk sample labels contain missing or empty values.",
      message_type = "error"
    )
  }
  group_vec <- if (!is.null(group.by)) {
    as.character(meta[[group.by]][meta_rows])
  } else {
    sample_vec
  }
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

  if (!is.null(features)) {
    feature_key <- toupper(rownames(expr))
    features_key <- toupper(features)
    keep <- feature_key %in% features_key
    expr <- expr[keep, , drop = FALSE]
  }
  if (!inherits(expr, c("dgCMatrix", "matrix", "Matrix"))) {
    expr <- as_matrix(expr)
  }
  sample_factor <- factor(sample_vec, levels = unique(sample_vec[order(meta_rows)]))
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
