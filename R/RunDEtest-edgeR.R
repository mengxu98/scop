# Explicit models are deliberately limited to the existing edgeR sample path.
RunDEtestCheckDesign <- function(design, contrast, test.use = "edgeR") {
  explicit <- !is.null(design) || !is.null(contrast)
  if (!explicit) {
    return(FALSE)
  }
  if (!identical(test.use, "edgeR")) {
    log_message("design and contrast are supported only for sample-level edgeR testing.", message_type = "error")
  }
  if (!inherits(design, "formula") || length(design) != 2L ||
    "." %in% all.vars(design)) {
    log_message("design must be a one-sided formula with explicit sample metadata columns.", message_type = "error")
  }
  uses_offset <- function(term) {
    if (!is.call(term)) {
      return(FALSE)
    }
    head <- term[[1L]]
    if (identical(head, as.name("offset")) ||
      (is.call(head) && length(head) == 3L &&
        as.character(head[[1L]]) %in% c("::", ":::") &&
        identical(head[[3L]], as.name("offset")))) {
      return(TRUE)
    }
    any(vapply(as.list(term)[-1L], uses_offset, logical(1)))
  }
  if (uses_offset(design)) {
    log_message("offset() terms are not supported in explicit edgeR designs; offsets would not be represented by the fitted model matrix.", message_type = "error")
  }
  if (!is.numeric(contrast) || !is.null(dim(contrast)) || !length(contrast) ||
    any(!is.finite(contrast)) || !any(contrast != 0) ||
    is.null(names(contrast)) || anyNA(names(contrast)) ||
    any(!nzchar(names(contrast))) || anyDuplicated(names(contrast))) {
    log_message("contrast must be a finite, nonzero named numeric vector matching every design column.", message_type = "error")
  }
  TRUE
}

RunDEtestEdgeRBackend <- function() {
  available <- check_r("edgeR", install = FALSE, verbose = FALSE)
  if (!isTRUE(all(unlist(available, use.names = FALSE)))) {
    log_message("edgeR is required for sample-level edgeR testing; install it before running RunDEtest().", message_type = "error")
  }
  invisible(TRUE)
}

RunDEtestModelSampleData <- function(metadata, sample_col, condition_col, design) {
  columns <- unique(c(sample_col, condition_col, all.vars(design)))
  missing <- setdiff(columns, colnames(metadata))
  if (length(missing)) {
    log_message("{paste(\"Missing design metadata columns:\", paste(missing, collapse = \", \"))}", message_type = "error")
  }
  metadata <- metadata[, columns, drop = FALSE]
  if (anyNA(metadata) || any(!nzchar(as.character(metadata[[sample_col]])))) {
    log_message("Explicit design metadata and sample IDs must not contain missing values or empty sample IDs.", message_type = "error")
  }
  for (column in setdiff(columns, sample_col)) {
    varying <- vapply(
      split(metadata[[column]], as.character(metadata[[sample_col]])),
      function(x) length(unique(x)) != 1L, logical(1)
    )
    if (any(varying)) {
      log_message("{paste0(\n        \"Design metadata '\", column, \"' must be constant within each aggregate sample. \",\n        \"For paired cell types, use distinct donor-by-cell-type sample IDs and include donor in design.\"\n      )}", message_type = "error")
    }
  }
  out <- metadata[!duplicated(as.character(metadata[[sample_col]])), , drop = FALSE]
  rownames(out) <- as.character(out[[sample_col]])
  out
}

RunDEtestEdgeRModel <- function(
  count_matrix, condition, condition1, condition2,
  sample_data, design, contrast, only.pos, logfc.threshold, p.adjust.method
) {
  RunDEtestCheckDesign(design, contrast)
  ids <- colnames(count_matrix)
  features <- rownames(count_matrix)
  valid_ids <- function(x, n) {
    !is.null(x) && length(x) == n &&
      !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x)
  }
  if (!valid_ids(ids, ncol(count_matrix)) || !valid_ids(features, nrow(count_matrix))) {
    log_message("Explicit edgeR models require unique, non-empty sample and feature IDs.", message_type = "error")
  }
  if (length(condition) != length(ids) || anyNA(condition)) {
    log_message("condition must contain one non-missing label per count-matrix sample.", message_type = "error")
  }
  pair <- resolve_condition_pair(condition, condition1, condition2, strict_two_levels = TRUE)
  keep <- as.character(condition) %in% unlist(pair, use.names = FALSE)
  condition <- factor(as.character(condition)[keep], levels = unlist(pair, use.names = FALSE))
  count_matrix <- as.matrix(count_matrix[, keep, drop = FALSE])
  ids <- colnames(count_matrix)
  if (any(table(condition) < 2L)) {
    log_message("Explicit edgeR models require at least 2 biological samples per condition; cells or technical libraries are not replicates.", message_type = "error")
  }
  if (any(!is.finite(count_matrix)) || any(count_matrix < 0) || any(count_matrix %% 1 != 0)) {
    log_message("Explicit edgeR models require finite, non-negative integer raw counts.", message_type = "error")
  }
  if (any(colSums(count_matrix) == 0)) {
    log_message("Explicit edgeR models cannot fit a sample with zero total counts.", message_type = "error")
  }
  if (!is.data.frame(sample_data) || !valid_ids(rownames(sample_data), nrow(sample_data)) ||
    !all(ids %in% rownames(sample_data))) {
    log_message("sample_data must have unique row names matching count-matrix sample IDs.", message_type = "error")
  }
  sample_data <- droplevels(sample_data[ids, , drop = FALSE])
  missing <- setdiff(all.vars(design), colnames(sample_data))
  if (length(missing)) {
    log_message("{paste(\"Missing design metadata columns:\", paste(missing, collapse = \", \"))}", message_type = "error")
  }
  model_frame <- stats::model.frame(design, data = sample_data, na.action = stats::na.fail)
  model_matrix <- stats::model.matrix(design, data = model_frame)
  if (nrow(model_matrix) != ncol(count_matrix) || any(!is.finite(model_matrix))) {
    log_message("design must produce one finite model-matrix row per sample.", message_type = "error")
  }
  rank <- qr(model_matrix)$rank
  if (rank < ncol(model_matrix)) {
    log_message("design is not full rank: remove confounded or redundant covariates before testing contrasts.", message_type = "error")
  }
  residual_df <- nrow(model_matrix) - rank
  if (residual_df < 1L) {
    log_message("design has no residual degrees of freedom; biological replication is required.", message_type = "error")
  }
  if (length(contrast) != ncol(model_matrix) || !setequal(names(contrast), colnames(model_matrix))) {
    log_message("{paste0(\n      \"contrast names must match every design column exactly: \",\n      paste(colnames(model_matrix), collapse = \", \")\n    )}", message_type = "error")
  }
  contrast <- contrast[colnames(model_matrix)]
  # Full column rank guarantees that every finite coefficient contrast is estimable.
  direction <- paste(paste0(
    format(contrast[contrast != 0], trim = TRUE), " * ",
    names(contrast)[contrast != 0]
  ), collapse = " + ")
  RunDEtestEdgeRBackend()
  if (!requireNamespace("statmod", quietly = TRUE)) {
    log_message("statmod is required for robust edgeR quasi-likelihood fitting; install it before running RunDEtest().", message_type = "error")
  }
  dge <- edgeR::DGEList(counts = count_matrix)
  keep_features <- edgeR::filterByExpr(dge, design = model_matrix)
  full <- data.frame(
    p_val = numeric(), avg_log2FC = numeric(), ave_expr = numeric(),
    pct.1 = numeric(), pct.2 = numeric(), F = numeric(), p_val_adj = numeric(),
    p_val_BH = numeric(), row.names = character()
  )
  normalized_samples <- NULL
  if (any(keep_features)) {
    dge <- dge[keep_features, , keep.lib.sizes = FALSE]
    if (any(colSums(dge$counts) == 0)) {
      log_message("A sample has zero library size after filterByExpr; the explicit model cannot be fitted.", message_type = "error")
    }
    dge <- edgeR::calcNormFactors(dge)
    dge <- edgeR::estimateDisp(dge, design = model_matrix)
    fit <- edgeR::glmQLFit(dge, design = model_matrix, robust = TRUE)
    qlf <- edgeR::glmQLFTest(fit, contrast = unname(contrast))
    tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
    detect <- dge$counts[rownames(tt), , drop = FALSE] > 0
    full <- data.frame(
      p_val = tt$PValue, avg_log2FC = tt$logFC, ave_expr = tt$logCPM,
      pct.1 = round(rowMeans(detect[, condition == pair$condition1, drop = FALSE]), 3),
      pct.2 = round(rowMeans(detect[, condition == pair$condition2, drop = FALSE]), 3),
      F = tt$F, row.names = rownames(tt)
    )
    # Both adjustments use the entire filterByExpr-tested family, before display filtering.
    full$p_val_adj <- stats::p.adjust(full$p_val, method = p.adjust.method)
    full$p_val_BH <- stats::p.adjust(full$p_val, method = "BH")
    normalized_samples <- dge$samples
  }
  out <- full[if (isTRUE(only.pos)) {
    full$avg_log2FC >= logfc.threshold
  } else {
    abs(full$avg_log2FC) >= logfc.threshold
  }, , drop = FALSE]
  attr(out, "edgeR_model") <- list(
    formula = paste(deparse(design), collapse = " "), design = model_matrix,
    contrast = contrast, positive_direction = direction, rank = rank, residual_df = residual_df,
    sample_data = sample_data, condition_pair = pair, sample_counts = table(condition),
    sample_unit = list(
      source = "user-declared count-matrix columns",
      biological_independence = "not verified from metadata"
    ),
    counts = count_matrix, library_sizes = colSums(count_matrix), normalized_samples = normalized_samples,
    input_features = features, tested_features = features[keep_features],
    filtered_features = features[!keep_features], full_results = full,
    p_adjust_method = p.adjust.method, bh_scope = "all filterByExpr-tested features per model",
    display_filter = list(only.pos = only.pos, logfc.threshold = logfc.threshold),
    backend = c(
      filter = "edgeR::filterByExpr", normalization = "edgeR::calcNormFactors",
      dispersion = "edgeR::estimateDisp", fit = "edgeR::glmQLFit", test = "edgeR::glmQLFTest"
    ),
    robust = TRUE, versions = c(R = as.character(getRversion()), edgeR = as.character(utils::packageVersion("edgeR")))
  )
  out
}

RunDEtestFeatureCoordinates <- function(object, assay, features, feature_type) {
  if (identical(feature_type, "gene")) {
    return(NULL)
  }
  if (inherits(object, "Seurat") && inherits(object[[assay]], "ChromatinAssay")) {
    ranges <- Signac::granges(object[[assay]])
    ids <- rownames(object[[assay]])
    if (!identical(atac_range_ids(ranges), ids)) {
      log_message("ChromatinAssay ranges must align exactly with peak-count rows before differential testing.", message_type = "error")
    }
  } else if (methods::is(object, "RangedSummarizedExperiment")) {
    ranges <- SummarizedExperiment::rowRanges(object)
    ids <- rownames(object)
    if (!inherits(ranges, "GRanges")) {
      log_message("Peak differential testing requires one GRanges interval per feature.", message_type = "error")
    }
  } else {
    log_message("Explicit peak models require a ChromatinAssay or RangedSummarizedExperiment with genomic coordinates.", message_type = "error")
  }
  ranges <- ranges[match(features, ids)]
  chromosomes <- as.character(GenomicRanges::seqnames(ranges))
  # GenomeInfoDb stores genome labels per seqlevel, not per peak. Resolve each
  # selected range by chromosome rather than recycling labels over count rows.
  genomes <- atac_backend_fun("GenomeInfoDb", "genome")(ranges)
  data.frame(
    feature = features, chromosome = chromosomes,
    start = GenomicRanges::start(ranges), end = GenomicRanges::end(ranges),
    genome = unname(as.character(genomes[chromosomes])),
    coordinate_system = "1-based closed", row.names = features, stringsAsFactors = FALSE
  )
}

RunDEtestAnnotateModel <- function(result, coordinates = NULL, feature_type = "gene") {
  result$gene <- rownames(result)
  result$feature <- rownames(result)
  result$feature_type <- rep(feature_type, nrow(result))
  if (!is.null(coordinates)) {
    result <- cbind(result, coordinates[rownames(result), setdiff(colnames(coordinates), "feature"), drop = FALSE])
  }
  result
}
