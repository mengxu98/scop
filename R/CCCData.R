#' @title Run LIANA cell-cell communication analysis
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams standard_scop
#' @param group.by Metadata column defining cell groups. Passed to
#' `liana::liana_wrap()` as `idents_col`.
#' @param method LIANA methods to run. Defaults to LIANA's internal methods.
#' @param resource LIANA ligand-receptor resource(s). Default is `"Consensus"`.
#' @param assay Assay used by LIANA. If `NULL`, LIANA uses the default assay.
#' @param min_cells Minimum cells per identity retained by LIANA.
#' @param return_all Whether LIANA should return all possible interactions.
#' @param ... Additional arguments passed to `liana::liana_wrap()`.
#'
#' @return A `Seurat` object with LIANA results stored in
#' `srt@tools[["LIANA"]]`.
#' @export
RunLIANA <- function(
  srt,
  group.by,
  method = c("natmi", "connectome", "logfc", "sca", "cellphonedb"),
  resource = "Consensus",
  assay = NULL,
  min_cells = 5,
  return_all = FALSE,
  verbose = TRUE,
  ...
) {
  check_r(c("saezlab/liana", "SingleCellExperiment"), verbose = FALSE)
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (!is.character(group.by) || length(group.by) != 1L || !group.by %in% colnames(srt[[]])) {
    log_message(
      "{.arg group.by} must be a valid metadata column in {.cls Seurat}",
      message_type = "error"
    )
  }

  log_message(
    "Running {.pkg LIANA} cell-cell communication analysis...",
    verbose = verbose
  )

  counts <- GetAssayData5(
    object = srt,
    assay = assay %||% SeuratObject::DefaultAssay(srt),
    layer = "counts"
  )
  logcounts <- GetAssayData5(
    object = srt,
    assay = assay %||% SeuratObject::DefaultAssay(srt),
    layer = "data"
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = counts,
      logcounts = logcounts
    ),
    colData = srt[[]]
  )

  dots <- list(...)
  if (is.null(dots$base)) {
    dots$base <- exp(1)
  }
  res <- do.call(
    getExportedValue("liana", "liana_wrap"),
    c(
      list(
        sce = sce,
        method = method,
        resource = resource,
        idents_col = group.by,
        assay = NULL,
        min_cells = min_cells,
        return_all = return_all,
        verbose = verbose
      ),
      dots
    )
  )

  long_table <- standardize_liana_result(res)
  liana_table <- ccc_long_to_liana(long_table)
  pair_table <- aggregate_ccc_long(long_table)

  srt@tools[["LIANA"]] <- list(
    method = "LIANA",
    results = res,
    long_table = long_table,
    liana_table = liana_table,
    pair_table = pair_table,
    parameters = list(
      group.by = group.by,
      method = method,
      resource = resource,
      assay = assay,
      min_cells = min_cells,
      return_all = return_all
    )
  )
  srt <- ccc_update_unified_bundle(
    srt = srt,
    method = "LIANA",
    bundle = srt@tools[["LIANA"]]
  )

  log_message(
    "{.pkg LIANA} analysis completed",
    message_type = "success",
    verbose = verbose
  )
  srt
}

#' @title Convert CCC results to a LIANA-like table
#'
#' @md
#' @inheritParams CCCHeatmap
#' @param aggregate Whether to aggregate duplicate
#' `source-target-ligand-receptor` rows. This should usually stay `TRUE` for
#' OmicVerse compatibility.
#' @param sample_col Optional column used as sample/context/dataset key. If
#' `NULL`, the first available column among `"sample"`, `"context"`,
#' `"condition"`, and `"dataset"` is used.
#' @param score_col Column used as the exported communication score.
#' @param pvalue_col Column used as the exported p-value/rank-like support.
#'
#' @return A data frame with LIANA/LIANA+ compatible columns:
#' `source`, `target`, `ligand_complex`, `receptor_complex`, `score`, and
#' `pvalue`, plus scop/OV-friendly metadata.
#' @export
ccc_to_liana <- function(
  srt,
  method = NULL,
  condition = NULL,
  dataset = 1,
  slot.name = "net",
  signaling = NULL,
  pairLR.use = NULL,
  sender.use = NULL,
  receiver.use = NULL,
  ligand.use = NULL,
  receptor.use = NULL,
  interaction.use = NULL,
  thresh = 0.05,
  aggregate = TRUE,
  sample_col = NULL,
  score_col = "score",
  pvalue_col = "pvalue"
) {
  df <- ccc_result_long_table(
    srt = srt,
    method = method,
    condition = condition,
    dataset = dataset,
    slot.name = slot.name,
    signaling = signaling,
    pairLR.use = pairLR.use,
    sender.use = sender.use,
    receiver.use = receiver.use,
    ligand.use = ligand.use,
    receptor.use = receptor.use,
    interaction.use = interaction.use,
    thresh = thresh
  )
  if (
    is.null(sample_col) &&
      "method" %in% colnames(df) &&
      length(unique(as.character(df$method))) > 1L
  ) {
    sample_col <- "method"
  }
  ccc_long_to_liana(
    df = df,
    aggregate = aggregate,
    sample_col = sample_col,
    score_col = score_col,
    pvalue_col = pvalue_col
  )
}

#' @title Convert CCC results to OmicVerse communication AnnData
#'
#' @md
#' @inheritParams ccc_to_liana
#' @param liana_res Optional precomputed LIANA-like data frame. If supplied,
#' `srt` is not required.
#' @param score_key Column in `liana_res` used for `layers[["means"]]`.
#' @param pvalue_key Column in `liana_res` used for `layers[["pvalues"]]`.
#' @param inverse_score Whether smaller `score_key` values should be converted
#' to larger communication strengths. Useful for rank-like metrics.
#' @param inverse_pvalue Whether smaller `pvalue_key` values should be inverted
#' before writing to `layers[["pvalues"]]`.
#' @param verbose Whether to print progress messages.
#' @param h5ad_path Optional output path. If provided, the AnnData object is
#' written to this file before being returned.
#'
#' @return A Python `anndata.AnnData` communication object compatible with
#' `ov.pl.ccc_*` plotting functions.
#' @export
ccc_to_adata <- function(
  srt = NULL,
  method = NULL,
  condition = NULL,
  dataset = 1,
  slot.name = "net",
  signaling = NULL,
  pairLR.use = NULL,
  sender.use = NULL,
  receiver.use = NULL,
  ligand.use = NULL,
  receptor.use = NULL,
  interaction.use = NULL,
  thresh = 0.05,
  liana_res = NULL,
  score_key = "score",
  pvalue_key = "pvalue",
  inverse_score = FALSE,
  inverse_pvalue = FALSE,
  sample_col = NULL,
  h5ad_path = NULL,
  verbose = TRUE
) {
  if (is.null(liana_res)) {
    if (is.null(srt)) {
      log_message(
        "Provide either {.arg srt} or {.arg liana_res}",
        message_type = "error"
      )
    }
    liana_res <- ccc_to_liana(
      srt = srt,
      method = method,
      condition = condition,
      dataset = dataset,
      slot.name = slot.name,
      signaling = signaling,
      pairLR.use = pairLR.use,
      sender.use = sender.use,
      receiver.use = receiver.use,
      ligand.use = ligand.use,
      receptor.use = receptor.use,
      interaction.use = interaction.use,
      thresh = thresh,
      sample_col = sample_col
    )
  } else {
    liana_res <- ccc_long_to_liana(
      df = liana_res,
      aggregate = TRUE,
      sample_col = sample_col,
      score_col = score_key,
      pvalue_col = pvalue_key
    )
    score_key <- "score"
    pvalue_key <- "pvalue"
  }

  if (!all(c("source", "target", "ligand_complex", "receptor_complex") %in% colnames(liana_res))) {
    log_message(
      "{.arg liana_res} must contain source, target, ligand_complex, and receptor_complex columns",
      message_type = "error"
    )
  }
  if (!score_key %in% colnames(liana_res)) {
    log_message(
      "{.arg score_key} ({.val {score_key}}) is not present in {.arg liana_res}",
      message_type = "error"
    )
  }
  if (!pvalue_key %in% colnames(liana_res)) {
    log_message(
      "{.arg pvalue_key} ({.val {pvalue_key}}) is not present in {.arg liana_res}",
      message_type = "error"
    )
  }

  if (!isTRUE(getOption("scop_skip_python_prepare", FALSE))) {
    PrepareEnv(modules = "scanpy")
    check_python(c("anndata", "numpy"), verbose = FALSE)
  }
  ad <- reticulate::import("anndata", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)

  sample_col <- ccc_resolve_sample_col(liana_res, sample_col)
  row_id <- paste(liana_res$source, liana_res$target, sep = "|")
  if (!is.null(sample_col)) {
    row_id <- paste0(row_id, "|", sample_col, "=", liana_res[[sample_col]])
  }
  col_id <- paste(liana_res$ligand_complex, liana_res$receptor_complex, sep = " -> ")

  obs_levels <- sort(unique(row_id))
  var_levels <- sort(unique(col_id))
  score_matrix <- matrix(
    0,
    nrow = length(obs_levels),
    ncol = length(var_levels),
    dimnames = list(obs_levels, var_levels)
  )
  pvalue_matrix <- matrix(
    1,
    nrow = length(obs_levels),
    ncol = length(var_levels),
    dimnames = list(obs_levels, var_levels)
  )

  score_values <- ccc_prepare_metric_values(
    liana_res[[score_key]],
    invert = inverse_score,
    fill = 0
  )
  pvalue_values <- ccc_prepare_metric_values(
    liana_res[[pvalue_key]],
    invert = inverse_pvalue,
    fill = 1
  )
  idx <- cbind(match(row_id, obs_levels), match(col_id, var_levels))
  score_matrix[idx] <- score_values
  pvalue_matrix[idx] <- pvalue_values

  obs <- unique(data.frame(
    row_id = row_id,
    sender = as.character(liana_res$source),
    receiver = as.character(liana_res$target),
    cell_type_pair = row_id,
    stringsAsFactors = FALSE
  ))
  if (!is.null(sample_col)) {
    obs[[sample_col]] <- as.character(liana_res[[sample_col]][match(obs$row_id, row_id)])
  }
  rownames(obs) <- obs$row_id
  obs <- obs[obs_levels, setdiff(colnames(obs), "row_id"), drop = FALSE]
  obs[] <- lapply(obs, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
  })

  var <- unique(data.frame(
    var_id = col_id,
    interacting_pair = as.character(liana_res$interacting_pair),
    pair_lr = as.character(liana_res$pair_lr),
    interaction_name = as.character(liana_res$interaction_name),
    interaction_name_2 = as.character(liana_res$interaction_name_2),
    classification = as.character(liana_res$classification),
    pathway_name = as.character(liana_res$classification),
    signaling = as.character(liana_res$classification),
    gene_a = as.character(liana_res$ligand_complex),
    gene_b = as.character(liana_res$receptor_complex),
    ligand = as.character(liana_res$ligand_complex),
    receptor = as.character(liana_res$receptor_complex),
    annotation_strategy = "scop_ccc_to_liana",
    classification_source = "scop",
    stringsAsFactors = FALSE
  ))
  var <- var[!duplicated(var$var_id), , drop = FALSE]
  rownames(var) <- var$var_id
  var <- var[var_levels, setdiff(colnames(var), "var_id"), drop = FALSE]
  var[] <- lapply(var, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
  })

  adata <- ad$AnnData(
    X = reticulate::np_array(score_matrix, dtype = np$float32),
    obs = obs,
    var = var
  )
  adata$obs_names <- reticulate::r_to_py(as.list(rownames(obs)))
  adata$var_names <- reticulate::r_to_py(as.list(rownames(var)))
  adata$layers$`__setitem__`(
    "means",
    reticulate::np_array(score_matrix, dtype = np$float32)
  )
  adata$layers$`__setitem__`(
    "pvalues",
    reticulate::np_array(pvalue_matrix, dtype = np$float32)
  )
  uns_liana_res <- liana_res
  uns_liana_res[] <- lapply(uns_liana_res, function(x) {
    if (is.factor(x)) {
      x <- as.character(x)
    }
    if (is.character(x)) {
      x[is.na(x)] <- ""
      return(x)
    }
    if (is.logical(x)) {
      x[is.na(x)] <- FALSE
      return(x)
    }
    if (is.numeric(x) || is.integer(x)) {
      return(x)
    }
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
  })
  adata$uns$`__setitem__`("liana_res", reticulate::r_to_py(uns_liana_res))
  adata$uns$`__setitem__`("liana_score_key", score_key)
  adata$uns$`__setitem__`("liana_pvalue_key", pvalue_key)
  adata$uns$`__setitem__`("liana_uns_key", "liana_res")
  adata$uns$`__setitem__`("liana_sample_key", sample_col %||% "none")
  adata$uns$`__setitem__`("liana_classification_reference", "scop")
  adata$uns$`__setitem__`("liana_classification_fallback", "Unclassified")

  if (!is.null(h5ad_path)) {
    h5ad_path <- normalizePath(h5ad_path, mustWork = FALSE)
    adata$write_h5ad(h5ad_path)
    log_message(
      "Wrote CCC communication AnnData to {.file {h5ad_path}}",
      verbose = verbose
    )
  }

  adata
}

ccc_result_long_table <- function(
  srt,
  method = NULL,
  condition = NULL,
  dataset = 1,
  slot.name = "net",
  signaling = NULL,
  pairLR.use = NULL,
  sender.use = NULL,
  receiver.use = NULL,
  ligand.use = NULL,
  receptor.use = NULL,
  interaction.use = NULL,
  thresh = 0.05
) {
  method <- detect_method(srt = srt, method = method)
  df <- ccc_long_table_for_method(
    srt = srt,
    method = method,
    condition = condition,
    dataset = dataset,
    slot.name = slot.name,
    signaling = signaling,
    pairLR.use = pairLR.use,
    sources.use = sender.use,
    targets.use = receiver.use,
    thresh = thresh
  )
  df <- standardize_long_df(df)
  if (nrow(df) == 0L) {
    return(df)
  }
  df <- filter_long_df(
    df = df,
    sender.use = sender.use,
    receiver.use = receiver.use,
    ligand.use = ligand.use,
    receptor.use = receptor.use,
    interaction.use = interaction.use,
    signaling = signaling,
    pairLR.use = pairLR.use
  )
  if (nrow(df) == 0L) {
    return(df)
  }
  if (!"method" %in% colnames(df)) {
    df$method <- method
  }
  df <- ccc_mark_significance(df, thresh = thresh)
  df
}

ccc_long_to_liana <- function(
  df,
  aggregate = TRUE,
  sample_col = NULL,
  score_col = "score",
  pvalue_col = "pvalue"
) {
  df <- standardize_long_df(df)
  if (is.null(df) || nrow(df) == 0L) {
    return(data.frame())
  }

  ligand <- as.character(df$ligand)
  receptor <- as.character(df$receptor)
  missing_lr <- is.na(ligand) | !nzchar(ligand) | is.na(receptor) | !nzchar(receptor)
  if (any(missing_lr)) {
    parsed <- ccc_parse_lr_label(df$pair_lr[missing_lr])
    empty_parsed <- !nzchar(parsed$ligand) | !nzchar(parsed$receptor)
    if (any(empty_parsed)) {
      parsed2 <- ccc_parse_lr_label(df$interaction_label[missing_lr][empty_parsed])
      parsed$ligand[empty_parsed] <- parsed2$ligand
      parsed$receptor[empty_parsed] <- parsed2$receptor
    }
    ligand[missing_lr] <- ifelse(nzchar(parsed$ligand), parsed$ligand, ligand[missing_lr])
    receptor[missing_lr] <- ifelse(nzchar(parsed$receptor), parsed$receptor, receptor[missing_lr])
  }

  sample_col <- ccc_resolve_sample_col(df, sample_col)
  score <- ccc_numeric_column(df, score_col)
  if (all(!is.finite(score))) {
    score <- ccc_numeric_column(
      df,
      c(
        "means",
        "prob",
        "interaction_score",
        "prioritization_score",
        "LRscore",
        "lrscore",
        "magnitude",
        "weight",
        "specificity"
      )
    )
  }
  pvalue <- ccc_numeric_column(df, pvalue_col)
  if (all(!is.finite(pvalue))) {
    pvalue <- ccc_numeric_column(
      df,
      c(
        "pvalue",
        "pval",
        "cellphone_pvals",
        "pvalues",
        "aggregate_rank",
        "specificity_rank",
        "magnitude_rank",
        "cellphonedb.pvalue"
      )
    )
  }
  if (all(!is.finite(score)) && any(is.finite(pvalue))) {
    score <- ccc_score_from_rank_like(pvalue)
  }
  score[!is.finite(score)] <- 0
  pvalue[!is.finite(pvalue)] <- 1

  out <- data.frame(
    source = as.character(df$sender),
    target = as.character(df$receiver),
    ligand_complex = ccc_clean_identifier(ligand, drop_receptor_suffix = FALSE),
    receptor_complex = ccc_clean_identifier(receptor, drop_receptor_suffix = FALSE),
    score = score,
    pvalue = pvalue,
    stringsAsFactors = FALSE
  )
  out$ligand_complex[is.na(out$ligand_complex)] <- ""
  out$receptor_complex[is.na(out$receptor_complex)] <- ""
  out$interaction_name_2 <- paste(out$ligand_complex, out$receptor_complex, sep = " - ")
  out$interacting_pair <- paste(out$ligand_complex, out$receptor_complex, sep = "_")
  out$pair_lr <- paste(out$ligand_complex, out$receptor_complex, sep = "-")
  out$interaction_name <- if ("interaction_name" %in% colnames(df)) {
    as.character(df$interaction_name)
  } else {
    out$interacting_pair
  }
  out$classification <- as.character(df$classification %||% df$pathway_name)
  out$classification[is.na(out$classification) | !nzchar(out$classification)] <- "Unclassified"
  out$method <- if ("method" %in% colnames(df)) as.character(df$method) else NA_character_
  if ("liana_method" %in% colnames(df)) {
    out$liana_method <- as.character(df$liana_method)
  }
  if ("resource" %in% colnames(df)) {
    out$resource <- as.character(df$resource)
  }
  if (!is.null(sample_col)) {
    out[[sample_col]] <- as.character(df[[sample_col]])
  }

  keep <- !is.na(out$source) & nzchar(out$source) &
    !is.na(out$target) & nzchar(out$target) &
    !is.na(out$ligand_complex) & nzchar(out$ligand_complex) &
    !is.na(out$receptor_complex) & nzchar(out$receptor_complex)
  out <- out[keep, , drop = FALSE]
  if (nrow(out) == 0L) {
    return(out)
  }

  if (isTRUE(aggregate)) {
    out <- ccc_aggregate_liana_table(out, sample_col = sample_col)
  }
  out$specificity_rank <- rank(-out$score, ties.method = "average", na.last = "keep")
  out$magnitude_rank <- out$specificity_rank
  out$aggregate_rank <- ifelse(
    is.finite(out$pvalue),
    out$pvalue,
    out$specificity_rank
  )
  rownames(out) <- NULL
  out
}

standardize_liana_result <- function(res) {
  raw <- ccc_flatten_liana_result(res)
  if (nrow(raw) == 0L) {
    return(data.frame())
  }
  if ("method" %in% colnames(raw)) {
    raw$liana_method <- as.character(raw$method)
  } else if (!"liana_method" %in% colnames(raw)) {
    raw$liana_method <- NA_character_
  }
  raw$method <- "LIANA"
  out <- standardize_long_df(raw)

  score <- ccc_numeric_column(
    out,
    c(
      "score",
      "sca.LRscore",
      "LRscore",
      "lrscore",
      "magnitude",
      "lr_means",
      "lr.mean",
      "lr_mean",
      "expr_prod",
      "ligand.expr",
      "receptor.expr",
      "weight",
      "specificity"
    )
  )
  pvalue <- ccc_numeric_column(
    out,
    c(
      "pvalue",
      "pval",
      "cellphone_pvals",
      "cellphonedb.pvalue",
      "aggregate_rank",
      "specificity_rank",
      "magnitude_rank"
    )
  )
  if (all(!is.finite(score)) && any(is.finite(pvalue))) {
    score <- ccc_score_from_rank_like(pvalue)
  }
  out$score <- score
  out$pvalue <- pvalue
  out$score[!is.finite(out$score)] <- 0
  out$pvalue[!is.finite(out$pvalue)] <- 1
  out <- ccc_mark_significance(out)
  out
}

ccc_flatten_liana_result <- function(x, path = character()) {
  if (is.null(x)) {
    return(data.frame())
  }
  if (inherits(x, "python.builtin.object")) {
    x <- py_to_r2(x)
  }
  if (is.data.frame(x)) {
    out <- standardize_df(x)
    out$liana_method <- path[1] %||% NA_character_
    out$resource <- path[2] %||% NA_character_
    return(out)
  }
  if (!is.list(x)) {
    return(data.frame())
  }
  nms <- names(x)
  pieces <- lapply(seq_along(x), function(i) {
    nm <- nms[i] %||% paste0("result_", i)
    if (is.na(nm) || !nzchar(nm)) {
      nm <- paste0("result_", i)
    }
    ccc_flatten_liana_result(x[[i]], path = c(path, nm))
  })
  pieces <- Filter(function(el) is.data.frame(el) && nrow(el) > 0L, pieces)
  if (length(pieces) == 0L) {
    return(data.frame())
  }
  common <- Reduce(union, lapply(pieces, colnames))
  pieces <- lapply(pieces, function(el) {
    missing_cols <- setdiff(common, colnames(el))
    for (nm in missing_cols) {
      el[[nm]] <- NA
    }
    el[, common, drop = FALSE]
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

ccc_parse_lr_label <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  parse_one <- function(value) {
    value <- trimws(value)
    if (!nzchar(value)) {
      return(c(ligand = "", receptor = ""))
    }
    for (pattern in c("\\s+-\\s+", "\\s+->\\s+", "\\|", "_", "-")) {
      parts <- strsplit(value, pattern, perl = TRUE)[[1]]
      parts <- trimws(parts)
      parts <- parts[nzchar(parts)]
      if (length(parts) >= 2L) {
        return(c(ligand = parts[1], receptor = paste(parts[-1], collapse = "_")))
      }
    }
    c(ligand = "", receptor = "")
  }
  parsed <- t(vapply(x, parse_one, character(2)))
  data.frame(
    ligand = parsed[, "ligand"],
    receptor = parsed[, "receptor"],
    stringsAsFactors = FALSE
  )
}

ccc_resolve_sample_col <- function(df, sample_col = NULL) {
  if (!is.null(sample_col)) {
    if (!sample_col %in% colnames(df)) {
      log_message(
        "{.arg sample_col} ({.val {sample_col}}) is not present in the CCC table",
        message_type = "error"
      )
    }
    return(sample_col)
  }
  candidates <- c("sample", "context", "condition", "dataset")
  hit <- candidates[candidates %in% colnames(df)][1]
  if (is.na(hit)) {
    NULL
  } else {
    hit
  }
}

ccc_numeric_column <- function(df, candidates) {
  for (candidate in candidates) {
    col <- ccc_pick_col(df, candidate)
    if (is.null(col)) {
      next
    }
    values <- suppressWarnings(as.numeric(df[[col]]))
    if (any(is.finite(values))) {
      return(values)
    }
  }
  rep(NA_real_, nrow(df))
}

ccc_score_from_rank_like <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(NA_real_, length(x))
  finite <- is.finite(x)
  if (!any(finite)) {
    return(out)
  }
  x_finite <- x[finite]
  if (all(x_finite >= 0 & x_finite <= 1)) {
    out[finite] <- 1 - x_finite
  } else {
    out[finite] <- 1 / (x_finite + 1e-12)
  }
  out
}

ccc_prepare_metric_values <- function(x, invert = FALSE, fill = 0) {
  x <- suppressWarnings(as.numeric(x))
  if (isTRUE(invert)) {
    x <- ccc_score_from_rank_like(x)
  }
  x[!is.finite(x)] <- fill
  x
}

ccc_aggregate_liana_table <- function(out, sample_col = NULL) {
  key_cols <- c("source", "target", "ligand_complex", "receptor_complex")
  if (!is.null(sample_col) && sample_col %in% colnames(out)) {
    key_cols <- c(key_cols, sample_col)
  }
  key <- do.call(paste, c(out[key_cols], sep = "\r"))
  split_idx <- split(seq_len(nrow(out)), key, drop = TRUE)
  pieces <- lapply(split_idx, function(idx) {
    x <- out[idx, , drop = FALSE]
    row <- x[1, , drop = FALSE]
    row$score <- sum(suppressWarnings(as.numeric(x$score)), na.rm = TRUE)
    pvalue <- suppressWarnings(as.numeric(x$pvalue))
    pvalue <- pvalue[is.finite(pvalue)]
    row$pvalue <- if (length(pvalue) > 0L) min(pvalue, na.rm = TRUE) else 1
    if ("classification" %in% colnames(row)) {
      row$classification <- ccc_first_nonempty(x$classification, "Unclassified")
    }
    if ("method" %in% colnames(row)) {
      row$method <- paste(unique(stats::na.omit(as.character(x$method))), collapse = ";")
    }
    if ("liana_method" %in% colnames(row)) {
      row$liana_method <- paste(unique(stats::na.omit(as.character(x$liana_method))), collapse = ";")
    }
    if ("resource" %in% colnames(row)) {
      row$resource <- paste(unique(stats::na.omit(as.character(x$resource))), collapse = ";")
    }
    row
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

ccc_first_nonempty <- function(x, default = NA_character_) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) {
    default
  } else {
    x[1]
  }
}
