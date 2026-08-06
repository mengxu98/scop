# CCC method dispatch --------------------------------------------------------

normalize_ccc_method <- function(method) {
  alias_map <- c(
    "CCC" = "CCC",
    "CellPhoneDB" = "CellphoneDB",
    "CellphoneDB" = "CellphoneDB",
    "Liana" = "LIANA",
    "liana" = "LIANA",
    "NicheNet" = "Nichenetr",
    "MultiNicheNet" = "MultiNichenetr",
    "SpatialCellChat" = "SpatialCellChat",
    "MDIC3" = "MDIC3"
  )
  alias_map_lower <- c(
    "ccc" = "CCC",
    "unified" = "CCC",
    "cellchat" = "CellChat",
    "spatialcellchat" = "SpatialCellChat",
    "spatial_cellchat" = "SpatialCellChat",
    "spatial cellchat" = "SpatialCellChat",
    "cellphonedb" = "CellphoneDB",
    "cellphone_db" = "CellphoneDB",
    "cellphone db" = "CellphoneDB",
    "liana" = "LIANA",
    "nichenet" = "Nichenetr",
    "nichenetr" = "Nichenetr",
    "multinichenet" = "MultiNichenetr",
    "multinichenetr" = "MultiNichenetr",
    "mdic3" = "MDIC3"
  )
  method_chr <- as.character(method)[1]
  if (is.na(method_chr) || !nzchar(method_chr)) {
    log_message(
      "{.arg method} must be a non-empty string or NULL to auto-detect",
      message_type = "error"
    )
  }
  method_chr <- trimws(method_chr)
  if (method_chr %in% names(alias_map)) {
    return(unname(alias_map[[method_chr]]))
  }
  key <- tolower(method_chr)
  if (key %in% names(alias_map_lower)) {
    return(unname(alias_map_lower[[key]]))
  }
  method_chr
}

ccc_method_runner <- function(method) {
  method <- normalize_ccc_method(method)
  switch(method,
    CellChat = RunCellChat,
    CellphoneDB = RunCellphoneDB,
    LIANA = RunLIANA,
    Nichenetr = RunNichenetr,
    MultiNichenetr = RunMultiNichenetr,
    SpatialCellChat = RunSpatialCellChat,
    MDIC3 = RunMDIC3,
    ccc_unsupported_method(method)
  )
}

ccc_unsupported_method <- function(method) {
  log_message(
    "Unsupported CCC method {.val {method}}",
    message_type = "error"
  )
}

ccc_semantic_long_table <- function(df, method = NULL) {
  out <- standardize_long_df(df)
  if (nrow(out) == 0L) {
    return(out)
  }
  if (!is.null(method)) out$method <- method
  if (!"method" %in% colnames(out)) out$method <- NA_character_
  for (nm in c("resource", "condition", "sample")) {
    if (!nm %in% colnames(out)) out[[nm]] <- NA_character_
  }
  for (nm in c("producer", "backend_version")) {
    if (!nm %in% colnames(out)) out[[nm]] <- NA_character_
  }
  if (!"score_type" %in% colnames(out)) out$score_type <- "backend_score"
  if (!"score_direction" %in% colnames(out)) out$score_direction <- "higher_better"
  if (!"pvalue_type" %in% colnames(out)) {
    out$pvalue_type <- ifelse(is.finite(suppressWarnings(as.numeric(out$pvalue))),
      "backend_support", "not_available"
    )
  }
  if (!"support_type" %in% colnames(out)) out$support_type <- "method_specific"
  if (!"priority_rank" %in% colnames(out)) out$priority_rank <- NA_real_
  if (!"priority_score" %in% colnames(out)) out$priority_score <- NA_real_

  method_values <- unique(as.character(out$method))
  for (method_value in method_values) {
    idx <- which(as.character(out$method) == method_value)
    score <- suppressWarnings(as.numeric(out$score[idx]))
    finite <- is.finite(score)
    if (!any(finite)) next
    direction <- unique(as.character(out$score_direction[idx]))
    decreasing <- !length(direction) || !identical(direction[1], "lower_better")
    rank_value <- rep(NA_real_, length(idx))
    rank_value[finite] <- rank(if (decreasing) -score[finite] else score[finite],
      ties.method = "average"
    )
    percentile <- rank_value / sum(finite)
    missing_rank <- !is.finite(out$priority_rank[idx])
    out$priority_rank[idx[missing_rank]] <- percentile[missing_rank]
    missing_score <- !is.finite(out$priority_score[idx])
    out$priority_score[idx[missing_score]] <- 1 - percentile[missing_score]
  }
  out
}

ccc_filter_table_context <- function(
  df,
  resource = NULL,
  condition = NULL,
  sample = NULL
) {
  if (!is.data.frame(df) || nrow(df) == 0L) {
    return(df)
  }
  if (!is.null(resource)) {
    available_resource <- if ("resource" %in% colnames(df)) {
      unique(as.character(df$resource))
    } else {
      character(0)
    }
    available_resource <- available_resource[
      !is.na(available_resource) & nzchar(available_resource)
    ]
    if (length(available_resource) == 0L) {
      log_message(
        "The selected CCC result does not provide {.arg resource} provenance",
        message_type = "error"
      )
    }
    df <- df[as.character(df$resource) %in% as.character(resource), , drop = FALSE]
  }
  if (!is.null(condition)) {
    available_condition <- if ("condition" %in% colnames(df)) {
      unique(as.character(df$condition))
    } else {
      character(0)
    }
    available_condition <- available_condition[
      !is.na(available_condition) & nzchar(available_condition)
    ]
    if (length(available_condition) == 0L) {
      log_message(
        "The selected CCC result does not provide {.arg condition} provenance",
        message_type = "error"
      )
    }
    df <- df[
      as.character(df$condition) %in% as.character(condition), ,
      drop = FALSE
    ]
  }
  if (!is.null(sample)) {
    sample_col <- c("sample", "context", "dataset")
    sample_col <- sample_col[sample_col %in% colnames(df)][1]
    available_sample <- if (!is.na(sample_col)) {
      unique(as.character(df[[sample_col]]))
    } else {
      character(0)
    }
    available_sample <- available_sample[
      !is.na(available_sample) & nzchar(available_sample)
    ]
    if (is.na(sample_col) || length(available_sample) == 0L) {
      log_message(
        "The selected CCC result does not provide sample or context provenance",
        message_type = "error"
      )
    }
    df <- df[as.character(df[[sample_col]]) %in% as.character(sample), , drop = FALSE]
  }
  df
}

ccc_prepare_filtered_object <- function(
  srt,
  method,
  resource = NULL,
  condition = NULL,
  sample = NULL
) {
  if (is.null(resource) && is.null(condition) && is.null(sample)) {
    return(srt)
  }
  method <- normalize_ccc_method(method)
  targets <- unique(c("CCC", if (!identical(method, "CCC")) method))
  targets <- targets[targets %in% names(srt@tools)]
  for (target in targets) {
    bundle <- srt@tools[[target]]
    for (field in c("long_table", "primary_table", "consensus_table")) {
      if (is.data.frame(bundle[[field]])) {
        bundle[[field]] <- ccc_filter_table_context(
          bundle[[field]],
          resource = resource,
          condition = condition,
          sample = sample
        )
      }
    }
    srt@tools[[target]] <- bundle
  }
  srt
}

ccc_bind_long_tables <- function(pieces) {
  pieces <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, pieces)
  if (length(pieces) == 0L) {
    return(data.frame())
  }
  common <- Reduce(union, lapply(pieces, colnames))
  pieces <- lapply(pieces, function(x) {
    missing <- setdiff(common, colnames(x))
    for (nm in missing) x[[nm]] <- NA
    x[, common, drop = FALSE]
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

ccc_combine_methods <- function(
  df,
  mode = c("separate", "support", "rank", "legacy")
) {
  mode <- match.arg(mode)
  df <- ccc_semantic_long_table(df)
  if (nrow(df) == 0L || identical(mode, "separate")) {
    return(df)
  }
  if (identical(mode, "legacy")) {
    log_message(
      "{.val combine_methods = 'legacy'} sums backend scores on incompatible scales and is deprecated",
      message_type = "warning"
    )
    return(df)
  }
  required <- c("sender", "receiver", "ligand", "receptor", "method")
  if (!all(required %in% colnames(df))) {
    log_message("Unified CCC data are missing method or interaction identifiers", message_type = "error")
  }
  has_interaction <- !is.na(df$ligand) & nzchar(trimws(as.character(df$ligand))) &
    !is.na(df$receptor) & nzchar(trimws(as.character(df$receptor)))
  df <- df[has_interaction, , drop = FALSE]
  if (nrow(df) == 0L) {
    return(df)
  }
  methods <- unique(as.character(df$method))
  methods <- methods[!is.na(methods) & nzchar(methods)]
  df$.ccc_method_percentile <- 1
  for (method_i in methods) {
    idx <- which(as.character(df$method) == method_i)
    raw_rank <- suppressWarnings(as.numeric(df$priority_rank[idx]))
    finite <- is.finite(raw_rank)
    if (any(finite)) {
      percentile <- (rank(raw_rank[finite], ties.method = "average") - 1) /
        max(1, sum(finite) - 1)
      df$.ccc_method_percentile[idx[finite]] <- percentile
    }
  }
  key_cols <- c("sender", "receiver", "ligand", "receptor")
  key <- do.call(paste, c(lapply(df[key_cols], as.character), sep = "\r"))
  split_idx <- split(seq_len(nrow(df)), key, drop = TRUE)
  rows <- lapply(split_idx, function(idx) {
    x <- df[idx, , drop = FALSE]
    row <- x[1, , drop = FALSE]
    support_methods <- unique(as.character(x$method))
    support_methods <- support_methods[!is.na(support_methods) & nzchar(support_methods)]
    row$method <- "CCC"
    row$support_methods <- paste(sort(support_methods), collapse = ";")
    row$support_count <- length(support_methods)
    row$support_fraction <- length(support_methods) / max(1L, length(methods))
    if (identical(mode, "support")) {
      row$score <- row$support_count
      row$priority_score <- row$support_fraction
      row$priority_rank <- 1 - row$support_fraction
      row$score_type <- "method_support_count"
      row$support_type <- "cross_method_support"
    } else {
      rank_values <- suppressWarnings(as.numeric(x$.ccc_method_percentile))
      ranks <- tapply(
        rank_values,
        as.character(x$method),
        function(value) {
          value <- value[is.finite(value)]
          if (length(value)) mean(value) else 1
        }
      )
      ranks <- as.numeric(ranks)
      missing_methods <- max(0L, length(methods) - length(ranks))
      ranks <- c(ranks, rep(1, missing_methods))
      mean_rank <- if (length(ranks)) mean(ranks) else 1
      row$score <- 1 - mean_rank
      row$priority_score <- 1 - mean_rank
      row$priority_rank <- mean_rank
      row$score_type <- "mean_within_method_percentile"
      row$support_type <- "scop_visualization_consensus"
    }
    row$.ccc_method_percentile <- NULL
    row$pvalue <- NA_real_
    row$pvalue_type <- "not_available"
    row$significant <- row$support_count > 0L
    row
  })
  ccc_bind_long_tables(rows)
}

ccc_prepare_combined_object <- function(
  srt,
  method,
  combine_methods = c("separate", "support", "rank", "legacy")
) {
  combine_methods <- match.arg(combine_methods)
  method <- normalize_ccc_method(method)
  if (!identical(method, "CCC") || identical(combine_methods, "separate")) {
    return(srt)
  }
  bundle <- srt@tools[["CCC"]]
  if (is.null(bundle)) {
    bundle <- ccc_build_unified_bundle(srt)
  }
  liana_methods <- tolower(as.character(
    srt@tools[["LIANA"]]$parameters$method %||% character(0)
  ))
  if (
    combine_methods %in% c("support", "rank") &&
      all(c("CellphoneDB", "LIANA") %in% (bundle$methods %||% character(0))) &&
      "cellphonedb" %in% liana_methods
  ) {
    log_message(
      paste0(
        "Standalone CellphoneDB and the LIANA consensus are not independent ",
        "evidence because LIANA includes its CellPhoneDB scoring method; ",
        "interpret this as backend agreement, not independent validation"
      ),
      message_type = "warning"
    )
  }
  bundle$long_table <- ccc_combine_methods(bundle$long_table, mode = combine_methods)
  bundle$pair_table <- aggregate_ccc_long(bundle$long_table, backend = "r")
  bundle$metadata$combine_methods <- combine_methods
  srt@tools[["CCC"]] <- bundle
  srt
}

ccc_plot_methods_separately <- function(call, srt, env) {
  unified <- ccc_semantic_long_table(srt@tools[["CCC"]]$long_table)
  methods <- unique(as.character(unified$method))
  methods <- methods[!is.na(methods) & nzchar(methods) & methods != "CCC"]
  plots <- lapply(methods, function(method) {
    method_table <- unified[as.character(unified$method) == method, , drop = FALSE]
    bundle <- srt@tools[[method]] %||% list(method = method)
    bundle$long_table <- method_table
    bundle$primary_table <- method_table
    bundle$pair_table <- aggregate_ccc_long(method_table, backend = "r")
    method_object <- srt
    method_object@tools[[method]] <- bundle
    next_call <- call
    next_call$srt <- method_object
    next_call$method <- method
    next_call$combine_methods <- "legacy"
    next_call$resource <- NULL
    next_call$condition <- NULL
    next_call$sample <- NULL
    eval(next_call, envir = env)
  })
  names(plots) <- methods
  plots
}
