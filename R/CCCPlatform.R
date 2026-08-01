#' @title List SCOP cell-cell communication methods
#'
#' @description
#' Returns the CCC producers integrated with SCOP's Seurat result, scheduler,
#' unified table, and plotting interfaces.
#'
#' @return A data frame describing registered CCC methods and capabilities.
#' @export
ListCCCMethods <- function() {
  registry <- ccc_method_registry()
  rows <- lapply(registry, function(x) {
    data.frame(
      method = x$method,
      producer = x$producer,
      default = isTRUE(x$default),
      requirements = paste(x$required, collapse = ", "),
      result_types = paste(
        c(
          "long", "pair", "primary", "raw",
          if (identical(x$method, "LIANA")) "consensus",
          if (isTRUE(x$native_result)) "native"
        ),
        collapse = ", "
      ),
      generic_plots = paste(x$generic_plots, collapse = ", "),
      native_plots = paste(x$native_plots, collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

ccc_method_registry <- function() {
  generic <- c("heatmap", "dot", "tile", "circle", "chord", "arrow",
    "sigmoid", "bipartite", "embedding_network", "bar", "sankey",
    "box", "violin")
  list(
    CellChat = list(
      method = "CellChat", producer = "RunCellChat", default = TRUE,
      required = character(0), object_arg = "srt",
      generic_plots = generic,
      native_plots = c("pathway", "role", "ranknet", "diff_network"),
      native_result = TRUE
    ),
    CellphoneDB = list(
      method = "CellphoneDB", producer = "RunCellphoneDB", default = TRUE,
      required = character(0), object_arg = "srt",
      generic_plots = generic, native_plots = character(0)
    ),
    LIANA = list(
      method = "LIANA", producer = "RunLIANA", default = TRUE,
      required = character(0), object_arg = "srt",
      generic_plots = generic, native_plots = c("consensus_rank")
    ),
    Nichenetr = list(
      method = "Nichenetr", producer = "RunNichenetr", default = FALSE,
      required = "receiver", object_arg = "srt",
      generic_plots = generic,
      native_plots = c("ligand_target", "ligand_activity", "gene")
    ),
    MultiNichenetr = list(
      method = "MultiNichenetr", producer = "RunMultiNichenetr", default = FALSE,
      required = c(
        "sample.by", "condition.by", "condition_oi",
        "condition_reference", "receiver_celltypes"
      ),
      object_arg = "srt", generic_plots = generic,
      native_plots = c("ligand_target", "ligand_activity", "gene")
    ),
    SpatialCellChat = list(
      method = "SpatialCellChat", producer = "RunSpatialCellChat", default = FALSE,
      required = character(0), object_arg = "srt",
      generic_plots = generic,
      native_plots = c(
        "spatial_network", "lr_spatial", "pathway", "incoming", "outgoing",
        "diffusion"
      ),
      native_result = TRUE
    ),
    MDIC3 = list(
      method = "MDIC3", producer = "RunMDIC3", default = FALSE,
      required = "grn or grn_method", object_arg = "object",
      generic_plots = generic, native_plots = character(0)
    )
  )
}

ccc_registry_entry <- function(method) {
  method <- normalize_ccc_method(method)
  entry <- ccc_method_registry()[[method]]
  if (is.null(entry)) {
    log_message(
      "Unsupported CCC method {.val {method}}. Use {.fn ListCCCMethods} to inspect registered methods.",
      message_type = "error"
    )
  }
  entry
}

#' @title Inspect CCC results stored in a Seurat object
#'
#' @param srt A Seurat object.
#'
#' @return A data frame with one row per registered CCC method.
#' @export
CCCResultInfo <- function(srt) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  run_status <- srt@tools[["RunCCC"]]$status %||% data.frame()
  unified_methods <- ccc_unified_methods(srt)
  has_unified <- !is.null(srt@tools[["CCC"]])
  rows <- lapply(ccc_method_registry(), function(entry) {
    bundle <- srt@tools[[entry$method]]
    status <- if (is.null(bundle)) "absent" else if (
      is.data.frame(bundle$primary_table %||% bundle$long_table) &&
        nrow(bundle$primary_table %||% bundle$long_table) > 0L
    ) {
      "completed"
    } else {
      "incomplete"
    }
    if (nrow(run_status) > 0L && entry$method %in% run_status$method) {
      run_value <- run_status$status[match(entry$method, run_status$method)]
      if (!is.na(run_value) && nzchar(run_value)) status <- run_value
    }
    if (
      identical(status, "completed") && has_unified &&
        !entry$method %in% unified_methods
    ) {
      status <- "stale"
    }
    data.frame(
      method = entry$method,
      status = status,
      rows = if (is.null(bundle)) 0L else nrow(bundle$primary_table %||% bundle$long_table %||% data.frame()),
      schema = as.character(bundle$metadata$schema %||% NA_character_),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @title Get a standardized or native CCC result
#'
#' @param srt A Seurat object.
#' @param method A registered CCC method.
#' @param type Result representation.
#' @param resource Optional LIANA resource when retrieving a consensus.
#' @param condition Stored CellChat or SpatialCellChat result name when
#' retrieving a native object.
#' @param sample Stored SpatialCellChat sample when retrieving a native object.
#'
#' @return The requested stored result.
#' @export
GetCCCResult <- function(
  srt,
  method,
  type = c("primary", "long", "pair", "consensus", "raw", "native"),
  resource = NULL,
  condition = NULL,
  sample = NULL
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  method <- normalize_ccc_method(method)
  type <- match.arg(type)
  bundle <- get_bundle(srt, method)
  semantic_method <- if (identical(method, "CCC")) NULL else method
  if (
    identical(method, "CellChat") && type %in% c("primary", "long", "pair")
  ) {
    cellchat_long <- bundle$primary_table %||% bundle$long_table
    if (!is.data.frame(cellchat_long) || nrow(cellchat_long) == 0L) {
      cellchat_long <- ccc_build_cellchat_long_table(srt)
    }
    if (identical(type, "pair")) {
      return(bundle$pair_table %||% aggregate_ccc_long(cellchat_long))
    }
    return(ccc_semantic_long_table(cellchat_long, method = method))
  }
  out <- switch(type,
    primary = ccc_semantic_long_table(
      bundle$primary_table %||% bundle$long_table,
      method = semantic_method
    ),
    long = ccc_semantic_long_table(bundle$long_table, method = semantic_method),
    pair = bundle$primary_pair_table %||% bundle$pair_table,
    consensus = {
      if (!is.null(resource)) {
        bundle$consensus_by_resource[[resource]]
      } else {
        bundle$consensus_table
      }
    },
    raw = bundle$results %||% bundle$raw_result %||% bundle$raw,
    native = {
      if (identical(method, "CellChat")) {
        GetCCCObject(
          object = srt,
          method = "CellChat",
          result.name = condition
        )
      } else if (identical(method, "SpatialCellChat")) {
        GetCCCObject(
          object = srt,
          method = "SpatialCellChat",
          result.name = condition,
          sample = sample
        )
      } else {
        bundle$native_object %||% bundle$cellchat_object
      }
    }
  )
  if (is.null(out)) {
    log_message(
      "Result type {.val {type}} is not available for {.val {method}}",
      message_type = "error"
    )
  }
  out
}

ccc_semantic_long_table <- function(df, method = NULL) {
  out <- standardize_long_df(df)
  if (nrow(out) == 0L) return(out)
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
      "backend_support", "not_available")
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
      ties.method = "average")
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
  if (!is.data.frame(df) || nrow(df) == 0L) return(df)
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
      as.character(df$condition) %in% as.character(condition),
      ,
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
  if (is.null(resource) && is.null(condition) && is.null(sample)) return(srt)
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
  if (length(pieces) == 0L) return(data.frame())
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
  if (nrow(df) == 0L || identical(mode, "separate")) return(df)
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
    warning(
      paste0(
        "Standalone CellphoneDB and the LIANA consensus are not independent ",
        "evidence because LIANA includes its CellPhoneDB scoring method; ",
        "interpret this as backend agreement, not independent validation"
      ),
      call. = FALSE
    )
  }
  bundle$long_table <- ccc_combine_methods(bundle$long_table, mode = combine_methods)
  bundle$pair_table <- aggregate_ccc_long(bundle$long_table, backend = "r")
  bundle$metadata$combine_methods <- combine_methods
  srt@tools[["CCC"]] <- bundle
  srt
}

ccc_plot_methods_separately <- function(call, srt, env) {
  methods <- ccc_unified_methods(srt)
  methods <- methods[methods != "CCC"]
  plots <- lapply(methods, function(method) {
    next_call <- call
    next_call$method <- method
    next_call$combine_methods <- "legacy"
    eval(next_call, envir = env)
  })
  names(plots) <- methods
  plots
}
