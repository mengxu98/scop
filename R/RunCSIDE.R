#' @title Run C-SIDE spatial differential expression
#'
#' @description
#' Run `spacexr` C-SIDE after RCTD to test cell type-specific spatial or
#' condition-aware differential expression.
#'
#' @md
#' @inheritParams RunRCTD
#' @param layer Assay layer used as the spatial expression source.
#' @param rctd_result Optional old-api `spacexr` RCTD object. If `NULL`,
#' `srt@tools[["RCTD"]]$object` from [RunRCTD()] is used.
#' @param explanatory.variable Named numeric vector used by
#' `spacexr::run.CSIDE.single()`. Names must match spatial spot names.
#' @param group.by Metadata column used to build C-SIDE regions when
#' `region_list` is not supplied.
#' @param condition.by Binary metadata column converted to a 0/1 explanatory
#' variable for `spacexr::run.CSIDE.single()`.
#' @param design Numeric design matrix used by `spacexr::run.CSIDE()`.
#' Row names must match `barcodes` or spatial spot names.
#' @param region_list Named list of barcode vectors used by
#' `spacexr::run.CSIDE.regions()`.
#' @param barcodes Barcodes used by `spacexr::run.CSIDE()` or
#' `spacexr::run.CSIDE.intercept()`. If `NULL`, row names of `design` or all
#' spatial spots are used where applicable.
#' @param mode C-SIDE mode. `"auto"` dispatches to `"general"`, `"regions"`,
#' `"single"`, or `"intercept"` based on the supplied design inputs.
#' @param celltypes Optional cell types passed to C-SIDE as `cell_types`.
#' @param features Optional features retained in the normalized result table.
#' C-SIDE itself still applies its own gene filtering through backend
#' parameters such as `gene_threshold`.
#' @param tool_name Name used to store detailed C-SIDE results in `srt@tools`.
#' @param ... Additional named parameters passed to the selected C-SIDE
#' backend, such as `cell_type_threshold`, `gene_threshold`,
#' `doublet_mode`, `cell_type_specific`, or `params_to_test`. When using the
#' stored result from [RunRCTD()] with `rctd_mode = "full"`, `doublet_mode`
#' defaults to `FALSE` unless explicitly supplied.
#'
#' @return A `Seurat` object with C-SIDE summary metadata and detailed results
#' stored in `srt@tools[[tool_name]]` when `store_results = TRUE`.
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- subset(
#'   visium_human_pancreas_sub,
#'   cells = colnames(visium_human_pancreas_sub)[1:120],
#'   features = rownames(visium_human_pancreas_sub)[1:400]
#' )
#' spatial$region <- ifelse(spatial$x > stats::median(spatial$x), "right", "left")
#' spatial$CSIDE_n_sig <- ifelse(spatial$region == "right", 12, 4)
#' spatial$CSIDE_mode <- "regions"
#' cside_result <- data.frame(
#'   feature = rownames(spatial)[1:4],
#'   celltype = rep(c("Ductal", "Endocrine"), each = 2),
#'   parameter = "right_vs_left",
#'   logFC = c(1.2, 0.8, -0.9, -1.1),
#'   statistic = c(4.1, 3.5, -3.2, -3.8),
#'   p_value = c(0.001, 0.004, 0.006, 0.002),
#'   q_value = c(0.004, 0.008, 0.010, 0.006),
#'   significant = TRUE,
#'   method = "regions"
#' )
#' spatial@tools$CSIDE <- list(
#'   result_table = cside_result,
#'   parameters = list(mode = "regions", group.by = "region")
#' )
#'
#' SpatialSpotPlot(
#'   spatial,
#'   group.by = "CSIDE_n_sig",
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#' SpatialSpotPlot(
#'   spatial,
#'   features = cside_result$feature[1:2],
#'   assay = "Spatial",
#'   layer = "counts",
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#'
#' if (
#'   isTRUE(check_r("dmcable/spacexr", verbose = FALSE)) &&
#'     !is.null(spatial@tools$RCTD) &&
#'     inherits(spatial@tools$RCTD$object, "RCTD")
#' ) {
#'   spatial <- RunCSIDE(
#'     spatial,
#'     group.by = "region",
#'     celltypes = c("Ductal", "Endocrine"),
#'     gene_threshold = 0.00005,
#'     cell_type_threshold = 125
#'   )
#' }
RunCSIDE <- function(
  srt,
  rctd_result = NULL,
  explanatory.variable = NULL,
  group.by = NULL,
  condition.by = NULL,
  design = NULL,
  region_list = NULL,
  barcodes = NULL,
  mode = c("auto", "single", "regions", "general", "intercept"),
  assay = NULL,
  layer = "counts",
  celltypes = NULL,
  features = NULL,
  prefix = "CSIDE",
  tool_name = "CSIDE",
  store_results = TRUE,
  verbose = TRUE,
  ...
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  cside_assert_scalar_string(prefix, "prefix")
  cside_assert_scalar_string(tool_name, "tool_name")
  cside_assert_null_or_scalar_string(group.by, "group.by")
  cside_assert_null_or_scalar_string(condition.by, "condition.by")
  mode <- match.arg(mode)
  assay <- assay %||% SeuratObject::DefaultAssay(srt)

  stored_rctd_mode <- cside_stored_rctd_mode(
    srt = srt,
    rctd_result = rctd_result
  )
  rctd_result <- cside_resolve_rctd_result(
    srt = srt,
    rctd_result = rctd_result
  )
  extra_args <- list(...)
  cside_validate_named_param_list(extra_args, "...")
  extra_args <- cside_apply_rctd_mode_defaults(
    extra_args = extra_args,
    stored_rctd_mode = stored_rctd_mode
  )
  mode <- cside_resolve_mode(
    mode = mode,
    explanatory.variable = explanatory.variable,
    group.by = group.by,
    condition.by = condition.by,
    design = design,
    region_list = region_list,
    barcodes = barcodes
  )
  backend_fun <- cside_get_backend_fun(mode)

  inputs <- cside_prepare_backend_inputs(
    srt = srt,
    mode = mode,
    explanatory.variable = explanatory.variable,
    group.by = group.by,
    condition.by = condition.by,
    design = design,
    region_list = region_list,
    barcodes = barcodes
  )

  log_message(
    "Run {.pkg spacexr} C-SIDE in {.val {mode}} mode",
    verbose = verbose
  )
  result <- cside_run_backend(
    backend_fun = backend_fun,
    mode = mode,
    rctd_result = rctd_result,
    inputs = inputs,
    celltypes = celltypes,
    extra_args = extra_args
  )

  extracted <- cside_extract_result(
    result = result,
    mode = mode,
    features = features
  )
  summary_meta <- cside_summary_metadata(
    srt = srt,
    prefix = prefix,
    mode = mode,
    result_table = extracted$result_table
  )
  srt <- Seurat::AddMetaData(srt, metadata = summary_meta)

  if (isTRUE(store_results)) {
    srt@tools[[tool_name]] <- list(
      result = result,
      result_table = extracted$result_table,
      sig_gene_list = extracted$sig_gene_list,
      all_gene_list = extracted$all_gene_list,
      gene_fits = extracted$gene_fits,
      design = inputs$design,
      region_list = inputs$region_list,
      barcodes = inputs$barcodes,
      parameters = list(
        assay = assay,
        layer = layer,
        mode = mode,
        group.by = group.by,
        condition.by = condition.by,
        celltypes = celltypes,
        features = features,
        prefix = prefix,
        tool_name = tool_name,
        backend_args = extra_args,
        condition_levels = inputs$condition_levels
      )
    )
  }

  log_message(
    "{.pkg C-SIDE} results stored in {.code srt@tools[[{tool_name}]]}",
    verbose = verbose
  )
  srt
}

cside_stored_rctd_mode <- function(srt, rctd_result = NULL) {
  if (!is.null(rctd_result)) {
    return(NULL)
  }
  stored <- srt@tools[["RCTD"]]
  if (is.null(stored) || is.null(stored$parameters$rctd_mode)) {
    return(NULL)
  }
  stored$parameters$rctd_mode
}

cside_apply_rctd_mode_defaults <- function(extra_args, stored_rctd_mode = NULL) {
  if (
    is.null(extra_args$doublet_mode) &&
      length(stored_rctd_mode) == 1L &&
      identical(as.character(stored_rctd_mode), "full")
  ) {
    extra_args$doublet_mode <- FALSE
  }
  extra_args
}

cside_resolve_rctd_result <- function(srt, rctd_result = NULL) {
  if (!is.null(rctd_result)) {
    return(cside_validate_rctd_result(rctd_result))
  }
  stored <- srt@tools[["RCTD"]]
  if (is.null(stored) || !is.list(stored) || is.null(stored$object)) {
    log_message(
      paste(
        "No RCTD object was supplied.",
        "Run {.fn RunRCTD} with an old-api {.pkg spacexr} backend first,",
        "or provide {.arg rctd_result} explicitly."
      ),
      message_type = "error"
    )
  }
  cside_validate_rctd_result(stored$object)
}

cside_validate_rctd_result <- function(rctd_result) {
  if (!inherits(rctd_result, "RCTD")) {
    log_message(
      paste(
        "{.arg rctd_result} must be an old-api {.pkg spacexr} {.cls RCTD}",
        "object. C-SIDE is not available for the Bioconductor",
        "{.pkg spacexr} createRctd/runRctd result object."
      ),
      message_type = "error"
    )
  }
  rctd_result
}

cside_resolve_mode <- function(
  mode,
  explanatory.variable,
  group.by,
  condition.by,
  design,
  region_list,
  barcodes
) {
  if (!identical(mode, "auto")) {
    return(mode)
  }
  if (!is.null(design)) {
    return("general")
  }
  if (!is.null(region_list) || !is.null(group.by)) {
    return("regions")
  }
  if (!is.null(explanatory.variable) || !is.null(condition.by)) {
    return("single")
  }
  if (!is.null(barcodes)) {
    return("intercept")
  }
  "intercept"
}

cside_prepare_backend_inputs <- function(
  srt,
  mode,
  explanatory.variable,
  group.by,
  condition.by,
  design,
  region_list,
  barcodes
) {
  switch(
    mode,
    single = cside_prepare_single_input(
      srt = srt,
      explanatory.variable = explanatory.variable,
      condition.by = condition.by
    ),
    regions = cside_prepare_regions_input(
      srt = srt,
      group.by = group.by,
      region_list = region_list
    ),
    general = cside_prepare_general_input(
      srt = srt,
      design = design,
      barcodes = barcodes
    ),
    intercept = cside_prepare_intercept_input(
      srt = srt,
      barcodes = barcodes
    )
  )
}

cside_prepare_single_input <- function(
  srt,
  explanatory.variable = NULL,
  condition.by = NULL
) {
  if (!is.null(explanatory.variable) && !is.null(condition.by)) {
    log_message(
      "Provide only one of {.arg explanatory.variable} or {.arg condition.by}",
      message_type = "error"
    )
  }
  condition_levels <- NULL
  if (!is.null(condition.by)) {
    if (!condition.by %in% colnames(srt@meta.data)) {
      log_message(
        "{.arg condition.by} {.val {condition.by}} is not in srt meta.data",
        message_type = "error"
      )
    }
    condition <- srt[[condition.by, drop = TRUE]]
    keep <- !is.na(condition)
    levels_use <- unique(as.character(condition[keep]))
    if (length(levels_use) != 2L) {
      log_message(
        "{.arg condition.by} must contain exactly two non-missing groups",
        message_type = "error"
      )
    }
    condition_levels <- levels_use
    explanatory.variable <- ifelse(
      as.character(condition) == levels_use[[2L]],
      1,
      0
    )
    explanatory.variable[!keep] <- NA_real_
    names(explanatory.variable) <- colnames(srt)
  }
  if (is.null(explanatory.variable)) {
    log_message(
      "{.arg explanatory.variable} or {.arg condition.by} is required for {.arg mode = 'single'}",
      message_type = "error"
    )
  }
  explanatory.variable <- cside_align_named_vector(
    explanatory.variable,
    cells = colnames(srt),
    arg = "explanatory.variable"
  )
  keep <- is.finite(explanatory.variable)
  if (sum(keep) < 2L || length(unique(explanatory.variable[keep])) < 2L) {
    log_message(
      "{.arg explanatory.variable} must contain at least two finite values",
      message_type = "error"
    )
  }
  list(
    explanatory.variable = explanatory.variable,
    barcodes = names(explanatory.variable)[keep],
    condition_levels = condition_levels,
    design = NULL,
    region_list = NULL
  )
}

cside_prepare_regions_input <- function(
  srt,
  group.by = NULL,
  region_list = NULL
) {
  if (!is.null(region_list)) {
    region_list <- cside_normalize_region_list(
      region_list = region_list,
      cells = colnames(srt)
    )
  } else {
    if (is.null(group.by) || !group.by %in% colnames(srt@meta.data)) {
      log_message(
        "{.arg group.by} must be a metadata column when {.arg region_list} is NULL",
        message_type = "error"
      )
    }
    groups <- srt[[group.by, drop = TRUE]]
    keep <- !is.na(groups)
    region_list <- split(colnames(srt)[keep], as.character(groups[keep]))
    region_list <- cside_normalize_region_list(
      region_list = region_list,
      cells = colnames(srt)
    )
  }
  list(
    explanatory.variable = NULL,
    barcodes = unique(unlist(region_list, use.names = FALSE)),
    condition_levels = NULL,
    design = NULL,
    region_list = region_list
  )
}

cside_prepare_general_input <- function(
  srt,
  design,
  barcodes = NULL
) {
  if (is.null(design)) {
    log_message(
      "{.arg design} is required for {.arg mode = 'general'}",
      message_type = "error"
    )
  }
  design <- as.data.frame(design, check.names = FALSE)
  is_numeric <- vapply(design, is.numeric, logical(1))
  if (!all(is_numeric)) {
    log_message(
      "{.arg design} must contain numeric columns only",
      message_type = "error"
    )
  }
  if (
    is.null(rownames(design)) ||
      any(!nzchar(rownames(design))) ||
      identical(rownames(design), as.character(seq_len(nrow(design))))
  ) {
    if (nrow(design) != ncol(srt)) {
      log_message(
        "{.arg design} must have row names or one row per spatial spot",
        message_type = "error"
      )
    }
    rownames(design) <- colnames(srt)
  }
  if (is.null(barcodes)) {
    barcodes <- rownames(design)
  }
  barcodes <- cside_validate_barcodes(barcodes, cells = rownames(design))
  design <- as.matrix(design[barcodes, , drop = FALSE])
  list(
    explanatory.variable = NULL,
    barcodes = barcodes,
    condition_levels = NULL,
    design = design,
    region_list = NULL
  )
}

cside_prepare_intercept_input <- function(srt, barcodes = NULL) {
  barcodes <- barcodes %||% colnames(srt)
  barcodes <- cside_validate_barcodes(barcodes, cells = colnames(srt))
  list(
    explanatory.variable = NULL,
    barcodes = barcodes,
    condition_levels = NULL,
    design = NULL,
    region_list = NULL
  )
}

cside_run_backend <- function(
  backend_fun,
  mode,
  rctd_result,
  inputs,
  celltypes = NULL,
  extra_args = list()
) {
  base_args <- switch(
    mode,
    single = list(
      myRCTD = rctd_result,
      explanatory.variable = inputs$explanatory.variable,
      cell_types = celltypes
    ),
    regions = list(
      myRCTD = rctd_result,
      region_list = inputs$region_list,
      cell_types = celltypes
    ),
    general = list(
      myRCTD = rctd_result,
      X = inputs$design,
      barcodes = inputs$barcodes,
      cell_types = celltypes
    ),
    intercept = list(
      myRCTD = rctd_result,
      barcodes = inputs$barcodes,
      cell_types = celltypes
    )
  )
  duplicate_args <- intersect(names(extra_args), names(base_args))
  if (length(duplicate_args) > 0L) {
    log_message(
      "{.arg ...} duplicates explicit C-SIDE argument(s): {.val {duplicate_args}}",
      message_type = "error"
    )
  }
  do.call(backend_fun, c(base_args, extra_args))
}

cside_get_backend_fun <- function(mode) {
  check_r("dmcable/spacexr", verbose = FALSE)
  fun <- switch(
    mode,
    single = "run.CSIDE.single",
    regions = "run.CSIDE.regions",
    general = "run.CSIDE",
    intercept = "run.CSIDE.intercept"
  )
  tryCatch(
    get_namespace_fun("spacexr", fun),
    error = function(e) {
      log_message(
        paste(
          "{.pkg spacexr} does not expose {.fn {fun}}.",
          "Install the original C-SIDE backend with",
          "{.code thisutils::check_r('dmcable/spacexr')}."
        ),
        message_type = "error"
      )
    }
  )
}

cside_extract_result <- function(result, mode, features = NULL) {
  if (is.list(result) && !is.null(result$result_table)) {
    result_table <- cside_normalize_result_table(result$result_table, mode = mode)
    if (!is.null(features)) {
      result_table <- result_table[result_table$feature %in% features, , drop = FALSE]
    }
    return(list(
      result_table = result_table,
      sig_gene_list = result$sig_gene_list %||% NULL,
      all_gene_list = result$all_gene_list %||% NULL,
      gene_fits = result$gene_fits %||% NULL
    ))
  }

  de_results <- cside_result_item(result, "de_results")
  internal_vars_de <- cside_result_item(result, "internal_vars_de")
  sig_gene_list <- de_results$sig_gene_list %||% NULL
  all_gene_list <- de_results$all_gene_list %||% NULL
  gene_fits <- de_results$gene_fits %||% NULL

  result_table <- cside_build_result_table(
    all_gene_list = all_gene_list,
    sig_gene_list = sig_gene_list,
    mode = mode
  )
  if (!is.null(features)) {
    result_table <- result_table[result_table$feature %in% features, , drop = FALSE]
  }
  if (nrow(result_table) == 0L && !is.null(internal_vars_de$cell_types)) {
    result_table <- cside_empty_result_table()
  }

  list(
    result_table = result_table,
    sig_gene_list = sig_gene_list,
    all_gene_list = all_gene_list,
    gene_fits = gene_fits
  )
}

cside_result_item <- function(result, item) {
  if (isS4(result) && methods::is(result, "RCTD")) {
    return(methods::slot(result, item))
  }
  if (is.list(result)) {
    return(result[[item]] %||% list())
  }
  list()
}

cside_build_result_table <- function(all_gene_list, sig_gene_list, mode) {
  if (is.null(all_gene_list)) {
    return(cside_empty_result_table())
  }
  if (is.data.frame(all_gene_list)) {
    all_gene_list <- list(CSIDE = all_gene_list)
  }
  if (!is.list(all_gene_list) || length(all_gene_list) == 0L) {
    return(cside_empty_result_table())
  }

  rows <- lapply(names(all_gene_list), function(celltype) {
    df <- all_gene_list[[celltype]]
    if (!is.data.frame(df) || nrow(df) == 0L) {
      return(NULL)
    }
    df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
    feature <- cside_result_features(df)
    sig_features <- cside_sig_features(sig_gene_list, celltype)
    p_value <- cside_pick_numeric(df, c("p_value", "p_val", "p", "p_val_best"))
    q_value <- cside_pick_numeric(df, c("q_value", "q_val", "padj", "fdr"))
    if (all(is.na(q_value)) && any(is.finite(p_value))) {
      q_value <- stats::p.adjust(p_value, method = "BH")
    }
    data.frame(
      feature = feature,
      celltype = celltype,
      parameter = cside_result_parameter(df),
      logFC = cside_pick_numeric(df, c("logFC", "log_fc", "log_fc_best", "log_fc_est")),
      statistic = cside_pick_numeric(df, c("statistic", "Z_score", "z_score", "Z_est")),
      p_value = p_value,
      q_value = q_value,
      significant = feature %in% sig_features,
      method = mode,
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) {
    return(cside_empty_result_table())
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

cside_normalize_result_table <- function(result_table, mode) {
  df <- as.data.frame(result_table, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"feature" %in% colnames(df)) {
    df$feature <- cside_result_features(df)
  }
  if (!"celltype" %in% colnames(df)) {
    df$celltype <- "CSIDE"
  }
  if (!"parameter" %in% colnames(df)) {
    df$parameter <- NA_character_
  }
  if (!"logFC" %in% colnames(df)) {
    df$logFC <- cside_pick_numeric(df, c("log_fc", "log_fc_best", "log_fc_est"))
  }
  if (!"statistic" %in% colnames(df)) {
    df$statistic <- cside_pick_numeric(df, c("Z_score", "z_score", "Z_est"))
  }
  if (!"p_value" %in% colnames(df)) {
    df$p_value <- cside_pick_numeric(df, c("p_val", "p", "p_val_best"))
  }
  if (!"q_value" %in% colnames(df)) {
    df$q_value <- cside_pick_numeric(df, c("q_val", "padj", "fdr"))
  }
  if (!"significant" %in% colnames(df)) {
    df$significant <- is.finite(df$q_value) & df$q_value < 0.05
  }
  if (!"method" %in% colnames(df)) {
    df$method <- mode
  }
  df[, c(
    "feature", "celltype", "parameter", "logFC", "statistic",
    "p_value", "q_value", "significant", "method"
  ), drop = FALSE]
}

cside_empty_result_table <- function() {
  data.frame(
    feature = character(),
    celltype = character(),
    parameter = character(),
    logFC = numeric(),
    statistic = numeric(),
    p_value = numeric(),
    q_value = numeric(),
    significant = logical(),
    method = character(),
    stringsAsFactors = FALSE
  )
}

cside_result_features <- function(df) {
  hit <- intersect(c("feature", "gene", "genes", "gene_id"), colnames(df))
  if (length(hit) > 0L) {
    return(as.character(df[[hit[[1L]]]]))
  }
  rn <- rownames(df)
  if (!is.null(rn) && length(rn) == nrow(df) && !all(grepl("^[0-9]+$", rn))) {
    return(as.character(rn))
  }
  paste0("Feature", seq_len(nrow(df)))
}

cside_result_parameter <- function(df) {
  if ("paramindex_best" %in% colnames(df)) {
    return(as.character(df$paramindex_best))
  }
  if (all(c("paramindex1_best", "paramindex2_best") %in% colnames(df))) {
    return(paste(df$paramindex1_best, df$paramindex2_best, sep = "-"))
  }
  rep(NA_character_, nrow(df))
}

cside_pick_numeric <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0L) {
    return(rep(NA_real_, nrow(df)))
  }
  suppressWarnings(as.numeric(df[[hit[[1L]]]]))
}

cside_sig_features <- function(sig_gene_list, celltype) {
  if (is.null(sig_gene_list) || is.null(sig_gene_list[[celltype]])) {
    return(character())
  }
  sig <- sig_gene_list[[celltype]]
  if (is.data.frame(sig)) {
    return(cside_result_features(sig))
  }
  as.character(sig)
}

cside_summary_metadata <- function(srt, prefix, mode, result_table) {
  n_sig <- if (nrow(result_table) == 0L) {
    0L
  } else {
    sum(result_table$significant, na.rm = TRUE)
  }
  data.frame(
    rep(n_sig, ncol(srt)),
    rep(mode, ncol(srt)),
    row.names = colnames(srt),
    stringsAsFactors = FALSE
  ) |>
    stats::setNames(c(paste0(prefix, "_n_sig"), paste0(prefix, "_mode")))
}

cside_align_named_vector <- function(x, cells, arg) {
  nms <- names(x)
  x <- suppressWarnings(as.numeric(x))
  names(x) <- nms
  if (is.null(names(x))) {
    if (length(x) != length(cells)) {
      log_message(
        "{.arg {arg}} must be named or have one value per spatial spot",
        message_type = "error"
      )
    }
    names(x) <- cells
  }
  missing_cells <- setdiff(names(x), cells)
  if (length(missing_cells) > 0L) {
    log_message(
      "{.arg {arg}} contains names not present in {.arg srt}: {.val {missing_cells}}",
      message_type = "error"
    )
  }
  x[cells]
}

cside_normalize_region_list <- function(region_list, cells) {
  if (!is.list(region_list) || length(region_list) < 3L) {
      log_message(
        "{.arg region_list} must contain at least three regions for spacexr run.CSIDE.regions. Use {.arg condition.by} for binary comparisons.",
        message_type = "error"
      )
  }
  if (is.null(names(region_list)) || any(!nzchar(names(region_list)))) {
    names(region_list) <- paste0("region", seq_along(region_list))
  }
  region_list <- lapply(region_list, function(x) {
    x <- unique(as.character(x))
    x[x %in% cells]
  })
  keep <- lengths(region_list) > 0L
  region_list <- region_list[keep]
  if (length(region_list) < 3L) {
      log_message(
        "{.arg region_list} must contain at least three non-empty regions matching spatial spots for spacexr run.CSIDE.regions.",
        message_type = "error"
      )
  }
  region_list
}

cside_validate_barcodes <- function(barcodes, cells) {
  barcodes <- unique(as.character(barcodes))
  if (length(barcodes) == 0L || any(is.na(barcodes) | !nzchar(barcodes))) {
    log_message(
      "{.arg barcodes} must contain non-empty spot names",
      message_type = "error"
    )
  }
  missing <- setdiff(barcodes, cells)
  if (length(missing) > 0L) {
    log_message(
      "{.arg barcodes} are not present in the expected input: {.val {missing}}",
      message_type = "error"
    )
  }
  barcodes
}

cside_validate_named_param_list <- function(x, arg_name) {
  if (length(x) == 0L) {
    return(invisible(TRUE))
  }
  nms <- names(x)
  if (is.null(nms) || any(is.na(nms) | !nzchar(nms))) {
    log_message(
      "{.arg {arg_name}} must contain named arguments only",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

cside_assert_scalar_string <- function(x, arg) {
  if (
    is.null(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !nzchar(x)
  ) {
    log_message(
      "{.arg {arg}} must be a single non-empty string",
      message_type = "error"
    )
  }
}

cside_assert_null_or_scalar_string <- function(x, arg) {
  if (is.null(x)) {
    return(invisible(TRUE))
  }
  cside_assert_scalar_string(x, arg)
}
