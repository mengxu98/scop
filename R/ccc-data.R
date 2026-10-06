get_cc_obj <- function(x) {
  if (!inherits(x, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  x@tools[["CellChat"]]
}

cc_names <- function(srt) {
  store <- get_cc_obj(srt)
  if (is.null(store) || is.null(store$results)) {
    return(character(0))
  }
  names(store$results)
}

comparison_cc_names <- function(srt) {
  store <- get_cc_obj(srt)
  if (is.null(store) || is.null(store$comparisons)) {
    return(character(0))
  }
  names(store$comparisons)
}

resolve_single_cc_condition <- function(srt, condition = NULL) {
  result_names <- cc_names(srt)
  if (!is.null(condition)) {
    if (!condition %in% result_names) {
      log_message(
        "Condition {.val {condition}} not found in CellChat results",
        message_type = "error"
      )
    }
    return(condition)
  }
  if (length(result_names) == 1L) {
    return(result_names[1])
  }
  if ("ALL" %in% result_names) {
    return("ALL")
  }
  log_message(
    "Multiple CellChat results found. Please specify {.arg condition}",
    message_type = "error"
  )
}

use_cc_single_condition <- function(srt, condition = NULL) {
  result_names <- cc_names(srt)
  cmp_names <- comparison_cc_names(srt)
  if (!is.null(condition)) {
    if (condition %in% result_names) {
      return(TRUE)
    }
    if (condition %in% cmp_names) {
      return(FALSE)
    }
    log_message(
      "{.arg condition} must be one of CellChat result names or comparison names",
      message_type = "error"
    )
  }
  if (length(result_names) == 1L && length(cmp_names) == 0L) {
    return(TRUE)
  }
  if ("ALL" %in% result_names && length(cmp_names) == 0L) {
    return(TRUE)
  }
  if (length(cmp_names) == 1L && length(result_names) == 0L) {
    return(FALSE)
  }
  log_message(
    "The requested plot is ambiguous because both single-condition results and/or comparison results are available. Please specify {.arg condition}",
    message_type = "error"
  )
}

resolve_group_index_single <- function(object, group.use = NULL) {
  if (is.null(group.use) || length(group.use) == 0L) {
    return(NULL)
  }

  groups <- character(0)
  if (inherits(object, "CellChat")) {
    groups <- tryCatch(
      {
        if (!is.null(object@idents)) {
          unique(as.character(object@idents))
        } else {
          character(0)
        }
      },
      error = function(e) character(0)
    )
  }
  if (length(groups) == 0L) {
    groups <- unique(as.character(group.use))
  }

  group.use_chr <- as.character(group.use)
  if (length(group.use_chr) == 1L && identical(group.use_chr, "all")) {
    return(groups)
  }

  group.use_num <- suppressWarnings(as.integer(group.use_chr))
  if (
    length(group.use_chr) > 0L &&
      all(!is.na(group.use_num)) &&
      all(group.use_chr == as.character(group.use_num))
  ) {
    if (any(group.use_num < 1L) || any(group.use_num > length(groups))) {
      log_message(
        "{.arg group.use} indices are out of range for available groups",
        message_type = "error"
      )
    }
    return(groups[group.use_num])
  }

  missing <- setdiff(group.use_chr, groups)
  if (length(missing) > 0L) {
    log_message(
      "Unknown group labels in {.arg group.use}: {.val {missing}}",
      message_type = "error"
    )
  }

  groups[groups %in% group.use_chr]
}

get_single_cc_obj <- function(srt, condition = NULL) {
  store <- get_cc_obj(srt)
  condition <- resolve_single_cc_condition(srt, condition = condition)
  store$results[[condition]]$cellchat_object
}

cc_get_cmp <- function(srt, condition = NULL) {
  store <- get_cc_obj(srt)
  cmp_names <- comparison_cc_names(srt)
  if (is.null(condition)) {
    if (length(cmp_names) == 1L) {
      condition <- cmp_names[1]
    } else if (length(cmp_names) == 0L) {
      log_message(
        "No CellChat comparisons found",
        message_type = "error"
      )
    } else {
      log_message(
        "Multiple CellChat comparisons found. Please specify {.arg condition}",
        message_type = "error"
      )
    }
  }
  if (!condition %in% cmp_names) {
    log_message(
      "Comparison {.val {condition}} not found in CellChat results",
      message_type = "error"
    )
  }
  store$comparisons[[condition]]
}

get_dataset_object <- function(srt, condition = NULL, dataset = 1) {
  store <- get_cc_obj(srt)
  result_names <- cc_names(srt)
  cmp_names <- comparison_cc_names(srt)

  if (!is.null(condition) && condition %in% result_names) {
    return(list(
      object = store$results[[condition]]$cellchat_object,
      seurat_object = store$results[[condition]]$seurat_object,
      label = condition,
      source = "single"
    ))
  }

  if (!is.null(condition) && condition %in% cmp_names) {
    cmp <- cc_get_cmp(srt, condition = condition)
    ds_name <- cc_pick_dataset_name(cmp, dataset)
    seu_object <- NULL
    if (
      !is.null(store$results[[ds_name]]) &&
        !is.null(store$results[[ds_name]]$seurat_object)
    ) {
      seu_object <- store$results[[ds_name]]$seurat_object
    }
    return(list(
      object = cmp$object.list[[ds_name]],
      seurat_object = seu_object,
      label = ds_name,
      source = "comparison"
    ))
  }

  if (is.null(condition) && length(result_names) == 1L) {
    return(list(
      object = store$results[[result_names[1]]]$cellchat_object,
      seurat_object = store$results[[result_names[1]]]$seurat_object,
      label = result_names[1],
      source = "single"
    ))
  }

  if (is.null(condition) && "ALL" %in% result_names) {
    return(list(
      object = store$results[["ALL"]]$cellchat_object,
      seurat_object = store$results[["ALL"]]$seurat_object,
      label = "ALL",
      source = "single"
    ))
  }

  if (is.null(condition) && length(cmp_names) == 1L) {
    cmp <- cc_get_cmp(srt, condition = cmp_names[1])
    ds_name <- cc_pick_dataset_name(cmp, dataset)
    seu_object <- NULL
    if (
      !is.null(store$results[[ds_name]]) &&
        !is.null(store$results[[ds_name]]$seurat_object)
    ) {
      seu_object <- store$results[[ds_name]]$seurat_object
    }
    return(list(
      object = cmp$object.list[[ds_name]],
      seurat_object = seu_object,
      label = ds_name,
      source = "comparison"
    ))
  }

  log_message(
    "Unable to determine which CellChat object to plot. Please specify {.arg condition}",
    message_type = "error"
  )
}

cc_resolve_dataset_index <- function(cmp, comparison = c(1, 2)) {
  ds_names <- names(cmp$object.list)
  if (is.numeric(comparison)) {
    idx <- as.integer(comparison)
    if (any(is.na(idx)) || any(idx < 1L) || any(idx > length(ds_names))) {
      log_message(
        "comparison indices out of range. Available datasets: {.val {seq_along(ds_names)}}",
        message_type = "error"
      )
    }
    return(idx)
  }
  idx <- match(as.character(comparison), ds_names)
  if (any(is.na(idx))) {
    log_message(
      "comparison names not found: {.val {as.character(comparison)[is.na(idx)]}}. Available datasets: {.val {ds_names}}",
      message_type = "error"
    )
  }
  idx
}

cc_pick_dataset_name <- function(cmp, dataset = 1) {
  nm <- names(cmp$object.list)
  if (is.character(dataset)) {
    if (!dataset %in% nm) {
      log_message(
        "dataset {.val {dataset}} not found. Available: {.val {nm}}",
        message_type = "error"
      )
    }
    return(dataset)
  }
  if (length(dataset) != 1L || dataset < 1L || dataset > length(nm)) {
    log_message(
      "dataset index out of range. Available: {.val {seq_along(nm)}}",
      message_type = "error"
    )
  }
  nm[dataset]
}

cmp_cc_label <- function(cmp, comp_idx = c(1, 2)) {
  nm <- names(cmp$object.list)[comp_idx]
  if (length(nm) <= 1L) {
    return(nm[1])
  }
  paste0(nm[1], "_vs_", nm[2])
}

subset_cc_table <- function(
  object,
  slot.name = "net",
  signaling = NULL,
  pairLR.use = NULL,
  sources.use = NULL,
  targets.use = NULL,
  thresh = 0.05,
  dataset = NULL
) {
  check_r("jinworks/CellChat", verbose = FALSE)
  if (!is.null(pairLR.use) && !is.data.frame(pairLR.use)) {
    pairLR.use <- data.frame(
      interaction_name = as.character(pairLR.use),
      stringsAsFactors = FALSE
    )
  }
  if (!is.null(pairLR.use)) {
    signaling <- NULL
  }
  df <- get_namespace_fun("CellChat", "subsetCommunication")(
    object = object,
    slot.name = slot.name,
    sources.use = sources.use,
    targets.use = targets.use,
    signaling = signaling,
    pairLR.use = pairLR.use,
    thresh = thresh
  )
  df <- as.data.frame(df)
  if (!is.null(dataset) && nrow(df) > 0L) {
    df$dataset <- dataset
  }
  df
}

detect_method <- function(srt, method = NULL) {
  if (!is.null(method)) {
    return(normalize_ccc_method(method))
  }
  available <- intersect(names(srt@tools), ccc_registered_methods())
  if (length(available) == 1L) {
    return(available[1])
  }
  if ("CCC" %in% names(srt@tools)) {
    return("CCC")
  }
  if (length(available) == 0L) {
    log_message(
      "No cell-cell communication results were found in {.cls Seurat}",
      message_type = "error"
    )
  }
  log_message(
    "Multiple CCC methods are available. Please specify {.arg method}. Candidates: {.val {available}}",
    message_type = "error"
  )
}

get_bundle <- function(srt, method) {
  x <- srt@tools[[method]]
  if (is.null(x)) {
    log_message(
      "{.pkg {method}} results not found in {.cls Seurat}",
      message_type = "error"
    )
  }
  x
}

ccc_unified_methods <- function(srt) {
  bundle <- srt@tools[["CCC"]]
  if (is.null(bundle)) {
    return(character(0))
  }
  methods <- bundle$methods %||% NULL
  if (is.null(methods) && is.data.frame(bundle$long_table) && "method" %in% colnames(bundle$long_table)) {
    methods <- unique(as.character(bundle$long_table$method))
  }
  methods <- unique(as.character(methods %||% character(0)))
  methods[!is.na(methods) & nzchar(methods)]
}

ccc_has_unified_method <- function(srt, method = NULL) {
  if (!inherits(srt, "Seurat") || is.null(srt@tools[["CCC"]])) {
    return(FALSE)
  }
  if (is.null(method)) {
    return(TRUE)
  }
  method <- normalize_ccc_method(method)
  identical(method, "CCC") || method %in% ccc_unified_methods(srt)
}

ccc_require_coordinate_contracts <- function(srt, methods) {
  producers <- c(
    SpatialCellChat = "RunSpatialCellChat()",
    COMMOT = "RunCOMMOT()",
    SpaTalk = "RunSpaTalk()"
  )
  methods <- unique(vapply(
    methods %||% character(0),
    normalize_ccc_method,
    character(1)
  ))
  methods <- intersect(methods, names(producers))
  for (method in methods) {
    spatial_require_coordinate_contract(
      srt@tools[[method]],
      producers[[method]]
    )
  }
  invisible(methods)
}

ccc_long_table_for_method <- function(
  srt,
  method,
  condition = NULL,
  dataset = 1,
  slot.name = "net",
  signaling = NULL,
  pairLR.use = NULL,
  sources.use = NULL,
  targets.use = NULL,
  thresh = 0.05
) {
  method <- detect_method(srt = srt, method = method)
  ccc_require_coordinate_contracts(srt, method)
  use_cellchat_direct <- identical(method, "CellChat") &&
    (
      !is.null(condition) ||
        !isTRUE(dataset == 1) ||
        !identical(slot.name, "net") ||
        !is.null(signaling) ||
        !is.null(pairLR.use) ||
        !is.null(sources.use) ||
        !is.null(targets.use)
    )
  if (
    !isTRUE(use_cellchat_direct) &&
      (identical(method, "CCC") || ccc_has_unified_method(srt, method = method))
  ) {
    filter_method <- if (identical(method, "CCC")) NULL else method
    return(ccc_get_unified_long_table(
      srt = srt,
      method = filter_method,
      thresh = thresh
    ))
  }
  if (identical(method, "CellChat")) {
    return(extract_long_table(
      srt = srt,
      condition = condition,
      dataset = dataset,
      slot.name = slot.name,
      signaling = signaling,
      pairLR.use = pairLR.use,
      sources.use = sources.use,
      targets.use = targets.use,
      thresh = thresh
    ))
  }
  bundle <- get_bundle(srt, method = method)
  bundle$long_table %||% data.frame()
}

ccc_available_methods <- function(srt) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  intersect(names(srt@tools), ccc_registered_methods())
}

ccc_build_cellchat_long_table <- function(srt, thresh = 0.05) {
  store <- get_cc_obj(srt)
  if (is.null(store) || is.null(store$results) || length(store$results) == 0L) {
    return(data.frame())
  }
  pieces <- lapply(names(store$results), function(condition) {
    obj <- store$results[[condition]]$cellchat_object
    if (is.null(obj)) {
      return(data.frame())
    }
    out <- tryCatch(
      subset_cc_table(
        object = obj,
        slot.name = "net",
        thresh = thresh,
        dataset = condition
      ),
      error = function(e) data.frame()
    )
    if (is.null(out) || nrow(out) == 0L) {
      return(data.frame())
    }
    out
  })
  pieces <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, pieces)
  if (length(pieces) == 0L) {
    return(data.frame())
  }
  df <- do.call(rbind, pieces)
  rownames(df) <- NULL
  df <- standardize_long_df(df)
  df <- ccc_mark_significance(df, thresh = thresh)
  df$method <- "CellChat"
  df
}

ccc_bundle_long_table <- function(srt, method, bundle = NULL, thresh = 0.05) {
  method <- normalize_ccc_method(method)
  spec <- ccc_method_spec(method)
  if (isFALSE(spec$supports_unified_edges)) {
    return(data.frame())
  }
  if (identical(method, "CellChat")) {
    bundle <- bundle %||% srt@tools[[method]]
    stored <- bundle$primary_table %||% bundle$long_table
    rebuilt <- ccc_build_cellchat_long_table(srt, thresh = thresh)
    source <- if (is.data.frame(rebuilt) && nrow(rebuilt) > 0L) rebuilt else stored
    df <- ccc_semantic_long_table(
      source,
      method = method
    )
    df <- ccc_mark_significance(df, thresh = thresh)
    if ("pvalue" %in% colnames(df)) {
      pvalue <- suppressWarnings(as.numeric(df$pvalue))
      df <- df[!is.finite(pvalue) | pvalue <= thresh, , drop = FALSE]
    }
    if (nrow(df) > 0L) df$producer <- "RunCellChat"
    return(df)
  }
  bundle <- bundle %||% get_bundle(srt, method = method)
  df <- ccc_semantic_long_table(
    bundle$primary_table %||% bundle$long_table %||% data.frame(),
    method = method
  )
  if (nrow(df) == 0L) {
    return(df)
  }
  if (!"method" %in% colnames(df)) {
    df$method <- method
  }
  df$method <- method
  provenance <- bundle$provenance %||% list()
  df$producer <- provenance$producer %||% paste0("Run", method)
  backend_version <- provenance$backend_version %||%
    provenance$backend_versions %||% NA_character_
  if (length(backend_version) > 1L) {
    backend_version <- paste(
      paste(names(backend_version), backend_version, sep = "="),
      collapse = ";"
    )
  }
  df$backend_version <- as.character(backend_version)[1]
  ccc_mark_significance(df, thresh = thresh)
}

ccc_build_unified_bundle <- function(
  srt,
  methods = NULL,
  thresh = 0.05,
  backend = c("cpp", "r")
) {
  backend <- match.arg(backend)
  methods <- methods %||% ccc_available_methods(srt)
  methods <- unique(vapply(methods, normalize_ccc_method, character(1)))
  methods <- setdiff(methods, "CCC")
  ccc_require_coordinate_contracts(srt, methods)
  edge_methods <- methods[vapply(methods, function(method) {
    isTRUE(ccc_method_spec(method)$supports_unified_edges)
  }, logical(1))]
  pieces <- lapply(edge_methods, function(method) {
    ccc_bundle_long_table(
      srt = srt,
      method = method,
      bundle = srt@tools[[method]],
      thresh = thresh
    )
  })
  long_table <- ccc_bind_long_tables(pieces)
  long_table <- ccc_semantic_long_table(long_table)
  pair_table <- aggregate_ccc_long(long_table, backend = backend)
  liana_sample_col <- if ("method" %in% colnames(long_table)) "method" else NULL
  liana_table <- ccc_long_to_liana(long_table, sample_col = liana_sample_col)
  list(
    method = "CCC",
    methods = sort(unique(as.character(long_table$method %||% edge_methods))),
    association_methods = setdiff(methods, edge_methods),
    long_table = long_table,
    pair_table = pair_table,
    liana_table = liana_table,
    metadata = list(
      methods = methods,
      updated_at = as.character(Sys.time()),
      backend = backend,
      backend_scope = "result aggregation and unified-table construction"
    )
  )
}

ccc_update_unified_bundle <- function(
  srt,
  method,
  bundle = NULL,
  thresh = 0.05,
  backend = c("cpp", "r")
) {
  backend <- match.arg(backend)
  method <- normalize_ccc_method(method)
  if (isFALSE(ccc_method_spec(method)$supports_unified_edges)) {
    return(srt)
  }
  new_long <- ccc_bundle_long_table(
    srt = srt,
    method = method,
    bundle = bundle,
    thresh = thresh
  )
  old <- srt@tools[["CCC"]]
  old_long <- if (!is.null(old)) {
    standardize_long_df(old$long_table %||% data.frame())
  } else {
    data.frame()
  }
  if (nrow(old_long) > 0L && "method" %in% colnames(old_long)) {
    old_long <- old_long[as.character(old_long$method) != method, , drop = FALSE]
  }
  long_table <- ccc_bind_long_tables(list(old_long, new_long))
  long_table <- ccc_semantic_long_table(long_table)
  old_methods <- if (!is.null(old$methods)) {
    setdiff(as.character(old$methods), method)
  } else {
    character(0)
  }
  methods <- sort(unique(c(
    old_methods,
    method,
    as.character(long_table$method %||% character(0))
  )))
  srt@tools[["CCC"]] <- list(
    method = "CCC",
    methods = methods,
    long_table = long_table,
    pair_table = aggregate_ccc_long(long_table, backend = backend),
    liana_table = ccc_long_to_liana(
      long_table,
      sample_col = if ("method" %in% colnames(long_table)) "method" else NULL
    ),
    metadata = list(
      methods = methods,
      updated_method = method,
      updated_at = as.character(Sys.time()),
      backend = backend,
      backend_scope = "result aggregation and unified-table construction"
    )
  )
  srt
}

ccc_get_unified_long_table <- function(srt, method = NULL, thresh = 0.05) {
  method <- if (is.null(method)) NULL else normalize_ccc_method(method)
  bundle <- srt@tools[["CCC"]]
  contract_methods <- if (is.null(method) || identical(method, "CCC")) {
    unique(c(
      ccc_unified_methods(srt),
      if (is.null(bundle)) ccc_available_methods(srt) else character(0)
    ))
  } else {
    method
  }
  ccc_require_coordinate_contracts(srt, contract_methods)
  if (is.null(bundle)) {
    bundle <- ccc_build_unified_bundle(srt = srt, thresh = thresh)
  }
  df <- ccc_semantic_long_table(bundle$long_table %||% data.frame())
  if (nrow(df) == 0L) {
    return(df)
  }
  if (!is.null(method) && !"CCC" %in% method && "method" %in% colnames(df)) {
    df <- df[as.character(df$method) %in% method, , drop = FALSE]
  }
  ccc_mark_significance(df, thresh = thresh)
}

extract_long_table <- function(
  srt,
  condition = NULL,
  dataset = 1,
  slot.name = "net",
  signaling = NULL,
  pairLR.use = NULL,
  sources.use = NULL,
  targets.use = NULL,
  thresh = 0.05
) {
  info <- get_dataset_object(srt, condition = condition, dataset = dataset)
  df <- subset_cc_table(
    object = info$object,
    slot.name = slot.name,
    signaling = signaling,
    pairLR.use = pairLR.use,
    sources.use = sources.use,
    targets.use = targets.use,
    thresh = thresh,
    dataset = info$label
  )
  df <- standardize_long_df(df)
  df <- ccc_mark_significance(df, thresh = thresh)
  df$method <- "CellChat"
  df
}

standardize_long_df <- function(df) {
  df <- standardize_df(df)
  if (is.null(df) || nrow(df) == 0L) {
    return(data.frame())
  }

  rename_map <- list(
    sender = c("sender", "source"),
    receiver = c("receiver", "target"),
    interaction_name = c("interaction_name", "interacting_pair", "interaction"),
    ligand = c("ligand", "ligand_complex", "gene_a", "partner_a", "from"),
    receptor = c("receptor", "receptor_complex", "gene_b", "partner_b", "to"),
    pathway_name = c("pathway_name", "classification", "signaling"),
    score = c("score", "prob", "means", "mean", "prioritization_score", "LRscore", "lrscore", "magnitude"),
    pvalue = c("pvalue", "pval", "pvalues", "cellphone_pvals", "aggregate_rank", "specificity_rank", "magnitude_rank")
  )
  out <- df
  for (nm in names(rename_map)) {
    col <- ccc_pick_col(out, rename_map[[nm]])
    if (!is.null(col) && !nm %in% colnames(out)) {
      colnames(out)[match(col, colnames(out))] <- nm
    }
  }
  for (nm in c(
    "sender",
    "receiver",
    "interaction_name",
    "ligand",
    "receptor",
    "pathway_name"
  )) {
    if (!nm %in% colnames(out)) {
      out[[nm]] <- NA_character_
    }
  }
  interaction_label_col <- ccc_pick_col(
    out,
    c("interaction_label", "interaction_name_2", "interaction_name", "interacting_pair")
  )
  if (!is.null(interaction_label_col) && !"interaction_label" %in% colnames(out)) {
    out[["interaction_label"]] <- out[[interaction_label_col]]
  }
  if (!"interaction_label" %in% colnames(out)) {
    out[["interaction_label"]] <- out[["interaction_name"]]
  }
  out[["interaction_label"]] <- ccc_display_interaction(out[["interaction_label"]])
  if (!"classification" %in% colnames(out)) {
    out[["classification"]] <- out[["pathway_name"]]
  }
  out[["classification"]] <- as.character(out[["classification"]])
  out[["classification"]][is.na(out[["classification"]]) | !nzchar(out[["classification"]])] <- "Unclassified"
  out[["pathway_name"]] <- out[["classification"]]
  pair_lr_col <- ccc_pick_col(
    out,
    c("pair_lr", "pairLR", "interacting_pair", "interaction_name_2", "interaction_name")
  )
  if (!is.null(pair_lr_col) && !"pair_lr" %in% colnames(out)) {
    out[["pair_lr"]] <- out[[pair_lr_col]]
  }
  if (!"pair_lr" %in% colnames(out)) {
    out[["pair_lr"]] <- paste(out[["ligand"]], out[["receptor"]], sep = "-")
  }
  if (!"score" %in% colnames(out)) {
    out$score <- NA_real_
  }
  if (!"pvalue" %in% colnames(out)) {
    out$pvalue <- NA_real_
  }
  if (!"ligand_display" %in% colnames(out)) {
    out$ligand_display <- ccc_display_gene(out$ligand)
  }
  if (!"receptor_display" %in% colnames(out)) {
    out$receptor_display <- ccc_display_gene(out$receptor)
  }
  if (!"interaction_display" %in% colnames(out)) {
    out$interaction_display <- out$interaction_label
  }
  out <- ccc_mark_significance(out)
  out
}

ccc_display_gene <- function(x) {
  x <- ccc_clean_identifier(x, drop_receptor_suffix = TRUE)
  x[is.na(x)] <- ""
  x
}

ccc_display_interaction <- function(x) {
  x <- ccc_clean_identifier(x, drop_receptor_suffix = FALSE)
  x[is.na(x)] <- ""
  x <- gsub("_", " - ", x, fixed = TRUE)
  x <- gsub("\\s*-\\s*", " - ", x)
  x
}

ccc_clean_identifier <- function(x, drop_receptor_suffix = FALSE) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(x)
  x <- sub("^(complex:|simple:)", "", x, ignore.case = TRUE)
  x <- sub("(?:_|\\s)complex$", "", x, ignore.case = TRUE)
  if (isTRUE(drop_receptor_suffix)) {
    x <- sub("^integrin[_\\s]+", "", x, ignore.case = TRUE)
    x <- sub("_receptor_inhibitor$", "", x, ignore.case = TRUE)
    x <- sub("_receptor$", "", x, ignore.case = TRUE)
    x <- sub("_ligand$", "", x, ignore.case = TRUE)
  }
  x
}

ccc_standardize_ligand_target_df <- function(
  df,
  top_n = 20,
  sender.use = NULL,
  receiver.use = NULL,
  context_df = NULL,
  sender_default = NULL,
  receiver_default = NULL
) {
  df <- standardize_df(df)
  if (is.null(df) || nrow(df) == 0L) {
    log_message(
      "No ligand-target table is available for plotting",
      message_type = "error"
    )
  }
  df <- rename_by_candidates(
    df,
    "ligand",
    c("ligand", "from", "test_ligand", "ligand_oi", "ligand_source")
  )
  df <- rename_by_candidates(
    df,
    "target",
    c("target", "to", "gene", "target_gene", "target_genes", "gene_oi")
  )
  df <- rename_by_candidates(
    df,
    "sender",
    c("sender", "sender_celltype", "celltype_sender", "source_celltype")
  )
  df <- rename_by_candidates(
    df,
    "receiver",
    c("receiver", "receiver_celltype", "celltype_receiver", "target_celltype")
  )
  df <- rename_by_candidates(
    df,
    "weight",
    c(
      "weight",
      "score",
      "regulatory_potential",
      "pearson",
      "aupr_corrected",
      "aupr",
      "activity",
      "prioritization_score"
    )
  )

  if (!"weight" %in% colnames(df)) {
    numeric_cols <- setdiff(
      colnames(df)[vapply(df, is.numeric, logical(1))],
      c("ligand", "target")
    )
    if (length(numeric_cols) > 0L) {
      colnames(df)[match(numeric_cols[1], colnames(df))] <- "weight"
    }
  }

  if (!all(c("ligand", "target", "weight") %in% colnames(df))) {
    log_message(
      "Ligand-target table must contain ligand, target and weight columns",
      message_type = "error"
    )
  }

  sender_default <- unique(stats::na.omit(as.character(sender_default)))
  receiver_default <- unique(stats::na.omit(as.character(receiver_default)))
  if (!"sender" %in% colnames(df) && length(sender_default) == 1L) {
    df$sender <- sender_default
  }
  if (!"receiver" %in% colnames(df) && length(receiver_default) == 1L) {
    df$receiver <- receiver_default
  }

  context_df <- standardize_df(context_df)
  if (nrow(context_df) > 0L) {
    context_df <- standardize_long_df(context_df)
  }

  sender_requested <- !is.null(sender.use) && length(sender.use) > 0L
  receiver_requested <- !is.null(receiver.use) && length(receiver.use) > 0L
  direct_sender <- "sender" %in% colnames(df)
  direct_receiver <- "receiver" %in% colnames(df)

  if (sender_requested && direct_sender) {
    df <- df[df$sender %in% sender.use, , drop = FALSE]
  }
  if (receiver_requested && direct_receiver) {
    df <- df[df$receiver %in% receiver.use, , drop = FALSE]
  }

  if (
    (sender_requested && !direct_sender) ||
      (receiver_requested && !direct_receiver)
  ) {
    if (
      nrow(context_df) > 0L &&
        all(c("ligand", "sender", "receiver") %in% colnames(context_df))
    ) {
      context_keep <- context_df
      if (sender_requested) {
        context_keep <- context_keep[
          context_keep$sender %in% sender.use, ,
          drop = FALSE
        ]
      }
      if (receiver_requested) {
        context_keep <- context_keep[
          context_keep$receiver %in% receiver.use, ,
          drop = FALSE
        ]
      }
      keep_ligands <- unique(as.character(context_keep$ligand))
      keep_ligands <- keep_ligands[!is.na(keep_ligands) & nzchar(keep_ligands)]
      if (length(keep_ligands) > 0L) {
        df <- df[df$ligand %in% keep_ligands, , drop = FALSE]
      } else {
        df <- df[0, , drop = FALSE]
      }
    } else {
      log_message(
        paste0(
          "{.arg sender.use}/{.arg receiver.use} were provided, but the ",
          "ligand-target table does not contain sender/receiver columns and ",
          "no compatible long-table context is available for filtering."
        ),
        message_type = "warning"
      )
    }
  }

  if (nrow(df) == 0L) {
    log_message(
      "No ligand-target records remain after applying the current filters",
      message_type = "error"
    )
  }

  df <- df[order(df$weight, decreasing = TRUE), , drop = FALSE]
  df <- utils::head(
    df,
    max(top_n, 1L) * max(3L, min(10L, length(unique(df$ligand))))
  )
  ligand_levels <- unique(df$ligand)
  target_levels <- unique(df$target)
  if (length(ligand_levels) > top_n) {
    ligand_levels <- ligand_levels[seq_len(top_n)]
    df <- df[df$ligand %in% ligand_levels, , drop = FALSE]
    target_levels <- unique(df$target)
  }
  df$weight <- suppressWarnings(as.numeric(df$weight))
  df$ligand <- factor(df$ligand, levels = rev(ligand_levels))
  df$target <- factor(df$target, levels = target_levels)
  df
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
