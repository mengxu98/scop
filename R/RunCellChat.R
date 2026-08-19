#' @title Run CellChat analysis
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams scop-params
#' @inheritParams CellDimPlot
#' @param species `"Homo_sapiens"`, `"Mus_musculus"`, or `"zebrafish"`.
#' @param annotation_selected Cell types to include. `NULL` uses all.
#' @param group_column Metadata column defining conditions or groups.
#' @param group_cmp Pairwise condition comparisons for differential CellChat.
#' @param thresh Threshold for centrality scores.
#' @param min.cells Minimum expressed cells required for genes used in CCC.
#' @param do.fast Use CellChat's fast Wilcoxon via `presto` (must be installed).
#' @param backend Post-processing / unified CCC table backend. Does not change
#' upstream CellChat inference.
#' @return A `Seurat` object with `CellChat` results stored in `srt@tools[["CellChat"]]`.
#'
#' @export
#'
#' @seealso
#' [CCCHeatmap], [CCCStatPlot], [CCCNetworkPlot]
#'
#' @references
#' [CellChat](https://github.com/jinworks/CellChat)
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#' pancreas_sub <- RunCellChat(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   species = "Mus_musculus"
#' )
#'
#' CCCNetworkPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   plot_type = "bipartite"
#' )
#'
#' CCCHeatmap(
#'   pancreas_sub,
#'   method = "CellChat",
#'   plot_type = "heatmap"
#' )
#'
#' CCCStatPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   plot_type = "violin",
#'   top_n = 50
#' )
RunCellChat <- function(
  srt,
  group.by,
  species = c("Homo_sapiens", "Mus_musculus", "zebrafish"),
  split.by = NULL,
  annotation_selected = NULL,
  group_column = NULL,
  group_cmp = NULL,
  thresh = 0.05,
  min.cells = 10,
  do.fast = FALSE,
  backend = c("cpp", "r"),
  assay = NULL,
  layer = "data",
  verbose = TRUE
) {
  backend <- match.arg(backend)
  log_message(
    "Start {.pkg CellChat} analysis",
    verbose = verbose
  )

  check_r("jinworks/CellChat", verbose = FALSE)
  identifyOverExpressedGenes <- get_namespace_fun("CellChat", "identifyOverExpressedGenes")
  identify_formals <- names(formals(identifyOverExpressedGenes))
  if (
    isTRUE(do.fast) &&
      ("..." %in% identify_formals || "do.fast" %in% identify_formals)
  ) {
    check_r("immunogenomics/presto", verbose = FALSE)
  }

  validate_cc_input(
    srt = srt,
    group.by = group.by,
    split.by = split.by,
    assay = assay %||% DefaultAssay(srt),
    layer = layer,
    annotation_selected = annotation_selected
  )

  assay <- assay %||% DefaultAssay(srt)
  species <- match.arg(species)
  srt_cc <- Seurat::DietSeurat(
    object = srt,
    assays = assay,
    dimreducs = NULL,
    graphs = NULL,
    misc = FALSE
  )
  SeuratObject::DefaultAssay(srt_cc) <- assay

  results <- list()
  comparisons <- list()

  if (is.null(group_column)) {
    res <- run_one_cc(
      seu = srt_cc,
      label = "ALL",
      group.by = group.by,
      split.by = split.by,
      annotation_selected = annotation_selected,
      assay = assay,
      layer = layer,
      species = species,
      thresh = thresh,
      min.cells = min.cells,
      do.fast = do.fast,
      verbose = verbose
    )
    if (!is.null(res)) {
      results[["ALL"]] <- res
    }
  } else {
    if (!group_column %in% colnames(srt@meta.data)) {
      log_message(
        "{.val {group_column}} does not exist in {.cls Seurat}",
        message_type = "error"
      )
    }

    cond_vec <- as.character(srt@meta.data[[group_column]])
    conditions <- unique(cond_vec)
    conditions <- conditions[!is.na(conditions)]

    for (condition in conditions) {
      log_message(
        "Processing condition: {.val {condition}}",
        verbose = verbose
      )
      cells_i <- colnames(srt_cc)[cond_vec == condition]
      seu_i <- srt_cc[, cells_i, drop = FALSE]
      res_i <- run_one_cc(
        seu = seu_i,
        label = condition,
        group.by = group.by,
        split.by = split.by,
        annotation_selected = annotation_selected,
        assay = assay,
        layer = layer,
        species = species,
        thresh = thresh,
        min.cells = min.cells,
        do.fast = do.fast,
        verbose = verbose
      )
      if (!is.null(res_i)) {
        results[[condition]] <- res_i
      }
    }

    pairwise_cmp <- expand_pairwise_spec(
      cmp_spec = group_cmp,
      available_names = names(results),
      label = "group_cmp",
      verbose = verbose
    )

    comparisons <- build_comparison_results(
      result_list = results,
      pairwise_cmp = pairwise_cmp,
      verbose = verbose
    )
  }

  bundle <- list(
    results = results,
    comparisons = comparisons,
    parameters = list(
      group.by = group.by,
      species = species,
      split.by = split.by,
      annotation_selected = annotation_selected,
      group_column = group_column,
      group_cmp = group_cmp,
      thresh = thresh,
      min.cells = min.cells,
      do.fast = do.fast,
      backend = backend,
      backend_scope = "result aggregation and unified-table construction",
      assay = assay,
      layer = layer
    )
  )

  srt@tools[["CellChat"]] <- bundle
  cellchat_long <- ccc_build_cellchat_long_table(srt, thresh = thresh)
  bundle$long_table <- cellchat_long
  bundle$primary_table <- cellchat_long
  bundle$pair_table <- aggregate_ccc_long(cellchat_long, backend = backend)
  bundle$metadata <- list(
    updated_at = as.character(Sys.time()),
    backend = backend
  )
  bundle$provenance <- list(
    producer = "RunCellChat",
    backend_version = as.character(utils::packageVersion("CellChat"))
  )
  srt@tools[["CellChat"]] <- bundle
  srt <- ccc_update_unified_bundle(
    srt = srt,
    method = "CellChat",
    bundle = bundle,
    thresh = thresh,
    backend = backend
  )

  log_message(
    "{.pkg CellChat} analysis completed",
    message_type = "success",
    verbose = verbose
  )
  srt
}
validate_cc_input <- function(
  srt,
  group.by,
  split.by = NULL,
  assay = NULL,
  layer = "data",
  annotation_selected = NULL
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (
    length(group.by) != 1L ||
      !is.character(group.by) ||
      !group.by %in% colnames(srt@meta.data)
  ) {
    log_message(
      "{.arg group.by} must be a valid metadata column in {.cls Seurat}",
      message_type = "error"
    )
  }
  if (!is.null(split.by) && !split.by %in% colnames(srt@meta.data)) {
    log_message(
      "{.val {split.by}} does not exist in {.cls Seurat}",
      message_type = "error"
    )
  }
  assay <- assay %||% DefaultAssay(srt)
  if (!assay %in% names(srt@assays)) {
    log_message(
      "{.val {assay}} does not exist in {.cls Seurat}",
      message_type = "error"
    )
  }
  if (!is.character(layer) || length(layer) != 1L || is.na(layer)) {
    log_message(
      "{.arg layer} must be a single non-missing character string",
      message_type = "error"
    )
  }

  if (!is.null(annotation_selected)) {
    available_annotations <- unique(as.character(srt@meta.data[[group.by]]))
    missing_annotations <- setdiff(annotation_selected, available_annotations)
    if (length(missing_annotations) > 0L) {
      log_message(
        "Missing annotations in {.val {group.by}}: {.val {missing_annotations}}",
        message_type = "error"
      )
    }
  }

  invisible(TRUE)
}

expand_pairwise_spec <- function(
  cmp_spec = NULL,
  available_names,
  label = "group_cmp",
  verbose = TRUE
) {
  available_names <- unique(as.character(stats::na.omit(available_names)))
  if (length(available_names) < 2L) {
    return(list())
  }

  if (is.null(cmp_spec)) {
    out <- utils::combn(available_names, 2, simplify = FALSE)
    names(out) <- vapply(
      out,
      function(x) paste0(x[1], "_vs_", x[2]),
      character(1)
    )
    return(out)
  }

  if (!is.list(cmp_spec)) {
    cmp_spec <- list(cmp_spec)
  }

  raw_names <- names(cmp_spec)
  out <- list()
  out_names <- character()

  for (i in seq_along(cmp_spec)) {
    cmp_i <- cmp_spec[[i]]

    if (is.numeric(cmp_i)) {
      idx <- as.integer(cmp_i)
      if (
        any(is.na(idx)) || any(idx < 1L) || any(idx > length(available_names))
      ) {
        log_message(
          "{.arg {label}} entry {.val {i}} contains indices outside the available range",
          message_type = "error"
        )
      }
      cmp_i <- available_names[idx]
    } else {
      cmp_i <- as.character(cmp_i)
    }

    cmp_i <- unique(stats::na.omit(cmp_i))
    if (length(cmp_i) < 2L) {
      log_message(
        "Each entry of {.arg {label}} must contain at least two groups. Entry {.val {i}} was skipped",
        message_type = "warning",
        verbose = verbose
      )
      next
    }

    missing_groups <- setdiff(cmp_i, available_names)
    if (length(missing_groups) > 0L) {
      log_message(
        "Missing groups in {.arg {label}}: {.val {missing_groups}}",
        message_type = "error"
      )
    }

    pairs_i <- if (length(cmp_i) == 2L) {
      list(cmp_i)
    } else {
      log_message(
        "Entry {.val {i}} of {.arg {label}} contains more than two groups; all pairwise comparisons will be created",
        verbose = verbose
      )
      utils::combn(cmp_i, 2, simplify = FALSE)
    }

    for (j in seq_along(pairs_i)) {
      pair_j <- as.character(pairs_i[[j]])
      nm <- if (
        !is.null(raw_names) && nzchar(raw_names[i]) && length(pairs_i) == 1L
      ) {
        raw_names[i]
      } else {
        paste0(pair_j[1], "_vs_", pair_j[2])
      }
      out[[length(out) + 1L]] <- pair_j
      out_names <- c(out_names, nm)
    }
  }

  if (length(out) == 0L) {
    return(list())
  }

  names(out) <- make.unique(out_names, sep = "_")
  out
}

run_one_cc <- function(
  seu,
  label,
  group.by,
  split.by = NULL,
  annotation_selected = NULL,
  assay = NULL,
  layer = "data",
  species,
  thresh = 0.05,
  min.cells = 10,
  do.fast = FALSE,
  verbose = TRUE
) {
  if (!inherits(seu, "Seurat")) {
    log_message(
      "{.arg seu} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }

  assay <- assay %||% DefaultAssay(seu)
  Idents(seu) <- seu@meta.data[[group.by]]

  if (!is.null(split.by)) {
    seu$samples <- as.factor(seu@meta.data[[split.by]])
  } else {
    seu$samples <- as.factor("All")
  }

  if (!is.null(annotation_selected)) {
    keep_cells <- colnames(seu)[
      as.character(seu@meta.data[[group.by]]) %in% annotation_selected
    ]
    if (length(keep_cells) == 0L) {
      log_message(
        "No cells retained for {.val {label}} after filtering {.arg annotation_selected}",
        message_type = "warning",
        verbose = verbose
      )
      return(NULL)
    }
    seu <- seu[, keep_cells, drop = FALSE]
    Idents(seu) <- seu@meta.data[[group.by]]
  }

  if (ncol(seu) == 0L) {
    log_message(
      "No cells found for {.val {label}}",
      message_type = "warning",
      verbose = verbose
    )
    return(NULL)
  }

  cc <- DoCellChat(
    object = seu,
    assay = assay,
    layer = layer,
    species = species,
    thresh = thresh,
    min.cells = min.cells,
    do.fast = do.fast
  )

  list(
    cellchat_object = cc,
    seurat_object = seu
  )
}

build_comparison_results <- function(
  result_list,
  pairwise_cmp,
  verbose = TRUE
) {
  if (length(pairwise_cmp) == 0L) {
    return(list())
  }

  out <- list()
  for (i in seq_along(pairwise_cmp)) {
    pair_i <- as.character(pairwise_cmp[[i]])
    nm_i <- names(pairwise_cmp)[i] %||% paste0(pair_i[1], "_vs_", pair_i[2])

    if (length(pair_i) != 2L) {
      log_message(
        "Comparison {.val {nm_i}} skipped because it does not contain exactly two datasets/groups",
        message_type = "warning",
        verbose = verbose
      )
      next
    }
    if (!all(pair_i %in% names(result_list))) {
      log_message(
        "Comparison {.val {nm_i}} skipped because one or more groups were not successfully analyzed",
        message_type = "warning",
        verbose = verbose
      )
      next
    }

    log_message(
      "Merging CellChat objects for comparison {.val {nm_i}}",
      verbose = verbose
    )

    object.list <- lapply(pair_i, function(x) result_list[[x]]$cellchat_object)
    names(object.list) <- pair_i

    merged_cc <- get_namespace_fun("CellChat", "mergeCellChat")(
      object.list = object.list,
      add.names = names(object.list)
    )

    out[[nm_i]] <- list(
      comparison_object = merged_cc,
      object.list = object.list,
      groups = pair_i
    )
  }

  out
}

DoCellChat <- function(
  object,
  assay = NULL,
  layer = "data",
  species,
  thresh = 0.05,
  min.cells = 10,
  do.fast = FALSE
) {
  assay <- assay %||% DefaultAssay(object)

  cell_labels <- as.character(Idents(object))
  if (any(cell_labels == "0")) {
    replacement <- "C0"
    while (replacement %in% cell_labels) {
      replacement <- paste0(replacement, "_")
    }
    cell_labels[cell_labels == "0"] <- replacement
    warning(
      sprintf(
        "Cell labels contain the CellChat-reserved value \"0\"; remapped to \"%s\"",
        replacement
      ),
      call. = FALSE
    )
  }
  metadata <- data.frame(label = cell_labels)
  metadata <- cbind(metadata, object@meta.data)

  expr_mat <- tryCatch(
    GetAssayData5(object = object, layer = layer, assay = assay),
    error = function(e) NULL
  )
  if (is.null(expr_mat)) {
    log_message(
      "Failed to extract expression data from assay {.val {assay}} and layer {.val {layer}}. Please ensure the layer exists and is normalized for CellChat",
      message_type = "error"
    )
  }

  object <- get_namespace_fun("CellChat", "createCellChat")(
    object = expr_mat,
    meta = metadata,
    group.by = "label"
  )

  object@DB <- cellchat_database(species = species)

  object <- get_namespace_fun("CellChat", "subsetData")(object)
  identifyOverExpressedGenes <- get_namespace_fun("CellChat", "identifyOverExpressedGenes")
  identify_args <- list(
    object = object,
    thresh.p = thresh,
    min.cells = min.cells,
    do.fast = isTRUE(do.fast)
  )
  identify_formals <- names(formals(identifyOverExpressedGenes))
  if (!"..." %in% identify_formals) {
    identify_args <- identify_args[names(identify_args) %in% identify_formals]
  }
  object <- do.call(identifyOverExpressedGenes, identify_args)
  object <- get_namespace_fun("CellChat", "identifyOverExpressedInteractions")(object)
  object <- get_namespace_fun("CellChat", "computeCommunProb")(object)
  object <- get_namespace_fun("CellChat", "filterCommunication")(
    object,
    min.cells = min.cells
  )
  object <- get_namespace_fun("CellChat", "computeCommunProbPathway")(object, thresh = thresh)
  object <- get_namespace_fun("CellChat", "aggregateNet")(object, thresh = thresh)
  object <- get_namespace_fun("CellChat", "netAnalysis_computeCentrality")(object, thresh = thresh)

  object
}

cellchat_database <- function(species) {
  data_name <- switch(species,
    "Mus_musculus" = "CellChatDB.mouse",
    "Homo_sapiens" = "CellChatDB.human",
    "zebrafish" = "CellChatDB.zebrafish",
    log_message(
      "Invalid species. Must be one of {.val {c('Homo_sapiens', 'Mus_musculus', 'zebrafish')}}",
      message_type = "error"
    )
  )

  db <- get_namespace_fun("CellChat", data_name)
  if (!is.null(db)) {
    return(db)
  }

  data_env <- new.env(parent = emptyenv())
  utils::data(list = data_name, package = "CellChat", envir = data_env)
  db <- get0(data_name, envir = data_env, inherits = FALSE)
  if (!is.null(db)) {
    return(db)
  }

  log_message(
    "Cannot load {.pkg CellChat} database object {.val {data_name}}. Update {.pkg CellChat} or report the output of {.code utils::data(package = 'CellChat')}.",
    message_type = "error"
  )
}
