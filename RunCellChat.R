
#' @title Run CellChat
#'
#' @description
#' [RunCellChat] runs CellChat on a Seurat object and stores the results in
#' `srt@tools[["CellChat"]]`. It supports either a single global analysis or
#' condition-specific analyses followed by pairwise merged comparisons.
#'
#'
#' @return
#'
#' @export
#'
#' @seealso
#' [CellChatPlot]
#'
#' @references
#' [CellChat](https://github.com/jinworks/CellChat),
#' [scDown::run_cellchatV2](https://htmlpreview.github.io/?https://raw.githubusercontent.com/BCH-RC/scDown/main/vignettes/scDown_CellChatV2.html)
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
    assay = NULL,
    layer = "data",
    verbose = TRUE
) {
  log_message(
    "Start {.pkg CellChat} analysis",
    verbose = verbose
  )

  check_r("jinworks/CellChat", verbose = FALSE)
  check_r("immunogenomics/presto", verbose = FALSE)

  .cc_validate_seurat_input(
    srt = srt,
    group.by = group.by,
    split.by = split.by,
    assay = assay %||% DefaultAssay(srt),
    layer = layer,
    annotation_selected = annotation_selected
  )

  assay <- assay %||% DefaultAssay(srt)
  species_use <- .cc_normalize_species(species)

  results <- list()
  comparisons <- list()

  if (is.null(group_column)) {
    res <- .cc_run_one_cellchat(
      seu = srt,
      label = "ALL",
      group.by = group.by,
      split.by = split.by,
      annotation_selected = annotation_selected,
      assay = assay,
      layer = layer,
      species = species_use,
      thresh = thresh,
      min.cells = min.cells,
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
      cells_i <- colnames(srt)[cond_vec == condition]
      seu_i <- srt[, cells_i, drop = FALSE]
      res_i <- .cc_run_one_cellchat(
        seu = seu_i,
        label = condition,
        group.by = group.by,
        split.by = split.by,
        annotation_selected = annotation_selected,
        assay = assay,
        layer = layer,
        species = species_use,
        thresh = thresh,
        min.cells = min.cells,
        verbose = verbose
      )
      if (!is.null(res_i)) {
        results[[condition]] <- res_i
      }
    }

    pairwise_cmp <- .cc_expand_pairwise_spec(
      cmp_spec = group_cmp,
      available_names = names(results),
      label = "group_cmp",
      verbose = verbose
    )

    comparisons <- .cc_build_comparison_results(
      result_list = results,
      pairwise_cmp = pairwise_cmp,
      verbose = verbose
    )
  }

  bundle <- .cc_make_cellchat_bundle(
    results = results,
    comparisons = comparisons,
    parameters = list(
      group.by = group.by,
      species = species_use,
      split.by = split.by,
      annotation_selected = annotation_selected,
      group_column = group_column,
      group_cmp = group_cmp,
      thresh = thresh,
      min.cells = min.cells,
      assay = assay,
      layer = layer
    )
  )

  srt@tools[["CellChat"]] <- bundle

  log_message(
    "CellChat analysis completed",
    message_type = "success",
    verbose = verbose
  )
  srt
}

.cc_validate_seurat_input <- function(
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
  if (length(group.by) != 1L || !is.character(group.by) || !group.by %in% colnames(srt@meta.data)) {
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

.cc_normalize_species <- function(species) {
  species <- match.arg(species, c("Homo_sapiens", "Mus_musculus", "zebrafish"))
  switch(
    species,
    "Mus_musculus" = "mouse",
    "Homo_sapiens" = "human",
    "zebrafish" = "zebrafish"
  )
}

.cc_expand_pairwise_spec <- function(
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
    names(out) <- vapply(out, function(x) paste0(x[1], "_vs_", x[2]), character(1))
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
      if (any(is.na(idx)) || any(idx < 1L) || any(idx > length(available_names))) {
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
      nm <- if (!is.null(raw_names) && nzchar(raw_names[i]) && length(pairs_i) == 1L) {
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

.cc_run_one_cellchat <- function(
    seu,
    label,
    group.by,
    split.by = NULL,
    annotation_selected = NULL,
    assay = NULL,
    layer = "data",
    species = c("human", "mouse", "zebrafish"),
    thresh = 0.05,
    min.cells = 10,
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
    keep_cells <- colnames(seu)[as.character(seu@meta.data[[group.by]]) %in% annotation_selected]
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

  cc <- .DoCellChat(
    object = seu,
    assay = assay,
    layer = layer,
    species = species,
    thresh = thresh,
    min.cells = min.cells
  )

  list(
    cellchat_object = cc,
    seurat_object = seu
  )
}

.cc_build_comparison_results <- function(
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

    merged_cc <- CellChat::mergeCellChat(
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

.cc_make_cellchat_bundle <- function(results, comparisons = list(), parameters = list()) {
  list(
    results = results,
    comparisons = comparisons,
    parameters = parameters
  )
}

.DoCellChat <- function(
    object,
    assay = NULL,
    layer = "data",
    species = c("human", "mouse", "zebrafish"),
    thresh = 0.05,
    min.cells = 10
) {
  assay <- assay %||% DefaultAssay(object)
  species <- match.arg(species)

  metadata <- data.frame(label = Idents(object))
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

  object <- CellChat::createCellChat(
    object = expr_mat,
    meta = metadata,
    group.by = "label"
  )

  object@DB <- switch(
    EXPR = species,
    "mouse" = CellChat::CellChatDB.mouse,
    "human" = CellChat::CellChatDB.human,
    "zebrafish" = CellChat::CellChatDB.zebrafish,
    log_message(
      "Invalid species. Must be one of {.val {c('human', 'mouse', 'zebrafish')}}",
      message_type = "error"
    )
  )

  object <- CellChat::subsetData(object)
  object <- CellChat::identifyOverExpressedGenes(
    object,
    thresh.p = thresh,
    min.cells = min.cells
  )
  object <- CellChat::identifyOverExpressedInteractions(object)
  object <- CellChat::computeCommunProb(object)
  object <- CellChat::filterCommunication(
    object,
    min.cells = min.cells
  )
  object <- CellChat::computeCommunProbPathway(object, thresh = thresh)
  object <- CellChat::aggregateNet(object, thresh = thresh)
  object <- CellChat::netAnalysis_computeCentrality(object, thresh = thresh)

  object
}

.SubsetCellChatMod <- function(
    object,
    idents.use,
    thresh = 0.05,
    verbose = TRUE
) {
  labels <- object@idents
  if (object@options$mode == "merged") {
    log_message(
      "Use the joint cell labels from the merged CellChat object",
      verbose = verbose
    )
    labels <- object@idents$joint
  }
  if (!is.factor(labels)) {
    labels <- factor(labels)
  }

  level_use0 <- levels(labels)
  level_use <- levels(labels)[levels(labels) %in% unique(labels)]
  level_use <- level_use[level_use %in% idents.use]
  level_use_index <- which(as.character(labels) %in% level_use)
  cells_use <- names(labels)[level_use_index]

  log_message(
    "The subset of cell groups used for CellChat analysis are {.val {level_use}}",
    verbose = verbose
  )

  data_subset <- object@data[, level_use_index, drop = FALSE]
  data_signaling_subset <- object@data.signaling[, level_use_index, drop = FALSE]
  meta_subset <- object@meta[level_use_index, , drop = FALSE]

  if (object@options$mode == "merged") {
    idents <- object@idents[1:(length(object@idents) - 1)]
    group_existing_index <- which(level_use0 %in% level_use)

    net_subset <- vector("list", length = length(object@net))
    netP_subset <- vector("list", length = length(object@netP))
    idents_subset <- vector("list", length = length(idents))
    names(net_subset) <- names(object@net)
    names(netP_subset) <- names(object@netP)
    names(idents_subset) <- names(object@idents[1:(length(object@idents) - 1)])

    images_subset <- vector("list", length = length(idents))
    names(images_subset) <- names(object@idents[1:(length(object@idents) - 1)])

    for (i in seq_along(idents)) {
      log_message(
        "Update slots object@images, object@net, object@netP, object@idents in dataset {.val {names(object@idents)[i]}}",
        verbose = verbose
      )
      images <- object@images[[i]]
      for (images_j in names(images)) {
        values <- images[[images_j]]
        if (images_j %in% c("coordinates")) {
          images[[images_j]] <- values[level_use_index, , drop = FALSE]
        }
        if (images_j %in% c("distance")) {
          images[[images_j]] <- values[group_existing_index, group_existing_index, drop = FALSE]
        }
      }
      images_subset[[i]] <- images

      net <- object@net[[i]]
      for (net.j in names(net)) {
        values <- net[[net.j]]
        if (net.j %in% c("prob", "pval")) {
          net[[net.j]] <- values[group_existing_index, group_existing_index, , drop = FALSE]
        }
        if (net.j %in% c("count", "sum", "weight")) {
          net[[net.j]] <- values[group_existing_index, group_existing_index, drop = FALSE]
        }
      }
      net_subset[[i]] <- net

      netP <- CellChat::computeCommunProbPathway(
        net = net_subset[[i]],
        pairLR.use = object@LR[[i]]$LRsig,
        thresh = thresh
      )
      netP$centr <- CellChat::netAnalysis_computeCentrality(
        net = net_subset[[i]]$prob
      )
      netP_subset[[i]] <- netP

      idents_subset[[i]] <- idents[[i]][names(idents[[i]]) %in% cells_use]
      idents_subset[[i]] <- factor(
        idents_subset[[i]],
        levels = levels(idents[[i]])[levels(idents[[i]]) %in% level_use]
      )
    }

    idents_subset$joint <- factor(
      object@idents$joint[level_use_index],
      levels = level_use
    )
  } else {
    group_existing_index <- which(level_use0 %in% level_use)

    images <- object@images
    for (images_j in names(images)) {
      values <- images[[images_j]]
      if (images_j %in% c("coordinates")) {
        images[[images_j]] <- values[level_use_index, , drop = FALSE]
      }
      if (images_j %in% c("distance")) {
        images[[images_j]] <- values[group_existing_index, group_existing_index, drop = FALSE]
      }
    }
    images_subset <- images

    net <- object@net
    for (net.j in names(net)) {
      values <- net[[net.j]]
      if (net.j %in% c("prob", "pval")) {
        net[[net.j]] <- values[group_existing_index, group_existing_index, , drop = FALSE]
      }
      if (net.j %in% c("count", "sum", "weight")) {
        net[[net.j]] <- values[group_existing_index, group_existing_index, drop = FALSE]
      }
    }
    net_subset <- net

    netP <- CellChat::computeCommunProbPathway(
      net = net_subset,
      pairLR.use = object@LR$LRsig,
      thresh = thresh
    )
    netP$centr <- CellChat::netAnalysis_computeCentrality(net = net_subset$prob)
    netP_subset <- netP

    idents_subset <- object@idents[level_use_index]
    idents_subset <- factor(idents_subset, levels = level_use)
  }

  methods::new(
    Class = "CellChat",
    data = data_subset,
    data.signaling = data_signaling_subset,
    images = images_subset,
    net = net_subset,
    netP = netP_subset,
    meta = meta_subset,
    idents = idents_subset,
    var.features = object@var.features,
    LR = object@LR,
    DB = object@DB,
    options = object@options
  )
}