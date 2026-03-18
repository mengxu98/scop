## =========================
## Enhanced self-contained CellChat module
## =========================

## ---------- helpers ----------
.msg <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) {
    message(sprintf(...))
  }
}

.stop2 <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

.check_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    .stop2("Package '%s' is required but not installed.", pkg)
  }
}

.get_assay_data <- function(object, assay = "RNA", layer = "data") {
  if (!inherits(object, "Seurat")) {
    .stop2("object must be a Seurat object")
  }
  
  out <- tryCatch(
    SeuratObject::GetAssayData(object = object, assay = assay, layer = layer),
    error = function(e) NULL
  )
  if (!is.null(out)) return(out)
  
  out <- tryCatch(
    Seurat::GetAssayData(object = object, assay = assay, slot = layer),
    error = function(e) NULL
  )
  if (!is.null(out)) return(out)
  
  .stop2("Failed to extract assay data from assay='%s', layer/slot='%s'.", assay, layer)
}

.save_plot_if_needed <- function(plot_obj, file = NULL, width = 8, height = 6, dpi = 300) {
  if (is.null(file)) return(invisible(plot_obj))
  
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  
  ext <- tolower(tools::file_ext(file))
  if (inherits(plot_obj, c("gg", "ggplot"))) {
    ggplot2::ggsave(filename = file, plot = plot_obj, width = width, height = height, dpi = dpi)
  } else if (ext %in% c("pdf")) {
    grDevices::pdf(file, width = width, height = height)
    print(plot_obj)
    grDevices::dev.off()
  } else if (ext %in% c("png")) {
    grDevices::png(file, width = width, height = height, units = "in", res = dpi)
    print(plot_obj)
    grDevices::dev.off()
  } else {
    warning("Unsupported auto-save format for non-ggplot object: ", file)
  }
  
  invisible(plot_obj)
}

.resolve_group_index <- function(object, groups = NULL) {
  if (is.null(groups)) return(NULL)
  
  lev <- levels(object@idents)
  if (is.numeric(groups)) return(groups)
  
  groups <- as.character(groups)
  idx <- match(groups, lev)
  if (any(is.na(idx))) {
    .stop2(
      "These cell groups were not found: %s. Available groups: %s",
      paste(groups[is.na(idx)], collapse = ", "),
      paste(lev, collapse = ", ")
    )
  }
  idx
}

.resolve_group_index_merged <- function(merged_object, groups = NULL) {
  if (is.null(groups)) return(NULL)
  
  if (!is.list(merged_object@idents) || is.null(merged_object@idents$joint)) {
    return(.resolve_group_index(merged_object, groups))
  }
  
  lev <- levels(merged_object@idents$joint)
  if (is.numeric(groups)) return(groups)
  
  groups <- as.character(groups)
  idx <- match(groups, lev)
  if (any(is.na(idx))) {
    .stop2(
      "These cell groups were not found in merged object: %s. Available groups: %s",
      paste(groups[is.na(idx)], collapse = ", "),
      paste(lev, collapse = ", ")
    )
  }
  idx
}

.resolve_comparison_dataset <- function(cmp_obj, comparison = c(1, 2)) {
  ds_names <- names(cmp_obj$object.list)
  if (is.numeric(comparison)) return(comparison)
  
  comparison <- as.character(comparison)
  idx <- match(comparison, ds_names)
  if (any(is.na(idx))) {
    .stop2(
      "comparison names not found: %s. Available datasets: %s",
      paste(comparison[is.na(idx)], collapse = ", "),
      paste(ds_names, collapse = ", ")
    )
  }
  idx
}

## ---------- core ----------
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
    assay = "RNA",
    verbose = TRUE
) {
  .check_pkg("Seurat")
  .check_pkg("CellChat")
  
  .msg("Start CellChat analysis", verbose = verbose)
  
  if (!inherits(srt, "Seurat")) .stop2("srt must be a Seurat object")
  if (!group.by %in% colnames(srt@meta.data)) .stop2("group.by '%s' does not exist in srt@meta.data", group.by)
  if (!is.null(split.by) && !split.by %in% colnames(srt@meta.data)) .stop2("split.by '%s' does not exist in srt@meta.data", split.by)
  if (!is.null(group_column) && !group_column %in% colnames(srt@meta.data)) .stop2("group_column '%s' does not exist in srt@meta.data", group_column)
  if (!assay %in% names(srt@assays)) .stop2("assay '%s' does not exist in Seurat object", assay)
  
  if (!is.null(annotation_selected)) {
    available_annotations <- unique(as.character(srt@meta.data[[group.by]]))
    missing_annotations <- setdiff(annotation_selected, available_annotations)
    if (length(missing_annotations) > 0) {
      .stop2("Missing annotations in '%s': %s", group.by, paste(missing_annotations, collapse = ", "))
    }
  }
  
  species <- match.arg(species)
  species <- switch(
    species,
    "Mus_musculus" = "mouse",
    "Homo_sapiens" = "human",
    "zebrafish" = "zebrafish"
  )
  
  Seurat::Idents(srt) <- srt@meta.data[[group.by]]
  
  if (!is.null(split.by)) {
    srt$samples <- as.factor(srt@meta.data[[split.by]])
  } else {
    srt$samples <- as.factor("All")
  }
  
  if (is.null(srt@tools[["CellChat"]])) srt@tools[["CellChat"]] <- list()
  
  cellchat_results <- list()
  comparison_results <- list()
  
  .run_one <- function(seu_use, label) {
    if (ncol(seu_use) == 0) {
      .msg("No cells found in %s", label, verbose = verbose)
      return(NULL)
    }
    
    Seurat::Idents(seu_use) <- seu_use@meta.data[[group.by]]
    
    cc <- .DoCellChat(
      object = seu_use,
      assay = assay,
      species = species,
      thresh = thresh,
      min.cells = min.cells
    )
    
    if (!is.null(annotation_selected)) {
      cc <- .SubsetCellChatMod(
        object = cc,
        idents.use = annotation_selected,
        thresh = thresh,
        verbose = verbose
      )
      cc <- CellChat::netAnalysis_computeCentrality(cc, thresh = thresh)
      
      keep_cells <- colnames(seu_use)[as.character(seu_use@meta.data[[group.by]]) %in% annotation_selected]
      seu_use <- seu_use[, keep_cells, drop = FALSE]
    }
    
    list(
      cellchat_object = cc,
      seurat_object = seu_use
    )
  }
  
  if (is.null(group_column)) {
    cellchat_results[["ALL"]] <- .run_one(srt, "ALL")
  } else {
    conditions <- unique(as.character(srt@meta.data[[group_column]]))
    conditions <- conditions[!is.na(conditions)]
    
    if (!is.null(group_cmp)) {
      available_groups <- unique(as.character(srt@meta.data[[group_column]]))
      all_groups_in_cmp <- unique(unlist(group_cmp))
      missing_groups <- setdiff(all_groups_in_cmp, available_groups)
      if (length(missing_groups) > 0) {
        .stop2("Missing groups in group_cmp: %s", paste(missing_groups, collapse = ", "))
      }
    }
    
    for (condition in conditions) {
      .msg("Processing group: %s", condition, verbose = verbose)
      cells_use <- rownames(srt@meta.data)[as.character(srt@meta.data[[group_column]]) == condition]
      srt_condition <- srt[, cells_use, drop = FALSE]
      res_condition <- .run_one(srt_condition, condition)
      if (!is.null(res_condition)) cellchat_results[[condition]] <- res_condition
    }
    
    if (!is.null(group_cmp)) {
      for (cmp in group_cmp) {
        if (length(cmp) != 2) {
          warning("Each element of group_cmp must contain exactly 2 groups. Skipped.")
          next
        }
        
        g1 <- cmp[1]
        g2 <- cmp[2]
        cmp_name <- paste0(g1, "_vs_", g2)
        
        if (!all(c(g1, g2) %in% names(cellchat_results))) {
          warning(sprintf("Comparison %s skipped because group is missing.", cmp_name))
          next
        }
        
        .msg("Merging comparison: %s", cmp_name, verbose = verbose)
        
        object.list <- list(
          cellchat_results[[g1]]$cellchat_object,
          cellchat_results[[g2]]$cellchat_object
        )
        names(object.list) <- c(g1, g2)
        
        merged_cc <- CellChat::mergeCellChat(object.list = object.list, add.names = names(object.list))
        
        comparison_results[[cmp_name]] <- list(
          comparison_object = merged_cc,
          groups = c(g1, g2),
          object.list = object.list
        )
      }
    }
  }
  
  srt@tools[["CellChat"]][["results"]] <- cellchat_results
  srt@tools[["CellChat"]][["comparisons"]] <- comparison_results
  srt@tools[["CellChat"]][["parameters"]] <- list(
    group.by = group.by,
    species = species,
    split.by = split.by,
    annotation_selected = annotation_selected,
    group_column = group_column,
    group_cmp = group_cmp,
    thresh = thresh,
    min.cells = min.cells,
    assay = assay
  )
  
  .msg("CellChat analysis completed", verbose = verbose)
  srt
}

RunCellChatMulti <- function(
    seurat_list,
    group.by,
    dataset_names = names(seurat_list),
    species = c("Homo_sapiens", "Mus_musculus", "zebrafish"),
    split.by = NULL,
    annotation_selected = NULL,
    thresh = 0.05,
    min.cells = 10,
    assay = "RNA",
    verbose = TRUE
) {
  .check_pkg("CellChat")
  if (!is.list(seurat_list) || length(seurat_list) < 2) {
    .stop2("seurat_list must be a list with at least 2 Seurat objects")
  }
  
  if (is.null(dataset_names)) {
    dataset_names <- paste0("dataset", seq_along(seurat_list))
  }
  if (length(dataset_names) != length(seurat_list)) {
    .stop2("dataset_names length must equal seurat_list length")
  }
  
  cc_results <- vector("list", length(seurat_list))
  names(cc_results) <- dataset_names
  
  for (i in seq_along(seurat_list)) {
    .msg("Running CellChat for dataset: %s", dataset_names[i], verbose = verbose)
    tmp <- RunCellChat(
      srt = seurat_list[[i]],
      group.by = group.by,
      species = species,
      split.by = split.by,
      annotation_selected = annotation_selected,
      group_column = NULL,
      group_cmp = NULL,
      thresh = thresh,
      min.cells = min.cells,
      assay = assay,
      verbose = verbose
    )
    cc_results[[i]] <- tmp@tools[["CellChat"]][["results"]][["ALL"]]$cellchat_object
  }
  
  merged_cc <- CellChat::mergeCellChat(object.list = cc_results, add.names = names(cc_results))
  
  list(
    object.list = cc_results,
    comparison_object = merged_cc,
    dataset_names = dataset_names,
    parameters = list(
      group.by = group.by,
      species = match.arg(species),
      split.by = split.by,
      annotation_selected = annotation_selected,
      thresh = thresh,
      min.cells = min.cells,
      assay = assay
    )
  )
}

.DoCellChat <- function(
    object,
    assay = "RNA",
    species = c("human", "mouse", "zebrafish"),
    thresh = 0.05,
    min.cells = 10
) {
  .check_pkg("CellChat")
  species <- match.arg(species)
  
  metadata <- data.frame(label = Seurat::Idents(object))
  metadata <- cbind(metadata, object@meta.data)
  
  expr_mat <- .get_assay_data(object, assay = assay, layer = "data")
  
  cellchat_object <- CellChat::createCellChat(
    object = expr_mat,
    meta = metadata,
    group.by = "label"
  )
  
  cellchat_object@DB <- switch(
    species,
    "mouse" = CellChat::CellChatDB.mouse,
    "human" = CellChat::CellChatDB.human,
    "zebrafish" = CellChat::CellChatDB.zebrafish
  )
  
  cellchat_object <- CellChat::subsetData(cellchat_object)
  cellchat_object <- CellChat::identifyOverExpressedGenes(cellchat_object, thresh.p = thresh, min.cells = min.cells)
  cellchat_object <- CellChat::identifyOverExpressedInteractions(cellchat_object)
  cellchat_object <- CellChat::computeCommunProb(cellchat_object)
  cellchat_object <- CellChat::filterCommunication(cellchat_object, min.cells = min.cells)
  cellchat_object <- CellChat::computeCommunProbPathway(cellchat_object, thresh = thresh)
  cellchat_object <- CellChat::aggregateNet(cellchat_object, thresh = thresh)
  cellchat_object <- CellChat::netAnalysis_computeCentrality(cellchat_object, thresh = thresh)
  
  cellchat_object
}

.SubsetCellChatMod <- function(
    object,
    idents.use,
    thresh = 0.05,
    verbose = TRUE
) {
  labels <- object@idents
  
  if (object@options$mode == "merged") {
    .msg("Use the joint cell labels from the merged CellChat object", verbose = verbose)
    labels <- object@idents$joint
  }
  
  if (!is.factor(labels)) labels <- factor(labels)
  
  level_use0 <- levels(labels)
  level_use <- levels(labels)[levels(labels) %in% unique(labels)]
  level_use <- level_use[level_use %in% idents.use]
  level_use_index <- which(as.character(labels) %in% level_use)
  cells_use <- names(labels)[level_use_index]
  
  .msg("Subset cell groups: %s", paste(level_use, collapse = ", "), verbose = verbose)
  
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
      netP$centr <- CellChat::netAnalysis_computeCentrality(net = net_subset[[i]]$prob)
      netP_subset[[i]] <- netP
      
      idents_subset[[i]] <- idents[[i]][names(idents[[i]]) %in% cells_use]
      idents_subset[[i]] <- factor(
        idents_subset[[i]],
        levels = levels(idents[[i]])[levels(idents[[i]]) %in% level_use]
      )
    }
    
    idents_subset$joint <- factor(object@idents$joint[level_use_index], levels = level_use)
    
  } else {
    group_existing_index <- which(level_use0 %in% level_use)
    
    images <- object@images
    for (images_j in names(images)) {
      values <- images[[images_j]]
      if (images_j %in% c("coordinates")) images[[images_j]] <- values[level_use_index, , drop = FALSE]
      if (images_j %in% c("distance")) images[[images_j]] <- values[group_existing_index, group_existing_index, drop = FALSE]
    }
    images_subset <- images
    
    net <- object@net
    for (net.j in names(net)) {
      values <- net[[net.j]]
      if (net.j %in% c("prob", "pval")) net[[net.j]] <- values[group_existing_index, group_existing_index, , drop = FALSE]
      if (net.j %in% c("count", "sum", "weight")) net[[net.j]] <- values[group_existing_index, group_existing_index, drop = FALSE]
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

ListCellChatResults <- function(srt) {
  if (is.null(srt@tools[["CellChat"]])) {
    return(list(results = character(0), comparisons = character(0)))
  }
  
  list(
    results = names(srt@tools[["CellChat"]][["results"]]),
    comparisons = names(srt@tools[["CellChat"]][["comparisons"]])
  )
}

GetCellChatResult <- function(srt, condition = "ALL") {
  if (is.null(srt@tools[["CellChat"]][["results"]][[condition]])) {
    stop("CellChat result not found for condition: ", condition, call. = FALSE)
  }
  srt@tools[["CellChat"]][["results"]][[condition]]
}

GetCellChatComparison <- function(srt, condition) {
  if (is.null(srt@tools[["CellChat"]][["comparisons"]][[condition]])) {
    stop("CellChat comparison not found: ", condition, call. = FALSE)
  }
  srt@tools[["CellChat"]][["comparisons"]][[condition]]
}

CellChatPlot <- function(
    srt,
    plot_type = c(
      "overview",
      "pathway",
      "bubble",
      "diff",
      "heatmap",
      "ranknet",
      "comparison",
      "lr_contribution",
      "individual_lr",
      "chord",
      "role_scatter",
      "role_heatmap",
      "role_change"
    ),
    condition = NULL,
    measure = c("count", "weight"),
    signaling = NULL,
    pairLR.use = NULL,
    sources.use = NULL,
    targets.use = NULL,
    comparison = c(1, 2),
    dataset = NULL,
    pattern = c("outgoing", "incoming", "all"),
    idents.use = NULL,
    signaling.exclude = NULL,
    layout = "circle",
    title.name = NULL,
    return.data = FALSE,
    file = NULL,
    width = 8,
    height = 6,
    dpi = 300,
    ...
) {
  plot_type <- match.arg(plot_type)
  measure <- match.arg(measure)
  pattern <- match.arg(pattern)
  
  ## ---------- single result modes ----------
  if (plot_type %in% c("overview", "pathway", "lr_contribution", "individual_lr", "chord")) {
    if (is.null(condition)) condition <- "ALL"
    obj <- GetCellChatResult(srt, condition = condition)$cellchat_object
    groupSize <- as.numeric(table(obj@idents))
    
    if (plot_type == "overview") {
      mat <- if (measure == "count") obj@net$count else obj@net$weight
      p <- CellChat::netVisual_circle(
        mat,
        vertex.weight = groupSize,
        weight.scale = TRUE,
        title.name = if (is.null(title.name)) paste(condition, measure) else title.name,
        ...
      )
      .save_plot_if_needed(p, file, width, height, dpi)
      return(p)
    }
    
    if (plot_type == "pathway") {
      if (is.null(signaling)) .stop2("signaling must be provided for plot_type='pathway'")
      p <- CellChat::netVisual_aggregate(
        obj,
        signaling = signaling,
        layout = layout,
        ...
      )
      .save_plot_if_needed(p, file, width, height, dpi)
      return(p)
    }
    
    if (plot_type == "lr_contribution") {
      if (is.null(signaling)) .stop2("signaling must be provided for plot_type='lr_contribution'")
      p <- CellChat::netAnalysis_contribution(
        obj,
        signaling = signaling,
        title = if (is.null(title.name)) signaling else title.name,
        ...
      )
      .save_plot_if_needed(p, file, width, height, dpi)
      return(p)
    }
    
    if (plot_type == "individual_lr") {
      if (is.null(signaling)) .stop2("signaling must be provided for plot_type='individual_lr'")
      if (is.null(pairLR.use)) .stop2("pairLR.use must be provided for plot_type='individual_lr'")
      p <- CellChat::netVisual_individual(
        obj,
        signaling = signaling,
        pairLR.use = pairLR.use,
        layout = layout,
        ...
      )
      .save_plot_if_needed(p, file, width, height, dpi)
      return(p)
    }
    
    if (plot_type == "chord") {
      if (is.null(signaling)) .stop2("signaling must be provided for plot_type='chord'")
      p <- CellChat::netVisual_aggregate(
        obj,
        signaling = signaling,
        layout = "chord",
        ...
      )
      .save_plot_if_needed(p, file, width, height, dpi)
      return(p)
    }
  }
  
  ## ---------- comparison / merged modes ----------
  cmp <- GetCellChatComparison(srt, condition = condition)
  merged_obj <- cmp$comparison_object
  obj.list <- cmp$object.list
  
  comp_idx <- .resolve_comparison_dataset(cmp, comparison = comparison)
  src_idx <- .resolve_group_index_merged(merged_obj, sources.use)
  tgt_idx <- .resolve_group_index_merged(merged_obj, targets.use)
  
  if (plot_type == "diff") {
    p <- CellChat::netVisual_diffInteraction(
      merged_obj,
      measure = if (measure == "count") "count" else "weight",
      ...
    )
    .save_plot_if_needed(p, file, width, height, dpi)
    return(p)
  }
  
  if (plot_type == "heatmap") {
    if (is.null(signaling)) {
      p <- CellChat::netVisual_heatmap(
        merged_obj,
        measure = if (measure == "count") "count" else "weight",
        ...
      )
      .save_plot_if_needed(p, file, width, height, dpi)
      return(p)
    } else {
      ht <- lapply(seq_along(obj.list), function(i) {
        CellChat::netVisual_heatmap(
          obj.list[[i]],
          signaling = signaling,
          title.name = names(obj.list)[i],
          ...
        )
      })
      return(ht)
    }
  }
  
  if (plot_type == "ranknet") {
    p <- CellChat::rankNet(
      merged_obj,
      mode = "comparison",
      stacked = TRUE,
      do.stat = TRUE,
      ...
    )
    .save_plot_if_needed(p, file, width, height, dpi)
    return(p)
  }
  
  if (plot_type %in% c("comparison", "bubble")) {
    ## comparison/bubble merged mode
    p <- CellChat::netVisual_bubble(
      merged_obj,
      sources.use = src_idx,
      targets.use = tgt_idx,
      signaling = signaling,
      pairLR.use = pairLR.use,
      comparison = comp_idx,
      return.data = return.data,
      title.name = title.name,
      ...
    )
    if (!isTRUE(return.data)) .save_plot_if_needed(p, file, width, height, dpi)
    return(p)
  }
  
  if (plot_type == "role_scatter") {
    if (is.null(dataset)) dataset <- 1
    if (is.character(dataset)) dataset <- match(dataset, names(obj.list))
    if (is.na(dataset) || dataset < 1 || dataset > length(obj.list)) {
      .stop2("dataset must be one of: %s", paste(names(obj.list), collapse = ", "))
    }
    p <- CellChat::netAnalysis_signalingRole_scatter(
      obj.list[[dataset]],
      title = if (is.null(title.name)) names(obj.list)[dataset] else title.name,
      ...
    )
    .save_plot_if_needed(p, file, width, height, dpi)
    return(p)
  }
  
  if (plot_type == "role_heatmap") {
    if (is.null(dataset)) dataset <- 1
    if (is.character(dataset)) dataset <- match(dataset, names(obj.list))
    if (is.na(dataset) || dataset < 1 || dataset > length(obj.list)) {
      .stop2("dataset must be one of: %s", paste(names(obj.list), collapse = ", "))
    }
    p <- CellChat::netAnalysis_signalingRole_heatmap(
      obj.list[[dataset]],
      pattern = pattern,
      signaling = signaling,
      title = if (is.null(title.name)) names(obj.list)[dataset] else title.name,
      ...
    )
    return(p)
  }
  
  if (plot_type == "role_change") {
    if (is.null(idents.use)) .stop2("idents.use must be provided for plot_type='role_change'")
    p <- CellChat::netAnalysis_signalingChanges_scatter(
      merged_obj,
      idents.use = idents.use,
      signaling.exclude = signaling.exclude,
      ...
    )
    .save_plot_if_needed(p, file, width, height, dpi)
    return(p)
  }
}
## =========================
## PATCH: robust bubble / chord / export / table
## append to the END of your current cellchat_module_fixed.R
## =========================

.check_pkg2 <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg), call. = FALSE)
  }
}

.cc_get_single_obj <- function(srt, condition = "ALL") {
  GetCellChatResult(srt, condition = condition)$cellchat_object
}

.cc_get_cmp <- function(srt, condition) {
  GetCellChatComparison(srt, condition = condition)
}

.cc_resolve_dataset_index <- function(cmp, comparison = c(1, 2)) {
  ds_names <- names(cmp$object.list)
  if (is.numeric(comparison)) return(comparison)
  idx <- match(as.character(comparison), ds_names)
  if (any(is.na(idx))) {
    stop(
      sprintf(
        "comparison names not found: %s. Available datasets: %s",
        paste(as.character(comparison)[is.na(idx)], collapse = ", "),
        paste(ds_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  idx
}

.cc_resolve_group_index_single <- function(object, groups = NULL) {
  if (is.null(groups)) return(NULL)
  lev <- levels(object@idents)
  if (is.numeric(groups)) return(groups)
  idx <- match(as.character(groups), lev)
  if (any(is.na(idx))) {
    stop(
      sprintf(
        "These cell groups were not found: %s. Available groups: %s",
        paste(as.character(groups)[is.na(idx)], collapse = ", "),
        paste(lev, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  idx
}

.cc_resolve_group_index_merged <- function(object, groups = NULL) {
  if (is.null(groups)) return(NULL)
  lev <- if (is.list(object@idents) && !is.null(object@idents$joint)) {
    levels(object@idents$joint)
  } else {
    levels(object@idents)
  }
  if (is.numeric(groups)) return(groups)
  idx <- match(as.character(groups), lev)
  if (any(is.na(idx))) {
    stop(
      sprintf(
        "These cell groups were not found in merged object: %s. Available groups: %s",
        paste(as.character(groups)[is.na(idx)], collapse = ", "),
        paste(lev, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  idx
}

.cc_safe_ggsave <- function(file, plot, width = 8, height = 6, dpi = 300) {
  if (is.null(file)) return(invisible(NULL))
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = file, plot = plot, width = width, height = height, dpi = dpi)
  invisible(file)
}

.cc_pick_dataset_name <- function(cmp, dataset = 1) {
  nm <- names(cmp$object.list)
  if (is.character(dataset)) {
    if (!dataset %in% nm) {
      stop(sprintf("dataset '%s' not found. Available: %s", dataset, paste(nm, collapse = ", ")), call. = FALSE)
    }
    return(dataset)
  }
  if (dataset < 1 || dataset > length(nm)) {
    stop(sprintf("dataset index out of range. Available: 1..%s", length(nm)), call. = FALSE)
  }
  nm[dataset]
}

## -------------------------
## unified table extractor
## -------------------------
CellChatTable <- function(
    srt,
    condition,
    mode = c("single", "comparison"),
    dataset = 1,
    signaling = NULL,
    sources.use = NULL,
    targets.use = NULL,
    pairLR.use = NULL,
    comparison = c(1, 2),
    slot.name = "net",
    thresh = 0.05,
    save.file = NULL
) {
  .check_pkg2("CellChat")
  mode <- match.arg(mode)
  
  if (mode == "single") {
    obj <- .cc_get_single_obj(srt, condition)
    src_idx <- .cc_resolve_group_index_single(obj, sources.use)
    tgt_idx <- .cc_resolve_group_index_single(obj, targets.use)
    
    df <- CellChat::subsetCommunication(
      object = obj,
      slot.name = slot.name,
      sources.use = src_idx,
      targets.use = tgt_idx,
      signaling = signaling,
      pairLR.use = pairLR.use,
      thresh = thresh
    )
    
    if (!is.null(save.file)) {
      dir.create(dirname(save.file), recursive = TRUE, showWarnings = FALSE)
      utils::write.csv(df, save.file, row.names = FALSE)
    }
    return(df)
  }
  
  cmp <- .cc_get_cmp(srt, condition)
  dataset_name <- .cc_pick_dataset_name(cmp, dataset)
  obj <- cmp$object.list[[dataset_name]]
  
  src_idx <- .cc_resolve_group_index_single(obj, sources.use)
  tgt_idx <- .cc_resolve_group_index_single(obj, targets.use)
  
  df <- CellChat::subsetCommunication(
    object = obj,
    slot.name = slot.name,
    sources.use = src_idx,
    targets.use = tgt_idx,
    signaling = signaling,
    pairLR.use = pairLR.use,
    thresh = thresh
  )
  
  df <- as.data.frame(df)
  df$dataset <- rep(dataset_name, nrow(df))
  
  ## try to add delta info from merged object if possible
  if (length(comparison) == 2) {
    df$comparison_name <- paste0(names(cmp$object.list)[comparison[1]], "_vs_", names(cmp$object.list)[comparison[2]])
  }
  
  if (!is.null(save.file)) {
    dir.create(dirname(save.file), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(df, save.file, row.names = FALSE)
  }
  
  df
}

## -------------------------
## robust custom bubble plot from table
## -------------------------
.CellChatBubblePlotFromTable <- function(
    df,
    title.name = NULL,
    color.by = c("prob", "pval"),
    remove.isolate = TRUE
) {
  .check_pkg2("ggplot2")
  color.by <- match.arg(color.by)
  
  if (is.null(df) || nrow(df) == 0) {
    stop("No communication records available for bubble plot.", call. = FALSE)
  }
  
  df <- as.data.frame(df)
  req_cols <- c("source", "target", "interaction_name", "prob", "pval")
  miss <- setdiff(req_cols, colnames(df))
  if (length(miss) > 0) {
    stop(sprintf("Missing required columns for bubble plot: %s", paste(miss, collapse = ", ")), call. = FALSE)
  }
  
  df$pair <- paste(df$source, "->", df$target)
  df$interaction_plot <- df$interaction_name
  
  if ("dataset" %in% colnames(df)) {
    df$pair <- paste0(df$pair, " [", df$dataset, "]")
  }
  
  if (isTRUE(remove.isolate)) {
    df <- df[df$prob > 0, , drop = FALSE]
  }
  
  if (nrow(df) == 0) {
    stop("No non-zero communication records remain after filtering.", call. = FALSE)
  }
  
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = pair,
      y = interaction_plot,
      size = prob,
      color = if (color.by == "pval") -log10(pval + 1e-300) else prob
    )
  ) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = title.name %||% "CellChat bubble plot",
      size = "prob",
      color = if (color.by == "pval") "-log10(pval)" else "prob"
    )
  
  p
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

## -------------------------
## optimized chord plot
## -------------------------
CellChatChord <- function(
    srt,
    condition = "ALL",
    signaling,
    dataset = NULL,
    mode = c("single", "comparison"),
    sources.use = NULL,
    targets.use = NULL,
    max.groups = 8,
    top.n = 30,
    reduce = TRUE,
    small.gap = 1,
    big.gap = 8,
    lab.cex = 0.5,
    title.name = NULL,
    file = NULL,
    width = 10,
    height = 10
) {
  .check_pkg2("CellChat")
  .check_pkg2("circlize")
  mode <- match.arg(mode)
  
  if (mode == "single") {
    obj <- .cc_get_single_obj(srt, condition)
  } else {
    cmp <- .cc_get_cmp(srt, condition)
    dataset_name <- .cc_pick_dataset_name(cmp, dataset %||% 1)
    obj <- cmp$object.list[[dataset_name]]
  }
  
  all_groups <- levels(obj@idents)
  
  src_idx <- .cc_resolve_group_index_single(obj, sources.use)
  tgt_idx <- .cc_resolve_group_index_single(obj, targets.use)
  
  ## use table first to reduce complexity
  df <- tryCatch(
    CellChat::subsetCommunication(
      object = obj,
      signaling = signaling,
      sources.use = src_idx,
      targets.use = tgt_idx,
      thresh = 0.05
    ),
    error = function(e) NULL
  )
  
  if (is.null(df) || nrow(df) == 0) {
    stop(sprintf("There is no significant communication of %s", paste(signaling, collapse = ", ")), call. = FALSE)
  }
  
  keep_groups <- unique(c(df$source, df$target))
  
  if (isTRUE(reduce)) {
    df <- df[order(df$prob, decreasing = TRUE), , drop = FALSE]
    if (nrow(df) > top.n) df <- df[seq_len(top.n), , drop = FALSE]
    keep_groups <- unique(c(df$source, df$target))
    
    if (length(keep_groups) > max.groups) {
      score <- sort(tapply(df$prob, c(df$source, df$target), sum, na.rm = TRUE), decreasing = TRUE)
      keep_groups <- names(score)[seq_len(min(max.groups, length(score)))]
      df <- df[df$source %in% keep_groups & df$target %in% keep_groups, , drop = FALSE]
    }
  }
  
  src_use2 <- match(unique(df$source), all_groups)
  tgt_use2 <- match(unique(df$target), all_groups)
  src_use2 <- src_use2[!is.na(src_use2)]
  tgt_use2 <- tgt_use2[!is.na(tgt_use2)]
  
  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    grDevices::pdf(file, width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  
  circlize::circos.clear()
  on.exit(try(circlize::circos.clear(), silent = TRUE), add = TRUE)
  
  CellChat::netVisual_chord_cell(
    object = obj,
    signaling = signaling,
    sources.use = src_use2,
    targets.use = tgt_use2,
    small.gap = small.gap,
    big.gap = big.gap,
    lab.cex = lab.cex,
    title.name = title.name %||% paste(signaling, collapse = ", ")
  )
}

## -------------------------
## overwrite CellChatPlot with robust bubble/chord behavior
## -------------------------
CellChatPlot <- function(
    srt,
    plot_type = c(
      "overview",
      "pathway",
      "bubble",
      "diff",
      "heatmap",
      "ranknet",
      "comparison",
      "lr_contribution",
      "individual_lr",
      "chord",
      "role_scatter",
      "role_heatmap",
      "role_change"
    ),
    condition = NULL,
    measure = c("count", "weight"),
    signaling = NULL,
    pairLR.use = NULL,
    sources.use = NULL,
    targets.use = NULL,
    comparison = c(1, 2),
    dataset = 1,
    pattern = c("outgoing", "incoming", "all"),
    idents.use = NULL,
    signaling.exclude = NULL,
    layout = "circle",
    title.name = NULL,
    return.data = FALSE,
    file = NULL,
    width = 8,
    height = 6,
    dpi = 300,
    ...
) {
  .check_pkg2("CellChat")
  .check_pkg2("ggplot2")
  
  plot_type <- match.arg(plot_type)
  measure <- match.arg(measure)
  pattern <- match.arg(pattern)
  
  ## single-result family
  if (plot_type %in% c("overview", "pathway", "lr_contribution", "individual_lr", "chord")) {
    if (is.null(condition)) condition <- "ALL"
    obj <- .cc_get_single_obj(srt, condition)
    groupSize <- as.numeric(table(obj@idents))
    
    if (plot_type == "overview") {
      mat <- if (measure == "count") obj@net$count else obj@net$weight
      p <- CellChat::netVisual_circle(
        mat,
        vertex.weight = groupSize,
        weight.scale = TRUE,
        title.name = title.name %||% paste(condition, measure),
        ...
      )
      if (!is.null(file) && inherits(p, c("gg", "ggplot"))) .cc_safe_ggsave(file, p, width, height, dpi)
      return(p)
    }
    
    if (plot_type == "pathway") {
      if (is.null(signaling)) stop("signaling must be provided for plot_type='pathway'", call. = FALSE)
      p <- CellChat::netVisual_aggregate(
        obj,
        signaling = signaling,
        layout = layout,
        ...
      )
      if (!is.null(file) && inherits(p, c("gg", "ggplot"))) .cc_safe_ggsave(file, p, width, height, dpi)
      return(p)
    }
    
    if (plot_type == "lr_contribution") {
      if (is.null(signaling)) stop("signaling must be provided for plot_type='lr_contribution'", call. = FALSE)
      p <- CellChat::netAnalysis_contribution(
        obj,
        signaling = signaling,
        title = title.name %||% signaling,
        ...
      )
      if (!is.null(file)) .cc_safe_ggsave(file, p, width, height, dpi)
      return(p)
    }
    
    if (plot_type == "individual_lr") {
      if (is.null(signaling)) stop("signaling must be provided for plot_type='individual_lr'", call. = FALSE)
      if (is.null(pairLR.use)) stop("pairLR.use must be provided for plot_type='individual_lr'", call. = FALSE)
      p <- CellChat::netVisual_individual(
        obj,
        signaling = signaling,
        pairLR.use = pairLR.use,
        layout = layout,
        ...
      )
      return(p)
    }
    
    if (plot_type == "chord") {
      if (is.null(signaling)) stop("signaling must be provided for plot_type='chord'", call. = FALSE)
      return(
        CellChatChord(
          srt = srt,
          condition = condition,
          signaling = signaling,
          mode = "single",
          sources.use = sources.use,
          targets.use = targets.use,
          file = file,
          width = width,
          height = height
        )
      )
    }
  }
  
  ## special: single-condition bubble
  if (plot_type == "bubble" && !is.null(condition) && condition %in% ListCellChatResults(srt)$results) {
    df <- CellChatTable(
      srt = srt,
      condition = condition,
      mode = "single",
      signaling = signaling,
      sources.use = sources.use,
      targets.use = targets.use,
      pairLR.use = pairLR.use
    )
    if (isTRUE(return.data)) return(df)
    
    p <- .CellChatBubblePlotFromTable(
      df,
      title.name = title.name %||% paste0(condition, " bubble")
    )
    if (!is.null(file)) .cc_safe_ggsave(file, p, width, height, dpi)
    return(p)
  }
  
  ## comparison family
  cmp <- .cc_get_cmp(srt, condition)
  merged_obj <- cmp$comparison_object
  obj.list <- cmp$object.list
  
  if (plot_type == "diff") {
    p <- CellChat::netVisual_diffInteraction(
      merged_obj,
      measure = if (measure == "count") "count" else "weight",
      ...
    )
    if (!is.null(file) && inherits(p, c("gg", "ggplot"))) .cc_safe_ggsave(file, p, width, height, dpi)
    return(p)
  }
  
  if (plot_type == "heatmap") {
    if (is.null(signaling)) {
      p <- CellChat::netVisual_heatmap(
        merged_obj,
        measure = if (measure == "count") "count" else "weight",
        ...
      )
      return(p)
    } else {
      ht <- lapply(seq_along(obj.list), function(i) {
        CellChat::netVisual_heatmap(
          obj.list[[i]],
          signaling = signaling,
          title.name = names(obj.list)[i],
          ...
        )
      })
      return(ht)
    }
  }
  
  if (plot_type == "ranknet") {
    p <- CellChat::rankNet(
      merged_obj,
      mode = "comparison",
      stacked = TRUE,
      do.stat = TRUE,
      ...
    )
    if (!is.null(file) && inherits(p, c("gg", "ggplot"))) .cc_safe_ggsave(file, p, width, height, dpi)
    return(p)
  }
  
  ## robust comparison bubble/comparison:
  ## no longer depend on netVisual_bubble when user filters source/target
  if (plot_type %in% c("comparison", "bubble")) {
    if (isTRUE(return.data)) {
      out <- list()
      for (ds in names(obj.list)) {
        out[[ds]] <- CellChatTable(
          srt = srt,
          condition = condition,
          mode = "comparison",
          dataset = ds,
          signaling = signaling,
          sources.use = sources.use,
          targets.use = targets.use,
          pairLR.use = pairLR.use,
          comparison = comparison
        )
      }
      return(out)
    }
    
    ## if no filtering at all, still allow original function
    if (is.null(sources.use) && is.null(targets.use) && is.null(pairLR.use) && is.null(signaling)) {
      p <- CellChat::netVisual_bubble(
        merged_obj,
        comparison = comparison,
        title.name = title.name,
        ...
      )
      return(p)
    }
    
    ## filtered comparison -> custom robust ggplot
    df_list <- lapply(names(obj.list), function(ds) {
      CellChatTable(
        srt = srt,
        condition = condition,
        mode = "comparison",
        dataset = ds,
        signaling = signaling,
        sources.use = sources.use,
        targets.use = targets.use,
        pairLR.use = pairLR.use,
        comparison = comparison
      )
    })
    
    df_list <- df_list[vapply(df_list, nrow, integer(1)) > 0]
    
    if (length(df_list) == 0) {
      stop("No communication records found under the current filtering conditions.", call. = FALSE)
    }
    
    df <- do.call(rbind, df_list)
    
    p <- .CellChatBubblePlotFromTable(
      df,
      title.name = title.name %||% paste0(condition, " comparison bubble")
    )
    if (!is.null(file)) .cc_safe_ggsave(file, p, width, height, dpi)
    return(p)
  }
  
  if (plot_type == "role_scatter") {
    ds_name <- .cc_pick_dataset_name(cmp, dataset)
    p <- CellChat::netAnalysis_signalingRole_scatter(
      obj.list[[ds_name]],
      title = title.name %||% ds_name,
      ...
    )
    if (!is.null(file)) .cc_safe_ggsave(file, p, width, height, dpi)
    return(p)
  }
  
  if (plot_type == "role_heatmap") {
    ds_name <- .cc_pick_dataset_name(cmp, dataset)
    p <- CellChat::netAnalysis_signalingRole_heatmap(
      obj.list[[ds_name]],
      pattern = pattern,
      signaling = signaling,
      title = title.name %||% ds_name,
      ...
    )
    return(p)
  }
  
  if (plot_type == "role_change") {
    if (is.null(idents.use)) stop("idents.use must be provided for plot_type='role_change'", call. = FALSE)
    p <- CellChat::netAnalysis_signalingChanges_scatter(
      merged_obj,
      idents.use = idents.use,
      signaling.exclude = signaling.exclude,
      ...
    )
    if (!is.null(file)) .cc_safe_ggsave(file, p, width, height, dpi)
    return(p)
  }
  
  stop("Unsupported plot_type.", call. = FALSE)
}