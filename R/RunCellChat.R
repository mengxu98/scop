#' @title Run CellChat
#'
#' @description
#' [RunCellChat] performs CellChat analysis on a Seurat object to investigate cell-to-cell communication.
#' The results are stored in the Seurat object and can be visualized using [CellChatPlot].
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams standard_scop
#' @inheritParams CellDimPlot
#' @param species The species of the data, either 'human', 'mouse' or 'zebrafish'.
#' @param annotation_selected A vector of cell annotations of interest for running the CellChat analysis.
#' If not provided, all cell types will be considered.
#' @param group_column Name of the metadata column in the Seurat object that defines conditions or groups.
#' @param group_cmp A list of pairwise condition comparisons for differential CellChat analysis.
#' @param thresh The threshold for computing centrality scores.
#' @param min.cells the minmum number of expressed cells required for the genes that are considered for cell-cell communication analysis.
#'
#' @export
#'
#' @seealso
#' [CellChatPlot]
#'
#' @references
#' [CellChat](https://github.com/jinworks/CellChat),
#' [scDown::run_cellchatV2](https://htmlpreview.github.io/?https://raw.githubusercontent.com/BCH-RC/scDown/main/vignettes/scDown_CellChatV2.html)
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunCellChat(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   species = "Mus_musculus"
#' )
#'
#' CellChatPlot(pancreas_sub)
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
    verbose = TRUE) {
  log_message(
    "Start {.pkg CellChat} analysis",
    verbose = verbose
  )

  check_r("jinworks/CellChat", verbose = FALSE)
  check_r("immunogenomics/presto", verbose = FALSE)
  if (!group.by %in% colnames(srt@meta.data)) {
    log_message(
      "{.val {group.by}} does not exist in {.cls Seurat}",
      message_type = "error"
    )
  }

  if (!is.null(split.by) && !split.by %in% colnames(srt@meta.data)) {
    log_message(
      "{.val {split.by}} does not exist in {.cls Seurat}",
      message_type = "error"
    )
  }

  if (!is.null(group_column) && !group_column %in% colnames(srt@meta.data)) {
    log_message(
      "{.val {group_column}} does not exist in {.cls Seurat}",
      message_type = "error"
    )
  }

  if (!is.null(annotation_selected)) {
    available_annotations <- unique(srt@meta.data[, group.by])
    missing_annotations <- setdiff(annotation_selected, available_annotations)
    if (length(missing_annotations) > 0) {
      log_message(
        "Missing annotations: {.val {missing_annotations}}",
        message_type = "error"
      )
    }
  }

  if (!is.null(group_cmp)) {
    available_groups <- unique(srt@meta.data[, group_column])
    all_groups_in_cmp <- unique(unlist(group_cmp))
    missing_groups <- setdiff(all_groups_in_cmp, available_groups)
    if (length(missing_groups) > 0) {
      log_message(
        "Missing groups in comparisons: {.val {missing_groups}}",
        message_type = "error"
      )
    }
  }
  species <- match.arg(species)
  if (species == "Mus_musculus") {
    species <- "mouse"
  }
  if (species == "Homo_sapiens") {
    species <- "human"
  }
  Idents(srt) <- srt@meta.data[, group.by]

  if (!is.null(split.by)) {
    srt$samples <- as.factor(srt@meta.data[, split.by])
  } else {
    srt$samples <- as.factor("All")
  }

  cellchat_results <- list()
  if (is.null(group_column)) {
    cellchat_object <- .DoCellChat(object = srt, species = species)

    if (is.null(annotation_selected)) {
      srt_condition <- srt
    } else {
      cellchat_object <- .SubsetCellChatMod(
        cellchat_object,
        idents.use = annotation_selected,
        thresh = thresh
      )
      cellchat_object <- CellChat::netAnalysis_computeCentrality(
        cellchat_object,
        thresh = thresh
      )
      srt_condition <- srt[, srt@meta.data[, group.by] %in% annotation_selected]
    }

    cellchat_results[["ALL"]] <- list(
      cellchat_object = cellchat_object,
      seurat_object = srt_condition
    )
  } else {
    conditions <- unique(as.vector(srt@meta.data[, group_column]))
    metadata_cond <- FetchData(object = srt, vars = group_column)

    for (condition in conditions) {
      log_message(
        "Processing condition: {.val {condition}}",
        verbose = verbose
      )

      srt_condition <- srt[, which(x = (metadata_cond == condition))]
      cellchat_object_condition <- .DoCellChat(
        object = srt_condition, species = species
      )

      if (!is.null(annotation_selected)) {
        cellchat_object_condition <- .SubsetCellChatMod(
          cellchat_object_condition,
          idents.use = annotation_selected,
          thresh = thresh
        )
        cellchat_object_condition <- CellChat::netAnalysis_computeCentrality(
          cellchat_object_condition,
          thresh = thresh
        )
        srt_condition <- srt_condition[, srt_condition@meta.data[, group.by] %in% annotation_selected]
      }

      cellchat_results[[condition]] <- list(
        cellchat_object = cellchat_object_condition,
        seurat_object = srt_condition
      )
    }
  }

  srt@tools[["CellChat"]][["results"]] <- cellchat_results
  srt@tools[["CellChat"]][["parameters"]] <- list(
    group.by = group.by,
    species = species,
    split.by = split.by,
    annotation_selected = annotation_selected,
    group_column = group_column,
    group_cmp = group_cmp
  )

  log_message(
    "CellChat analysis completed",
    message_type = "success",
    verbose = verbose
  )

  return(srt)
}

.DoCellChat <- function(
    object,
    assay = "RNA",
    species = c("human", "mouse", "zebrafish"),
    thresh = 0.05,
    min.cells = 10) {
  metadata <- data.frame(label = Idents(object))
  metadata <- cbind(metadata, object@meta.data)
  object <- CellChat::createCellChat(
    GetAssayData5(object, layer = "data", assay = assay),
    meta = metadata, group.by = "label"
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

  return(object)
}

.SubsetCellChatMod <- function(
    object,
    idents.use,
    thresh = 0.05,
    verbose = TRUE) {
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

  data_subset <- object@data[, level_use_index]
  data_signaling_subset <- object@data.signaling[, level_use_index]
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
          values_new <- values[level_use_index, ]
          images[[images_j]] <- values_new
        }
        if (images_j %in% c("distance")) {
          values_new <- values[group_existing_index, group_existing_index, drop = FALSE]
          images[[images_j]] <- values_new
        }
      }
      images_subset[[i]] <- images

      net <- object@net[[i]]
      for (net.j in names(net)) {
        values <- net[[net.j]]
        if (net.j %in% c("prob", "pval")) {
          values_new <- values[group_existing_index, group_existing_index, ]
          net[[net.j]] <- values_new
        }
        if (net.j %in% c("count", "sum", "weight")) {
          values_new <- values[group_existing_index, group_existing_index]
          net[[net.j]] <- values_new
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
    log_message(
      "Update slots object@images, object@net, object@netP in a single dataset...",
      verbose = verbose
    )

    group_existing_index <- which(level_use0 %in% level_use)

    images <- object@images
    for (images_j in names(images)) {
      values <- images[[images_j]]
      if (images_j %in% c("coordinates")) {
        values_new <- values[level_use_index, ]
        images[[images_j]] <- values_new
      }
      if (images_j %in% c("distance")) {
        values_new <- values[group_existing_index, group_existing_index, drop = FALSE]
        images[[images_j]] <- values_new
      }
    }
    images_subset <- images

    net <- object@net
    for (net.j in names(net)) {
      values <- net[[net.j]]
      if (net.j %in% c("prob", "pval")) {
        values_new <- values[group_existing_index, group_existing_index, , drop = FALSE]
        net[[net.j]] <- values_new
      }
      if (net.j %in% c("count", "sum", "weight")) {
        values_new <- values[group_existing_index, group_existing_index, drop = FALSE]
        net[[net.j]] <- values_new
      }
    }
    net_subset <- net

    netP <- CellChat::computeCommunProbPathway(
      net = net_subset,
      pairLR.use = object@LR$LRsig,
      thresh = thresh
    )
    netP$centr <- CellChat::netAnalysis_computeCentrality(
      net = net_subset$prob
    )
    netP_subset <- netP
    idents_subset <- object@idents[level_use_index]
    idents_subset <- factor(idents_subset, levels = level_use)
  }

  object_subset <- methods::new(
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
  return(object_subset)
}
