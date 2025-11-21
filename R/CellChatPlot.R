#' @title Plot CellChat analysis results
#'
#' @description
#' [CellChatPlot] creates various visualizations for CellChat analysis results stored in a Seurat object.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt A Seurat object that has been processed with [RunCellChat].
#' @param plot_type Type of plot to create.
#' Options: `"aggregate"`, `"pathway"`, `"comparison"`,
#' `"heatmap"`, `"circle"`, `"bubble"`, `"gene"`.
#' @param condition Condition to plot (if multiple conditions exist).
#' @param pathway Specific pathway to visualize (for pathway, bubble, and gene plots).
#' If `NULL`, uses top pathways.
#' @param dirpath Directory to save plots.
#' @param output_format Format of output figure: `"png"` or `"pdf"`.
#' Default is `"png"`.
#' @param top_n Number of top pathways to use for plotting.
#' Default is `10`.
#' @param base_height Base height multiplier for all plots.
#' Default is `1`.
#' @param base_width Base width multiplier for all plots.
#' Default is `1`.
#' @param res Resolution for PNG output.
#' Default is `300`.
#'
#' @export
#'
#' @seealso
#' [RunCellChat]
#'
#' @examples
#' options(log_message.verbose = FALSE)
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunCellChat(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   species = "mouse"
#' )
#'
#' CellChatPlot(pancreas_sub, plot_type = "aggregate")
#'
#' CellChatPlot(pancreas_sub, plot_type = "pathway")
#'
#' CellChatPlot(pancreas_sub, plot_type = "bubble")
#'
#' CellChatPlot(pancreas_sub, plot_type = "gene")
#'
#' CellChatPlot(pancreas_sub, plot_type = "heatmap")
CellChatPlot <- function(
    srt,
    plot_type = "aggregate",
    condition = NULL,
    pathway = NULL,
    dirpath = NULL,
    output_format = "pdf",
    top_n = 10,
    base_height = 1,
    base_width = 1,
    res = 300,
    verbose = TRUE) {
  if (is.null(srt@tools[["CellChat"]])) {
    log_message(
      "No CellChat results found in the Seurat object. Please run {.fn RunCellChat} first.",
      message_type = "error"
    )
  }

  cellchat_results <- srt@tools[["CellChat"]][["results"]]
  parameters <- srt@tools[["CellChat"]][["parameters"]]

  if (is.null(condition)) {
    if (length(cellchat_results) == 1) {
      condition <- names(cellchat_results)[1]
    } else {
      log_message(
        "Multiple conditions found. Please specify condition parameter",
        message_type = "error"
      )
    }
  }

  if (!condition %in% names(cellchat_results)) {
    log_message(
      "Condition {.val {condition}} not found in CellChat results",
      message_type = "error"
    )
  }

  cellchat_object <- cellchat_results[[condition]]$cellchat_object
  seurat_object <- cellchat_results[[condition]]$seurat_object

  if (!is.null(dirpath)) {
    if (!dir.exists(dirpath)) {
      dir.create(dirpath, recursive = TRUE)
    }
  }

  log_message(
    "Creating {.val {plot_type}} plot for condition {.val {condition}}",
    verbose = verbose
  )

  if (!is.null(dirpath)) {
    if (!dir.exists(dirpath)) {
      dir.create(dirpath, recursive = TRUE)
    }
  }

  if (plot_type == "aggregate") {
    p <- .create_aggregate_plots(
      cellchat_object = cellchat_object,
      condition = condition,
      dirpath = dirpath,
      output_format = output_format,
      base_height = base_height,
      base_width = base_width,
      res = res
    )
  } else if (plot_type == "pathway") {
    if (is.null(pathway)) {
      pathways_top <- cellchat_object@netP$pathways[1:top_n]
      pathways_top <- pathways_top[!is.na(pathways_top)]
      log_message(
        "Using top {.val {length(pathways_top)}} pathways",
        verbose = verbose
      )
    } else {
      pathways_top <- pathway
    }

    p <- .create_pathway_plots(
      cellchat_object = cellchat_object,
      seurat_object = seurat_object,
      pathways_to_show = pathways_top,
      condition = condition,
      dirpath = dirpath,
      output_format = output_format,
      base_height = base_height,
      base_width = base_width,
      res = res
    )
  } else if (plot_type == "heatmap") {
    p <- .create_heatmap_plots(
      cellchat_object = cellchat_object,
      condition = condition,
      dirpath = dirpath,
      output_format = output_format,
      base_height = base_height,
      base_width = base_width,
      res = res
    )
  } else if (plot_type == "circle") {
    p <- .create_circle_plots(
      cellchat_object = cellchat_object,
      condition = condition,
      dirpath = dirpath,
      output_format = output_format,
      base_height = base_height,
      base_width = base_width,
      res = res
    )
  } else if (plot_type == "bubble") {
    if (is.null(pathway)) {
      pathways_top <- cellchat_object@netP$pathways[1:min(3, top_n)]
      pathways_top <- pathways_top[!is.na(pathways_top)]
      log_message(
        "No pathway specified. Using top {.val {length(pathways_top)}} pathways for bubble plots",
        verbose = verbose
      )
      p2_list <- lapply(
        pathways_top, function(path) {
          .create_bubble_plots(
            cellchat_object = cellchat_object,
            pathway = path,
            condition = condition,
            dirpath = dirpath,
            output_format = output_format,
            base_height = base_height,
            base_width = base_width,
            res = res
          )
        }
      )
      p <- patchwork::wrap_plots(p2_list)
    } else {
      p <- .create_bubble_plots(
        cellchat_object = cellchat_object,
        pathway = pathway,
        condition = condition,
        dirpath = dirpath,
        output_format = output_format,
        base_height = base_height,
        base_width = base_width,
        res = res
      )
    }
  } else if (plot_type == "gene") {
    if (is.null(pathway)) {
      pathways_top <- cellchat_object@netP$pathways[1:min(3, top_n)]
      pathways_top <- pathways_top[!is.na(pathways_top)]
      log_message(
        "No pathway specified. Using top {.val {length(pathways_top)}} pathways for gene expression plots",
        verbose = verbose
      )

      p2_list <- lapply(
        pathways_top, function(path) {
          .create_gene_plots(
            cellchat_object = cellchat_object,
            seurat_object = seurat_object,
            pathway = path,
            condition = condition,
            dirpath = dirpath,
            output_format = output_format,
            species = parameters$species,
            base_height = base_height,
            base_width = base_width,
            res = res
          )
        }
      )
      p <- patchwork::wrap_plots(p2_list)
    } else {
      p <- .create_gene_plots(
        cellchat_object = cellchat_object,
        seurat_object = seurat_object,
        pathway = pathway,
        condition = condition,
        dirpath = dirpath,
        output_format = output_format,
        species = parameters$species,
        base_height = base_height,
        base_width = base_width,
        res = res
      )
    }
  } else if (plot_type == "comparison") {
    if (is.null(parameters$group_cmp)) {
      log_message(
        "No group comparisons found in CellChat results",
        message_type = "error"
      )
    }
    log_message(
      "Comparison plots not yet implemented",
      message_type = "warning"
    )
  } else {
    log_message(
      "Unknown plot_type: {.val {plot_type}}. Available options: aggregate, pathway, heatmap, circle, bubble, gene, comparison",
      message_type = "error"
    )
  }

  log_message(
    "Plot creation completed",
    message_type = "success",
    verbose = verbose
  )
  return(p)
}

.create_output_device <- function(
    filename,
    output_format = "pdf",
    height = 6,
    width = 6,
    res = 300,
    units = "in") {
  if (output_format == "png") {
    grDevices::png(
      filename,
      height = height,
      width = width,
      res = res,
      units = units
    )
  } else if (output_format == "pdf") {
    grDevices::pdf(
      filename,
      height = height,
      width = width
    )
  }
}

.save_ggplot <- function(
    plot,
    filename,
    output_format = "pdf",
    height = 6,
    width = 6,
    res = 300) {
  if (output_format == "png") {
    ggsave(
      filename,
      plot = plot,
      height = height,
      width = width,
      dpi = res
    )
  } else if (output_format == "pdf") {
    ggsave(filename, plot = plot, height = height, width = width)
  }
}

.create_aggregate_plots <- function(
    cellchat_object,
    condition,
    dirpath = NULL,
    output_format = "pdf",
    base_height = 1,
    base_width = 1,
    res = 600) {
  group_size <- as.numeric(table(cellchat_object@idents))
  numofcelltypes <- length(group_size)
  height_in <- max(4, 2 * (numofcelltypes / 4)) * base_height
  width_in <- max(6, 8 / 3 * (numofcelltypes / 4 + 1)) * base_width

  if (!is.null(dirpath)) {
    circle_filename <- paste0(
      dirpath, "/", condition, "_net_interaction_and_weight.", output_format
    )
    .create_output_device(
      circle_filename,
      output_format = output_format,
      height = height_in,
      width = width_in,
      res = res
    )
  }

  graphics::par(mfrow = c(1, 2), xpd = TRUE)
  CellChat::netVisual_circle(
    cellchat_object@net$count,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = "Number of interactions"
  )
  CellChat::netVisual_circle(
    cellchat_object@net$weight,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = "Interaction weights/strength"
  )

  if (!is.null(dirpath)) {
    grDevices::dev.off()
  }

  p_scatter <- CellChat::netAnalysis_signalingRole_scatter(cellchat_object)
  print(p_scatter)

  if (!is.null(dirpath)) {
    p_scatter_filename <- paste0(
      dirpath, "/", condition, "_signaling_role.", output_format
    )
    .save_ggplot(
      p_scatter,
      p_scatter_filename,
      output_format = output_format,
      height = 6 * base_height,
      width = 6 * base_width
    )
  }
  return(p_scatter)
}

.create_pathway_plots <- function(
    cellchat_object,
    seurat_object,
    pathways_to_show,
    condition,
    dirpath = NULL,
    output_format = "pdf",
    base_height = 1,
    base_width = 1,
    res = 300) {
  for (pathway in pathways_to_show) {
    if (!is.null(dirpath)) {
      chord_filename <- paste0(
        dirpath, "/", pathway, "_",
        condition, "_signaling_strength_chord.", output_format
      )
      .create_output_device(
        chord_filename,
        output_format = output_format,
        height = 5 * base_height,
        width = 5 * base_width,
        res = res
      )
    }

    CellChat::netVisual_aggregate(
      cellchat_object,
      signaling = pathway,
      title.space = 3,
      layout = "chord"
    )

    if (!is.null(dirpath)) {
      grDevices::dev.off()
    }

    p <- CellChat::netVisual_bubble(
      cellchat_object,
      signaling = pathway,
      remove.isolate = FALSE,
      font.size = 6
    )

    if (!is.null(dirpath)) {
      bubble_filename <- paste0(
        dirpath, "/", pathway, "_",
        condition, "_LR_bubble_plot.", output_format
      )
      bubble_height <- 1 + 0.3 * length(unique(p$data$interaction_name)) * base_height
      bubble_width <- 1 + 0.3 * length(unique(p$data$source.target)) * base_width

      .create_output_device(
        bubble_filename,
        output_format = output_format,
        height = bubble_height,
        width = bubble_width,
        res = res
      )
      print(p)
      grDevices::dev.off()
    }
  }
}

.create_heatmap_plots <- function(
    cellchat_object,
    condition,
    dirpath = NULL,
    output_format = "png",
    base_height = 1,
    base_width = 1,
    res = 300) {
  group_size <- as.numeric(table(cellchat_object@idents))
  pathway_num <- length(cellchat_object@netP$pathways)

  height_in <- 2 * 3 * (ceiling(pathway_num / 50)) * base_height
  width_in <- 8 / 3 * 3.5 * ceiling(length(group_size) / 30) * base_width
  heatmap_height <- 10 * ceiling(pathway_num / 50) * base_height
  heatmap_width <- 10 * ceiling(length(group_size) / 30) * base_width

  if (!is.null(dirpath)) {
    heatmap_filename <- paste0(
      dirpath, "/", condition, "_outgoing_incoming_signal.", output_format
    )
    .create_output_device(
      heatmap_filename,
      output_format = output_format,
      height = height_in,
      width = width_in,
      res = res
    )
  }

  ht1 <- CellChat::netAnalysis_signalingRole_heatmap(
    cellchat_object,
    pattern = "outgoing",
    height = heatmap_height,
    width = heatmap_width,
    font.size = 6
  )
  ht2 <- CellChat::netAnalysis_signalingRole_heatmap(
    cellchat_object,
    pattern = "incoming",
    height = heatmap_height,
    width = heatmap_width,
    font.size = 6
  )
  ComplexHeatmap::draw(ht1 + ht2)

  if (!is.null(dirpath)) {
    grDevices::dev.off()
  }
}

.create_circle_plots <- function(
    cellchat_object,
    condition,
    dirpath = NULL,
    output_format = "png",
    base_height = 1,
    base_width = 1,
    res = 300) {
  group_size <- as.numeric(table(cellchat_object@idents))
  numofcelltypes <- length(group_size)
  mat <- cellchat_object@net$weight

  height_in <- 2 * 3 * ceiling(numofcelltypes / 4) * base_height
  width_in <- 2 * 4 * 3 * base_width

  if (!is.null(dirpath)) {
    circle_filename <- paste0(
      dirpath, "/",
      condition, "_net_weight_per_celltype.", output_format
    )
    .create_output_device(
      circle_filename,
      output_format = output_format,
      height = height_in,
      width = width_in,
      res = res
    )
  }

  graphics::par(mfrow = c(ceiling(length(group_size) / 4), 4), xpd = TRUE)
  for (i in seq_len(nrow(mat))) {
    mat2 <- matrix(
      0,
      nrow = nrow(mat),
      ncol = ncol(mat),
      dimnames = dimnames(mat)
    )
    mat2[i, ] <- mat[i, ]
    CellChat::netVisual_circle(
      mat2,
      vertex.weight = group_size,
      weight.scale = TRUE,
      edge.weight.max = max(mat),
      vertex.label.cex = 1 + 4 / numofcelltypes,
      title.name = rownames(mat)[i]
    )
  }

  if (!is.null(dirpath)) {
    grDevices::dev.off()
  }
}

.create_bubble_plots <- function(
    cellchat_object,
    pathway,
    condition,
    dirpath = NULL,
    output_format = "png",
    base_height = 1,
    base_width = 1,
    res = 300) {
  p <- CellChat::netVisual_bubble(
    cellchat_object,
    signaling = pathway,
    remove.isolate = FALSE,
    font.size = 7
  )

  bubble_height <- (2 + 0.4 * length(unique(p$data$interaction_name))) * base_height
  bubble_width <- (2 + 25 / 300 * length(unique(p$data$source.target))) * base_width

  if (!is.null(dirpath)) {
    bubble_filename <- paste0(
      dirpath, "/", pathway, "_",
      condition, "_LR_bubble_plot.", output_format
    )
    .create_output_device(
      bubble_filename,
      output_format = output_format,
      height = bubble_height,
      width = bubble_width,
      res = res
    )
    print(p)
    grDevices::dev.off()
  }
  return(p)
}

.create_gene_plots <- function(
    cellchat_object,
    seurat_object,
    pathway,
    condition,
    dirpath = NULL,
    output_format = "pdf",
    species,
    base_height = 1,
    base_width = 1,
    res = 300) {
  pairLR <- CellChat::extractEnrichedLR(
    cellchat_object,
    signaling = pathway,
    geneLR.return = FALSE
  )
  LRs_genes <- unique(
    unlist(strsplit(split = "_", x = pairLR$interaction_name))
  )

  genes1 <- LRs_genes[LRs_genes %in% toupper(rownames(seurat_object))]
  genes2 <- LRs_genes[!(LRs_genes %in% toupper(rownames(seurat_object)))]
  genes22 <- sub(".?.?", "", genes2)
  LRs_genes <- c(genes1, genes22[genes22 %in% toupper(rownames(seurat_object))])

  if (species == "mouse") {
    genes <- rownames(cellchat_object@data)
    indices <- match(LRs_genes, toupper(genes))
    LRs_genes <- genes[indices]
  }
  p2 <- VlnPlot(
    object = seurat_object,
    features = LRs_genes,
    pt.size = -1,
    stack = TRUE
  )
  if (length(LRs_genes) == 1) {
    gene_width <- 4 * base_width
    gene_height <- (1 + .5 * length(levels(seurat_object))) * base_height
  } else {
    gene_width <- (2 + length(LRs_genes)) * base_width
    gene_height <- (1 + .5 * length(levels(seurat_object))) * base_height
  }

  if (!is.null(dirpath)) {
    gene_filename <- paste0(
      dirpath, "/", pathway, "_",
      condition, "_signaling_gene.", output_format
    )
    .create_output_device(
      gene_filename,
      output_format = output_format,
      height = gene_height,
      width = gene_width,
      res = res
    )
    print(p2)
    grDevices::dev.off()
  }
  return(p2)
}
