#' @title Visualize metacell partitions on a dimensionality reduction
#'
#' @md
#' @inheritParams CellDimPlot
#' @param srt A `Seurat` object with metacell results from `RunMetaCell()`.
#' @param show_cells Logical. If `TRUE`, the original single-cell points are
#' drawn as a semi-transparent background layer behind the metacell centroids.
#' @param color.by Metadata column in the metacell Seurat used to color
#' centroids. If `NULL`, metacell centroids use one fixed color.
#' @param palette_metacell Color palette for the metacell centroid layer.
#' Default is "Chinese".
#' @param palcolor_metacell Custom colors for the metacell centroid layer.
#' @param cell.alpha Alpha value for the original single-cell background layer.
#' @param cell.size Point size for the original single-cell background layer.
#' @param stroke Point border stroke width for metacell centroids.
#' @param show_metacell_size Whether to map `metacell_size` to centroid size.
#' @param metacell_size_range Point-size range used when
#' `show_metacell_size = TRUE`.
#' @param cell_param A named list of extra arguments passed to [CellDimPlot()]
#' for the original single-cell background layer.
#' @param metacell_param A named list of extra arguments passed to [CellDimPlot()]
#' for the metacell centroid/query layer.
#' @param return_layer Logical. If `TRUE`, returns a named list of ggplot2
#' layers/scales (`cells`, `fill_scale`, `centroids`, `scale_fill`,
#' `scale_size`, `labels`) instead of a complete plot.
#' @param ... Additional arguments passed to [geom_point()] for metacell
#' centroids.
#'
#' @return A `ggplot` object, or a named list of ggplot2 layers when
#' `return_layer = TRUE`.
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub, verbose = FALSE)
#' mc <- RunMetaCell(
#'   pancreas_sub,
#'   method = "supercell",
#'   gamma = 20
#' )
#'
#' MetaCellPlot(
#'   mc,
#'   group.by = "CellType",
#'   palette_metacell = "ChineseSet8"
#' )
#'
#' MetaCellPlot(
#'   mc,
#'   group.by = "CellType",
#'   reduction = "umap",
#'   palette = "ChineseSet8",
#'   show_cells = TRUE
#' )
#'
#' CellDimPlot(
#'   pancreas_sub, group.by = "CellType"
#' ) +
#'   MetaCellPlot(
#'     mc,
#'     group.by = "CellType",
#'     return_layer = TRUE,
#'     palette_metacell = "ChineseSet8"
#'   )
MetaCellPlot <- function(
  srt,
  reduction = NULL,
  show_cells = FALSE,
  group.by = NULL,
  color.by = NULL,
  dims = c(1, 2),
  label = FALSE,
  palette = "Chinese",
  palcolor = NULL,
  palette_metacell = "Chinese",
  palcolor_metacell = NULL,
  pt.size = 1.2,
  pt.alpha = 1,
  cell.alpha = 1,
  cell.size = 0.7,
  stroke = 0.5,
  show_metacell_size = TRUE,
  metacell_size_range = NULL,
  cell_param = list(),
  metacell_param = list(),
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  return_layer = FALSE,
  ...
) {
  original_srt <- srt@misc[["original_srt"]]
  membership <- srt@misc[["cell_membership"]]

  if (!is.null(original_srt) && !is.null(membership)) {
    embedding_src <- original_srt
    membership <- as.character(membership)
    mc_meta <- srt@meta.data
  } else {
    mc_tool <- srt@tools[["Metacell"]]
    if (is.null(mc_tool)) {
      log_message(
        "No metacell data found. Run {.fn RunMetaCell} first.",
        message_type = "error"
      )
    }
    membership <- mc_tool[["membership"]]
    if (is.null(membership) || length(membership) != ncol(srt)) {
      log_message(
        "{.code srt@tools[[\"Metacell\"]]$membership} is missing or does not match cell count",
        message_type = "error"
      )
    }
    membership <- as.character(membership)
    names(membership) <- colnames(srt)
    embedding_src <- srt
    mc_meta <- srt@meta.data
  }
  names(membership) <- colnames(embedding_src)

  if (is.null(reduction)) {
    reduction <- DefaultReduction(embedding_src)
  } else {
    reduction <- DefaultReduction(embedding_src, pattern = reduction)
  }
  if (!reduction %in% names(embedding_src@reductions)) {
    log_message(
      "{.val {reduction}} is not in the reduction names of the data source. Run {.fn standard_scop} on the original Seurat first.",
      message_type = "error"
    )
  }
  embedding <- Seurat::Embeddings(embedding_src, reduction = reduction)
  if (max(dims) > ncol(embedding)) {
    dims <- seq_len(min(2, ncol(embedding)))
  }
  if (length(dims) < 2 || max(dims) > ncol(embedding)) {
    log_message(
      "Reduction {.val {reduction}} does not have enough dimensions",
      message_type = "error"
    )
  }
  dim_labels <- colnames(embedding)[dims]
  if (any(is.na(dim_labels)) || any(!nzchar(dim_labels))) {
    dim_labels <- paste0("Dim", dims)
  }

  embedding <- embedding[names(membership), , drop = FALSE]
  centroids <- do.call(
    rbind,
    lapply(
      split(as.data.frame(embedding), membership),
      function(df) colMeans(df[, dims, drop = FALSE])
    )
  )
  centroids <- as.data.frame(centroids)
  colnames(centroids) <- c("dim1", "dim2")
  centroids[["metacell"]] <- rownames(centroids)
  centroids[["size"]] <- as.integer(
    table(membership)[centroids[["metacell"]]]
  )

  if (!is.null(group.by)) {
    color.by <- group.by
  }
  metacell_size_range <- metacell_size_range %||%
    c(pt.size, pt.size * 2.5)

  centroid_counts <- Matrix::Matrix(
    0,
    nrow = 1,
    ncol = nrow(centroids),
    sparse = TRUE
  )
  rownames(centroid_counts) <- "MetaCellPlot"
  colnames(centroid_counts) <- centroids[["metacell"]]
  centroid_srt <- Seurat::CreateSeuratObject(counts = centroid_counts)
  centroid_srt[["metacell_size"]] <- centroids[["size"]]
  centroid_srt[["metacell"]] <- centroids[["metacell"]]
  if (!is.null(color.by) && color.by %in% colnames(mc_meta)) {
    centroid_srt@meta.data[[color.by]] <- mc_meta[
      colnames(centroid_srt),
      color.by,
      drop = TRUE
    ]
  } else if (
    !is.null(color.by) &&
      color.by %in% c("Metacell_id", "metacell", "metacell_id")
  ) {
    centroid_srt@meta.data[[color.by]] <- colnames(centroid_srt)
  }
  centroid_embeddings <- as.matrix(centroids[, c("dim1", "dim2"), drop = FALSE])
  rownames(centroid_embeddings) <- centroids[["metacell"]]
  colnames(centroid_embeddings) <- paste0("MetaCellPlot_", seq_len(2))
  centroid_srt[["MetaCellPlot"]] <- SeuratObject::CreateDimReducObject(
    embeddings = centroid_embeddings,
    key = "MetaCellPlot_",
    assay = SeuratObject::DefaultAssay(centroid_srt)
  )

  cell_group <- color.by
  centroid_group <- color.by
  cell_palcolor <- palcolor
  centroid_palcolor <- palcolor_metacell
  embedding_plot <- embedding_src
  if (is.null(color.by) || !color.by %in% colnames(embedding_plot@meta.data)) {
    cell_group <- "MetaCellPlot_cells"
    embedding_plot@meta.data[[cell_group]] <- factor("Cells")
    cell_palcolor <- "grey65"
  }
  if (is.null(color.by) || !color.by %in% colnames(centroid_srt@meta.data)) {
    centroid_group <- "MetaCellPlot_metacells"
    centroid_srt@meta.data[[centroid_group]] <- factor("Metacells")
    if (is.null(centroid_palcolor)) {
      centroid_palcolor <- unname(
        palette_colors(palette = palette_metacell, n = 1L)[1]
      )
    }
  }
  if (
    !is.null(color.by) &&
      color.by %in% colnames(embedding_plot@meta.data) &&
      color.by %in% colnames(centroid_srt@meta.data)
  ) {
    color_levels <- unique(c(
      as.character(embedding_plot@meta.data[[color.by]]),
      as.character(centroid_srt@meta.data[[color.by]])
    ))
    color_levels <- color_levels[!is.na(color_levels)]
    embedding_plot@meta.data[[color.by]] <- factor(
      as.character(embedding_plot@meta.data[[color.by]]),
      levels = color_levels
    )
    centroid_srt@meta.data[[color.by]] <- factor(
      as.character(centroid_srt@meta.data[[color.by]]),
      levels = color_levels
    )
  }

  cell_args <- utils::modifyList(
    list(
      srt = embedding_plot,
      group.by = cell_group,
      reduction = reduction,
      dims = dims,
      pt.size = cell.size,
      pt.alpha = cell.alpha,
      palette = palette,
      palcolor = cell_palcolor,
      show_stat = FALSE,
      label = FALSE,
      legend.position = if (identical(cell_group, "MetaCellPlot_cells")) {
        "none"
      } else {
        legend.position
      },
      legend.direction = legend.direction,
      theme_use = theme_use,
      theme_args = theme_args,
      force = TRUE,
      verbose = FALSE
    ),
    cell_param
  )
  metacell_args <- utils::modifyList(
    list(
      srt = centroid_srt,
      group.by = centroid_group,
      reduction = "MetaCellPlot",
      dims = c(1, 2),
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      palette = palette_metacell,
      palcolor = centroid_palcolor,
      show_stat = FALSE,
      label = FALSE,
      legend.position = legend.position,
      legend.direction = legend.direction,
      theme_use = theme_use,
      theme_args = theme_args,
      force = TRUE,
      verbose = FALSE
    ),
    metacell_param
  )

  p_metacell <- do.call(CellDimPlot, metacell_args)
  p_metacell <- p_metacell + ggplot2::guides(color = "none", fill = "none")
  metacell_data <- p_metacell$data
  metacell_data[["metacell_size"]] <- centroid_srt@meta.data[
    rownames(metacell_data),
    "metacell_size",
    drop = TRUE
  ]
  metacell_pt_size <- metacell_args[["pt.size"]]
  metacell_pt_alpha <- metacell_args[["pt.alpha"]]
  metacell_groups <- levels(metacell_data[["group.by"]])
  metacell_colors <- palette_colors(
    metacell_groups,
    palette = palette_metacell,
    palcolor = centroid_palcolor,
    NA_keep = TRUE
  )
  centroid_layer <- if (isTRUE(show_metacell_size)) {
    ggplot2::geom_point(
      data = metacell_data,
      ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        fill = .data[["group.by"]],
        size = .data[["metacell_size"]]
      ),
      shape = 21,
      color = "black",
      stroke = stroke,
      alpha = metacell_pt_alpha,
      ...
    )
  } else {
    ggplot2::geom_point(
      data = metacell_data,
      ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        fill = .data[["group.by"]]
      ),
      shape = 21,
      color = "black",
      size = metacell_pt_size,
      stroke = stroke,
      alpha = metacell_pt_alpha,
      ...
    )
  }
  metacell_fill_scale <- ggplot2::scale_fill_manual(
    name = if (!is.null(color.by)) {
      paste0("Metacell: ", color.by)
    } else {
      "Metacell"
    },
    values = metacell_colors,
    na.value = "grey80"
  )
  metacell_size_scale <- if (isTRUE(show_metacell_size)) {
    ggplot2::scale_size_continuous(
      name = "Metacell size",
      range = metacell_size_range
    )
  } else {
    NULL
  }
  cell_layer <- NULL

  label_layer <- NULL
  if (isTRUE(label)) {
    label_layer <- ggrepel::geom_text_repel(
      data = centroids,
      mapping = ggplot2::aes(
        x = .data[["dim1"]],
        y = .data[["dim2"]],
        label = .data[["metacell"]]
      ),
      size = 3,
      max.overlaps = 50,
      show.legend = FALSE
    )
  }

  if (isTRUE(return_layer)) {
    return(list(
      cells = cell_layer,
      fill_scale = ggnewscale::new_scale_fill(),
      centroids = centroid_layer,
      scale_fill = metacell_fill_scale,
      scale_size = metacell_size_scale,
      labels = label_layer
    ))
  }

  if (isTRUE(show_cells)) {
    p <- do.call(CellDimPlot, cell_args)
  } else {
    p <- p_metacell
  }

  p <- p +
    ggnewscale::new_scale_fill() +
    centroid_layer +
    metacell_fill_scale
  if (!is.null(metacell_size_scale)) {
    p <- p + metacell_size_scale
  }

  if (!is.null(label_layer)) {
    p <- p + label_layer
  }

  p
}
