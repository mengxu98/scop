#' @title Heatmap plot for dynamic features along lineages
#'
#' @md
#' @inheritParams GroupHeatmap
#' @param srt A Seurat object.
#' @param lineages A character vector specifying the lineages to plot.
#' @param features A character vector specifying the features to plot.
#' By default, this parameter is set to NULL, and the dynamic features will be determined by the parameters
#' \code{min_expcells}, \code{r.sq}, \code{dev.expl}, \code{padjust} and \code{num_intersections}.
#' @param use_fitted A logical indicating whether to use fitted values. Default is FALSE.
#' @param border A logical indicating whether to add a border to the heatmap. Default is TRUE.
#' @param flip A logical indicating whether to flip the heatmap. Default is FALSE.
#' @param min_expcells A numeric value specifying the minimum number of expected cells. Default is 20.
#' @param r.sq A numeric value specifying the R-squared threshold. Default is 0.2.
#' @param dev.expl A numeric value specifying the deviance explained threshold. Default is 0.2.
#' @param padjust A numeric value specifying the p-value adjustment threshold. Default is 0.05.
#' @param num_intersections This parameter is a numeric vector used to determine the number of intersections among lineages. It helps in selecting which dynamic features will be used. By default, when this parameter is set to NULL, all dynamic features that pass the specified threshold will be used for each lineage.
#' @param cell_density A numeric value is used to define the cell density within each cell bin. By default, this parameter is set to 1, which means that all cells will be included within each cell bin.
#' @param cell_bins A numeric value specifying the number of cell bins. Default is 100.
#' @param order_by A character vector specifying the order of the heatmap. Default is "peaktime".
#' @param family A character specifying the model used to calculate the dynamic features if needed. By default, this parameter is set to NULL, and the appropriate family will be automatically determined.
#' @param cluster_features_by A character vector specifying which lineage to use when clustering features. By default, this parameter is set to NULL, which means that all lineages will be used.
#' @param pseudotime_label A numeric vector specifying the pseudotime label. Default is NULL.
#' @param pseudotime_label_color A character string specifying the pseudotime label color. Default is "black".
#' @param pseudotime_label_linetype A numeric value specifying the pseudotime label line type. Default is 2.
#' @param pseudotime_label_linewidth A numeric value specifying the pseudotime label line width. Default is 3.
#' @param pseudotime_palette A character vector specifying the color palette to use for pseudotime.
#' @param pseudotime_palcolor A list specifying the colors to use for the pseudotime in the heatmap.
#' @param separate_annotation A character vector of names of annotations to be displayed in separate annotation blocks. Each name should match a column name in the metadata of the Seurat object.
#' @param separate_annotation_palette A character vector specifying the color palette to use for separate annotations.
#' @param separate_annotation_palcolor A list specifying the colors to use for each level of the separate annotations.
#' @param separate_annotation_params A list of other parameters to be passed to the HeatmapAnnotation function when creating the separate annotation blocks.
#' @param reverse_ht A logical indicating whether to reverse the heatmap. Default is NULL.
#'
#' @seealso \code{\link{RunDynamicFeatures}} \code{\link{RunDynamicEnrichment}}
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' # pancreas_sub <- AnnotateFeatures(
#' #   srt = pancreas_sub,
#' #   species = "Mus_musculus",
#' #   db = c("TF", "CSPA")
#' # )
#'
#' pancreas_sub <- RunSlingshot(
#'   srt = pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP"
#' )
#' pancreas_sub <- RunDynamicFeatures(
#'   srt = pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   n_candidates = 200
#' )
#'
#' ht1 <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   lineages = "Lineage1",
#'   n_split = 5,
#'   split_method = "kmeans-peaktime",
#'   cell_annotation = "SubCellType"
#' )
#' ht1$plot
#'
#' panel_fix(ht1$plot, raster = TRUE, dpi = 50)
#'
#' ht2 <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   lineages = "Lineage1",
#'   features = c(
#'     "Sox9",
#'     "Neurod2",
#'     "Isl1",
#'     "Rbp4",
#'     "Pyy", "S_score", "G2M_score"
#'   ),
#'   cell_annotation = "SubCellType"
#' )
#' ht2$plot
#'
#' panel_fix(
#'   ht2$plot,
#'   height = 5,
#'   width = 5,
#'   raster = TRUE,
#'   dpi = 50
#' )
#'
#' ht3 <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   n_split = 5,
#'   split_method = "kmeans",
#'   cluster_rows = TRUE,
#'   cell_annotation = "SubCellType"
#' )
#' ht3$plot
#'
#' ht4 <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   reverse_ht = "Lineage1",
#'   use_fitted = TRUE,
#'   n_split = 6,
#'   split_method = "mfuzz",
#'   heatmap_palette = "viridis",
#'   cell_annotation = c(
#'     "SubCellType", "Phase", "G2M_score"
#'   ),
#'   cell_annotation_palette = c(
#'     "Paired", "simspec", "Purples"
#'   ),
#'   separate_annotation = list(
#'     "SubCellType", c("Arxes1", "Ncoa2")
#'   ),
#'   separate_annotation_palette = c(
#'     "Paired", "Set1"
#'   ),
#'   separate_annotation_params = list(
#'     height = grid::unit(10, "mm")
#'   ),
#'   feature_annotation = c("TF", "CSPA"),
#'   feature_annotation_palcolor = list(
#'     c("gold", "steelblue"),
#'     c("forestgreen")
#'   ),
#'   pseudotime_label = 25,
#'   pseudotime_label_color = "red"
#' )
#' ht4$plot
#'
#' ht5 <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   reverse_ht = "Lineage1",
#'   use_fitted = TRUE,
#'   n_split = 6,
#'   split_method = "mfuzz",
#'   heatmap_palette = "viridis",
#'   cell_annotation = c(
#'     "SubCellType", "Phase", "G2M_score"
#'   ),
#'   cell_annotation_palette = c(
#'     "Paired", "simspec", "Purples"
#'   ),
#'   separate_annotation = list(
#'     "SubCellType", c("Arxes1", "Ncoa2")
#'   ),
#'   separate_annotation_palette = c("Paired", "Set1"),
#'   separate_annotation_params = list(width = grid::unit(10, "mm")),
#'   feature_annotation = c("TF", "CSPA"),
#'   feature_annotation_palcolor = list(
#'     c("gold", "steelblue"),
#'     c("forestgreen")
#'   ),
#'   pseudotime_label = 25,
#'   pseudotime_label_color = "red",
#'   flip = TRUE, column_title_rot = 45
#' )
#' ht5$plot
#'
#' ht6 <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   reverse_ht = "Lineage1",
#'   cell_annotation = "SubCellType",
#'   n_split = 5,
#'   split_method = "mfuzz",
#'   species = "Mus_musculus",
#'   db = "GO_BP",
#'   anno_terms = TRUE,
#'   anno_keys = TRUE,
#'   anno_features = TRUE
#' )
#' ht6$plot
DynamicHeatmap <- function(
    srt,
    lineages,
    features = NULL,
    use_fitted = FALSE,
    border = TRUE,
    flip = FALSE,
    min_expcells = 20,
    r.sq = 0.2,
    dev.expl = 0.2,
    padjust = 0.05,
    num_intersections = NULL,
    cell_density = 1,
    cell_bins = 100,
    order_by = c("peaktime", "valleytime"),
    layer = "counts",
    assay = NULL,
    exp_method = c("zscore", "raw", "fc", "log2fc", "log1p"),
    exp_legend_title = NULL,
    limits = NULL,
    lib_normalize = identical(layer, "counts"),
    libsize = NULL,
    family = NULL,
    cluster_features_by = NULL,
    cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_names_side = ifelse(flip, "left", "right"),
    column_names_side = ifelse(flip, "bottom", "top"),
    row_names_rot = 0,
    column_names_rot = 90,
    row_title = NULL,
    column_title = NULL,
    row_title_side = "left",
    column_title_side = "top",
    row_title_rot = 0,
    column_title_rot = ifelse(flip, 90, 0),
    feature_split = NULL,
    feature_split_by = NULL,
    n_split = NULL,
    split_order = NULL,
    split_method = c(
      "mfuzz",
      "kmeans",
      "kmeans-peaktime",
      "hclust",
      "hclust-peaktime"
    ),
    decreasing = FALSE,
    fuzzification = NULL,
    anno_terms = FALSE,
    anno_keys = FALSE,
    anno_features = FALSE,
    terms_width = grid::unit(4, "in"),
    terms_fontsize = 8,
    keys_width = grid::unit(2, "in"),
    keys_fontsize = c(6, 10),
    features_width = grid::unit(2, "in"),
    features_fontsize = c(6, 10),
    IDtype = "symbol",
    species = "Homo_sapiens",
    db_update = FALSE,
    db_version = "latest",
    db_combine = FALSE,
    convert_species = FALSE,
    Ensembl_version = 103,
    mirror = NULL,
    db = "GO_BP",
    TERM2GENE = NULL,
    TERM2NAME = NULL,
    minGSSize = 10,
    maxGSSize = 500,
    GO_simplify = FALSE,
    GO_simplify_cutoff = "p.adjust < 0.05",
    simplify_method = "Wang",
    simplify_similarityCutoff = 0.7,
    pvalueCutoff = NULL,
    padjustCutoff = 0.05,
    topTerm = 5,
    show_termid = FALSE,
    topWord = 20,
    words_excluded = NULL,
    nlabel = 20,
    features_label = NULL,
    label_size = 10,
    label_color = "black",
    pseudotime_label = NULL,
    pseudotime_label_color = "black",
    pseudotime_label_linetype = 2,
    pseudotime_label_linewidth = 3,
    heatmap_palette = "viridis",
    heatmap_palcolor = NULL,
    pseudotime_palette = "cividis",
    pseudotime_palcolor = NULL,
    feature_split_palette = "simspec",
    feature_split_palcolor = NULL,
    cell_annotation = NULL,
    cell_annotation_palette = "Paired",
    cell_annotation_palcolor = NULL,
    cell_annotation_params = if (flip) {
      list(width = grid::unit(5, "mm"))
    } else {
      list(height = grid::unit(5, "mm"))
    },
    feature_annotation = NULL,
    feature_annotation_palette = "Dark2",
    feature_annotation_palcolor = NULL,
    feature_annotation_params = if (flip) {
      list(height = grid::unit(5, "mm"))
    } else {
      list(width = grid::unit(5, "mm"))
    },
    separate_annotation = NULL,
    separate_annotation_palette = "Paired",
    separate_annotation_palcolor = NULL,
    separate_annotation_params = if (flip) {
      list(width = grid::unit(10, "mm"))
    } else {
      list(height = grid::unit(10, "mm"))
    },
    reverse_ht = NULL,
    use_raster = NULL,
    raster_device = "png",
    raster_by_magick = FALSE,
    height = NULL,
    width = NULL,
    units = "inch",
    seed = 11,
    ht_params = list()) {
  set.seed(seed)
  if (isTRUE(raster_by_magick)) {
    check_r("magick")
  }

  split_method <- match.arg(split_method)
  order_by <- match.arg(order_by)
  data_nm <- c(
    ifelse(isTRUE(lib_normalize), "normalized", ""), layer
  )
  data_nm <- paste(data_nm[data_nm != ""], collapse = " ")
  if (length(exp_method) == 1 && is.function(exp_method)) {
    exp_name <- paste0(
      as.character(x = formals()$exp_method),
      "(",
      data_nm,
      ")"
    )
  } else {
    exp_method <- match.arg(exp_method)
    exp_name <- paste0(exp_method, "(", data_nm, ")")
  }
  exp_name <- exp_legend_title %||% exp_name

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (any(!lineages %in% colnames(srt@meta.data))) {
    lineages_missing <- lineages[!lineages %in% colnames(srt@meta.data)]
    for (l in lineages_missing) {
      if (paste0("DynamicFeatures_", l) %in% names(srt@tools)) {
        pseudotime <- srt@tools[[paste0("DynamicFeatures_", l)]][["lineages"]]
        srt@meta.data[[l]] <- srt@meta.data[[pseudotime]]
      } else {
        log_message(
          "lineages: ", l, " is not in the meta data of the Seurat object",
          message_type = "error"
        )
      }
    }
  }

  if (is.null(feature_split_by)) {
    feature_split_by <- lineages
  }

  if (any(!feature_split_by %in% lineages)) {
    log_message(
      "'feature_split_by' must be a subset of lineages.",
      message_type = "error"
    )
  }
  if (
    !split_method %in%
      c("mfuzz", "kmeans", "kmeans-peaktime", "hclust", "hclust-peaktime")
  ) {
    log_message(
      "'split_method' must be one of 'mfuzz', 'kmeans', 'kmeans-peaktime', 'hclust', 'hclust-peaktime'.",
      message_type = "error"
    )
  }

  if (!is.null(feature_split) && is.null(names(feature_split))) {
    log_message(
      "'feature_split' must be named.",
      message_type = "error"
    )
  }
  if (!is.null(feature_split) && !is.factor(feature_split)) {
    feature_split <- factor(feature_split, levels = unique(feature_split))
  }
  if (is.numeric(pseudotime_label)) {
    if (length(pseudotime_label_color) == 1) {
      pseudotime_label_color <- rep(
        pseudotime_label_color,
        length(pseudotime_label)
      )
    }
    if (length(pseudotime_label_linetype) == 1) {
      pseudotime_label_linetype <- rep(
        pseudotime_label_linetype,
        length(pseudotime_label)
      )
    }
    if (length(pseudotime_label_linewidth) == 1) {
      pseudotime_label_linewidth <- rep(
        pseudotime_label_linewidth,
        length(pseudotime_label)
      )
    }
    npal <- unique(c(
      length(pseudotime_label),
      length(pseudotime_label_color),
      length(pseudotime_label_linetype),
      length(pseudotime_label_linewidth)
    ))
    if (length(npal[npal != 0]) > 1) {
      log_message(
        "Parameters for the pseudotime_label must be the same length!",
        message_type = "error"
      )
    }
  }

  if (!is.null(cell_annotation)) {
    if (length(cell_annotation_palette) == 1) {
      cell_annotation_palette <- rep(
        cell_annotation_palette,
        length(cell_annotation)
      )
    }
    if (length(cell_annotation_palcolor) == 1) {
      cell_annotation_palcolor <- rep(
        cell_annotation_palcolor,
        length(cell_annotation)
      )
    }
    npal <- unique(c(
      length(cell_annotation_palette),
      length(cell_annotation_palcolor),
      length(cell_annotation)
    ))
    if (length(npal[npal != 0]) > 1) {
      log_message(
        "cell_annotation_palette and cell_annotation_palcolor must be the same length as cell_annotation",
        message_type = "error"
      )
    }
    if (
      any(
        !cell_annotation %in%
          c(colnames(srt@meta.data), rownames(srt@assays[[assay]]))
      )
    ) {
      log_message(
        "cell_annotation: ",
        paste0(
          cell_annotation[
            !cell_annotation %in%
              c(colnames(srt@meta.data), rownames(srt@assays[[assay]]))
          ],
          collapse = ","
        ),
        " is not in the Seurat object."
      )
    }
  }

  if (!is.null(feature_annotation)) {
    if (length(feature_annotation_palette) == 1) {
      feature_annotation_palette <- rep(
        feature_annotation_palette,
        length(feature_annotation)
      )
    }
    if (length(feature_annotation_palcolor) == 1) {
      feature_annotation_palcolor <- rep(
        feature_annotation_palcolor,
        length(feature_annotation)
      )
    }
    npal <- unique(c(
      length(feature_annotation_palette),
      length(feature_annotation_palcolor),
      length(feature_annotation)
    ))
    if (length(npal[npal != 0]) > 1) {
      log_message(
        "feature_annotation_palette and feature_annotation_palcolor must be the same length as feature_annotation",
        message_type = "error"
      )
    }
    srt_assay_features <- GetFeaturesData(srt@assays[[assay]])
    if (
      any(!feature_annotation %in% colnames(srt_assay_features))
    ) {
      log_message(
        "feature_annotation: ",
        paste0(
          feature_annotation[
            !feature_annotation %in% colnames(srt_assay_features)
          ],
          collapse = ","
        ),
        " is not in the meta data of the ",
        assay,
        " assay in the Seurat object.",
        message_type = "error"
      )
    }
  }

  if (!is.null(separate_annotation)) {
    if (length(separate_annotation_palette) == 1) {
      separate_annotation_palette <- rep(
        separate_annotation_palette,
        length(separate_annotation)
      )
    }
    if (length(separate_annotation_palcolor) == 1) {
      separate_annotation_palcolor <- rep(
        separate_annotation_palcolor,
        length(separate_annotation)
      )
    }
    npal <- unique(c(
      length(separate_annotation_palette),
      length(separate_annotation_palcolor),
      length(separate_annotation)
    ))
    if (length(npal[npal != 0]) > 1) {
      log_message(
        "separate_annotation_palette and separate_annotation_palcolor must be the same length as separate_annotation",
        message_type = "error"
      )
    }
    if (
      any(
        !unique(unlist(separate_annotation)) %in%
          c(colnames(srt@meta.data), rownames(srt@assays[[assay]]))
      )
    ) {
      log_message(
        "separate_annotation: ",
        paste0(
          unique(unlist(separate_annotation))[
            !unique(unlist(separate_annotation)) %in%
              c(colnames(srt@meta.data), rownames(srt@assays[[assay]]))
          ],
          collapse = ","
        ),
        " is not in the Seurat object.",
        message_type = "error"
      )
    }
  }

  if (length(width) == 1) {
    width <- rep(width, length(lineages))
  }
  if (length(height) == 1) {
    height <- rep(height, length(lineages))
  }
  if (length(width) >= 1) {
    names(width) <- lineages
  }
  if (length(height) >= 1) {
    names(height) <- lineages
  }

  column_split <- NULL
  if (isTRUE(flip)) {
    cluster_rows_raw <- cluster_rows
    cluster_columns_raw <- cluster_columns
    cluster_row_slices_raw <- cluster_row_slices
    cluster_column_slices_raw <- cluster_column_slices
    cluster_rows <- cluster_columns_raw
    cluster_columns <- cluster_rows_raw
    cluster_row_slices <- cluster_column_slices_raw
    cluster_column_slices <- cluster_row_slices_raw
  }

  cell_union <- unique(colnames(srt@assays[[1]])[apply(
    srt@meta.data[, lineages, drop = FALSE],
    1,
    function(x) !all(is.na(x))
  )])
  Pseudotime_assign <- Matrix::rowMeans(
    srt@meta.data[cell_union, lineages, drop = FALSE],
    na.rm = TRUE
  )
  cell_metadata <- cbind.data.frame(
    data.frame(row.names = cell_union, cells = cell_union),
    Pseudotime_assign = Pseudotime_assign,
    srt@meta.data[cell_union, lineages, drop = FALSE]
  )
  if (cell_density != 1) {
    bins <- cut(
      Pseudotime_assign,
      breaks = seq(
        min(Pseudotime_assign, na.rm = TRUE),
        max(Pseudotime_assign, na.rm = TRUE),
        length.out = cell_bins
      ),
      include.lowest = TRUE
    )
    ncells <- ceiling(max(table(bins), na.rm = TRUE) * cell_density)
    log_message("ncell/bin=", ncells, "(", cell_bins, "bins)")
    cell_keep <- unlist(sapply(levels(bins), function(x) {
      cells <- names(Pseudotime_assign)[bins == x]
      out <- sample(cells, size = min(length(cells), ncells))
      return(out)
    }))
    cell_metadata <- cell_metadata[cell_keep, , drop = FALSE]
  }
  cell_order_list <- list()
  for (l in lineages) {
    if (
      !is.null(reverse_ht) && (l %in% lineages[reverse_ht] || l %in% reverse_ht)
    ) {
      cell_metadata_sub <- stats::na.omit(cell_metadata[, l, drop = FALSE])
      cell_metadata_sub <- cell_metadata_sub[
        order(cell_metadata_sub[[l]], decreasing = TRUE), ,
        drop = FALSE
      ]
    } else {
      cell_metadata_sub <- stats::na.omit(cell_metadata[, l, drop = FALSE])
      cell_metadata_sub <- cell_metadata_sub[
        order(cell_metadata_sub[[l]], decreasing = FALSE), ,
        drop = FALSE
      ]
    }
    cell_order_list[[l]] <- paste0(rownames(cell_metadata_sub), l)
  }

  if (!is.null(cell_annotation)) {
    cell_metadata <- cbind.data.frame(
      cell_metadata,
      cbind.data.frame(
        srt@meta.data[
          rownames(cell_metadata),
          c(intersect(cell_annotation, colnames(srt@meta.data))),
          drop = FALSE
        ],
        Matrix::t(GetAssayData5(
          srt,
          assay = assay,
          layer = "data"
        )[
          intersect(cell_annotation, rownames(srt@assays[[assay]])) %||%
            integer(),
          rownames(cell_metadata),
          drop = FALSE
        ])
      )
    )
  }

  dynamic <- list()
  if (is.null(features)) {
    for (l in lineages) {
      DynamicFeatures <- srt@tools[[paste0("DynamicFeatures_", l)]][[
        "DynamicFeatures"
      ]]
      if (is.null(DynamicFeatures)) {
        log_message(
          "DynamicFeatures result for ",
          l,
          " found in the srt object. Should perform RunDynamicFeatures first!",
          message_type = "error"
        )
      }
      DynamicFeatures <- DynamicFeatures[
        DynamicFeatures$exp_ncells > min_expcells &
          DynamicFeatures$r.sq > r.sq &
          DynamicFeatures$dev.expl > dev.expl &
          DynamicFeatures$padjust < padjust, ,
        drop = FALSE
      ]
      dynamic[[l]] <- DynamicFeatures
      features <- c(features, DynamicFeatures[["features"]])
    }
    log_message(
      length(unique(features)),
      " features from ",
      paste0(lineages, collapse = ","),
      " passed the threshold (exp_ncells>",
      min_expcells,
      " & r.sq>",
      r.sq,
      " & dev.expl>",
      dev.expl,
      " & padjust<",
      padjust,
      "): \n",
      paste0(utils::head(features, 10), collapse = ","),
      "..."
    )
  } else {
    for (l in lineages) {
      DynamicFeatures <- srt@tools[[paste0("DynamicFeatures_", l)]][[
        "DynamicFeatures"
      ]]
      if (is.null(DynamicFeatures)) {
        srt <- RunDynamicFeatures(
          srt,
          lineages = l,
          features = features,
          assay = assay,
          layer = layer,
          family = family,
          libsize = libsize
        )
        DynamicFeatures <- srt@tools[[paste0("DynamicFeatures_", l)]][[
          "DynamicFeatures"
        ]]
      }
      if (any(!features %in% rownames(DynamicFeatures))) {
        srt <- RunDynamicFeatures(
          srt,
          lineages = l,
          features = features[!features %in% rownames(DynamicFeatures)],
          assay = assay,
          layer = layer,
          family = family,
          libsize = libsize
        )
        DynamicFeatures <- rbind(
          DynamicFeatures,
          srt@tools[[paste0("DynamicFeatures_", l)]][["DynamicFeatures"]]
        )
      }
      DynamicFeatures <- DynamicFeatures[features, , drop = FALSE]
      dynamic[[l]] <- DynamicFeatures
    }
  }

  all_calculated <- Reduce(
    intersect,
    lapply(
      lineages,
      function(l) {
        rownames(
          srt@tools[[paste0("DynamicFeatures_", l)]][["DynamicFeatures"]]
        )
      }
    )
  )

  features_tab <- table(features)
  features <- names(features_tab)[which(
    features_tab %in% (num_intersections %||% seq_along(lineages))
  )]
  if (!all(features %in% all_calculated)) {
    log_message(
      "Some features were missing in at least one lineage: \n",
      paste0(utils::head(features[!features %in% all_calculated], 10), collapse = ","),
      "..."
    )
    features <- intersect(features, all_calculated)
  }

  gene <- features[features %in% rownames(srt@assays[[assay]])]
  meta <- features[features %in% colnames(srt@meta.data)]
  if (length(gene) == 0 && length(meta) == 0) {
    log_message(
      "No dynamic features found in the meta.data or in the assay: ", assay,
      message_type = "error"
    )
  }
  feature_metadata <- data.frame(
    row.names = features,
    features = features
  )
  for (l in lineages) {
    feature_metadata[rownames(dynamic[[l]]), paste0(l, order_by)] <- dynamic[[
      l
    ]][, order_by]
    feature_metadata[
      rownames(dynamic[[l]]),
      paste0(l, "exp_ncells")
    ] <- dynamic[[l]][, "exp_ncells"]
    feature_metadata[rownames(dynamic[[l]]), paste0(l, "r.sq")] <- dynamic[[
      l
    ]][, "r.sq"]
    feature_metadata[rownames(dynamic[[l]]), paste0(l, "dev.expl")] <- dynamic[[
      l
    ]][, "dev.expl"]
    feature_metadata[rownames(dynamic[[l]]), paste0(l, "padjust")] <- dynamic[[
      l
    ]][, "padjust"]
  }
  feature_metadata[, order_by] <- apply(
    feature_metadata[, paste0(lineages, order_by), drop = FALSE],
    1,
    max,
    na.rm = TRUE
  )
  feature_metadata <- feature_metadata[
    order(feature_metadata[, order_by], decreasing = decreasing), ,
    drop = FALSE
  ]
  feature_metadata <- feature_metadata[
    rownames(feature_metadata) %in% features, ,
    drop = FALSE
  ]
  features <- rownames(feature_metadata)
  if (!is.null(feature_annotation)) {
    feature_metadata <- cbind.data.frame(
      feature_metadata,
      GetFeaturesData(srt, assay = assay)[
        rownames(feature_metadata),
        feature_annotation,
        drop = FALSE
      ]
    )
  }

  if (isTRUE(use_fitted)) {
    mat_list <- list()
    for (l in lineages) {
      fitted_matrix <- srt@tools[[paste0("DynamicFeatures_", l)]][[
        "fitted_matrix"
      ]][, -1]
      rownames(fitted_matrix) <- paste0(rownames(fitted_matrix), l)
      mat_list[[l]] <- Matrix::t(fitted_matrix[, features])
    }
    mat_raw <- do.call(cbind, mat_list)
  } else {
    mat_list <- list()
    y_libsize <- Matrix::colSums(
      GetAssayData5(
        srt,
        assay = assay,
        layer = "counts"
      )
    )

    for (l in lineages) {
      cells <- gsub(pattern = l, replacement = "", x = cell_order_list[[l]])
      mat_tmp <- Matrix::as.matrix(
        rbind(
          GetAssayData5(
            srt,
            assay = assay,
            layer = layer
          )[gene, cells, drop = FALSE],
          Matrix::t(
            srt@meta.data[cells, meta, drop = FALSE]
          )
        )
      )[features, , drop = FALSE]
      if (isTRUE(lib_normalize) && min(mat_tmp, na.rm = TRUE) >= 0) {
        if (!is.null(libsize)) {
          libsize_use <- libsize
        } else {
          libsize_use <- y_libsize[colnames(mat_tmp)]
          isfloat <- any(libsize_use %% 1 != 0, na.rm = TRUE)
          if (isTRUE(isfloat)) {
            libsize_use <- rep(1, length(libsize_use))
            log_message(
              "The values in the 'counts' layer are non-integer. Set the library size to 1.",
              message_type = "warning"
            )
          }
        }
        mat_tmp[gene, ] <- Matrix::t(
          Matrix::t(
            mat_tmp[gene, , drop = FALSE]
          ) /
            libsize_use *
            stats::median(y_libsize)
        )
      }
      colnames(mat_tmp) <- paste0(colnames(mat_tmp), l)
      mat_list[[l]] <- mat_tmp
    }
    mat_raw <- do.call(cbind, mat_list)
  }

  mat <- matrix_process(mat_raw, method = exp_method)
  mat[is.infinite(mat)] <- max(abs(mat[!is.infinite(mat)]), na.rm = TRUE) *
    ifelse(mat[is.infinite(mat)] > 0, 1, -1)
  mat[is.na(mat)] <- mean(mat, na.rm = TRUE)

  mat_split <- mat[, unlist(cell_order_list[feature_split_by]), drop = FALSE]

  if (is.null(limits)) {
    if (!is.function(exp_method) && exp_method %in% c("zscore", "log2fc")) {
      b <- ceiling(
        min(abs(stats::quantile(mat, c(0.01, 0.99), na.rm = TRUE)), na.rm = TRUE) * 2
      ) /
        2
      colors <- circlize::colorRamp2(
        seq(-b, b, length = 100),
        palette_scop(
          palette = heatmap_palette,
          palcolor = heatmap_palcolor
        )
      )
    } else {
      b <- stats::quantile(mat, c(0.01, 0.99), na.rm = TRUE)
      colors <- circlize::colorRamp2(
        seq(b[1], b[2], length = 100),
        palette_scop(
          palette = heatmap_palette,
          palcolor = heatmap_palcolor
        )
      )
    }
  } else {
    colors <- circlize::colorRamp2(
      seq(limits[1], limits[2], length = 100),
      palette_scop(
        palette = heatmap_palette,
        palcolor = heatmap_palcolor
      )
    )
  }

  lgd <- list()
  lgd[["ht"]] <- ComplexHeatmap::Legend(
    title = exp_name,
    col_fun = colors,
    border = TRUE
  )

  ha_top_list <- list()
  pseudotime <- stats::na.omit(unlist(cell_metadata[, lineages]))
  pseudotime_col <- circlize::colorRamp2(
    breaks = seq(
      min(pseudotime, na.rm = TRUE),
      max(pseudotime, na.rm = TRUE),
      length = 100
    ),
    colors = palette_scop(
      palette = pseudotime_palette,
      palcolor = pseudotime_palcolor
    )
  )

  for (l in lineages) {
    ha_top_list[[l]] <- ComplexHeatmap::HeatmapAnnotation(
      Pseudotime = ComplexHeatmap::anno_simple(
        x = cell_metadata[
          gsub(pattern = l, replacement = "", x = cell_order_list[[l]]),
          l
        ],
        col = pseudotime_col,
        which = ifelse(flip, "row", "column"),
        na_col = "transparent",
        border = TRUE
      ),
      which = ifelse(flip, "row", "column"),
      show_annotation_name = l == lineages[1],
      annotation_name_side = ifelse(flip, "top", "left")
    )
  }

  lgd[["pseudotime"]] <- ComplexHeatmap::Legend(
    title = "Pseudotime",
    col_fun = pseudotime_col,
    border = TRUE
  )

  if (!is.null(cell_annotation)) {
    for (i in seq_along(cell_annotation)) {
      cellan <- cell_annotation[i]
      palette <- cell_annotation_palette[i]
      palcolor <- cell_annotation_palcolor[[i]]
      cell_anno <- cell_metadata[, cellan]
      names(cell_anno) <- rownames(cell_metadata)
      if (!is.numeric(cell_anno)) {
        if (is.logical(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = c(TRUE, FALSE))
        } else if (!is.factor(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = unique(cell_anno))
        }
        for (l in lineages) {
          lineage_cells <- gsub(
            pattern = l,
            replacement = "",
            x = cell_order_list[[l]]
          )
          ha_cell <- list()
          ha_cell[[cellan]] <- ComplexHeatmap::anno_simple(
            x = as.character(cell_anno[lineage_cells]),
            col = palette_scop(
              cell_anno,
              palette = palette,
              palcolor = palcolor
            ),
            which = ifelse(flip, "row", "column"),
            na_col = "transparent",
            border = TRUE
          )
          anno_args <- c(
            ha_cell,
            which = ifelse(flip, "row", "column"),
            show_annotation_name = l == lineages[1],
            annotation_name_side = ifelse(flip, "top", "left")
          )
          anno_args <- c(
            anno_args,
            cell_annotation_params[setdiff(
              names(cell_annotation_params),
              names(anno_args)
            )]
          )
          ha_top <- do.call(ComplexHeatmap::HeatmapAnnotation, args = anno_args)
          if (is.null(ha_top_list[[l]])) {
            ha_top_list[[l]] <- ha_top
          } else {
            ha_top_list[[l]] <- c(ha_top_list[[l]], ha_top)
          }
        }
        lgd[[cellan]] <- ComplexHeatmap::Legend(
          title = cellan,
          labels = levels(cell_anno),
          legend_gp = grid::gpar(
            fill = palette_scop(
              cell_anno,
              palette = palette,
              palcolor = palcolor
            )
          ),
          border = TRUE
        )
      } else {
        col_fun <- circlize::colorRamp2(
          breaks = seq(
            min(cell_anno, na.rm = TRUE),
            max(cell_anno, na.rm = TRUE),
            length = 100
          ),
          colors = palette_scop(palette = palette, palcolor = palcolor)
        )
        for (l in lineages) {
          lineage_cells <- gsub(
            pattern = l,
            replacement = "",
            x = cell_order_list[[l]]
          )
          ha_cell <- list()
          ha_cell[[cellan]] <- ComplexHeatmap::anno_simple(
            x = cell_anno[lineage_cells],
            col = col_fun,
            which = ifelse(flip, "row", "column"),
            na_col = "transparent",
            border = TRUE
          )
          anno_args <- c(
            ha_cell,
            which = ifelse(flip, "row", "column"),
            show_annotation_name = l == lineages[1],
            annotation_name_side = ifelse(flip, "top", "left")
          )
          anno_args <- c(
            anno_args,
            cell_annotation_params[setdiff(
              names(cell_annotation_params),
              names(anno_args)
            )]
          )
          ha_top <- do.call(ComplexHeatmap::HeatmapAnnotation, args = anno_args)
          if (is.null(ha_top_list[[l]])) {
            ha_top_list[[l]] <- ha_top
          } else {
            ha_top_list[[l]] <- c(ha_top_list[[l]], ha_top)
          }
        }
        lgd[[cellan]] <- ComplexHeatmap::Legend(
          title = cellan,
          col_fun = col_fun,
          border = TRUE
        )
      }
    }
  }

  if (!is.null(separate_annotation)) {
    subplots_list <- list()
    for (i in seq_along(separate_annotation)) {
      cellan <- separate_annotation[[i]]
      palette <- separate_annotation_palette[i]
      palcolor <- separate_annotation_palcolor[[i]]
      if (length(cellan) == 1 && cellan %in% colnames(srt@meta.data)) {
        cell_anno <- srt@meta.data[[cellan]]
      } else {
        cell_anno <- numeric()
      }
      if (!is.numeric(cell_anno)) {
        if (is.logical(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = c(TRUE, FALSE))
        } else if (!is.factor(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = unique(cell_anno))
        }
        for (l in lineages) {
          lineage_cells <- gsub(
            pattern = l,
            replacement = "",
            x = cell_order_list[[l]]
          )
          subplots <- CellDensityPlot(
            srt = srt,
            cells = lineage_cells,
            group.by = cellan,
            features = l,
            decreasing = TRUE,
            x_order = "rank",
            palette = palette,
            palcolor = palcolor,
            flip = flip,
            reverse = l %in% lineages[reverse_ht] || l %in% reverse_ht
          ) +
            theme_void()
          subplots_list[[paste0(cellan, ":", l)]] <- subplots
          graphics <- list()
          nm <- paste0(cellan, ":", l)
          funbody <- paste0(
            "
            g <- as_grob(subplots_list[['",
            nm,
            "']] + theme_void() + theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = 'none'));
            g$name <- '",
            nm,
            "';
            grid::grid.draw(g);
            grid::grid.rect(gp = grid::gpar(fill = 'transparent', col = 'black'));
            "
          )
          funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
          eval(
            parse(
              text = paste(
                "block_graphics <- function(index, levels) {",
                funbody,
                "}",
                sep = ""
              )
            ),
            envir = environment()
          )

          ha_cell <- list()
          ha_cell[[paste0(cellan, "\n(separate)")]] <- ComplexHeatmap::anno_block(
            panel_fun = block_graphics,
            which = ifelse(flip, "row", "column"),
            show_name = l == lineages[1]
          )
          anno_args <- c(
            ha_cell,
            which = ifelse(flip, "row", "column"),
            show_annotation_name = l == lineages[1],
            annotation_name_side = ifelse(flip, "top", "left")
          )
          anno_args <- c(
            anno_args,
            separate_annotation_params[setdiff(
              names(separate_annotation_params),
              names(anno_args)
            )]
          )
          ha_top <- do.call(ComplexHeatmap::HeatmapAnnotation, args = anno_args)
          if (is.null(ha_top_list[[l]])) {
            ha_top_list[[l]] <- ha_top
          } else {
            ha_top_list[[l]] <- c(ha_top_list[[l]], ha_top)
          }
        }
        lgd[[paste0("separate:", cellan)]] <- ComplexHeatmap::Legend(
          title = paste0(cellan, "\n(separate)"),
          labels = levels(cell_anno),
          legend_gp = grid::gpar(
            fill = palette_scop(
              cell_anno,
              palette = palette,
              palcolor = palcolor
            )
          ),
          border = TRUE
        )
      } else {
        for (l in lineages) {
          lineage_cells <- gsub(
            pattern = l,
            replacement = "",
            x = cell_order_list[[l]]
          )
          subplots <- DynamicPlot(
            srt = srt,
            cells = lineage_cells,
            lineages = l,
            group.by = NULL,
            features = cellan,
            line_palette = palette,
            line_palcolor = palcolor,
            add_rug = FALSE,
            legend.position = "none",
            compare_features = TRUE,
            x_order = "rank",
            flip = flip,
            reverse = l %in% lineages[reverse_ht] || l %in% reverse_ht
          ) +
            theme_void()
          subplots_list[[paste0(
            paste0(cellan, collapse = ","),
            ":",
            l
          )]] <- subplots
          graphics <- list()
          nm <- paste0(paste0(cellan, collapse = ","), ":", l)
          funbody <- paste0(
            "
            g <- as_grob(subplots_list[['",
            nm,
            "']] + theme_void() + theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = 'none'));
            g$name <- '",
            nm,
            "';
            grid::grid.draw(g);
            grid::grid.rect(gp = grid::gpar(fill = 'transparent', col = 'black'));
            "
          )
          funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
          eval(
            parse(
              text = paste(
                "block_graphics <- function(index, levels) {",
                funbody,
                "}",
                sep = ""
              )
            ),
            envir = environment()
          )

          ha_cell <- list()
          ha_cell[[paste0(
            paste0(cellan, collapse = ","),
            "\n(separate)"
          )]] <- ComplexHeatmap::anno_block(
            panel_fun = block_graphics,
            which = ifelse(flip, "row", "column"),
            show_name = l == lineages[1]
          )
          anno_args <- c(
            ha_cell,
            which = ifelse(flip, "row", "column"),
            show_annotation_name = l == lineages[1],
            annotation_name_side = ifelse(flip, "top", "left")
          )
          anno_args <- c(
            anno_args,
            separate_annotation_params[setdiff(
              names(separate_annotation_params),
              names(anno_args)
            )]
          )
          ha_top <- do.call(ComplexHeatmap::HeatmapAnnotation, args = anno_args)
          if (is.null(ha_top_list[[l]])) {
            ha_top_list[[l]] <- ha_top
          } else {
            ha_top_list[[l]] <- c(ha_top_list[[l]], ha_top)
          }
        }
        lgd[[paste0("separate:", paste0(cellan, collapse = ","))]] <- ComplexHeatmap::Legend(
          title = "Features\n(separate)",
          labels = cellan,
          legend_gp = grid::gpar(
            fill = palette_scop(cellan, palette = palette, palcolor = palcolor)
          ),
          border = TRUE
        )
      }
    }
  }

  if (is.null(feature_split)) {
    if (is.null(n_split) || isTRUE(nrow(mat) <= n_split)) {
      row_split_raw <- row_split <- feature_split <- NULL
    } else {
      if (n_split == 1) {
        row_split_raw <- row_split <- feature_split <- stats::setNames(
          rep(1, nrow(mat_split)),
          rownames(mat_split)
        )
      } else {
        if (split_method == "mfuzz") {
          status <- tryCatch(check_r("e1071"), error = identity)
          if (inherits(status, "error")) {
            log_message(
              "The e1071 package was not found. Switch split_method to 'kmeans'",
              message_type = "warning"
            )
            split_method <- "kmeans"
          } else {
            mat_split_tmp <- mat_split
            colnames(mat_split_tmp) <- make.unique(colnames(mat_split_tmp))
            mat_split_tmp <- standardise(mat_split_tmp)
            min_fuzzification <- mestimate(mat_split_tmp)
            if (is.null(fuzzification)) {
              fuzzification <- min_fuzzification + 0.1
            } else {
              if (fuzzification <= min_fuzzification) {
                log_message(
                  "fuzzification value is samller than estimated:",
                  round(min_fuzzification, 2),
                  message_type = "warning"
                )
              }
            }
            cl <- e1071::cmeans(
              mat_split_tmp,
              centers = n_split,
              method = "cmeans",
              m = fuzzification
            )
            if (length(cl$cluster) == 0) {
              log_message(
                "Clustering with mfuzz failed (fuzzification=",
                round(fuzzification, 2),
                "). Please set a larger fuzzification parameter manually.",
                message_type = "error"
              )
            }
            row_split <- feature_split <- cl$cluster
          }
        }
        if (split_method == "kmeans") {
          km <- stats::kmeans(
            mat_split,
            centers = n_split,
            iter.max = 1e4,
            nstart = 20
          )
          row_split <- feature_split <- km$cluster
        }
        if (split_method == "kmeans-peaktime") {
          feature_y <- feature_metadata[rownames(mat_split), order_by]
          names(feature_y) <- rownames(mat_split)
          km <- stats::kmeans(
            feature_y,
            centers = n_split,
            iter.max = 1e4,
            nstart = 20
          )
          row_split <- feature_split <- km$cluster
        }
        if (split_method == "hclust") {
          hc <- stats::hclust(
            stats::as.dist(
              proxyC::dist(mat_split)
            )
          )
          row_split <- feature_split <- stats::cutree(hc, k = n_split)
        }
        if (split_method == "hclust-peaktime") {
          feature_y <- feature_metadata[rownames(mat_split), order_by]
          names(feature_y) <- rownames(mat_split)
          hc <- stats::hclust(
            stats::as.dist(
              proxyC::dist(feature_y)
            )
          )
          row_split <- feature_split <- stats::cutree(hc, k = n_split)
        }
      }
      df <- data.frame(
        row_split = row_split,
        order_by = feature_metadata[names(row_split), order_by]
      )
      df_order <- stats::aggregate(df, by = list(row_split), FUN = mean)
      df_order <- df_order[
        order(df_order[["order_by"]], decreasing = decreasing), ,
        drop = FALSE
      ]
      if (!is.null(split_order)) {
        df_order <- df_order[split_order, , drop = FALSE]
      }
      split_levels <- c()
      for (i in seq_len(nrow(df_order))) {
        raw_nm <- df_order[i, "row_split"]
        feature_split[feature_split == raw_nm] <- paste0("C", i)
        level <- paste0("C", i, "(", sum(row_split == raw_nm), ")")
        row_split[row_split == raw_nm] <- level
        split_levels <- c(split_levels, level)
      }
      row_split_raw <- row_split <- factor(row_split, levels = split_levels)
      feature_split <- factor(
        feature_split,
        levels = paste0("C", seq_len(nrow(df_order)))
      )
    }
  } else {
    row_split_raw <- row_split <- feature_split <- feature_split[row.names(mat)]
  }
  if (!is.null(feature_split)) {
    feature_metadata[["feature_split"]] <- feature_split
  } else {
    feature_metadata[["feature_split"]] <- NA
  }

  ha_left <- NULL
  if (!is.null(row_split)) {
    if (isTRUE(cluster_row_slices)) {
      if (!isTRUE(cluster_rows)) {
        dend <- ComplexHeatmap::cluster_within_group(
          Matrix::t(mat_split),
          row_split_raw
        )
        cluster_rows <- dend
        row_split <- length(unique(row_split_raw))
      }
    }
    funbody <- paste0(
      "
      grid::grid.rect(gp = grid::gpar(fill = palette_scop(",
      paste0("c('", paste0(levels(row_split_raw), collapse = "','"), "')"),
      ",palette = '",
      feature_split_palette,
      "',palcolor=c(",
      paste0(
        "'",
        paste0(unlist(feature_split_palcolor), collapse = "','"),
        "'"
      ),
      "))[nm]))
    "
    )
    funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
    eval(
      parse(
        text = paste(
          "panel_fun <- function(index, nm) {",
          funbody,
          "}",
          sep = ""
        )
      ),
      envir = environment()
    )
    ha_clusters <- ComplexHeatmap::HeatmapAnnotation(
      features_split = ComplexHeatmap::anno_block(
        align_to = split(seq_along(row_split_raw), row_split_raw),
        panel_fun = methods::getFunction("panel_fun", where = environment()),
        width = grid::unit(0.1, "in"),
        height = grid::unit(0.1, "in"),
        show_name = FALSE,
        which = ifelse(flip, "column", "row")
      ),
      which = ifelse(flip, "column", "row"),
      border = TRUE
    )
    if (is.null(ha_left)) {
      ha_left <- ha_clusters
    } else {
      ha_left <- c(ha_left, ha_clusters)
    }
    lgd[["Cluster"]] <- ComplexHeatmap::Legend(
      title = "Cluster",
      labels = intersect(levels(row_split_raw), row_split_raw),
      legend_gp = grid::gpar(
        fill = palette_scop(
          intersect(levels(row_split_raw), row_split_raw),
          type = "discrete",
          palette = feature_split_palette,
          palcolor = feature_split_palcolor,
          matched = TRUE
        )
      ),
      border = TRUE
    )
  }

  if (isTRUE(cluster_rows) && !is.null(cluster_features_by)) {
    mat_cluster <- mat[,
      unlist(cell_order_list[cluster_features_by]),
      drop = FALSE
    ]
    if (is.null(row_split)) {
      dend <- stats::as.dendrogram(
        stats::hclust(
          stats::as.dist(
            proxyC::dist(mat_cluster)
          )
        )
      )
      dend_ordered <- stats::reorder(
        dend,
        wts = colMeans(mat_cluster),
        agglo.FUN = mean
      )
      cluster_rows <- dend_ordered
    } else {
      row_split <- length(unique(row_split_raw))
      dend <- cluster_within_group2(
        Matrix::t(mat_cluster),
        row_split_raw
      )
      cluster_rows <- dend
    }
  }

  l <- lineages[1]
  ht_args <- list(
    matrix = mat[, cell_order_list[[l]], drop = FALSE],
    col = colors,
    row_split = row_split,
    cluster_rows = cluster_rows,
    cluster_row_slices = cluster_row_slices,
    column_split = column_split,
    cluster_columns = cluster_columns,
    cluster_column_slices = cluster_column_slices,
    use_raster = TRUE
  )
  ht_args <- c(ht_args, ht_params[setdiff(names(ht_params), names(ht_args))])
  ht_list <- do.call(ComplexHeatmap::Heatmap, args = ht_args)
  features_ordered <- rownames(mat)[unlist(
    suppressWarnings(
      ComplexHeatmap::row_order(
        ht_list
      )
    )
  )]
  feature_metadata[["index"]] <- stats::setNames(
    object = seq_along(features_ordered),
    nm = features_ordered
  )[rownames(feature_metadata)]

  if (is.null(features_label)) {
    if (nlabel > 0) {
      if (length(features) > nlabel) {
        index_from <- ceiling((length(features_ordered) / nlabel) / 2)
        index_to <- length(features_ordered)
        index <- unique(
          round(
            seq(
              from = index_from,
              to = index_to,
              length.out = nlabel
            )
          )
        )
      } else {
        index <- seq_along(features_ordered)
      }
    } else {
      index <- NULL
    }
  } else {
    index <- which(features_ordered %in% features_label)
    drop <- setdiff(features_label, features_ordered)
    if (length(drop) > 0) {
      log_message(
        paste0(paste0(drop, collapse = ","), "was not found in the features"),
        message_type = "warning"
      )
    }
  }
  if (length(index) > 0) {
    ha_mark <- ComplexHeatmap::HeatmapAnnotation(
      gene = ComplexHeatmap::anno_mark(
        at = which(rownames(feature_metadata) %in% features_ordered[index]),
        labels = feature_metadata[
          which(rownames(feature_metadata) %in% features_ordered[index]),
          "features"
        ],
        side = ifelse(flip, "top", "left"),
        labels_gp = grid::gpar(fontsize = label_size, col = label_color),
        link_gp = grid::gpar(fontsize = label_size, col = label_color),
        which = ifelse(flip, "column", "row")
      ),
      which = ifelse(flip, "column", "row"),
      show_annotation_name = FALSE
    )
    if (is.null(ha_left)) {
      ha_left <- ha_mark
    } else {
      ha_left <- c(ha_mark, ha_left)
    }
  }

  ha_right <- NULL
  if (length(lineages) > 1) {
    ha_list <- list()
    for (l in lineages) {
      ha_list[[l]] <- ComplexHeatmap::anno_simple(
        x = is.na(feature_metadata[, paste0(l, order_by)]) + 0,
        col = c("0" = "#181830", "1" = "transparent"),
        width = grid::unit(5, "mm"),
        height = grid::unit(5, "mm"),
        which = ifelse(flip, "column", "row")
      )
    }
    ha_lineage <- do.call(
      ComplexHeatmap::HeatmapAnnotation,
      args = c(
        ha_list,
        which = ifelse(flip, "column", "row"),
        annotation_name_side = ifelse(flip, "left", "top"),
        border = TRUE
      )
    )
    if (is.null(ha_right)) {
      ha_right <- ha_lineage
    } else {
      ha_right <- c(ha_right, ha_lineage)
    }
  }
  if (!is.null(feature_annotation)) {
    for (i in seq_along(feature_annotation)) {
      featan <- feature_annotation[i]
      palette <- feature_annotation_palette[i]
      palcolor <- feature_annotation_palcolor[[i]]
      featan_values <- feature_metadata[, featan]
      if (!is.numeric(featan_values)) {
        if (is.logical(featan_values)) {
          featan_values <- factor(featan_values, levels = c(TRUE, FALSE))
        } else if (!is.factor(featan_values)) {
          featan_values <- factor(featan_values, levels = unique(featan_values))
        }
        ha_feature <- list()
        ha_feature[[featan]] <- ComplexHeatmap::anno_simple(
          x = as.character(featan_values),
          col = palette_scop(
            featan_values,
            palette = palette,
            palcolor = palcolor
          ),
          which = ifelse(flip, "column", "row"),
          na_col = "transparent",
          border = TRUE
        )
        anno_args <- c(
          ha_feature,
          which = ifelse(flip, "column", "row"),
          show_annotation_name = TRUE,
          annotation_name_side = ifelse(flip, "left", "top"),
          border = TRUE
        )
        anno_args <- c(
          anno_args,
          feature_annotation_params[setdiff(
            names(feature_annotation_params),
            names(anno_args)
          )]
        )
        ha_feature <- do.call(
          ComplexHeatmap::HeatmapAnnotation,
          args = anno_args
        )
        if (is.null(ha_right)) {
          ha_right <- ha_feature
        } else {
          ha_right <- c(ha_right, ha_feature)
        }
        lgd[[featan]] <- ComplexHeatmap::Legend(
          title = featan,
          labels = levels(featan_values),
          legend_gp = grid::gpar(
            fill = palette_scop(
              featan_values,
              palette = palette,
              palcolor = palcolor
            )
          ),
          border = TRUE
        )
      } else {
        col_fun <- circlize::colorRamp2(
          breaks = seq(
            min(featan_values, na.rm = TRUE),
            max(featan_values, na.rm = TRUE),
            length = 100
          ),
          colors = palette_scop(
            palette = palette, palcolor = palcolor
          )
        )
        ha_feature <- list()
        ha_feature[[featan]] <- ComplexHeatmap::anno_simple(
          x = featan_values,
          col = col_fun,
          which = ifelse(flip, "column", "row"),
          na_col = "transparent",
          border = TRUE
        )
        anno_args <- c(
          ha_feature,
          which = ifelse(flip, "column", "row"),
          show_annotation_name = TRUE,
          annotation_name_side = ifelse(flip, "left", "top"),
          border = TRUE
        )
        anno_args <- c(
          anno_args,
          feature_annotation_params[setdiff(
            names(feature_annotation_params),
            names(anno_args)
          )]
        )
        ha_feature <- do.call(ComplexHeatmap::HeatmapAnnotation, args = anno_args)
        if (is.null(ha_right)) {
          ha_right <- ha_feature
        } else {
          ha_right <- c(ha_right, ha_feature)
        }
        lgd[[featan]] <- ComplexHeatmap::Legend(
          title = featan,
          col_fun = col_fun,
          border = TRUE
        )
      }
    }
  }

  enrichment <- heatmap_enrichment(
    geneID = feature_metadata[["features"]],
    geneID_groups = feature_metadata[["feature_split"]],
    feature_split_palette = feature_split_palette,
    feature_split_palcolor = feature_split_palcolor,
    ha_right = ha_right,
    flip = flip,
    anno_terms = anno_terms,
    anno_keys = anno_keys,
    anno_features = anno_features,
    terms_width = terms_width,
    terms_fontsize = terms_fontsize,
    keys_width = keys_width,
    keys_fontsize = keys_fontsize,
    features_width = features_width,
    features_fontsize = features_fontsize,
    IDtype = IDtype,
    species = species,
    db_update = db_update,
    db_version = db_version,
    db_combine = db_combine,
    convert_species = convert_species,
    Ensembl_version = Ensembl_version,
    mirror = mirror,
    db = db,
    TERM2GENE = TERM2GENE,
    TERM2NAME = TERM2NAME,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    GO_simplify = GO_simplify,
    GO_simplify_cutoff = GO_simplify_cutoff,
    simplify_method = simplify_method,
    simplify_similarityCutoff = simplify_similarityCutoff,
    pvalueCutoff = pvalueCutoff,
    padjustCutoff = padjustCutoff,
    topTerm = topTerm,
    show_termid = show_termid,
    topWord = topWord,
    words_excluded = words_excluded
  )
  res <- enrichment$res
  ha_right <- enrichment$ha_right

  ht_list <- NULL
  for (l in lineages) {
    if (l == lineages[1]) {
      left_annotation <- ha_left
    } else {
      left_annotation <- NULL
    }
    if (l == lineages[length(lineages)]) {
      right_annotation <- ha_right
    } else {
      right_annotation <- NULL
    }
    ht_args <- list(
      name = l,
      matrix = if (flip) {
        Matrix::t(mat[, cell_order_list[[l]], drop = FALSE])
      } else {
        mat[, cell_order_list[[l]], drop = FALSE]
      },
      col = colors,
      row_title = row_title %||% if (flip) l else character(0),
      row_title_side = row_title_side,
      column_title = column_title %||% if (flip) character(0) else l,
      column_title_side = if (flip) "top" else column_title_side,
      row_title_rot = row_title_rot,
      column_title_rot = column_title_rot,
      row_split = if (flip) column_split else row_split,
      column_split = if (flip) row_split else column_split,
      cluster_rows = if (flip) cluster_columns else cluster_rows,
      cluster_columns = if (flip) cluster_rows else cluster_columns,
      cluster_row_slices = if (flip) {
        cluster_column_slices
      } else {
        cluster_row_slices
      },
      cluster_column_slices = if (flip) {
        cluster_row_slices
      } else {
        cluster_column_slices
      },
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      row_names_side = row_names_side,
      column_names_side = column_names_side,
      row_names_rot = row_names_rot,
      column_names_rot = column_names_rot,
      top_annotation = if (flip) left_annotation else ha_top_list[[l]],
      left_annotation = if (flip) ha_top_list[[l]] else left_annotation,
      bottom_annotation = if (flip) right_annotation else NULL,
      right_annotation = if (flip) NULL else right_annotation,
      show_heatmap_legend = FALSE,
      border = border,
      use_raster = use_raster,
      raster_device = raster_device,
      raster_by_magick = raster_by_magick,
      width = if (is.numeric(width[l])) grid::unit(width[l], units = units) else NULL,
      height = if (is.numeric(height[l])) {
        grid::unit(height[l], units = units)
      } else {
        NULL
      }
    )
    if (any(names(ht_params) %in% names(ht_args))) {
      log_message(
        "ht_params: ",
        paste0(intersect(names(ht_params), names(ht_args)), collapse = ","),
        " were duplicated and will not be used.",
        message_type = "warning"
      )
    }
    ht_args <- c(ht_args, ht_params[setdiff(names(ht_params), names(ht_args))])
    if (isTRUE(flip)) {
      if (is.null(ht_list)) {
        ht_list <- do.call(ComplexHeatmap::Heatmap, args = ht_args)
      } else {
        ht_list <- ht_list %v% do.call(ComplexHeatmap::Heatmap, args = ht_args)
      }
    } else {
      ht_list <- ht_list + do.call(ComplexHeatmap::Heatmap, args = ht_args)
    }
  }

  if (
    (!is.null(row_split) && length(index) > 0) ||
      any(c(anno_terms, anno_keys, anno_features)) ||
      !is.null(width) ||
      !is.null(height)
  ) {
    fix <- TRUE
    if (is.null(width) || is.null(height)) {
      log_message(
        "\nThe size of the heatmap is fixed because certain elements are not scalable.\nThe width and height of the heatmap are determined by the size of the current viewport.\nIf you want to have more control over the size, you can manually set the parameters 'width' and 'height'.\n"
      )
    }
  } else {
    fix <- FALSE
  }
  rendersize <- heatmap_rendersize(
    width = width,
    height = height,
    units = units,
    ha_top_list = ha_top_list,
    ha_left = ha_left,
    ha_right = ha_right,
    ht_list = ht_list,
    legend_list = lgd,
    flip = flip
  )
  width_sum <- rendersize[["width_sum"]]
  height_sum <- rendersize[["height_sum"]]

  if (isTRUE(fix)) {
    fixsize <- heatmap_fixsize(
      width = width,
      width_sum = width_sum,
      height = height,
      height_sum = height_sum,
      units = units,
      ht_list = ht_list,
      legend_list = lgd
    )
    ht_width <- fixsize[["ht_width"]]
    ht_height <- fixsize[["ht_height"]]

    g_tree <- grid::grid.grabExpr(
      {
        ComplexHeatmap::draw(ht_list, annotation_legend_list = lgd)
        for (enrich in db) {
          enrich_anno <- names(ha_right)[grep(
            paste0("_split_", enrich),
            names(ha_right)
          )]
          if (length(enrich_anno) > 0) {
            for (enrich_anno_element in enrich_anno) {
              enrich_obj <- strsplit(enrich_anno_element, "_split_")[[1]][1]
              ComplexHeatmap::decorate_annotation(
                enrich_anno_element,
                slice = 1,
                {
                  grid::grid.text(
                    paste0(enrich, " (", enrich_obj, ")"),
                    x = grid::unit(1, "npc"),
                    y = grid::unit(1, "npc") + grid::unit(2.5, "mm"),
                    just = c("left", "bottom")
                  )
                }
              )
            }
          }
        }
        if (is.numeric(pseudotime_label)) {
          for (n in seq_along(pseudotime_label)) {
            pse <- pseudotime_label[n]
            col <- pseudotime_label_color[n]
            lty <- pseudotime_label_linetype[n]
            lwd <- pseudotime_label_linewidth[n]
            for (l in lineages) {
              for (slice in 1:max(nlevels(row_split), 1)) {
                ComplexHeatmap::decorate_heatmap_body(
                  l,
                  {
                    pseudotime <- cell_metadata[
                      gsub(
                        pattern = l,
                        replacement = "",
                        x = cell_order_list[[l]]
                      ),
                      l
                    ]
                    i <- which.min(abs(pseudotime - pse))
                    if (flip) {
                      x <- 1 - (i / length(pseudotime))
                      grid::grid.lines(
                        c(0, 1),
                        c(x, x),
                        gp = grid::gpar(lty = lty, lwd = lwd, col = col)
                      )
                    } else {
                      i <- which.min(abs(pseudotime - pse))
                      x <- i / length(pseudotime)
                      grid::grid.lines(
                        c(x, x),
                        c(0, 1),
                        gp = grid::gpar(lty = lty, lwd = lwd, col = col)
                      )
                    }
                  },
                  row_slice = ifelse(flip, 1, slice),
                  column_slice = ifelse(flip, slice, 1)
                )
              }
            }
          }
        }
      },
      width = ht_width,
      height = ht_height,
      wrap = TRUE,
      wrap.grobs = TRUE
    )
  } else {
    ht_width <- grid::unit(width_sum, units = units)
    ht_height <- grid::unit(height_sum, units = units)
    g_tree <- grid::grid.grabExpr(
      {
        ComplexHeatmap::draw(ht_list, annotation_legend_list = lgd)
        for (enrich in db) {
          enrich_anno <- names(ha_right)[grep(
            paste0("_split_", enrich),
            names(ha_right)
          )]
          if (length(enrich_anno) > 0) {
            for (enrich_anno_element in enrich_anno) {
              enrich_obj <- strsplit(enrich_anno_element, "_split_")[[1]][1]
              ComplexHeatmap::decorate_annotation(
                enrich_anno_element,
                slice = 1,
                {
                  grid::grid.text(
                    paste0(enrich, " (", enrich_obj, ")"),
                    x = grid::unit(1, "npc"),
                    y = grid::unit(1, "npc") + grid::unit(2.5, "mm"),
                    just = c("left", "bottom")
                  )
                }
              )
            }
          }
        }
        if (is.numeric(pseudotime_label)) {
          for (n in seq_along(pseudotime_label)) {
            pse <- pseudotime_label[n]
            col <- pseudotime_label_color[n]
            lty <- pseudotime_label_linetype[n]
            lwd <- pseudotime_label_linewidth[n]
            for (l in lineages) {
              for (slice in 1:max(nlevels(row_split), 1)) {
                ComplexHeatmap::decorate_heatmap_body(
                  l,
                  {
                    pseudotime <- cell_metadata[
                      gsub(
                        pattern = l,
                        replacement = "",
                        x = cell_order_list[[l]]
                      ),
                      l
                    ]
                    i <- which.min(abs(pseudotime - pse))
                    if (flip) {
                      x <- 1 - (i / length(pseudotime))
                      grid::grid.lines(
                        c(0, 1),
                        c(x, x),
                        gp = grid::gpar(lty = lty, lwd = lwd, col = col)
                      )
                    } else {
                      i <- which.min(abs(pseudotime - pse))
                      x <- i / length(pseudotime)
                      grid::grid.lines(
                        c(x, x),
                        c(0, 1),
                        gp = grid::gpar(lty = lty, lwd = lwd, col = col)
                      )
                    }
                  },
                  row_slice = ifelse(flip, 1, slice),
                  column_slice = ifelse(flip, slice, 1)
                )
              }
            }
          }
        }
      },
      width = ht_width,
      height = ht_height,
      wrap = TRUE,
      wrap.grobs = TRUE
    )
  }

  if (isTRUE(fix)) {
    p <- panel_fix_overall(
      g_tree,
      width = as.numeric(ht_width),
      height = as.numeric(ht_height),
      units = units
    )
  } else {
    p <- patchwork::wrap_plots(g_tree)
  }

  return(
    list(
      plot = p,
      matrix = mat,
      cell_order = cell_order_list,
      feature_split = feature_split,
      cell_metadata = cell_metadata,
      feature_metadata = feature_metadata,
      enrichment = res
    )
  )
}
