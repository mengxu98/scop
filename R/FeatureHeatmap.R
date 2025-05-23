#' FeatureHeatmap
#'
#' @inheritParams GroupHeatmap
#' @param max_cells An integer, maximum number of cells to sample per group, default is 100.
#' @param cell_order A vector of cell names defining the order of cells, default is NULL.
#'
#' @seealso \code{\link{RunDEtest}}
#'
#' @importFrom ComplexHeatmap %v%
#' @export
#'
#' @examples
#' library(dplyr)
#' data("pancreas_sub")
#' pancreas_sub <- RunDEtest(
#'   pancreas_sub,
#'   group_by = "CellType"
#' )
#' de_filter <- filter(
#'   pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox,
#'   p_val_adj < 0.05 & avg_log2FC > 1
#' )
#' ht1 <- FeatureHeatmap(
#'   srt = pancreas_sub,
#'   features = de_filter$gene,
#'   group.by = "CellType",
#'   split.by = "Phase",
#'   cell_split_palette = "Dark2"
#' )
#' ht1$plot
#' panel_fix(
#'   ht1$plot,
#'   height = 4,
#'   width = 6,
#'   raster = TRUE,
#'   dpi = 50
#' )
#'
#' ht2 <- FeatureHeatmap(
#'   srt = pancreas_sub,
#'   features = de_filter$gene,
#'   group.by = c("CellType", "SubCellType"),
#'   n_split = 4,
#'   cluster_rows = TRUE,
#'   cluster_row_slices = TRUE,
#'   cluster_columns = TRUE,
#'   cluster_column_slices = TRUE,
#'   ht_params = list(row_gap = grid::unit(0, "mm"))
#' )
#' ht2$plot
#'
#' ht3 <- FeatureHeatmap(
#'   srt = pancreas_sub,
#'   features = de_filter$gene,
#'   feature_split = de_filter$group1,
#'   group.by = "CellType",
#'   species = "Mus_musculus",
#'   db = "GO_BP",
#'   anno_terms = TRUE,
#'   anno_keys = TRUE,
#'   anno_features = TRUE
#' )
#' ht3$plot
#'
#' pancreas_sub <- AnnotateFeatures(
#'   pancreas_sub,
#'   species = "Mus_musculus",
#'   db = c("TF", "CSPA")
#' )
#' ht4 <- FeatureHeatmap(
#'   srt = pancreas_sub,
#'   features = de_filter$gene,
#'   n_split = 4,
#'   group.by = "CellType",
#'   heatmap_palette = "viridis",
#'   feature_annotation = c("TF", "CSPA"),
#'   feature_annotation_palcolor = list(
#'     c("gold", "steelblue"), c("forestgreen")
#'   ),
#'   cell_annotation = c("Phase", "G2M_score"),
#'   cell_annotation_palette = c("Dark2", "Purples")
#' )
#' ht4$plot
#'
#' ht5 <- FeatureHeatmap(
#'   srt = pancreas_sub,
#'   features = de_filter$gene,
#'   n_split = 4,
#'   group.by = "CellType",
#'   heatmap_palette = "viridis",
#'   feature_annotation = c("TF", "CSPA"),
#'   feature_annotation_palcolor = list(
#'     c("gold", "steelblue"), c("forestgreen")
#'   ),
#'   cell_annotation = c("Phase", "G2M_score"),
#'   cell_annotation_palette = c("Dark2", "Purples"),
#'   flip = TRUE,
#'   column_title_rot = 45
#' )
#' ht5$plot
#'
#' pancreas_sub <- RunSlingshot(
#'   srt = pancreas_sub,
#'   group_by = "SubCellType",
#'   reduction = "UMAP"
#' )
#' ht6 <- FeatureHeatmap(
#'   srt = pancreas_sub,
#'   features = de_filter$gene,
#'   nlabel = 10,
#'   cell_order = names(sort(pancreas_sub$Lineage1)),
#'   cell_annotation = c("SubCellType", "Lineage1"),
#'   cell_annotation_palette = c("Paired", "cividis")
#' )
#' ht6$plot
FeatureHeatmap <- function(
    srt,
    features = NULL,
    cells = NULL,
    group.by = NULL,
    split.by = NULL,
    within_groups = FALSE,
    max_cells = 100,
    cell_order = NULL,
    border = TRUE,
    flip = FALSE,
    layer = "counts",
    assay = NULL,
    exp_method = c("zscore", "raw", "fc", "log2fc", "log1p"),
    exp_legend_title = NULL,
    limits = NULL,
    lib_normalize = identical(layer, "counts"),
    libsize = NULL,
    feature_split = NULL,
    feature_split_by = NULL,
    n_split = NULL,
    split_order = NULL,
    split_method = c("kmeans", "hclust", "mfuzz"),
    decreasing = FALSE,
    fuzzification = NULL,
    cluster_features_by = NULL,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cluster_row_slices = FALSE,
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
    heatmap_palette = "RdBu",
    heatmap_palcolor = NULL,
    group_palette = "Paired",
    group_palcolor = NULL,
    cell_split_palette = "simspec",
    cell_split_palcolor = NULL,
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
  data_nm <- c(ifelse(isTRUE(lib_normalize), "normalized", ""), layer)
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

  assay <- assay %||% DefaultAssay(srt)
  if (length(feature_split) != 0 && length(feature_split) != length(features)) {
    stop("feature_split must be the same length as features")
  }
  if (is.null(group.by)) {
    srt@meta.data[["All.groups"]] <- factor("")
    group.by <- "All.groups"
  }
  if (any(!group.by %in% colnames(srt@meta.data))) {
    stop(
      group.by[!group.by %in% colnames(srt@meta.data)],
      " is not in the meta data of the Seurat object."
    )
  }
  if (!is.null(group.by)) {
    for (g in group.by) {
      if (!is.factor(srt@meta.data[[g]])) {
        srt@meta.data[[g]] <- factor(
          srt@meta.data[[g]],
          levels = unique(srt@meta.data[[g]])
        )
      }
    }
  }
  if (length(split.by) > 1) {
    stop("'split.by' only support one variable.")
  }
  if (any(!split.by %in% colnames(srt@meta.data))) {
    stop(
      split.by[!split.by %in% colnames(srt@meta.data)],
      " is not in the meta data of the Seurat object."
    )
  }
  if (!is.null(split.by)) {
    if (!is.factor(srt@meta.data[[split.by]])) {
      srt@meta.data[[split.by]] <- factor(
        srt@meta.data[[split.by]],
        levels = unique(srt@meta.data[[split.by]])
      )
    }
  }
  if (length(group_palette) == 1) {
    group_palette <- rep(group_palette, length(group.by))
  }
  if (length(group_palette) != length(group.by)) {
    stop("'group_palette' must be the same length as 'group.by'")
  }
  group_palette <- stats::setNames(group_palette, nm = group.by)
  raw.group.by <- group.by
  raw.group_palette <- group_palette
  if (isTRUE(within_groups)) {
    new.group.by <- c()
    new.group_palette <- group_palette
    for (g in group.by) {
      groups <- split(colnames(srt), srt[[g, drop = TRUE]])
      new.group_palette[g] <- list(rep(new.group_palette[g], length(groups)))
      for (nm in names(groups)) {
        srt[[make.names(nm)]] <- factor(
          NA,
          levels = levels(srt[[g, drop = TRUE]])
        )
        srt[[make.names(nm)]][colnames(srt) %in% groups[[nm]], ] <- nm
        new.group.by <- c(new.group.by, make.names(nm))
      }
    }
    group.by <- new.group.by
    group_palette <- unlist(new.group_palette)
  }
  if (is.null(feature_split_by)) {
    feature_split_by <- group.by
  }
  if (any(!feature_split_by %in% group.by)) {
    stop("feature_split_by must be a subset of group.by")
  }
  if (!is.null(feature_split) && !is.factor(feature_split)) {
    feature_split <- factor(feature_split, levels = unique(feature_split))
  }
  if (length(feature_split) != 0 && length(feature_split) != length(features)) {
    stop("feature_split must be the same length as features")
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
      stop(
        "cell_annotation_palette and cell_annotation_palcolor must be the same length as cell_annotation"
      )
    }
    if (
      any(
        !cell_annotation %in%
          c(colnames(srt@meta.data), rownames(Seurat::GetAssay(
            srt,
            assay = assay
          )))
      )
    ) {
      stop(
        "cell_annotation: ",
        paste0(
          cell_annotation[
            !cell_annotation %in%
              c(colnames(srt@meta.data), rownames(Seurat::GetAssay(
                srt,
                assay = assay
              )))
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
      stop(
        "feature_annotation_palette and feature_annotation_palcolor must be the same length as feature_annotation"
      )
    }
    if (
      any(!feature_annotation %in% colnames(GetFeaturesData(srt, assay = assay)))
    ) {
      stop(
        "feature_annotation: ",
        paste0(
          feature_annotation[
            !feature_annotation %in% colnames(GetFeaturesData(srt, assay = assay))
          ],
          collapse = ","
        ),
        " is not in the meta data of the ",
        assay,
        " assay in the Seurat object."
      )
    }
  }
  if (length(width) == 1) {
    width <- rep(width, length(group.by))
  }
  if (length(height) == 1) {
    height <- rep(height, length(group.by))
  }
  if (length(width) >= 1) {
    names(width) <- group.by
  }
  if (length(height) >= 1) {
    names(height) <- group.by
  }

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

  if (is.null(cells)) {
    cells <- colnames(srt@assays[[1]])
  }
  if (all(!cells %in% colnames(srt@assays[[1]]))) {
    stop("No cells found.")
  }
  if (!all(cells %in% colnames(srt@assays[[1]]))) {
    warning("Some cells not found.", immediate. = TRUE)
  }
  cells <- intersect(cells, colnames(srt@assays[[1]]))

  if (is.null(features)) {
    features <- VariableFeatures(srt, assay = assay)
  }
  index <- features %in%
    c(rownames(Seurat::GetAssay(
      srt,
      assay = assay
    )), colnames(srt@meta.data))
  features <- features[index]
  features_unique <- make.unique(features)
  if (!is.null(feature_split)) {
    feature_split <- feature_split[index]
    names(feature_split) <- features_unique
  }

  cell_groups <- list()
  for (cell_group in group.by) {
    if (!is.factor(srt@meta.data[[cell_group]])) {
      srt@meta.data[[cell_group]] <- factor(
        srt@meta.data[[cell_group]],
        levels = unique(srt@meta.data[[cell_group]])
      )
    }
    if (is.null(split.by)) {
      cell_groups[[cell_group]] <- unlist(
        lapply(levels(srt@meta.data[[cell_group]]), function(x) {
          cells_sub <- colnames(srt@assays[[1]])[which(
            srt@meta.data[[cell_group]] == x
          )]
          cells_sub <- intersect(cells, cells_sub)
          size <- ifelse(
            length(cells_sub) > max_cells,
            max_cells,
            length(cells_sub)
          )
          cells_sample <- sample(cells_sub, size)
          out <- stats::setNames(rep(x, size), cells_sample)
          return(out)
        }),
        use.names = TRUE
      )
      levels <- levels(srt@meta.data[[cell_group]])
      cell_groups[[cell_group]] <- factor(
        cell_groups[[cell_group]],
        levels = levels[levels %in% cell_groups[[cell_group]]]
      )
    } else {
      if (!is.factor(srt@meta.data[[split.by]])) {
        srt@meta.data[[split.by]] <- factor(
          srt@meta.data[[split.by]],
          levels = unique(srt@meta.data[[split.by]])
        )
      }
      cell_groups[[cell_group]] <- unlist(
        lapply(levels(srt@meta.data[[cell_group]]), function(x) {
          cells_sub <- colnames(srt@assays[[1]])[
            srt@meta.data[[cell_group]] == x
          ]
          cells_sub <- intersect(cells, cells_sub)
          cells_tmp <- NULL
          for (sp in levels(srt@meta.data[[split.by]])) {
            cells_sp <- cells_sub[srt@meta.data[cells_sub, split.by] == sp]
            size <- ifelse(
              length(cells_sp) > max_cells,
              max_cells,
              length(cells_sp)
            )
            cells_sample <- sample(cells_sp, size)
            cells_tmp <- c(
              cells_tmp,
              stats::setNames(rep(paste0(x, " : ", sp), size), cells_sample)
            )
          }
          size <- ifelse(
            length(cells_tmp) > max_cells,
            max_cells,
            length(cells_tmp)
          )
          out <- sample(cells_tmp, size)
          return(out)
        }),
        use.names = TRUE
      )
      levels <- apply(
        expand.grid(
          levels(srt@meta.data[[split.by]]),
          levels(srt@meta.data[[cell_group]])
        ),
        1,
        function(x) paste0(x[2:1], collapse = " : ")
      )
      cell_groups[[cell_group]] <- factor(
        cell_groups[[cell_group]],
        levels = levels[levels %in% cell_groups[[cell_group]]]
      )
    }
    if (!is.null(cell_order)) {
      cell_groups[[cell_group]] <- cell_groups[[cell_group]][intersect(
        cell_order,
        names(cell_groups[[cell_group]])
      )]
    }
  }

  gene <- features[features %in% rownames(Seurat::GetAssay(
    srt,
    assay = assay
  ))]
  gene_unique <- features_unique[features %in% rownames(Seurat::GetAssay(
    srt,
    assay = assay
  ))]
  meta <- features[features %in% colnames(srt@meta.data)]
  all_cells <- unique(unlist(lapply(cell_groups, names)))
  mat_raw <- Matrix::as.matrix(
    rbind(
      SeuratObject::GetAssayData(
        srt,
        assay = assay,
        layer = layer
      )[gene, all_cells, drop = FALSE],
      Matrix::t(srt@meta.data[all_cells, meta, drop = FALSE])
    )
  )[features, , drop = FALSE]
  rownames(mat_raw) <- features_unique
  if (isTRUE(lib_normalize) && min(mat_raw, na.rm = TRUE) >= 0) {
    if (!is.null(libsize)) {
      libsize_use <- libsize
    } else {
      libsize_use <- colSums(
        SeuratObject::GetAssayData(
          srt,
          assay = assay,
          layer = "counts"
        )[, colnames(mat_raw), drop = FALSE]
      )
      isfloat <- any(libsize_use %% 1 != 0, na.rm = TRUE)
      if (isTRUE(isfloat)) {
        libsize_use <- rep(1, length(libsize_use))
        warning(
          "The values in the 'counts' layer are non-integer. Set the library size to 1.",
          immediate. = TRUE
        )
      }
    }
    mat_raw[gene_unique, ] <- Matrix::t(
      Matrix::t(mat_raw[gene_unique, , drop = FALSE]) /
        libsize_use *
        median(libsize_use)
    )
  }

  # data used to plot heatmap
  mat_list <- list()
  for (cell_group in group.by) {
    mat_tmp <- mat_raw[, names(cell_groups[[cell_group]])]
    mat_tmp <- matrix_process(mat_tmp, method = exp_method)
    mat_tmp[is.infinite(mat_tmp)] <- max(
      abs(mat_tmp[!is.infinite(mat_tmp)]),
      na.rm = TRUE
    ) *
      ifelse(mat_tmp[is.infinite(mat_tmp)] > 0, 1, -1)
    mat_tmp[is.na(mat_tmp)] <- mean(mat_tmp, na.rm = TRUE)
    mat_list[[cell_group]] <- mat_tmp
  }

  # data used to do clustering
  # if (length(feature_split_by) == 1) {
  #   mat_split <- mat_list[[feature_split_by]]
  # } else {
  #   # mat_split <- mat_list[, unlist(lapply(cell_groups[feature_split_by], names))]
  #   # mat_split <- matrix_process(mat_split, method = exp_method)
  #   mat_split <- do.call(cbind, mat_list[feature_split_by])
  #   mat_split[is.infinite(mat_split)] <- max(abs(mat_split[!is.infinite(mat_split)])) * ifelse(mat_split[is.infinite(mat_split)] > 0, 1, -1)
  #   mat_split[is.na(mat_split)] <- mean(mat_split, na.rm = TRUE)
  # }
  mat_split <- do.call(cbind, mat_list[feature_split_by])

  if (is.null(limits)) {
    if (!is.function(exp_method) && exp_method %in% c("zscore", "log2fc")) {
      b <- ceiling(
        min(
          abs(quantile(do.call(cbind, mat_list), c(0.01, 0.99), na.rm = TRUE)),
          na.rm = TRUE
        ) *
          2
      ) /
        2
      colors <- circlize::colorRamp2(
        seq(-b, b, length = 100),
        palette_scop(palette = heatmap_palette, palcolor = heatmap_palcolor)
      )
    } else {
      b <- quantile(do.call(cbind, mat_list), c(0.01, 0.99), na.rm = TRUE)
      colors <- circlize::colorRamp2(
        seq(b[1], b[2], length = 100),
        palette_scop(palette = heatmap_palette, palcolor = heatmap_palcolor)
      )
    }
  } else {
    colors <- circlize::colorRamp2(
      seq(limits[1], limits[2], length = 100),
      palette_scop(palette = heatmap_palette, palcolor = heatmap_palcolor)
    )
  }

  cell_metadata <- cbind.data.frame(
    data.frame(row.names = colnames(mat_raw), cells = colnames(mat_raw)),
    cbind.data.frame(
      srt@meta.data[
        colnames(mat_raw),
        c(group.by, intersect(cell_annotation, colnames(srt@meta.data))),
        drop = FALSE
      ],
      Matrix::t(
        SeuratObject::GetAssayData(
          srt,
          assay = assay,
          layer = "data"
        )[
          intersect(cell_annotation, rownames(Seurat::GetAssay(
            srt,
            assay = assay
          ))),
          colnames(mat_raw),
          drop = FALSE
        ]
      )
    )
  )
  feature_metadata <- cbind.data.frame(
    data.frame(
      row.names = features_unique,
      features = features,
      features_uique = features_unique
    ),
    GetFeaturesData(srt, assay = assay)[
      features,
      c(feature_annotation),
      drop = FALSE
    ]
  )
  feature_metadata[, "duplicated"] <- feature_metadata[["features"]] %in%
    features[duplicated(features)]

  lgd <- list()
  lgd[["ht"]] <- ComplexHeatmap::Legend(title = exp_name, col_fun = colors, border = TRUE)

  ha_top_list <- NULL
  cluster_columns_list <- list()
  column_split_list <- list()
  for (i in seq_along(group.by)) {
    cell_group <- group.by[i]
    cluster_columns_list[[cell_group]] <- cluster_columns
    column_split_list[[cell_group]] <- cell_groups[[cell_group]]
    if (isTRUE(cluster_column_slices)) {
      if (!isTRUE(cluster_columns)) {
        if (nlevels(column_split_list[[cell_group]]) == 1) {
          stop(
            "cluster_column_slices=TRUE can not be used when there is only one group."
          )
        }
        dend <- ComplexHeatmap::cluster_within_group(
          mat_list[[cell_group]],
          column_split_list[[cell_group]]
        )
        cluster_columns_list[[cell_group]] <- dend
        column_split_list[[cell_group]] <- length(unique(column_split_list[[
          cell_group
        ]]))
      }
    }
    if (cell_group != "All.groups") {
      funbody <- paste0(
        "
        grid::grid.rect(gp = grid::gpar(fill = palette_scop(",
        paste0(
          "c('",
          paste0(levels(srt@meta.data[[cell_group]]), collapse = "','"),
          "')"
        ),
        ",palette = '",
        group_palette[i],
        "',palcolor=c(",
        paste0("'", paste0(group_palcolor[[i]], collapse = "','"), "'"),
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

      anno <- list()
      anno[[cell_group]] <- anno_block(
        align_to = split(
          seq_along(cell_groups[[cell_group]]),
          gsub(
            pattern = " : .*",
            replacement = "",
            x = cell_groups[[cell_group]]
          )
        ),
        panel_fun = getFunction("panel_fun", where = environment()),
        which = ifelse(flip, "row", "column"),
        show_name = FALSE
      )
      ha_cell_group <- do.call(
        ComplexHeatmap::HeatmapAnnotation,
        args = c(
          anno,
          which = ifelse(flip, "row", "column"),
          show_annotation_name = TRUE,
          annotation_name_side = ifelse(flip, "top", "left"),
          border = TRUE
        )
      )
      ha_top_list[[cell_group]] <- ha_cell_group
    }

    if (!is.null(split.by)) {
      funbody <- paste0(
        "
      grid::grid.rect(gp = grid::gpar(fill = palette_scop(",
        paste0(
          "c('",
          paste0(levels(srt@meta.data[[split.by]]), collapse = "','"),
          "')"
        ),
        ",palette = '",
        cell_split_palette,
        "',palcolor=c(",
        paste0("'", paste0(unlist(cell_split_palcolor), collapse = "','"), "'"),
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

      anno <- list()
      anno[[split.by]] <- anno_block(
        align_to = split(
          seq_along(cell_groups[[cell_group]]),
          gsub(
            pattern = ".* : ",
            replacement = "",
            x = cell_groups[[cell_group]]
          )
        ),
        panel_fun = getFunction("panel_fun", where = environment()),
        which = ifelse(flip, "row", "column"),
        show_name = i == 1
      )
      ha_split_by <- do.call(
        ComplexHeatmap::HeatmapAnnotation,
        args = c(
          anno,
          which = ifelse(flip, "row", "column"),
          show_annotation_name = TRUE,
          annotation_name_side = ifelse(flip, "top", "left"),
          border = TRUE
        )
      )
      if (is.null(ha_top_list[[cell_group]])) {
        ha_top_list[[cell_group]] <- ha_split_by
      } else {
        ha_top_list[[cell_group]] <- c(ha_top_list[[cell_group]], ha_split_by)
      }
    }
  }
  for (i in seq_along(raw.group.by)) {
    cell_group <- raw.group.by[i]
    if (cell_group != "All.groups") {
      lgd[[cell_group]] <- ComplexHeatmap::Legend(
        title = cell_group,
        labels = levels(srt@meta.data[[cell_group]]),
        legend_gp = grid::gpar(
          fill = palette_scop(
            levels(srt@meta.data[[cell_group]]),
            palette = raw.group_palette[i],
            palcolor = group_palcolor[[i]]
          )
        ),
        border = TRUE
      )
    }
  }
  if (!is.null(split.by)) {
    lgd[[split.by]] <- ComplexHeatmap::Legend(
      title = split.by,
      labels = levels(srt@meta.data[[split.by]]),
      legend_gp = grid::gpar(
        fill = palette_scop(
          levels(srt@meta.data[[split.by]]),
          palette = cell_split_palette,
          palcolor = cell_split_palcolor
        )
      ),
      border = TRUE
    )
  }

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
        for (cell_group in group.by) {
          ha_cell <- list()
          ha_cell[[cellan]] <- ComplexHeatmap::anno_simple(
            x = as.character(cell_anno[names(cell_groups[[cell_group]])]),
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
            show_annotation_name = cell_group == group.by[1],
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
          if (is.null(ha_top_list[[cell_group]])) {
            ha_top_list[[cell_group]] <- ha_top
          } else {
            ha_top_list[[cell_group]] <- c(ha_top_list[[cell_group]], ha_top)
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
        for (cell_group in group.by) {
          ha_cell <- list()
          ha_cell[[cellan]] <- ComplexHeatmap::anno_simple(
            x = cell_anno[names(cell_groups[[cell_group]])],
            col = col_fun,
            which = ifelse(flip, "row", "column"),
            na_col = "transparent",
            border = TRUE
          )
          anno_args <- c(
            ha_cell,
            which = ifelse(flip, "row", "column"),
            show_annotation_name = cell_group == group.by[1],
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
          if (is.null(ha_top_list[[cell_group]])) {
            ha_top_list[[cell_group]] <- ha_top
          } else {
            ha_top_list[[cell_group]] <- c(ha_top_list[[cell_group]], ha_top)
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

  if (is.null(feature_split)) {
    if (is.null(n_split) || isTRUE(nrow(mat_split) <= n_split)) {
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
            warning(
              "The e1071 package was not found. Switch split_method to 'kmeans'",
              immediate. = TRUE
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
                warning(
                  "fuzzification value is samller than estimated:",
                  round(min_fuzzification, 2),
                  immediate. = TRUE
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
              stop(
                "Clustering with mfuzz failed (fuzzification=",
                round(fuzzification, 2),
                "). Please set a larger fuzzification parameter manually."
              )
            }
            # mfuzz.plot(eset, cl,new.window = FALSE)
            row_split <- feature_split <- cl$cluster
          }
        }
        if (split_method == "kmeans") {
          km <- kmeans(
            mat_split,
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
          row_split <- feature_split <- cutree(hc, k = n_split)
        }
      }
      groupmean <- aggregate(
        Matrix::t(mat_split),
        by = list(unlist(cell_groups[feature_split_by])),
        mean
      )
      maxgroup <- groupmean[, 1][apply(
        groupmean[, names(row_split)],
        2,
        which.max
      )]
      df <- data.frame(row_split = row_split, order_by = maxgroup)
      df_order <- aggregate(
        df[["order_by"]],
        by = list(df[, "row_split"]),
        FUN = function(x) names(sort(table(x), decreasing = TRUE))[1]
      )
      df_order[, "row_split"] <- df_order[, "Group.1"]
      df_order[["order_by"]] <- as.numeric(factor(
        df_order[["x"]],
        levels = levels(maxgroup)
      ))
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
    row_split_raw <- row_split <- feature_split <- feature_split[row.names(
      mat_split
    )]
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
        dend <- ComplexHeatmap::cluster_within_group(Matrix::t(mat_split), row_split_raw)
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
      features_split = anno_block(
        align_to = split(seq_along(row_split_raw), row_split_raw),
        panel_fun = getFunction("panel_fun", where = environment()),
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
    mat_cluster <- do.call(cbind, mat_list[cluster_features_by])
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
      dend <- cluster_within_group2(Matrix::t(mat_cluster), row_split_raw)
      cluster_rows <- dend
    }
  }

  cell_group <- group.by[1]
  ht_args <- list(
    name = cell_group,
    matrix = mat_list[[cell_group]],
    col = colors,
    row_split = row_split,
    column_split = column_split_list[[cell_group]],
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns_list[[cell_group]],
    cluster_row_slices = cluster_row_slices,
    cluster_column_slices = cluster_column_slices,
    use_raster = TRUE
  )
  ht_args <- c(ht_args, ht_params[setdiff(names(ht_params), names(ht_args))])
  ht_list <- do.call(ComplexHeatmap::Heatmap, args = ht_args)
  features_ordered <- rownames(
    mat_list[[1]]
  )[unlist(
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
        index <- unique(round(seq(
          from = index_from,
          to = index_to,
          length.out = nlabel
        )))
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
      warning(
        paste0(paste0(drop, collapse = ","), "was not found in the features"),
        immediate. = TRUE
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
        ha_feature <- do.call(ComplexHeatmap::HeatmapAnnotation, args = anno_args)
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
          colors = palette_scop(palette = palette, palcolor = palcolor)
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
  for (cell_group in group.by) {
    if (cell_group == group.by[1]) {
      left_annotation <- ha_left
    } else {
      left_annotation <- NULL
    }
    if (cell_group == group.by[length(group.by)]) {
      right_annotation <- ha_right
    } else {
      right_annotation <- NULL
    }

    ht_args <- list(
      name = cell_group,
      matrix = if (flip) Matrix::t(mat_list[[cell_group]]) else mat_list[[cell_group]],
      col = colors,
      row_title = row_title %||%
        if (flip) {
          ifelse(cell_group != "All.groups", cell_group, "")
        } else {
          character(0)
        },
      row_title_side = row_title_side,
      column_title = column_title %||%
        if (flip) {
          character(0)
        } else {
          ifelse(cell_group != "All.groups", cell_group, "")
        },
      column_title_side = column_title_side,
      row_title_rot = row_title_rot,
      column_title_rot = column_title_rot,
      row_split = if (flip) column_split_list[[cell_group]] else row_split,
      column_split = if (flip) row_split else column_split_list[[cell_group]],
      cluster_rows = if (flip) {
        cluster_columns_list[[cell_group]]
      } else {
        cluster_rows
      },
      cluster_columns = if (flip) {
        cluster_rows
      } else {
        cluster_columns_list[[cell_group]]
      },
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
      top_annotation = if (flip) left_annotation else ha_top_list[[cell_group]],
      left_annotation = if (flip) {
        ha_top_list[[cell_group]]
      } else {
        left_annotation
      },
      bottom_annotation = if (flip) right_annotation else NULL,
      right_annotation = if (flip) NULL else right_annotation,
      show_heatmap_legend = FALSE,
      border = border,
      use_raster = use_raster,
      raster_device = raster_device,
      raster_by_magick = raster_by_magick,
      width = if (is.numeric(width[cell_group])) {
        grid::unit(width[cell_group], units = units)
      } else {
        NULL
      },
      height = if (is.numeric(height[cell_group])) {
        grid::unit(height[cell_group], units = units)
      } else {
        NULL
      }
    )
    if (!is.null(split.by) && !isTRUE(cluster_column_slices)) {
      groups_order <- sapply(
        strsplit(levels(column_split_list[[cell_group]]), " : "),
        function(x) x[[1]]
      )
      gaps_order <- paste(
        groups_order[2:length(groups_order)],
        groups_order[1:(length(groups_order) - 1)],
        sep = "->"
      )
      gaps <- rep(grid::unit(1, "mm"), length(gaps_order))
      gaps[
        groups_order[2:length(groups_order)] ==
          groups_order[1:(length(groups_order) - 1)]
      ] <- grid::unit(0, "mm")
      if (isTRUE(flip)) {
        ht_args[["row_gap"]] <- gaps
      } else {
        ht_args[["column_gap"]] <- gaps
      }
    }
    if (any(names(ht_params) %in% names(ht_args))) {
      warning(
        "ht_params: ",
        paste0(intersect(names(ht_params), names(ht_args)), collapse = ","),
        " were duplicated and will not be used.",
        immediate. = TRUE
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
      message(
        "The size of the heatmap is fixed because certain elements are not scalable."
      )
      message(
        "The width and height of the heatmap are determined by the size of the current viewport."
      )
      message(
        "If you want to have more control over the size, you can manually set the parameters 'width' and 'height'."
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
      matrix_list = mat_list,
      feature_split = feature_split,
      cell_metadata = cell_metadata,
      feature_metadata = feature_metadata,
      enrichment = res
    )
  )
}

heatmap_fixsize <- function(
    width,
    width_sum,
    height,
    height_sum,
    units,
    ht_list,
    legend_list) {
  ht <- ComplexHeatmap::draw(ht_list, annotation_legend_list = legend_list)
  ht_width <- ComplexHeatmap:::width(ht)
  ht_height <- ComplexHeatmap:::height(ht)
  # g_tree <- grid::grid.grabExpr(
  #   {
  #     ht <- ComplexHeatmap::draw(ht_list, annotation_legend_list = legend_list)
  #     ht_width <- ComplexHeatmap:::width(ht)
  #     ht_height <- ComplexHeatmap:::height(ht)
  #     if (inherits(ht_list, "HeatmapList")) {
  #       for (nm in names(ht_list@ht_list)) {
  #         if (is.null(names(width))) {
  #           width_fix <- width[1]
  #         } else {
  #           width_fix <- width[nm]
  #         }
  #         if (is.null(names(height))) {
  #           height_fix <- height[1]
  #         } else {
  #           height_fix <- height[nm]
  #         }
  #         ht_list@ht_list[[nm]]@matrix_param$width <- grid::unit(
  #           width_fix %||% dim(ht_list@ht_list[[nm]]@matrix)[1],
  #           units = "null"
  #         )
  #         ht_list@ht_list[[nm]]@matrix_param$height <- grid::unit(
  #           height_fix %||% dim(ht_list@ht_list[[nm]]@matrix)[2],
  #           units = "null"
  #         )
  #       }
  #     } else if (inherits(ht_list, "Heatmap")) {
  #       ht_list@matrix_param$width <- grid::unit(
  #         width[1] %||% dim(ht_list@matrix)[1],
  #         units = "null"
  #       )
  #       ht_list@matrix_param$height <- grid::unit(
  #         height[1] %||% dim(ht_list@matrix)[2],
  #         units = "null"
  #       )
  #     } else {
  #       stop("ht_list is not a class of HeatmapList or Heatmap.")
  #     }
  #   },
  #   width = grid::unit(width_sum, units = units),
  #   height = grid::unit(height_sum, units = units),
  #   wrap = TRUE,
  #   wrap.grobs = TRUE
  # )
  if (grid::unitType(ht_width) == "npc") {
    ht_width <- grid::unit(width_sum, units = units)
  }
  if (grid::unitType(ht_height) == "npc") {
    ht_height <- grid::unit(height_sum, units = units)
  }
  if (is.null(width)) {
    ht_width <- max(
      grid::convertWidth(
        ht@layout$max_left_component_width,
        units,
        valueOnly = TRUE
      ) +
        grid::convertWidth(
          ht@layout$max_right_component_width,
          units,
          valueOnly = TRUE
        ) +
        grid::convertWidth(
          sum(ht@layout$max_title_component_width),
          units,
          valueOnly = TRUE
        ) +
        grid::convertWidth(
          ht@annotation_legend_param$size[1],
          units,
          valueOnly = TRUE
        ) +
        grid::convertWidth(grid::unit(1, "in"), units, valueOnly = TRUE),
      grid::convertWidth(grid::unit(0.95, "npc"), units, valueOnly = TRUE)
    )
    ht_width <- grid::unit(ht_width, units)
  }
  if (is.null(height)) {
    ht_height <- max(
      grid::convertHeight(
        ht@layout$max_top_component_height,
        units,
        valueOnly = TRUE
      ) +
        grid::convertHeight(
          ht@layout$max_bottom_component_height,
          units,
          valueOnly = TRUE
        ) +
        grid::convertHeight(
          sum(ht@layout$max_title_component_height),
          units,
          valueOnly = TRUE
        ) +
        grid::convertHeight(grid::unit(1, "in"), units, valueOnly = TRUE),
      grid::convertHeight(
        ht@annotation_legend_param$size[2],
        units,
        valueOnly = TRUE
      ),
      grid::convertHeight(grid::unit(0.95, "npc"), units, valueOnly = TRUE)
    )
    ht_height <- grid::unit(ht_height, units)
  }
  ht_width <- convertUnit(ht_width, unitTo = units)
  ht_height <- convertUnit(ht_height, unitTo = units)

  return(list(ht_width = ht_width, ht_height = ht_height))
}
