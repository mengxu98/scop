#' @title The Cell Correlation Heatmap
#'
#' @description
#' This function generates a heatmap to visualize the similarity between different cell types or conditions.
#' It takes in Seurat objects or expression matrices as input and calculates pairwise similarities or distance.
#'
#' @md
#' @param srt_query A Seurat object or count matrix representing the query dataset.
#' This dataset will be used to calculate the similarities between cells.
#' @param srt_ref A Seurat object or count matrix representing the reference dataset.
#' If provided, the similarities will be calculated between cells from the query and reference datasets.
#' If not provided, the similarities will be calculated within the query dataset.
#' @param bulk_ref A count matrix representing bulk data.
#' If provided, the similarities will be calculated between cells from the query dataset and bulk data.
#' @param query_group The grouping variable in the query dataset.
#' This variable will be used to group cells in the heatmap rows.
#' If not provided, all cells will be treated as one group.
#' @param ref_group The grouping variable in the reference dataset.
#' This variable will be used to group cells in the heatmap columns.
#' If not provided, all cells will be treated as one group.
#' @param query_assay The assay to use for the query dataset.
#' If not provided, the default assay of the query dataset will be used.
#' @param ref_assay The assay to use for the reference dataset.
#' If not provided, the default assay of the reference dataset will be used.
#' @param query_reduction The dimensionality reduction method to use for the query dataset.
#' If not provided, no dimensionality reduction will be applied to the query dataset.
#' @param ref_reduction The dimensionality reduction method to use for the reference dataset.
#' If not provided, no dimensionality reduction will be applied to the reference dataset.
#' @param query_dims The dimensions to use for the query dataset.
#' If not provided, the first 30 dimensions will be used.
#' @param ref_dims The dimensions to use for the reference dataset.
#' If not provided, the first 30 dimensions will be used.
#' @param query_collapsing Whether to collapse cells within each query group before calculating similarities.
#' If set to TRUE, the similarities will be calculated between query groups rather than individual cells.
#' @param ref_collapsing Detail description of the `query_collapsing` argument.
#' @param features A vector of feature names to include in the heatmap.
#' If not provided, a default set of highly variable features (HVF) will be used.
#' @param features_type The type of features to use.
#' Options are `"HVF"` for highly variable features,
#' `"DE"` for differentially expressed features between query and reference groups.
#' @param feature_source The source of features to use.
#' Options are `"query"` to use only features from the query dataset,
#' `"ref"` to use only features from the reference dataset, or `"both"` to use features from both datasets.
#' If not provided or set to "both", features will be selected from both datasets.
#' @param nfeatures The maximum number of features to include in the heatmap.
#' Default is 2000.
#' @param DEtest_param The parameters to use for differential expression testing.
#' This should be a list with two elements:
#' `"max.cells.per.ident"` specifying the maximum number of cells per group for differential expression testing,
#' and `"test.use"` specifying the statistical test to use for differential expression testing.
#' Default parameters will be used.
#' @param DE_threshold The threshold for differential expression.
#' Only features with adjusted p-values below this threshold will be considered differentially expressed.
#' @param distance_metric The distance metric to use for calculating similarities between cells.
#' This can be any of the following:
#' `"cosine"`, `"pearson"`, `"spearman"`, `"correlation"`, `"jaccard"`, `"ejaccard"`, `"dice"`, `"edice"`, `"hamman"`, `"simple matching"`, or `"faith"`.
#' Dhe default is `"cosine"`.
#' @param k The number of nearest neighbors to use for calculating similarities.
#' Default is 30.
#' @param filter_lowfreq The minimum frequency threshold for selecting query dataset features.
#' Features with a frequency below this threshold will be excluded from the heatmap.
#' Default is 0.
#' @param prefix The prefix to use for the KNNPredict tool layer in the query dataset.
#' This can be used to avoid conflicts with other tools in the Seurat object.
#' The default is `"KNNPredict"`.
#' @param exp_legend_title The title for the color legend in the heatmap.
#' If not provided, a default title based on the similarity metric will be used.
#' @param border Whether to add a border around each heatmap cell.
#' The default is `TRUE`.
#' @param flip Whether to flip the orientation of the heatmap.
#' If set to TRUE, the rows and columns of the heatmap will be swapped.
#' This can be useful for visualizing large datasets in a more compact form.
#' The default is `FALSE`.
#' @param limits The limits for the color scale in the heatmap.
#' If not provided, the default is to use the range of similarity values.
#' @param cluster_rows Whether to cluster the rows of the heatmap.
#' If set to TRUE, the rows will be rearranged based on hierarchical clustering.
#' The default is `FALSE`.
#' @param cluster_columns Whether to cluster the columns of the heatmap.
#' If set to TRUE, the columns will be rearranged based on hierarchical clustering.
#' The default is `FALSE`.
#' @param show_row_names Whether to show the row names in the heatmap.
#' The default is `FALSE`.
#' @param show_column_names Whether to show the column names in the heatmap.
#' The default is `FALSE`.
#' @param row_names_side The side of the heatmap to show the row names.
#' Options are `"left"` or `"right"`.
#' Default is `"left"`.
#' @param column_names_side The side of the heatmap to show the column names.
#' Options are `"top"` or `"bottom"`.
#' If not provided, Default is `"top"`.
#' @param row_names_rot The rotation angle of the row names.
#' If not provided, Default is `0` degrees.
#' @param column_names_rot The rotation angle of the column names.
#' If not provided, the default is 90 degrees.
#' @param row_title The title for the row names in the heatmap.
#' If not provided, the default is to use the query grouping variable.
#' @param column_title The title for the column names in the heatmap.
#' Default is to use the reference grouping variable.
#' @param row_title_side The side of the heatmap to show the row title.
#' Options are `"top"` or `"bottom"`.
#' Default is `"left"`.
#' @param column_title_side The side of the heatmap to show the column title.
#' Options are `"left"` or `"right"`.
#' Default is `"top"`.
#' @param row_title_rot The rotation angle of the row title.
#' Default is `90` degrees.
#' @param column_title_rot The rotation angle of the column title.
#' Default is `0` degrees.
#' @param nlabel The maximum number of labels to show on each side of the heatmap.
#' If set to 0, no labels will be shown.
#' This can be useful for reducing clutter in large heatmaps.
#' Default is 0.
#' @param label_cutoff The similarity cutoff for showing labels.
#' Only cells with similarity values above this cutoff will have labels.
#' Default is 0.
#' @param label_by The dimension to use for labeling cells.
#' Options are `"row"` to label cells by row, `"column"` to label cells by column, or `"both"` to label cells by both row and column.
#' Default is `"row"`.
#' @param label_size The size of the labels.
#' Default is `10`.
#' @param heatmap_palette The color palette to use for the heatmap.
#' This can be any of the palettes available in the circlize package.
#' Default is `"RdBu"`.
#' @param heatmap_palcolor The specific colors to use for the heatmap palette.
#' This should be a vector of color names or RGB values.
#' Default is `NULL`.
#' @param query_group_palette The color palette to use for the query group legend.
#' This can be any of the palettes available in the circlize package.
#' Default is `"Paired"`.
#' @param query_group_palcolor The specific colors to use for the query group palette.
#' This should be a vector of color names or RGB values.
#' Default is `NULL`.
#' @param ref_group_palette The color palette to use for the reference group legend.
#' This can be any of the palettes available in the circlize package.
#' Default is `"simspec"`.
#' @param ref_group_palcolor The specific colors to use for the reference group palette.
#' This should be a vector of color names or RGB values.
#' Default is `NULL`.
#' @param query_cell_annotation A vector of cell metadata column names or assay feature names to use for highlighting specific cells in the heatmap.
#' Each element of the vector will create a separate cell annotation track in the heatmap.
#' If not provided, no cell annotations will be shown.
#' @param query_cell_annotation_palette The color palette to use for the query cell annotation tracks.
#' This can be any of the palettes available in the circlize package.
#' If a single color palette is provided, it will be used for all cell annotation tracks.
#' If multiple color palettes are provided, each track will be assigned a separate palette.
#' Default is `"Paired"`.
#' @param query_cell_annotation_palcolor The specific colors to use for the query cell annotation palettes.
#' This should be a list of vectors, where each vector contains the colors for a specific cell annotation track.
#' If a single color vector is provided, it will be used for all cell annotation tracks.
#' Default is `NULL`.
#' @param query_cell_annotation_params Additional parameters to customize the appearance of the query cell annotation tracks.
#' This should be a list with named elements, where the names correspond to parameter names in the [ComplexHeatmap::Heatmap] function.
#' Any conflicting parameters will override the defaults set by this function.
#' @param ref_cell_annotation_params Detail description of the `query_cell_annotation_params` argument.
#' @param ref_cell_annotation A vector of cell metadata column names or assay feature names to use for highlighting specific cells in the heatmap.
#' Each element of the vector will create a separate cell annotation track in the heatmap.
#' If not provided, no cell annotations will be shown.
#' @param ref_cell_annotation_palette The color palette to use for the reference cell annotation tracks.
#' This can be any of the palettes available in the circlize package.
#' If a single color palette is provided, it will be used for all cell annotation tracks.
#' If multiple color palettes are provided, each track will be assigned a separate palette.
#' Default is `"Paired"`.
#' @param ref_cell_annotation_palcolor The specific colors to use for the reference cell annotation palettes.
#' This should be a list of vectors, where each vector contains the colors for a specific cell annotation track.
#' If a single color vector is provided, it will be used for all cell annotation tracks.
#' If multiple color vectors are provided, each track will be assigned a separate color vector.
#' Default is `NULL`.
#' @param use_raster Whether to use raster images for rendering the heatmap.
#' If set to `TRUE`, the heatmap will be rendered as a raster image using the raster_device argument.
#' Default is determined based on the number of rows and columns in the heatmap.
#' @param raster_device The raster device to use for rendering the heatmap.
#' This should be a character string specifying the device name, such as `"png"`, `"jpeg"`, or `"pdf"`.
#' Default is `"png"`.
#' @param raster_by_magick Whether to use the magick package for rendering rasters.
#' If set to `TRUE`, the magick package will be used instead of the raster package.
#' This can be useful for rendering large heatmaps more efficiently.
#' If the magick package is not installed, this argument will be ignored.
#' @param width The width of the heatmap in the specified units.
#' If not provided, the width will be automatically determined based on the number of columns in the heatmap and the default unit.
#' @param height The height of the heatmap in the specified units.
#' If not provided, the height will be automatically determined based on the number of rows in the heatmap and the default unit.
#' @param units The units to use for the width and height of the heatmap.
#' Default is `"inch"`, Options are `"mm"`, `"cm"`, or `"inch"`.
#' @param seed The random seed to use for reproducible results.
#' Default is `11`.
#' @param ht_params Additional parameters to customize the appearance of the heatmap.
#' This should be a list with named elements, where the names correspond to parameter names in the [ComplexHeatmap::Heatmap] function.
#' Any conflicting parameters will override the defaults set by this function.
#'
#' @return
#' A list with the following elements:
#' \itemize{
#'   \item \code{plot:} The heatmap plot as a ggplot object.
#'   \item \code{features:} The features used in the heatmap.
#'   \item \code{simil_matrix:} The similarity matrix used to generate the heatmap.
#'   \item \code{simil_name:} The name of the similarity metric used to generate the heatmap.
#'   \item \code{cell_metadata:} The cell metadata used to generate the heatmap.
#' }
#'
#' @seealso
#' [RunKNNMap], [RunKNNPredict]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' ht1 <- CellCorHeatmap(
#'   srt_query = pancreas_sub,
#'   query_group = "SubCellType"
#' )
#' ht1$plot
#'
#' data(panc8_sub)
#' # Simply convert genes from human to mouse and preprocess the data
#' genenames <- make.unique(
#'   thisutils::capitalize(
#'     rownames(panc8_sub),
#'     force_tolower = TRUE
#'   )
#' )
#' names(genenames) <- rownames(panc8_sub)
#' panc8_sub <- RenameFeatures(
#'   panc8_sub,
#'   newnames = genenames
#' )
#' panc8_sub <- CheckDataMerge(
#'   panc8_sub,
#'   batch = "tech"
#' )[["srt_merge"]]
#'
#' ht2 <- CellCorHeatmap(
#'   srt_query = pancreas_sub,
#'   srt_ref = panc8_sub,
#'   nlabel = 3,
#'   label_cutoff = 0.6,
#'   query_group = "SubCellType",
#'   ref_group = "celltype",
#'   query_cell_annotation = "Phase",
#'   query_cell_annotation_palette = "Set2",
#'   ref_cell_annotation = "tech",
#'   ref_cell_annotation_palette = "Set3",
#'   width = 4,
#'   height = 4
#' )
#' ht2$plot
#'
#' ht3 <- CellCorHeatmap(
#'   srt_query = pancreas_sub,
#'   srt_ref = panc8_sub,
#'   query_group = "SubCellType",
#'   query_collapsing = FALSE,
#'   cluster_rows = TRUE,
#'   ref_group = "celltype",
#'   ref_collapsing = FALSE,
#'   cluster_columns = TRUE
#' )
#' ht3$plot
#'
#' ht4 <- CellCorHeatmap(
#'   srt_query = pancreas_sub,
#'   srt_ref = panc8_sub,
#'   show_row_names = TRUE,
#'   show_column_names = TRUE,
#'   query_group = "SubCellType",
#'   ref_group = "celltype",
#'   query_cell_annotation = c(
#'     "Sox9", "Rbp4", "Gcg", "Nap1l2", "Xist"
#'   ),
#'   ref_cell_annotation = c(
#'     "Sox9", "Rbp4", "Gcg", "Nap1l2", "Xist"
#'   ),
#'   height = 2.5,
#'   width = 3.5
#' )
#' ht4$plot
CellCorHeatmap <- function(
    srt_query,
    srt_ref = NULL,
    bulk_ref = NULL,
    query_group = NULL,
    ref_group = NULL,
    query_assay = NULL,
    ref_assay = NULL,
    query_reduction = NULL,
    ref_reduction = NULL,
    query_dims = 1:30,
    ref_dims = 1:30,
    query_collapsing = !is.null(query_group),
    ref_collapsing = TRUE,
    features = NULL,
    features_type = c("HVF", "DE"),
    feature_source = "both",
    nfeatures = 2000,
    DEtest_param = list(
      max.cells.per.ident = 200,
      test.use = "wilcox"
    ),
    DE_threshold = "p_val_adj < 0.05",
    distance_metric = "cosine",
    k = 30,
    filter_lowfreq = 0,
    prefix = "KNNPredict",
    exp_legend_title = NULL,
    border = TRUE,
    flip = FALSE,
    limits = NULL,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_names_side = "left",
    column_names_side = "top",
    row_names_rot = 0,
    column_names_rot = 90,
    row_title = NULL,
    column_title = NULL,
    row_title_side = "left",
    column_title_side = "top",
    row_title_rot = 90,
    column_title_rot = 0,
    nlabel = 0,
    label_cutoff = 0,
    label_by = "row",
    label_size = 10,
    heatmap_palette = "RdBu",
    heatmap_palcolor = NULL,
    query_group_palette = "Paired",
    query_group_palcolor = NULL,
    ref_group_palette = "simspec",
    ref_group_palcolor = NULL,
    query_cell_annotation = NULL,
    query_cell_annotation_palette = "Paired",
    query_cell_annotation_palcolor = NULL,
    query_cell_annotation_params = if (flip) {
      list(height = grid::unit(10, "mm"))
    } else {
      list(width = grid::unit(10, "mm"))
    },
    ref_cell_annotation = NULL,
    ref_cell_annotation_palette = "Paired",
    ref_cell_annotation_palcolor = NULL,
    ref_cell_annotation_params = if (flip) {
      list(width = grid::unit(10, "mm"))
    } else {
      list(height = grid::unit(10, "mm"))
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

  ref_legend <- TRUE
  if (is.null(srt_ref) && is.null(bulk_ref)) {
    srt_ref <- srt_query
    ref_group <- query_group
    ref_assay <- query_assay
    ref_reduction <- query_reduction
    ref_dims <- query_dims
    ref_collapsing <- query_collapsing
    ref_group_palette <- query_group_palette
    ref_group_palcolor <- query_group_palcolor
    ref_cell_annotation <- query_cell_annotation
    ref_cell_annotation_palette <- query_cell_annotation_palette
    ref_cell_annotation_palcolor <- query_cell_annotation_palcolor
    ref_legend <- FALSE
  }
  if (!is.null(bulk_ref)) {
    srt_ref <- Seurat::CreateSeuratObject(
      counts = bulk_ref,
      meta.data = data.frame(celltype = colnames(bulk_ref)),
      assay = "RNA"
    )
    ref_group <- "CellType"
    ref_assay <- "RNA"
    ref_reduction <- NULL
    ref_collapsing <- FALSE
  }
  query_assay <- query_assay %||% SeuratObject::DefaultAssay(srt_query)
  ref_assay <- ref_assay %||% SeuratObject::DefaultAssay(srt_ref)
  other_params <- list(
    query_group = query_group,
    query_reduction = query_reduction,
    query_assay = query_assay,
    query_dims = query_dims,
    query_collapsing = query_collapsing,
    ref_group = ref_group,
    ref_reduction = ref_reduction,
    ref_assay = ref_assay,
    ref_dims = ref_dims,
    ref_collapsing = ref_collapsing
  )

  srt_query_tools <- srt_query@tools[[paste0(prefix, "_classification")]]
  if (is.null(srt_query_tools[["distance_matrix"]]) || !identical(other_params, srt_query_tools[["other_params"]])) {
    srt_query <- RunKNNPredict(
      srt_query = srt_query,
      srt_ref = srt_ref,
      query_group = query_group,
      ref_group = ref_group,
      query_assay = query_assay,
      ref_assay = ref_assay,
      query_reduction = query_reduction,
      ref_reduction = ref_reduction,
      query_dims = query_dims,
      ref_dims = ref_dims,
      query_collapsing = query_collapsing,
      ref_collapsing = ref_collapsing,
      features = features,
      features_type = features_type,
      feature_source = feature_source,
      nfeatures = nfeatures,
      DEtest_param = DEtest_param,
      DE_threshold = DE_threshold,
      distance_metric = distance_metric,
      k = k,
      filter_lowfreq = filter_lowfreq,
      prefix = prefix,
      nn_method = "raw",
      return_full_distance_matrix = TRUE
    )
  }

  srt_query_tools <- srt_query@tools[[paste0(prefix, "_classification")]]
  distance_matrix <- srt_query_tools[["distance_matrix"]]
  distance_metric <- srt_query_tools[["distance_metric"]]
  simil_methods <- c(
    "cosine",
    "pearson",
    "spearman",
    "correlation",
    "jaccard",
    "ejaccard",
    "dice",
    "edice",
    "hamman",
    "simple matching",
    "faith"
  )
  dist_methods <- c(
    "euclidean",
    "chisquared",
    "kullback",
    "manhattan",
    "maximum",
    "canberra",
    "minkowski",
    "hamming"
  )
  if (distance_metric %in% simil_methods) {
    simil_matrix <- Matrix::t(
      as_matrix(1 - distance_matrix)
    )
    simil_name <- paste0(capitalize(distance_metric), " similarity")
  } else if (distance_metric %in% dist_methods) {
    simil_matrix <- Matrix::t(
      as_matrix(
        1 - distance_matrix / max(distance_matrix, na.rm = TRUE)
      )
    )
    simil_name <- paste0(
      "1-dist[",
      distance_metric,
      "]/max(dist[",
      distance_metric,
      "])"
    )
  }
  simil_matrix[is.infinite(simil_matrix)] <- max(
    abs(simil_matrix[!is.infinite(simil_matrix)]),
    na.rm = TRUE
  ) *
    ifelse(simil_matrix[is.infinite(simil_matrix)] > 0, 1, -1)
  simil_matrix[is.na(simil_matrix)] <- 0
  exp_name <- exp_legend_title %||% simil_name

  cell_groups <- list()
  if (is.null(query_group)) {
    srt_query@meta.data[["All.groups"]] <- factor("")
    query_group <- "All.groups"
  }
  if (is.null(ref_group)) {
    srt_ref@meta.data[["All.groups"]] <- factor("")
    ref_group <- "All.groups"
  }
  if (!is.factor(srt_query[[query_group, drop = TRUE]])) {
    srt_query@meta.data[[query_group]] <- factor(
      srt_query[[query_group, drop = TRUE]],
      levels = unique(srt_query[[query_group, drop = TRUE]])
    )
  }
  cell_groups[["query_group"]] <- unlist(
    lapply(
      levels(
        srt_query[[query_group, drop = TRUE]]
      ),
      function(x) {
        cells_sub <- colnames(srt_query)[which(
          srt_query[[query_group, drop = TRUE]] == x
        )]
        stats::setNames(
          object = rep(x, length(cells_sub)),
          nm = cells_sub
        )
      }
    ),
    use.names = TRUE
  )

  levels <- levels(srt_query[[query_group, drop = TRUE]])
  cell_groups[["query_group"]] <- factor(
    cell_groups[["query_group"]],
    levels = levels[levels %in% cell_groups[["query_group"]]]
  )

  if (!is.factor(srt_ref[[ref_group, drop = TRUE]])) {
    srt_ref@meta.data[[ref_group]] <- factor(
      srt_ref[[ref_group, drop = TRUE]],
      levels = unique(srt_ref[[ref_group, drop = TRUE]])
    )
  }
  cell_groups[["ref_group"]] <- unlist(
    lapply(levels(srt_ref[[ref_group, drop = TRUE]]), function(x) {
      cells_sub <- colnames(srt_ref)[which(
        srt_ref[[ref_group, drop = TRUE]] == x
      )]
      stats::setNames(
        object = rep(x, length(cells_sub)),
        nm = cells_sub
      )
    }),
    use.names = TRUE
  )

  levels <- levels(srt_ref[[ref_group, drop = TRUE]])
  cell_groups[["ref_group"]] <- factor(
    cell_groups[["ref_group"]],
    levels = levels[levels %in% cell_groups[["ref_group"]]]
  )

  if (isTRUE(query_collapsing)) {
    simil_matrix <- simil_matrix[
      intersect(
        levels(cell_groups[["query_group"]]),
        rownames(simil_matrix)
      ), ,
      drop = FALSE
    ]
  } else {
    simil_matrix <- simil_matrix[
      intersect(
        names(cell_groups[["query_group"]]),
        rownames(simil_matrix)
      ), ,
      drop = FALSE
    ]
  }

  if (isTRUE(ref_collapsing)) {
    simil_matrix <- simil_matrix[,
      intersect(
        levels(cell_groups[["ref_group"]]),
        colnames(simil_matrix)
      ),
      drop = FALSE
    ]
  } else {
    simil_matrix <- simil_matrix[,
      intersect(
        names(cell_groups[["ref_group"]]),
        colnames(simil_matrix)
      ),
      drop = FALSE
    ]
  }

  if (!is.null(query_cell_annotation)) {
    if (length(query_cell_annotation_palette) == 1) {
      query_cell_annotation_palette <- rep(
        query_cell_annotation_palette,
        length(query_cell_annotation)
      )
    }
    if (length(query_cell_annotation_palcolor) == 1) {
      query_cell_annotation_palcolor <- rep(
        query_cell_annotation_palcolor,
        length(query_cell_annotation)
      )
    }
    npal <- unique(c(
      length(query_cell_annotation_palette),
      length(query_cell_annotation_palcolor),
      length(query_cell_annotation)
    ))
    if (length(npal[npal != 0]) > 1) {
      log_message(
        "query_cell_annotation_palette and query_cell_annotation_palcolor must be the same length as query_cell_annotation",
        message_type = "error"
      )
    }
    if (any(!query_cell_annotation %in% c(colnames(srt_query@meta.data), rownames(srt_query[[query_assay]])))) {
      log_message(
        "query_cell_annotation: ",
        paste0(
          query_cell_annotation[
            !query_cell_annotation %in%
              c(
                colnames(srt_query@meta.data),
                rownames(srt_query[[query_assay]])
              )
          ],
          collapse = ","
        ),
        " is not in the Seurat object.",
        message_type = "error"
      )
    }
  }
  if (!is.null(ref_cell_annotation)) {
    if (length(ref_cell_annotation_palette) == 1) {
      ref_cell_annotation_palette <- rep(
        ref_cell_annotation_palette,
        length(ref_cell_annotation)
      )
    }
    if (length(ref_cell_annotation_palcolor) == 1) {
      ref_cell_annotation_palcolor <- rep(
        ref_cell_annotation_palcolor,
        length(ref_cell_annotation)
      )
    }
    npal <- unique(c(
      length(ref_cell_annotation_palette),
      length(ref_cell_annotation_palcolor),
      length(ref_cell_annotation)
    ))
    if (length(npal[npal != 0]) > 1) {
      log_message(
        "ref_cell_annotation_palette and ref_cell_annotation_palcolor must be the same length as ref_cell_annotation",
        message_type = "error"
      )
    }
    if (
      any(
        !ref_cell_annotation %in%
          c(colnames(srt_ref@meta.data), rownames(srt_ref[[ref_assay]]))
      )
    ) {
      log_message(
        "ref_cell_annotation: ",
        paste0(
          ref_cell_annotation[
            !ref_cell_annotation %in%
              c(colnames(srt_ref@meta.data), rownames(srt_ref[[ref_assay]]))
          ],
          collapse = ","
        ),
        " is not in the Seurat object.",
        message_type = "error"
      )
    }
  }

  if (isTRUE(flip)) {
    cluster_rows_raw <- cluster_rows
    cluster_columns_raw <- cluster_columns
    cluster_rows <- cluster_columns_raw
    cluster_columns <- cluster_rows_raw
  }
  if (is.null(limits)) {
    colors <- circlize::colorRamp2(
      seq(
        min(simil_matrix, na.rm = TRUE),
        max(simil_matrix, na.rm = TRUE),
        length = 100
      ),
      palette_colors(palette = heatmap_palette, palcolor = heatmap_palcolor)
    )
  } else {
    colors <- circlize::colorRamp2(
      seq(limits[1], limits[2], length = 100),
      palette_colors(palette = heatmap_palette, palcolor = heatmap_palcolor)
    )
  }

  cell_metadata <- data.frame(
    row.names = c(
      paste0("query_", colnames(srt_query)),
      paste0("ref_", colnames(srt_ref))
    ),
    cells = c(colnames(srt_query), colnames(srt_ref))
  )
  query_metadata <- cbind.data.frame(
    srt_query@meta.data[
      cell_metadata[["cells"]],
      c(
        query_group,
        intersect(query_cell_annotation, colnames(srt_query@meta.data))
      ),
      drop = FALSE
    ],
    as.data.frame(
      Matrix::t(
        GetAssayData5(
          object = srt_ref,
          assay = ref_assay,
          layer = "data",
          verbose = FALSE
        )[
          intersect(
            query_cell_annotation,
            rownames(srt_query[[query_assay]])
          ) %||%
            integer(), ,
          drop = FALSE
        ]
      )
    )[cell_metadata[["cells"]], , drop = FALSE]
  )
  colnames(query_metadata) <- paste0("query_", colnames(query_metadata))
  ref_metadata <- cbind.data.frame(
    srt_ref@meta.data[
      cell_metadata[["cells"]],
      c(ref_group, intersect(ref_cell_annotation, colnames(srt_ref@meta.data))),
      drop = FALSE
    ],
    as.data.frame(
      Matrix::t(
        GetAssayData5(
          object = srt_ref,
          assay = ref_assay,
          layer = "data",
          verbose = FALSE
        )[
          intersect(
            ref_cell_annotation,
            rownames(srt_ref[[ref_assay]])
          ) %||%
            integer(), ,
          drop = FALSE
        ]
      )
    )[cell_metadata[["cells"]], , drop = FALSE]
  )
  colnames(ref_metadata) <- paste0("ref_", colnames(ref_metadata))
  cell_metadata <- cbind.data.frame(
    cell_metadata,
    cbind.data.frame(query_metadata, ref_metadata)
  )

  lgd <- list()
  lgd[["ht"]] <- ComplexHeatmap::Legend(
    title = exp_name,
    col_fun = colors,
    border = TRUE
  )

  ha_query_list <- NULL
  if (query_group != "All.groups") {
    if (
      isFALSE(query_collapsing) &&
        ((isFALSE(flip) & isTRUE(cluster_rows)) ||
          (isTRUE(flip) & isTRUE(cluster_columns)))
    ) {
      query_cell_annotation <- c(query_group, query_cell_annotation)
      query_cell_annotation_palette <- c(
        query_group_palette,
        query_cell_annotation_palette
      )
      query_cell_annotation_palcolor <- c(
        list(query_group_palcolor),
        query_cell_annotation_palcolor %||% list(NULL)
      )
    } else {
      funbody <- paste0(
        "
        grid::grid.rect(gp = grid::gpar(fill = palette_colors(",
        paste0(
          "c('",
          paste0(
            levels(srt_query[[query_group, drop = TRUE]]),
            collapse = "','"
          ),
          "')"
        ),
        ",palette = '",
        query_group_palette,
        "',palcolor=c(",
        paste0("'", paste0(query_group_palcolor, collapse = "','"), "'"),
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
      if (isTRUE(query_collapsing)) {
        anno[[paste0(c("Query", query_group), collapse = ":")]] <- ComplexHeatmap::anno_block(
          align_to = split(
            seq_along(levels(cell_groups[["query_group"]])),
            levels(cell_groups[["query_group"]])
          ),
          panel_fun = methods::getFunction("panel_fun", where = environment()),
          which = ifelse(flip, "column", "row"),
          show_name = FALSE
        )
      } else {
        anno[[paste0(c("Query", query_group), collapse = ":")]] <- ComplexHeatmap::anno_block(
          align_to = split(
            seq_along(cell_groups[["query_group"]]),
            cell_groups[["query_group"]]
          ),
          panel_fun = methods::getFunction("panel_fun", where = environment()),
          which = ifelse(flip, "column", "row"),
          show_name = FALSE
        )
      }
      ha_cell_group <- do.call(
        ComplexHeatmap::HeatmapAnnotation,
        args = c(
          anno,
          which = ifelse(flip, "column", "row"),
          show_annotation_name = TRUE,
          annotation_name_side = ifelse(flip, "left", "bottom"),
          border = TRUE
        )
      )
      ha_query_list[[paste0(
        c("Query", query_group),
        collapse = ":"
      )]] <- ha_cell_group
      lgd[[paste0(c("Query", query_group), collapse = ":")]] <- ComplexHeatmap::Legend(
        title = paste0(c("Query", query_group), collapse = ":"),
        labels = levels(srt_query[[query_group, drop = TRUE]]),
        legend_gp = grid::gpar(
          fill = palette_colors(
            levels(srt_query[[query_group, drop = TRUE]]),
            palette = query_group_palette,
            palcolor = query_group_palcolor
          )
        ),
        border = TRUE
      )
    }
  }

  ha_ref_list <- NULL
  if (ref_group != "All.groups") {
    if (
      isFALSE(ref_collapsing) &&
        ((isFALSE(flip) & isTRUE(cluster_columns)) ||
          (isTRUE(flip) & isTRUE(cluster_rows)))
    ) {
      ref_cell_annotation <- c(ref_group, ref_cell_annotation)
      ref_cell_annotation_palette <- c(
        ref_group_palette,
        ref_cell_annotation_palette
      )
      ref_cell_annotation_palcolor <- c(
        list(ref_group_palcolor),
        ref_cell_annotation_palcolor %||% list(NULL)
      )
    } else {
      funbody <- paste0(
        "
        grid::grid.rect(gp = grid::gpar(fill = palette_colors(",
        paste0(
          "c('",
          paste0(levels(srt_ref[[ref_group, drop = TRUE]]), collapse = "','"),
          "')"
        ),
        ",palette = '",
        ref_group_palette,
        "',palcolor=c(",
        paste0("'", paste0(ref_group_palcolor, collapse = "','"), "'"),
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
      if (isTRUE(ref_collapsing)) {
        anno[[paste0(c("Ref", ref_group), collapse = ":")]] <- ComplexHeatmap::anno_block(
          align_to = split(
            seq_along(levels(cell_groups[["ref_group"]])),
            levels(cell_groups[["ref_group"]])
          ),
          panel_fun = methods::getFunction("panel_fun", where = environment()),
          which = ifelse(!flip, "column", "row"),
          show_name = FALSE
        )
      } else {
        anno[[paste0(c("Ref", ref_group), collapse = ":")]] <- ComplexHeatmap::anno_block(
          align_to = split(
            seq_along(cell_groups[["ref_group"]]),
            cell_groups[["ref_group"]]
          ),
          panel_fun = methods::getFunction("panel_fun", where = environment()),
          which = ifelse(!flip, "column", "row"),
          show_name = FALSE
        )
      }

      ha_cell_group <- do.call(
        ComplexHeatmap::HeatmapAnnotation,
        args = c(
          anno,
          which = ifelse(!flip, "column", "row"),
          show_annotation_name = TRUE,
          annotation_name_side = ifelse(!flip, "left", "bottom"),
          border = TRUE
        )
      )
      ha_ref_list[[paste0(
        c("Ref", ref_group),
        collapse = ":"
      )]] <- ha_cell_group
      lgd[[paste0(c("Ref", ref_group), collapse = ":")]] <- ComplexHeatmap::Legend(
        title = paste0(c("Ref", ref_group), collapse = ":"),
        labels = levels(srt_ref[[ref_group, drop = TRUE]]),
        legend_gp = grid::gpar(
          fill = palette_colors(
            levels(srt_ref[[ref_group, drop = TRUE]]),
            palette = ref_group_palette,
            palcolor = ref_group_palcolor
          )
        ),
        border = TRUE
      )
    }
  }

  if (!is.null(query_cell_annotation)) {
    query_subplots_list <- list()
    for (i in seq_along(query_cell_annotation)) {
      cellan <- query_cell_annotation[i]
      palette <- query_cell_annotation_palette[i]
      palcolor <- query_cell_annotation_palcolor[[i]]
      cell_anno <- cell_metadata[, paste0("query_", cellan)]
      names(cell_anno) <- rownames(cell_metadata)
      if (!is.numeric(cell_anno)) {
        if (is.logical(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = c(TRUE, FALSE))
        } else if (!is.factor(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = unique(cell_anno))
        }
        if (isTRUE(query_collapsing)) {
          subplots <- CellStatPlot(
            srt_query,
            flip = !flip,
            cells = gsub("query_", "", names(cell_groups[["query_group"]])),
            plot_type = "pie",
            group.by = query_group,
            stat.by = cellan,
            palette = palette,
            palcolor = palcolor,
            individual = TRUE,
            combine = FALSE
          )
          query_subplots_list[[paste0(cellan, ":", query_group)]] <- subplots
          graphics <- list()
          for (nm in names(subplots)) {
            funbody <- paste0(
              "
              g <- as_grob(query_subplots_list[['",
              cellan,
              ":",
              query_group,
              "']]",
              "[['",
              nm,
              "']] + theme_void() + theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = 'none'));
              g$name <- '",
              paste0(cellan, ":", query_group, "-", nm),
              "';
              grid::grid.draw(g)
              "
            )
            funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
            eval(
              parse(
                text = paste(
                  "graphics[[nm]] <- function(x, y, w, h) {",
                  funbody,
                  "}",
                  sep = ""
                )
              ),
              envir = environment()
            )
          }
          x_nm <- sapply(
            strsplit(levels(cell_groups[["query_group"]]), " : "),
            function(x) {
              if (length(x) == 2) {
                paste0(c(query_group, x[1], x[2]), collapse = ":")
              } else {
                paste0(c(query_group, x[1], ""), collapse = ":")
              }
            }
          )

          ha_cell <- list()
          ha_cell[[cellan]] <- ComplexHeatmap::anno_customize(
            x = x_nm,
            graphics = graphics,
            which = ifelse(flip, "column", "row"),
            border = TRUE,
            verbose = FALSE
          )
          anno_args <- c(
            ha_cell,
            which = ifelse(flip, "column", "row"),
            show_annotation_name = TRUE,
            annotation_name_side = ifelse(flip, "left", "bottom")
          )
          anno_args <- c(
            anno_args,
            query_cell_annotation_params[setdiff(
              names(query_cell_annotation_params),
              names(anno_args)
            )]
          )
          ha_query <- do.call(ComplexHeatmap::HeatmapAnnotation, args = anno_args)
          if (
            is.null(ha_query_list[[paste0(
              c("Query", query_group),
              collapse = ":"
            )]])
          ) {
            ha_query_list[[paste0(
              c("Query", query_group),
              collapse = ":"
            )]] <- ha_query
          } else {
            ha_query_list[[paste0(
              c("Query", query_group),
              collapse = ":"
            )]] <- c(
              ha_query_list[[paste0(c("Query", query_group), collapse = ":")]],
              ha_query
            )
          }
        } else {
          ha_cell <- list()
          ha_cell[[cellan]] <- ComplexHeatmap::anno_simple(
            x = as.character(cell_anno[paste0(
              "query_",
              names(cell_groups[["query_group"]])
            )]),
            col = palette_colors(
              cell_anno,
              palette = palette,
              palcolor = palcolor
            ),
            which = ifelse(flip, "column", "row"),
            na_col = "transparent",
            border = TRUE
          )
          anno_args <- c(
            ha_cell,
            which = ifelse(flip, "column", "row"),
            show_annotation_name = TRUE,
            annotation_name_side = ifelse(flip, "left", "bottom")
          )
          anno_args <- c(
            anno_args,
            query_cell_annotation_params[setdiff(
              names(query_cell_annotation_params),
              names(anno_args)
            )]
          )
          ha_query <- do.call(ComplexHeatmap::HeatmapAnnotation, args = anno_args)
          if (
            is.null(ha_query_list[[paste0(
              c("Query", query_group),
              collapse = ":"
            )]])
          ) {
            ha_query_list[[paste0(
              c("Query", query_group),
              collapse = ":"
            )]] <- ha_query
          } else {
            ha_query_list[[paste0(
              c("Query", query_group),
              collapse = ":"
            )]] <- c(
              ha_query_list[[paste0(c("Query", query_group), collapse = ":")]],
              ha_query
            )
          }
        }
        lgd[[paste0(c("Query", cellan), collapse = ":")]] <- ComplexHeatmap::Legend(
          title = paste0(c("Query", cellan), collapse = ":"),
          labels = levels(cell_anno),
          legend_gp = grid::gpar(
            fill = palette_colors(
              cell_anno,
              palette = palette,
              palcolor = palcolor
            )
          ),
          border = TRUE
        )
      } else {
        if (isTRUE(query_collapsing)) {
          subplots <- FeatureStatPlot(
            srt_query,
            assay = query_assay,
            layer = "data",
            flip = !flip,
            stat.by = cellan,
            cells = gsub("query_", "", names(cell_groups[["query_group"]])),
            group.by = query_group,
            palette = query_group_palette,
            palcolor = query_group_palcolor,
            fill.by = "group",
            same.y.lims = TRUE,
            individual = TRUE,
            combine = FALSE
          )
          query_subplots_list[[paste0(cellan, ":", query_group)]] <- subplots
          graphics <- list()
          for (nm in names(subplots)) {
            funbody <- paste0(
              "
              g <- as_grob(query_subplots_list[['",
              cellan,
              ":",
              query_group,
              "']]",
              "[['",
              nm,
              "']]  + facet_null() + theme_void() + theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = 'none'));
              g$name <- '",
              paste0(cellan, ":", query_group, "-", nm),
              "';
              grid::grid.draw(g)
              "
            )
            funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
            eval(
              parse(
                text = paste(
                  "graphics[[nm]] <- function(x, y, w, h) {",
                  funbody,
                  "}",
                  sep = ""
                )
              ),
              envir = environment()
            )
          }
          x_nm <- sapply(
            strsplit(levels(cell_groups[["query_group"]]), " : "),
            function(x) {
              if (length(x) == 2) {
                paste0(c(cellan, query_group, x[1], x[2]), collapse = ":")
              } else {
                paste0(c(cellan, query_group, x[1], ""), collapse = ":")
              }
            }
          )
          ha_cell <- list()
          ha_cell[[cellan]] <- ComplexHeatmap::anno_customize(
            x = x_nm,
            graphics = graphics,
            which = ifelse(flip, "column", "row"),
            border = TRUE,
            verbose = FALSE
          )
          anno_args <- c(
            ha_cell,
            which = ifelse(flip, "column", "row"),
            show_annotation_name = TRUE,
            annotation_name_side = ifelse(flip, "left", "bottom")
          )
          anno_args <- c(
            anno_args,
            query_cell_annotation_params[setdiff(
              names(query_cell_annotation_params),
              names(anno_args)
            )]
          )
          ha_query <- do.call(ComplexHeatmap::HeatmapAnnotation, args = anno_args)
          if (
            is.null(ha_query_list[[paste0(
              c("Query", query_group),
              collapse = ":"
            )]])
          ) {
            ha_query_list[[paste0(
              c("Query", query_group),
              collapse = ":"
            )]] <- ha_query
          } else {
            ha_query_list[[paste0(
              c("Query", query_group),
              collapse = ":"
            )]] <- c(
              ha_query_list[[paste0(c("Query", query_group), collapse = ":")]],
              ha_query
            )
          }
        } else {
          col_fun <- circlize::colorRamp2(
            breaks = seq(
              min(cell_anno, na.rm = TRUE),
              max(cell_anno, na.rm = TRUE),
              length = 100
            ),
            colors = palette_colors(palette = palette, palcolor = palcolor)
          )
          ha_cell <- list()
          ha_cell[[cellan]] <- ComplexHeatmap::anno_simple(
            x = cell_anno[paste0(
              "query_",
              names(cell_groups[["query_group"]])
            )],
            col = col_fun,
            which = ifelse(flip, "column", "row"),
            na_col = "transparent",
            border = TRUE
          )
          anno_args <- c(
            ha_cell,
            which = ifelse(flip, "column", "row"),
            show_annotation_name = TRUE,
            annotation_name_side = ifelse(flip, "left", "bottom")
          )
          anno_args <- c(
            anno_args,
            query_cell_annotation_params[setdiff(
              names(query_cell_annotation_params),
              names(anno_args)
            )]
          )
          ha_query <- do.call(ComplexHeatmap::HeatmapAnnotation, args = anno_args)
          if (
            is.null(ha_query_list[[paste0(
              c("Query", query_group),
              collapse = ":"
            )]])
          ) {
            ha_query_list[[paste0(
              c("Query", query_group),
              collapse = ":"
            )]] <- ha_query
          } else {
            ha_query_list[[paste0(
              c("Query", query_group),
              collapse = ":"
            )]] <- c(
              ha_query_list[[paste0(c("Query", query_group), collapse = ":")]],
              ha_query
            )
          }
          lgd[[paste0(c("Query", cellan), collapse = ":")]] <- ComplexHeatmap::Legend(
            title = paste0(c("Query", cellan), collapse = ":"),
            col_fun = col_fun,
            border = TRUE
          )
        }
      }
    }
  }

  if (!is.null(ref_cell_annotation)) {
    ref_subplots_list <- list()
    for (i in seq_along(ref_cell_annotation)) {
      cellan <- ref_cell_annotation[i]
      palette <- ref_cell_annotation_palette[i]
      palcolor <- ref_cell_annotation_palcolor[[i]]
      cell_anno <- cell_metadata[, paste0("ref_", cellan)]
      names(cell_anno) <- rownames(cell_metadata)
      if (!is.numeric(cell_anno)) {
        if (is.logical(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = c(TRUE, FALSE))
        } else if (!is.factor(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = unique(cell_anno))
        }
        if (isTRUE(ref_collapsing)) {
          subplots <- CellStatPlot(
            srt_ref,
            flip = flip,
            cells = gsub("ref_", "", names(cell_groups[["ref_group"]])),
            plot_type = "pie",
            group.by = ref_group,
            stat.by = cellan,
            palette = palette,
            palcolor = palcolor,
            individual = TRUE,
            combine = FALSE
          )
          ref_subplots_list[[paste0(cellan, ":", ref_group)]] <- subplots
          graphics <- list()
          for (nm in names(subplots)) {
            funbody <- paste0(
              "
              g <- as_grob(ref_subplots_list[['",
              cellan,
              ":",
              ref_group,
              "']]",
              "[['",
              nm,
              "']] + theme_void() + theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = 'none'));
              g$name <- '",
              paste0(cellan, ":", ref_group, "-", nm),
              "';
              grid::grid.draw(g)
              "
            )
            funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
            eval(
              parse(
                text = paste(
                  "graphics[[nm]] <- function(x, y, w, h) {",
                  funbody,
                  "}",
                  sep = ""
                )
              ),
              envir = environment()
            )
          }
          x_nm <- sapply(
            strsplit(levels(cell_groups[["ref_group"]]), " : "),
            function(x) {
              if (length(x) == 2) {
                paste0(c(ref_group, x[1], x[2]), collapse = ":")
              } else {
                paste0(c(ref_group, x[1], ""), collapse = ":")
              }
            }
          )

          ha_cell <- list()
          ha_cell[[cellan]] <- ComplexHeatmap::anno_customize(
            x = x_nm,
            graphics = graphics,
            which = ifelse(!flip, "column", "row"),
            border = TRUE,
            verbose = FALSE
          )
          anno_args <- c(
            ha_cell,
            which = ifelse(!flip, "column", "row"),
            show_annotation_name = TRUE,
            annotation_name_side = ifelse(!flip, "left", "bottom")
          )
          anno_args <- c(
            anno_args,
            ref_cell_annotation_params[setdiff(
              names(ref_cell_annotation_params),
              names(anno_args)
            )]
          )
          ha_ref <- do.call(ComplexHeatmap::HeatmapAnnotation, args = anno_args)
          if (
            is.null(ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]])
          ) {
            ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]] <- ha_ref
          } else {
            ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]] <- c(
              ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]],
              ha_ref
            )
          }
        } else {
          ha_cell <- list()
          ha_cell[[cellan]] <- ComplexHeatmap::anno_simple(
            x = as.character(cell_anno[paste0(
              "ref_",
              names(cell_groups[["ref_group"]])
            )]),
            col = palette_colors(
              cell_anno,
              palette = palette,
              palcolor = palcolor
            ),
            which = ifelse(!flip, "column", "row"),
            na_col = "transparent",
            border = TRUE
          )
          anno_args <- c(
            ha_cell,
            which = ifelse(!flip, "column", "row"),
            show_annotation_name = TRUE,
            annotation_name_side = ifelse(!flip, "left", "bottom")
          )
          anno_args <- c(
            anno_args,
            ref_cell_annotation_params[setdiff(
              names(ref_cell_annotation_params),
              names(anno_args)
            )]
          )
          ha_ref <- do.call(ComplexHeatmap::HeatmapAnnotation, args = anno_args)
          if (
            is.null(ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]])
          ) {
            ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]] <- ha_ref
          } else {
            ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]] <- c(
              ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]],
              ha_ref
            )
          }
        }
        lgd[[paste0(c("Ref", cellan), collapse = ":")]] <- ComplexHeatmap::Legend(
          title = paste0(c("Ref", cellan), collapse = ":"),
          labels = levels(cell_anno),
          legend_gp = grid::gpar(
            fill = palette_colors(
              cell_anno,
              palette = palette,
              palcolor = palcolor
            )
          ),
          border = TRUE
        )
      } else {
        if (isTRUE(ref_collapsing)) {
          subplots <- FeatureStatPlot(
            srt_ref,
            assay = ref_assay,
            layer = "data",
            flip = flip,
            stat.by = cellan,
            cells = gsub("ref_", "", names(cell_groups[["ref_group"]])),
            group.by = ref_group,
            palette = ref_group_palette,
            palcolor = ref_group_palcolor,
            fill.by = "group",
            same.y.lims = TRUE,
            individual = TRUE,
            combine = FALSE
          )
          ref_subplots_list[[paste0(cellan, ":", ref_group)]] <- subplots
          graphics <- list()
          for (nm in names(subplots)) {
            funbody <- paste0(
              "
              g <- as_grob(ref_subplots_list[['",
              cellan,
              ":",
              ref_group,
              "']]",
              "[['",
              nm,
              "']]  + facet_null() + theme_void() + theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = 'none'));
              g$name <- '",
              paste0(cellan, ":", ref_group, "-", nm),
              "';
              grid::grid.draw(g)
              "
            )
            funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
            eval(
              parse(
                text = paste(
                  "graphics[[nm]] <- function(x, y, w, h) {",
                  funbody,
                  "}",
                  sep = ""
                )
              ),
              envir = environment()
            )
          }
          x_nm <- sapply(
            strsplit(levels(cell_groups[["ref_group"]]), " : "),
            function(x) {
              if (length(x) == 2) {
                paste0(c(cellan, ref_group, x[1], x[2]), collapse = ":")
              } else {
                paste0(c(cellan, ref_group, x[1], ""), collapse = ":")
              }
            }
          )
          ha_cell <- list()
          ha_cell[[cellan]] <- ComplexHeatmap::anno_customize(
            x = x_nm,
            graphics = graphics,
            which = ifelse(!flip, "column", "row"),
            border = TRUE,
            verbose = FALSE
          )
          anno_args <- c(
            ha_cell,
            which = ifelse(!flip, "column", "row"),
            show_annotation_name = TRUE,
            annotation_name_side = ifelse(!flip, "left", "bottom")
          )
          anno_args <- c(
            anno_args,
            ref_cell_annotation_params[setdiff(
              names(ref_cell_annotation_params),
              names(anno_args)
            )]
          )
          ha_ref <- do.call(ComplexHeatmap::HeatmapAnnotation, args = anno_args)
          if (
            is.null(ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]])
          ) {
            ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]] <- ha_ref
          } else {
            ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]] <- c(
              ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]],
              ha_ref
            )
          }
        } else {
          col_fun <- circlize::colorRamp2(
            breaks = seq(
              min(cell_anno, na.rm = TRUE),
              max(cell_anno, na.rm = TRUE),
              length = 100
            ),
            colors = palette_colors(palette = palette, palcolor = palcolor)
          )
          ha_cell <- list()
          ha_cell[[cellan]] <- ComplexHeatmap::anno_simple(
            x = cell_anno[paste0("ref_", names(cell_groups[["ref_group"]]))],
            col = col_fun,
            which = ifelse(!flip, "column", "row"),
            na_col = "transparent",
            border = TRUE
          )
          anno_args <- c(
            ha_cell,
            which = ifelse(!flip, "column", "row"),
            show_annotation_name = TRUE,
            annotation_name_side = ifelse(!flip, "left", "bottom")
          )
          anno_args <- c(
            anno_args,
            ref_cell_annotation_params[setdiff(
              names(ref_cell_annotation_params),
              names(anno_args)
            )]
          )
          ha_ref <- do.call(ComplexHeatmap::HeatmapAnnotation, args = anno_args)
          if (
            is.null(ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]])
          ) {
            ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]] <- ha_ref
          } else {
            ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]] <- c(
              ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]],
              ha_ref
            )
          }
          lgd[[paste0(c("Ref", cellan), collapse = ":")]] <- ComplexHeatmap::Legend(
            title = paste0(c("Ref", cellan), collapse = ":"),
            col_fun = col_fun,
            border = TRUE
          )
        }
      }
    }
  }

  if (isFALSE(ref_legend)) {
    lgd <- lgd[grep("Ref", names(lgd), invert = TRUE)]
  }

  layer_fun <- function(j, i, x, y, w, h, fill) {
    if (nlabel > 0) {
      if (flip) {
        mat <- Matrix::t(simil_matrix)
      } else {
        mat <- simil_matrix
      }
      value <- ComplexHeatmap::pindex(mat, i, j)
      ind_mat <- ComplexHeatmap::restore_matrix(j, i, x, y)

      inds <- NULL
      if (label_by %in% c("row", "both")) {
        for (row in 1:nrow(ind_mat)) {
          ind <- ind_mat[row, ]
          ind <- ind[which(
            value[ind] >=
              max(
                c(sort(value[ind], decreasing = TRUE)[nlabel]),
                na.rm = TRUE
              ) &
              value[ind] >= label_cutoff
          )]
          inds <- c(inds, ind)
        }
      }
      if (label_by %in% c("column", "both")) {
        for (column in 1:ncol(ind_mat)) {
          ind <- ind_mat[, column]
          ind <- ind[which(
            value[ind] >=
              max(
                c(sort(value[ind], decreasing = TRUE)[nlabel]),
                na.rm = TRUE
              ) &
              value[ind] >= label_cutoff
          )]
          inds <- c(inds, ind)
        }
      }
      if (label_by == "both") {
        inds <- inds[duplicated(inds)]
      }
      if (length(inds) > 0) {
        theta <- seq(pi / 8, 2 * pi, length.out = 16)
        lapply(theta, function(i) {
          x_out <- x[inds] + grid::unit(cos(i) * label_size / 30, "mm")
          y_out <- y[inds] + grid::unit(sin(i) * label_size / 30, "mm")
          grid::grid.text(
            round(value[inds], 2),
            x = x_out,
            y = y_out,
            gp = grid::gpar(fontsize = label_size, col = "white")
          )
        })
        grid::grid.text(
          round(value[inds], 2),
          x[inds],
          y[inds],
          gp = grid::gpar(fontsize = label_size, col = "black")
        )
      }
    }
  }

  ht_list <- NULL
  ht_args <- list(
    name = exp_name,
    matrix = if (flip) Matrix::t(simil_matrix) else simil_matrix,
    col = colors,
    layer_fun = layer_fun,
    row_title = row_title %||%
      if (flip) {
        paste0(c("Ref", ref_group), collapse = ":")
      } else {
        paste0(c("Query", query_group), collapse = ":")
      },
    row_title_side = row_title_side,
    column_title = column_title %||%
      if (flip) {
        paste0(c("Query", query_group), collapse = ":")
      } else {
        paste0(c("Ref", ref_group), collapse = ":")
      },
    column_title_side = column_title_side,
    row_title_rot = row_title_rot,
    column_title_rot = column_title_rot,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    row_names_side = row_names_side,
    column_names_side = column_names_side,
    row_names_rot = row_names_rot,
    column_names_rot = column_names_rot,
    top_annotation = if (flip) ha_query_list[[1]] else ha_ref_list[[1]],
    left_annotation = if (flip) ha_ref_list[[1]] else ha_query_list[[1]],
    show_heatmap_legend = FALSE,
    border = border,
    use_raster = use_raster,
    raster_device = raster_device,
    raster_by_magick = raster_by_magick,
    width = if (is.numeric(width)) grid::unit(width, units = units) else NULL,
    height = if (is.numeric(height)) grid::unit(height, units = units) else NULL
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

  if (!is.null(width) || !is.null(height)) {
    fix <- TRUE
  } else {
    fix <- FALSE
  }
  rendersize <- heatmap_rendersize(
    width = width,
    height = height,
    units = units,
    ha_top_list = ha_ref_list,
    ha_left = ha_query_list[[1]],
    ha_right = NULL,
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
      features = features,
      simil_matrix = simil_matrix,
      simil_name = simil_name,
      cell_metadata = cell_metadata
    )
  )
}
