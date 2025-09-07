#' @title The Group Heatmap
#'
#' @md
#' @param srt A Seurat object.
#' @param features The features to include in the heatmap.
#' Default is NULL.
#' @param group.by A character vector specifying the groups to group by.
#' Default is NULL.
#' @param split.by A character vector specifying the variable to split the heatmap by.
#'  Default is NULL.
#' @param within_groups A logical value indicating whether to create separate heatmap scales for each group or within each group.
#' Default is FALSE.
#' @param grouping.var A character vector that specifies another variable for grouping, such as certain conditions.
#' Default is NULL.
#' @param numerator A character vector specifying the value to use as the numerator in the grouping.var grouping.
#' Default is NULL.
#' @param cells A character vector specifying the cells to include in the heatmap.
#' Default is NULL.
#' @param aggregate_fun A function to use for aggregating data within groups.
#' Default is [base::mean].
#' @param exp_cutoff A numeric value specifying the threshold for cell counting if \code{add_dot} is TRUE.
#' Default is 0.
#' @param border A logical value indicating whether to add a border to the heatmap.
#' Default is TRUE.
#' @param flip A logical value indicating whether to flip the heatmap.
#' Default is FALSE.
#' @param layer A character vector specifying the layer in the Seurat object to use.
#' Default is "counts".
#' @param assay A character vector specifying the assay in the Seurat object to use.
#' Default is NULL.
#' @param exp_method A character vector specifying the method for calculating expression values.
#' Default is "zscore" with options "zscore", "raw", "fc", "log2fc", "log1p".
#' @param exp_legend_title A character vector specifying the title for the legend of expression value.
#' Default is NULL.
#' @param limits A two-length numeric vector specifying the limits for the color scale.
#' Default is NULL.
#' @param lib_normalize A logical value indicating whether to normalize the data by library size.
#' @param libsize A numeric vector specifying the library size for each cell.
#' Default is NULL.
#' @param feature_split A factor specifying how to split the features.
#' Default is NULL.
#' @param feature_split_by A character vector specifying which group.by to use when splitting features (into n_split feature clusters).
#' Default is NULL.
#' @param n_split An integer specifying the number of feature splits (feature clusters) to create.
#' Default is NULL.
#' @param split_order A numeric vector specifying the order of splits. Default is NULL.
#' @param split_method A character vector specifying the method for splitting features.
#' Default is "kmeans" with options "kmeans", "hclust", "mfuzz").
#' @param decreasing A logical value indicating whether to sort feature splits in decreasing order.
#' Default is FALSE.
#' @param fuzzification A numeric value specifying the fuzzification coefficient.
#' Default is NULL.
#' @param cluster_features_by A character vector specifying which group.by to use when clustering features.
#' Default is NULL. By default, this parameter is set to NULL, which means that all groups will be used.
#' @param cluster_rows A logical value indicating whether to cluster rows in the heatmap.
#' Default is FALSE.
#' @param cluster_columns A logical value indicating whether to cluster columns in the heatmap.
#' Default is FALSE.
#' @param cluster_row_slices A logical value indicating whether to cluster row slices in the heatmap.
#' Default is FALSE.
#' @param cluster_column_slices A logical value indicating whether to cluster column slices in the heatmap.
#' Default is FALSE.
#' @param show_row_names A logical value indicating whether to show row names in the heatmap.
#' Default is FALSE.
#' @param show_column_names A logical value indicating whether to show column names in the heatmap.
#' Default is FALSE.
#' @param row_names_side A character vector specifying the side to place row names.
#' @param column_names_side A character vector specifying the side to place column names.
#' @param row_names_rot A numeric value specifying the rotation angle for row names.
#' Default is 0.
#' @param column_names_rot A numeric value specifying the rotation angle for column names.
#' Default is 90.
#' @param row_title A character vector specifying the title for rows.
#' Default is NULL.
#' @param column_title A character vector specifying the title for columns.
#' Default is NULL.
#' @param row_title_side A character vector specifying the side to place row title.
#' Default is "left".
#' @param column_title_side A character vector specifying the side to place column title.
#' Default is "top".
#' @param row_title_rot A numeric value specifying the rotation angle for row title.
#' Default is 0.
#' @param column_title_rot A numeric value specifying the rotation angle for column title.
#' @param anno_terms A logical value indicating whether to include term annotations.
#' Default is FALSE.
#' @param anno_keys A logical value indicating whether to include key annotations.
#' Default is FALSE.
#' @param anno_features A logical value indicating whether to include feature annotations.
#' Default is FALSE.
#' @param terms_width A unit specifying the width of term annotations.
#' Default is unit(4, "in").
#' @param terms_fontsize A numeric vector specifying the font size(s) for term annotations.
#' Default is 8.
#' @param keys_width A unit specifying the width of key annotations.
#' Default is unit(2, "in").
#' @param keys_fontsize A two-length numeric vector specifying the minimum and maximum font size(s) for key annotations.
#' Default is c(6, 10).
#' @param features_width A unit specifying the width of feature annotations.
#' Default is unit(2, "in").
#' @param features_fontsize A two-length numeric vector specifying the minimum and maximum font size(s) for feature annotations.
#' Default is c(6, 10).
#' @param IDtype A character vector specifying the type of IDs for features.
#' Default is "symbol".
#' @param species A character vector specifying the species for features.
#' Default is "Homo_sapiens".
#' @param db_update A logical value indicating whether to update the database.
#' Default is FALSE.
#' @param db_version A character vector specifying the version of the database.
#' Default is "latest".
#' @param db_combine A logical value indicating whether to use a combined database.
#' Default is FALSE.
#' @param convert_species A logical value indicating whether to use a species-converted database if annotation is missing for `species`.
#' Default is FALSE.
#' @param Ensembl_version An integer specifying the Ensembl version.
#' Default is 103.
#' @param mirror A character vector specifying the mirror for the Ensembl database.
#' Default is NULL.
#' @param db A character vector specifying the database to use.
#' Default is "GO_BP".
#' @param TERM2GENE A data.frame specifying the TERM2GENE mapping for the database.
#' Default is NULL.
#' @param TERM2NAME A data.frame specifying the TERM2NAME mapping for the database.
#' Default is NULL.
#' @param minGSSize An integer specifying the minimum gene set size for the database.
#' Default is 10.
#' @param maxGSSize An integer specifying the maximum gene set size for the database.
#' Default is 500.
#' @param GO_simplify A logical value indicating whether to simplify gene ontology terms.
#' Default is FALSE.
#' @param GO_simplify_cutoff A character vector specifying the cutoff for GO simplification.
#' Default is "p.adjust < 0.05".
#' @param simplify_method A character vector specifying the method for GO simplification.
#' Default is "Wang".
#' @param simplify_similarityCutoff A numeric value specifying the similarity cutoff for GO simplification.
#' Default is 0.7.
#' @param pvalueCutoff A numeric vector specifying the p-value cutoff(s) for significance.
#' Default is NULL.
#' @param padjustCutoff A numeric value specifying the adjusted p-value cutoff for significance.
#' Default is 0.05.
#' @param topTerm An integer specifying the number of top terms to include.
#' Default is 5.
#' @param show_termid A logical value indicating whether to show term IDs.
#' Default is FALSE.
#' @param topWord An integer specifying the number of top words to include.
#' Default is 20.
#' @param words_excluded A character vector specifying the words to exclude.
#' Default is NULL.
#' @param nlabel An integer specifying the number of labels to include.
#' Default is 0.
#' @param features_label A character vector specifying the features to label.
#' Default is NULL.
#' @param label_size A numeric value specifying the size of labels.
#' Default is 10.
#' @param label_color A character vector specifying the color of labels.
#' Default is "black".
#' @param add_bg A logical value indicating whether to add a background to the heatmap.
#' Default is FALSE.
#' @param bg_alpha A numeric value specifying the alpha value for the background color.
#' Default is 0.5.
#' @param add_dot A logical value indicating whether to add dots to the heatmap.
#' The size of dot represents percentage of expressed cells based on the specified `exp_cutoff`.
#' Default is FALSE.
#' @param dot_size A unit specifying the base size of the dots.
#' Default is unit(8, "mm").
#' @param add_reticle A logical value indicating whether to add reticles to the heatmap.
#' Default is FALSE.
#' @param reticle_color A character vector specifying the color of the reticles.
#' Default is "grey".
#' @param add_violin A logical value indicating whether to add violins to the heatmap.
#' Default is FALSE.
#' @param fill.by A character vector specifying what to fill the violin.
#' Possible values are "group", "feature", or "expression".
#' Default is "feature".
#' @param fill_palette A character vector specifying the palette to use for fill.
#' Default is "Dark2".
#' @param fill_palcolor A character vector specifying the fill color to use.
#' Default is NULL.
#' @param heatmap_palette A character vector specifying the palette to use for the heatmap.
#' Default is "RdBu".
#' @param heatmap_palcolor A character vector specifying the heatmap color to use.
#' Default is NULL.
#' @param group_palette A character vector specifying the palette to use for groups.
#' Default is "Paired".
#' @param group_palcolor A character vector specifying the group color to use.
#' Default is NULL.
#' @param cell_split_palette A character vector specifying the palette to use for cell splits.
#' Default is "simspec".
#' @param cell_split_palcolor A character vector specifying the cell split color to use.
#' Default is NULL.
#' @param feature_split_palette A character vector specifying the palette to use for feature splits.
#' Default is "simspec".
#' @param feature_split_palcolor A character vector specifying the feature split color to use.
#' Default is NULL.
#' @param cell_annotation A character vector specifying the cell annotation(s) to include.
#' Default is NULL.
#' @param cell_annotation_palette A character vector specifying the palette to use for cell annotations.
#' The length of the vector should match the number of cell_annotation. Default is "Paired".
#' @param cell_annotation_palcolor A list of character vector specifying the cell annotation color(s) to use.
#' The length of the list should match the number of cell_annotation. Default is NULL.
#' @param cell_annotation_params A list specifying additional parameters for cell annotations.
#' Default is a list with width = unit(1, "cm") if flip is TRUE, else a list with height = unit(1, "cm").
#' @param feature_annotation A character vector specifying the feature annotation(s) to include.
#' Default is NULL.
#' @param feature_annotation_palette A character vector specifying the palette to use for feature annotations.
#' The length of the vector should match the number of feature_annotation.
#' Default is "Dark2".
#' @param feature_annotation_palcolor A list of character vector specifying the feature annotation color to use.
#' The length of the list should match the number of feature_annotation.
#' Default is NULL.
#' @param feature_annotation_params A list specifying additional parameters for feature annotations.
#' Default is an empty list.
#' @param use_raster A logical value indicating whether to use a raster device for plotting.
#' Default is NULL.
#' @param raster_device A character vector specifying the raster device to use.
#' Default is "png".
#' @param raster_by_magick A logical value indicating whether to use the 'magick' package for raster.
#' Default is FALSE.
#' @param height A numeric vector specifying the height(s) of the heatmap body.
#' Default is NULL.
#' @param width A numeric vector specifying the width(s) of the heatmap body.
#' Default is NULL.
#' @param units A character vector specifying the units for the height and width.
#' Default is "inch".
#' @param seed An integer specifying the random seed. Default is 11.
#' @param ht_params A list specifying additional parameters passed to the [ComplexHeatmap::Heatmap] function.
#' Default is an empty list.
#' @param verbose A logical value indicating whether to print messages.
#' Default is TRUE.
#' @param ... Additional arguments passed to the [ComplexHeatmap::Heatmap] function.
#'
#' @seealso [RunDEtest]
#'
#' @return
#' A list with the following elements:
#' \itemize{
#'   \item \code{plot:} The heatmap plot.
#'   \item \code{matrix_list:} A list of matrix for each `group.by` used in the heatmap.
#'   \item \code{feature_split:} NULL or a factor if splitting is performed in the heatmap.
#'   \item \code{cell_metadata:} Meta data of cells used to generate the heatmap.
#'   \item \code{feature_metadata:} Meta data of features used to generate the heatmap.
#'   \item \code{enrichment:} NULL or a enrichment result generated by [RunEnrichment] when any of the parameters `anno_terms`, `anno_keys`, or `anno_features` is set to TRUE.
#' }
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' data(pancreas_sub)
#'
#' pancreas_sub <- AnnotateFeatures(
#'   pancreas_sub,
#'   species = "Mus_musculus",
#'   db = c("CSPA", "TF")
#' )
#'
#' ht1 <- GroupHeatmap(
#'   pancreas_sub,
#'   features = c(
#'     "Sox9", "Anxa2", "Bicc1", # Ductal
#'     "Neurog3", "Hes6", # EPs
#'     "Fev", "Neurod1", # Pre-endocrine
#'     "Rbp4", "Pyy", # Endocrine
#'     "Ins1", "Gcg", "Sst", "Ghrl"
#'     # Beta, Alpha, Delta, Epsilon
#'   ),
#'   group.by = c("CellType", "SubCellType")
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
#' # pancreas_sub <- RunDEtest(
#' #  pancreas_sub,
#' #   group_by = "CellType"
#' # )
#' de_filter <- filter(
#'   pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox,
#'   p_val_adj < 0.05 & avg_log2FC > 1
#' )
#'
#' ht2 <- GroupHeatmap(
#'   pancreas_sub,
#'   features = de_filter$gene,
#'   group.by = "CellType",
#'   split.by = "Phase",
#'   cell_split_palette = "Dark2",
#'   cluster_rows = TRUE,
#'   cluster_columns = TRUE
#' )
#' ht2$plot
#'
#' ht3 <- GroupHeatmap(
#'   pancreas_sub,
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
#' de_top <- de_filter %>%
#'   group_by(gene) %>%
#'   top_n(1, avg_log2FC) %>%
#'   group_by(group1) %>%
#'   top_n(3, avg_log2FC)
#' ht4 <- GroupHeatmap(
#'   pancreas_sub,
#'   features = de_top$gene,
#'   feature_split = de_top$group1,
#'   group.by = "CellType",
#'   heatmap_palette = "YlOrRd",
#'   cell_annotation = c(
#'     "Phase", "G2M_score", "Neurod2"
#'   ),
#'   cell_annotation_palette = c(
#'     "Dark2", "Paired", "Paired"
#'   ),
#'   cell_annotation_params = list(
#'     height = grid::unit(10, "mm")
#'   ),
#'   feature_annotation = c("TF", "CSPA"),
#'   feature_annotation_palcolor = list(
#'     c("gold", "steelblue"),
#'     c("forestgreen")
#'   ),
#'   add_dot = TRUE,
#'   add_bg = TRUE,
#'   nlabel = 0,
#'   show_row_names = TRUE
#' )
#' ht4$plot
#'
#' ht5 <- GroupHeatmap(
#'   pancreas_sub,
#'   features = de_top$gene,
#'   feature_split = de_top$group1,
#'   group.by = "CellType",
#'   heatmap_palette = "YlOrRd",
#'   cell_annotation = c(
#'     "Phase", "G2M_score", "Neurod2"
#'   ),
#'   cell_annotation_palette = c(
#'     "Dark2", "Paired", "Paired"
#'   ),
#'   cell_annotation_params = list(
#'     width = grid::unit(10, "mm")
#'   ),
#'   feature_annotation = c("TF", "CSPA"),
#'   feature_annotation_palcolor = list(
#'     c("gold", "steelblue"), c("forestgreen")
#'   ),
#'   add_dot = TRUE,
#'   add_bg = TRUE,
#'   flip = TRUE,
#'   column_title_rot = 45,
#'   nlabel = 0,
#'   show_row_names = TRUE
#' )
#' ht5$plot
#'
#' ht6 <- GroupHeatmap(
#'   pancreas_sub,
#'   features = de_top$gene,
#'   feature_split = de_top$group1,
#'   group.by = "CellType",
#'   add_violin = TRUE,
#'   cluster_rows = TRUE,
#'   nlabel = 0,
#'   show_row_names = TRUE
#' )
#' ht6$plot
#'
#' ht7 <- GroupHeatmap(
#'   pancreas_sub,
#'   features = de_top$gene,
#'   feature_split = de_top$group1,
#'   group.by = "CellType",
#'   add_violin = TRUE,
#'   fill.by = "expression",
#'   fill_palette = "Blues",
#'   cluster_rows = TRUE,
#'   nlabel = 0,
#'   show_row_names = TRUE
#' )
#' ht7$plot
#'
#' ht8 <- GroupHeatmap(
#'   pancreas_sub,
#'   features = de_top$gene,
#'   group.by = "CellType",
#'   split.by = "Phase",
#'   n_split = 4,
#'   cluster_rows = TRUE,
#'   cluster_columns = TRUE,
#'   cluster_row_slices = TRUE,
#'   cluster_column_slices = TRUE,
#'   add_dot = TRUE,
#'   add_reticle = TRUE,
#'   heatmap_palette = "viridis",
#'   nlabel = 0,
#'   show_row_names = TRUE,
#'   ht_params = list(
#'     row_gap = grid::unit(0, "mm"),
#'     row_names_gp = grid::gpar(fontsize = 10)
#'   )
#' )
#' ht8$plot
GroupHeatmap <- function(
    srt,
    features = NULL,
    group.by = NULL,
    split.by = NULL,
    within_groups = FALSE,
    grouping.var = NULL,
    numerator = NULL,
    cells = NULL,
    aggregate_fun = base::mean,
    exp_cutoff = 0,
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
    add_bg = FALSE,
    bg_alpha = 0.5,
    add_dot = FALSE,
    dot_size = grid::unit(8, "mm"),
    add_reticle = FALSE,
    reticle_color = "grey",
    add_violin = FALSE,
    fill.by = "feature",
    fill_palette = "Dark2",
    fill_palcolor = NULL,
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
      list(width = grid::unit(10, "mm"))
    } else {
      list(height = grid::unit(10, "mm"))
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
    ht_params = list(),
    verbose = TRUE,
    ...) {
  set.seed(seed)

  if (isTRUE(raster_by_magick)) {
    check_r("magick")
  }
  if (is.null(features)) {
    log_message(
      "{.arg features} is not provided.",
      message_type = "error"
    )
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

  if (!is.null(grouping.var) && exp_method != "log2fc") {
    log_message(
      "When 'grouping.var' is specified, 'exp_method' can only be 'log2fc'",
      message_type = "warning",
      verbose = verbose
    )
    exp_method <- "log2fc"
  }
  exp_name <- exp_legend_title %||% exp_name

  if (!is.null(grouping.var)) {
    if (identical(split.by, grouping.var)) {
      log_message(
        "'grouping.var' must be different from 'split.by'",
        message_type = "error"
      )
    }
    if (!is.factor(srt@meta.data[[grouping.var]])) {
      srt@meta.data[[grouping.var]] <- factor(
        srt@meta.data[[grouping.var]],
        levels = unique(srt@meta.data[[grouping.var]])
      )
    }
    if (is.null(numerator)) {
      numerator <- levels(srt@meta.data[[grouping.var]])[1]
      log_message(
        "'{.arg numerator}' is not specified. Use the first level in '{.arg grouping.var}': ",
        numerator,
        message_type = "warning",
        verbose = verbose
      )
    } else {
      if (!numerator %in% levels(srt@meta.data[, grouping.var])) {
        log_message(
          "'{.arg numerator}' is not an element of the '{.arg grouping.var}'",
          message_type = "error"
        )
      }
    }
    srt@meta.data[["grouping.var.use"]] <- srt@meta.data[[grouping.var]] ==
      numerator
    add_dot <- FALSE
    exp_name <- paste0(numerator, "/", "other\n", exp_method, "(", data_nm, ")")
  }

  assay <- assay %||% DefaultAssay(srt)
  assay_use <- Seurat::GetAssay(srt, assay = assay)

  if (is.null(group.by)) {
    srt@meta.data[["All.groups"]] <- factor("")
    group.by <- "All.groups"
  }

  if (any(!group.by %in% colnames(srt@meta.data))) {
    log_message(
      group.by[!group.by %in% colnames(srt@meta.data)],
      " is not in the meta data of the Seurat object.",
      message_type = "error"
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
    log_message(
      "'{.arg split.by}' only support one variable.",
      message_type = "error"
    )
  }
  if (any(!split.by %in% colnames(srt@meta.data))) {
    log_message(
      split.by[!split.by %in% colnames(srt@meta.data)],
      " is not in the meta data of the Seurat object.",
      message_type = "error"
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
  group_elements <- unlist(lapply(
    srt@meta.data[, group.by, drop = FALSE],
    function(x) length(unique(x))
  ))

  if (any(group_elements == 1) && exp_method == "zscore") {
    log_message(
      "{.arg exp_method} {.val zscore} cannot be applied to the group(s) consisting of one element: ",
      paste0(names(group_elements)[group_elements == 1], collapse = ","),
      message_type = "error"
    )
  }
  if (length(group_palette) == 1) {
    group_palette <- rep(group_palette, length(group.by))
  }
  if (length(group_palette) != length(group.by)) {
    log_message(
      "{.arg group_palette} must be the same length as {.arg group.by}",
      message_type = "error"
    )
  }
  group_palette <- stats::setNames(group_palette, nm = group.by)
  raw_group_by <- group.by
  raw_group_palette <- group_palette
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

  if (!is.null(feature_split) && !is.factor(feature_split)) {
    feature_split <- factor(feature_split, levels = unique(feature_split))
  }
  if (length(feature_split) != 0 && length(feature_split) != length(features)) {
    log_message(
      "{.arg feature_split} must be the same length as {.arg features}",
      message_type = "error"
    )
  }
  if (is.null(feature_split_by)) {
    feature_split_by <- group.by
  }
  if (any(!feature_split_by %in% group.by)) {
    log_message(
      "{.arg feature_split_by} must be a subset of {.arg group.by}",
      message_type = "error"
    )
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
        "{.arg cell_annotation_palette} and {.arg cell_annotation_palcolor} must be the same length as {.arg cell_annotation}",
        message_type = "error"
      )
    }
    if (
      any(
        !cell_annotation %in%
          c(colnames(srt@meta.data), rownames(assay_use))
      )
    ) {
      log_message(
        "cell_annotation: ",
        paste0(
          cell_annotation[
            !cell_annotation %in%
              c(colnames(srt@meta.data), rownames(assay_use))
          ],
          collapse = ","
        ),
        " is not in the Seurat object.",
        message_type = "error"
      )
    }
  }
  feature_annotation_data <- GetFeaturesData(
    srt,
    assay = assay
  )
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
        "{.arg feature_annotation_palette} and {.arg feature_annotation_palcolor} must be the same length as {.arg feature_annotation}",
        message_type = "error"
      )
    }
    if (any(!feature_annotation %in% colnames(feature_annotation_data))) {
      log_message(
        "{.arg feature_annotation}: ",
        paste0(
          feature_annotation[
            !feature_annotation %in% colnames(feature_annotation_data)
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
    cells <- colnames(assay_use)
  }
  if (all(!cells %in% colnames(assay_use))) {
    log_message(
      "No cells found",
      message_type = "error"
    )
  }
  if (!all(cells %in% colnames(assay_use))) {
    log_message(
      "Some cells not found",
      message_type = "warning",
      verbose = verbose
    )
  }
  cells <- intersect(cells, colnames(assay_use))

  if (is.null(features)) {
    features <- SeuratObject::VariableFeatures(
      srt,
      assay = assay
    )
  }
  index <- features %in% c(rownames(assay_use), colnames(srt@meta.data))
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
    cell_groups[[cell_group]] <- stats::setNames(
      srt@meta.data[cells, cell_group],
      cells
    )
    cell_groups[[cell_group]] <- stats::na.omit(cell_groups[[cell_group]])
    cell_groups[[cell_group]] <- factor(
      cell_groups[[cell_group]],
      levels = levels(cell_groups[[cell_group]])[
        levels(cell_groups[[cell_group]]) %in% cell_groups[[cell_group]]
      ]
    )

    if (!is.null(split.by)) {
      if (!is.factor(srt@meta.data[[split.by]])) {
        srt@meta.data[[split.by]] <- factor(
          srt@meta.data[[split.by]],
          levels = unique(srt@meta.data[[split.by]])
        )
      }
      levels <- apply(
        expand.grid(
          levels(srt@meta.data[[split.by]]),
          levels(cell_groups[[cell_group]])
        ),
        1,
        function(x) paste0(x[2:1], collapse = " : ")
      )
      cell_groups[[cell_group]] <- stats::setNames(
        paste0(
          cell_groups[[cell_group]][cells],
          " : ",
          srt@meta.data[cells, split.by]
        ),
        cells
      )
      cell_groups[[cell_group]] <- factor(
        cell_groups[[cell_group]],
        levels = levels[levels %in% cell_groups[[cell_group]]]
      )
    }

    if (!is.null(grouping.var)) {
      levels <- apply(
        expand.grid(c("TRUE", "FALSE"), levels(cell_groups[[cell_group]])),
        1,
        function(x) paste0(x[2:1], collapse = " ; ")
      )
      cell_groups[[cell_group]] <- stats::setNames(
        paste0(
          cell_groups[[cell_group]][cells],
          " ; ",
          srt@meta.data[cells, "grouping.var.use"]
        ),
        cells
      )
      cell_groups[[cell_group]] <- factor(
        cell_groups[[cell_group]],
        levels = levels[levels %in% cell_groups[[cell_group]]]
      )
    }
  }

  gene <- features[features %in% rownames(assay_use)]
  gene_unique <- features_unique[features %in% rownames(assay_use)]
  meta <- features[features %in% colnames(srt@meta.data)]

  mat_raw <- as_matrix(
    rbind(
      GetAssayData5(
        srt,
        assay = assay,
        layer = layer
      )[gene, cells, drop = FALSE],
      Matrix::t(srt@meta.data[cells, meta, drop = FALSE])
    )
  )[features, , drop = FALSE]
  rownames(mat_raw) <- features_unique
  if (isTRUE(lib_normalize) && min(mat_raw, na.rm = TRUE) >= 0) {
    if (!is.null(libsize)) {
      libsize_use <- libsize
    } else {
      libsize_use <- Matrix::colSums(
        GetAssayData5(
          srt,
          assay = assay,
          layer = "counts"
        )[, colnames(mat_raw), drop = FALSE]
      )
      isfloat <- any(libsize_use %% 1 != 0, na.rm = TRUE)
      if (isTRUE(isfloat)) {
        libsize_use <- rep(1, length(libsize_use))
        log_message(
          "The values in the {.val counts} layer are non-integer. Set the library size to {.val 1}",
          message_type = "warning",
          verbose = verbose
        )
        if (!is.null(grouping.var)) {
          exp_name <- paste0(
            numerator,
            "/",
            "other\n",
            exp_method,
            "(",
            layer,
            ")"
          )
        } else {
          exp_name <- paste0(exp_method, "(", layer, ")")
        }
      }
    }
    mat_raw[gene_unique, ] <- Matrix::t(
      Matrix::t(
        mat_raw[gene_unique, , drop = FALSE]
      ) /
        libsize_use *
        stats::median(libsize_use)
    )
  }

  mat_raw_list <- list()
  mat_perc_list <- list()
  for (cell_group in names(cell_groups)) {
    mat_tmp <- Matrix::t(
      stats::aggregate(
        Matrix::t(mat_raw[features_unique, , drop = FALSE]),
        by = list(cell_groups[[cell_group]][colnames(mat_raw)]),
        FUN = aggregate_fun
      )
    )
    colnames(mat_tmp) <- mat_tmp[1, , drop = FALSE]
    mat_tmp <- mat_tmp[-1, , drop = FALSE]
    class(mat_tmp) <- "numeric"
    mat_raw_list[[cell_group]] <- mat_tmp

    mat_perc <- Matrix::t(
      stats::aggregate(
        Matrix::t(mat_raw[features_unique, , drop = FALSE]),
        by = list(cell_groups[[cell_group]][colnames(mat_raw)]),
        FUN = function(x) {
          sum(x > exp_cutoff) / length(x)
        }
      )
    )
    colnames(mat_perc) <- mat_perc[1, , drop = FALSE]
    mat_perc <- mat_perc[-1, , drop = FALSE]
    class(mat_perc) <- "numeric"
    if (isTRUE(flip)) {
      mat_perc <- Matrix::t(mat_perc)
    }
    mat_perc_list[[cell_group]] <- mat_perc
  }

  # data used to plot heatmap
  mat_list <- list()
  for (cell_group in group.by) {
    mat_tmp <- mat_raw_list[[cell_group]]
    if (is.null(grouping.var)) {
      mat_tmp <- matrix_process(mat_tmp, method = exp_method)
      mat_tmp[is.infinite(mat_tmp)] <- max(
        abs(mat_tmp[!is.infinite(mat_tmp)]),
        na.rm = TRUE
      ) *
        ifelse(mat_tmp[is.infinite(mat_tmp)] > 0, 1, -1)
      mat_tmp[is.na(mat_tmp)] <- mean(mat_tmp, na.rm = TRUE)
      mat_list[[cell_group]] <- mat_tmp
    } else {
      compare_groups <- strsplit(colnames(mat_tmp), " ; ")
      names_keep <- names(which(
        table(sapply(compare_groups, function(x) x[[1]])) == 2
      ))
      group_keep <- which(sapply(
        compare_groups,
        function(x) x[[1]] %in% names_keep
      ))
      group_TRUE <- intersect(
        group_keep,
        which(sapply(compare_groups, function(x) x[[2]]) == "TRUE")
      )
      group_FALSE <- intersect(
        group_keep,
        which(sapply(compare_groups, function(x) x[[2]]) == "FALSE")
      )
      mat_tmp <- log2(mat_tmp[, group_TRUE] / mat_tmp[, group_FALSE])
      colnames(mat_tmp) <- gsub(" ; .*", "", colnames(mat_tmp))
      mat_tmp[is.infinite(mat_tmp)] <- max(
        abs(mat_tmp[!is.infinite(mat_tmp)]),
        na.rm = TRUE
      ) *
        ifelse(mat_tmp[is.infinite(mat_tmp)] > 0, 1, -1)
      mat_tmp[is.na(mat_tmp)] <- 0
      mat_list[[cell_group]] <- mat_tmp
      cell_groups[[cell_group]] <- factor(
        gsub(" ; .*", "", cell_groups[[cell_group]]),
        levels = unique(gsub(" ; .*", "", levels(cell_groups[[cell_group]])))
      )
    }
  }

  # data used to do clustering
  mat_split <- do.call(cbind, mat_list[feature_split_by])

  if (is.null(limits)) {
    if (!is.function(exp_method) && exp_method %in% c("zscore", "log2fc")) {
      b <- ceiling(
        min(
          abs(
            stats::quantile(do.call(cbind, mat_list), c(0.01, 0.99), na.rm = TRUE)
          ),
          na.rm = TRUE
        ) * 2
      ) / 2
      colors <- circlize::colorRamp2(
        seq(-b, b, length = 100),
        palette_colors(palette = heatmap_palette, palcolor = heatmap_palcolor)
      )
    } else {
      b <- stats::quantile(do.call(cbind, mat_list), c(0.01, 0.99), na.rm = TRUE)
      colors <- circlize::colorRamp2(
        seq(b[1], b[2], length = 100),
        palette_colors(palette = heatmap_palette, palcolor = heatmap_palcolor)
      )
    }
  } else {
    colors <- circlize::colorRamp2(
      seq(limits[1], limits[2], length = 100),
      palette_colors(palette = heatmap_palette, palcolor = heatmap_palcolor)
    )
  }

  cell_metadata <- cbind.data.frame(
    data.frame(row.names = cells, cells = cells),
    cbind.data.frame(
      srt@meta.data[
        cells,
        c(group.by, intersect(cell_annotation, colnames(srt@meta.data))),
        drop = FALSE
      ],
      Matrix::t(
        GetAssayData5(
          srt,
          assay = assay,
          layer = "data"
        )[intersect(cell_annotation, rownames(assay_use)) %||% integer(), cells, drop = FALSE]
      )
    )
  )
  feature_metadata <- cbind.data.frame(
    data.frame(
      row.names = features_unique,
      features = features,
      features_uique = features_unique
    ),
    feature_annotation_data[
      features,
      intersect(
        feature_annotation,
        colnames(feature_annotation_data)
      ),
      drop = FALSE
    ]
  )
  feature_metadata[, "duplicated"] <- feature_metadata[["features"]] %in%
    features[duplicated(features)]

  lgd <- list()
  lgd[["ht"]] <- ComplexHeatmap::Legend(
    title = exp_name,
    col_fun = colors,
    border = TRUE
  )
  if (isTRUE(add_dot)) {
    lgd[["point"]] <- ComplexHeatmap::Legend(
      labels = paste0(seq(20, 100, length.out = 5), "%"),
      title = "Percent",
      type = "points",
      pch = 21,
      size = dot_size * seq(0.2, 1, length.out = 5),
      grid_height = dot_size * seq(0.2, 1, length.out = 5) * 0.8,
      grid_width = dot_size,
      legend_gp = grid::gpar(fill = "grey30"),
      border = FALSE,
      background = "transparent",
      direction = "vertical"
    )
  }

  ha_top_list <- NULL
  cluster_columns_list <- list()
  column_split_list <- list()
  for (i in seq_along(group.by)) {
    cell_group <- group.by[i]
    cluster_columns_list[[cell_group]] <- cluster_columns
    if (is.null(split.by)) {
      column_split_list[[cell_group]] <- NULL
    } else {
      column_split_list[[cell_group]] <- factor(
        gsub(" : .*", "", levels(cell_groups[[cell_group]])),
        levels = levels(srt@meta.data[[cell_group]])
      )
    }
    if (isTRUE(cluster_column_slices) && !is.null(split.by)) {
      if (isFALSE(cluster_columns)) {
        if (nlevels(column_split_list[[cell_group]]) == 1) {
          log_message(
            "cluster_column_slices=TRUE can not be used when there is only one group.",
            message_type = "error"
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
        grid::grid.rect(gp = grid::gpar(fill = palette_colors(",
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
      anno[[cell_group]] <- ComplexHeatmap::anno_block(
        align_to = split(
          seq_along(levels(cell_groups[[cell_group]])),
          gsub(
            pattern = " : .*",
            replacement = "",
            x = levels(cell_groups[[cell_group]])
          )
        ),
        panel_fun = methods::getFunction("panel_fun", where = environment()),
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
      grid::grid.rect(gp = grid::gpar(fill = palette_colors(",
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
      anno[[split.by]] <- ComplexHeatmap::anno_block(
        align_to = split(
          seq_along(levels(cell_groups[[cell_group]])),
          gsub(
            pattern = ".* : ",
            replacement = "",
            x = levels(cell_groups[[cell_group]])
          )
        ),
        panel_fun = methods::getFunction("panel_fun", where = environment()),
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

  for (i in seq_along(raw_group_by)) {
    cell_group <- raw_group_by[i]
    if (cell_group != "All.groups") {
      lgd[[cell_group]] <- ComplexHeatmap::Legend(
        title = cell_group,
        labels = levels(srt@meta.data[[cell_group]]),
        legend_gp = grid::gpar(
          fill = palette_colors(
            levels(srt@meta.data[[cell_group]]),
            palette = raw_group_palette[i],
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
        fill = palette_colors(
          levels(srt@meta.data[[split.by]]),
          palette = cell_split_palette,
          palcolor = cell_split_palcolor
        )
      ),
      border = TRUE
    )
  }

  if (!is.null(cell_annotation)) {
    subplots_list <- list()
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
          subplots <- CellStatPlot(
            srt,
            flip = flip,
            cells = names(cell_groups[[cell_group]]),
            plot_type = "pie",
            stat.by = cellan,
            group.by = cell_group,
            split.by = split.by,
            palette = palette,
            palcolor = palcolor,
            title = NULL,
            individual = TRUE,
            combine = FALSE
          )
          subplots_list[[paste0(cellan, ":", cell_group)]] <- subplots
          graphics <- list()
          for (nm in names(subplots)) {
            funbody <- paste0(
              "
              g <- as_grob(subplots_list[['",
              cellan,
              ":",
              cell_group,
              "']]",
              "[['",
              nm,
              "']]  + facet_null() + theme_void() + theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = 'none'));
              g$name <- '",
              paste0(cellan, ":", cell_group, "-", nm),
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
            strsplit(levels(cell_groups[[cell_group]]), " : "),
            function(x) {
              if (length(x) == 2) {
                paste0(c(cell_group, x[1], x[2]), collapse = ":")
              } else {
                paste0(c(cell_group, x[1], ""), collapse = ":")
              }
            }
          )

          ha_cell <- list()
          ha_cell[[cellan]] <- ComplexHeatmap::anno_customize(
            x = x_nm,
            graphics = graphics,
            which = ifelse(flip, "row", "column"),
            border = TRUE,
            verbose = FALSE
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
            fill = palette_colors(
              cell_anno,
              palette = palette,
              palcolor = palcolor
            )
          ),
          border = TRUE
        )
      } else {
        for (cell_group in group.by) {
          subplots <- FeatureStatPlot(
            srt,
            assay = assay,
            layer = "data",
            flip = flip,
            stat.by = cellan,
            cells = names(cell_groups[[cell_group]]),
            group.by = cell_group,
            split.by = split.by,
            palette = palette,
            palcolor = palcolor,
            fill.by = "group",
            same.y.lims = TRUE,
            individual = TRUE,
            combine = FALSE
          )
          subplots_list[[paste0(cellan, ":", cell_group)]] <- subplots
          graphics <- list()
          for (nm in names(subplots)) {
            funbody <- paste0(
              "
              g <- as_grob(subplots_list[['",
              cellan,
              ":",
              cell_group,
              "']]",
              "[['",
              nm,
              "']]  + facet_null() + theme_void() + theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = 'none'));
              g$name <- '",
              paste0(cellan, ":", cell_group, "-", nm),
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
            strsplit(levels(cell_groups[[cell_group]]), " : "),
            function(x) {
              if (length(x) == 2) {
                paste0(c(cellan, cell_group, x[1], x[2]), collapse = ":")
              } else {
                paste0(c(cellan, cell_group, x[1], ""), collapse = ":")
              }
            }
          )
          ha_cell <- list()
          ha_cell[[cellan]] <- ComplexHeatmap::anno_customize(
            x = x_nm,
            graphics = graphics,
            which = ifelse(flip, "row", "column"),
            border = TRUE,
            verbose = FALSE
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
            log_message(
              "The {.pkg e1071} package was not found. Switch {.arg split_method} to {.val kmeans}",
              message_type = "warning",
              verbose = verbose
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
                  message_type = "warning",
                  verbose = verbose
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
        if (split_method == "hclust") {
          hc <- stats::hclust(
            stats::as.dist(proxyC::dist(mat_split))
          )
          row_split <- feature_split <- stats::cutree(hc, k = n_split)
        }
      }
      groupmean <- stats::aggregate(
        Matrix::t(mat_split),
        by = list(unlist(lapply(cell_groups[feature_split_by], levels))),
        mean
      )
      maxgroup <- groupmean[, 1][apply(
        groupmean[, names(row_split)],
        2,
        which.max
      )]
      maxgroup <- factor(
        maxgroup,
        levels = levels(unlist(cell_groups[feature_split_by]))
      )
      df <- data.frame(row_split = row_split, order_by = maxgroup)
      df_order <- stats::aggregate(
        df[["order_by"]],
        by = list(df[, "row_split"]),
        FUN = function(x) names(sort(table(x), decreasing = TRUE))[1]
      )
      df_order[, "row_split"] <- df_order[, "Group.1"]
      df_order[["order_by"]] <- as.numeric(
        factor(
          df_order[["x"]],
          levels = levels(maxgroup)
        )
      )
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
      if (isFALSE(cluster_rows)) {
        dend <- ComplexHeatmap::cluster_within_group(
          Matrix::t(mat_split), row_split_raw
        )
        cluster_rows <- dend
        row_split <- length(unique(row_split_raw))
      }
    }
    funbody <- paste0(
      "
      grid::grid.rect(gp = grid::gpar(fill = palette_colors(",
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
      which = ifelse(flip, "column", "row")
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
        fill = palette_colors(
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
          stats::as.dist(proxyC::dist(mat_cluster))
        )
      )
      dend_ordered <- stats::reorder(
        dend,
        wts = Matrix::colMeans(mat_cluster),
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
  ht_order <- unlist(
    suppressWarnings(
      ComplexHeatmap::row_order(
        ht_list
      )
    )
  )
  features_ordered <- rownames(mat_list[[1]])[ht_order]
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
        # ha_feature[[featan]] <- ComplexHeatmap::anno_block(
        #   align_to = split(seq_along(featan_values), featan_values),
        #   panel_fun = function(index, nm) {
        #     grid::grid.rect(gp = grid::gpar(
        #       fill = palette_colors(
        #         featan_values,
        #         palette = palette[i],
        #         palcolor = palcolor
        #       )[nm],
        #       col = NA
        #     ))
        #   },
        #   which = ifelse(flip, "column", "row")
        # )
        ha_feature[[featan]] <- ComplexHeatmap::anno_simple(
          x = as.character(featan_values),
          col = palette_colors(
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
            fill = palette_colors(
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
          colors = palette_colors(palette = palette[i], palcolor = palcolor)
        )
        ha_feature <- list()
        # ha_feature[[featan]] <- ComplexHeatmap::anno_block(
        #   align_to = split(seq_along(featan_values), featan_values),
        #   panel_fun = function(index, nm) {
        #     grid::grid.rect(gp = grid::gpar(
        #       fill = col_fun(featan_values[nm]),
        #       col = NA
        #     ))
        #   },
        #   which = ifelse(flip, "column", "row")
        # )
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
  vlnplots_list <- NULL
  x_nm_list <- NULL
  y_nm_list <- NULL
  if (fill.by == "group") {
    palette <- group_palette
    palcolor <- group_palcolor
  } else {
    palette <- feature_annotation_palette
    palcolor <- feature_annotation_palcolor
  }

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
    if (isTRUE(add_violin)) {
      vlnplots <- FeatureStatPlot(
        srt,
        assay = assay,
        layer = "data",
        flip = flip,
        stat.by = rownames(mat_list[[cell_group]]),
        cells = names(cell_groups[[cell_group]]),
        group.by = cell_group,
        split.by = split.by,
        palette = fill_palette,
        palcolor = fill_palcolor,
        fill.by = fill.by,
        same.y.lims = TRUE,
        individual = TRUE,
        combine = FALSE
      )
      lgd[["ht"]] <- NULL

      for (nm in names(vlnplots)) {
        gtable <- as_grob(
          vlnplots[[nm]] +
            ggplot2::facet_null() +
            ggplot2::theme_void() +
            ggplot2::theme(legend.position = "none")
        )
        gtable$name <- paste0(cell_group, "-", nm)
        vlnplots[[nm]] <- gtable
      }
      vlnplots_list[[paste0("heatmap_group:", cell_group)]] <- vlnplots
      x_nm <- rownames(mat_list[[cell_group]])
      x_nm_list[[paste0("heatmap_group:", cell_group)]] <- x_nm
      y_nm <- sapply(
        strsplit(levels(cell_groups[[cell_group]]), " : "),
        function(x) {
          if (length(x) == 2) {
            paste0(c(cell_group, x[1], x[2]), collapse = ":")
          } else {
            paste0(c(cell_group, x[1], ""), collapse = ":")
          }
        }
      )
      y_nm_list[[paste0("heatmap_group:", cell_group)]] <- y_nm
    }

    funbody <- paste0(
      if (isTRUE(add_dot) || isTRUE(add_violin)) {
        "grid::grid.rect(x, y,
          width = width, height = height,
          gp = grid::gpar(col = 'white', lwd = 1, fill = 'white')
        );"
      },
      if (isTRUE(add_bg)) {
        paste0(
          "
        grid::grid.rect(x, y,
          width = width, height = height,
          gp = grid::gpar(col = fill, lwd = 1, fill = adjcolors(fill, ",
          bg_alpha,
          "))
        );
        "
        )
      },
      if (isTRUE(add_reticle)) {
        paste0(
          "
        ind_mat = ComplexHeatmap::restore_matrix(j, i, x, y);
        ind_top = ind_mat[1,];
        ind_left = ind_mat[,1];
        for(col in seq_len(ncol(ind_mat))){
          grid::grid.lines(
            x = grid::unit(rep(x[ind_top[col]],each=2),'npc'),
            y = grid::unit(c(0,1),'npc'),
            gp = grid::gpar(col = '",
          reticle_color,
          "', lwd = 1.5));
        };
        for(row in seq_len(nrow(ind_mat))){
          grid::grid.lines(
            x = grid::unit(c(0,1),'npc'),
            y = grid::unit(rep(y[ind_left[row]],each=2),'npc'),
            gp = grid::gpar(col = '",
          reticle_color,
          "', lwd = 1.5));
        };
        "
        )
      },
      if (isTRUE(add_dot)) {
        paste0(
          "perc <- ComplexHeatmap::pindex(mat_perc_list[['",
          cell_group,
          "']]",
          ", i, j);
        grid::grid.points(x, y,
          pch = 21,
          size = dot_size*perc,
          gp = grid::gpar(col = 'black', lwd = 1, fill = fill)
        );
        "
        )
      },
      if (isTRUE(add_violin)) {
        if (isTRUE(flip)) {
          paste0(
            "
        groblist <- extractgrobs(vlnplots = vlnplots_list[[paste0('heatmap_group:', '",
            cell_group,
            "')]],
               x_nm =  x_nm_list[[paste0('heatmap_group:', '",
            cell_group,
            "')]],
               y_nm= y_nm_list[[paste0('heatmap_group:', '",
            cell_group,
            "')]],
               x = j,y = i);
        grid_draw(groblist, x = x, y = y, width = width, height = height);
        "
          )
        } else {
          paste0(
            "
        groblist <- extractgrobs(vlnplots = vlnplots_list[[paste0('heatmap_group:', '",
            cell_group,
            "')]],
               x_nm =  x_nm_list[[paste0('heatmap_group:', '",
            cell_group,
            "')]],
               y_nm= y_nm_list[[paste0('heatmap_group:', '",
            cell_group,
            "')]],
               x = i,y = j);
        grid_draw(groblist, x = x, y = y, width = width, height = height);
        "
          )
        }
      }
    )

    funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
    eval(
      parse(
        text = paste(
          "layer_fun <- function(j, i, x, y, width, height, fill) {",
          funbody,
          "}",
          sep = ""
        )
      ),
      envir = environment()
    )
    ht_args <- list(
      name = cell_group,
      matrix = if (flip) {
        Matrix::t(mat_list[[cell_group]])
      } else {
        mat_list[[cell_group]]
      },
      col = colors,
      layer_fun = methods::getFunction(
        "layer_fun",
        where = environment()
      ),
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
      row_split = if (flip) {
        column_split_list[[cell_group]]
      } else {
        row_split
      },
      column_split = if (flip) {
        row_split
      } else {
        column_split_list[[cell_group]]
      },
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
        "The size of the heatmap is fixed because certain elements are not scalable."
      )
      log_message(
        "The width and height of the heatmap are determined by the size of the current viewport."
      )
      log_message(
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
  } else {
    ht_width <- grid::unit(width_sum, units = units)
    ht_height <- grid::unit(height_sum, units = units)
  }
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
      g_tree = g_tree,
      matrix_list = mat_list,
      feature_split = feature_split,
      cell_metadata = cell_metadata,
      feature_metadata = feature_metadata,
      enrichment = res
    )
  )
}
