#' @title NMF similarity heatmap
#'
#' @md
#' @inheritParams GroupHeatmap
#' @param srt A Seurat object containing an NMF dimensional reduction.
#' @param plot_type Plot type. `"cells"` plots cell/spot similarity from NMF
#' embeddings. `"features"` plots feature similarity from NMF loadings.
#' @param reduction Name of the NMF reduction. Default is `"nmf"`.
#' @param dims Dimensions/components from the NMF reduction to use. If `NULL`,
#' all available dimensions are used.
#' @param cells Cells/spots to include when `plot_type = "cells"`.
#' @param features Features to include when `plot_type = "features"`. If `NULL`,
#' variable features shared with the loading matrix are used; if none are found,
#' all features in the loading matrix are used.
#' @param similarity_metric Similarity metric.
#' @param cell_annotation Metadata columns to show as column annotations in
#' cell mode.
#' @param feature_annotation Feature metadata columns to show as column
#' annotations in feature mode.
#' @param cluster_palette Palette used for NMF cluster/program annotations.
#' @param cluster_palcolor Optional custom colors for NMF cluster/program
#' annotations.
#' @param heatmap_limits Numeric breaks for the heatmap color scale. If `NULL`,
#' defaults to `c(0, 0.35, 0.75, 1)`.
#'
#' @return
#' A list with the following elements:
#' \itemize{
#'   \item `plot`: The heatmap plot as a patchwork/ggplot object.
#'   \item `similarity_matrix`: The ordered similarity matrix used for plotting.
#'   \item `nmf_cluster`: The ordered NMF cluster/program assignment.
#'   \item `order`: The ordered row/column names.
#'   \item `metadata`: Ordered cell or feature metadata used for annotations.
#'   \item `enrichment`: Enrichment results for feature mode when requested,
#'   otherwise `NULL`.
#' }
#'
#' @seealso [RunNMF]
#'
#' @export
#'
#' @examples
#' library(Matrix)
#' data(pancreas_sub)
#' pancreas_sub <- Seurat::NormalizeData(pancreas_sub, verbose = FALSE)
#' pancreas_sub <- Seurat::FindVariableFeatures(
#'   pancreas_sub,
#'   nfeatures = 1000,
#'   verbose = FALSE
#' )
#' pancreas_sub <- RunNMF(
#'   pancreas_sub,
#'   features = SeuratObject::VariableFeatures(pancreas_sub),
#'   nbes = 5,
#'   maxit = 50,
#'   verbose = FALSE
#' )
#' ht_cells <- NMFHeatmap(
#'   pancreas_sub,
#'   plot_type = "cells",
#'   cell_annotation = "CellType",
#'   width = 3,
#'   height = 0.5
#' )
#' ht_cells$plot
#'
#' ht_features <- NMFHeatmap(
#'   pancreas_sub,
#'   plot_type = "features",
#'   width = 3,
#'   height = 0.5
#' )
#' ht_features$plot
NMFHeatmap <- function(
  srt,
  plot_type = c("cells", "features"),
  reduction = "nmf",
  dims = NULL,
  cells = NULL,
  features = NULL,
  similarity_metric = "cosine",
  cell_annotation = NULL,
  feature_annotation = NULL,
  assay = NULL,
  border = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_names_side = "left",
  column_names_side = "top",
  row_names_rot = 0,
  column_names_rot = 90,
  row_title = NULL,
  column_title = NULL,
  anno_terms = FALSE,
  anno_keys = FALSE,
  anno_features = FALSE,
  terms_width = grid::unit(4, "in"),
  terms_fontsize = 8,
  terms_stat = "none",
  terms_stat_digits = 2,
  terms_stat_test = TRUE,
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
  Ensembl_version = NULL,
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
  heatmap_palette = "simspec",
  heatmap_palcolor = c("#ffffe5", "#d9f0d3", "#74add1", "#2166ac"),
  heatmap_limits = NULL,
  cluster_palette = "simspec",
  cluster_palcolor = NULL,
  cell_annotation_palette = "Chinese",
  cell_annotation_palcolor = NULL,
  feature_annotation_palette = "Dark2",
  feature_annotation_palcolor = NULL,
  use_raster = NULL,
  raster_device = "png",
  raster_by_magick = FALSE,
  height = NULL,
  width = NULL,
  units = "inch",
  cores = 1,
  seed = 11,
  legend.position = "right",
  ht_params = list(),
  verbose = TRUE
) {
  set.seed(seed)
  plot_type <- match.arg(plot_type)
  similarity_metric <- match.arg(similarity_metric, choices = "cosine")
  if (isTRUE(raster_by_magick)) {
    check_r("magick", verbose = FALSE)
  }
  if (!reduction %in% names(srt@reductions)) {
    log_message(
      "{.arg reduction} {.val {reduction}} is not present in {.arg srt}.",
      message_type = "error"
    )
  }

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  nmf_data <- nmf_heatmap_input(
    srt = srt,
    plot_type = plot_type,
    reduction = reduction,
    dims = dims,
    cells = cells,
    features = features
  )
  mat <- nmf_data[["matrix"]]
  dims_use <- nmf_data[["dims"]]

  n_items <- nrow(mat)
  similarity_size_gib <- sprintf("%.2f", n_items^2 * 8 / 1024^3)
  log_message(
    "{.fn NMFHeatmap} input: {.val {n_items}} {plot_type} x {.val {length(dims_use)}} NMF dimensions. Computing a {.val {n_items}} x {.val {n_items}} similarity matrix (~{similarity_size_gib} GiB dense numeric matrix).",
    verbose = verbose
  )
  if (n_items > 5000) {
    log_message(
      "{.fn NMFHeatmap} similarity calculation scales quadratically with the number of rows; use {.arg cells}, {.arg features}, or fewer {.arg dims} if this step is slow on the current machine.",
      message_type = "warning",
      verbose = verbose
    )
  }
  similarity <- nmf_heatmap_cosine(mat)
  log_message("Ordering {.fn NMFHeatmap} rows and columns ...", verbose = verbose)
  nmf_strength <- apply(mat, 1, max, na.rm = TRUE)
  nmf_cluster <- max.col(mat, ties.method = "first")
  nmf_cluster <- factor(nmf_cluster, levels = sort(unique(nmf_cluster)))
  names(nmf_cluster) <- rownames(mat)
  names(nmf_strength) <- rownames(mat)

  item_order <- order(nmf_cluster, -nmf_strength, rownames(mat))
  similarity <- similarity[item_order, item_order, drop = FALSE]
  nmf_cluster <- nmf_cluster[rownames(similarity)]
  nmf_strength <- nmf_strength[rownames(similarity)]

  cluster_title <- if (identical(plot_type, "cells")) {
    "NMF cluster"
  } else {
    "NMF program"
  }
  cluster_cols <- palette_colors(
    levels(nmf_cluster),
    palette = cluster_palette,
    palcolor = cluster_palcolor
  )
  cluster_cols <- cluster_cols[levels(nmf_cluster)]
  metadata <- nmf_heatmap_metadata(
    srt = srt,
    assay = assay,
    plot_type = plot_type,
    ids = rownames(similarity),
    annotations = if (identical(plot_type, "cells")) {
      cell_annotation
    } else {
      feature_annotation
    },
    nmf_cluster = nmf_cluster,
    nmf_strength = nmf_strength,
    cluster_title = cluster_title
  )

  top_annotation <- nmf_heatmap_top_annotation(
    metadata = metadata,
    annotations = if (identical(plot_type, "cells")) {
      cell_annotation
    } else {
      feature_annotation
    },
    cluster_title = cluster_title,
    cluster_cols = cluster_cols,
    annotation_palette = if (identical(plot_type, "cells")) {
      cell_annotation_palette
    } else {
      feature_annotation_palette
    },
    annotation_palcolor = if (identical(plot_type, "cells")) {
      cell_annotation_palcolor
    } else {
      feature_annotation_palcolor
    },
    border = border
  )
  left_annotation <- do.call(
    ComplexHeatmap::rowAnnotation,
    args = c(
      stats::setNames(list(nmf_cluster), cluster_title),
      list(
        col = stats::setNames(list(cluster_cols), cluster_title),
        annotation_name_side = "bottom",
        annotation_width = grid::unit(5, "mm"),
        border = border
      )
    )
  )

  enrichment <- list(res = NULL, ha_right = NULL)
  if (
    identical(plot_type, "features") &&
      any(c(anno_terms, anno_keys, anno_features))
  ) {
    log_message(
      "Running enrichment annotation for {.fn NMFHeatmap} feature clusters ...",
      verbose = verbose
    )
    enrichment <- heatmap_enrichment(
      geneID = rownames(similarity),
      geneID_groups = nmf_cluster,
      feature_split_palette = cluster_palette,
      feature_split_palcolor = cluster_palcolor,
      ha_right = NULL,
      flip = FALSE,
      anno_terms = anno_terms,
      anno_keys = anno_keys,
      anno_features = anno_features,
      terms_width = terms_width,
      terms_fontsize = terms_fontsize,
      terms_stat = terms_stat,
      terms_stat_digits = terms_stat_digits,
      terms_stat_test = terms_stat_test,
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
      words_excluded = words_excluded,
      cores = cores
    )
  } else if (
    !identical(plot_type, "features") &&
      any(c(anno_terms, anno_keys, anno_features))
  ) {
    log_message(
      "Enrichment annotations are only available when {.arg plot_type = 'features'}.",
      message_type = "warning",
      verbose = verbose
    )
  }
  lgd <- enrichment[["lgd"]]

  heatmap_col <- nmf_heatmap_col_fun(
    values = similarity,
    limits = heatmap_limits,
    palette = heatmap_palette,
    palcolor = heatmap_palcolor
  )
  ht_name <- if (identical(plot_type, "cells")) {
    "Similarity\nindex"
  } else {
    "Loading\nsimilarity"
  }
  if (is.null(use_raster)) {
    use_raster <- nrow(similarity) * ncol(similarity) > 1e6
  }

  log_message(
    "Building {.pkg ComplexHeatmap} object for {.fn NMFHeatmap} ...",
    verbose = verbose
  )
  ht_args <- c(
    list(
      matrix = similarity,
      name = ht_name,
      col = heatmap_col,
      top_annotation = top_annotation,
      left_annotation = left_annotation,
      right_annotation = enrichment[["ha_right"]],
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      row_names_side = row_names_side,
      column_names_side = column_names_side,
      row_names_rot = row_names_rot,
      column_names_rot = column_names_rot,
      row_title = row_title,
      column_title = column_title,
      border = border,
      use_raster = use_raster,
      raster_device = raster_device,
      raster_by_magick = raster_by_magick,
      heatmap_legend_param = list(
        at = c(0, 0.5, 1),
        labels = c("0", "0.5", "1.0")
      )
    ),
    ht_params
  )
  ht_list <- do.call(ComplexHeatmap::Heatmap, args = ht_args)

  plot <- nmf_heatmap_render_plot(
    ht_list = ht_list,
    ht_name = ht_name,
    nmf_cluster = nmf_cluster,
    top_annotation = top_annotation,
    left_annotation = left_annotation,
    right_annotation = enrichment[["ha_right"]],
    db = db,
    anno_terms = anno_terms,
    anno_keys = anno_keys,
    anno_features = anno_features,
    width = width,
    height = height,
    units = units,
    legend_list = lgd,
    legend.position = legend.position,
    verbose = verbose
  )

  list(
    plot = plot,
    similarity_matrix = similarity,
    nmf_cluster = nmf_cluster,
    order = rownames(similarity),
    metadata = metadata,
    enrichment = enrichment[["res"]],
    dims = dims_use,
    similarity_metric = similarity_metric
  )
}

nmf_heatmap_input <- function(
  srt,
  plot_type,
  reduction,
  dims,
  cells,
  features
) {
  if (identical(plot_type, "cells")) {
    mat <- Seurat::Embeddings(srt, reduction = reduction)
    if (!is.null(cells)) {
      cells_use <- intersect(cells, rownames(mat))
      if (length(cells_use) == 0) {
        log_message(
          "No requested {.arg cells} were found in the NMF embedding.",
          message_type = "error"
        )
      }
      mat <- mat[cells_use, , drop = FALSE]
    }
  } else {
    mat <- SeuratObject::Loadings(srt[[reduction]])
    if (is.null(features)) {
      features <- intersect(SeuratObject::VariableFeatures(srt), rownames(mat))
      if (length(features) == 0) {
        features <- rownames(mat)
      }
    } else {
      features <- intersect(features, rownames(mat))
      if (length(features) == 0) {
        log_message(
          "No requested {.arg features} were found in the NMF loading matrix.",
          message_type = "error"
        )
      }
    }
    mat <- mat[features, , drop = FALSE]
  }
  dims <- nmf_heatmap_dims(dims = dims, n_dims = ncol(mat))
  mat <- as_matrix(mat[, dims, drop = FALSE])
  keep <- Matrix::rowSums(abs(mat)) > 0
  if (!all(keep)) {
    mat <- mat[keep, , drop = FALSE]
  }
  if (nrow(mat) < 2) {
    log_message(
      "At least two rows with non-zero NMF values are required.",
      message_type = "error"
    )
  }
  list(matrix = mat, dims = dims)
}

nmf_heatmap_dims <- function(dims, n_dims) {
  if (is.null(dims)) {
    return(seq_len(n_dims))
  }
  dims <- unique(as.integer(dims))
  dims <- dims[!is.na(dims)]
  dims <- dims[dims >= 1 & dims <= n_dims]
  if (length(dims) == 0) {
    log_message(
      "{.arg dims} does not select any available NMF dimensions.",
      message_type = "error"
    )
  }
  dims
}

nmf_heatmap_cosine <- function(mat) {
  norm <- sqrt(Matrix::rowSums(mat^2))
  norm[norm == 0 | !is.finite(norm)] <- 1
  mat_norm <- mat / norm
  mat_norm[!is.finite(mat_norm)] <- 0
  sim <- tcrossprod(mat_norm)
  sim[!is.finite(sim)] <- 0
  sim
}

nmf_heatmap_metadata <- function(
  srt,
  assay,
  plot_type,
  ids,
  annotations,
  nmf_cluster,
  nmf_strength,
  cluster_title
) {
  if (identical(plot_type, "cells")) {
    metadata <- srt@meta.data[
      ids,
      intersect(annotations, colnames(srt@meta.data)),
      drop = FALSE
    ]
  } else {
    feature_data <- GetFeaturesData(srt, assay = assay)
    metadata <- feature_data[
      ids,
      intersect(annotations, colnames(feature_data)),
      drop = FALSE
    ]
  }
  metadata <- as.data.frame(metadata)
  metadata[[cluster_title]] <- nmf_cluster
  metadata[["Max loading"]] <- nmf_strength
  metadata <- metadata[,
    c(
      cluster_title,
      "Max loading",
      setdiff(colnames(metadata), c(cluster_title, "Max loading"))
    ),
    drop = FALSE
  ]
  metadata
}

nmf_heatmap_top_annotation <- function(
  metadata,
  annotations,
  cluster_title,
  cluster_cols,
  annotation_palette,
  annotation_palcolor,
  border
) {
  annotations <- intersect(annotations, colnames(metadata))
  show_columns <- c(cluster_title, annotations)
  if (length(show_columns) == 0) {
    return(NULL)
  }
  ann_list <- list()
  col_list <- list()
  for (i in seq_along(show_columns)) {
    ann <- show_columns[i]
    values <- metadata[[ann]]
    names(values) <- rownames(metadata)
    ann_list[[ann]] <- values
    if (identical(ann, cluster_title)) {
      col_list[[ann]] <- cluster_cols
    } else {
      palette <- annotation_palette[min(i - 1, length(annotation_palette))]
      palcolor <- NULL
      if (!is.null(annotation_palcolor)) {
        if (is.list(annotation_palcolor)) {
          palcolor <- annotation_palcolor[[min(
            i - 1,
            length(annotation_palcolor)
          )]]
        } else {
          palcolor <- annotation_palcolor
        }
      }
      col_list[[ann]] <- nmf_heatmap_annotation_col(values, palette, palcolor)
    }
  }
  do.call(
    ComplexHeatmap::HeatmapAnnotation,
    args = c(
      ann_list,
      list(
        col = col_list,
        which = "column",
        annotation_name_side = "left",
        border = border
      )
    )
  )
}

nmf_heatmap_annotation_col <- function(values, palette, palcolor) {
  if (is.numeric(values)) {
    vals <- values[is.finite(values)]
    if (length(vals) == 0) {
      vals <- c(0, 1)
    }
    breaks <- unique(as.numeric(stats::quantile(
      vals,
      probs = c(0.02, 0.5, 0.98),
      na.rm = TRUE
    )))
    if (length(breaks) == 1) {
      breaks <- breaks + c(-1, 0, 1)
    }
    colors <- palette_colors(palette = palette, palcolor = palcolor)
    return(circlize::colorRamp2(
      breaks = seq(min(breaks), max(breaks), length.out = length(colors)),
      colors = colors
    ))
  }
  if (is.logical(values)) {
    values <- factor(values, levels = c(TRUE, FALSE))
  } else if (!is.factor(values)) {
    values <- factor(values, levels = unique(values))
  }
  palette_colors(levels(values), palette = palette, palcolor = palcolor)
}

nmf_heatmap_col_fun <- function(values, limits, palette, palcolor) {
  colors <- if (!is.null(palcolor)) {
    palcolor
  } else {
    palette_colors(palette = palette)
  }
  if (is.null(limits)) {
    limits <- if (length(colors) == 4) {
      c(0, 0.35, 0.75, 1)
    } else {
      seq(0, 1, length.out = length(colors))
    }
  }
  if (length(limits) == 2) {
    limits <- seq(limits[1], limits[2], length.out = length(colors))
  }
  if (length(limits) != length(colors)) {
    log_message(
      "{.arg heatmap_limits} must have length 2 or the same length as heatmap colors.",
      message_type = "error"
    )
  }
  circlize::colorRamp2(limits, colors)
}

nmf_heatmap_render_plot <- function(
  ht_list,
  ht_name,
  nmf_cluster,
  top_annotation,
  left_annotation,
  right_annotation,
  db,
  anno_terms,
  anno_keys,
  anno_features,
  width,
  height,
  units,
  legend_list = NULL,
  legend.position = "right",
  verbose = TRUE
) {
  fix <- !is.null(width) || !is.null(height) || !is.null(right_annotation)
  log_message(
    "Calculating {.fn NMFHeatmap} render size ...",
    verbose = verbose
  )
  rendersize <- heatmap_rendersize(
    width = width,
    height = height,
    units = units,
    ha_top_list = list(top_annotation),
    ha_left = left_annotation,
    ha_right = right_annotation,
    ht_list = ht_list,
    legend_list = legend_list,
    flip = FALSE
  )
  width_sum <- rendersize[["width_sum"]]
  height_sum <- rendersize[["height_sum"]]
  if (isTRUE(fix)) {
    log_message(
      "Fixing {.fn NMFHeatmap} panel size ...",
      verbose = verbose
    )
    fixsize_env <- new.env(parent = emptyenv())
    invisible(grid::grid.grabExpr(
      {
        fixsize_env[["value"]] <- heatmap_fixsize(
          width = width,
          width_sum = width_sum,
          height = height,
          height_sum = height_sum,
          units = units,
          ht_list = ht_list,
          legend_list = legend_list
        )
      }
    ))
    fixsize <- fixsize_env[["value"]]
    rm(fixsize_env)
    ht_width <- fixsize[["ht_width"]]
    ht_height <- fixsize[["ht_height"]]
  } else {
    ht_width <- grid::unit(width_sum, units = units)
    ht_height <- grid::unit(height_sum, units = units)
  }
  log_message(
    "Drawing {.fn NMFHeatmap}; this can take time for large similarity matrices ...",
    verbose = verbose
  )
  g_tree <- grid::grid.grabExpr(
    {
      ComplexHeatmap::draw(
        ht_list,
        heatmap_legend_side = legend.position,
        annotation_legend_list = legend_list,
        annotation_legend_side = legend.position
      )
      nmf_heatmap_decorate_blocks(ht_name, nmf_cluster)
      if (
        !is.null(right_annotation) &&
          any(c(anno_terms, anno_keys, anno_features))
      ) {
        nmf_heatmap_decorate_enrichment(right_annotation, db)
      }
    },
    width = ht_width,
    height = ht_height,
    wrap = TRUE,
    wrap.grobs = TRUE
  )
  log_message(
    "Assembling {.fn NMFHeatmap} plot object ...",
    verbose = verbose
  )
  if (isTRUE(fix)) {
    panel_fix_overall(
      g_tree,
      width = as.numeric(ht_width),
      height = as.numeric(ht_height),
      units = units
    )
  } else {
    patchwork::wrap_plots(g_tree)
  }
}

nmf_heatmap_decorate_blocks <- function(ht_name, nmf_cluster) {
  cluster_runs <- split(seq_along(nmf_cluster), nmf_cluster)
  ComplexHeatmap::decorate_heatmap_body(ht_name, {
    n <- length(nmf_cluster)
    for (idx in cluster_runs) {
      start <- min(idx)
      end <- max(idx)
      center <- ((start + end) / 2 - 0.5) / n
      size <- length(idx) / n
      grid::grid.rect(
        x = grid::unit(center, "npc"),
        y = grid::unit(1 - center, "npc"),
        width = grid::unit(size, "npc"),
        height = grid::unit(size, "npc"),
        gp = grid::gpar(fill = NA, col = "black", lwd = 1.4)
      )
    }
  })
}

nmf_heatmap_decorate_enrichment <- function(right_annotation, db) {
  for (enrich in db) {
    enrich_anno <- names(right_annotation)[grep(
      paste0("_split_", enrich),
      names(right_annotation)
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
}
