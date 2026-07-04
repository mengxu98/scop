#' @title Run scFEA flux estimation for a Seurat object
#'
#' @param srt A Seurat object.
#' @param assay Assay to use as expression matrix. Default is
#' `DefaultAssay(srt)`.
#' @param layer Assay layer to use. Default is `"data"`.
#' @param species One of `"human"` or `"mouse"`, selecting the M168 scFEA files.
#' @param n_epoch Number of scFEA training epochs.
#' @param sc_imputation Whether to run MAGIC imputation inside the scFEA
#' backend.
#' @param assay_flux Name of the assay storing module flux scores.
#' @param assay_balance Name of the assay storing metabolite balance scores.
#' @param store_metadata Whether to also append flux and balance values to
#' `srt@meta.data`.
#' @param data_dir Optional directory containing scFEA M168 CSV resources. If
#' `NULL`, files are downloaded from `mengxu98/datasets` and cached with
#' `tools::R_user_dir("scop", "data")`.
#' @param seed Random seed passed to R and the Python scFEA backend.
#' @param max_cells Maximum number of cells used for GNN training. When the
#' input has more cells, a random subset is sampled for training and the
#' trained model predicts fluxes for all cells in batches. This drastically
#' reduces peak memory for large datasets. Set `NULL` (default) to train on
#' all cells, matching the original upstream behaviour. A value such as
#' `20000` can be useful on machines with 16-32 GiB RAM.
#' @param verbose Whether to print progress messages.
#'
#' @return A Seurat object with `assay_flux`, `assay_balance`, and
#' `srt@tools[["scFEA"]]`.
#'
#' @details
#' The scFEA GNN architecture, loss function, training loop, and
#' M168 human/mouse model files are downloaded from the `mengxu98/datasets`
#' repository and cached outside the package; the R side prepares a Seurat
#' expression matrix and stores flux / balance outputs back into Seurat.
#'
#' scFEA is licensed for academic, non-commercial use. See the bundled
#' `inst/python/scfea/LICENSE` file for details. Data resources are downloaded
#' from \url{https://github.com/mengxu98/datasets/tree/main/scFEA}.
#'
#' @export
RunscFEA <- function(
  srt,
  assay = NULL,
  layer = "data",
  species = c("human", "mouse"),
  n_epoch = 100,
  sc_imputation = FALSE,
  assay_flux = "scFEAflux",
  assay_balance = "scFEAbalance",
  store_metadata = FALSE,
  data_dir = NULL,
  seed = 16,
  max_cells = NULL,
  verbose = TRUE
) {
  species <- match.arg(species)
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (!is.numeric(n_epoch) || length(n_epoch) != 1 || n_epoch < 1) {
    log_message(
      "{.arg n_epoch} must be a positive number",
      message_type = "error"
    )
  }
  if (!is.null(max_cells) && (
    !is.numeric(max_cells) || length(max_cells) != 1 || is.na(max_cells) ||
      max_cells < 1
  )) {
    log_message(
      "{.arg max_cells} must be {.code NULL} or a positive number",
      message_type = "error"
    )
  }

  set.seed(seed)
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  log_message(
    "Running {.pkg scFEA} on assay {.val {assay}} layer {.val {layer}}",
    verbose = verbose
  )

  configure_python_thread_env()
  if (isTRUE(reticulate::py_available(initialize = FALSE))) {
    tryCatch(
      reticulate::py_run_string(
        "
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['KMP_WARNINGS'] = '0'
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
os.environ['NUMBA_NUM_THREADS'] = '1'
"
      ),
      error = function(...) NULL
    )
  }

  ok <- check_python(
    c("torch", "numpy", "pandas", "tqdm"),
    verbose = verbose
  )
  if (!isTRUE(ok)) {
    log_message(
      paste0(
        "Missing Python dependencies for {.pkg scFEA}. ",
        "Run {.fn check_python} with {.val torch}, {.val numpy}, ",
        "{.val pandas}, and {.val tqdm}, or run {.fn PrepareEnv}."
      ),
      message_type = "error"
    )
  }
  if (isTRUE(sc_imputation)) {
    ok_magic <- check_python("magic-impute", verbose = verbose)
    if (!isTRUE(ok_magic)) {
      log_message(
        "{.arg sc_imputation = TRUE} requires Python package {.pkg magic-impute}",
        message_type = "error"
      )
    }
  }

  expr_mat <- GetAssayData5(srt, assay = assay, layer = layer)
  if (is.null(expr_mat)) {
    log_message(
      "Layer {.val {layer}} was not found in assay {.val {assay}}",
      message_type = "error"
    )
  }
  if (!inherits(expr_mat, c("dgCMatrix", "matrix", "Matrix"))) {
    expr_mat <- as_matrix(expr_mat)
  }
  if (nrow(expr_mat) == 0 || ncol(expr_mat) == 0) {
    log_message(
      "The selected expression matrix is empty",
      message_type = "error"
    )
  }
  n_input_genes <- nrow(expr_mat)
  n_input_cells <- ncol(expr_mat)

  data_dir <- resolve_scfea_dir(data_dir = data_dir, verbose = verbose)
  module_genes <- scfea_module_genes(species, data_dir = data_dir)
  overlap_genes <- intersect(rownames(expr_mat), module_genes)
  if (length(overlap_genes) == 0) {
    log_message(
      "No genes overlap the {.pkg scFEA} M168 module gene file",
      message_type = "error"
    )
  }
  log_message(
    "Input matrix: {.val {nrow(expr_mat)}} genes x {.val {ncol(expr_mat)}} cells; ",
    "{.val {length(overlap_genes)}} genes overlap scFEA M168 modules",
    verbose = verbose
  )

  keep_genes <- Matrix::rowSums(expr_mat) > 0
  n_nonzero_genes <- sum(keep_genes)
  scfea_genes <- intersect(rownames(expr_mat)[keep_genes], module_genes)
  if (length(scfea_genes) == 0) {
    log_message(
      "No non-zero genes overlap the {.pkg scFEA} M168 module gene file",
      message_type = "error"
    )
  }
  expr_mat <- expr_mat[scfea_genes, , drop = FALSE]

  dense_bytes <- prod(dim(expr_mat)) * 8
  dense_gib <- dense_bytes / (1024^3)
  if (dense_gib > 2) {
    log_message(
      "Dense coercion of {.val {nrow(expr_mat)}} genes x {.val {ncol(expr_mat)}} cells ",
      "would allocate approximately {.val {round(dense_gib, 1)}} GiB. ",
      "Consider subsetting cells or using a machine with more memory.",
      message_type = "warning"
    )
  }
  log_message(
    "Using {.val {nrow(expr_mat)}} non-zero overlapping genes for scFEA input",
    verbose = verbose
  )

  expr_cells_genes <- as_matrix(Matrix::t(expr_mat))
  storage.mode(expr_cells_genes) <- "double"

  conda <- resolve_conda("auto")
  envname <- get_envname(NULL)
  python_path <- conda_python(conda = conda, envname = envname)
  configure_python_runtime(python_path)

  scfea <- reticulate::import_from_path(
    "scfea",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = FALSE
  )
  pd <- reticulate::import("pandas", convert = FALSE)

  py_expr <- pd$DataFrame(
    data = reticulate::r_to_py(expr_cells_genes, convert = FALSE),
    index = reticulate::r_to_py(rownames(expr_cells_genes), convert = FALSE),
    columns = reticulate::r_to_py(colnames(expr_cells_genes), convert = FALSE)
  )

  result_py <- scfea$run_scfea(
    py_expr,
    species = species,
    sc_imputation = isTRUE(sc_imputation),
    n_epoch = as.integer(n_epoch),
    seed = as.integer(seed),
    verbose = isTRUE(verbose),
    data_dir = data_dir,
    max_cells = if (is.null(max_cells)) NULL else as.integer(max_cells)
  )

  flux_df <- scfea_py_dataframe_to_r(
    result_py[["flux"]],
    expected_rows = colnames(expr_mat)
  )
  balance_df <- scfea_py_dataframe_to_r(
    result_py[["balance"]],
    expected_rows = colnames(expr_mat)
  )

  missing_flux_cells <- setdiff(colnames(expr_mat), rownames(flux_df))
  missing_balance_cells <- setdiff(colnames(expr_mat), rownames(balance_df))
  if (length(missing_flux_cells) > 0 || length(missing_balance_cells) > 0) {
    log_message(
      "scFEA output cell names do not match the input Seurat cells",
      message_type = "error"
    )
  }
  flux_df <- flux_df[colnames(expr_mat), , drop = FALSE]
  balance_df <- balance_df[colnames(expr_mat), , drop = FALSE]

  flux_mat <- t(as_matrix(flux_df))
  storage.mode(flux_mat) <- "double"
  rownames(flux_mat) <- scfea_module_feature_id(rownames(flux_mat))
  colnames(flux_mat) <- colnames(expr_mat)

  balance_mat <- t(as_matrix(balance_df))
  storage.mode(balance_mat) <- "double"
  colnames(balance_mat) <- colnames(expr_mat)

  module_info <- scfea_module_info(data_dir = data_dir)
  module_info <- module_info[match(rownames(flux_mat), module_info$feature_id), , drop = FALSE]
  rownames(module_info) <- module_info$feature_id

  compound_info <- scfea_compound_info(species, data_dir = data_dir)
  compound_info <- compound_info[match(rownames(balance_mat), compound_info$compound_name), , drop = FALSE]
  rownames(compound_info) <- compound_info$compound_name

  srt[[assay_flux]] <- Seurat::CreateAssayObject(data = flux_mat)
  srt[[assay_flux]] <- Seurat::AddMetaData(
    srt[[assay_flux]],
    metadata = module_info
  )
  srt[[assay_balance]] <- Seurat::CreateAssayObject(data = balance_mat)
  srt[[assay_balance]] <- Seurat::AddMetaData(
    srt[[assay_balance]],
    metadata = compound_info
  )

  if (isTRUE(store_metadata)) {
    flux_meta <- as.data.frame(t(flux_mat), check.names = FALSE)
    colnames(flux_meta) <- paste0(assay_flux, "_", make.names(rownames(flux_mat)))
    balance_meta <- as.data.frame(t(balance_mat), check.names = FALSE)
    colnames(balance_meta) <- paste0(assay_balance, "_", make.names(rownames(balance_mat)))
    srt <- Seurat::AddMetaData(srt, metadata = flux_meta)
    srt <- Seurat::AddMetaData(srt, metadata = balance_meta)
  }

  loss <- reticulate::py_to_r(result_py[["predictions"]][["loss"]])
  srt@tools[["scFEA"]] <- list(
    flux = flux_mat,
    balance = balance_mat,
    module_info = module_info,
    compound_info = compound_info,
    loss = loss,
    species = species,
    n_epoch = as.integer(n_epoch),
    seed = as.integer(seed),
    assay = assay,
    layer = layer,
    data_dir = data_dir,
    assay_flux = assay_flux,
    assay_balance = assay_balance,
    sc_imputation = isTRUE(sc_imputation),
    max_cells = max_cells,
    input_summary = list(
      n_input_genes = n_input_genes,
      n_nonzero_genes = n_nonzero_genes,
      n_input_cells = n_input_cells,
      n_module_genes = length(module_genes),
      n_overlap_genes = length(overlap_genes),
      n_scfea_genes = nrow(expr_mat),
      overlap_genes = overlap_genes
    )
  )

  log_message(
    "Stored scFEA flux assay {.val {assay_flux}} ({.val {nrow(flux_mat)}} modules) ",
    "and balance assay {.val {assay_balance}} ({.val {nrow(balance_mat)}} metabolites)",
    verbose = verbose
  )
  srt
}

#' Plot scFEA module flux heatmap
#'
#' @param srt A Seurat object returned by [RunscFEA()].
#' @param assay Flux assay name.
#' @param layer Flux assay layer.
#' @param group.by Metadata column used to aggregate cells.
#' @param features Optional scFEA modules to plot. Supports Seurat-safe feature
#' ids (`"M-1"`), scFEA ids (`"M_1"`), labels (`"M1"`), or reaction labels.
#' @param modules Optional scFEA modules to plot. This is a user-facing alias
#' for `features` that also accepts combined labels such as
#' `"M150: PRPP -> UMP"`.
#' @param label_by Row label mode: module id, reaction, or both.
#' @param add_sm_anno Whether to add `SM_anno` row annotation and split rows by
#' pathway class.
#' @param scale_rows Whether to z-score each row after group aggregation.
#' @param cluster_rows,cluster_columns Passed to [ComplexHeatmap::Heatmap()].
#' @param sm_anno_label_rot Rotation angle for `SM_anno` row-split labels.
#' @param show_row_names Whether to show row labels.
#' @param show_column_names Whether to show aggregated group labels.
#' @param heatmap_limit Numeric clipping limit for row z-scores.
#' @param heatmap_column_width Width for each aggregated group column. Numeric
#' values are interpreted as millimeters.
#' @param column_names_rot Rotation angle for heatmap column labels.
#' @param heatmap_palette,heatmap_palcolor Continuous heatmap palette passed to
#' `palette_colors()`.
#' @param group_palette,group_palcolor Palette for the top group annotation.
#' @param feature_split_palette,feature_split_palcolor Palette for `SM_anno`
#' row annotation.
#' @param border Whether to draw borders around the heatmap body and
#' annotations. Kept for backward compatibility. The more specific
#' `heatmap_border`, `cell_annotation_border`, and `feature_annotation_border`
#' arguments inherit from this value when left as `NULL`.
#' @param heatmap_border,cell_annotation_border,feature_annotation_border
#' Whether to draw borders for the heatmap body, cell annotations, and feature
#' annotations, respectively. Defaults inherit from `border`.
#' @param heatmap_border_palcolor,cell_annotation_border_palcolor,feature_annotation_border_palcolor
#' Border colors for the heatmap body, cell annotations, and feature
#' annotations when their matching border argument is `TRUE`.
#' @param heatmap_border_size,cell_annotation_border_size,feature_annotation_border_size
#' Border line widths for the heatmap body, cell annotations, and feature
#' annotations when their matching border argument is `TRUE`.
#' @param use_raster,raster_by_magick Raster settings passed to
#' [ComplexHeatmap::Heatmap()].
#' @param column_names_gp,column_title_gp,row_title_gp Font settings passed to
#' [grid::gpar()]-aware ComplexHeatmap arguments.
#' @param legend_title_gp,legend_labels_gp Font settings for heatmap and
#' annotation legends.
#' @param row_gap Gap between `SM_anno` row splits.
#' @param sm_anno_size Size of the left `SM_anno` annotation strip.
#' @param group_anno_size Size of the top group annotation strip.
#' @param mark_features Optional scFEA modules to mark with
#' [ComplexHeatmap::anno_mark()]. Accepts the same module ids and labels as
#' `features`.
#' @param mark_label_by Label mode for marked rows. Defaults to `label_by`.
#' @param mark_labels_gp Font settings for marked row labels. If `NULL`, a
#' compact default is chosen from `mark_label_by`.
#' @param mark_annotation_width Width of the right mark annotation. Numeric
#' values are interpreted as millimeters.
#' @param mark_link_width Link width used by [ComplexHeatmap::anno_mark()].
#' @param name Heatmap legend title.
#' @param ht_params Additional parameters passed to
#' [ComplexHeatmap::Heatmap()], overriding defaults when names overlap.
#' @param width,height Suggested export size. Stored as attributes on the
#' returned heatmap object.
#' @param verbose Whether to print messages.
#'
#' @return A `ComplexHeatmap` heatmap object.
#' @export
scFEAHeatmap <- function(
  srt,
  assay = "scFEAflux",
  layer = "data",
  group.by = NULL,
  features = NULL,
  modules = NULL,
  label_by = c("module", "reaction", "module_reaction"),
  add_sm_anno = TRUE,
  scale_rows = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  sm_anno_label_rot = 0,
  show_row_names = FALSE,
  show_column_names = FALSE,
  heatmap_limit = 2,
  heatmap_column_width = NULL,
  column_names_rot = 90,
  heatmap_palette = "RdBu",
  heatmap_palcolor = NULL,
  group_palette = "Chinese",
  group_palcolor = NULL,
  feature_split_palette = "simspec",
  feature_split_palcolor = NULL,
  border = TRUE,
  heatmap_border = NULL,
  cell_annotation_border = NULL,
  feature_annotation_border = NULL,
  heatmap_border_palcolor = "black",
  cell_annotation_border_palcolor = "black",
  feature_annotation_border_palcolor = "black",
  heatmap_border_size = 1,
  cell_annotation_border_size = 1,
  feature_annotation_border_size = 1,
  use_raster = TRUE,
  raster_by_magick = FALSE,
  column_names_gp = grid::gpar(fontsize = 8),
  column_title_gp = grid::gpar(fontsize = 11),
  row_title_gp = grid::gpar(fontsize = 7.5),
  legend_title_gp = grid::gpar(fontsize = 8),
  legend_labels_gp = grid::gpar(fontsize = 7),
  row_gap = grid::unit(1.2, "mm"),
  sm_anno_size = grid::unit(2.5, "mm"),
  group_anno_size = grid::unit(3, "mm"),
  mark_features = NULL,
  mark_label_by = NULL,
  mark_labels_gp = NULL,
  mark_annotation_width = NULL,
  mark_link_width = grid::unit(5, "mm"),
  name = NULL,
  ht_params = list(),
  width = NULL,
  height = NULL,
  verbose = TRUE
) {
  heatmap_border <- heatmap_border %||% border
  cell_annotation_border <- cell_annotation_border %||% border
  feature_annotation_border <- feature_annotation_border %||% border
  heatmap_border_color <- heatmap_border_color(
    heatmap_border,
    heatmap_border_palcolor
  )
  cell_annotation_border_color <- heatmap_border_color(
    cell_annotation_border,
    cell_annotation_border_palcolor
  )
  feature_annotation_border_color <- heatmap_border_color(
    feature_annotation_border,
    feature_annotation_border_palcolor
  )
  heatmap_border_size <- heatmap_border_size_value(heatmap_border_size)
  cell_annotation_border_size <- heatmap_border_size_value(cell_annotation_border_size)
  feature_annotation_border_size <- heatmap_border_size_value(feature_annotation_border_size)
  heatmap_border <- heatmap_border_enabled(heatmap_border)
  cell_annotation_border <- heatmap_border_enabled(cell_annotation_border)
  feature_annotation_border <- heatmap_border_enabled(feature_annotation_border)

  label_by <- match.arg(label_by)
  mark_label_by <- mark_label_by %||% label_by
  mark_label_by <- match.arg(mark_label_by, c("module", "reaction", "module_reaction"))
  if (!is.null(features) && !is.null(modules)) {
    log_message(
      "{.arg features} and {.arg modules} cannot both be supplied",
      message_type = "error"
    )
  }
  features <- features %||% modules
  if (is.null(width)) {
    width <- if (label_by == "module") 7.2 else 10.8
  }
  if (is.null(height)) {
    height <- 11.5
  }
  mat <- scfea_get_assay_matrix(srt, assay = assay, layer = layer)
  available_features <- rownames(mat)
  module_info_all <- scfea_get_module_info(srt, assay = assay)
  features_use <- scfea_resolve_module_features(
    features = features,
    module_info = module_info_all,
    available_features = available_features,
    verbose = verbose
  )
  mat <- mat[features_use, , drop = FALSE]
  module_info <- module_info_all[features_use, , drop = FALSE]

  group_df <- scfea_group_cells(srt, group.by = group.by, cells = colnames(mat))
  group_levels <- unique(as.character(group_df$group))
  agg_mat <- vapply(
    group_levels,
    function(group_i) {
      cells_i <- group_df$cell[group_df$group == group_i]
      Matrix::rowMeans(mat[, cells_i, drop = FALSE])
    },
    numeric(nrow(mat))
  )
  if (is.null(dim(agg_mat))) {
    agg_mat <- matrix(
      agg_mat,
      nrow = nrow(mat),
      dimnames = list(rownames(mat), group_levels)
    )
  } else {
    rownames(agg_mat) <- rownames(mat)
  }

  plot_mat <- agg_mat
  if (isTRUE(scale_rows)) {
    plot_mat <- t(scale(t(plot_mat)))
    plot_mat[is.na(plot_mat)] <- 0
    if (!is.null(heatmap_limit) && is.finite(heatmap_limit)) {
      plot_mat[plot_mat > heatmap_limit] <- heatmap_limit
      plot_mat[plot_mat < -heatmap_limit] <- -heatmap_limit
    }
  }
  color_limit <- heatmap_limit
  if (is.null(color_limit) || !is.finite(color_limit)) {
    color_limit <- max(abs(plot_mat), na.rm = TRUE)
  }
  if (!is.finite(color_limit) || color_limit == 0) {
    color_limit <- 1
  }

  row_labels <- scfea_module_labels(module_info, label_by = label_by)
  names(row_labels) <- rownames(plot_mat)
  name <- name %||% if (isTRUE(scale_rows)) "z-score" else "flux"

  if (is.numeric(heatmap_column_width)) {
    heatmap_column_width <- grid::unit(heatmap_column_width, "mm")
  }
  heatmap_body_width <- NULL
  if (!is.null(heatmap_column_width)) {
    heatmap_body_width <- heatmap_column_width * ncol(plot_mat)
  }

  left_annotation <- NULL
  row_split <- NULL
  if (isTRUE(add_sm_anno)) {
    sm_anno <- module_info$SM_anno
    sm_anno[is.na(sm_anno)] <- "Unknown"
    sm_anno <- factor(sm_anno, levels = unique(sm_anno))
    sm_cols <- palette_colors(
      levels(sm_anno),
      palette = feature_split_palette,
      palcolor = feature_split_palcolor
    )
    left_annotation <- ComplexHeatmap::rowAnnotation(
      SM_anno = sm_anno,
      col = list(SM_anno = sm_cols),
      show_annotation_name = FALSE,
      simple_anno_size = sm_anno_size,
      border = feature_annotation_border,
      gp = heatmap_border_gp(
        feature_annotation_border,
        feature_annotation_border_color
      ),
      annotation_legend_param = list(
        SM_anno = list(
          title = "SM_anno",
          labels_gp = legend_labels_gp,
          title_gp = legend_title_gp
        )
      )
    )
    row_split <- sm_anno
  }

  top_annotation <- NULL
  if (!is.null(group.by)) {
    group_cols <- palette_colors(
      group_levels,
      palette = group_palette,
      palcolor = group_palcolor
    )
    group_df <- data.frame(group = colnames(plot_mat), check.names = FALSE)
    colnames(group_df) <- group.by
    top_annotation <- ComplexHeatmap::HeatmapAnnotation(
      df = group_df,
      col = stats::setNames(list(group_cols), group.by),
      show_annotation_name = FALSE,
      simple_anno_size = group_anno_size,
      border = cell_annotation_border,
      gp = heatmap_border_gp(
        cell_annotation_border,
        cell_annotation_border_color
      ),
      annotation_legend_param = stats::setNames(
        list(list(
          title = group.by,
          labels_gp = legend_labels_gp,
          title_gp = legend_title_gp
        )),
        group.by
      )
    )
  }

  right_annotation <- NULL
  if (!is.null(mark_features)) {
    mark_features_use <- scfea_resolve_module_features(
      features = mark_features,
      module_info = module_info_all,
      available_features = available_features,
      verbose = verbose
    )
    mark_features_use <- intersect(mark_features_use, rownames(plot_mat))
    if (length(mark_features_use) == 0) {
      log_message(
        "No requested {.arg mark_features} are present in the plotted scFEA modules",
        message_type = "warning",
        verbose = verbose
      )
    } else {
      mark_row_labels <- scfea_module_labels(module_info, label_by = mark_label_by)
      names(mark_row_labels) <- rownames(plot_mat)
      mark_labels_gp <- mark_labels_gp %||% grid::gpar(
        fontsize = if (mark_label_by == "module") 8 else 6.5
      )
      if (is.null(mark_annotation_width)) {
        mark_annotation_width <- grid::unit(if (mark_label_by == "module") 20 else 65, "mm")
      } else if (is.numeric(mark_annotation_width)) {
        mark_annotation_width <- grid::unit(mark_annotation_width, "mm")
      }
      right_annotation <- ComplexHeatmap::rowAnnotation(
        mark = ComplexHeatmap::anno_mark(
          at = match(mark_features_use, rownames(plot_mat)),
          labels = unname(mark_row_labels[mark_features_use]),
          labels_gp = mark_labels_gp,
          lines_gp = grid::gpar(col = "black", lwd = 0.7),
          padding = grid::unit(0.8, "mm"),
          link_width = mark_link_width
        ),
        show_annotation_name = FALSE,
        width = mark_annotation_width
      )
    }
  }

  colors <- circlize::colorRamp2(
    seq(-color_limit, color_limit, length = 100),
    palette_colors(palette = heatmap_palette, palcolor = heatmap_palcolor)
  )
  ht_args <- list(
    matrix = plot_mat,
    name = name,
    col = colors,
    width = heatmap_body_width,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    row_labels = row_labels,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    use_raster = use_raster,
    raster_by_magick = raster_by_magick,
    column_names_rot = column_names_rot,
    column_names_gp = column_names_gp,
    column_title_gp = column_title_gp,
    left_annotation = left_annotation,
    top_annotation = top_annotation,
    right_annotation = right_annotation,
    row_split = row_split,
    row_title_rot = sm_anno_label_rot,
    row_title_gp = row_title_gp,
    row_gap = row_gap,
    border = heatmap_border,
    column_title = group.by %||% "All cells",
    heatmap_legend_param = list(
      legend_height = grid::unit(5, "cm"),
      title_position = "leftcenter-rot",
      title_gp = legend_title_gp,
      labels_gp = legend_labels_gp
    )
  )
  ht_args <- utils::modifyList(ht_args, ht_params)
  ht <- do.call(ComplexHeatmap::Heatmap, ht_args)
  attr(ht, "width") <- width
  attr(ht, "height") <- height
  ht
}

#' Plot scFEA flux Cohen's d volcano plots
#'
#' @param srt A Seurat object returned by [RunscFEA()].
#' @param group.by Metadata column defining groups. If `ident.1` and `ident.2`
#' are both `NULL`, each group is compared against all remaining cells.
#' @param ident.1,ident.2 Group names to compare. Cohen's d is
#' `mean(ident.1) - mean(ident.2)` divided by pooled standard deviation. If
#' both are `NULL`, one-vs-rest contrasts are generated automatically.
#' @param assay Flux assay name.
#' @param layer Flux assay layer.
#' @param label_by Label mode for significant pathway modules.
#' @param pathways Optional `SM_anno` classes to plot.
#' @param p_adj_cutoff Adjusted p-value cutoff for highlighting, labels, and
#' the horizontal threshold line.
#' @param cohen_cutoff Absolute Cohen's d cutoff for highlighting, labels, and
#' the vertical threshold lines.
#' @param combine Whether to combine pathway plots into 2 x 2 paged panels.
#' If `FALSE`, returns one plot per pathway.
#' @param width,height Suggested export size. Stored as attributes on the
#' returned plot object.
#' @param verbose Whether to print messages.
#'
#' @return For a single contrast, a paged list of `ggplot` objects when
#' `combine = TRUE`, otherwise a named list of pathway `ggplot` objects.
#' Statistics are stored in the `"data"` attribute. For automatic one-vs-rest
#' mode, a named list of single-contrast results is returned and combined
#' statistics are stored in the outer `"data"` attribute.
#' @export
scFEAVolcanoPlot <- function(
  srt,
  group.by,
  ident.1 = NULL,
  ident.2 = NULL,
  assay = "scFEAflux",
  layer = "data",
  label_by = c("reaction", "module", "module_reaction"),
  pathways = NULL,
  p_adj_cutoff = 1e-3,
  cohen_cutoff = 0.2,
  combine = TRUE,
  width = 12,
  height = 10.4,
  verbose = TRUE
) {
  label_by <- match.arg(label_by)
  mat <- scfea_get_assay_matrix(srt, assay = assay, layer = layer)
  module_info <- scfea_get_module_info(srt, assay = assay)
  module_info <- module_info[rownames(mat), , drop = FALSE]
  contrasts <- scfea_prepare_contrasts(
    srt = srt,
    group.by = group.by,
    cells = colnames(mat),
    ident.1 = ident.1,
    ident.2 = ident.2
  )

  if (isTRUE(contrasts$automatic)) {
    results <- lapply(
      contrasts$items,
      function(contrast_i) {
        scfea_volcano_plot_one(
          mat = mat,
          module_info = module_info,
          group_df = contrast_i$group_df,
          ident.1 = contrast_i$ident.1,
          ident.2 = contrast_i$ident.2,
          label_by = label_by,
          pathways = pathways,
          p_adj_cutoff = p_adj_cutoff,
          cohen_cutoff = cohen_cutoff,
          combine = combine,
          width = width,
          height = height,
          verbose = verbose
        )
      }
    )
    names(results) <- names(contrasts$items)
    attr(results, "data") <- scfea_bind_contrast_data(results, contrasts$items)
    attr(results, "width") <- width
    attr(results, "height") <- height
    return(results)
  }

  contrast_i <- contrasts$items[[1]]
  scfea_volcano_plot_one(
    mat = mat,
    module_info = module_info,
    group_df = contrast_i$group_df,
    ident.1 = contrast_i$ident.1,
    ident.2 = contrast_i$ident.2,
    label_by = label_by,
    pathways = pathways,
    p_adj_cutoff = p_adj_cutoff,
    cohen_cutoff = cohen_cutoff,
    combine = combine,
    width = width,
    height = height,
    verbose = verbose
  )
}

scfea_volcano_plot_one <- function(
  mat,
  module_info,
  group_df,
  ident.1,
  ident.2,
  label_by,
  pathways,
  p_adj_cutoff,
  cohen_cutoff,
  combine,
  width,
  height,
  verbose
) {
  stats <- scfea_compare_features(
    mat = mat,
    group_df = group_df,
    ident.1 = ident.1,
    ident.2 = ident.2
  )
  stats <- cbind(
    stats,
    module_info[
      stats$feature_id,
      setdiff(colnames(module_info), "feature_id"),
      drop = FALSE
    ]
  )
  stats$plot_label <- scfea_module_labels(stats, label_by = label_by)
  stats$neg_log10_padj <- -log10(pmax(stats$p_val_adj, .Machine$double.xmin))

  padj_line <- -log10(p_adj_cutoff)
  finite_cohen <- stats$cohen_d[is.finite(stats$cohen_d)]
  finite_padj <- stats$neg_log10_padj[is.finite(stats$neg_log10_padj)]
  x_breaks <- sort(unique(c(
    if (length(finite_cohen) > 0) pretty(finite_cohen) else numeric(),
    -cohen_cutoff,
    cohen_cutoff
  )))
  y_breaks <- sort(unique(c(
    if (length(finite_padj) > 0) pretty(c(finite_padj, padj_line)) else numeric(),
    padj_line
  )))
  threshold_label <- function(x) {
    formatC(x, format = "fg", digits = 4)
  }

  pathways <- pathways %||% unique(stats$SM_anno)
  pathways <- intersect(pathways, unique(stats$SM_anno))
  if (length(pathways) == 0) {
    log_message(
      "No {.arg pathways} match the scFEA module annotation",
      message_type = "error"
    )
  }

  check_r("cowplot", verbose = FALSE)
  plots <- lapply(
    pathways,
    function(pathway_i) {
      plot_df <- stats
      in_pathway <- !is.na(plot_df$SM_anno) & plot_df$SM_anno == pathway_i
      pathway_up_df <- plot_df[
        in_pathway & !is.na(plot_df$cohen_d) & plot_df$cohen_d > 0,
        ,
        drop = FALSE
      ]
      pathway_down_df <- plot_df[
        in_pathway & !is.na(plot_df$cohen_d) & plot_df$cohen_d < 0,
        ,
        drop = FALSE
      ]
      label_df <- plot_df[
        in_pathway &
          !is.na(plot_df$p_val_adj) &
          plot_df$p_val_adj <= p_adj_cutoff &
          abs(plot_df$cohen_d) >= cohen_cutoff,
        ,
        drop = FALSE
      ]

      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = .data$cohen_d, y = .data$neg_log10_padj)
      ) +
        ggplot2::geom_hline(
          yintercept = padj_line,
          color = "#999999",
          linetype = "dashed",
          linewidth = 0.8
        ) +
        ggplot2::geom_vline(
          xintercept = c(-cohen_cutoff, cohen_cutoff),
          color = "#999999",
          linetype = "dashed",
          linewidth = 0.8
        ) +
        ggplot2::geom_point(
          size = 1.2,
          color = "grey50",
          alpha = 0.8,
          na.rm = TRUE
        ) +
        ggplot2::geom_point(
          data = pathway_up_df,
          size = 1.2,
          shape = 16,
          color = "red",
          alpha = 0.8,
          na.rm = TRUE
        ) +
        ggplot2::geom_point(
          data = pathway_down_df,
          size = 1.2,
          shape = 16,
          color = "blue4",
          alpha = 0.8,
          na.rm = TRUE
        ) +
        ggplot2::geom_point(
          data = label_df,
          size = 3,
          shape = 16,
          color = "green4",
          alpha = 0.8,
          na.rm = TRUE
        ) +
        ggrepel::geom_text_repel(
          data = label_df,
          ggplot2::aes(label = .data$plot_label),
          color = "black",
          size = 4,
          arrow = grid::arrow(
            ends = "first",
            length = grid::unit(0.01, "npc")
          ),
          box.padding = 0.2,
          point.padding = 0.3,
          segment.color = "black",
          segment.size = 0.3,
          force = 1,
          max.iter = 3000,
          max.overlaps = Inf,
          na.rm = TRUE
        ) +
        ggplot2::labs(
          x = "Cohen's D",
          y = "-Log10(P-value)",
          title = paste0(ident.1, " vs ", ident.2, "\n", pathway_i)
        ) +
        ggplot2::scale_x_continuous(
          breaks = x_breaks,
          labels = threshold_label
        ) +
        ggplot2::scale_y_continuous(
          breaks = y_breaks,
          labels = threshold_label
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(
            linewidth = 0.8,
            colour = "black"
          ),
          axis.title = ggplot2::element_text(size = 14),
          axis.text = ggplot2::element_text(size = 8, color = "black"),
          plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
          plot.margin = grid::unit(c(0, 1, 2, 1), "cm")
        )

      cowplot::ggdraw(xlim = c(0, 1), ylim = c(0, 1.1)) +
        cowplot::draw_plot(p, x = 0, y = 0, width = 1, height = 1) +
        cowplot::draw_line(
          x = c(0.55, 0.75),
          y = c(0.1, 0.1),
          lineend = "round",
          size = 1,
          col = "#B51F29",
          arrow = grid::arrow(
            angle = 15,
            length = grid::unit(0.1, "inches"),
            type = "closed"
          )
        ) +
        cowplot::draw_line(
          x = c(0.25, 0.45),
          y = c(0.1, 0.1),
          lineend = "round",
          size = 1,
          col = "#006699",
          arrow = grid::arrow(
            angle = 15,
            length = grid::unit(0.1, "inches"),
            type = "closed",
            ends = "first"
          )
        ) +
        cowplot::draw_text(
          paste0("Activate in\n", ident.1),
          size = 12,
          x = 0.88,
          y = 0.1,
          color = "black",
          fontface = "italic"
        ) +
        cowplot::draw_text(
          paste0("Activate in\n", ident.2),
          size = 12,
          x = 0.15,
          y = 0.1,
          color = "black",
          fontface = "italic"
        )
    }
  )
  names(plots) <- pathways

  if (isTRUE(combine)) {
    split_four <- split(plots, ceiling(seq_along(plots) / 4))
    pages <- lapply(
      split_four,
      function(plot_group) {
        cowplot::plot_grid(
          plotlist = plot_group,
          ncol = 2,
          nrow = 2,
          align = "hv"
        )
      }
    )
    attr(pages, "data") <- stats
    attr(pages, "plots") <- plots
    attr(pages, "width") <- width
    attr(pages, "height") <- height
    return(pages)
  }

  attr(plots, "data") <- stats
  attr(plots, "width") <- width
  attr(plots, "height") <- height
  plots
}

#' Plot scFEA metabolite balance changes
#'
#' @param srt A Seurat object returned by [RunscFEA()].
#' @param group.by Metadata column defining groups. If `ident.1` and `ident.2`
#' are both `NULL`, each group is compared against all remaining cells.
#' @param ident.1,ident.2 Group names to compare. Difference is
#' `mean(ident.1) - mean(ident.2)`. If both are `NULL`, one-vs-rest contrasts
#' are generated automatically.
#' @param assay Balance assay name.
#' @param layer Balance assay layer.
#' @param top_n If `NULL`, plot all metabolites. If positive, plot the top
#' `top_n` increased and top `top_n` decreased metabolites after the
#' adjusted-p-value filter.
#' @param p_adj_cutoff Adjusted p-value cutoff used when `top_n` is not `NULL`.
#' @param title Optional plot title.
#'
#' @return For a single contrast, a list with `plot` and plotted `data`. For
#' automatic one-vs-rest mode, a named list of single-contrast results is
#' returned and combined statistics are stored in the outer `"data"` attribute.
#' @export
scFEABalanceBarPlot <- function(
  srt,
  group.by,
  ident.1 = NULL,
  ident.2 = NULL,
  assay = "scFEAbalance",
  layer = "data",
  top_n = NULL,
  p_adj_cutoff = 0.01,
  title = NULL
) {
  mat <- scfea_get_assay_matrix(srt, assay = assay, layer = layer)
  contrasts <- scfea_prepare_contrasts(
    srt = srt,
    group.by = group.by,
    cells = colnames(mat),
    ident.1 = ident.1,
    ident.2 = ident.2
  )

  if (isTRUE(contrasts$automatic)) {
    results <- lapply(
      contrasts$items,
      function(contrast_i) {
        scfea_balance_plot_one(
          mat = mat,
          group_df = contrast_i$group_df,
          ident.1 = contrast_i$ident.1,
          ident.2 = contrast_i$ident.2,
          top_n = top_n,
          p_adj_cutoff = p_adj_cutoff,
          title = title
        )
      }
    )
    names(results) <- names(contrasts$items)
    attr(results, "data") <- scfea_bind_contrast_data(
      results,
      contrasts$items,
      data_fun = function(x) x$data
    )
    return(results)
  }

  contrast_i <- contrasts$items[[1]]
  scfea_balance_plot_one(
    mat = mat,
    group_df = contrast_i$group_df,
    ident.1 = contrast_i$ident.1,
    ident.2 = contrast_i$ident.2,
    top_n = top_n,
    p_adj_cutoff = p_adj_cutoff,
    title = title
  )
}

scfea_balance_plot_one <- function(
  mat,
  group_df,
  ident.1,
  ident.2,
  top_n,
  p_adj_cutoff,
  title
) {
  stats <- scfea_compare_features(
    mat = mat,
    group_df = group_df,
    ident.1 = ident.1,
    ident.2 = ident.2
  )
  stats$compound_name <- stats$feature_id
  stats$direction <- ifelse(stats$mean_diff >= 0, ident.1, ident.2)

  plot_df <- stats
  if (!is.null(top_n)) {
    if (!is.numeric(top_n) || length(top_n) != 1 || top_n < 1) {
      log_message(
        "{.arg top_n} must be NULL or a positive number",
        message_type = "error"
      )
    }
    sig_df <- plot_df[!is.na(plot_df$p_val_adj) & plot_df$p_val_adj <= p_adj_cutoff, , drop = FALSE]
    if (nrow(sig_df) == 0) {
      log_message(
        "No metabolites passed {.arg p_adj_cutoff}; using all metabolites for top_n selection",
        message_type = "warning"
      )
      sig_df <- plot_df
    }
    up_df <- sig_df[sig_df$mean_diff > 0, , drop = FALSE]
    up_df <- up_df[order(up_df$mean_diff, decreasing = TRUE), , drop = FALSE]
    down_df <- sig_df[sig_df$mean_diff < 0, , drop = FALSE]
    down_df <- down_df[order(down_df$mean_diff, decreasing = FALSE), , drop = FALSE]
    plot_df <- rbind(
      utils::head(up_df, top_n),
      utils::head(down_df, top_n)
    )
  }

  if (nrow(plot_df) == 0) {
    log_message(
      "No metabolites available to plot after {.arg top_n} selection",
      message_type = "warning"
    )
  }

  plot_df <- plot_df[order(plot_df$mean_diff, decreasing = FALSE), , drop = FALSE]
  plot_df$compound_name <- factor(plot_df$compound_name, levels = plot_df$compound_name)
  title <- title %||% paste0(ident.1, " vs ", ident.2, " scFEA balance")

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$compound_name,
      y = .data$mean_diff,
      fill = .data$direction
    )
  ) +
    ggplot2::geom_col(width = 0.78) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, color = "grey30") +
    ggplot2::scale_fill_manual(
      values = stats::setNames(c("#D73027", "#4575B4"), c(ident.1, ident.2))
    ) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = paste0("Mean balance difference (", ident.1, " - ", ident.2, ")"),
      fill = "Higher in"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "top"
    )

  list(plot = p, data = plot_df)
}

scfea_get_assay_matrix <- function(srt, assay, layer) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (!assay %in% names(srt@assays)) {
    log_message(
      "Assay {.val {assay}} was not found in the Seurat object",
      message_type = "error"
    )
  }
  mat <- GetAssayData5(srt, assay = assay, layer = layer)
  if (is.null(mat)) {
    log_message(
      "Layer {.val {layer}} was not found in assay {.val {assay}}",
      message_type = "error"
    )
  }
  if (!inherits(mat, c("dgCMatrix", "matrix", "Matrix"))) {
    mat <- as_matrix(mat)
  }
  mat
}

scfea_group_cells <- function(srt, group.by, cells) {
  if (is.null(group.by)) {
    return(data.frame(cell = cells, group = "All cells", stringsAsFactors = FALSE))
  }
  if (length(group.by) != 1 || !group.by %in% colnames(srt@meta.data)) {
    log_message(
      "{.arg group.by} must be one metadata column in {.code srt@meta.data}",
      message_type = "error"
    )
  }
  group <- srt@meta.data[cells, group.by, drop = TRUE]
  if (any(is.na(group))) {
    log_message(
      "{.arg group.by} contains NA values for selected cells",
      message_type = "error"
    )
  }
  data.frame(cell = cells, group = as.character(group), stringsAsFactors = FALSE)
}

scfea_prepare_contrasts <- function(srt, group.by, cells, ident.1, ident.2) {
  if (missing(group.by) || is.null(group.by)) {
    log_message(
      "{.arg group.by} is required for this scFEA comparison",
      message_type = "error"
    )
  }
  if (is.null(ident.1) != is.null(ident.2)) {
    log_message(
      "{.arg ident.1} and {.arg ident.2} must be supplied together, or both left NULL for one-vs-rest",
      message_type = "error"
    )
  }

  group_df <- scfea_group_cells(srt, group.by = group.by, cells = cells)
  if (!is.null(ident.1) && !is.null(ident.2)) {
    ident.1 <- as.character(ident.1)
    ident.2 <- as.character(ident.2)
    contrast <- paste0(ident.1, "_vs_", ident.2)
    return(list(
      automatic = FALSE,
      items = list(
        contrast = list(
          contrast = contrast,
          ident.1 = ident.1,
          ident.2 = ident.2,
          group_df = group_df
        )
      )
    ))
  }

  group_raw <- srt@meta.data[cells, group.by, drop = TRUE]
  group_levels <- if (is.factor(group_raw)) {
    levels(group_raw)[levels(group_raw) %in% as.character(group_raw)]
  } else {
    unique(as.character(group_raw))
  }
  if (length(group_levels) < 2) {
    log_message(
      "{.arg group.by} must contain at least two groups for automatic one-vs-rest comparisons",
      message_type = "error"
    )
  }

  items <- lapply(
    group_levels,
    function(group_i) {
      one_group_df <- group_df
      one_group_df$group <- ifelse(one_group_df$group == group_i, group_i, "rest")
      list(
        contrast = paste0(group_i, "_vs_rest"),
        ident.1 = group_i,
        ident.2 = "rest",
        group_df = one_group_df
      )
    }
  )
  names(items) <- group_levels
  list(automatic = TRUE, items = items)
}

scfea_bind_contrast_data <- function(results, contrasts, data_fun = function(x) attr(x, "data")) {
  rows <- lapply(
    names(results),
    function(contrast_name) {
      data_i <- data_fun(results[[contrast_name]])
      if (is.null(data_i)) {
        return(NULL)
      }
      data_i <- as.data.frame(data_i, stringsAsFactors = FALSE)
      data.frame(
        contrast = contrasts[[contrast_name]]$contrast,
        ident.1 = contrasts[[contrast_name]]$ident.1,
        ident.2 = contrasts[[contrast_name]]$ident.2,
        data_i,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }
  )
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

scfea_compare_features <- function(mat, srt = NULL, group.by = NULL, ident.1, ident.2, group_df = NULL) {
  if (is.null(group_df)) {
    if (missing(group.by) || is.null(group.by)) {
      log_message(
        "{.arg group.by} is required for this scFEA comparison",
        message_type = "error"
      )
    }
    group_df <- scfea_group_cells(srt, group.by = group.by, cells = colnames(mat))
  } else {
    if (!all(c("cell", "group") %in% colnames(group_df))) {
      log_message(
        "{.arg group_df} must contain {.val cell} and {.val group} columns",
        message_type = "error"
      )
    }
    group_df <- group_df[match(colnames(mat), group_df$cell), , drop = FALSE]
    if (any(is.na(group_df$cell))) {
      log_message(
        "{.arg group_df} must contain every cell in the scFEA matrix",
        message_type = "error"
      )
    }
    group_df$group <- as.character(group_df$group)
  }
  cells1 <- group_df$cell[group_df$group == as.character(ident.1)]
  cells2 <- group_df$cell[group_df$group == as.character(ident.2)]
  if (length(cells1) == 0 || length(cells2) == 0) {
    log_message(
      "{.arg ident.1} and {.arg ident.2} must both be present in {.arg group.by}",
      message_type = "error"
    )
  }

  x1 <- mat[, cells1, drop = FALSE]
  x2 <- mat[, cells2, drop = FALSE]
  mean1 <- Matrix::rowMeans(x1)
  mean2 <- Matrix::rowMeans(x2)
  sd1 <- MatrixGenerics::rowSds(x1)
  sd2 <- MatrixGenerics::rowSds(x2)
  pooled_sd <- sqrt(
    ((length(cells1) - 1) * sd1^2 + (length(cells2) - 1) * sd2^2) /
      (length(cells1) + length(cells2) - 2)
  )
  cohen_d <- (mean1 - mean2) / pooled_sd
  cohen_d[!is.finite(cohen_d)] <- 0

  p_val <- vapply(
    seq_len(nrow(mat)),
    function(i) {
      tryCatch(
        stats::wilcox.test(
          x = as.numeric(x1[i, ]),
          y = as.numeric(x2[i, ]),
          exact = FALSE
        )$p.value,
        error = function(e) NA_real_
      )
    },
    numeric(1)
  )

  data.frame(
    feature_id = rownames(mat),
    mean_1 = mean1,
    mean_2 = mean2,
    mean_diff = mean1 - mean2,
    cohen_d = cohen_d,
    p_val = p_val,
    p_val_adj = stats::p.adjust(p_val, method = "BH"),
    stringsAsFactors = FALSE,
    row.names = rownames(mat)
  )
}

scfea_files <- c(
  "Human_M168_information.symbols.csv",
  "cName_c70_m168.csv",
  "cName_complete_mouse_c70_m168.csv",
  "cmMat_c70_m168.csv",
  "cmMat_complete_mouse_c70_m168.csv",
  "module_gene_complete_mouse_m168.csv",
  "module_gene_m168.csv"
)

resolve_scfea_dir <- function(data_dir = NULL, verbose = TRUE) {
  if (!is.null(data_dir)) {
    if (!dir.exists(data_dir)) {
      log_message(
        "Directory {.path {data_dir}} does not exist",
        message_type = "error"
      )
    }
    missing_files <- scfea_files[!file.exists(file.path(data_dir, scfea_files))]
    if (length(missing_files) > 0) {
      log_message(
        "Missing scFEA resource files in {.path {data_dir}}: {.val {missing_files}}",
        message_type = "error"
      )
    }
    return(normalizePath(data_dir, mustWork = TRUE))
  }

  cache_dir <- file.path(
    tools::R_user_dir("scop", "data"),
    "scFEA"
  )
  if (dir.exists(cache_dir) && all(file.exists(file.path(cache_dir, scfea_files)))) {
    return(cache_dir)
  }

  cache_key <- list("scFEA", "M168", "2026-06-01")
  if (requireNamespace("R.cache", quietly = TRUE)) {
    cached <- R.cache::loadCache(key = cache_key)
    cached_dir <- if (!is.null(cached$data_dir)) {
      normalizePath(cached$data_dir, mustWork = FALSE)
    } else {
      NULL
    }
    cache_dir_norm <- normalizePath(cache_dir, mustWork = FALSE)
    if (
      !is.null(cached) &&
        identical(cached_dir, cache_dir_norm) &&
        dir.exists(cached_dir) &&
        all(file.exists(file.path(cached_dir, scfea_files)))
    ) {
      log_message(
        "Using scFEA data from datasets cache: {.path {cached_dir}}",
        verbose = verbose
      )
      return(cached_dir)
    }
  }

  log_message(
    "Downloading scFEA model data from datasets GitHub repository...",
    verbose = verbose
  )
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  old_timeout <- getOption("timeout")
  options(timeout = max(600, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  for (fname in scfea_files) {
    url <- paste0(
      "https://raw.githubusercontent.com/mengxu98/datasets/main/scFEA/",
      fname
    )
    dest <- file.path(cache_dir, fname)
    if (!file.exists(dest)) {
      log_message(
        "  Downloading {.path {fname}} ...",
        verbose = verbose
      )
      utils::download.file(
        url = url,
        destfile = dest,
        mode = "wb",
        quiet = !verbose
      )
    }
  }

  if (requireNamespace("R.cache", quietly = TRUE)) {
    R.cache::saveCache(
      list(
        data_dir = cache_dir,
        files = scfea_files,
        version = "2026-06-01"
      ),
      key = cache_key,
      comment = paste0(
        "2026-06-01 nterm:",
        length(scfea_files),
        "|scFEA-M168"
      )
    )
  }

  cache_dir
}

scfea_module_info <- function(data_dir = NULL, verbose = TRUE) {
  data_dir <- resolve_scfea_dir(data_dir = data_dir, verbose = verbose)
  info <- utils::read.csv(
    file.path(data_dir, "Human_M168_information.symbols.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  module_id <- info[[1]]
  if (!grepl("^M_", module_id[1])) {
    module_id <- rownames(info)
  }
  info$module_id <- module_id
  info <- info[grepl("^M_", info$module_id), , drop = FALSE]
  info$feature_id <- scfea_module_feature_id(info$module_id)
  info$module_label <- sub("_", "", info$module_id, fixed = TRUE)
  info$reaction_label <- paste0(
    info$Compound_IN_name,
    " -> ",
    info$Compound_OUT_name
  )
  info$SM_anno <- scfea_supermodule_labels()[as.character(info$Supermodule_id)]
  info <- info[, c(
    "module_id", "feature_id", "module_label", "reaction_label",
    "SM_anno", "Module_id", "Compound_IN_name", "Compound_IN_ID",
    "Compound_OUT_name", "Compound_OUT_ID", "Supermodule_id"
  )]
  rownames(info) <- info$feature_id
  info
}

scfea_get_module_info <- function(srt, assay = "scFEAflux") {
  module_info <- NULL
  if ("scFEA" %in% names(srt@tools) && !is.null(srt@tools[["scFEA"]]$module_info)) {
    module_info <- srt@tools[["scFEA"]]$module_info
  } else if ("ScFEA" %in% names(srt@tools) && !is.null(srt@tools[["ScFEA"]]$module_info)) {
    module_info <- srt@tools[["ScFEA"]]$module_info
  }
  if (is.null(module_info)) {
    module_info <- scfea_module_info(verbose = FALSE)
  }
  if (!all(rownames(module_info) %in% rownames(srt[[assay]]))) {
    module_info <- module_info[intersect(rownames(module_info), rownames(srt[[assay]])), , drop = FALSE]
  }
  rownames(module_info) <- module_info$feature_id
  module_info
}

scfea_compound_info <- function(species = c("human", "mouse"), data_dir = NULL, verbose = TRUE) {
  species <- match.arg(species)
  file <- scfea_species_files(species, data_dir = data_dir, verbose = verbose)[["cName_file"]]
  cname <- utils::read.csv(
    file,
    header = FALSE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  compound_info <- data.frame(
    compound_name = as.character(unlist(cname[1, ], use.names = FALSE)),
    compound_id = as.character(unlist(cname[2, ], use.names = FALSE)),
    stringsAsFactors = FALSE
  )
  rownames(compound_info) <- compound_info$compound_name
  compound_info
}

scfea_module_genes <- function(species = c("human", "mouse"), data_dir = NULL, verbose = TRUE) {
  species <- match.arg(species)
  file <- scfea_species_files(species, data_dir = data_dir, verbose = verbose)[["moduleGene_file"]]
  module_gene <- utils::read.csv(
    file,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  genes <- unique(as.character(unlist(module_gene, use.names = FALSE)))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  genes
}

scfea_species_files <- function(species = c("human", "mouse"), data_dir = NULL, verbose = TRUE) {
  species <- match.arg(species)
  data_dir <- resolve_scfea_dir(data_dir = data_dir, verbose = verbose)
  files <- list(
    human = list(
      moduleGene_file = file.path(data_dir, "module_gene_m168.csv"),
      stoichiometry_matrix = file.path(data_dir, "cmMat_c70_m168.csv"),
      cName_file = file.path(data_dir, "cName_c70_m168.csv")
    ),
    mouse = list(
      moduleGene_file = file.path(data_dir, "module_gene_complete_mouse_m168.csv"),
      stoichiometry_matrix = file.path(data_dir, "cmMat_complete_mouse_c70_m168.csv"),
      cName_file = file.path(data_dir, "cName_complete_mouse_c70_m168.csv")
    )
  )
  files[[species]]
}

scfea_module_feature_id <- function(module_id) {
  gsub("_", "-", module_id, fixed = TRUE)
}

scfea_normalize_module_query <- function(x) {
  x <- as.character(x)
  x <- gsub("\u2212|\u2013|\u2014", "-", x)
  x <- gsub("\u2192|\u21d2|\u27f6|\u27f9", "->", x)
  x <- gsub("\\s*->\\s*", " -> ", x)
  x <- gsub("\\s*:\\s*", ": ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

scfea_resolve_module_features <- function(
  features,
  module_info,
  available_features,
  verbose = TRUE
) {
  if (is.null(features)) {
    return(available_features)
  }
  features <- scfea_normalize_module_query(features)
  module_reaction_label <- paste0(module_info$module_label, ": ", module_info$reaction_label)
  maps <- list(
    feature_id = stats::setNames(module_info$feature_id, module_info$feature_id),
    module_id = stats::setNames(module_info$feature_id, module_info$module_id),
    module_label = stats::setNames(module_info$feature_id, module_info$module_label),
    reaction_label = stats::setNames(module_info$feature_id, module_info$reaction_label),
    module_reaction_label = stats::setNames(module_info$feature_id, module_reaction_label)
  )
  maps <- c(
    maps,
    lapply(
      maps,
      function(map_i) {
        stats::setNames(unname(map_i), scfea_normalize_module_query(names(map_i)))
      }
    )
  )
  resolved <- unlist(
    lapply(
      features,
      function(feature_i) {
        hit <- NULL
        for (map_i in maps) {
          if (feature_i %in% names(map_i)) {
            hit <- unname(map_i[[feature_i]])
            break
          }
        }
        hit
      }
    ),
    use.names = FALSE
  )
  valid_features <- unique(unlist(lapply(maps, names), use.names = FALSE))
  missing_features <- setdiff(features, valid_features)
  if (length(missing_features) > 0) {
    log_message(
      "Dropping unknown scFEA module features: {.val {missing_features}}",
      message_type = "warning",
      verbose = verbose
    )
  }
  resolved <- intersect(unique(resolved), available_features)
  if (length(resolved) == 0) {
    log_message(
      "No requested scFEA module features are available in the assay",
      message_type = "error"
    )
  }
  resolved
}

scfea_module_labels <- function(module_info, label_by = c("module", "reaction", "module_reaction")) {
  label_by <- match.arg(label_by)
  if (label_by == "module") {
    return(module_info$module_label)
  }
  if (label_by == "reaction") {
    return(module_info$reaction_label)
  }
  paste0(module_info$module_label, ": ", module_info$reaction_label)
}

scfea_py_dataframe_to_r <- function(py_df, expected_rows = NULL) {
  df <- reticulate::py_to_r(py_df)
  if (!is.data.frame(df)) {
    df <- as.data.frame(df, check.names = FALSE)
  }
  default_rows <- list(
    as.character(seq_len(nrow(df))),
    as.character(seq_len(nrow(df)) - 1L)
  )
  if (ncol(df) > 1) {
    first_col <- df[[1]]
    first_is_index <- !is.numeric(first_col) &&
      length(unique(as.character(first_col))) == nrow(df)
    if (
      first_is_index &&
        (identical(rownames(df), default_rows[[1]]) ||
          identical(rownames(df), default_rows[[2]]))
    ) {
      rownames(df) <- as.character(first_col)
      df <- df[, -1, drop = FALSE]
    }
  }
  if (!is.null(expected_rows) && (
    is.null(rownames(df)) ||
      any(rownames(df) %in% c("", NA)) ||
      identical(rownames(df), default_rows[[1]]) ||
      identical(rownames(df), default_rows[[2]])
  )) {
    rownames(df) <- expected_rows
  }
  row_ids <- rownames(df)
  df <- as.data.frame(
    lapply(df, function(x) as.numeric(as.character(x))),
    check.names = FALSE
  )
  rownames(df) <- row_ids
  if (!is.null(expected_rows) && length(expected_rows) == nrow(df)) {
    rownames(df) <- rownames(df) %||% expected_rows
    if (!all(expected_rows %in% rownames(df))) {
      rownames(df) <- expected_rows
    }
  }
  df
}

scfea_supermodule_labels <- function() {
  c(
    "1" = "Glycolysis_TCA_cycle",
    "2" = "Serine_metabolism",
    "3" = "Pentose_phosphate_pathway",
    "4" = "Fatty_acid_metabolism",
    "5" = "Aspartate_metabolism",
    "6" = "Beta_Alanine_metabolism",
    "7" = "Propionyl_CoA_metabolism",
    "8" = "Glutamate_metabolism",
    "9" = "BCAA_metabolism",
    "10" = "Urea_cycle",
    "11" = "Spermine_metabolism",
    "12" = "Transporters",
    "13" = "Hyaluronic_acid_synthesis",
    "14" = "Glycogen_synthesis",
    "15" = "Glycosaminoglycan_synthesis",
    "16" = "N_linked_Glycan_synthesis",
    "17" = "O_linked_Glycan_synthesis",
    "18" = "Sialic_acid_synthesis",
    "19" = "Glycan_synthesis",
    "20" = "Purine_synthesis",
    "21" = "Pyrimidine_synthesis",
    "22" = "Steroid_hormone_synthesis"
  )
}
