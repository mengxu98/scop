#' @title Plot dynamic features across pseudotime
#'
#' @md
#' @inheritParams RunStandardWorkflow
#' @inheritParams CellDimPlot
#' @inheritParams GroupHeatmap
#' @inheritParams RunDynamicFeatures
#' @inheritParams scop-params
#' @param features Features to use.
#' @param lineages Lineages to plot.
#' @param group_use Groups from `group.by` to
#' keep. If both `group_use` and `cells` are provided, their intersection will
#' be used.
#' @param fit.by Optional metadata column used to fit an independent trajectory
#' for each group over that group's observed pseudotime range. The fit is not
#' extrapolated beyond observed cells. This does not change the existing
#' point-coloring behavior of `group.by`. Groups with fewer than 10 cells or 10
#' unique pseudotime values are shown as raw points or rugs but are not fitted.
#' @param cells Cell names to use.
#' @param family GLM family. `NULL` chooses automatically.
#' @param exp_method Expression transform: `"log1p"`, `"raw"`, `"zscore"`, `"fc"`,
#' or `"log2fc"`. With `fit.by`, transform parameters are shared across groups
#' within each lineage.
#' @param lib_normalize Library-size normalize. Defaults to `TRUE` when `layer`
#' is counts.
#' @param libsize Library size for each cell.
#' @param compare_lineages,compare_features Compare lineages or features on one
#' plot. Grouped line or interval plots support up to six combined
#' lineage-feature series.
#' @param add_line,add_interval,line.size,line_palette,line_palcolor Fitted line and CI.
#' @param add_point,pt.size,point_palette,point_palcolor Overlay points.
#' @param add_rug Draw rugs.
#' @param flip,reverse Flip or reverse the x-axis.
#' @param x_order Order of x-axis values.
#' @param aspect.ratio Panel aspect ratio.
#'
#' @seealso [RunDynamicFeatures]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#' pancreas_sub <- RunSlingshot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP"
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   lineages = paste0("Lineage", 1:2),
#'   lineages_span = 0.1
#' )
#'
#' DynamicPlot(
#'   pancreas_sub,
#'   lineages = "Lineage1",
#'   features = c("Arxes1", "Ncoa2", "G2M_score"),
#'   group.by = "SubCellType",
#'   group_use = c("Ductal", "Beta"),
#'   compare_features = TRUE
#' )
#'
#' # Demonstration only: create two conditions spanning the same pseudotime.
#' lineage1_cells <- rownames(pancreas_sub@meta.data)[
#'   order(pancreas_sub$Lineage1, na.last = NA)
#' ]
#' pancreas_sub$condition <- NA_character_
#' pancreas_sub$condition[lineage1_cells] <- rep(
#'   c("control", "HUA"),
#'   length.out = length(lineage1_cells)
#' )
#' DynamicPlot(
#'   pancreas_sub,
#'   lineages = "Lineage1",
#'   features = "Vim",
#'   group.by = "condition",
#'   fit.by = "condition"
#' )
#'
#' DynamicPlot(
#'   pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   features = c("Arxes1", "Ncoa2", "G2M_score"),
#'   group.by = "SubCellType",
#'   compare_lineages = TRUE,
#'   compare_features = FALSE
#' )
#'
#' DynamicPlot(
#'   pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   features = c("Arxes1", "Ncoa2", "G2M_score"),
#'   group.by = "SubCellType",
#'   compare_lineages = FALSE,
#'   compare_features = FALSE
#' )
DynamicPlot <- function(
  srt,
  lineages,
  features,
  group.by = NULL,
  group_use = NULL,
  cells = NULL,
  layer = "counts",
  assay = NULL,
  family = NULL,
  exp_method = c(
    "log1p",
    "raw",
    "zscore",
    "fc",
    "log2fc"
  ),
  lib_normalize = identical(layer, "counts"),
  libsize = NULL,
  compare_lineages = TRUE,
  compare_features = FALSE,
  add_line = TRUE,
  add_interval = TRUE,
  line.size = 1,
  line_palette = "Dark2",
  line_palcolor = NULL,
  add_point = TRUE,
  pt.size = 1,
  point_palette = "Chinese",
  point_palcolor = NULL,
  add_rug = TRUE,
  flip = FALSE,
  reverse = FALSE,
  x_order = c("value", "rank"),
  aspect.ratio = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  cores = 1,
  verbose = TRUE,
  seed = 11,
  fit.by = NULL
) {
  set.seed(seed)

  check_r("MatrixGenerics", verbose = FALSE)
  x_order <- match.arg(x_order)
  if (!is.null(group.by)) {
    group_missing <- setdiff(group.by, colnames(srt@meta.data))
    if (length(group_missing) > 0) {
      log_message(
        "{.val {group_missing}} is not in the meta.data of srt object",
        message_type = "error"
      )
    }
  }
  if (!is.null(fit.by)) {
    if (length(fit.by) != 1) {
      log_message(
        "{.arg fit.by} must be one metadata column.",
        message_type = "error"
      )
    }
    if (!fit.by %in% colnames(srt@meta.data)) {
      log_message(
        "{.val {fit.by}} is not in the meta.data of srt object",
        message_type = "error"
      )
    }
  }
  if (!is.null(group_use)) {
    if (is.null(group.by)) {
      log_message(
        "{.arg group_use} requires {.arg group.by}.",
        message_type = "error"
      )
    }
    if (length(group.by) > 1) {
      log_message(
        "{.arg group_use} supports one {.arg group.by} column.",
        message_type = "error"
      )
    }
    group_use <- unique(as.character(group_use))
    group_use <- group_use[nzchar(group_use)]
    group_cells <- rownames(srt@meta.data)[
      as.character(srt@meta.data[[group.by]]) %in% group_use
    ]
    if (length(group_cells) == 0) {
      log_message(
        "No cells found for {.arg group_use} in {.val {group.by}}.",
        message_type = "error"
      )
    }
    cells <- if (is.null(cells)) {
      group_cells
    } else {
      intersect(cells, group_cells)
    }
    if (length(cells) == 0) {
      log_message(
        "No cells remain after intersecting {.arg cells} and {.arg group_use}.",
        message_type = "error"
      )
    }
  }

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

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!is.null(family) && length(family) > 1 && is.null(names(family))) {
    if (length(family) != length(features)) {
      log_message(
        "{.arg family} must be one character or a vector of the same length as features",
        message_type = "error"
      )
    }
    names(family) <- features
  }
  internal_cols <- make.unique(c(
    colnames(srt@meta.data),
    ".scop_fit_group",
    ".scop_fit_series",
    ".scop_fit_path"
  ))
  internal_cols <- internal_cols[
    ncol(srt@meta.data) + seq_len(3L)
  ]
  fit_group_col <- internal_cols[[1]]
  fit_series_col <- internal_cols[[2]]
  fit_path_col <- internal_cols[[3]]
  gene <- features[features %in% rownames(srt@assays[[assay]])]
  meta <- features[features %in% colnames(srt@meta.data)]
  features <- c(gene, meta)
  if (length(features) == 0) {
    log_message(
      "No feature found in the srt object.",
      message_type = "error"
    )
  }

  cell_union <- c()
  raw_matrix_list <- list()
  fitted_matrix_list <- list()
  upr_matrix_list <- list()
  lwr_matrix_list <- list()
  fit_lineage_list <- list()
  fit_group_list <- list()
  if (is.null(fit.by)) {
    fit_groups <- NA_character_
    fit_group_levels <- "all"
  } else {
    fit_cells <- rownames(srt@meta.data)
    if (!is.null(cells)) {
      fit_cells <- intersect(fit_cells, cells)
    }
    fit_values <- srt@meta.data[fit_cells, fit.by, drop = TRUE]
    fit_values_present <- as.character(fit_values[!is.na(fit_values)])
    fit_values_present <- fit_values_present[nzchar(fit_values_present)]
    if (is.factor(srt@meta.data[[fit.by]])) {
      fit_group_levels <- levels(srt@meta.data[[fit.by]])
      fit_groups <- fit_group_levels[fit_group_levels %in% fit_values_present]
    } else {
      fit_groups <- unique(fit_values_present)
      fit_group_levels <- fit_groups
    }
    if (length(fit_groups) == 0) {
      log_message(
        "No non-missing groups found for {.arg fit.by} in the selected cells. Raw points and rugs will be retained without fitted trajectories.",
        message_type = "warning",
        verbose = verbose
      )
    }
    raw_feature_matrix <- GetAssayData5(
      srt,
      assay = assay,
      layer = layer
    )[gene, , drop = FALSE]
    if (length(meta) > 0) {
      raw_feature_matrix <- rbind(
        raw_feature_matrix,
        Matrix::t(srt@meta.data[, meta, drop = FALSE])
      )
    }
    raw_feature_matrix <- raw_feature_matrix[features, , drop = FALSE]
  }
  add_group_matrices <- function(
    raw_matrix,
    fitted_matrix,
    upr_matrix,
    lwr_matrix,
    lineage,
    fit_group
  ) {
    fit_id <- paste0("fit", length(raw_matrix_list) + 1L)
    raw_matrix_list[[fit_id]] <<- raw_matrix
    fitted_matrix_list[[fit_id]] <<- fitted_matrix
    upr_matrix_list[[fit_id]] <<- upr_matrix
    lwr_matrix_list[[fit_id]] <<- lwr_matrix
    fit_lineage_list[[fit_id]] <<- lineage
    fit_group_list[[fit_id]] <<- fit_group
    cell_union <<- unique(c(cell_union, rownames(raw_matrix)))
  }
  empty_fit_matrix <- function(cell_names) {
    matrix(
      NA_real_,
      nrow = length(cell_names),
      ncol = length(features),
      dimnames = list(cell_names, features)
    )
  }
  subset_family <- function(features_fit) {
    if (is.null(family) || length(family) == 1) {
      return(family)
    }
    if (length(family) != length(features)) {
      log_message(
        "{.arg family} must be one character or a vector of the same length as features",
        message_type = "error"
      )
    }
    family_all <- family
    if (is.null(names(family_all))) {
      names(family_all) <- features
    }
    family_all[features_fit]
  }
  fit_libsize_all <- NULL
  feature_uses_libsize <- stats::setNames(
    rep(FALSE, length(features)),
    features
  )
  if (!is.null(fit.by) && layer == "counts" && length(gene) > 0) {
    counts_status <- CheckDataType(
      srt,
      assay = assay,
      layer = "counts",
      verbose = FALSE
    )
    family_all <- subset_family(features)
    if (is.null(family_all)) {
      family_all <- stats::setNames(
        rep("gaussian", length(features)),
        features
      )
      if (counts_status == "raw_counts") {
        family_all[gene] <- "nb"
      }
    } else if (length(family_all) == 1) {
      family_all <- stats::setNames(
        rep(family_all, length(features)),
        features
      )
    }
    feature_uses_libsize[gene] <- family_all[gene] != "gaussian"
    if (any(feature_uses_libsize)) {
      if (is.null(libsize)) {
        fit_libsize_all <- if (counts_status == "raw_counts") {
          Matrix::colSums(
            GetAssayData5(srt, assay = assay, layer = "counts")
          )
        } else {
          stats::setNames(rep(1, ncol(srt)), colnames(srt))
        }
      } else {
        if (!length(libsize) %in% c(1, ncol(srt))) {
          log_message(
            "{.arg libsize} must be length of 1 or the number of cells",
            message_type = "error"
          )
        }
        fit_libsize_all <- if (length(libsize) == 1) {
          stats::setNames(rep(libsize, ncol(srt)), colnames(srt))
        } else {
          stats::setNames(as.numeric(libsize), colnames(srt))
        }
      }
    }
  }
  for (l in lineages) {
    if (!is.null(fit.by)) {
      for (g in fit_groups) {
        group_cells <- fit_cells[
          !is.na(fit_values) & as.character(fit_values) == g
        ]
        lineage_values <- srt@meta.data[group_cells, l, drop = TRUE]
        valid_cells <- group_cells[is.finite(lineage_values)]
        valid_pseudotime <- lineage_values[is.finite(lineage_values)]
        if (length(valid_cells) == 0) {
          log_message(
            "Skip {.val {g}} in {.val {l}}: no cells have finite pseudotime values.",
            message_type = "warning",
            verbose = verbose
          )
          next
        }
        if (length(valid_cells) < 10 || length(unique(valid_pseudotime)) < 10) {
          ordered_cells <- valid_cells[order(valid_pseudotime)]
          raw_matrix <- t(as_matrix(
            raw_feature_matrix[, ordered_cells, drop = FALSE]
          ))
          empty_fit <- empty_fit_matrix(ordered_cells)
          add_group_matrices(
            raw_matrix,
            empty_fit,
            empty_fit,
            empty_fit,
            l,
            g
          )
          log_message(
            "Skip the fitted trajectory for {.val {g}} in {.val {l}}: the default GAM requires at least ten cells and ten unique pseudotime values. Raw points and rugs are retained.",
            message_type = "warning",
            verbose = verbose
          )
          next
        }

        feature_values <- as_matrix(
          raw_feature_matrix[, valid_cells, drop = FALSE]
        )
        feature_supported <- vapply(
          seq_along(features),
          function(i) {
            valid_feature <- is.finite(feature_values[i, ])
            if (feature_uses_libsize[[i]]) {
              libsize_group <- fit_libsize_all[valid_cells]
              valid_feature <- valid_feature &
                is.finite(libsize_group) & libsize_group > 0
            }
            sum(valid_feature) >= 10 &&
              length(unique(valid_pseudotime[valid_feature])) >= 10
          },
          logical(1)
        )
        features_fit <- features[feature_supported]
        features_raw_only <- features[!feature_supported]
        if (length(features_raw_only) > 0) {
          log_message(
            "Skip the fitted trajectory for {.val {features_raw_only}} in {.val {g}} / {.val {l}}: each feature requires at least ten finite usable observations and ten unique pseudotime values. Raw points and rugs are retained.",
            message_type = "warning",
            verbose = verbose
          )
        }
        ordered_cells <- valid_cells[order(valid_pseudotime)]
        if (length(features_fit) == 0) {
          raw_matrix <- t(as_matrix(
            raw_feature_matrix[, ordered_cells, drop = FALSE]
          ))
          empty_fit <- empty_fit_matrix(ordered_cells)
          add_group_matrices(
            raw_matrix,
            empty_fit,
            empty_fit,
            empty_fit,
            l,
            g
          )
          next
        }

        srt_tmp <- srt
        exclude_cells <- setdiff(rownames(srt_tmp@meta.data), valid_cells)
        srt_tmp@meta.data[exclude_cells, l] <- NA_real_
        srt_tmp <- RunDynamicFeatures(
          srt_tmp,
          lineages = l,
          features = features_fit,
          assay = assay,
          layer = layer,
          family = subset_family(features_fit),
          libsize = libsize,
          cores = cores,
          verbose = verbose
        )
        fit_result <- srt_tmp@tools[[paste0("DynamicFeatures_", l)]]
        fitted_cells <- rownames(fit_result[["raw_matrix"]])
        raw_matrix <- t(as_matrix(
          raw_feature_matrix[, fitted_cells, drop = FALSE]
        ))
        fitted_matrix <- empty_fit_matrix(fitted_cells)
        upr_matrix <- empty_fit_matrix(fitted_cells)
        lwr_matrix <- empty_fit_matrix(fitted_cells)
        fitted_matrix[, features_fit] <- as_matrix(
          fit_result[["fitted_matrix"]][, features_fit, drop = FALSE]
        )
        upr_matrix[, features_fit] <- as_matrix(
          fit_result[["upr_matrix"]][, features_fit, drop = FALSE]
        )
        lwr_matrix[, features_fit] <- as_matrix(
          fit_result[["lwr_matrix"]][, features_fit, drop = FALSE]
        )
        add_group_matrices(
          raw_matrix,
          fitted_matrix,
          upr_matrix,
          lwr_matrix,
          l,
          g
        )
      }
      unassigned_cells <- fit_cells[
        is.na(fit_values) | !nzchar(as.character(fit_values))
      ]
      unassigned_pseudotime <- srt@meta.data[
        unassigned_cells,
        l,
        drop = TRUE
      ]
      unassigned_cells <- unassigned_cells[is.finite(unassigned_pseudotime)]
      unassigned_pseudotime <- unassigned_pseudotime[
        is.finite(unassigned_pseudotime)
      ]
      if (length(unassigned_cells) > 0) {
        ordered_cells <- unassigned_cells[order(unassigned_pseudotime)]
        raw_matrix <- t(as_matrix(
          raw_feature_matrix[, ordered_cells, drop = FALSE]
        ))
        empty_fit <- empty_fit_matrix(ordered_cells)
        add_group_matrices(
          raw_matrix,
          empty_fit,
          empty_fit,
          empty_fit,
          l,
          NA_character_
        )
      }
      next
    }

    features_exist <- c()
    raw_matrix <- NULL
    if (paste0("DynamicFeatures_", l) %in% names(srt@tools)) {
      raw_matrix <- srt@tools[[paste0("DynamicFeatures_", l)]][["raw_matrix"]][
        ,
        -1,
        drop = FALSE
      ]
      fitted_matrix <- srt@tools[[paste0("DynamicFeatures_", l)]][[
        "fitted_matrix"
      ]][, -1, drop = FALSE]
      upr_matrix <- srt@tools[[paste0("DynamicFeatures_", l)]][["upr_matrix"]][
        ,
        -1,
        drop = FALSE
      ]
      lwr_matrix <- srt@tools[[paste0("DynamicFeatures_", l)]][["lwr_matrix"]][
        ,
        -1,
        drop = FALSE
      ]
      features_exist <- colnames(raw_matrix)
    }
    feature_calcu <- features[!features %in% features_exist]
    if (length(feature_calcu) > 0) {
      srt_tmp <- RunDynamicFeatures(
        srt,
        lineages = l,
        features = feature_calcu,
        assay = assay,
        layer = layer,
        family = family,
        libsize = libsize,
        cores = cores,
        verbose = verbose
      )
      if (is.null(raw_matrix)) {
        raw_matrix <- fitted_matrix <- upr_matrix <- lwr_matrix <- matrix(
          NA,
          nrow = nrow(srt_tmp@tools[[paste0("DynamicFeatures_", l)]][[
            "raw_matrix"
          ]]),
          ncol = 0
        )
      }
      raw_matrix <- cbind(
        raw_matrix,
        srt_tmp@tools[[paste0("DynamicFeatures_", l)]][["raw_matrix"]][,
          feature_calcu,
          drop = FALSE
        ]
      )
      fitted_matrix <- cbind(
        fitted_matrix,
        srt_tmp@tools[[paste0("DynamicFeatures_", l)]][["fitted_matrix"]][,
          feature_calcu,
          drop = FALSE
        ]
      )
      upr_matrix <- cbind(
        upr_matrix,
        srt_tmp@tools[[paste0("DynamicFeatures_", l)]][["upr_matrix"]][,
          feature_calcu,
          drop = FALSE
        ]
      )
      lwr_matrix <- cbind(
        lwr_matrix,
        srt_tmp@tools[[paste0("DynamicFeatures_", l)]][["lwr_matrix"]][,
          feature_calcu,
          drop = FALSE
        ]
      )
    }
    raw_matrix_list[[l]] <- as_matrix(
      raw_matrix[, features, drop = FALSE]
    )
    fitted_matrix_list[[l]] <- as_matrix(
      fitted_matrix[, features, drop = FALSE]
    )
    upr_matrix_list[[l]] <- as_matrix(
      upr_matrix[, features, drop = FALSE]
    )
    lwr_matrix_list[[l]] <- as_matrix(
      lwr_matrix[, features, drop = FALSE]
    )
    fit_lineage_list[[l]] <- l
    fit_group_list[[l]] <- NA_character_
    cell_union <- unique(c(cell_union, rownames(raw_matrix)))
  }

  x_assign <- rowMeans(
    srt@meta.data[cell_union, lineages, drop = FALSE],
    na.rm = TRUE
  )
  cell_metadata <- cbind.data.frame(
    data.frame(row.names = cell_union),
    x_assign = x_assign,
    srt@meta.data[cell_union, lineages, drop = FALSE]
  )

  cell_order_list <- list()
  for (l in lineages) {
    cell_metadata_sub <- stats::na.omit(
      cell_metadata[, l, drop = FALSE]
    )
    cell_metadata_sub <- cell_metadata_sub[
      order(
        cell_metadata_sub[[l]],
        decreasing = FALSE
      ), ,
      drop = FALSE
    ]
    cell_order_list[[l]] <- paste0(rownames(cell_metadata_sub), l)
  }

  df_list <- list()
  y_libsize <- Matrix::colSums(
    GetAssayData5(
      srt,
      assay = assay,
      layer = "counts"
    )
  )
  if (!is.null(libsize)) {
    if (!length(libsize) %in% c(1, ncol(srt))) {
      log_message(
        "{.arg libsize} must be length of 1 or the number of cells",
        message_type = "error"
      )
    }
    libsize_all <- if (length(libsize) == 1) {
      stats::setNames(rep(libsize, ncol(srt)), colnames(srt))
    } else {
      stats::setNames(as.numeric(libsize), colnames(srt))
    }
    normalization_center <- stats::median(libsize_all)
  } else {
    libsize_all <- NULL
    normalization_center <- stats::median(y_libsize)
  }

  for (fit_id in names(raw_matrix_list)) {
    raw_matrix <- raw_matrix_list[[fit_id]]
    if (
      isTRUE(lib_normalize) &&
        length(gene) > 0 &&
        min(raw_matrix[, gene, drop = FALSE], na.rm = TRUE) >= 0
    ) {
      if (!is.null(libsize_all)) {
        libsize_use <- libsize_all[rownames(raw_matrix)]
      } else {
        libsize_use <- y_libsize[rownames(raw_matrix)]
        isfloat <- any(libsize_use %% 1 != 0, na.rm = TRUE)
        if (isTRUE(isfloat)) {
          libsize_use <- rep(1, length(libsize_use))
          log_message(
            "The values in the 'counts' layer are non-integer. Set the library size to 1.",
            message_type = "warning"
          )
        }
      }
      raw_matrix[, gene] <- raw_matrix[, gene, drop = FALSE] /
        libsize_use *
        normalization_center
    }
    raw_matrix_list[[fit_id]] <- raw_matrix
  }

  transform_reference_list <- stats::setNames(
    lapply(lineages, function(l) {
      fit_ids <- names(fit_lineage_list)[
        unlist(fit_lineage_list, use.names = FALSE) == l
      ]
      do.call(rbind, raw_matrix_list[fit_ids])
    }),
    lineages
  )
  map_lineage_matrices <- function(matrix_list, transform) {
    transformed_list <- matrix_list
    for (l in lineages) {
      fit_ids <- names(fit_lineage_list)[
        unlist(fit_lineage_list, use.names = FALSE) == l
      ]
      if (length(fit_ids) == 0) {
        next
      }
      row_counts <- vapply(matrix_list[fit_ids], nrow, integer(1))
      combined <- do.call(rbind, matrix_list[fit_ids])
      transformed <- as_matrix(transform(combined))
      if (!identical(dim(transformed), dim(combined))) {
        log_message(
          "{.arg exp_method} must preserve the dimensions of the expression matrix.",
          message_type = "error"
        )
      }
      dimnames(transformed) <- dimnames(combined)
      row_ends <- cumsum(row_counts)
      row_starts <- row_ends - row_counts + 1L
      for (i in seq_along(fit_ids)) {
        rows <- if (row_counts[[i]] == 0) {
          integer()
        } else {
          seq.int(row_starts[[i]], row_ends[[i]])
        }
        transformed_list[[fit_ids[[i]]]] <- transformed[
          rows, ,
          drop = FALSE
        ]
      }
    }
    transformed_list
  }

  if (is.function(exp_method)) {
    apply_exp_method <- function(x) t(exp_method(t(x)))
    raw_matrix_list <- map_lineage_matrices(
      raw_matrix_list,
      apply_exp_method
    )
    fitted_matrix_list <- map_lineage_matrices(
      fitted_matrix_list,
      apply_exp_method
    )
    upr_matrix_list <- map_lineage_matrices(upr_matrix_list, apply_exp_method)
    lwr_matrix_list <- map_lineage_matrices(lwr_matrix_list, apply_exp_method)
  } else {
    for (fit_id in names(raw_matrix_list)) {
      l <- fit_lineage_list[[fit_id]]
      transform_reference <- transform_reference_list[[l]]
      raw_matrix <- raw_matrix_list[[fit_id]]
      fitted_matrix <- fitted_matrix_list[[fit_id]]
      upr_matrix <- upr_matrix_list[[fit_id]]
      lwr_matrix <- lwr_matrix_list[[fit_id]]

      if (exp_method == "zscore") {
        center <- Matrix::colMeans(transform_reference)
        scale_sd <- MatrixGenerics::colSds(transform_reference)
        raw_matrix <- scale(raw_matrix, center = center, scale = scale_sd)
        fitted_matrix <- scale(fitted_matrix, center = center, scale = scale_sd)
        upr_matrix <- scale(upr_matrix, center = center, scale = scale_sd)
        lwr_matrix <- scale(lwr_matrix, center = center, scale = scale_sd)
      } else if (exp_method == "fc") {
        colm <- Matrix::colMeans(transform_reference)
        raw_matrix <- t(t(raw_matrix) / colm)
        fitted_matrix <- t(t(fitted_matrix) / colm)
        upr_matrix <- t(t(upr_matrix) / colm)
        lwr_matrix <- t(t(lwr_matrix) / colm)
      } else if (exp_method == "log2fc") {
        colm <- Matrix::colMeans(transform_reference)
        raw_matrix <- t(log2(t(raw_matrix) / colm))
        fitted_matrix <- t(log2(t(fitted_matrix) / colm))
        upr_matrix <- t(log2(t(upr_matrix) / colm))
        lwr_matrix <- t(log2(t(lwr_matrix) / colm))
      } else if (exp_method == "log1p") {
        raw_matrix <- log1p(raw_matrix)
        fitted_matrix <- log1p(fitted_matrix)
        upr_matrix <- log1p(upr_matrix)
        lwr_matrix <- log1p(lwr_matrix)
      }
      raw_matrix_list[[fit_id]] <- raw_matrix
      fitted_matrix_list[[fit_id]] <- fitted_matrix
      upr_matrix_list[[fit_id]] <- upr_matrix
      lwr_matrix_list[[fit_id]] <- lwr_matrix
    }
  }

  replace_infinite <- function(x) {
    if (is.null(fit.by)) {
      infinite <- is.infinite(x)
      if (!any(infinite)) {
        return(x)
      }
      finite_values <- x[is.finite(x)]
      cap <- if (length(finite_values) > 0) {
        max(abs(finite_values))
      } else {
        NA_real_
      }
      x[infinite] <- cap * ifelse(x[infinite] > 0, 1, -1)
      return(x)
    }
    for (feature in colnames(x)) {
      infinite <- is.infinite(x[, feature])
      if (!any(infinite)) {
        next
      }
      finite_values <- x[is.finite(x[, feature]), feature]
      cap <- if (length(finite_values) > 0) {
        max(abs(finite_values))
      } else {
        NA_real_
      }
      x[infinite, feature] <- cap * ifelse(
        x[infinite, feature] > 0,
        1,
        -1
      )
    }
    x
  }
  raw_matrix_list <- map_lineage_matrices(raw_matrix_list, replace_infinite)
  fitted_matrix_list <- map_lineage_matrices(
    fitted_matrix_list,
    replace_infinite
  )
  upr_matrix_list <- map_lineage_matrices(upr_matrix_list, replace_infinite)
  lwr_matrix_list <- map_lineage_matrices(lwr_matrix_list, replace_infinite)

  for (fit_id in names(raw_matrix_list)) {
    l <- fit_lineage_list[[fit_id]]
    fit_group <- fit_group_list[[fit_id]]
    raw_matrix <- raw_matrix_list[[fit_id]]
    fitted_matrix <- fitted_matrix_list[[fit_id]]
    upr_matrix <- upr_matrix_list[[fit_id]]
    lwr_matrix <- lwr_matrix_list[[fit_id]]

    raw <- as.data.frame(cbind(
      cell_metadata[rownames(raw_matrix), c(l, "x_assign")],
      raw_matrix
    ))
    colnames(raw)[1] <- "Pseudotime"
    raw[["Cell"]] <- rownames(raw)
    raw[["Value"]] <- "raw"
    raw <- reshape2::melt(
      raw,
      id.vars = c("Cell", "Pseudotime", "x_assign", "Value"),
      value.name = "exp",
      variable.name = "Features"
    )

    fitted <- as.data.frame(cbind(
      cell_metadata[rownames(fitted_matrix), c(l, "x_assign")],
      fitted_matrix
    ))
    colnames(fitted)[1] <- "Pseudotime"
    fitted[["Cell"]] <- rownames(fitted)
    fitted[["Value"]] <- "fitted"
    fitted <- reshape2::melt(
      fitted,
      id.vars = c("Cell", "Pseudotime", "x_assign", "Value"),
      value.name = "exp",
      variable.name = "Features"
    )

    upr <- as.data.frame(
      cbind(
        cell_metadata[rownames(upr_matrix), c(l, "x_assign")],
        upr_matrix
      )
    )
    colnames(upr)[1] <- "Pseudotime"
    upr[["Cell"]] <- rownames(upr)
    upr[["Value"]] <- "upr"
    upr <- reshape2::melt(
      upr,
      id.vars = c("Cell", "Pseudotime", "x_assign", "Value"),
      value.name = "exp",
      variable.name = "Features"
    )

    lwr <- as.data.frame(cbind(
      cell_metadata[rownames(lwr_matrix), c(l, "x_assign")],
      lwr_matrix
    ))
    colnames(lwr)[1] <- "Pseudotime"
    lwr[["Cell"]] <- rownames(lwr)
    lwr[["Value"]] <- "lwr"
    lwr <- reshape2::melt(
      lwr,
      id.vars = c("Cell", "Pseudotime", "x_assign", "Value"),
      value.name = "exp",
      variable.name = "Features"
    )

    raw[["upr"]] <- NA
    raw[["lwr"]] <- NA
    fitted[["upr"]] <- upr[["exp"]]
    fitted[["lwr"]] <- lwr[["exp"]]

    df_tmp <- rbind(raw, fitted)
    df_tmp[["Lineages"]] <- factor(l, levels = lineages)
    df_tmp[[fit_group_col]] <- if (is.null(fit.by)) {
      factor("all")
    } else {
      factor(fit_group, levels = fit_group_levels)
    }
    df_list[[fit_id]] <- df_tmp
  }
  df_all <- do.call(rbind, df_list)
  rownames(df_all) <- NULL

  if (!is.null(group.by)) {
    cell_group <- srt@meta.data[df_all[["Cell"]], group.by, drop = FALSE]
    if (!is.factor(cell_group[, group.by])) {
      cell_group[, group.by] <- factor(
        cell_group[, group.by],
        levels = unique(cell_group[, group.by])
      )
    }
    df_all <- cbind(df_all, cell_group)
  }
  lineage_code <- match(as.character(df_all[["Lineages"]]), lineages)
  feature_code <- match(as.character(df_all[["Features"]]), features)
  fit_group_code <- match(
    as.character(df_all[[fit_group_col]]),
    fit_group_levels
  )
  if (isTRUE(compare_lineages) && isTRUE(compare_features)) {
    fit_series_index <- lineage_code +
      (feature_code - 1L) * length(lineages)
    fit_series_labels <- as.vector(outer(
      lineages,
      features,
      paste,
      sep = " - "
    ))
  } else if (isTRUE(compare_features)) {
    fit_series_index <- feature_code
    fit_series_labels <- features
  } else if (isTRUE(compare_lineages)) {
    fit_series_index <- lineage_code
    fit_series_labels <- lineages
  } else {
    fit_series_index <- rep(1L, nrow(df_all))
    fit_series_labels <- "Series"
  }
  fit_series_levels <- paste0("series", seq_along(fit_series_labels))
  names(fit_series_labels) <- fit_series_levels
  df_all[[fit_series_col]] <- factor(
    paste0("series", fit_series_index),
    levels = fit_series_levels
  )
  fit_path_index <- lineage_code +
    (feature_code - 1L) * length(lineages) +
    (fit_group_code - 1L) * length(lineages) * length(features)
  df_all[[fit_path_col]] <- factor(paste0("path", fit_path_index))

  fitted_group_levels <- if (is.null(fit.by)) {
    character()
  } else {
    fitted_groups_all <- unique(as.character(df_all[[fit_group_col]][
      df_all[["Value"]] == "fitted" & is.finite(df_all[["exp"]])
    ]))
    fit_group_levels[fit_group_levels %in% fitted_groups_all]
  }
  fitted_series_all <- unique(as.character(df_all[[fit_series_col]][
    df_all[["Value"]] == "fitted" & is.finite(df_all[["exp"]])
  ]))
  fit_series <- fit_series_levels[fit_series_levels %in% fitted_series_all]
  linetype_values <- c(
    "solid", "dashed", "dotted", "dotdash", "longdash", "twodash"
  )[seq_along(fit_series_levels)]
  names(linetype_values) <- fit_series_levels

  if (!is.null(cells)) {
    df_all <- df_all[df_all[["Cell"]] %in% cells, , drop = FALSE]
  }
  df_all <- df_all[sample(seq_len(nrow(df_all))), , drop = FALSE]

  plist <- list()
  legend <- NULL
  if (isTRUE(compare_lineages)) {
    lineages_use <- list(lineages)
    lineages_formula <- "."
  } else {
    lineages_use <- lineages
    lineages_formula <- "Lineages"
  }
  if (isTRUE(compare_features)) {
    features_use <- list(features)
    features_formula <- "."
  } else {
    features_use <- features
    features_formula <- "Features"
  }
  formula <- paste(lineages_formula, "~", features_formula)
  fill_by <- "Lineages"
  if (lineages_formula == "." && length(lineages) > 1) {
    lineages_guide <- TRUE
  } else {
    lineages_guide <- FALSE
    if (isTRUE(compare_features)) {
      fill_by <- "Features"
    }
  }
  if (features_formula == "." && length(features) > 1) {
    features_guide <- TRUE
  } else {
    features_guide <- FALSE
  }
  if (!is.null(fit.by)) {
    fill_by <- fit_group_col
  }

  for (l in lineages_use) {
    for (f in features_use) {
      df <- subset(
        df_all,
        df_all[["Lineages"]] %in% l & df_all[["Features"]] %in% f
      )
      if (x_order == "rank") {
        df[, "x_assign"] <- rank(df[, "x_assign"])
        df[, "Pseudotime"] <- rank(df[, "Pseudotime"])
      }
      point_columns <- c("Cell", "x_assign", "exp", group.by)
      if (!isTRUE(compare_features)) {
        point_columns <- c(point_columns, "Features")
      }
      if (!isTRUE(compare_lineages)) {
        point_columns <- c(point_columns, "Lineages")
      }
      df_point <- unique(df[
        df[["Value"]] == "raw",
        unique(point_columns)
      ])
      fitted_groups <- unique(as.character(df[[fit_group_col]][
        df[["Value"]] == "fitted" & is.finite(df[["exp"]])
      ]))
      point_groups <- if (is.null(fit.by)) {
        character()
      } else {
        unique(as.character(df[[fit_group_col]][df[["Value"]] == "raw"]))
      }
      hide_point_guide <- identical(group.by, fit.by) &&
        isTRUE(add_line) &&
        all(point_groups %in% fitted_groups)
      distinguish_fit_series <- !is.null(fit.by) && length(fit_series) > 1
      if (
        (isTRUE(add_line) || isTRUE(add_interval)) &&
          isTRUE(distinguish_fit_series) &&
          length(fit_series_levels) > 6
      ) {
        log_message(
          "Grouped combined plots support at most six lineage-feature series. Set {.arg compare_lineages} or {.arg compare_features} to {.val FALSE}.",
          message_type = "error"
        )
      }
      if (isTRUE(compare_features)) {
        raw_point <- NULL
      } else {
        if (isTRUE(add_point)) {
          if (is.null(group.by)) {
            raw_point <- geom_point(
              data = df_point,
              mapping = aes(x = .data[["x_assign"]], y = .data[["exp"]]),
              size = pt.size,
              alpha = 0.8
            )
          } else {
            raw_point <- list(
              geom_point(
                data = df_point,
                mapping = aes(
                  x = .data[["x_assign"]],
                  y = .data[["exp"]],
                  color = .data[[group.by]]
                ),
                size = pt.size,
                alpha = 0.8
              ),
              scale_color_manual(
                values = palette_colors(
                  df[[group.by]],
                  palette = point_palette,
                  palcolor = point_palcolor
                ),
                guide = if (isTRUE(hide_point_guide)) {
                  "none"
                } else {
                  guide_legend(
                    override.aes = list(alpha = 1, size = 3),
                    order = 1
                  )
                }
              ),
              ggnewscale::new_scale_color()
            )
          }
        } else {
          raw_point <- NULL
        }
      }
      if (isTRUE(add_rug)) {
        if (is.null(group.by)) {
          rug <- list(geom_rug(
            data = df_point,
            mapping = aes(x = .data[["x_assign"]]),
            alpha = 1,
            length = grid::unit(0.05, "npc"),
            show.legend = FALSE
          ))
        } else {
          rug <- list(
            geom_rug(
              data = df_point,
              mapping = aes(
                x = .data[["x_assign"]],
                color = .data[[group.by]]
              ),
              alpha = 1,
              length = grid::unit(0.05, "npc"),
              show.legend = isTRUE(compare_features)
            ),
            scale_color_manual(
              values = palette_colors(
                df[[group.by]],
                palette = point_palette,
                palcolor = point_palcolor
              )
            ),
            ggnewscale::new_scale_color()
          )
        }
      } else {
        rug <- NULL
      }

      if (isTRUE(add_interval)) {
        interval_data <- subset(
          df,
          df[["Value"]] == "fitted" &
            is.finite(df[["exp"]]) &
            is.finite(df[["lwr"]]) &
            is.finite(df[["upr"]])
        )
        distinguish_interval_series <- !is.null(fit.by) &&
          !isTRUE(add_line) &&
          isTRUE(distinguish_fit_series)
        interval_mapping <- if (isTRUE(distinguish_interval_series)) {
          aes(
            x = .data[["Pseudotime"]],
            y = .data[["exp"]],
            ymin = .data[["lwr"]],
            ymax = .data[["upr"]],
            fill = .data[[fill_by]],
            linetype = .data[[fit_series_col]],
            group = .data[[fit_path_col]]
          )
        } else {
          aes(
            x = .data[["Pseudotime"]],
            y = .data[["exp"]],
            ymin = .data[["lwr"]],
            ymax = .data[["upr"]],
            fill = .data[[fill_by]],
            group = .data[[fit_path_col]]
          )
        }
        interval_geom <- if (isTRUE(distinguish_interval_series)) {
          geom_ribbon(
            data = interval_data,
            mapping = interval_mapping,
            alpha = 0.4,
            color = "grey30",
            linewidth = 0.4
          )
        } else {
          geom_ribbon(
            data = interval_data,
            mapping = interval_mapping,
            alpha = 0.4,
            color = "grey90"
          )
        }
        interval_linetype_scale <- if (isTRUE(distinguish_interval_series)) {
          scale_linetype_manual(
            name = "Lineage - feature",
            values = linetype_values,
            limits = fit_series,
            labels = fit_series_labels[fit_series],
            drop = FALSE,
            guide = guide_legend(order = 3)
          )
        } else {
          NULL
        }
        interval <- list(
          interval_geom,
          scale_fill_manual(
            values = palette_colors(
              df[[fill_by]],
              palette = if (!is.null(fit.by) && identical(group.by, fit.by)) {
                point_palette
              } else {
                line_palette
              },
              palcolor = if (!is.null(fit.by) && identical(group.by, fit.by)) {
                point_palcolor
              } else {
                line_palcolor
              }
            ),
            name = if (is.null(fit.by)) NULL else fit.by,
            limits = if (is.null(fit.by)) NULL else fitted_group_levels,
            drop = is.null(fit.by),
            guide = if (!is.null(fit.by) && isTRUE(add_line)) {
              "none"
            } else if (!is.null(fit.by)) {
              guide_legend(order = 2)
            } else if (
              fill_by == "Features" || lineages_guide || length(l) == 1
            ) {
              "none"
            } else {
              guide_legend()
            }
          ),
          interval_linetype_scale,
          ggnewscale::new_scale_fill()
        )
      } else {
        interval <- NULL
      }
      if (!is.null(fit.by)) {
        if (isTRUE(add_line)) {
          line_mapping <- if (isTRUE(distinguish_fit_series)) {
            aes(
              x = .data[["Pseudotime"]],
              y = .data[["exp"]],
              color = .data[[fit_group_col]],
              linetype = .data[[fit_series_col]],
              group = .data[[fit_path_col]]
            )
          } else {
            aes(
              x = .data[["Pseudotime"]],
              y = .data[["exp"]],
              color = .data[[fit_group_col]],
              group = .data[[fit_path_col]]
            )
          }
          linetype_scale <- if (isTRUE(distinguish_fit_series)) {
            scale_linetype_manual(
              name = "Lineage - feature",
              values = linetype_values,
              limits = fit_series,
              labels = fit_series_labels[fit_series],
              guide = guide_legend(order = 3)
            )
          } else {
            NULL
          }
          line <- list(
            geom_line(
              data = subset(
                df,
                df[["Value"]] == "fitted" & is.finite(df[["exp"]])
              ),
              mapping = line_mapping,
              linewidth = line.size,
              alpha = 0.8
            ),
            scale_color_manual(
              name = fit.by,
              values = palette_colors(
                df[[fit_group_col]],
                palette = if (identical(group.by, fit.by)) {
                  point_palette
                } else {
                  line_palette
                },
                palcolor = if (identical(group.by, fit.by)) {
                  point_palcolor
                } else {
                  line_palcolor
                }
              ),
              limits = fitted_group_levels,
              drop = FALSE,
              guide = guide_legend(
                override.aes = list(alpha = 1, linewidth = 2),
                order = 2
              )
            ),
            linetype_scale,
            ggnewscale::new_scale_color()
          )
        } else {
          line <- NULL
        }
      } else if (isTRUE(compare_features)) {
        line <- list(
          geom_line(
            data = subset(df, df[["Value"]] == "fitted"),
            mapping = aes(
              x = .data[["Pseudotime"]],
              y = .data[["exp"]],
              color = .data[["Features"]],
              group = .data[[fit_path_col]]
            ),
            linewidth = line.size,
            alpha = 0.8
          ),
          scale_color_manual(
            values = palette_colors(
              df[["Features"]],
              palette = line_palette,
              palcolor = line_palcolor
            ),
            guide = if (features_guide) {
              guide_legend(
                override.aes = list(alpha = 1, size = 2),
                order = 2
              )
            } else {
              "none"
            }
          ),
          ggnewscale::new_scale_color()
        )
      } else {
        if (isTRUE(add_line)) {
          line <- list(
            geom_line(
              data = subset(df, df[["Value"]] == "fitted"),
              mapping = aes(
                x = .data[["Pseudotime"]],
                y = .data[["exp"]],
                color = .data[["Lineages"]],
                group = .data[[fit_path_col]]
              ),
              linewidth = line.size,
              alpha = 0.8
            ),
            scale_color_manual(
              values = palette_colors(
                df[["Lineages"]],
                palette = line_palette,
                palcolor = line_palcolor
              ),
              guide = if (lineages_guide) {
                guide_legend(
                  override.aes = list(alpha = 1, size = 2),
                  order = 2
                )
              } else {
                "none"
              }
            ),
            ggnewscale::new_scale_color()
          )
        } else {
          line <- NULL
        }
      }

      x_trans <- ifelse(flip, "reverse", "identity")
      x_trans <- ifelse(
        reverse,
        setdiff(c("reverse", "identity"), x_trans),
        x_trans
      )
      p <- ggplot() +
        scale_x_continuous(trans = x_trans, expand = expansion(c(0, 0))) +
        scale_y_continuous(expand = expansion(c(0.1, 0.05))) +
        raw_point +
        rug +
        interval +
        line +
        labs(
          x = ifelse(x_order == "rank", "Pseudotime(rank)", "Pseudotime"),
          y = exp_name
        ) +
        facet_grid(
          stats::formula(formula),
          scales = "free"
        ) +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction
        )

      if (isTRUE(flip)) {
        p <- p + coord_flip()
      }
      if (is.null(legend)) {
        legend <- get_legend(
          p +
            theme(legend.position = "bottom")
        )
      }
      plist[[paste(
        paste0(l, collapse = "_"),
        paste0(f, collapse = "_"),
        sep = "."
      )]] <- p + theme(legend.position = "none")
    }
  }

  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- patchwork::wrap_plots(
        plotlist = plist,
        nrow = nrow,
        ncol = ncol,
        byrow = byrow
      )
    } else {
      plot <- plist[[1]]
    }
    if (legend.position != "none") {
      gtable <- as_grob(plot)
      gtable <- add_grob(gtable, legend, legend.position)
      plot <- patchwork::wrap_plots(gtable)
    }
    return(plot)
  } else {
    return(plist)
  }
}
