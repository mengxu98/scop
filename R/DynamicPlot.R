#' @title Plot dynamic features across pseudotime
#'
#' @md
#' @param srt A Seurat object.
#' @param features A character vector specifying the features to plot.
#' @param lineages A character vector specifying the lineages to plot.
#' @param group.by A character specifying a metadata column to group the cells by. Default is NULL.
#' @param cells A character vector specifying the cells to include in the plot. Default is NULL.
#' @param layer A character string specifying the layer to use for the analysis. Default is "counts".
#' @param assay A character string specifying the assay to use for the analysis. Default is NULL.
#' @param family A character specifying the model used to calculate the dynamic features if needed.
#' By default, this parameter is set to NULL, and the appropriate family will be automatically determined.
#' @param exp_method A character specifying the method to transform the expression values.
#' Default is "log1p" with options "log1p", "raw", "zscore", "fc", "log2fc".
#' @param lib_normalize A boolean specifying whether to normalize the expression values using library size.
#' Default the `layer` is counts, this parameter is set to TRUE.
#' Otherwise, it is set to FALSE.
#' @param libsize A numeric vector specifying the library size for each cell. Default is NULL.
#' @param compare_lineages A boolean specifying whether to compare the lineages in the plot. Default is TRUE.
#' @param compare_features A boolean specifying whether to compare the features in the plot. Default is FALSE.
#' @param add_line A boolean specifying whether to add lines to the plot. Default is TRUE.
#' @param add_interval A boolean specifying whether to add confidence intervals to the plot. Default is TRUE.
#' @param line.size A numeric specifying the size of the lines. Default is 1.
#' @param line_palette A character string specifying the name of the palette to use for the line colors. Default is "Dark2".
#' @param line_palcolor A vector specifying the colors to use for the line palette. Default is NULL.
#' @param add_point A boolean specifying whether to add points to the plot. Default is TRUE.
#' @param pt.size A numeric specifying the size of the points. Default is 1.
#' @param point_palette A character string specifying the name of the palette to use for the point colors. Default is "Paired".
#' @param point_palcolor A vector specifying the colors to use for the point palette. Default is NULL.
#' @param add_rug A boolean specifying whether to add rugs to the plot. Default is TRUE.
#' @param flip A boolean specifying whether to flip the x-axis. Default is FALSE.
#' @param reverse A boolean specifying whether to reverse the x-axis. Default is FALSE.
#' @param x_order A character specifying the order of the x-axis values. Default is c("value", "rank").
#' @param aspect.ratio A numeric specifying the aspect ratio of the plot. Default is NULL.
#' @param legend.position A character string specifying the position of the legend in the plot. Default is "right".
#' @param legend.direction A character string specifying the direction of the legend in the plot. Default is "vertical".
#' @param theme_use A character string specifying the name of the theme to use for the plot. Default is "theme_scop".
#' @param theme_args A list specifying the arguments to pass to the theme function. Default is list().
#' @param combine A boolean specifying whether to combine multiple plots into a single plot. Default is TRUE.
#' @param nrow A numeric specifying the number of rows in the combined plot. Default is NULL.
#' @param ncol A numeric specifying the number of columns in the combined plot. Default is NULL.
#' @param byrow A boolean specifying whether to fill plots by row in the combined plot. Default is TRUE.
#' @param seed A numeric specifying the random seed. Default is 11.
#'
#' @seealso [RunDynamicFeatures]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
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
#'   compare_features = TRUE
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
    point_palette = "Paired",
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
    seed = 11) {
  set.seed(seed)

  check_r("MatrixGenerics")
  x_order <- match.arg(x_order)
  if (!is.null(group.by) && !group.by %in% colnames(srt@meta.data)) {
    log_message(
      paste0(group.by, " is not in the meta.data of srt object."),
      message_type = "error"
    )
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
  for (l in lineages) {
    features_exist <- c()
    raw_matrix <- NULL
    if (paste0("DynamicFeatures_", l) %in% names(srt@tools)) {
      raw_matrix <- srt@tools[[paste0("DynamicFeatures_", l)]][["raw_matrix"]][
        ,
        -1
      ]
      fitted_matrix <- srt@tools[[paste0("DynamicFeatures_", l)]][[
        "fitted_matrix"
      ]][, -1]
      upr_matrix <- srt@tools[[paste0("DynamicFeatures_", l)]][["upr_matrix"]][
        ,
        -1
      ]
      lwr_matrix <- srt@tools[[paste0("DynamicFeatures_", l)]][["lwr_matrix"]][
        ,
        -1
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
        libsize = libsize
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
  for (l in lineages) {
    raw_matrix <- raw_matrix_list[[l]]
    fitted_matrix <- fitted_matrix_list[[l]]
    upr_matrix <- upr_matrix_list[[l]]
    lwr_matrix <- lwr_matrix_list[[l]]
    if (isTRUE(lib_normalize) && min(raw_matrix[, gene], na.rm = TRUE) >= 0) {
      if (!is.null(libsize)) {
        libsize_use <- libsize
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
        stats::median(y_libsize)
    }

    if (is.function(exp_method)) {
      raw_matrix <- t(exp_method(t(raw_matrix)))
      fitted_matrix <- t(exp_method(t(fitted_matrix)))
      upr_matrix <- t(exp_method(t(upr_matrix)))
      lwr_matrix <- t(exp_method(t(lwr_matrix)))
    } else if (exp_method == "raw") {
      raw_matrix <- raw_matrix
      fitted_matrix <- fitted_matrix
      upr_matrix <- upr_matrix
      lwr_matrix <- lwr_matrix
    } else if (exp_method == "zscore") {
      center <- colMeans(raw_matrix)
      sd <- MatrixGenerics::colSds(raw_matrix)
      raw_matrix <- scale(raw_matrix, center = center, scale = sd)
      fitted_matrix <- scale(fitted_matrix, center = center, scale = sd)
      upr_matrix <- scale(upr_matrix, center = center, scale = sd)
      lwr_matrix <- scale(lwr_matrix, center = center, scale = sd)
    } else if (exp_method == "fc") {
      colm <- colMeans(raw_matrix)
      raw_matrix <- t(t(raw_matrix) / colm)
      fitted_matrix <- t(t(fitted_matrix) / colm)
      upr_matrix <- t(t(upr_matrix) / colm)
      lwr_matrix <- t(t(lwr_matrix) / colm)
    } else if (exp_method == "log2fc") {
      colm <- colMeans(raw_matrix)
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
    raw_matrix[is.infinite(raw_matrix)] <- max(
      abs(raw_matrix[!is.infinite(raw_matrix)]),
      na.rm = TRUE
    ) *
      ifelse(raw_matrix[is.infinite(raw_matrix)] > 0, 1, -1)
    fitted_matrix[is.infinite(fitted_matrix)] <- max(abs(fitted_matrix[
      !is.infinite(fitted_matrix)
    ])) *
      ifelse(fitted_matrix[is.infinite(fitted_matrix)] > 0, 1, -1)
    upr_matrix[is.infinite(upr_matrix)] <- max(
      abs(upr_matrix[!is.infinite(upr_matrix)]),
      na.rm = TRUE
    ) *
      ifelse(upr_matrix[is.infinite(upr_matrix)] > 0, 1, -1)
    lwr_matrix[is.infinite(lwr_matrix)] <- max(
      abs(lwr_matrix[!is.infinite(lwr_matrix)]),
      na.rm = TRUE
    ) *
      ifelse(lwr_matrix[is.infinite(lwr_matrix)] > 0, 1, -1)

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
    df_list[[l]] <- df_tmp
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
  df_all[["LineagesFeatures"]] <- paste(
    df_all[["Lineages"]],
    df_all[["Features"]],
    sep = "-"
  )

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
      df_point <- unique(df[
        df[["Value"]] == "raw",
        c("Cell", "x_assign", "exp", group.by)
      ])
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
                )
              ),
              scale_fill_manual(
                values = palette_colors(
                  df[[group.by]],
                  palette = point_palette,
                  palcolor = point_palcolor
                ),
                guide = guide_legend(
                  override.aes = list(alpha = 1, size = 3),
                  order = 1
                )
              ),
              ggnewscale::new_scale_color(),
              ggnewscale::new_scale_fill()
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
        interval <- list(
          geom_ribbon(
            data = subset(df, df[["Value"]] == "fitted"),
            mapping = aes(
              x = .data[["Pseudotime"]],
              y = .data[["exp"]],
              ymin = .data[["lwr"]],
              ymax = .data[["upr"]],
              fill = .data[[fill_by]],
              group = .data[["LineagesFeatures"]]
            ),
            alpha = 0.4,
            color = "grey90"
          ),
          scale_fill_manual(
            values = palette_colors(
              df[[fill_by]],
              palette = line_palette,
              palcolor = line_palcolor
            ),
            guide = if (
              fill_by == "Features" || lineages_guide || length(l) == 1
            ) {
              "none"
            } else {
              guide_legend()
            }
          ),
          ggnewscale::new_scale_fill()
        )
      } else {
        interval <- NULL
      }
      if (isTRUE(compare_features)) {
        line <- list(
          geom_line(
            data = subset(df, df[["Value"]] == "fitted"),
            mapping = aes(
              x = .data[["Pseudotime"]],
              y = .data[["exp"]],
              color = .data[["Features"]],
              group = .data[["LineagesFeatures"]]
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
                group = .data[["LineagesFeatures"]]
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
