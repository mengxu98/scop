#' Features correlation plot
#' This function creates a correlation plot to visualize the pairwise correlations between selected features in a Seurat object.
#'
#' @param srt A Seurat object.
#' @param features A character vector specifying the features to compare. Should be present in both the assay data and the metadata of the Seurat object.
#' @param group.by A character string specifying the column in the metadata to group cells by.
#' @param split.by A character string specifying the column in the metadata to split the plot by.
#' @param cells A character vector specifying the cells to include in the plot. If NULL (default), all cells will be included.
#' @param layer A character string specifying the layer in the Seurat object to use. Defaults to "data".
#' @param assay A character string specifying the assay to use. Defaults to the default assay in the Seurat object.
#' @param cor_method A character string specifying the correlation method to use. Can be "pearson" (default) or "spearman".
#' @param adjust A numeric value specifying the adjustment factor for the width of the violin plots. Defaults to 1.
#' @param margin A numeric value specifying the margin size for the plot. Defaults to 1.
#' @param reverse A logical value indicating whether to reverse the order of the features in the plot. Defaults to FALSE.
#' @param add_equation A logical value indicating whether to add the equation of the linear regression line to each scatter plot. Defaults to FALSE.
#' @param add_r2 A logical value indicating whether to add the R-squared value of the linear regression line to each scatter plot. Defaults to TRUE.
#' @param add_pvalue A logical value indicating whether to add the p-value of the linear regression line to each scatter plot. Defaults to TRUE.
#' @param add_smooth A logical value indicating whether to add a smoothed line to each scatter plot. Defaults to TRUE.
#' @param palette A character string specifying the name of the color palette to use for the groups. Defaults to "Paired".
#' @param palcolor A character string specifying the color for the groups. Defaults to NULL.
#' @param cor_palette A character string specifying the name of the color palette to use for the correlation. Defaults to "RuBu".
#' @param cor_palcolor A character string specifying the color for the correlation. Defaults to "RuBu".
#' @param cor_range A two-length numeric vector specifying the range for the correlation.
#' @param pt.size A numeric value specifying the size of the points in the scatter plots. If NULL (default), the size will be automatically determined based on the number of cells.
#' @param pt.alpha A numeric value between 0 and 1 specifying the transparency of the points in the scatter plots. Defaults to 1.
#' @param cells.highlight A logical value or a character vector specifying the cells to highlight in the scatter plots. If TRUE, all cells will be highlighted. Defaults to NULL.
#' @param cols.highlight A character string specifying the color for the highlighted cells. Defaults to "black".
#' @param sizes.highlight A numeric value specifying the size of the highlighted cells in the scatter plots. Defaults to 1.
#' @param alpha.highlight A numeric value between 0 and 1 specifying the transparency of the highlighted cells in the scatter plots. Defaults to 1.
#' @param stroke.highlight A numeric value specifying the stroke size of the highlighted cells in the scatter plots. Defaults to 0.5.
#' @param calculate_coexp A logical value indicating whether to calculate the co-expression of selected features. Defaults to FALSE.
#' @param raster A logical value indicating whether to use raster graphics for scatter plots. Defaults to NULL.
#' @param raster.dpi A numeric vector specifying the dpi (dots per inch) resolution for raster graphics in the scatter plots. Defaults to c(512, 512).
#' @param aspect.ratio A numeric value specifying the aspect ratio of the scatter plots. Defaults to 1.
#' @param title A character string specifying the title for the correlation plot. Defaults to NULL.
#' @param subtitle A character string specifying the subtitle for the correlation plot. Defaults to NULL.
#' @param legend.position A character string specifying the position of the legend. Can be "right" (default), "left", "top", or "bottom".
#' @param legend.direction A character string specifying the direction of the legend. Can be "vertical" (default) or "horizontal".
#' @param theme_use A character string specifying the name of the theme to use for the plot. Defaults to "theme_scop".
#' @param theme_args A list of arguments to pass to the theme function. Defaults to an empty list.
#' @param combine A logical value indicating whether to combine the plots into a single plot. Defaults to TRUE.
#' @param nrow A numeric value specifying the number of rows in the combined plot. If NULL (default), the number of rows will be automatically determined.
#' @param ncol A numeric value specifying the number of columns in the combined plot. If NULL (default), the number of columns will be automatically determined.
#' @param byrow A logical value indicating whether to fill the combined plot byrow (top to bottom, left to right). Defaults to TRUE.
#' @param force A logical value indicating whether to force the creation of the plot, even if it contains more than 50 subplots. Defaults to FALSE.
#' @param seed A numeric value specifying the random seed for reproducibility. Defaults to 11.
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- Seurat::NormalizeData(pancreas_sub)
#' FeatureCorPlot(
#'   pancreas_sub,
#'   features = rownames(pancreas_sub)[1:5],
#'   group.by = "SubCellType"
#' )
#' FeatureCorPlot(
#'   pancreas_sub,
#'   features = c(
#'     "nFeature_RNA",
#'     "nCount_RNA",
#'     "nFeature_spliced",
#'     "nCount_spliced",
#'     "nFeature_unspliced",
#'     "nCount_unspliced"
#'   ),
#'   group.by = "SubCellType",
#'   cor_palette = "Greys",
#'   cor_range = c(0, 1)
#' )
#' FeatureCorPlot(
#'   pancreas_sub,
#'   features = c("nFeature_RNA", "nCount_RNA"),
#'   group.by = "SubCellType",
#'   add_equation = TRUE
#' )
FeatureCorPlot <- function(
    srt,
    features,
    group.by = NULL,
    split.by = NULL,
    cells = NULL,
    layer = "data",
    assay = NULL,
    cor_method = "pearson",
    adjust = 1,
    margin = 1,
    reverse = FALSE,
    add_equation = FALSE,
    add_r2 = TRUE,
    add_pvalue = TRUE,
    add_smooth = TRUE,
    palette = "Paired",
    palcolor = NULL,
    cor_palette = "RdBu",
    cor_palcolor = NULL,
    cor_range = c(-1, 1),
    pt.size = NULL,
    pt.alpha = 1,
    cells.highlight = NULL,
    cols.highlight = "black",
    sizes.highlight = 1,
    alpha.highlight = 1,
    stroke.highlight = 0.5,
    calculate_coexp = FALSE,
    raster = NULL,
    raster.dpi = c(512, 512),
    aspect.ratio = 1,
    title = NULL,
    subtitle = NULL,
    legend.position = "right",
    legend.direction = "vertical",
    theme_use = "theme_scop",
    theme_args = list(),
    combine = TRUE,
    nrow = NULL,
    ncol = NULL,
    byrow = TRUE,
    force = FALSE,
    seed = 11) {
  set.seed(seed)

  if (is.null(features)) {
    log_message(
      "'features' must be provided.",
      message_type = "error"
    )
  }
  if (!inherits(features, "character")) {
    log_message(
      "'features' is not a character vectors",
      message_type = "error"
    )
  }
  assay <- assay %||% DefaultAssay(srt)
  if (is.null(split.by)) {
    split.by <- "All.groups"
    srt@meta.data[[split.by]] <- factor("")
  }
  if (is.null(group.by)) {
    group.by <- "All.groups"
    srt@meta.data[[group.by]] <- factor("")
  }
  for (i in c(split.by, group.by)) {
    if (!i %in% colnames(srt@meta.data)) {
      log_message(
        paste0(i, " is not in the meta.data of srt object."),
        message_type = "error"
      )
    }
    if (!is.factor(srt@meta.data[[i]])) {
      srt@meta.data[[i]] <- factor(
        srt@meta.data[[i]],
        levels = unique(srt@meta.data[[i]])
      )
    }
  }
  if (!is.null(cells.highlight) & isFALSE(cells.highlight)) {
    if (!any(cells.highlight %in% colnames(srt@assays[[1]]))) {
      log_message("No cells in 'cells.highlight' found in srt.",
        message_type = "error"
      )
    }
    if (!all(cells.highlight %in% colnames(srt@assays[[1]]))) {
      log_message(
        "Some cells in 'cells.highlight' not found in srt.",
        message_type = "warning"
      )
    }
    cells.highlight <- intersect(cells.highlight, colnames(srt@assays[[1]]))
  }
  if (isTRUE(cells.highlight)) {
    cells.highlight <- colnames(srt@assays[[1]])
  }

  features_drop <- features[
    !features %in% c(rownames(srt@assays[[assay]]), colnames(srt@meta.data))
  ]
  if (length(features_drop) > 0) {
    log_message(
      paste0(features_drop, collapse = ","),
      " are not in the features of srt.",
      message_type = "warning"
    )
    features <- features[!features %in% features_drop]
  }

  features_gene <- features[features %in% rownames(srt@assays[[assay]])]
  features_meta <- features[features %in% colnames(srt@meta.data)]
  if (length(intersect(features_gene, features_meta)) > 0) {
    log_message(
      "Features appear in both gene names and metadata names: ",
      paste0(intersect(features_gene, features_meta), collapse = ","),
      message_type = "warning"
    )
  }

  if (isTRUE(calculate_coexp) && length(features_gene) > 0) {
    if (length(features_meta) > 0) {
      log_message(
        paste(features_meta, collapse = ","),
        "is not used when calculating co-expression",
        message_type = "warning"
      )
    }
    if (status %in% c("raw_counts", "raw_normalized_counts")) {
      srt@meta.data[["CoExp"]] <- apply(
        GetAssayData5(
          srt,
          assay = assay,
          layer = layer
        )[features_gene, , drop = FALSE],
        2,
        function(x) exp(mean(log(x)))
      )
    } else if (status == "log_normalized_counts") {
      srt@meta.data[["CoExp"]] <- apply(
        expm1(
          GetAssayData5(
            srt,
            assay = assay,
            layer = layer
          )[features_gene, , drop = FALSE]
        ),
        2,
        function(x) log1p(exp(mean(log(x))))
      )
    } else {
      log_message(
        "Can not determine the data type.",
        message_type = "error"
      )
    }
    features <- c(features, "CoExp")
    features_meta <- c(features_meta, "CoExp")
  }
  if (length(features_gene) > 0) {
    dat_gene <- Matrix::t(
      GetAssayData5(
        srt,
        assay = assay,
        layer = layer
      )[features_gene, , drop = FALSE]
    )
  } else {
    dat_gene <- matrix(nrow = ncol(srt@assays[[1]]), ncol = 0)
  }
  if (length(features_meta) > 0) {
    dat_meta <- as_matrix(srt@meta.data[, features_meta, drop = FALSE])
  } else {
    dat_meta <- matrix(nrow = ncol(srt@assays[[1]]), ncol = 0)
  }
  dat_exp <- cbind(dat_gene, dat_meta)
  features <- unique(features[features %in% c(features_gene, features_meta)])
  if (length(features) < 2) {
    log_message(
      "features must be a vector of length at least 2.",
      message_type = "error"
    )
  }

  if (!is.numeric(dat_exp) && !inherits(dat_exp, "Matrix")) {
    log_message(
      "'features' must be type of numeric variable.",
      message_type = "error"
    )
  }
  if (!inherits(dat_exp, "dgCMatrix")) {
    dat_exp <- SeuratObject::as.sparse(
      as_matrix(dat_exp)
    )
  }
  if (length(features) > 10 && isFALSE(force)) {
    log_message(
      "More than 10 features to be paired compared which will generate more than 50 plots.",
      message_type = "warning"
    )
    answer <- utils::askYesNo("Are you sure to continue?", default = FALSE)
    if (isFALSE(answer)) {
      return(invisible(NULL))
    }
  }
  dat_use <- srt@meta.data[, unique(c(split.by, group.by)), drop = FALSE]
  dat_use <- cbind(dat_use, dat_exp[row.names(dat_use), , drop = FALSE])
  if (!is.null(cells)) {
    dat_use <- dat_use[intersect(rownames(dat_use), cells), , drop = FALSE]
  }

  if (is.null(pt.size)) {
    pt.size <- min(3000 / nrow(dat_use), 0.5)
  }
  raster <- raster %||% (nrow(dat_use) * ncol(combn(features, m = 2)) > 1e5)
  if (isTRUE(raster)) {
    check_r("scattermore")
  }
  if (!is.null(x = raster.dpi)) {
    if (!is.numeric(x = raster.dpi) || length(x = raster.dpi) != 2) {
      log_message(
        "'raster.dpi' must be a two-length numeric vector",
        message_type = "error"
      )
    }
  }

  plist <- list()
  colors <- palette_colors(
    levels(dat_use[[group.by]]),
    palette = palette,
    palcolor = palcolor
  )
  cor_colors <- palette_colors(
    x = seq(cor_range[1], cor_range[2], length.out = 200),
    palette = cor_palette,
    palcolor = cor_palcolor
  )
  bound <- strsplit(gsub("\\(|\\)|\\[|\\]", "", names(cor_colors)), ",")
  bound <- lapply(bound, as.numeric)
  df_bound <- do.call(rbind, bound)
  rownames(df_bound) <- cor_colors
  df_bound[1, 1] <- df_bound[1, 1] - 0.01

  pair <- as.data.frame(Matrix::t(combn(features, m = 2)))
  colnames(pair) <- c("feature1", "feature2")
  pair_expand <- expand.grid(features, features, stringsAsFactors = TRUE)
  colnames(pair_expand) <- c("feature1", "feature2")
  pair_expand[["feature1"]] <- factor(
    pair_expand[["feature1"]],
    levels = levels(pair_expand[["feature2"]])
  )

  for (s in levels(dat_use[[split.by]])) {
    dat <- dat_use[dat_use[[split.by]] == s, , drop = FALSE]
    feature_mat <- Matrix::t(dat_exp[rownames(dat), features])
    if (cor_method %in% c("pearson", "spearman")) {
      if (cor_method == "spearman") {
        feature_mat <- Matrix::t(apply(feature_mat, 1, rank))
      }
      cor_method <- "correlation"
    }
    pair_sim <- proxyC::simil(
      x = feature_mat,
      method = cor_method
    )
    if (isTRUE(reverse)) {
      order1 <- rev(pair_expand[, 1])
      order2 <- rev(pair_expand[, 2])
      levels(order1) <- rev(levels(order1))
      levels(order2) <- rev(levels(order2))
    } else {
      order1 <- pair_expand[, 1]
      order2 <- pair_expand[, 2]
    }
    plotlist <- mapply(
      FUN = function(x, y) {
        f1 <- as.character(x)
        f2 <- as.character(y)
        f1_index <- as.numeric(x)
        f2_index <- as.numeric(y)
        p <- ggplot(data = dat) +
          do.call(theme_use, theme_args) +
          theme(
            aspect.ratio = aspect.ratio,
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            plot.margin = margin(margin, margin, margin, margin),
            legend.position = "none"
          )
        if (f1_index == f2_index) {
          p <- p +
            geom_violin(
              aes(
                x = .data[[group.by]],
                y = .data[[f1]],
                fill = .data[[group.by]]
              ),
              scale = "width",
              adjust = adjust,
              trim = TRUE,
              na.rm = TRUE
            ) +
            scale_x_discrete(
              position = ifelse(isTRUE(reverse), "top", "bottom")
            ) +
            scale_y_continuous(
              position = ifelse(isTRUE(reverse), "right", "left")
            )
        } else {
          p <- p +
            scale_x_continuous(
              n.breaks = 3,
              labels = scales::number_format(),
              position = ifelse(isTRUE(reverse), "top", "bottom")
            ) +
            scale_y_continuous(
              n.breaks = 3,
              labels = scales::number_format(),
              position = ifelse(isTRUE(reverse), "right", "left")
            )
        }
        if (f1_index < f2_index) {
          if (isTRUE(raster)) {
            p <- p +
              scattermore::geom_scattermore(
                mapping = aes(
                  x = .data[[f1]],
                  y = .data[[f2]],
                  color = .data[[group.by]]
                ),
                pointsize = ceiling(pt.size),
                alpha = pt.alpha,
                pixels = raster.dpi
              )
          } else {
            p <- p +
              geom_point(
                aes(
                  x = .data[[f1]],
                  y = .data[[f2]],
                  color = .data[[group.by]]
                ),
                alpha = pt.alpha,
                size = pt.size
              )
          }
          if (isTRUE(add_smooth)) {
            p <- p +
              geom_smooth(
                aes(x = .data[[f1]], y = .data[[f2]]),
                alpha = 0.5,
                method = "lm",
                color = "red",
                formula = y ~ x,
                na.rm = TRUE
              )
          }
          if (any(isTRUE(add_equation), isTRUE(add_r2), isTRUE(add_pvalue))) {
            m <- stats::lm(dat[[f2]] ~ dat[[f1]])
            if (stats::coef(m)[2] >= 0) {
              eq1 <- substitute(
                italic(y) == a + b %.% italic(x),
                list(
                  a = format(as.numeric(stats::coef(m)[1]), digits = 2),
                  b = format(as.numeric(stats::coef(m)[2]), digits = 2)
                )
              )
            } else {
              eq1 <- substitute(
                italic(y) == a - b %.% italic(x),
                list(
                  a = format(as.numeric(stats::coef(m)[1]), digits = 2),
                  b = format(-as.numeric(stats::coef(m)[2]), digits = 2)
                )
              )
            }
            eq1 <- as.character(as.expression(eq1))
            eq2 <- substitute(
              italic(r)^2 ~ "=" ~ r2,
              list(
                r2 = format(summary(m)$r.squared, digits = 2)
              )
            )
            eq2 <- as.character(as.expression(eq2))
            eq3 <- substitute(
              italic(p) ~ "=" ~ pvalue,
              list(
                pvalue = format(summary(m)$coefficients[2, 4], digits = 2)
              )
            )
            eq3 <- as.character(as.expression(eq3))
            eqs <- c(eq1, eq2, eq3)
            vjusts <- c(1.3, 1.3 * 2, 1.3 * 2^2)
            i <- c(isTRUE(add_equation), isTRUE(add_r2), isTRUE(add_pvalue))
            p <- p +
              annotate(
                geom = GeomTextRepel,
                x = -Inf,
                y = Inf,
                label = eqs[i],
                color = "black",
                bg.color = "white",
                bg.r = 0.1,
                size = 3.5,
                point.size = NA,
                max.overlaps = 100,
                force = 0,
                min.segment.length = Inf,
                hjust = -0.05,
                vjust = vjusts[1:sum(i)],
                parse = TRUE
              )
          }
          if (!is.null(cells.highlight)) {
            cell_df <- subset(p$data, rownames(p$data) %in% cells.highlight)
            if (nrow(cell_df) > 0) {
              # point_size <- p$layers[[1]]$aes_params$size
              if (isTRUE(raster)) {
                p <- p +
                  scattermore::geom_scattermore(
                    data = cell_df,
                    aes(x = .data[[f1]], y = .data[[f2]]),
                    color = cols.highlight,
                    pointsize = floor(sizes.highlight) + stroke.highlight,
                    alpha = alpha.highlight,
                    pixels = raster.dpi
                  ) +
                  scattermore::geom_scattermore(
                    data = cell_df,
                    aes(
                      x = .data[[f1]],
                      y = .data[[f2]],
                      color = .data[[group.by]]
                    ),
                    pointsize = floor(sizes.highlight),
                    alpha = alpha.highlight,
                    pixels = raster.dpi
                  )
              } else {
                p <- p +
                  suppressWarnings(geom_point(
                    data = cell_df,
                    aes(x = .data[[f1]], y = .data[[f2]]),
                    color = cols.highlight,
                    size = sizes.highlight + stroke.highlight,
                    alpha = alpha.highlight
                  )) +
                  suppressWarnings(geom_point(
                    data = cell_df,
                    aes(
                      x = .data[[f1]],
                      y = .data[[f2]],
                      color = .data[[group.by]]
                    ),
                    size = sizes.highlight,
                    alpha = alpha.highlight
                  ))
              }
            }
          }
        }
        if (f1_index > f2_index) {
          label <- paste0(f1, "\n", f2, "\nCor: ", round(pair_sim[f1, f2], 3)) # "\n","f1_index:",f1_index," ","f2_index:",f2_index
          label_pos <- (max(dat_exp[rownames(dat), ], na.rm = TRUE) +
            min(dat_exp[rownames(dat), ], na.rm = TRUE)) /
            2
          fill <- rownames(df_bound)[
            df_bound[, 1] < pair_sim[f1, f2] & df_bound[, 2] >= pair_sim[f1, f2]
          ]
          p <- p +
            annotate(
              geom = "rect",
              xmin = -Inf,
              xmax = Inf,
              ymin = -Inf,
              ymax = Inf,
              fill = fill
            ) +
            annotate(
              geom = GeomTextRepel,
              x = label_pos,
              y = label_pos,
              label = label,
              fontface = "bold",
              color = "black",
              bg.color = "white",
              bg.r = 0.1,
              size = 3.5,
              point.size = NA
            )
        }

        if (f1_index == 1 & f2_index != 1) {
          p <- p +
            theme(
              axis.ticks.y = element_line(),
              axis.text.y = element_text(size = 10)
            )
        }
        if (f2_index == length(features) & f1_index != length(features)) {
          p <- p +
            theme(
              axis.ticks.x = element_line(),
              axis.text.x = element_text(size = 10)
            )
        }
        if (f1_index == 1) {
          p <- p + labs(y = f2) + theme(axis.title.y = element_text(size = 12))
        }
        if (f2_index == length(features)) {
          p <- p + labs(x = f1) + theme(axis.title.x = element_text(size = 12))
        }
        p <- p +
          scale_color_manual(
            name = paste0(group.by, ":"),
            values = colors,
            labels = names(colors)
          ) +
          scale_fill_manual(
            name = paste0(group.by, ":"),
            values = colors,
            labels = names(colors)
          )
        return(p)
      },
      x = order1,
      y = order2,
      SIMPLIFY = FALSE
    )

    legend_list <- NULL
    if (length(features) > 1) {
      legend_list[["correlation"]] <- get_legend(
        ggplot(
          data.frame(range = cor_range, x = 1, y = 1),
          aes(x = x, y = y, fill = range)
        ) +
          geom_point() +
          scale_fill_gradientn(
            name = paste0("Correlation"),
            limits = cor_range,
            n.breaks = 3,
            colors = cor_colors,
            guide = guide_colorbar(
              frame.colour = "black",
              ticks.colour = "black",
              title.hjust = 0
            )
          ) +
          do.call(theme_use, theme_args) +
          theme(
            aspect.ratio = aspect.ratio,
            legend.position = "bottom",
            legend.direction = legend.direction
          )
      )
    }
    if (nlevels(dat[[group.by]]) > 1) {
      legend_list[["group.by"]] <- suppressWarnings(
        get_legend(
          plotlist[[1]] +
            guides(
              fill = guide_legend(
                title.hjust = 0,
                order = 1,
                override.aes = list(size = 4, color = "black", alpha = 1)
              )
            ) +
            do.call(theme_use, theme_args) +
            theme(
              aspect.ratio = aspect.ratio,
              legend.position = "bottom",
              legend.direction = legend.direction
            )
        )
      )
    }

    grob_row <- list()
    plotlist <- suppressWarnings(lapply(plotlist, as_grob))
    for (i in seq(1, length(plotlist), length(features))) {
      grob_row[[paste0(
        i:(i + length(features) - 1),
        collapse = "-"
      )]] <- do.call(cbind, plotlist[i:(i + length(features) - 1)])
    }
    gtable <- do.call(rbind, grob_row)
    if (length(legend_list) > 0) {
      legend_list <- legend_list[!sapply(legend_list, is.null)]
      if (legend.direction == "vertical") {
        legend <- do.call(cbind, legend_list)
      } else {
        legend <- do.call(rbind, legend_list)
      }
      gtable <- add_grob(gtable, legend, legend.position)
    }
    if (nlevels(dat_use[[split.by]]) > 1) {
      split_grob <- grid::textGrob(
        s,
        just = "center",
        gp = grid::gpar(fontface = "bold", fontsize = 13)
      )
      gtable <- add_grob(gtable, split_grob, "top")
    }
    if (!is.null(subtitle)) {
      subtitle_grob <- grid::textGrob(
        subtitle,
        x = 0,
        hjust = 0,
        gp = grid::gpar(fontface = "italic", fontsize = 13)
      )
      gtable <- add_grob(gtable, subtitle_grob, "top")
    }
    if (!is.null(title)) {
      title_grob <- grid::textGrob(title, x = 0, hjust = 0, gp = grid::gpar(fontsize = 14))
      gtable <- add_grob(gtable, title_grob, "top", 2 * grid::grobHeight(title_grob))
    }
    p <- patchwork::wrap_plots(gtable)
    plist[[paste0(s)]] <- p
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
    return(plot)
  } else {
    return(plist)
  }
}
