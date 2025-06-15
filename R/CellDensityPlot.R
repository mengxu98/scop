#' CellDensityPlot
#'
#' Plots the density of specified features in a single or multiple groups,
#' grouped by specified variables.
#'
#' @param srt A Seurat object.
#' @param features A character vector specifying the features to plot.
#' @param group.by A character vector specifying the variables to group the data by.
#' @param split.by A character vector specifying the variables to split the data by.
#' Default is NULL, which means no splitting is performed.
#' @param assay A character specifying the assay to use from the Seurat object.
#'   Default is NULL, which means the default assay will be used.
#' @param layer A character specifying the layer to use from the assay. Default is "data".
#' @param flip A logical indicating whether to flip the x-axis. Default is FALSE.
#' @param reverse A logical indicating whether to reverse the y-axis. Default is FALSE.
#' @param x_order A character specifying how to order the x-axis. Can be "value" or "rank". Default is "value".
#' @param decreasing A logical indicating whether to order the groups in decreasing order. Default is NULL.
#' @param palette A character specifying the color palette to use for grouping variables. Default is "Paired".
#' @param palcolor A character specifying the color to use for each group. Default is NULL.
#' @param cells A character vector specifying the cells to plot. Default is NULL, which means all cells are included.
#' @param keep_empty A logical indicating whether to keep empty groups. Default is FALSE.
#' @param y.nbreaks An integer specifying the number of breaks on the y-axis. Default is 4.
#' @param y.min A numeric specifying the minimum value on the y-axis. Default is NULL, which means the minimum value will be automatically determined.
#' @param y.max A numeric specifying the maximum value on the y-axis. Default is NULL, which means the maximum value will be automatically determined.
#' @param same.y.lims A logical indicating whether to use the same y-axis limits for all plots. Default is FALSE.
#' @param aspect.ratio A numeric specifying the aspect ratio of the plot. Default is NULL, which means the aspect ratio will be automatically determined.
#' @param title A character specifying the title of the plot. Default is NULL.
#' @param subtitle A character specifying the subtitle of the plot. Default is NULL.
#' @param legend.position A character specifying the position of the legend. Default is "right".
#' @param legend.direction A character specifying the direction of the legend. Default is "vertical".
#' @param theme_use A character specifying the theme to use. Default is "theme_scop".
#' @param theme_args A list of arguments to pass to the theme function.
#' @param combine A logical indicating whether to combine multiple plots into a single plot. Default is TRUE.
#' @param nrow An integer specifying the number of rows in the combined plot.
#'   Default is NULL, which means determined automatically based on the number of plots.
#' @param ncol An integer specifying the number of columns in the combined plot.
#'   Default is NULL, which means determined automatically based on the number of plots.
#' @param byrow A logical indicating whether to add plots by row or by column in the combined plot. Default is TRUE.
#' @param force A logical indicating whether to continue plotting if there are more than 50 features. Default is FALSE.
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' data("pancreas_sub")
#' CellDensityPlot(
#'   pancreas_sub,
#'   features = "Sox9",
#'   group.by = "SubCellType"
#' )
#'
#' pancreas_sub <- RunSlingshot( # bug in methods::validObject(out)
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP"
#' )
#' CellDensityPlot(
#'   pancreas_sub,
#'   features = "Lineage1",
#'   group.by = "SubCellType",
#'   aspect.ratio = 1
#' )
#' CellDensityPlot(
#'   pancreas_sub,
#'   features = "Lineage1",
#'   group.by = "SubCellType",
#'   flip = TRUE
#' )
#' }
CellDensityPlot <- function(
  srt,
  features,
  group.by = NULL,
  split.by = NULL,
  assay = NULL,
  layer = "data",
  flip = FALSE,
  reverse = FALSE,
  x_order = c("value", "rank"),
  decreasing = NULL,
  palette = "Paired",
  palcolor = NULL,
  cells = NULL,
  keep_empty = FALSE,
  y.nbreaks = 4,
  y.min = NULL,
  y.max = NULL,
  same.y.lims = FALSE,
  aspect.ratio = NULL,
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
  force = FALSE
) {
  check_r("ggridges")
  assay <- assay %||% DefaultAssay(srt)
  x_order <- match.arg(x_order)
  if (is.null(features)) {
    stop("'features' must be provided.")
  }
  if (!inherits(features, "character")) {
    stop("'features' is not a character vectors")
  }
  if (is.null(group.by)) {
    group.by <- "All.groups"
    srt@meta.data[[group.by]] <- factor("")
  }
  if (is.null(split.by)) {
    split.by <- "All.groups"
    srt@meta.data[[split.by]] <- factor("")
  }
  if (group.by == split.by & group.by == "All.groups") {
    legend.position <- "none"
  }
  for (i in c(group.by, split.by)) {
    if (!i %in% colnames(srt@meta.data)) {
      stop(paste0(i, " is not in the meta.data of srt object."))
    }
    if (!is.factor(srt@meta.data[[i]])) {
      srt@meta.data[[i]] <- factor(
        srt@meta.data[[i]],
        levels = unique(srt@meta.data[[i]])
      )
    }
  }

  features <- unique(features)
  features_drop <- features[
    !features %in% c(rownames(srt@assays[[assay]]), colnames(srt@meta.data))
  ]
  # print(colnames(srt@meta.data))
  if (length(features_drop) > 0) {
    warning(
      paste0(features_drop, collapse = ","),
      " are not in the features of srt.",
      immediate. = TRUE
    )
    features <- features[!features %in% features_drop]
  }

  features_gene <- features[features %in% rownames(srt@assays[[assay]])]
  features_meta <- features[features %in% colnames(srt@meta.data)]
  if (length(intersect(features_gene, features_meta)) > 0) {
    warning(
      "Features appear in both gene names and metadata names: ",
      paste0(intersect(features_gene, features_meta), collapse = ",")
    )
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
    dat_meta <- Matrix::as.matrix(srt@meta.data[, features_meta, drop = FALSE])
  } else {
    dat_meta <- matrix(nrow = ncol(srt@assays[[1]]), ncol = 0)
  }
  dat_exp <- cbind(dat_gene, dat_meta)
  features <- unique(features[features %in% c(features_gene, features_meta)])

  if (!is.numeric(dat_exp) && !inherits(dat_exp, "Matrix")) {
    stop("'features' must be type of numeric variable.")
  }
  if (length(features) > 50 && !isTRUE(force)) {
    warning("More than 50 features to be plotted", immediate. = TRUE)
    answer <- utils::askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
      return(invisible(NULL))
    }
  }

  dat_use <- cbind(
    dat_exp,
    srt@meta.data[row.names(dat_exp), c(group.by, split.by), drop = FALSE]
  )
  if (!is.null(cells)) {
    dat_use <- dat_use[intersect(rownames(dat_use), cells), , drop = FALSE]
  }

  if (isTRUE(same.y.lims) && is.null(y.max)) {
    y.max <- max(
      Matrix::as.matrix(dat_exp[,
        features
      ])[is.finite(Matrix::as.matrix(dat_exp[, features]))],
      na.rm = TRUE
    )
  }
  if (isTRUE(same.y.lims) && is.null(y.min)) {
    y.min <- min(
      Matrix::as.matrix(dat_exp[,
        features
      ])[is.finite(Matrix::as.matrix(dat_exp[, features]))],
      na.rm = TRUE
    )
  }

  plist <- list()
  for (f in features) {
    for (g in group.by) {
      colors <- palette_scop(
        levels(dat_use[[g]]),
        palette = palette,
        palcolor = palcolor
      )
      for (s in levels(dat_use[[split.by]])) {
        dat <- dat_use[dat_use[[split.by]] == s, , drop = FALSE]
        if (any(is.infinite(dat[, f]))) {
          dat[, f][dat[, f] == max(dat[, f], na.rm = TRUE)] <- max(
            dat[, f][is.finite(dat[, f])],
            na.rm = TRUE
          )
          dat[, f][dat[, f] == min(dat[, f], na.rm = TRUE)] <- min(
            dat[, f][is.finite(dat[, f])],
            na.rm = TRUE
          )
        }
        dat[, "cell"] <- rownames(dat)
        if (x_order == "value") {
          dat[, "value"] <- dat[, f]
        } else {
          dat[, "value"] <- rank(dat[, f])
        }
        dat[, "features"] <- f
        dat[, "split.by"] <- s
        dat <- dat[!is.na(dat[[f]]), , drop = FALSE]

        y_max_use <- y.max %||%
          suppressWarnings(max(
            dat[, "value"][is.finite(x = dat[, "value"])],
            na.rm = TRUE
          ))
        y_min_use <- y.min %||%
          suppressWarnings(min(
            dat[, "value"][is.finite(x = dat[, "value"])],
            na.rm = TRUE
          ))

        if (!is.null(decreasing)) {
          levels <- dat |>
            dplyr::group_by_at(g) |>
            dplyr::summarise_at(
              .funs = stats::median,
              .vars = f,
              na.rm = TRUE
            ) |>
            dplyr::arrange_at(
              .vars = f,
              .funs = if (decreasing) dplyr::desc else list()
            ) |>
            dplyr::pull(g) |>
            as.character()
          dat[["order"]] <- factor(dat[[g]], levels = levels)
        } else {
          dat[["order"]] <- factor(dat[[g]], levels = rev(levels(dat[[g]])))
        }
        if (flip) {
          dat[["order"]] <- factor(dat[[g]], levels = levels(dat[[g]]))
          aspect.ratio <- 1 / aspect.ratio
          if (length(aspect.ratio) == 0 || is.na(aspect.ratio)) {
            aspect.ratio <- NULL
          }
        }
        p <- ggplot(
          dat,
          aes(x = .data[["value"]], y = .data[["order"]], fill = .data[[g]])
        ) +
          ggridges::geom_density_ridges()
        p <- p +
          scale_fill_manual(
            name = paste0(g, ":"),
            values = colors
          )
        y.trans <- ifelse(flip, "reverse", "identity")
        y.trans <- ifelse(
          reverse,
          setdiff(c("reverse", "identity"), y.trans),
          y.trans
        )

        limits <- if (y.trans == "reverse") {
          c(y_max_use, y_min_use)
        } else {
          c(y_min_use, y_max_use)
        }
        p <- p +
          scale_y_discrete(drop = !keep_empty, expand = c(0, 0)) +
          scale_x_continuous(
            limits = limits,
            trans = y.trans,
            n.breaks = y.nbreaks,
            expand = c(0, 0)
          )
        if (split.by != "All.groups") {
          p <- p + facet_grid(. ~ split.by)
        }
        p <- p + labs(title = title, subtitle = subtitle, x = f, y = g)
        if (isTRUE(flip)) {
          p <- p +
            do.call(theme_use, theme_args) +
            theme(
              aspect.ratio = aspect.ratio,
              strip.text.x = element_text(angle = 0),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.ticks.x = element_line(),
              panel.grid.major.x = element_line(color = "grey", linetype = 2),
              legend.position = legend.position,
              legend.direction = legend.direction
            ) +
            coord_flip()
        } else {
          p <- p +
            do.call(theme_use, theme_args) +
            theme(
              aspect.ratio = aspect.ratio,
              strip.text.y = element_text(angle = 0),
              axis.text.x = element_text(),
              axis.text.y = element_text(hjust = 1),
              axis.ticks.y = element_line(),
              panel.grid.major.y = element_line(color = "grey", linetype = 2),
              legend.position = legend.position,
              legend.direction = legend.direction
            )
        }
        plist[[paste0(f, ":", g, ":", paste0(s, collapse = ","))]] <- p
      }
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
    return(plot)
  } else {
    return(plist)
  }
  return(p)
}
