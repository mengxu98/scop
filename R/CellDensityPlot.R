#' @title Cell density plot
#'
#' @description
#' Feature density by group.
#'
#' @md
#' @inheritParams CellDimPlot
#' @inheritParams FeatureDimPlot
#' @inheritParams scop-params
#' @param features Features to plot.
#' @param flip,reverse Flip the x-axis or reverse the y-axis.
#' @param x_order `"value"` or `"rank"`.
#' @param decreasing Order groups decreasingly.
#' @param cells Cell names to include. `NULL` uses all cells.
#' @param keep_empty Keep empty groups.
#' @param y.nbreaks,y.min,y.max Y-axis breaks and limits. `NULL` limits are automatic.
#' @param same.y.lims Share y-axis limits across panels.
#' @param aspect.ratio Panel aspect ratio.
#' @param force Draw even when there are more than 50 features.
#'
#' @seealso
#' [CellStatPlot]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#' CellDensityPlot(
#'   pancreas_sub,
#'   features = "Sox9",
#'   group.by = "SubCellType"
#' )
#'
#' pancreas_sub <- RunSlingshot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP"
#' )
#'
#' CellDensityPlot(
#'   pancreas_sub,
#'   features = "Lineage1",
#'   group.by = "SubCellType",
#'   aspect.ratio = 1
#' )
#'
#' CellDensityPlot(
#'   pancreas_sub,
#'   features = "Lineage1",
#'   group.by = "SubCellType",
#'   flip = TRUE
#' )
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
  palette = "Chinese",
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
  force = FALSE,
  verbose = TRUE
) {
  check_r("ggridges", verbose = FALSE)
  assay <- assay %||% DefaultAssay(srt)
  x_order <- match.arg(x_order)
  if (is.null(features)) {
    log_message(
      "{.arg features} must be provided",
      message_type = "error"
    )
  }
  if (!inherits(features, "character")) {
    log_message(
      "{.arg features} is not a character vectors",
      message_type = "error"
    )
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
      log_message(
        "{.val {i}} is not in the meta.data of object",
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

  features <- unique(features)
  features_drop <- features[
    !features %in% c(rownames(srt@assays[[assay]]), colnames(srt@meta.data))
  ]
  if (length(features_drop) > 0) {
    log_message(
      "{.val {features_drop}} are not in the features of srt",
      message_type = "warning",
      verbose = verbose
    )
    features <- features[!features %in% features_drop]
  }

  features_gene <- features[features %in% rownames(srt@assays[[assay]])]
  features_meta <- features[features %in% colnames(srt@meta.data)]
  if (length(intersect(features_gene, features_meta)) > 0) {
    log_message(
      "Features appear in both gene names and metadata names: {.val {intersect(features_gene, features_meta)}}",
      message_type = "warning",
      verbose = verbose
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
    dat_meta <- as_matrix(srt@meta.data[, features_meta, drop = FALSE])
  } else {
    dat_meta <- matrix(nrow = ncol(srt@assays[[1]]), ncol = 0)
  }
  dat_exp <- cbind(dat_gene, dat_meta)
  features <- unique(features[features %in% c(features_gene, features_meta)])

  if (!is.numeric(dat_exp) && !inherits(dat_exp, "Matrix")) {
    log_message(
      "{.arg features} must be type of numeric variable",
      message_type = "error"
    )
  }
  if (length(features) > 50 && isFALSE(force)) {
    log_message(
      "More than 50 {.arg features} to be plotted",
      message_type = "warning",
      verbose = verbose
    )
    answer <- utils::askYesNo("Are you sure to continue?", default = FALSE)
    if (isFALSE(answer)) {
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
      as_matrix(dat_exp[
        ,
        features
      ])[is.finite(as_matrix(dat_exp[, features]))],
      na.rm = TRUE
    )
  }
  if (isTRUE(same.y.lims) && is.null(y.min)) {
    y.min <- min(
      as_matrix(dat_exp[
        ,
        features
      ])[is.finite(as_matrix(dat_exp[, features]))],
      na.rm = TRUE
    )
  }

  plist <- list()
  for (f in features) {
    for (g in group.by) {
      colors <- palette_colors(
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

  combine_plot_list(
    plist,
    combine = combine, nrow = nrow, ncol = ncol, byrow = byrow
  )
}
