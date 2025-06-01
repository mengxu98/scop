#' Statistical plot of cells
#'
#' @inheritParams StatPlot
#' @param srt A Seurat object.
#' @param cells A character vector specifying the cells to include in the plot. Default is NULL.
#'
#' @seealso \code{\link{StatPlot}}
#'
#' @export
#'
#' @examples
#' data("pancreas_sub")
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "SubCellType",
#'   label = TRUE
#' )
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "SubCellType",
#'   label = TRUE
#' ) %>%
#'   panel_fix(height = 2, width = 3)
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "SubCellType",
#'   stat_type = "count",
#'   position = "dodge",
#'   label = TRUE
#' )
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "SubCellType",
#'   bg.by = "CellType",
#'   palette = "Set1",
#'   stat_type = "count",
#'   position = "dodge"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   plot_type = "bar"
#' )
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   plot_type = "rose"
#' )
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   plot_type = "ring"
#' )
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   plot_type = "pie"
#' )
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   plot_type = "dot"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   plot_type = "bar"
#' )
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   plot_type = "rose"
#' )
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   plot_type = "ring"
#' )
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   plot_type = "area"
#' )
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   plot_type = "dot"
#' )
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   plot_type = "trend"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   plot_type = "bar",
#'   individual = TRUE
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   stat_type = "count",
#'   plot_type = "bar"
#' )
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   stat_type = "count",
#'   plot_type = "rose"
#' )
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   stat_type = "count",
#'   plot_type = "ring"
#' )
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   stat_type = "count",
#'   plot_type = "area"
#' )
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   stat_type = "count",
#'   plot_type = "dot"
#' )
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   stat_type = "count",
#'   plot_type = "trend"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   stat_type = "count",
#'   plot_type = "bar",
#'   position = "dodge",
#'   label = TRUE
#' )
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   stat_type = "count",
#'   plot_type = "rose",
#'   position = "dodge",
#'   label = TRUE
#' )
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   stat_type = "count",
#'   plot_type = "ring",
#'   position = "dodge",
#'   label = TRUE
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = c("CellType", "Phase"),
#'   plot_type = "sankey"
#' )
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = c("CellType", "Phase"),
#'   plot_type = "chord"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = c("CellType", "Phase"),
#'   plot_type = "venn",
#'   stat_level = list(
#'     CellType = c("Ductal", "Ngn3 low EP"),
#'     Phase = "S"
#'   )
#' )
#' pancreas_sub$Progenitor <- pancreas_sub$CellType %in% c("Ngn3 low EP", "Ngn3 high EP")
#' pancreas_sub$G2M <- pancreas_sub$Phase == "G2M"
#' pancreas_sub$Sox9_Expressed <- pancreas_sub[["RNA"]]@counts["Sox9", ] > 0
#' pancreas_sub$Neurog3_Expressed <- pancreas_sub[["RNA"]]@counts["Neurog3", ] > 0
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = c(
#'     "Progenitor", "G2M", "Sox9_Expressed", "Neurog3_Expressed"
#'   ),
#'   plot_type = "venn",
#'   stat_level = "TRUE"
#' )
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = c(
#'     "Progenitor", "G2M", "Sox9_Expressed", "Neurog3_Expressed"
#'   ),
#'   plot_type = "upset",
#'   stat_level = "TRUE"
#' )
#' sum(
#'   pancreas_sub$Progenitor == "FALSE" &
#'     pancreas_sub$G2M == "FALSE" &
#'     pancreas_sub$Sox9_Expressed == "TRUE" &
#'     pancreas_sub$Neurog3_Expressed == "FALSE"
#' )
CellStatPlot <- function(
    srt,
    stat.by,
    group.by = NULL,
    split.by = NULL,
    bg.by = NULL,
    cells = NULL,
    flip = FALSE,
    NA_color = "grey",
    NA_stat = TRUE,
    keep_empty = FALSE,
    individual = FALSE,
    stat_level = NULL,
    plot_type = c(
      "bar",
      "rose",
      "ring",
      "pie",
      "trend",
      "area",
      "dot",
      "sankey",
      "chord",
      "venn",
      "upset"
    ),
    stat_type = c("percent", "count"),
    position = c("stack", "dodge"),
    palette = "Paired",
    palcolor = NULL,
    alpha = 1,
    bg_palette = "Paired",
    bg_palcolor = NULL,
    bg_alpha = 0.2,
    label = FALSE,
    label.size = 3.5,
    label.fg = "black",
    label.bg = "white",
    label.bg.r = 0.1,
    aspect.ratio = NULL,
    title = NULL,
    subtitle = NULL,
    xlab = NULL,
    ylab = NULL,
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
  cells <- cells %||% colnames(srt@assays[[1]])
  meta.data <- srt@meta.data[cells, , drop = FALSE]

  plot <- StatPlot(
    meta.data = meta.data,
    stat.by = stat.by,
    group.by = group.by,
    split.by = split.by,
    bg.by = bg.by,
    flip = flip,
    NA_color = NA_color,
    NA_stat = NA_stat,
    keep_empty = keep_empty,
    individual = individual,
    stat_level = stat_level,
    plot_type = plot_type,
    stat_type = stat_type,
    position = position,
    palette = palette,
    palcolor = palcolor,
    alpha = alpha,
    bg_palette = bg_palette,
    bg_palcolor = bg_palcolor,
    bg_alpha = bg_alpha,
    label = label,
    label.size = label.size,
    label.fg = label.fg,
    label.bg = label.bg,
    label.bg.r = label.bg.r,
    aspect.ratio = aspect.ratio,
    title = title,
    subtitle = subtitle,
    xlab = xlab,
    ylab = ylab,
    legend.position = legend.position,
    legend.direction = legend.direction,
    theme_use = theme_use,
    theme_args = theme_args,
    combine = combine,
    nrow = nrow,
    ncol = ncol,
    byrow = byrow,
    force = force,
    seed = seed
  )

  return(plot)
}

#' StatPlot
#'
#' Visualizes data using various plot types such as bar plots,
#' rose plots, ring plots, pie charts, trend plots, area plots,
#' dot plots, sankey plots, chord plots, venn diagrams, and upset plots.
#'
#' @param meta.data The data frame containing the data to be plotted.
#' @param stat.by The column name(s) in \code{meta.data} specifying the variable(s) to be plotted.
#' @param group.by The column name in \code{meta.data} specifying the grouping variable.
#' @param split.by The column name in \code{meta.data} specifying the splitting variable.
#' @param bg.by The column name in \code{meta.data} specifying the background variable for bar plots.
#' @param flip Logical indicating whether to flip the plot.
#' @param NA_color The color to use for missing values.
#' @param NA_stat Logical indicating whether to include missing values in the plot.
#' @param keep_empty Logical indicating whether to keep empty groups in the plot.
#' @param individual Logical indicating whether to plot individual groups separately.
#' @param stat_level The level(s) of the variable(s) specified in \code{stat.by} to include in the plot.
#' @param plot_type The type of plot to create.
#' Can be one of "bar", "rose", "ring", "pie", "trend", "area", "dot", "sankey", "chord", "venn", or "upset".
#' @param stat_type The type of statistic to compute for the plot.
#' Can be one of "percent" or "count".
#' @param position The position adjustment for the plot.
#' Can be one of "stack" or "dodge".
#' @param palette The name of the color palette to use for the plot.
#' @param palcolor The color to use in the color palette.
#' @param alpha The transparency level for the plot.
#' @param bg_palette The name of the background color palette to use for bar plots.
#' @param bg_palcolor The color to use in the background color palette.
#' @param bg_alpha The transparency level for the background color in bar plots.
#' @param label Logical indicating whether to add labels on the plot.
#' @param label.size The size of the labels.
#' @param label.fg The foreground color of the labels.
#' @param label.bg The background color of the labels.
#' @param label.bg.r The radius of the rounded corners of the label background.
#' @param aspect.ratio The aspect ratio of the plot.
#' @param title The main title of the plot.
#' @param subtitle The subtitle of the plot.
#' @param xlab The x-axis label of the plot.
#' @param ylab The y-axis label of the plot.
#' @param legend.position The position of the legend in the plot.
#' Can be one of "right", "left", "bottom", "top", or "none".
#' @param legend.direction The direction of the legend in the plot. Can be one of "vertical" or "horizontal".
#' @param theme_use The name of the theme to use for the plot. Can be one of the predefined themes or a custom theme.
#' @param theme_args A list of arguments to be passed to the theme function.
#' @param combine Logical indicating whether to combine multiple plots into a single plot.
#' @param nrow The number of rows in the combined plot.
#' @param ncol The number of columns in the combined plot.
#' @param byrow Logical indicating whether to fill the plot by row or by column.
#' @param force Logical indicating whether to force the plot even if some variables have more than 100 levels.
#' @param seed The random seed to use for reproducible results.
#'
#' @seealso \code{\link{CellStatPlot}}
#'
#' @export
#'
#' @examples
#' data("pancreas_sub")
#' head(pancreas_sub@meta.data)
#' StatPlot(
#'   pancreas_sub@meta.data,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   plot_type = "bar",
#'   label = TRUE
#' )
#'
#' head(pancreas_sub[["RNA"]]@meta.features)
#' StatPlot(
#'   pancreas_sub[["RNA"]]@meta.features,
#'   stat.by = "highly_variable_genes",
#'   plot_type = "ring",
#'   label = TRUE
#' )
#'
#' pancreas_sub <- AnnotateFeatures(
#'   pancreas_sub,
#'   species = "Mus_musculus",
#'   IDtype = "symbol",
#'   db = "GeneType"
#' )
#' head(pancreas_sub[["RNA"]]@meta.features)
#' StatPlot(
#'   pancreas_sub[["RNA"]]@meta.features,
#'   stat.by = "highly_variable_genes",
#'   group.by = "GeneType",
#'   stat_type = "count",
#'   plot_type = "bar",
#'   position = "dodge",
#'   label = TRUE,
#'   NA_stat = FALSE
#' )
StatPlot <- function(
    meta.data,
    stat.by,
    group.by = NULL,
    split.by = NULL,
    bg.by = NULL,
    flip = FALSE,
    NA_color = "grey",
    NA_stat = TRUE,
    keep_empty = FALSE,
    individual = FALSE,
    stat_level = NULL,
    plot_type = c(
      "bar",
      "rose",
      "ring",
      "pie",
      "trend",
      "area",
      "dot",
      "sankey",
      "chord",
      "venn",
      "upset"
    ),
    stat_type = c("percent", "count"),
    position = c("stack", "dodge"),
    palette = "Paired",
    palcolor = NULL,
    alpha = 1,
    bg_palette = "Paired",
    bg_palcolor = NULL,
    bg_alpha = 0.2,
    label = FALSE,
    label.size = 3.5,
    label.fg = "black",
    label.bg = "white",
    label.bg.r = 0.1,
    aspect.ratio = NULL,
    title = NULL,
    subtitle = NULL,
    xlab = NULL,
    ylab = NULL,
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

  stat_type <- match.arg(stat_type)
  plot_type <- match.arg(plot_type)
  position <- match.arg(position)

  if (nrow(meta.data) == 0) {
    stop("meta.data is empty.")
  }
  if (is.null(group.by)) {
    group.by <- "All.groups"
    xlab <- ""
    meta.data[[group.by]] <- factor("")
  }
  if (is.null(split.by)) {
    split.by <- "All.groups"
    meta.data[[split.by]] <- factor("")
  }

  for (i in unique(c(group.by, split.by, bg.by))) {
    if (!i %in% colnames(meta.data)) {
      stop(paste0(i, " is not in the meta.data."))
    }
    if (!is.factor(meta.data[[i]])) {
      meta.data[[i]] <- factor(meta.data[[i]], levels = unique(meta.data[[i]]))
    }
  }

  bg_map <- NULL
  if (!is.null(bg.by)) {
    for (g in group.by) {
      df_table <- table(meta.data[[g]], meta.data[[bg.by]])
      if (max(rowSums(df_table > 0), na.rm = TRUE) > 1) {
        stop("'group.by' must be a part of 'bg.by'")
      } else {
        bg_map[[g]] <- stats::setNames(
          colnames(df_table)[apply(df_table, 1, function(x) which(x > 0))],
          rownames(df_table)
        )
      }
    }
  } else {
    for (g in group.by) {
      bg_map[[g]] <- stats::setNames(levels(meta.data[[g]]), levels(meta.data[[g]]))
    }
  }

  for (i in unique(stat.by)) {
    if (!i %in% colnames(meta.data)) {
      stop(paste0(i, " is not in the meta.data."))
    }
    if (plot_type %in% c("venn", "upset")) {
      if (!is.factor(meta.data[[i]]) && !is.logical(meta.data[[i]])) {
        meta.data[[i]] <- factor(
          meta.data[[i]],
          levels = unique(meta.data[[i]])
        )
      }
    } else if (!is.factor(meta.data[[i]])) {
      meta.data[[i]] <- factor(meta.data[[i]], levels = unique(meta.data[[i]]))
    }
  }

  if (length(stat.by) >= 2) {
    if (!plot_type %in% c("sankey", "chord", "venn", "upset")) {
      stop(
        "plot_type must be one of 'sankey', 'chord', 'venn' and 'upset' whtn multiple 'stat.by' provided."
      )
    }
    if (length(stat.by) > 2 && plot_type == "chord") {
      stop(
        "'stat.by' can only be a vector of length 2 when 'plot_type' is 'chord'."
      )
    }
    if (length(stat.by) > 7 && plot_type == "venn") {
      stop(
        "'stat.by' can only be a vector of length <= 7 when 'plot_type' is 'venn'."
      )
    }
  }

  levels <- unique(
    unlist(
      lapply(
        meta.data[, stat.by, drop = FALSE],
        function(x) {
          if (is.factor(x)) {
            return(levels(x))
          }
          if (is.logical(x)) {
            return(as.character(unique(x)))
          }
        }
      )
    )
  )

  if (plot_type %in% c("venn", "upset")) {
    if (is.null(stat_level)) {
      stat_level <- lapply(stat.by, function(stat) {
        levels(meta.data[[stat]])[1] %||% sort(unique(meta.data[[stat]]))[1]
      })
      message("stat_level is set to ", paste0(stat_level, collapse = ","))
    } else {
      if (length(stat_level) == 1) {
        stat_level <- rep(stat_level, length(stat.by))
      }
      if (length(stat_level) != length(stat.by)) {
        stop("'stat_level' must be of length 1 or the same length as 'stat.by'")
      }
    }
    if (is.null(names(stat_level))) {
      names(stat_level) <- stat.by
    }
    for (i in stat.by) {
      meta.data[[i]] <- meta.data[[i]] %in% stat_level[[i]]
    }
  }

  if (plot_type %in% c("rose", "ring", "pie")) {
    aspect.ratio <- 1
  }

  if (
    any(group.by != "All.groups") &&
      plot_type %in% c("sankey", "chord", "venn", "upset")
  ) {
    warning(
      "group.by is not used when plot sankey, chord, venn or upset",
      immediate. = TRUE
    )
  }
  if (
    stat_type == "percent" &&
      plot_type %in% c("sankey", "chord", "venn", "upset")
  ) {
    warning(
      "stat_type is forcibly set to 'count' when plot sankey, chord, venn or upset",
      immediate. = TRUE
    )
    stat_type <- "count"
  }

  dat_all <- meta.data[,
    unique(c(stat.by, group.by, split.by, bg.by)),
    drop = FALSE
  ]
  nlev <- sapply(dat_all, nlevels)
  nlev <- nlev[nlev > 100]
  if (length(nlev) > 0 && !isTRUE(force)) {
    warning(
      paste(names(nlev), sep = ","),
      " have more than 100 levels.",
      immediate. = TRUE
    )
    answer <- utils::askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
      return(invisible(NULL))
    }
  }
  dat_split <- split.data.frame(dat_all, dat_all[[split.by]])

  plist <- list()
  if (plot_type %in% c("bar", "rose", "ring", "pie", "trend", "area", "dot")) {
    xlab <- xlab %||% group.by
    ylab <- ylab %||% ifelse(stat_type == "count", "Count", "Percentage")
    if (identical(theme_use, "theme_blank")) {
      theme_args[["xlab"]] <- xlab
      theme_args[["ylab"]] <- ylab
      if (plot_type %in% c("rose", "ring", "pie")) {
        theme_args[["add_coord"]] <- FALSE
      }
    }
    colors <- palette_scop(
      dat_all[[stat.by]],
      palette = palette,
      palcolor = palcolor,
      NA_color = NA_color,
      NA_keep = TRUE
    )

    comb_list <- list()
    comb <- expand.grid(
      stat_name = stat.by,
      group_name = group.by,
      stringsAsFactors = FALSE
    )
    if (isTRUE(individual)) {
      for (g in group.by) {
        comb_list[[g]] <- merge(
          comb,
          expand.grid(
            group_name = g,
            group_element = levels(dat_all[[g]]),
            split_name = levels(dat_all[[split.by]]),
            stringsAsFactors = FALSE
          ),
          by = "group_name"
        )
      }
    } else {
      for (g in group.by) {
        comb_list[[g]] <- merge(
          comb,
          expand.grid(
            group_name = g,
            group_element = list(levels(dat_all[[g]])),
            split_name = levels(dat_all[[split.by]]),
            stringsAsFactors = FALSE
          ),
          by = "group_name"
        )
      }
    }
    comb <- do.call(rbind, comb_list)
    rownames(comb) <- paste0(
      comb[["group_name"]],
      ":",
      sapply(comb[["group_element"]], function(x) paste0(x, collapse = ",")),
      ":",
      comb[["split_name"]]
    )

    plist <- lapply(
      stats::setNames(rownames(comb), rownames(comb)), function(i) {
        stat.by <- comb[i, "stat_name"]
        sp <- comb[i, "split_name"]
        g <- comb[i, "group_name"]
        single_group <- comb[[i, "group_element"]]
        colors_use <- colors[
          names(colors) %in%
            dat_split[[ifelse(split.by == "All.groups", 1, sp)]][[stat.by]]
        ]
        if (
          any(is.na(dat_split[[ifelse(split.by == "All.groups", 1, sp)]][[
            stat.by
          ]])) &&
            isTRUE(NA_stat)
        ) {
          colors_use <- c(colors_use, colors["NA"])
        }
        if (stat_type == "percent") {
          dat_use <- dat_split[[ifelse(split.by == "All.groups", 1, sp)]] %>%
            stats::xtabs(
              formula = paste0("~", stat.by, "+", g),
              addNA = NA_stat
            ) %>%
            as.data.frame() %>%
            dplyr::group_by(
              dplyr::across(
                dplyr::all_of(g)
              ),
              .drop = FALSE
            ) %>%
            dplyr::mutate(groupn = sum(Freq)) %>%
            dplyr::group_by(
              dplyr::across(
                dplyr::all_of(
                  c(stat.by, g)
                )
              ),
              .drop = FALSE
            ) %>%
            dplyr::mutate(value = Freq / groupn) %>%
            as.data.frame()
        } else {
          dat_use <- dat_split[[ifelse(split.by == "All.groups", 1, sp)]] %>%
            stats::xtabs(
              formula = paste0("~", stat.by, "+", g),
              addNA = NA_stat
            ) %>%
            as.data.frame() %>%
            dplyr::mutate(value = Freq)
        }
        dat <- dat_use[dat_use[[g]] %in% single_group, , drop = FALSE]
        dat[[g]] <- factor(
          dat[[g]],
          levels = levels(dat[[g]])[levels(dat[[g]]) %in% dat[[g]]]
        )
        dat <- dat[!is.na(dat[["value"]]), , drop = FALSE]
        if (!is.null(bg.by)) {
          bg <- bg.by
          bg_color <- palette_scop(
            levels(dat_all[[bg]]),
            palette = bg_palette,
            palcolor = bg_palcolor
          )
        } else {
          bg <- g
          bg_color <- palette_scop(
            levels(dat_all[[bg]]),
            palcolor = bg_palcolor %||%
              rep(c("transparent", "grey85"), nlevels(dat_all[[bg]]))
          )
        }

        if (isTRUE(flip)) {
          dat[[g]] <- factor(dat[[g]], levels = rev(levels(dat[[g]])))
          aspect.ratio <- 1 / aspect.ratio
          if (length(aspect.ratio) == 0 || is.na(aspect.ratio)) {
            aspect.ratio <- NULL
          }
        }
        if (plot_type == "ring") {
          dat[[g]] <- factor(dat[[g]], levels = c("   ", levels(dat[[g]])))
          dat <- rbind(dat, dat[nrow(dat) + 1, , drop = FALSE])
          dat[nrow(dat), g] <- "   "
        }
        if (plot_type == "dot") {
          position_use <- position_identity()
          scalex <- scale_x_discrete(drop = !keep_empty)
        } else {
          if (position == "stack") {
            position_use <- position_stack(vjust = 0.5)
            scalex <- scale_x_discrete(drop = !keep_empty, expand = c(0, 0))
            scaley <- scale_y_continuous(
              labels = if (stat_type == "count") {
                scales::number
              } else {
                scales::percent
              },
              expand = c(0, 0)
            )
          } else if (position == "dodge") {
            if (plot_type == "area") {
              position_use <- position_dodge2(width = 0.9, preserve = "total")
            } else {
              position_use <- position_dodge2(width = 0.9, preserve = "single")
            }
            scalex <- scale_x_discrete(drop = !keep_empty)
            scaley <- scale_y_continuous(
              limits = c(0, max(dat[["value"]], na.rm = TRUE) * 1.1),
              labels = if (stat_type == "count") {
                scales::number
              } else {
                scales::percent
              },
              expand = c(0, 0)
            )
          }
        }
        if (position == "stack") {
          bg_layer <- NULL
        } else {
          bg_data <- stats::na.omit(unique(dat[, g, drop = FALSE]))
          bg_data[["x"]] <- as.numeric(bg_data[[g]])
          bg_data[["xmin"]] <- ifelse(
            bg_data[["x"]] == min(bg_data[["x"]]),
            -Inf,
            bg_data[["x"]] - 0.5
          )
          bg_data[["xmax"]] <- ifelse(
            bg_data[["x"]] == max(bg_data[["x"]]),
            Inf,
            bg_data[["x"]] + 0.5
          )
          bg_data[["ymin"]] <- -Inf
          bg_data[["ymax"]] <- Inf
          bg_data[["fill"]] <- bg_color[bg_map[[g]][as.character(bg_data[[g]])]]
          bg_layer <- geom_rect(
            data = bg_data,
            xmin = bg_data[["xmin"]],
            xmax = bg_data[["xmax"]],
            ymin = bg_data[["ymin"]],
            ymax = bg_data[["ymax"]],
            fill = bg_data[["fill"]],
            alpha = bg_alpha,
            inherit.aes = FALSE
          )
        }

        if (plot_type == "bar") {
          p <- ggplot(
            dat,
            aes(x = .data[[g]], y = value, group = .data[[stat.by]])
          ) +
            bg_layer +
            geom_col(
              aes(fill = .data[[stat.by]]),
              width = 0.8,
              color = "black",
              alpha = alpha,
              position = position_use
            ) +
            scalex +
            scaley
        }
        if (plot_type == "trend") {
          dat_area <- dat[rep(seq_len(nrow(dat)), each = 2), , drop = FALSE]
          dat_area[[g]] <- as.numeric(dat_area[[g]])
          dat_area[seq(1, nrow(dat_area), 2), g] <- dat_area[
            seq(1, nrow(dat_area), 2),
            g
          ] -
            0.3
          dat_area[seq(2, nrow(dat_area), 2), g] <- dat_area[
            seq(2, nrow(dat_area), 2),
            g
          ] +
            0.3
          p <- ggplot(
            dat,
            aes(x = .data[[g]], y = value, fill = .data[[stat.by]])
          ) +
            bg_layer +
            geom_area(
              data = dat_area,
              mapping = aes(x = .data[[g]], fill = .data[[stat.by]]),
              alpha = alpha / 2,
              color = "grey50",
              position = position_use
            ) +
            geom_col(
              aes(fill = .data[[stat.by]]),
              width = 0.6,
              color = "black",
              alpha = alpha,
              position = position_use
            ) +
            scalex +
            scaley
        }
        if (plot_type == "rose") {
          p <- ggplot(
            dat,
            aes(x = .data[[g]], y = value, group = .data[[stat.by]])
          ) +
            bg_layer +
            geom_col(
              aes(fill = .data[[stat.by]]),
              width = 0.8,
              color = "black",
              alpha = alpha,
              position = position_use
            ) +
            scalex +
            scaley +
            coord_polar(theta = "x", start = ifelse(flip, pi / 2, 0))
        }
        if (plot_type == "ring" || plot_type == "pie") {
          p <- ggplot(
            dat,
            aes(x = .data[[g]], y = value, group = .data[[stat.by]])
          ) +
            bg_layer +
            geom_col(
              aes(fill = .data[[stat.by]]),
              width = 0.8,
              color = "black",
              alpha = alpha,
              position = position_use
            ) +
            scalex +
            scaley +
            coord_polar(theta = "y", start = ifelse(flip, pi / 2, 0))
        }
        if (plot_type == "area") {
          p <- ggplot(
            dat,
            aes(x = .data[[g]], y = value, group = .data[[stat.by]])
          ) +
            bg_layer +
            geom_area(
              aes(fill = .data[[stat.by]]),
              color = "black",
              alpha = alpha,
              position = position_use
            ) +
            scalex +
            scaley
        }
        if (plot_type == "dot") {
          p <- ggplot(dat, aes(x = .data[[g]], y = .data[[stat.by]])) +
            bg_layer +
            geom_point(
              aes(fill = .data[[stat.by]], size = value),
              color = "black",
              alpha = alpha,
              shape = 21,
              position = position_use
            ) +
            scalex +
            scale_size_area(name = capitalize(stat_type), max_size = 12) +
            guides(size = guide_legend(override.aes = list(fill = "grey30")))
        }
        if (isTRUE(label)) {
          if (plot_type == "dot") {
            p <- p +
              ggrepel::geom_text_repel(
                aes(
                  x = .data[[g]],
                  y = .data[[stat.by]],
                  label = if (stat_type == "count") {
                    value
                  } else {
                    paste0(round(value * 100, 1), "%")
                  },
                ),
                colour = label.fg,
                size = label.size,
                bg.color = label.bg,
                bg.r = label.bg.r,
                point.size = NA,
                max.overlaps = 100,
                min.segment.length = 0,
                force = 0,
                position = position_use
              )
          } else {
            p <- p +
              ggrepel::geom_text_repel(
                aes(
                  label = if (stat_type == "count") {
                    value
                  } else {
                    paste0(round(value * 100, 1), "%")
                  },
                  y = value
                ),
                colour = label.fg,
                size = label.size,
                bg.color = label.bg,
                bg.r = label.bg.r,
                point.size = NA,
                max.overlaps = 100,
                min.segment.length = 0,
                force = 0,
                position = position_use
              )
          }
        }
        if (plot_type %in% c("rose")) {
          axis_text_x <- element_text()
        } else if (plot_type %in% c("ring", "pie")) {
          axis_text_x <- element_text()
        } else {
          axis_text_x <- element_text(
            angle = 45, hjust = 1, vjust = 1
          )
        }
        title <- title %||% sp
        p <- p +
          labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
          scale_fill_manual(
            name = paste0(stat.by, ":"),
            values = colors_use,
            na.value = colors_use["NA"],
            drop = FALSE,
            limits = names(colors_use),
            na.translate = TRUE
          ) +
          do.call(theme_use, theme_args) +
          theme(
            aspect.ratio = aspect.ratio,
            axis.text.x = axis_text_x,
            legend.position = legend.position,
            legend.direction = legend.direction,
            panel.grid.major = if (plot_type == "trend" & stat_type == "percent") {
              element_blank()
            } else {
              element_line(colour = "grey80", linetype = 2)
            }
          ) +
          guides(
            fill = guide_legend(
              title.hjust = 0,
              order = 1,
              override.aes = list(size = 4, color = "black", alpha = 1)
            )
          )
        if (isTRUE(flip) && !plot_type %in% c("pie", "rose")) {
          p <- p + coord_flip()
        }
        return(p)
      }
    )
  } else if (plot_type %in% c("chord", "sankey", "venn", "upset")) {
    colors <- palette_scop(stat.by, palette = palette, palcolor = palcolor)
    if (plot_type == "chord" && isTRUE(combine)) {
      temp <- tempfile(fileext = "png")
      grDevices::png(temp)
      grDevices::dev.control("enable")
      nlev <- nlevels(dat_all[[split.by]])
      if (is.null(nrow) && is.null(ncol)) {
        nrow <- ceiling(sqrt(nlev))
        ncol <- ceiling(nlev / nrow)
      }
      if (is.null(nrow)) {
        nrow <- ceiling(sqrt(ncol))
      }
      if (is.null(ncol)) {
        ncol <- ceiling(sqrt(nrow))
      }
      graphics::par(mfrow = c(nrow, ncol))
    }

    for (sp in levels(dat_all[[split.by]])) {
      dat_use <- dat_split[[ifelse(split.by == "All.groups", 1, sp)]]
      if (plot_type == "venn") {
        check_r("ggVennDiagram")
        dat_list <- as.list(dat_use[, stat.by])
        dat_list <- lapply(
          stats::setNames(
            names(dat_list), names(dat_list)
          ),
          function(x) {
            lg <- dat_list[[x]]
            names(lg) <- rownames(dat_use)
            cellkeep <- names(lg)[lg]
            return(cellkeep)
          }
        )
        venn <- ggVennDiagram::Venn(dat_list)
        data <- ggVennDiagram::process_data(venn)
        dat_venn_region <- ggVennDiagram::venn_region(data)
        idname <- dat_venn_region[["name"]][
          dat_venn_region[["name"]] %in% stat.by
        ]
        names(idname) <- dat_venn_region[["id"]][
          dat_venn_region[["name"]] %in% stat.by
        ]
        idcomb <- strsplit(dat_venn_region[["id"]], split = "")
        colorcomb <- lapply(idcomb, function(x) colors[idname[as.character(x)]])
        dat_venn_region[["colors"]] <- sapply(
          colorcomb,
          function(x) blendcolors(x, mode = "blend")
        )
        dat_venn_region[["label"]] <- paste0(
          dat_venn_region[["count"]],
          "\n",
          round(
            dat_venn_region[["count"]] / sum(dat_venn_region[["count"]]) * 100,
            1
          ),
          "%"
        )
        dat_venn_setedge <- ggVennDiagram::venn_setedge(data)
        dat_venn_setedge[["colors"]] <- colors[stat.by[as.numeric(
          dat_venn_setedge[["id"]]
        )]]

        venn_regionedge_data <- ggVennDiagram::venn_regionedge(data)
        venn_regionedge_data[["colors"]] <- dat_venn_region[["colors"]][match(
          venn_regionedge_data[["id"]], dat_venn_region[["id"]]
        )]

        p <- ggplot() +
          geom_polygon(
            data = venn_regionedge_data,
            aes(X, Y, fill = colors, group = id),
            alpha = alpha
          ) +
          geom_path(
            data = dat_venn_setedge,
            aes(X, Y, group = id),
            color = "black",
            linewidth = 1,
            show.legend = FALSE
          ) +
          ggrepel::geom_text_repel(
            data = ggVennDiagram::venn_setlabel(data),
            aes(X, Y, label = paste0(
              name, "\n(", count, ")"
            )),
            fontface = "bold",
            colour = label.fg,
            size = label.size + 0.5,
            bg.color = label.bg,
            bg.r = label.bg.r,
            point.size = NA,
            max.overlaps = 100,
            force = 0,
            min.segment.length = 0,
            segment.colour = "black"
          ) +
          ggrepel::geom_text_repel(
            data = ggVennDiagram::venn_regionlabel(data),
            aes(X, Y, label = count),
            colour = label.fg,
            size = label.size,
            bg.color = label.bg,
            bg.r = label.bg.r,
            point.size = NA,
            max.overlaps = 100,
            force = 0,
            min.segment.length = 0,
            segment.colour = "black"
          ) +
          scale_fill_identity() +
          coord_equal() +
          theme(
            plot.title = element_text(hjust = 0.5),
            plot.background = element_blank(),
            panel.background = element_blank(),
            axis.title.y = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank()
          )
        p <- p + labs(x = sp, title = title, subtitle = subtitle)
      }

      if (plot_type == "upset") {
        check_r("ggupset")
        for (n in seq_len(nrow(dat_use))) {
          dat_use[["intersection"]][n] <- list(stat.by[unlist(dat_use[
            n,
            stat.by
          ])])
        }
        dat_use <- dat_use[
          sapply(dat_use[["intersection"]], length) > 0, ,
          drop = FALSE
        ]
        p <- ggplot(dat_use, aes(x = intersection)) +
          geom_bar(
            aes(fill = after_stat(count)),
            color = "black",
            width = 0.5,
            show.legend = FALSE
          ) +
          ggrepel::geom_text_repel(
            aes(label = after_stat(count)),
            stat = "count",
            colour = label.fg,
            size = label.size,
            bg.color = label.bg,
            bg.r = label.bg.r,
            point.size = NA,
            max.overlaps = 100,
            force = 0,
            min.segment.length = 0,
            segment.colour = "black"
          ) +
          labs(
            title = title,
            subtitle = subtitle,
            x = sp,
            y = "Intersection size"
          ) +
          ggupset::scale_x_upset(sets = stat.by, n_intersections = 20) +
          scale_fill_gradientn(
            colors = palette_scop(palette = "material-indigo")
          ) +
          theme_scop(
            aspect.ratio = 0.6,
            panel.grid.major = element_line(colour = "grey80", linetype = 2)
          ) +
          ggupset::theme_combmatrix(
            combmatrix.label.text = element_text(size = 12, color = "black"),
            combmatrix.label.extra_spacing = 6
          )
        p <- p + labs(title = title, subtitle = subtitle)
      }

      if (plot_type == "sankey") {
        colors <- palette_scop(
          c(
            unique(
              unlist(
                lapply(
                  dat_all[, stat.by, drop = FALSE],
                  levels
                )
              )
            ),
            NA
          ),
          palette = palette,
          palcolor = palcolor,
          NA_keep = TRUE,
          NA_color = NA_color
        )

        legend_list <- list()
        for (l in stat.by) {
          df <- data.frame(
            factor(levels(dat_use[[l]]), levels = levels(dat_use[[l]]))
          )
          colnames(df) <- l

          legend_list[[l]] <- get_legend(
            ggplot(data = df) +
              geom_col(
                aes(x = 1, y = 1, fill = .data[[l]]),
                color = "black"
              ) +
              scale_fill_manual(
                values = colors[levels(dat_use[[l]])]
              ) +
              guides(
                fill = guide_legend(
                  title.hjust = 0,
                  title.vjust = 0,
                  order = 1,
                  override.aes = list(size = 4, color = "black", alpha = 1)
                )
              ) +
              theme_scop(
                legend.position = "bottom",
                legend.direction = legend.direction
              )
          )

          if (any(is.na(dat_use[[l]]))) {
            raw_levels <- levels(dat_use[[l]])
            dat_use[[l]] <- as.character(dat_use[[l]])
            dat_use[[l]][is.na(dat_use[[l]])] <- "NA"
            dat_use[[l]] <- factor(dat_use[[l]], levels = c(raw_levels, "NA"))
          }
        }

        if (legend.direction == "vertical") {
          legend <- do.call(cbind, legend_list)
        } else {
          legend <- do.call(rbind, legend_list)
        }

        dat <- suppressWarnings(
          make_long(
            dat_use,
            dplyr::all_of(stat.by)
          )
        )
        dat$node <- factor(dat$node, levels = rev(names(colors)))
        p0 <- ggplot(
          dat,
          aes(
            x = x,
            next_x = next_x,
            node = node,
            next_node = next_node,
            fill = node
          )
        ) +
          geom_sankey(
            color = "black",
            flow.alpha = alpha,
            show.legend = FALSE,
            na.rm = FALSE
          ) +
          scale_fill_manual(values = colors, drop = FALSE) +
          scale_x_discrete(expand = c(0, 0.2)) +
          theme_void() +
          theme(axis.text.x = element_text())
        gtable <- as_grob(p0)
        gtable <- add_grob(
          gtable = gtable,
          grob = legend,
          position = legend.position
        )
        p <- patchwork::wrap_plots(gtable)
      }

      if (plot_type == "chord") {
        colors <- palette_scop(
          c(
            unique(
              unlist(
                lapply(
                  dat_all[, stat.by, drop = FALSE],
                  levels
                )
              )
            ),
            NA
          ),
          palette = palette,
          palcolor = palcolor,
          NA_keep = TRUE,
          NA_color = NA_color
        )
        M <- table(
          dat_use[[stat.by[1]]],
          dat_use[[stat.by[2]]],
          useNA = "ifany"
        )
        m <- matrix(M, ncol = ncol(M), dimnames = dimnames(M))
        colnames(m)[is.na(colnames(m))] <- "NA"
        circlize::chordDiagram(
          m,
          grid.col = colors,
          transparency = 0.2,
          link.lwd = 1,
          link.lty = 1,
          link.border = 1
        )
        circlize::circos.clear()
        p <- grDevices::recordPlot()
      }

      plist[[sp]] <- p
    }
  }
  if (isTRUE(combine) && plot_type == "chord") {
    plot <- grDevices::recordPlot()
    grDevices::dev.off()
    unlink(temp)
    return(plot)
  }
  if (isTRUE(combine) && plot_type != "chord") {
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
