#' @title The default theme for scop plot function.
#'
#' @md
#' @param aspect.ratio Aspect ratio of the panel.
#' @param base_size Base font size
#' @param ... Arguments passed to the [ggplot2::theme].
#'
#' @export
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(x = wt, y = mpg, colour = factor(cyl))) +
#'   geom_point()
#' p + theme_scop()
theme_scop <- function(
    aspect.ratio = NULL,
    base_size = 12,
    ...) {
  text_size_scale <- base_size / 12
  args1 <- list(
    aspect.ratio = aspect.ratio,
    text = element_text(
      size = 12 * text_size_scale,
      color = "black"
    ),
    plot.title = element_text(
      size = 14 * text_size_scale,
      colour = "black", vjust = 1
    ),
    plot.subtitle = element_text(
      size = 13 * text_size_scale,
      hjust = 0,
      margin = margin(b = 3)
    ),
    plot.background = element_rect(
      fill = "white",
      color = "white"
    ),
    axis.line = element_blank(),
    axis.title = element_text(
      size = 13 * text_size_scale,
      colour = "black"
    ),
    axis.text = element_text(
      size = 12 * text_size_scale,
      colour = "black"
    ),
    strip.text = element_text(
      size = 12.5 * text_size_scale,
      colour = "black",
      hjust = 0.5,
      margin = margin(3, 3, 3, 3)
    ),
    strip.background = element_rect(
      fill = "transparent", linetype = 0
    ),
    strip.switch.pad.grid = grid::unit(-1, "pt"),
    strip.switch.pad.wrap = grid::unit(-1, "pt"),
    strip.placement = "outside",
    legend.title = element_text(
      size = 12 * text_size_scale,
      colour = "black",
      hjust = 0
    ),
    legend.text = element_text(
      size = 11 * text_size_scale,
      colour = "black"
    ),
    legend.key = element_rect(
      fill = "transparent",
      color = "transparent"
    ),
    legend.key.size = grid::unit(10, "pt"),
    legend.background = element_blank(),
    panel.background = element_rect(
      fill = "white",
      color = "white"
    ),
    panel.border = element_rect(
      fill = "transparent",
      colour = "black",
      linewidth = 1
    )
  )
  args2 <- as.list(match.call())[-1]
  call_envir <- parent.frame(1)
  args2 <- lapply(
    args2, function(arg) {
      if (is.symbol(arg)) {
        eval(arg, envir = call_envir)
      } else if (is.call(arg)) {
        eval(arg, envir = call_envir)
      } else {
        arg
      }
    }
  )
  for (n in names(args2)) {
    args1[[n]] <- args2[[n]]
  }
  args <- args1[names(args1) %in% methods::formalArgs(theme)]
  out <- do.call(
    what = theme,
    args = args
  )
  return(out)
}

#' @title Blank theme
#'
#' @description This function creates a theme with all elements blank except for axis lines and labels.
#' It can optionally add coordinate axes in the plot.
#'
#' @md
#' @param add_coord Whether to add coordinate arrows. Default is `TRUE`.
#' @param xlen_npc The length of the x-axis arrow in "npc".
#' @param ylen_npc The length of the y-axis arrow in "npc".
#' @param xlab The label of the x-axis.
#' @param ylab The label of the y-axis.
#' @param lab_size The size of the axis labels.
#' @param ... Arguments passed to the [ggplot2::theme].
#'
#' @export
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(x = wt, y = mpg, colour = factor(cyl))) +
#'   geom_point()
#' p + theme_blank()
#' p + theme_blank(xlab = "x-axis", ylab = "y-axis", lab_size = 16)
theme_blank <- function(
    add_coord = TRUE,
    xlen_npc = 0.15,
    ylen_npc = 0.15,
    xlab = "",
    ylab = "",
    lab_size = 12, ...) {
  args1 <- list(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.background = element_blank(),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.size = grid::unit(10, "pt"),
    plot.margin = margin(
      lab_size + 2,
      lab_size + 2,
      lab_size + 2,
      lab_size + 2,
      unit = "points"
    )
  )
  args2 <- as.list(match.call())[-1]
  call_envir <- parent.frame(1)
  args2 <- lapply(
    args2, function(arg) {
      if (is.symbol(arg)) {
        eval(arg, envir = call_envir)
      } else if (is.call(arg)) {
        eval(arg, envir = call_envir)
      } else {
        arg
      }
    }
  )
  for (n in names(args2)) {
    args1[[n]] <- args2[[n]]
  }
  args <- args1[names(args1) %in% methods::formalArgs(theme)]
  out <- do.call(
    what = theme,
    args = args
  )
  if (isTRUE(add_coord)) {
    g <- grid::grobTree(
      grid::gList(
        grid::linesGrob(
          x = grid::unit(c(0, xlen_npc), "npc"),
          y = grid::unit(c(0, 0), "npc"),
          arrow = grid::arrow(
            length = grid::unit(0.02, "npc")
          ),
          gp = grid::gpar(lwd = 2)
        ),
        grid::textGrob(
          label = xlab,
          x = grid::unit(0, "npc"),
          y = grid::unit(0, "npc"),
          vjust = 4 / 3,
          hjust = 0,
          gp = grid::gpar(fontsize = lab_size)
        ),
        grid::linesGrob(
          x = grid::unit(c(0, 0), "npc"),
          y = grid::unit(c(0, ylen_npc), "npc"),
          arrow = grid::arrow(length = grid::unit(0.02, "npc")),
          gp = grid::gpar(lwd = 2)
        ),
        grid::textGrob(
          label = ylab, x = grid::unit(0, "npc"),
          y = grid::unit(0, "npc"),
          vjust = -2 / 3,
          hjust = 0,
          rot = 90,
          gp = grid::gpar(fontsize = lab_size)
        )
      )
    )
    return(list(
      list(ggplot2::annotation_custom(g)),
      list(theme_scop() + out),
      list(ggplot2::coord_cartesian(clip = "off"))
    ))
  } else {
    return(list(
      list(theme_scop() + out)
    ))
  }
}

#' @title Color palettes collected
#'
#' @description This function creates a color palette for a given vector of values.
#'
#' @md
#' @param x A vector of character/factor or numeric values.
#' If missing, numeric values 1:n will be used as x.
#' @param n The number of colors to return for numeric values.
#' @param palette Palette name. All available palette names can be queried with [show_palettes].
#' @param palcolor Custom colors used to create a color palette.
#' @param type Type of `x`.
#' Can be one of `"auto"`, `"discrete"` or `"continuous"`.
#' The default is `"auto"`, which automatically detects if `x` is a numeric value.
#' @param matched Whether to return a color vector of the same length as `x`.
#' Default is `FALSE`.
#' @param reverse Whether to invert the colors.
#' Default is `FALSE`.
#' @param NA_keep Whether to keep the color assignment to NA in `x`.
#' Default is `FALSE`.
#' @param NA_color Color assigned to NA if `NA_keep` is `TRUE`.
#' Default is `"grey80"`.
#'
#' @seealso [show_palettes], [palette_list]
#'
#' @export
#'
#' @examples
#' x <- c(1:3, NA, 3:5)
#' (pal1 <- palette_colors(
#'   x,
#'   palette = "Spectral"
#' ))
#' (pal2 <- palette_colors(
#'   x,
#'   palcolor = c("red", "white", "blue")
#' ))
#' (pal3 <- palette_colors(
#'   x,
#'   palette = "Spectral",
#'   n = 10
#' ))
#' (pal4 <- palette_colors(
#'   x,
#'   palette = "Spectral",
#'   n = 10,
#'   reverse = TRUE
#' ))
#' (pal5 <- palette_colors(
#'   x,
#'   palette = "Spectral",
#'   matched = TRUE
#' ))
#' (pal6 <- palette_colors(
#'   x,
#'   palette = "Spectral",
#'   matched = TRUE,
#'   NA_keep = TRUE
#' ))
#' show_palettes(
#'   list(pal1, pal2, pal3, pal4, pal5, pal6)
#' )
#'
#' all_palettes <- show_palettes(return_palettes = TRUE)
#' names(all_palettes)
palette_colors <- function(
    x,
    n = 100,
    palette = "Paired",
    palcolor = NULL,
    type = "auto",
    matched = FALSE,
    reverse = FALSE,
    NA_keep = FALSE,
    NA_color = "grey80") {
  palette_list <- scop::palette_list
  if (missing(x)) {
    x <- 1:n
    type <- "continuous"
  }
  if (!palette %in% names(palette_list)) {
    log_message(
      "The palette ({.val {palette}}) is invalid! Check the available palette names with {.fn show_palettes}. Or pass palette colors via the {.arg palcolor} parameter",
      message_type = "error"
    )
  }
  if (is.list(palcolor)) {
    palcolor <- unlist(palcolor)
  }
  if (all(palcolor == "")) {
    palcolor <- palette_list[[palette]]
  }
  if (is.null(palcolor) || length(palcolor) == 0) {
    palcolor <- palette_list[[palette]]
  }
  if (!is.null(names(palcolor))) {
    if (all(x %in% names(palcolor))) {
      palcolor <- palcolor[intersect(names(palcolor), x)]
    }
  }
  pal_n <- length(palcolor)

  if (!type %in% c("auto", "discrete", "continuous")) {
    log_message(
      "'type' must be one of 'auto','discrete' and 'continuous'.",
      message_type = "error"
    )
  }
  if (type == "auto") {
    if (is.numeric(x)) {
      type <- "continuous"
    } else {
      type <- "discrete"
    }
  }

  if (type == "discrete") {
    if (!is.factor(x)) {
      x <- factor(x, levels = unique(x))
    }
    n_x <- nlevels(x)
    if (isTRUE(attr(palcolor, "type") == "continuous")) {
      color <- grDevices::colorRampPalette(palcolor)(n_x)
    } else {
      color <- ifelse(rep(n_x, n_x) <= pal_n,
        palcolor[1:n_x],
        grDevices::colorRampPalette(palcolor)(n_x)
      )
    }
    names(color) <- levels(x)
    if (any(is.na(x))) {
      color <- c(color, stats::setNames(NA_color, "NA"))
    }
    if (isTRUE(matched)) {
      color <- color[x]
      color[is.na(color)] <- NA_color
    }
  } else if (type == "continuous") {
    if (!is.numeric(x) && all(!is.na(x))) {
      log_message(
        "'x' must be type of numeric when use continuous color palettes.",
        message_type = "error"
      )
    }
    if (all(is.na(x))) {
      values <- as.factor(rep(0, n))
    } else if (length(unique(stats::na.omit(as.numeric(x)))) == 1) {
      values <- as.factor(rep(unique(stats::na.omit(as.numeric(x))), n))
    } else {
      if (isTRUE(matched)) {
        values <- cut(
          x,
          breaks = seq(min(x, na.rm = TRUE),
            max(x, na.rm = TRUE),
            length.out = n + 1
          ),
          include.lowest = TRUE
        )
      } else {
        values <- cut(
          1:100,
          breaks = seq(min(x, na.rm = TRUE),
            max(x, na.rm = TRUE),
            length.out = n + 1
          ),
          include.lowest = TRUE
        )
      }
    }

    n_x <- nlevels(values)
    color <- ifelse(rep(n_x, n_x) <= pal_n,
      palcolor[1:n_x],
      grDevices::colorRampPalette(palcolor)(n_x)
    )
    names(color) <- levels(values)
    if (any(is.na(x))) {
      color <- c(color, stats::setNames(NA_color, "NA"))
    }
    if (isTRUE(matched)) {
      if (all(is.na(x))) {
        color <- NA_color
      } else if (length(unique(stats::na.omit(x))) == 1) {
        color <- color[as.character(unique(stats::na.omit(x)))]
        color[is.na(color)] <- NA_color
      } else {
        color <- color[as.character(values)]
        color[is.na(color)] <- NA_color
      }
    }
  }

  if (isTRUE(reverse)) {
    color <- rev(color)
  }
  if (isFALSE(NA_keep)) {
    color <- color[names(color) != "NA"]
  }
  return(color)
}

#' @title Show the color palettes
#'
#' @description This function displays color palettes using ggplot2.
#'
#' @md
#' @param palettes A list of color palettes.
#' If `NULL`, uses default palettes.
#' @param type A character vector specifying the type of palettes to include.
#' Default is "discrete".
#' @param index A numeric vector specifying the indices of the palettes to include.
#' Default is `NULL`.
#' @param palette_names A character vector specifying the names of the scop palettes to include.
#' Default is `NULL`.
#' @param return_names A logical value indicating whether to return the names of the selected palettes.
#' Default is `TRUE`.
#' @param return_palettes A logical value indicating whether to return the colors of selected palettes.
#' Default is `FALSE`.
#'
#' @seealso [palette_colors], [palette_list]
#'
#' @export
#'
#' @examples
#' show_palettes(
#'   palettes = list(
#'     c("red", "blue", "green"),
#'     c("yellow", "purple", "orange")
#'   )
#' )
#' all_palettes <- show_palettes(return_palettes = TRUE)
#' names(all_palettes)
#' all_palettes[["simspec"]]
#' show_palettes(index = 1:10)
#' show_palettes(
#'   type = "discrete",
#'   index = 1:10
#' )
#' show_palettes(
#'   type = "continuous",
#'   index = 1:10
#' )
#' show_palettes(
#'   palette_names = c(
#'     "Paired", "nejm", "simspec", "Spectral", "jet"
#'   ),
#'   return_palettes = TRUE
#' )
show_palettes <- function(
    palettes = NULL,
    type = c("discrete", "continuous"),
    index = NULL,
    palette_names = NULL,
    return_names = TRUE,
    return_palettes = FALSE) {
  palette_list <- scop::palette_list
  if (!is.null(palettes)) {
    palette_list <- palettes
  } else {
    palette_list <- palette_list[unlist(lapply(palette_list, function(x) isTRUE(attr(x, "type") %in% type)))]
  }
  index <- index[index %in% seq_along(palette_list)]
  if (!is.null(index)) {
    palette_list <- palette_list[index]
  }
  if (is.null(names(palette_list))) {
    names(palette_list) <- seq_along(palette_list)
  }
  if (is.null(palette_names)) {
    palette_names <- palette_names %||% names(palette_list)
  }
  if (any(!palette_names %in% names(palette_list))) {
    palette_not_found <- palette_names[!palette_names %in% names(palette_list)]
    log_message(
      "Can not find the palettes: {.val {palette_not_found}}",
      message_type = "error"
    )
  }
  palette_list <- palette_list[palette_names]

  df <- data.frame(
    palette = rep(
      names(palette_list), sapply(palette_list, length)
    ), color = unlist(palette_list)
  )
  df[["palette"]] <- factor(df[["palette"]], levels = rev(unique(df[["palette"]])))
  df[["color_order"]] <- factor(seq_len(nrow(df)), levels = seq_len(nrow(df)))
  df[["proportion"]] <- as.numeric(1 / table(df$palette)[df$palette])

  p <- ggplot(
    data = df, aes(y = .data[["palette"]], x = .data[["proportion"]], fill = .data[["color_order"]])
  ) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = df[["color"]]) +
    scale_x_continuous(expand = c(0, 0), trans = "reverse") +
    theme_scop(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      panel.border = element_blank()
    )
  print(p)

  if (isTRUE(return_palettes)) {
    return(palette_list)
  }
  if (isTRUE(return_names)) {
    return(palette_names)
  }
}

#' @title Set the panel width/height of a plot to a fixed value
#'
#' @description
#' The ggplot object, when stored, can only specify the height and width of the entire plot, not the panel.
#' The latter is obviously more important to control the final result of a plot.
#' This function can set the panel width/height of plot to a fixed value and rasterize it.
#'
#' @md
#' @param x A ggplot object, a grob object, or a combined plot made by patchwork or cowplot package.
#' @param panel_index Specify the panel to be fixed.
#' If `NULL`, will fix all panels.
#' @param respect Whether row heights and column widths should respect each other.
#' @param width The desired width of the fixed panels.
#' @param height The desired height of the fixed panels.
#' @param margin The margin to add around each panel, in inches.
#' Default is `1`.
#' @param padding The padding to add around each panel, in inches.
#' Default is `0`.
#' @param units The units in which `height`, `width` and `margin` are given.
#' Can be `mm`, `cm`, `in`, etc. See [grid::unit].
#' @param raster Whether to rasterize the panel.
#' @param dpi Plot resolution.
#' @param return_grob Whether to return a grob object instead of a wrapped `patchwork` object.
#' Default is `FALSE`.
#' @param save `NULL` or the file name used to save the plot.
#' @param bg_color The background color of the plot.
#' @param verbose Whether to print messages.
#' Default is `FALSE`.
#' @param ... Additional arguments passed to other functions.
#'
#' @export
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(
#'   data = mtcars, aes(x = mpg, y = wt, colour = cyl)
#' ) +
#'   geom_point() +
#'   facet_wrap(~gear, nrow = 2)
#' # fix the size of panel
#' panel_fix(
#'   p,
#'   width = 5,
#'   height = 3,
#'   units = "cm"
#' )
#' # rasterize the panel
#' panel_fix(
#'   p,
#'   width = 5,
#'   height = 3,
#'   units = "cm",
#'   raster = TRUE,
#'   dpi = 90
#' )
#'
#' # `panel_fix` will build and render the plot when the input is a ggplot object.
#' # so after `panel_fix`, the size of the object will be changed.
#' object.size(p)
#' object.size(
#'   panel_fix(
#'     p,
#'     width = 5,
#'     height = 3,
#'     units = "cm"
#'   )
#' )
#'
#' ## save the plot with appropriate size
#' # p_fix <- panel_fix(
#' #   p,
#' #   width = 5,
#' #   height = 3,
#' #   units = "cm"
#' # )
#' # plot_size <- attr(p_fix, "size")
#' # ggsave(
#' #   filename = "p_fix.png",
#' #   plot = p_fix,
#' #   units = plot_size$units,
#' #   width = plot_size$width,
#' #   height = plot_size$height
#' # )
#' ## or save the plot directly
#' # p_fix <- panel_fix(
#' #   p,
#' #   width = 5,
#' #   height = 3,
#' #   units = "cm",
#' #   save = "p_fix.png"
#' # )
#'
#' # fix the panel of the plot combined by `patchwork`
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' p1 <- CellDimPlot(
#'   pancreas_sub,
#'   "Phase",
#'   aspect.ratio = 1
#' )
#' p2 <- FeatureDimPlot(
#'   pancreas_sub,
#'   "Ins1",
#'   aspect.ratio = 0.5
#' )
#' p <- p1 / p2
#' # fix the panel size for each plot,
#' # the width will be calculated automatically based on `aspect.ratio`
#' panel_fix(p, height = 1)
#'
#' # fix the panel of the plot combined by plot_grid
#' if (requireNamespace("cowplot", quietly = TRUE)) {
#'   p1 <- CellDimPlot(
#'     pancreas_sub,
#'     c("Phase", "SubCellType"),
#'     label = TRUE
#'   )
#'   p2 <- FeatureDimPlot(
#'     pancreas_sub,
#'     c("Ins1", "Gcg"),
#'     label = TRUE
#'   )
#'   p <- cowplot::plot_grid(
#'     p1,
#'     p2,
#'     nrow = 2
#'   )
#'   # plot is combined by plot_grid
#'   # fix the size of panel for each plot
#'   panel_fix(p, height = 1)
#'   # rasterize the panel while keeping all labels and text in vector format
#'   panel_fix(p, height = 1, raster = TRUE, dpi = 30)
#' }
#'
#' # fix the panel of the heatmap
#' ht <- GroupHeatmap(pancreas_sub,
#'   features = c(
#'     "Sox9", "Anxa2", "Bicc1", # Ductal
#'     "Neurog3", "Hes6", # EPs
#'     "Fev", "Neurod1", # Pre-endocrine
#'     "Rbp4", "Pyy", # Endocrine
#'     "Ins1", "Gcg", "Sst", "Ghrl"
#'     # Beta, Alpha, Delta, Epsilon
#'   ),
#'   group.by = c("CellType", "SubCellType"),
#'   show_row_names = TRUE
#' )
#' # the size of the heatmap is not fixed and can be resized by zooming the viewport
#' ht$plot
#' # fix the size of the heatmap according the current viewport
#' panel_fix(ht$plot)
#' # rasterize the heatmap body
#' panel_fix(ht$plot, raster = TRUE, dpi = 30)
#' # fix the size of overall heatmap including annotation and legend
#' panel_fix(ht$plot, height = 4, width = 6)
panel_fix <- function(
    x = NULL,
    panel_index = NULL,
    respect = NULL,
    width = NULL,
    height = NULL,
    margin = 1,
    padding = 0,
    units = "in",
    raster = FALSE,
    dpi = 300,
    return_grob = FALSE,
    bg_color = "white",
    save = NULL,
    verbose = FALSE,
    ...) {
  if (!inherits(x, "gtable")) {
    tryCatch(
      {
        gtable <- as_gtable(x)
      },
      error = function(error) {
        log_message(
          error, "\nCannot convert the x to a gtable object",
          message_type = "error"
        )
      }
    )
  } else {
    gtable <- x
  }
  args <- as.list(match.call())[-1]
  depth <- args[["depth"]]
  if (is.null(depth)) {
    depth <- 1
  }

  if (is.null(panel_index)) {
    non_zero <- grep(
      pattern = "zeroGrob",
      vapply(gtable$grobs, as.character, character(1)),
      invert = TRUE
    )
    panel_index <- grep(
      pattern = "panel|full",
      gtable[["layout"]][["name"]]
    )
    panel_index <- intersect(panel_index, non_zero)
  }
  if (length(panel_index) == 0 && length(gtable$grobs) == 1) {
    panel_index <- 1
  }
  add_margin <- TRUE
  for (i in panel_index) {
    geom_index <- grep(
      pattern = "GeomDrawGrob",
      names(gtable$grobs[[i]][["children"]])
    )
    if (length(geom_index) > 0) {
      log_message(
        "panel {.val {i}} is detected as generated by plot_grid",
        verbose = verbose
      )
      for (j in geom_index) {
        subgrob <- gtable$grobs[[i]][["children"]][[j]][["children"]][[1]][["children"]][[1]]

        if (length(subgrob$grobs[[1]][["children"]]) > 0 &&
          all(sapply(subgrob$grobs[[1]][["children"]], function(x) inherits(x, "recordedGrob")))) {
          subgrob <- panel_fix_overall(
            x = subgrob$grobs[[1]][["children"]],
            width = width,
            height = height,
            margin = padding,
            units = units,
            raster = raster,
            dpi = dpi,
            return_grob = TRUE
          )
        } else {
          subgrob <- panel_fix(
            x = subgrob,
            width = width,
            height = height,
            margin = padding,
            units = units,
            raster = raster,
            dpi = dpi,
            return_grob = TRUE,
            verbose = verbose,
            depth = depth + 1
          )
        }
        gtable$grobs[[i]][["children"]][[j]][["children"]][[1]][["children"]][[1]] <- subgrob
      }
      sum_width <- grid::convertWidth(
        sum(subgrob[["widths"]]),
        unitTo = units,
        valueOnly = TRUE
      ) / as.numeric(gtable$grobs[[i]][["children"]][[j]]$vp$width)
      sum_height <- grid::convertHeight(
        sum(subgrob[["heights"]]),
        unitTo = units,
        valueOnly = TRUE
      ) / as.numeric(gtable$grobs[[i]][["children"]][[j]]$vp$height)
      gtable <- panel_fix_overall(
        gtable,
        panel_index = i,
        width = sum_width,
        height = sum_height,
        margin = ifelse(depth == 1, margin, 0),
        units = units,
        raster = FALSE,
        return_grob = TRUE
      )
    } else if (gtable$grobs[[i]]$name == "layout" || inherits(x, "patchwork")) {
      log_message(
        "panel {.val {i}} is detected as generated by patchwork",
        verbose = verbose
      )
      # if (i == panel_index[1] && length(panel_index) > 1 && isTRUE(verbose)) {
      #   log_message("More than 2 panels detected. panel_fix may not work as expected.")
      # }
      subgrob <- gtable$grobs[[i]]
      if (length(subgrob[["children"]]) > 0 &&
        all(sapply(subgrob[["children"]], function(x) inherits(x, "recordedGrob")))) {
        subgrob <- panel_fix_overall(
          subgrob[["children"]],
          width = width,
          height = height,
          margin = 0,
          units = units,
          raster = raster,
          dpi = dpi,
          return_grob = TRUE
        )
      } else {
        subgrob <- panel_fix(
          subgrob,
          width = width,
          height = height,
          margin = 0,
          units = units,
          raster = raster,
          dpi = dpi,
          return_grob = TRUE,
          verbose = verbose,
          depth = depth + 1
        )
      }
      gtable$grobs[[i]] <- subgrob
      layout <- gtable$layout
      layout[["rowranges"]] <- lapply(
        seq_len(nrow(layout)),
        function(n) layout$t[n]:layout$b[n]
      )
      layout[["colranges"]] <- lapply(
        seq_len(nrow(layout)),
        function(n) layout$l[n]:layout$r[n]
      )
      p_row <- c(layout$t[i], layout$b[i])
      p_col <- c(layout$l[i], layout$r[i])
      background_index <- grep(
        pattern = "background", layout$name
      )
      background_index <- background_index[order(layout$z[background_index], decreasing = TRUE)]
      for (bgi in background_index) {
        if (all(p_row %in% layout[["rowranges"]][[bgi]]) && all(p_col %in% layout[["colranges"]][[bgi]])) {
          p_background_index <- bgi
          break
        }
      }
      gtable <- gtable::gtable_add_rows(
        gtable,
        heights = grid::unit(padding, units),
        pos = layout$t[p_background_index] - 1
      )
      gtable <- gtable::gtable_add_rows(
        gtable,
        heights = grid::unit(padding, units),
        pos = layout$b[p_background_index]
      )
      gtable <- gtable::gtable_add_cols(
        gtable,
        widths = grid::unit(padding, units),
        pos = layout$l[p_background_index] - 1
      )
      gtable <- gtable::gtable_add_cols(
        gtable,
        widths = grid::unit(padding, units),
        pos = layout$r[p_background_index]
      )
      sum_width <- grid::convertWidth(
        sum(subgrob[["widths"]]),
        unitTo = units,
        valueOnly = TRUE
      )
      sum_height <- grid::convertHeight(
        sum(subgrob[["heights"]]),
        unitTo = units,
        valueOnly = TRUE
      )

      gtable <- panel_fix_overall(
        gtable,
        panel_index = i,
        width = sum_width,
        height = sum_height,
        margin = ifelse(depth == 1 & add_margin, margin, 0),
        units = units,
        raster = FALSE,
        respect = TRUE,
        return_grob = TRUE
      )
      if (depth == 1 & add_margin) {
        add_margin <- FALSE
      }
    } else {
      gtable <- panel_fix_overall(
        gtable,
        panel_index = i,
        width = width,
        height = height,
        margin = margin,
        units = units,
        raster = raster,
        dpi = dpi,
        return_grob = TRUE
      )
    }
  }

  if (!is.null(respect)) {
    gtable$respect <- respect
  }

  if (isTRUE(return_grob)) {
    return(gtable)
  } else {
    p <- patchwork::wrap_plots(gtable) +
      theme(
        plot.background = element_rect(
          fill = bg_color, color = bg_color
        )
      )
    if (units != "null") {
      plot_width <- grid::convertWidth(
        sum(gtable[["widths"]]),
        unitTo = units,
        valueOnly = TRUE
      )
      plot_height <- grid::convertHeight(
        sum(gtable[["heights"]]),
        unitTo = units,
        valueOnly = TRUE
      )
      attr(p, "size") <- list(
        width = plot_width,
        height = plot_height,
        units = units
      )
    }

    if (!is.null(save) && is.character(save) && nchar(save) > 0) {
      if (units == "null") {
        log_message(
          "{.arg units} can not be 'null' if want to save the plot",
          message_type = "error"
        )
      }
      filename <- normalizePath(save)
      log_message(
        "Save the plot to the file: {.file {filename}}",
        verbose = verbose
      )
      if (!dir.exists(dirname(filename))) {
        dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
      }
      ggplot2::ggsave(
        plot = p,
        filename = filename,
        width = plot_width,
        height = plot_height,
        units = units,
        dpi = dpi,
        limitsize = FALSE
      )
    }
    return(p)
  }
}

#' @rdname panel_fix
#' @export
panel_fix_overall <- function(
    x,
    panel_index = NULL,
    respect = NULL,
    width = NULL,
    height = NULL,
    margin = 1,
    units = "in",
    raster = FALSE,
    dpi = 300,
    return_grob = FALSE,
    bg_color = "white",
    save = NULL,
    verbose = TRUE) {
  if (!inherits(x, "gtable")) {
    if (inherits(x, "gTree")) {
      x <- x[["children"]]
    }
    tryCatch(
      {
        gtable <- as_gtable(x)
      },
      error = function(error) {
        log_message(
          error, "\nCannot convert the x to a gtable object",
          message_type = "error"
        )
      }
    )
  } else {
    gtable <- x
  }

  if (is.null(panel_index)) {
    non_zero <- grep(
      pattern = "zeroGrob",
      vapply(
        gtable$grobs, as.character, character(1)
      ), invert = TRUE
    )
    panel_index <- grep("panel|full", gtable[["layout"]][["name"]])
    panel_index <- intersect(panel_index, non_zero)
  }
  if (length(panel_index) == 0 && length(gtable$grobs) == 1) {
    panel_index <- 1
  }
  if (!length(width) %in% c(0, 1, length(panel_index)) || !length(height) %in% c(0, 1, length(panel_index))) {
    log_message(
      "The length of 'width' and 'height' must be 1 or the length of panels.",
      message_type = "error"
    )
  }

  if (inherits(x, "gList")) {
    panel_index <- 1
    panel_index_h <- panel_index_w <- list(1)
    w_comp <- h_comp <- list(grid::unit(1, "null"))
    w <- h <- list(grid::unit(1, "null"))
  } else if (length(panel_index) > 0) {
    panel_index_w <- panel_index_h <- list()
    w_comp <- h_comp <- list()
    w <- h <- list()
    for (i in seq_along(panel_index)) {
      index <- panel_index[i]
      panel_index_h[[i]] <- sort(
        unique(c(
          gtable[["layout"]][["t"]][index],
          gtable[["layout"]][["b"]][index]
        ))
      )
      panel_index_w[[i]] <- sort(
        unique(c(
          gtable[["layout"]][["l"]][index],
          gtable[["layout"]][["r"]][index]
        ))
      )
      w_comp[[i]] <- gtable[["widths"]][seq(min(panel_index_w[[i]]), max(panel_index_w[[i]]))]
      h_comp[[i]] <- gtable[["heights"]][seq(min(panel_index_h[[i]]), max(panel_index_h[[i]]))]

      if (length(w_comp[[i]]) == 1) {
        w[[i]] <- w_comp[[i]]
      } else if (length(w_comp[[i]]) > 1 && any(grid::unitType(w_comp[[i]]) == "null")) {
        w[[i]] <- grid::unit(1, units = "null")
      } else {
        w[[i]] <- sum(w_comp[[i]])
      }
      if (length(h_comp[[i]]) == 1) {
        h[[i]] <- h_comp[[i]]
      } else if (length(h_comp[[i]]) > 1 && any(grid::unitType(h_comp[[i]]) == "null")) {
        h[[i]] <- grid::unit(1, units = "null")
      } else {
        h[[i]] <- sum(h_comp[[i]])
      }
    }
  } else {
    log_message(
      "No panel detected",
      message_type = "error"
    )
  }

  if (units != "null") {
    raw_w <- sapply(
      w, function(x) {
        grid::convertWidth(x, unitTo = units, valueOnly = TRUE)
      }
    )
    raw_h <- sapply(
      h, function(x) {
        grid::convertHeight(x, unitTo = units, valueOnly = TRUE)
      }
    )
    for (i in seq_along(w)) {
      if (grid::unitType(w[[i]]) == "null" || grid::convertUnit(w[[i]], unitTo = "cm", valueOnly = TRUE) < 1e-10) {
        raw_w[i] <- 0
      }
    }
    for (i in seq_along(h)) {
      if (grid::unitType(h[[i]]) == "null" || grid::convertUnit(h[[i]], unitTo = "cm", valueOnly = TRUE) < 1e-10) {
        raw_h[i] <- 0
      }
    }
    if (isTRUE(gtable$respect)) {
      raw_aspect <- sapply(h, as.vector) / sapply(w, as.vector)
    } else {
      if (all(raw_w != 0) && all(raw_h != 0)) {
        raw_aspect <- raw_h / raw_w
      } else {
        raw_aspect <- grid::convertHeight(
          grid::unit(1, "npc"), "cm",
          valueOnly = TRUE
        ) / grid::convertWidth(grid::unit(1, "npc"), "cm", valueOnly = TRUE)
      }
    }

    if (is.null(width) && is.null(height)) {
      width <- raw_w
      height <- raw_h
      if (all(width == 0) && all(height == 0)) {
        width <- grid::convertWidth(
          grid::unit(1, "npc"), units,
          valueOnly = TRUE
        )
        height <- grid::convertHeight(
          grid::unit(1, "npc"), units,
          valueOnly = TRUE
        )
        if (isTRUE(gtable$respect)) {
          if (raw_aspect <= 1) {
            height <- width * raw_aspect
          } else {
            width <- height / raw_aspect
          }
        }
      }
    }

    for (i in seq_along(raw_aspect)) {
      if (is.finite(raw_aspect[i]) && raw_aspect[i] != 0) {
        if (is.null(width[i]) || is.na(width[i]) || width[i] == 0) {
          width[i] <- height[i] / raw_aspect[i]
        }
        if (is.null(height[i]) || is.na(height[i]) || height[i] == 0) {
          height[i] <- width[i] * raw_aspect[i]
        }
      }
    }

    for (i in seq_along(width)) {
      if (inherits(width[i], "unit")) {
        width[i] <- grid::convertWidth(
          width[i],
          unitTo = units,
          valueOnly = TRUE
        )
      }
    }
    for (i in seq_along(height)) {
      if (inherits(height[i], "unit")) {
        height[i] <- grid::convertHeight(
          height[i],
          unitTo = units,
          valueOnly = TRUE
        )
      }
    }
  }

  if (length(width) == 1) {
    width <- rep(width, length(panel_index))
  }
  if (length(height) == 1) {
    height <- rep(height, length(panel_index))
  }
  for (i in seq_along(panel_index)) {
    if (!is.null(width)) {
      width_unit <- width[i] / length(w_comp[[i]])
      gtable[["widths"]][seq(min(panel_index_w[[i]]), max(panel_index_w[[i]]))] <- rep(grid::unit(width_unit, units = units), length(w_comp[[i]]))
    }
    if (!is.null(height)) {
      height_unit <- height[i] / length(h_comp[[i]])
      gtable[["heights"]][seq(min(panel_index_h[[i]]), max(panel_index_h[[i]]))] <- rep(grid::unit(height_unit, units = units), length(h_comp[[i]]))
    }
  }
  gtable <- gtable::gtable_add_padding(
    gtable,
    padding = grid::unit(margin, units = units)
  )

  if (isTRUE(raster)) {
    check_r(c("png", "ragg"))
    for (i in seq_along(panel_index)) {
      index <- panel_index[i]
      g <- g_new <- gtable$grobs[[index]]
      vp <- g$vp
      children_order <- g$childrenOrder
      if (is.null(g$vp)) {
        g$vp <- grid::viewport()
      }

      for (j in seq_along(g[["children"]])) {
        child <- g[["children"]][[j]]
        child_nm <- names(g[["children"]])[j]
        if (!is.null(child$vp) ||
          any(grepl("(text)|(label)", child_nm)) ||
          any(grepl("(text)|(segments)|(legend)", class(child$list[[1]])))) {
          zero <- ggplot2::zeroGrob()
          zero$name <- g[["children"]][[j]]$name
          g[["children"]][[j]] <- zero
        } else if (inherits(child$list[[1]], "grob") || is.null(child$list[[1]])) {
          g_new[["children"]][[j]] <- ggplot2::zeroGrob()
        }
      }
      temp <- tempfile(fileext = "png")
      ragg::agg_png(
        temp,
        width = width[i],
        height = height[i],
        bg = "transparent",
        res = dpi,
        units = units
      )
      grid::grid.draw(g)
      grDevices::dev.off()
      g_ras <- grid::rasterGrob(png::readPNG(temp, native = TRUE))
      unlink(temp)
      g <- grid::addGrob(g_new, g_ras)
      g$vp <- vp
      g$childrenOrder <- c(g_ras$name, children_order)
      gtable$grobs[[index]] <- g
    }
  }

  if (!is.null(respect)) {
    gtable$respect <- respect
  }

  if (isTRUE(return_grob)) {
    return(gtable)
  } else {
    p <- patchwork::wrap_plots(gtable) +
      theme(plot.background = element_rect(fill = bg_color, color = bg_color))
    if (units != "null") {
      plot_width <- grid::convertWidth(
        sum(gtable[["widths"]]),
        unitTo = units,
        valueOnly = TRUE
      )
      plot_height <- grid::convertHeight(
        sum(gtable[["heights"]]),
        unitTo = units,
        valueOnly = TRUE
      )
      attr(p, "size") <- list(
        width = plot_width,
        height = plot_height,
        units = units
      )
    }

    if (!is.null(save) && is.character(save) && nchar(save) > 0) {
      if (units == "null") {
        log_message(
          "{.arg units} can not be 'null' if want to save the plot",
          message_type = "error"
        )
      }
      filename <- normalizePath(save)
      log_message(
        "Save the plot to the file: {.file {filename}}",
        verbose = verbose
      )
      if (!dir.exists(dirname(filename))) {
        dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
      }
      ggplot2::ggsave(
        plot = p,
        filename = filename,
        width = plot_width,
        height = plot_height,
        units = units,
        dpi = dpi,
        limitsize = FALSE
      )
    }
    return(p)
  }
}

#' @title Drop unused data in the plot
#'
#' @md
#' @param p A `ggplot` object or a `patchwork` object.
#'
#' @export
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(data = mtcars, aes(x = mpg, y = wt, colour = cyl)) +
#'   geom_point() +
#'   scale_x_continuous(limits = c(10, 30)) +
#'   scale_y_continuous(limits = c(1, 6)) +
#'   theme_scop()
#' object.size(p)
#'
#' p_drop <- drop_data(p)
#' object.size(p_drop)
#'
#' p / p_drop
drop_data <- function(p) {
  UseMethod(generic = "drop_data", object = p)
}

#' @export
#' @rdname drop_data
#' @method drop_data ggplot
drop_data.ggplot <- function(p) {
  p <- ggplot2:::plot_clone(p)

  # fix the scales for x/y axis and 'fill', 'color', 'shape',...
  for (i in seq_along(p$scales$scales)) {
    if (inherits(p$scales$scales[[i]], "ScaleDiscrete")) {
      p$scales$scales[[i]][["drop"]] <- FALSE
    }
    if (inherits(p$scales$scales[[i]], "ScaleContinuous")) {
      limits <- p$scales$scales[[i]]$get_limits()
      if (p$scales$scales[[i]]$aesthetics[1] == "x") {
        p$coordinates$limits$x <- limits
      }
      if (p$scales$scales[[i]]$aesthetics[1] == "y") {
        p$coordinates$limits$y <- limits
      }
    }
  }

  vars <- get_vars(p)
  # drop main data
  if (length(p$data) > 0) {
    vars_modified <- names(
      which(
        sapply(
          p$data[, intersect(colnames(p$data), vars), drop = FALSE], class
        ) == "character"
      )
    )
    for (v in vars_modified) {
      p$data[[v]] <- as.factor(p$data[[v]])
    }
    p$data <- p$data[1, , drop = FALSE]
  }

  # drop layer data
  for (i in seq_along(p$layers)) {
    if (length(p$layers[[i]]$data) > 0) {
      vars_modified <- names(
        which(
          sapply(
            p$layers[[i]]$data[, intersect(colnames(p$layers[[i]]$data), vars), drop = FALSE], class
          ) == "character"
        )
      )
      for (v in vars_modified) {
        p$layers[[i]]$data[[v]] <- as.factor(p$layers[[i]]$data[[v]])
      }
      p$layers[[i]]$data <- p$layers[[i]]$data[1, , drop = FALSE]
    }
  }

  return(p)
}

#' @export
#' @rdname drop_data
#' @method drop_data patchwork
drop_data.patchwork <- function(p) {
  for (i in seq_along(p$patches$plots)) {
    p$patches$plots[[i]] <- drop_data(p$patches$plots[[i]])
  }
  drop_data.ggplot(p)
}

#' @export
#' @rdname drop_data
#' @method drop_data default
drop_data.default <- function(p) {
  p
}

#' @title Slim unused data in the plot
#'
#' @md
#' @param p A `ggplot` object or a `patchwork` object.
#' @export
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(data = mtcars, aes(x = mpg, y = wt, colour = cyl)) +
#'   geom_point()
#' object.size(p)
#' colnames(p$data)
#'
#' p_slim <- slim_data(p)
#' object.size(p_slim)
#' colnames(p_slim$data)
slim_data <- function(p) {
  UseMethod(generic = "slim_data", object = p)
}

#' @export
#' @rdname slim_data
#' @method slim_data ggplot
slim_data.ggplot <- function(p) {
  vars <- get_vars(p)
  if (length(vars) > 0) {
    p$data <- p$data[, intersect(colnames(p$data), vars), drop = FALSE]
    for (i in seq_along(p$layers)) {
      if (length(p$layers[[i]]$data) > 0) {
        p$layers[[i]]$data <- p$layers[[i]]$data[, intersect(colnames(p$layers[[i]]$data), vars), drop = FALSE]
      }
    }
  }
  return(p)
}

#' @export
#' @rdname slim_data
#' @method slim_data patchwork
slim_data.patchwork <- function(p) {
  for (i in seq_along(p$patches$plots)) {
    p$patches$plots[[i]] <- slim_data(p$patches$plots[[i]])
  }
  p <- slim_data.ggplot(p)
  return(p)
}

#' @export
#' @method slim_data default
slim_data.default <- function(p) {
  return(p)
}


#' @title Get used vars in a ggplot object
#'
#' @md
#' @param p A `ggplot` object.
#' @param reverse If `TRUE` then will return unused vars.
#' @param verbose Whether to print messages.
#'
#' @export
get_vars <- function(p, reverse, verbose = FALSE) {
  mappings <- c(
    as.character(p$mapping),
    unlist(
      lapply(p$layers, function(x) as.character(x$mapping))
    ),
    unlist(
      lapply(p$layers, function(x) names(p$layers[[1]]$aes_params))
    ),
    names(
      p$facet$params$facets
    ), names(p$facet$params$rows), names(p$facet$params$cols)
  )
  vars <- unique(
    unlist(
      strsplit(
        gsub(
          "[~\\[\\]\\\"\\(\\)]", " ", unique(mappings),
          perl = TRUE
        ), " "
      )
    )
  )
  vars_used <- intersect(
    unique(
      c(
        colnames(p$data), unlist(lapply(p$layers, function(x) colnames(x$data)))
      )
    ), vars
  )

  log_message(
    "vars_used: ", paste0(vars_used, collapse = ","), "\n",
    "vars_notused: ", paste0(setdiff(names(p$data), vars), collapse = ","),
    verbose = verbose
  )
  return(vars_used)
}

#' @title Convert a color with specified alpha level
#'
#' @md
#' @param colors Color vectors.
#' @param alpha Alpha level in `[0,1]`.
#'
#' @export
#'
#' @examples
#' colors <- c("red", "blue", "green")
#' adjcolors(colors, 0.5)
#' ggplot2::alpha(colors, 0.5)
#'
#' show_palettes(
#'   list(
#'     "raw" = colors,
#'     "adjcolors" = adjcolors(colors, 0.5),
#'     "ggplot2::alpha" = ggplot2::alpha(colors, 0.5)
#'   )
#' )
adjcolors <- function(colors, alpha) {
  color_df <- as.data.frame(
    grDevices::col2rgb(colors) / 255
  )
  colors_out <- sapply(
    color_df, function(color) {
      color_rgb <- RGBA2RGB(list(color, alpha))
      grDevices::rgb(color_rgb[1], color_rgb[2], color_rgb[3])
    }
  )
  colors_out
}

#' @title Blends a list of colors using the specified blend mode
#'
#' @md
#' @param colors Color vectors.
#' @param mode Blend mode.
#' One of `"blend"`, `"average"`, `"screen"`, or `"multiply"`.
#'
#' @export
#'
#' @examples
#' blend <- c(
#'   "red",
#'   "green",
#'   blendcolors(c("red", "green"),
#'     mode = "blend"
#'   )
#' )
#' average <- c(
#'   "red",
#'   "green",
#'   blendcolors(c("red", "green"),
#'     mode = "average"
#'   )
#' )
#' screen <- c(
#'   "red",
#'   "green",
#'   blendcolors(c("red", "green"),
#'     mode = "screen"
#'   )
#' )
#' multiply <- c(
#'   "red",
#'   "green",
#'   blendcolors(c("red", "green"),
#'     mode = "multiply"
#'   )
#' )
#' show_palettes(
#'   list(
#'     "blend" = blend,
#'     "average" = average,
#'     "screen" = screen,
#'     "multiply" = multiply
#'   )
#' )
blendcolors <- function(
    colors,
    mode = c("blend", "average", "screen", "multiply")) {
  mode <- match.arg(mode)
  colors <- colors[!is.na(colors)]
  if (length(colors) == 0) {
    return(NA)
  }
  if (length(colors) == 1) {
    return(colors)
  }
  rgb <- as.list(
    as.data.frame(
      grDevices::col2rgb(colors) / 255
    )
  )
  Clist <- lapply(rgb, function(x) {
    list(x, 1)
  })
  blend_color <- BlendRGBList(Clist, mode = mode)
  blend_color <- grDevices::rgb(blend_color[1], blend_color[2], blend_color[3])
  return(blend_color)
}

RGBA2RGB <- function(RGBA, BackGround = c(1, 1, 1)) {
  A <- RGBA[[length(RGBA)]]
  RGB <- RGBA[[-length(RGBA)]] * A + BackGround * (1 - A)
  return(RGB)
}

Blend2Color <- function(C1, C2, mode = "blend") {
  c1 <- C1[[1]]
  c1a <- C1[[2]]
  c2 <- C2[[1]]
  c2a <- C2[[2]]
  A <- 1 - (1 - c1a) * (1 - c2a)
  if (A < 1.0e-6) {
    return(list(c(0, 0, 0), 1))
  }
  if (mode == "blend") {
    out <- (c1 * c1a + c2 * c2a * (1 - c1a)) / A
    A <- 1
  }
  if (mode == "average") {
    out <- (c1 + c2) / 2
    out[out > 1] <- 1
  }
  if (mode == "screen") {
    out <- 1 - (1 - c1) * (1 - c2)
  }
  if (mode == "multiply") {
    out <- c1 * c2
  }
  return(list(out, A))
}

BlendRGBList <- function(
    Clist,
    mode = "blend",
    RGB_BackGround = c(1, 1, 1)) {
  N <- length(Clist)
  ClistUse <- Clist
  while (N != 1) {
    temp <- ClistUse
    ClistUse <- list()
    for (C in temp[1:(length(temp) - 1)]) {
      c1 <- C[[1]]
      a1 <- C[[2]]
      c2 <- temp[[length(temp)]][[1]]
      a2 <- temp[[length(temp)]][[2]]
      ClistUse <- append(
        ClistUse,
        list(
          Blend2Color(
            C1 = list(c1, a1 * (1 - 1 / N)),
            C2 = list(c2, a2 * 1 / N),
            mode = mode
          )
        )
      )
    }
    N <- length(ClistUse)
  }
  Result <- list(ClistUse[[1]][[1]], ClistUse[[1]][[2]])
  Result <- RGBA2RGB(Result, BackGround = RGB_BackGround)
  return(Result)
}

build_patchwork <- function(
    x,
    guides = "auto",
    table_rows = 18,
    table_cols = 15,
    panel_row = 10,
    panel_col = 8) {
  x$layout <- utils::modifyList(
    patchwork:::default_layout,
    x$layout[!vapply(x$layout, is.null, logical(1))]
  )

  guides <- if (guides == "collect" && x$layout$guides != "keep") {
    "collect"
  } else {
    x$layout$guides
  }
  gt <- lapply(
    x$plots,
    patchwork:::plot_table,
    guides = guides
  )
  fixed_asp <- vapply(
    gt, function(x) isTRUE(x$respect), logical(1)
  )
  guide_grobs <- unlist(
    lapply(gt, `[[`, "collected_guides"),
    recursive = FALSE
  )
  gt <- lapply(
    gt,
    patchwork:::simplify_gt
  )
  gt <- patchwork:::add_insets(gt)
  if (is.null(x$layout$design)) {
    if (is.null(x$layout$ncol) && !is.null(x$layout$widths) && length(x$layout$widths) > 1) {
      x$layout$ncol <- length(x$layout$widths)
    }
    if (is.null(x$layout$nrow) && !is.null(x$layout$heights) && length(x$layout$heights) > 1) {
      x$layout$nrow <- length(x$layout$heights)
    }
    dims <- ggplot2::wrap_dims(
      length(gt),
      nrow = x$layout$nrow,
      ncol = x$layout$ncol
    )
    x$layout$design <- patchwork:::create_design(
      dims[2],
      dims[1],
      x$layout$byrow
    )
  } else {
    dims <- c(
      max(x$layout$design$b),
      max(x$layout$design$r)
    )
  }

  gt_new <- gtable::gtable(
    grid::unit(rep(0, table_cols * dims[2]), "null"),
    grid::unit(rep(0, table_rows * dims[1]), "null")
  )
  design <- as.data.frame(unclass(x$layout$design))
  if (nrow(design) < length(gt)) {
    log_message(
      "Too few patch areas to hold all plots. Dropping plots",
      message_type = "warning"
    )
    gt <- gt[seq_len(nrow(design))]
    fixed_asp <- fixed_asp[seq_len(nrow(design))]
  } else {
    design <- design[seq_along(gt), ]
  }
  if (any(design$t < 1)) design$t[design$t < 1] <- 1
  if (any(design$l < 1)) design$l[design$l < 1] <- 1
  if (any(design$b > dims[1])) design$b[design$b > dims[1]] <- dims[1]
  if (any(design$r > dims[2])) design$r[design$r > dims[2]] <- dims[2]
  max_z <- lapply(gt, function(x) max(x$layout$z))
  max_z <- c(0, cumsum(max_z))
  gt_new$layout <- do.call(
    rbind,
    lapply(
      seq_along(gt), function(i) {
        loc <- design[i, ]
        lay <- gt[[i]]$layout
        lay$name <- paste0(lay$name, "-", i)
        lay$t <- lay$t +
          ifelse(
            lay$t <= panel_row, (loc$t - 1) * table_rows,
            (loc$b - 1) * table_rows
          )
        lay$l <- lay$l +
          ifelse(
            lay$l <= panel_col,
            (loc$l - 1) * table_cols,
            (loc$r - 1) * table_cols
          )
        lay$b <- lay$b +
          ifelse(lay$b < panel_row,
            (loc$t - 1) * table_rows,
            (loc$b - 1) * table_rows
          )
        lay$r <- lay$r +
          ifelse(lay$r < panel_col,
            (loc$l - 1) * table_cols,
            (loc$r - 1) * table_cols
          )
        lay$z <- lay$z + max_z[i]
        lay
      }
    )
  )
  table_dimensions <- patchwork:::table_dims(
    lapply(gt, `[[`, "widths"),
    lapply(gt, `[[`, "heights"),
    design,
    dims[2],
    dims[1]
  )
  gt_new$grobs <- patchwork:::set_grob_sizes(
    gt,
    table_dimensions$widths,
    table_dimensions$heights, design
  )
  gt_new$widths <- table_dimensions$widths
  gt_new$heights <- table_dimensions$heights
  widths <- rep(x$layout$widths, length.out = dims[2])
  heights <- rep(x$layout$heights, length.out = dims[1])
  gt_new <- patchwork:::set_panel_dimensions(
    gt_new,
    gt,
    widths,
    heights,
    fixed_asp,
    design
  )
  if (x$layout$guides == "collect") {
    guide_grobs <- patchwork:::collapse_guides(guide_grobs)
    if (length(guide_grobs) != 0) {
      theme <- x$annotation$theme
      if (!attr(theme, "complete")) {
        theme <- ggplot2::theme_get() + theme
      }
      guide_grobs <- patchwork:::assemble_guides(guide_grobs, theme)
      gt_new <- patchwork:::attach_guides(gt_new, guide_grobs, theme)
    }
  } else {
    gt_new$collected_guides <- guide_grobs
  }

  class(gt_new) <- c("gtable_patchwork", class(gt_new))
  gt_new
}

patchwork_grob <- function(x, ...) {
  annotation <- utils::modifyList(
    patchwork:::default_annotation,
    x$patches$annotation[!vapply(x$patches$annotation, is.null, logical(1))]
  )
  x <- patchwork:::recurse_tags(
    x,
    annotation$tag_levels,
    annotation$tag_prefix,
    annotation$tag_suffix,
    annotation$tag_sep
  )$patches
  plot <- patchwork:::get_patches(x)
  gtable <- build_patchwork(plot)
  gtable <- patchwork:::annotate_table(gtable, annotation)
  class(gtable) <- setdiff(class(gtable), "gtable_patchwork")
  gtable
}

as_grob <- function(plot, ...) {
  if (inherits(plot, "gList")) {
    grid::grobTree(plot)
  } else if (inherits(plot, "patchwork")) {
    patchwork_grob(plot, ...)
  } else if (inherits(plot, "ggplot")) {
    ggplot2::ggplotGrob(plot)
  } else {
    cli::cli_alert_warning(
      "Cannot convert object of {.cls {class(plot)}} into a grob"
    )
  }
}

as_gtable <- function(plot, ...) {
  if (inherits(plot, "gtable")) {
    return(plot)
  }
  if (inherits(plot, "grob")) {
    u <- grid::unit(1, "null")
    gt <- gtable::gtable_col(NULL, list(plot), u, u)
    gt$layout$clip <- "inherit"
    return(gt)
  } else {
    grob <- as_grob(plot, ...)
    if (inherits(grob, "gtable")) {
      return(grob)
    } else {
      return(cowplot::as_gtable(grob, ...))
    }
  }
}

get_legend <- function(plot) {
  plot <- as_gtable(plot)
  grob_names <- plot$layout$name
  grobs <- plot$grobs
  grob_index <- which(
    grepl(
      "guide-box-bottom",
      grob_names
    )
  )
  grob_index <- grob_index[1]
  matched_grobs <- grobs[[grob_index]]
  matched_grobs
}

add_grob <- function(
    gtable,
    grob,
    position = c("top", "bottom", "left", "right", "none"),
    space = NULL,
    clip = "on") {
  position <- match.arg(position)
  if (position == "none" || is.null(grob)) {
    return(gtable)
  }

  if (is.null(space)) {
    if (gtable::is.gtable(grob)) {
      if (position %in% c("top", "bottom")) {
        space <- sum(grob$heights)
      } else {
        space <- sum(grob$widths)
      }
    } else if (grid::is.grob(grob)) {
      if (position %in% c("top", "bottom")) {
        space <- grid::grobHeight(grob)
      } else {
        space <- grid::grobWidth(grob)
      }
    }
  }

  if (position == "top") {
    gtable <- gtable::gtable_add_rows(
      gtable,
      space,
      0
    )
    gtable <- gtable::gtable_add_grob(
      gtable, grob,
      t = 1,
      l = mean(
        gtable$layout[grepl(pattern = "panel", x = gtable$layout$name), "l"]
      ),
      clip = clip
    )
  }
  if (position == "bottom") {
    gtable <- gtable::gtable_add_rows(
      gtable,
      space,
      -1
    )
    gtable <- gtable::gtable_add_grob(
      gtable, grob,
      t = dim(gtable)[1],
      l = mean(
        gtable$layout[grepl(pattern = "panel", x = gtable$layout$name), "l"]
      ), clip = clip
    )
  }
  if (position == "left") {
    gtable <- gtable::gtable_add_cols(
      gtable,
      space,
      0
    )
    gtable <- gtable::gtable_add_grob(
      gtable,
      grob,
      t = mean(
        gtable$layout[grep("panel", gtable$layout$name), "t"]
      ),
      l = 1,
      clip = clip
    )
  }
  if (position == "right") {
    gtable <- gtable::gtable_add_cols(
      gtable,
      space,
      -1
    )
    gtable <- gtable::gtable_add_grob(
      gtable,
      grob,
      t = mean(
        gtable$layout[grep("panel", gtable$layout$name), "t"]
      ),
      l = dim(gtable)[2],
      clip = clip
    )
  }
  return(gtable)
}
