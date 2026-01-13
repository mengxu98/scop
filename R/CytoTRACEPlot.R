#' @title Plot CytoTRACE 2 Results
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams CellDimPlot
#' @param ... Additional arguments to be passed to [CellDimPlot] and [FeatureDimPlot].
#'
#' @return
#' If `combine = TRUE`, returns a `patchwork` object combining all plots.
#' If `combine = FALSE`, returns a named list of ggplot objects:
#' \itemize{
#'   \item \code{Score}: UMAP plot colored by score computed by CytoTRACE2;
#'   \item \code{Potency}: UMAP plot colored by potency category computed by CytoTRACE2;
#'   \item \code{Relative}: UMAP plot colored by relative score computed by CytoTRACE2;
#'   \item \code{Phenotype}: UMAP plot colored by phenotype (if \code{group.by} is provided);
#'   \item \code{Boxplot}: Boxplot of score computed by CytoTRACE2 corresponding to phenotype (if \code{group.by} is provided).
#' }
#'
#' @export
#' @seealso
#' [RunCytoTRACE], [CellDimPlot], [FeatureDimPlot]
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunCytoTRACE(pancreas_sub)
#'
#' CytoTRACEPlot(pancreas_sub, group.by = "CellType")
#'
#' plots <- CytoTRACEPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   combine = FALSE
#' )
#' plots$Boxplot
CytoTRACEPlot <- function(
    srt,
    reduction = NULL,
    group.by = NULL,
    combine = TRUE,
    nrow = NULL,
    ncol = NULL,
    byrow = TRUE,
    pt.size = NULL,
    pt.alpha = 1,
    palette = "Paired",
    palcolor = NULL,
    theme_use = "theme_scop",
    theme_args = list(),
    verbose = TRUE,
    ...) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat}",
      message_type = "error"
    )
  }

  required_cols <- c(
    "CytoTRACE2_Potency",
    "CytoTRACE2_Score",
    "CytoTRACE2_Relative"
  )
  missing_cols <- required_cols[!required_cols %in% colnames(srt@meta.data)]
  if (length(missing_cols) > 0) {
    log_message(
      "Missing required CytoTRACE2 results: {.val {missing_cols}}, please run {.fn RunCytoTRACE} first",
      message_type = "error"
    )
  }

  if (!is.null(group.by)) {
    if (!group.by %in% colnames(srt@meta.data)) {
      log_message(
        "Column {.val {group.by}} not found in {.cls Seurat}",
        message_type = "error"
      )
    }
  }

  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  }

  colors_list <- c(
    "#D70440", "#ED5736", "#F9BD10",
    "#FFF799", "#519673", "#806D9E"
  )
  potency_order <- c(
    "Totipotent", "Pluripotent", "Multipotent",
    "Oligopotent", "Unipotent", "Differentiated"
  )
  names(colors_list) <- potency_order

  potency_levels_raw <- unique(
    as.character(srt@meta.data[["CytoTRACE2_Potency"]])
  )
  potency_levels <- intersect(potency_order, potency_levels_raw)
  srt@meta.data[["CytoTRACE2_Potency"]] <- factor(
    srt@meta.data[["CytoTRACE2_Potency"]],
    levels = potency_levels
  )
  potency_colors <- colors_list[potency_levels]

  plist <- list()
  plist[["Score"]] <- FeatureDimPlot(
    srt,
    features = "CytoTRACE2_Score",
    reduction = reduction,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    combine = FALSE,
    ...
  )[[1]] +
    ggplot2::guides(
      color = ggplot2::guide_colorbar(title = "Score")
    )

  plist[["Potency"]] <- CellDimPlot(
    srt,
    group.by = "CytoTRACE2_Potency",
    reduction = reduction,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    palcolor = potency_colors,
    combine = FALSE,
    ...
  )[[1]] +
    ggplot2::guides(
      color = ggplot2::guide_legend(title = "Potency"),
      fill = ggplot2::guide_legend(title = "Potency")
    )

  plist[["Relative"]] <- FeatureDimPlot(
    srt,
    features = "CytoTRACE2_Relative",
    reduction = reduction,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    combine = FALSE,
    ...
  )[[1]] +
    ggplot2::guides(
      color = ggplot2::guide_colorbar(title = "Relative")
    )

  if (!is.null(group.by)) {
    plist[["Phenotype"]] <- CellDimPlot(
      srt,
      group.by = group.by,
      reduction = reduction,
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      palette = palette,
      palcolor = palcolor,
      combine = FALSE,
      ...
    )[[1]]

    plist[["Boxplot"]] <- potency_boxplot(
      srt = srt,
      group.by = group.by,
      pt.alpha = pt.alpha,
      verbose = verbose
    )
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

potency_boxplot <- function(
    srt,
    group.by,
    pt.alpha = 1,
    theme_use = "theme_scop",
    theme_args = list(
      axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1, vjust = 1
      ),
      aspect.ratio = 1
    ),
    verbose = TRUE) {
  mtd <- srt@meta.data[, c(group.by, "CytoTRACE2_Score"), drop = FALSE]
  colnames(mtd)[1] <- "Phenotype_CT"

  mtd <- mtd[!is.na(mtd[["Phenotype_CT"]]), , drop = FALSE]
  mtd <- mtd[!is.na(mtd[["CytoTRACE2_Score"]]), , drop = FALSE]

  if (nrow(mtd) == 0) {
    log_message(
      "No valid data for boxplot after removing NA values",
      message_type = "warning",
      verbose = verbose
    )
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void()
    )
  }

  medians <- mtd |>
    dplyr::group_by(.data[["Phenotype_CT"]]) |>
    dplyr::summarise(
      Potency = stats::median(
        .data[["CytoTRACE2_Score"]],
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      dplyr::desc(.data[["Potency"]])
    )

  mtd <- mtd |>
    dplyr::inner_join(medians, by = "Phenotype_CT")
  mtd[["Phenotype_CT"]] <- factor(
    mtd[["Phenotype_CT"]],
    levels = medians[["Phenotype_CT"]]
  )

  labels <- c(
    "Differentiated", "Unipotent", "Oligopotent",
    "Multipotent", "Pluripotent", "Totipotent"
  )

  colors <- c(
    "#D70440", "#ED5736", "#F9BD10",
    "#FFF799", "#519673", "#806D9E"
  )
  p <- ggplot2::ggplot(
    mtd,
    ggplot2::aes(
      x = .data[["Phenotype_CT"]],
      y = .data[["CytoTRACE2_Score"]]
    )
  ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = .data[["Potency"]]),
      width = 0.8,
      alpha = pt.alpha,
      outlier.shape = NA
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(fill = .data[["Potency"]]),
      width = 0.05,
      height = 0,
      alpha = pt.alpha,
      shape = 21,
      stroke = 0.1,
      size = 1
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, by = 0.2),
      limits = c(0, 1),
      sec.axis = ggplot2::sec_axis(
        trans = ~.,
        breaks = seq(0, 1, by = 1 / 12),
        labels = c(
          "", "Differentiated", "", "Unipotent", "",
          "Oligopotent", "", "Multipotent", "",
          "Pluripotent", "", "Totipotent", ""
        )
      )
    ) +
    ggplot2::scale_fill_gradientn(
      colors = rev(colors),
      breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
      limits = c(0, 1),
      labels = labels
    ) +
    ggplot2::scale_color_gradientn(
      colors = rev(colors),
      breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
      limits = c(0, 1),
      labels = labels
    ) +
    ggplot2::scale_x_discrete(
      labels = function(x) stringr::str_wrap(x, width = 10)
    ) +
    ggplot2::labs(
      x = "Phenotype",
      y = "Potency score"
    ) +
    do.call(theme_use, theme_args)

  return(p)
}
