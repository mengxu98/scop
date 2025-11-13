#' @title Proportion Test Plot
#'
#' @description
#' Generate proportion test plots based on the results from [RunProportionTest].
#'
#' @md
#' @param srt A Seurat object containing proportion test results.
#' @param comparison A character string specifying which comparison to plot.
#' If NULL, plots all comparisons.
#' @param FDR_threshold FDR value cutoff for significance.
#' @param log2FD_threshold Absolute value of log2FD cutoff for significance.
#' @param order_by Method to order clusters.
#' Options: "name" (alphabetical), "value" (by log2FD value).
#' @param pt.size The size of the points.
#' @param pt.alpha The transparency of the points.
#' @param cols.sig Color for significant points and confidence intervals.
#' @param cols.ns Color for non-significant points and confidence intervals.
#' @param aspect.ratio The aspect ratio of the plot.
#' @param xlab A character string specifying the x-axis label.
#' @param ylab A character string specifying the y-axis label.
#' @param theme_use A character string specifying the theme to use for the plot.
#' @param theme_args A list of theme arguments to pass to the `theme_use` function.
#' @param legend.position Position of the legend.
#' @param legend.direction Direction of the legend.
#' @param legend.title Title of the legend.
#' @param combine Whether to combine the plots for each comparison into a single plot.
#' @param nrow An integer value specifying the number of rows in the combined plot.
#' @param ncol An integer value specifying the number of columns in the combined plot.
#' @param byrow Whether to arrange the plots by row in the combined plot.
#'
#' @seealso [RunProportionTest]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- RunProportionTest(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   split.by = "Phase"
#' )
#'
#' ProportionTestPlot(pancreas_sub)
#'
#' # Plot specific comparisons
#' ProportionTestPlot(
#'   pancreas_sub,
#'   comparison = c("G2M_vs_G1", "G2M_vs_S")
#' )
#'
#' # Plot paired comparisons using list format
#' ProportionTestPlot(
#'   pancreas_sub,
#'   comparison = list(c("G2M", "G1"))
#' )
#'
#' ProportionTestPlot(
#'   pancreas_sub,
#'   cols.sig = "blue",
#'   comparison = list(c("G2M", "G1"))
#' )
ProportionTestPlot <- function(
    srt,
    comparison = NULL,
    FDR_threshold = 0.05,
    log2FD_threshold = log2(1.5),
    order_by = c("value", "name"),
    pt.size = 1,
    pt.alpha = 1,
    cols.sig = "red",
    cols.ns = "grey",
    aspect.ratio = NULL,
    xlab = "Cell Type",
    ylab = "log2 (FD)",
    theme_use = "theme_scop",
    theme_args = list(),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.title = "Significance",
    combine = TRUE,
    nrow = NULL,
    ncol = NULL,
    byrow = TRUE) {
  if (!"ProportionTest" %in% names(srt@tools)) {
    log_message(
      "Cannot find the ProportionTest result. Perform {.fn RunProportionTest} first",
      message_type = "error"
    )
  }

  results <- srt@tools[["ProportionTest"]][["results"]]
  if (is.null(results) || length(results) == 0) {
    log_message(
      "No proportion test results found",
      message_type = "error"
    )
  }

  if (!is.null(comparison)) {
    if (is.list(comparison)) {
      if (all(sapply(comparison, function(x) is.character(x) && length(x) == 2))) {
        comparison_names <- c()
        for (pair in comparison) {
          comparison_names <- c(
            comparison_names,
            paste0(pair[1], "_vs_", pair[2]),
            paste0(pair[2], "_vs_", pair[1])
          )
        }
        comparison <- comparison_names
      } else {
        comparison <- unlist(comparison)
      }
    }

    missing_comparisons <- comparison[!comparison %in% names(results)]
    if (length(missing_comparisons) > 0) {
      log_message(
        "Specified comparisons not found in results: {.val {missing_comparisons}}",
        message_type = "error"
      )
    }
    results <- results[comparison]
  }

  plist <- list()
  for (comp_name in names(results)) {
    plot_data <- results[[comp_name]]

    plot_data$significance <- ifelse(
      plot_data$FDR < FDR_threshold & abs(plot_data$obs_log2FD) > log2FD_threshold,
      paste("FDR <", FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2)),
      "n.s."
    )

    plot_data$significance <- factor(
      plot_data$significance,
      levels = c(
        paste("FDR <", FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2)),
        "n.s."
      )
    )

    order_by <- match.arg(order_by)
    if (order_by == "value") {
      clusters_factor <- factor(plot_data$clusters)
      cluster_means <- tapply(plot_data$obs_log2FD, clusters_factor, mean, na.rm = TRUE)
      cluster_order <- names(sort(cluster_means, decreasing = TRUE))
      plot_data$clusters <- factor(plot_data$clusters, levels = cluster_order)
    } else if (order_by == "name") {
      plot_data$clusters <- factor(
        plot_data$clusters,
        levels = sort(unique(plot_data$clusters))
      )
    }

    p <- ggplot(
      plot_data, aes(x = clusters, y = obs_log2FD)
    ) +
      geom_pointrange(
        aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = significance),
        size = pt.size,
        alpha = pt.alpha
      ) +
      geom_hline(
        yintercept = log2FD_threshold, lty = 2, color = "grey"
      ) +
      geom_hline(
        yintercept = -log2FD_threshold, lty = 2, color = "grey"
      ) +
      geom_hline(
        yintercept = 0, color = "black"
      ) +
      scale_color_manual(
        name = legend.title,
        labels = c(
          paste(
            "FDR <", FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2)
          ),
          "n.s."
        ),
        values = c(cols.sig, cols.ns)
      ) +
      labs(
        title = paste0(plot_data$group2[1], " vs ", plot_data$group1[1]),
        x = xlab,
        y = ylab
      ) +
      coord_flip() +
      do.call(theme_use, theme_args) +
      theme(
        aspect.ratio = aspect.ratio,
        legend.position = legend.position,
        legend.direction = if (legend.position %in% c("top", "bottom")) {
          "horizontal"
        } else {
          "vertical"
        }
      )

    plist[[comp_name]] <- p
  }

  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot_names <- names(plist)

      paired_groups <- list()
      unpaired_plots <- c()

      for (i in seq_along(plot_names)) {
        if (i %in% unlist(paired_groups)) next

        current_name <- plot_names[i]
        groups <- strsplit(current_name, "_vs_")[[1]]
        if (length(groups) == 2) {
          reverse_name <- paste0(groups[2], "_vs_", groups[1])
          if (reverse_name %in% plot_names) {
            reverse_idx <- which(plot_names == reverse_name)
            paired_groups[[length(paired_groups) + 1]] <- c(i, reverse_idx)
          } else {
            unpaired_plots <- c(unpaired_plots, i)
          }
        } else {
          unpaired_plots <- c(unpaired_plots, i)
        }
      }

      reordered_plots <- list()
      reordered_names <- c()

      for (pair in paired_groups) {
        reordered_plots <- c(reordered_plots, plist[pair])
        reordered_names <- c(reordered_names, plot_names[pair])
      }

      if (length(unpaired_plots) > 0) {
        reordered_plots <- c(reordered_plots, plist[unpaired_plots])
        reordered_names <- c(reordered_names, plot_names[unpaired_plots])
      }

      if (is.null(nrow) && is.null(ncol)) {
        n_plots <- length(reordered_plots)
        if (length(paired_groups) > 0) {
          ncol <- 2
          nrow <- ceiling(n_plots / ncol)
        } else {
          ncol <- ceiling(sqrt(n_plots))
          nrow <- ceiling(n_plots / ncol)
        }
      }

      plot <- patchwork::wrap_plots(
        plotlist = reordered_plots,
        nrow = nrow,
        ncol = ncol,
        byrow = byrow
      ) +
        patchwork::plot_layout(guides = "collect")

      plot <- plot +
        patchwork::plot_annotation(
          theme = theme(plot.margin = margin(10, 10, 10, 10))
        ) +
        patchwork::plot_annotation(
          theme = theme(
            legend.position = legend.position,
            legend.direction = if (legend.position %in% c("top", "bottom")) {
              "horizontal"
            } else {
              "vertical"
            }
          )
        )
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
}
