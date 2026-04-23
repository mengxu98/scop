% Generated manually because roxygen2 is unavailable in this environment
% Please edit documentation in R/ProportionTestPlot.R
\name{ProportionTestPlot}
\alias{ProportionTestPlot}
\title{Proportion Test Plot}
\usage{
ProportionTestPlot(
  srt,
  comparison = NULL,
  proportion_method = NULL,
  result_level = c("group", "neighborhood"),
  plot_type = c(
    "effect",
    "umap",
    "volcano",
    "manhattan",
    "ring",
    "milo_beeswarm",
    "milo_graph",
    "sccoda_forest"
  ),
  umap_mode = c("discrete", "continuous"),
  reduction = "UMAP",
  projection_args = list(),
  FDR_threshold = 0.05,
  log2FD_threshold = log2(1.5),
  order_by = c("value", "name"),
  palette = "RdBu",
  palcolor = NULL,
  group_palette = "Chinese",
  group_palcolor = NULL,
  pt.size = 1,
  pt.alpha = 1,
  cols.sig = "red",
  cols.ns = "grey",
  cols.increase = "#d7301f",
  cols.decrease = "#2b8cbe",
  effect_color_mode = c("directional", "classic"),
  nlabel = 5,
  features_label = NULL,
  label = FALSE,
  label.fg = "black",
  label.bg = "white",
  label.bg.r = 0.1,
  label.size = 4,
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
  byrow = TRUE,
  seed = 11
)
}
\arguments{
\item{srt}{A Seurat object containing proportion test results.}

\item{comparison}{A character vector specifying which comparisons to plot.
If \code{NULL}, all stored comparisons are plotted.}

\item{proportion_method}{Optional method to select from method-layer results.
If \code{NULL}, active/most recent method is used.}

\item{result_level}{Result level for plotting: \code{"group"} or \code{"neighborhood"}.}

\item{plot_type}{Plot type. Main recommended values are \code{"effect"} and
\code{"umap"}. Legacy values are retained for compatibility:
\code{"volcano"}, \code{"manhattan"}, \code{"ring"},
\code{"milo_beeswarm"}, \code{"milo_graph"}, and \code{"sccoda_forest"}.}

\item{umap_mode}{UMAP mapping mode when \code{plot_type = "umap"}.
\code{"discrete"} maps cells to \code{Increased/Decreased/NS};
\code{"continuous"} maps group-level \code{obs_log2FD}.}

\item{reduction}{Reduction name used for UMAP projection.}

\item{projection_args}{Additional arguments passed through to
\code{CellDimPlot} (discrete mode) or \code{FeatureDimPlot} (continuous mode).}

\item{FDR_threshold}{FDR value cutoff for significance.}

\item{log2FD_threshold}{Absolute value of log2FD cutoff for significance.}

\item{order_by}{Method to order clusters. Options: \code{"name"}, \code{"value"}.}

\item{palette}{Palette name for continuous effect coloring.}

\item{palcolor}{Custom colors for \code{palette}.}

\item{group_palette}{Palette name for cluster/group coloring.}

\item{group_palcolor}{Custom colors for \code{group_palette}.}

\item{pt.size}{Point size.}

\item{pt.alpha}{Point transparency.}

\item{cols.sig}{Legacy color for significant/credible points and intervals.}

\item{cols.ns}{Color for non-significant points and intervals.}

\item{cols.increase}{Color for increased DA groups in directional effect mode.}

\item{cols.decrease}{Color for decreased DA groups in directional effect mode.}

\item{effect_color_mode}{Effect color style: \code{"directional"} (default)
or \code{"classic"}.}

\item{nlabel}{Number of labels when \code{label = TRUE} and
\code{features_label = NULL}.}

\item{features_label}{Character vector specifying labels to draw.}

\item{label}{Whether to draw labels.}

\item{label.fg}{Label foreground color.}

\item{label.bg}{Label background color.}

\item{label.bg.r}{Label background rounding radius.}

\item{label.size}{Label size.}

\item{aspect.ratio}{Aspect ratio of each panel.}

\item{xlab}{X-axis label.}

\item{ylab}{Y-axis label.}

\item{theme_use}{Theme used. Can be a character string or a theme function.}

\item{theme_args}{Additional arguments passed to \code{theme_use}.}

\item{legend.position}{Legend position.}

\item{legend.direction}{Legend direction.}

\item{legend.title}{Legend title.}

\item{combine}{Whether to combine into one patchwork object.}

\item{nrow}{Number of rows for combined plot.}

\item{ncol}{Number of columns for combined plot.}

\item{byrow}{Arrange panels by row.}

\item{seed}{Random seed used for jitter and label ranking ties.}
}
\description{
Generate differential abundance visualizations from outputs of
\link{RunProportionTest}. Main display styles are effect-size view and UMAP
projection view, while legacy plot types remain callable for compatibility.
}
\seealso{
\link{RunProportionTest}
}
