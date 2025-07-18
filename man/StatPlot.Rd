% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CellStatPlot.R
\name{StatPlot}
\alias{StatPlot}
\title{StatPlot}
\usage{
StatPlot(
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
  plot_type = c("bar", "rose", "ring", "pie", "trend", "area", "dot", "sankey", "chord",
    "venn", "upset"),
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
  seed = 11
)
}
\arguments{
\item{meta.data}{The data frame containing the data to be plotted.}

\item{stat.by}{The column name(s) in \code{meta.data} specifying the variable(s) to be plotted.}

\item{group.by}{The column name in \code{meta.data} specifying the grouping variable.}

\item{split.by}{The column name in \code{meta.data} specifying the splitting variable.}

\item{bg.by}{The column name in \code{meta.data} specifying the background variable for bar plots.}

\item{flip}{Logical indicating whether to flip the plot.}

\item{NA_color}{The color to use for missing values.}

\item{NA_stat}{Logical indicating whether to include missing values in the plot.}

\item{keep_empty}{Logical indicating whether to keep empty groups in the plot.}

\item{individual}{Logical indicating whether to plot individual groups separately.}

\item{stat_level}{The level(s) of the variable(s) specified in \code{stat.by} to include in the plot.}

\item{plot_type}{The type of plot to create.
Can be one of "bar", "rose", "ring", "pie", "trend", "area", "dot", "sankey", "chord", "venn", or "upset".}

\item{stat_type}{The type of statistic to compute for the plot.
Can be one of "percent" or "count".}

\item{position}{The position adjustment for the plot.
Can be one of "stack" or "dodge".}

\item{palette}{The name of the color palette to use for the plot.}

\item{palcolor}{The color to use in the color palette.}

\item{alpha}{The transparency level for the plot.}

\item{bg_palette}{The name of the background color palette to use for bar plots.}

\item{bg_palcolor}{The color to use in the background color palette.}

\item{bg_alpha}{The transparency level for the background color in bar plots.}

\item{label}{Logical indicating whether to add labels on the plot.}

\item{label.size}{The size of the labels.}

\item{label.fg}{The foreground color of the labels.}

\item{label.bg}{The background color of the labels.}

\item{label.bg.r}{The radius of the rounded corners of the label background.}

\item{aspect.ratio}{The aspect ratio of the plot.}

\item{title}{The main title of the plot.}

\item{subtitle}{The subtitle of the plot.}

\item{xlab}{The x-axis label of the plot.}

\item{ylab}{The y-axis label of the plot.}

\item{legend.position}{The position of the legend in the plot.
Can be one of "right", "left", "bottom", "top", or "none".}

\item{legend.direction}{The direction of the legend in the plot. Can be one of "vertical" or "horizontal".}

\item{theme_use}{The name of the theme to use for the plot. Can be one of the predefined themes or a custom theme.}

\item{theme_args}{A list of arguments to be passed to the theme function.}

\item{combine}{Logical indicating whether to combine multiple plots into a single plot.}

\item{nrow}{The number of rows in the combined plot.}

\item{ncol}{The number of columns in the combined plot.}

\item{byrow}{Logical indicating whether to fill the plot by row or by column.}

\item{force}{Logical indicating whether to force the plot even if some variables have more than 100 levels.}

\item{seed}{The random seed to use for reproducible results.}
}
\description{
Visualizes data using various plot types such as bar plots,
rose plots, ring plots, pie charts, trend plots, area plots,
dot plots, sankey plots, chord plots, venn diagrams, and upset plots.
}
\examples{
data(pancreas_sub)
head(pancreas_sub@meta.data)
StatPlot(
  pancreas_sub@meta.data,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "bar",
  label = TRUE
)

StatPlot(
  pancreas_sub[["RNA"]]@meta.data,
  stat.by = "highly_variable_genes",
  plot_type = "ring",
  label = TRUE,
  NA_stat = FALSE
)

pancreas_sub <- AnnotateFeatures(
  pancreas_sub,
  species = "Mus_musculus",
  IDtype = "symbol",
  # db = c("VerSeDa", "TF")
  db = "VerSeDa"
)
StatPlot(
  GetFeaturesData(pancreas_sub, "RNA"),
  stat.by = "TF",
  group.by = "VerSeDa",
  stat_type = "count",
  plot_type = "bar",
  position = "dodge",
  label = TRUE,
  NA_stat = FALSE
)
}
\seealso{
\link{CellStatPlot}
}
