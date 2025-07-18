% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EnrichmentPlot.R
\name{EnrichmentPlot}
\alias{EnrichmentPlot}
\title{EnrichmentPlot}
\usage{
EnrichmentPlot(
  srt,
  db = "GO_BP",
  group_by = NULL,
  test.use = "wilcox",
  res = NULL,
  plot_type = c("bar", "dot", "lollipop", "network", "enrichmap", "wordcloud",
    "comparison"),
  split_by = c("Database", "Groups"),
  color_by = "Database",
  group_use = NULL,
  id_use = NULL,
  pvalueCutoff = NULL,
  padjustCutoff = 0.05,
  topTerm = ifelse(plot_type == "enrichmap", 100, 6),
  compare_only_sig = FALSE,
  topWord = 100,
  word_type = c("term", "feature"),
  word_size = c(2, 8),
  words_excluded = NULL,
  network_layout = "fr",
  network_labelsize = 5,
  network_blendmode = "blend",
  network_layoutadjust = TRUE,
  network_adjscale = 60,
  network_adjiter = 100,
  enrichmap_layout = "fr",
  enrichmap_cluster = "fast_greedy",
  enrichmap_label = c("term", "feature"),
  enrichmap_labelsize = 5,
  enrlichmap_nlabel = 4,
  enrichmap_show_keyword = FALSE,
  enrichmap_mark = c("ellipse", "hull"),
  enrichmap_expand = c(0.5, 0.5),
  character_width = 50,
  lineheight = 0.5,
  palette = "Spectral",
  palcolor = NULL,
  aspect.ratio = 1,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  seed = 11
)
}
\arguments{
\item{srt}{A Seurat object containing the results of RunDEtest and RunEnrichment.
If specified, enrichment results will be extracted from the Seurat object automatically.
If not specified, the \code{res} arguments must be provided.}

\item{db}{The database to use for enrichment plot. Default is "GO_BP".}

\item{group_by}{A character vector specifying the grouping variable in the Seurat object.
This argument is only used if \code{srt} is specified.}

\item{test.use}{A character vector specifying the test to be used in differential expression analysis.
This argument is only used if \code{srt} is specified.}

\item{res}{Enrichment results generated by RunEnrichment function.
If provided, 'srt', 'test.use' and 'group_by' are ignored.}

\item{plot_type}{The type of plot to generate.
Options are: "bar", "dot", "lollipop", "network", "enrichmap", "wordcloud", "comparison".
Default is "bar".}

\item{split_by}{The splitting variable(s) for the plot.
Can be "Database", "Groups", or both.
Default is c("Database", "Groups") for plots.}

\item{color_by}{The variable used for coloring.
Default is "Database".}

\item{group_use}{The group(s) to be used for enrichment plot.
Default is NULL.}

\item{id_use}{List of IDs to be used to display specific terms in the enrichment plot.
Default value is NULL.}

\item{pvalueCutoff}{The p-value cutoff.
Only work when \code{padjustCutoff} is NULL.
Default is NULL.}

\item{padjustCutoff}{The p-adjusted cutoff.
Default is 0.05.}

\item{topTerm}{The number of top terms to display. Default is 6, or 100 if 'plot_type' is "enrichmap".}

\item{compare_only_sig}{Whether to compare only significant terms.
Default is FALSE.}

\item{topWord}{The number of top words to display for wordcloud.
Default is 100.}

\item{word_type}{The type of words to display in wordcloud.
Options are "term" and "feature".
Default is "term".}

\item{word_size}{The size range for words in wordcloud. Default is c(2, 8).}

\item{words_excluded}{Words to be excluded from the wordcloud.
The default value is NULL, which means that the built-in words (\code{\link{words_excluded}}) will be used.}

\item{network_layout}{The layout algorithm to use for network plot.
Options are "fr", "kk","random", "circle", "tree", "grid", or other algorithm from 'igraph' package.
Default is "fr".}

\item{network_labelsize}{The label size for network plot.
Default is 5.}

\item{network_blendmode}{The blend mode for network plot.
Default is "blend".}

\item{network_layoutadjust}{Whether to adjust the layout of the network plot to avoid overlapping words.
Default is TRUE.}

\item{network_adjscale}{The scale for adjusting network plot layout.
Default is 60.}

\item{network_adjiter}{The number of iterations for adjusting network plot layout.
Default is 100.}

\item{enrichmap_layout}{The layout algorithm to use for enrichmap plot.
Options are "fr", "kk","random", "circle", "tree", "grid", or other algorithm from 'igraph' package.
Default is "fr".}

\item{enrichmap_cluster}{The clustering algorithm to use for enrichmap plot.
Options are "walktrap", "fast_greedy", or other algorithm from 'igraph' package.
Default is "fast_greedy".}

\item{enrichmap_label}{The label type for enrichmap plot.
Options are "term" and "feature". Default is "term".}

\item{enrichmap_labelsize}{The label size for enrichmap plot.
Default is 5.}

\item{enrlichmap_nlabel}{The number of labels to display for each cluster in enrichmap plot.
Default is 4.}

\item{enrichmap_show_keyword}{Whether to show the keyword of terms or features in enrichmap plot.
Default is FALSE.}

\item{enrichmap_mark}{The mark shape for enrichmap plot.
Options are "ellipse" and "hull".
Default is "ellipse".}

\item{enrichmap_expand}{The expansion factor for enrichmap plot.
Default is c(0.5, 0.5).}

\item{character_width}{The maximum width of character of descriptions.
Default is 50.}

\item{lineheight}{The line height for y-axis labels.
Default is 0.5.}

\item{palette}{The color palette to use.
Default is "Spectral".}

\item{palcolor}{Custom colors for palette.
Default is NULL.}

\item{aspect.ratio}{The aspect ratio of the plot.
Default is 1.}

\item{legend.position}{The position of the legend.
Default is "right".}

\item{legend.direction}{The direction of the legend.
Default is "vertical".}

\item{theme_use}{The theme to use for the plot.
Default is "theme_scop".}

\item{theme_args}{The arguments to pass to the theme.
Default is an empty list.}

\item{combine}{Whether to combine multiple plots into a single plot.
Default is TRUE.}

\item{nrow}{The number of rows in the combined plot.
Default is NULL, calculated based on the number of plots.}

\item{ncol}{The number of columns in the combined plot.
Default is NULL, calculated based on the number of plots.}

\item{byrow}{Whether to fill the combined plot by row.
Default is TRUE.}

\item{seed}{The random seed to use. Default is 11.}
}
\description{
This function generates various types of plots for enrichment (over-representation) analysis.
}
\examples{
data(pancreas_sub)
# pancreas_sub <- RunDEtest(
#   pancreas_sub,
#   group_by = "CellType"
# )

pancreas_sub <- RunEnrichment(
  srt = pancreas_sub,
  db = c("GO_BP", "GO_CC"),
  group_by = "CellType",
  species = "Mus_musculus"
)

EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP",
  group_by = "CellType",
  group_use = "Ductal",
  plot_type = "bar"
)

EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP",
  group_by = "CellType",
  plot_type = "bar",
  color_by = "Groups",
  ncol = 2
)

EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP",
  group_by = "CellType",
  plot_type = "bar",
  id_use = list(
    "Ductal" = c(
      "GO:0002181", "GO:0045787",
      "GO:0006260", "GO:0050679"
    ),
    "Ngn3-low-EP" = c(
      "GO:0050678", "GO:0051101",
      "GO:0072091", "GO:0006631"
    ),
    "Ngn3-high-EP" = c(
      "GO:0035270", "GO:0030325",
      "GO:0008637", "GO:0030856"
    ),
    "Pre-endocrine" = c(
      "GO:0090276", "GO:0031018",
      "GO:0030073", "GO:1903532"
    )
  )
)

EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP",
  group_by = "CellType",
  topTerm = 3,
  plot_type = "comparison"
)

EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP",
  group_by = "CellType",
  topTerm = 3,
  plot_type = "comparison",
  compare_only_sig = TRUE
)

EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP",
  group_by = "CellType",
  group_use = c("Ductal", "Endocrine"),
  plot_type = "comparison"
)

EnrichmentPlot(
  pancreas_sub,
  db = c("GO_BP", "GO_CC"),
  group_by = "CellType",
  group_use = c("Ductal", "Endocrine"),
  plot_type = "bar",
  split_by = "Groups"
)

EnrichmentPlot(
  pancreas_sub,
  db = c("GO_BP", "GO_CC"),
  group_by = "CellType",
  group_use = c("Ductal", "Endocrine"),
  plot_type = "bar",
  split_by = "Database",
  color_by = "Groups"
)

EnrichmentPlot(
  pancreas_sub,
  db = c("GO_BP", "GO_CC"),
  group_by = "CellType",
  group_use = c("Ductal", "Endocrine"),
  plot_type = "bar",
  split_by = c("Database", "Groups")
)

EnrichmentPlot(
  pancreas_sub,
  db = c("GO_BP", "GO_CC"),
  group_by = "CellType",
  group_use = c("Ductal", "Endocrine"),
  plot_type = "bar",
  split_by = c("Groups", "Database")
)
EnrichmentPlot(
  pancreas_sub,
  db = c("GO_BP", "GO_CC"),
  group_by = "CellType",
  plot_type = "bar",
  split_by = "Database",
  color_by = "Groups",
  palette = "Set1"
)

EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP",
  group_by = "CellType",
  group_use = "Ductal",
  plot_type = "dot",
  palette = "GdRd"
)

EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP",
  group_by = "CellType",
  group_use = "Ductal",
  plot_type = "lollipop",
  palette = "GdRd"
)
EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP",
  group_by = "CellType",
  group_use = "Ductal",
  plot_type = "wordcloud"
)
EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP",
  group_by = "CellType",
  group_use = "Ductal",
  plot_type = "wordcloud",
  word_type = "feature"
)
EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP",
  group_by = "CellType",
  group_use = "Ductal",
  plot_type = "network"
)
EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP",
  group_by = "CellType",
  group_use = "Ductal",
  plot_type = "network",
  id_use = c(
    "GO:0050678",
    "GO:0035270",
    "GO:0090276",
    "GO:0030073"
  )
)
EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP",
  group_by = "CellType",
  group_use = "Ductal",
  plot_type = "network",
  network_layoutadjust = FALSE
)

EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP",
  group_by = "CellType",
  group_use = "Ductal",
  plot_type = "network",
  topTerm = 4,
  network_blendmode = "average",
  theme_use = "theme_blank",
  theme_args = list(add_coord = FALSE)
) |> panel_fix(height = 5)

EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP",
  group_by = "CellType",
  group_use = "Ductal",
  plot_type = "enrichmap"
)

EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP",
  group_by = "CellType",
  group_use = "Ductal",
  plot_type = "enrichmap",
  enrichmap_expand = c(2, 1)
)

EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP",
  group_by = "CellType",
  group_use = "Ductal",
  plot_type = "enrichmap",
  enrichmap_show_keyword = TRUE,
  character_width = 10
)

EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP",
  group_by = "CellType",
  group_use = "Ductal",
  plot_type = "enrichmap",
  topTerm = 200,
  enrichmap_mark = "hull",
  enrichmap_label = "feature",
  enrlichmap_nlabel = 3,
  character_width = 10,
  theme_use = "theme_blank",
  theme_args = list(add_coord = FALSE)
) |> panel_fix(height = 4)

pancreas_sub <- RunEnrichment(
  srt = pancreas_sub,
  db = c("MP", "DO"),
  group_by = "CellType",
  convert_species = TRUE,
  species = "Mus_musculus"
)

EnrichmentPlot(
  pancreas_sub,
  db = c("MP", "DO"),
  group_by = "CellType",
  group_use = "Ductal",
  ncol = 1
)
}
\seealso{
\code{\link{RunEnrichment}}
}
