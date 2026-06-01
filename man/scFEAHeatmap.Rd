\name{scFEAHeatmap}
\alias{scFEAHeatmap}
\title{Plot scFEA module flux heatmap}
\usage{
scFEAHeatmap(
  srt,
  assay = "scFEAflux",
  layer = "data",
  group.by = NULL,
  features = NULL,
  modules = NULL,
  label_by = c("module", "reaction", "module_reaction"),
  add_sm_anno = TRUE,
  scale_rows = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  sm_anno_label_rot = 0,
  show_row_names = FALSE,
  heatmap_limit = 2,
  heatmap_column_width = grid::unit(10, "mm"),
  column_names_rot = 45,
  column_names_gp = grid::gpar(fontsize = 8),
  column_title_gp = grid::gpar(fontsize = 11),
  row_title_gp = grid::gpar(fontsize = 7.5),
  legend_title_gp = grid::gpar(fontsize = 8),
  legend_labels_gp = grid::gpar(fontsize = 7),
  row_gap = grid::unit(1.2, "mm"),
  sm_anno_size = grid::unit(2.5, "mm"),
  group_anno_size = grid::unit(3, "mm"),
  mark_features = NULL,
  mark_label_by = NULL,
  mark_labels_gp = NULL,
  mark_annotation_width = NULL,
  mark_link_width = grid::unit(5, "mm"),
  name = NULL,
  width = NULL,
  height = NULL,
  verbose = TRUE
)
}
\arguments{
\item{srt}{A Seurat object returned by \code{\link[=RunscFEA]{RunscFEA()}}.}
\item{assay}{Flux assay name.}
\item{layer}{Flux assay layer.}
\item{group.by}{Metadata column used to aggregate cells. If \code{NULL}, all cells are averaged together.}
\item{features}{Optional scFEA modules to plot. Supports Seurat-safe feature ids such as \code{"M-1"}, scFEA ids such as \code{"M_1"}, labels such as \code{"M1"}, or reaction labels.}
\item{modules}{Optional scFEA modules to plot. This is a user-facing alias for \code{features} that also accepts combined labels such as \code{"M150: PRPP -> UMP"}.}
\item{label_by}{Row label mode: module id, reaction, or both.}
\item{add_sm_anno}{Whether to add \code{SM_anno} row annotation and split rows by pathway class.}
\item{scale_rows}{Whether to z-score each row after group aggregation.}
\item{cluster_rows, cluster_columns}{Passed to \code{\link[ComplexHeatmap:Heatmap]{ComplexHeatmap::Heatmap()}}.}
\item{sm_anno_label_rot}{Rotation angle for \code{SM_anno} row-split labels.}
\item{show_row_names}{Whether to show row labels.}
\item{heatmap_limit}{Numeric clipping limit for row z-scores.}
\item{heatmap_column_width}{Width for each aggregated group column. Numeric values are interpreted as millimeters.}
\item{column_names_rot}{Rotation angle for heatmap column labels.}
\item{column_names_gp, column_title_gp, row_title_gp}{Font settings passed to \code{\link[grid:gpar]{grid::gpar()}}-aware ComplexHeatmap arguments.}
\item{legend_title_gp, legend_labels_gp}{Font settings for heatmap and annotation legends.}
\item{row_gap}{Gap between \code{SM_anno} row splits.}
\item{sm_anno_size}{Size of the left \code{SM_anno} annotation strip.}
\item{group_anno_size}{Size of the top group annotation strip.}
\item{mark_features}{Optional scFEA modules to mark with \code{\link[ComplexHeatmap:anno_mark]{ComplexHeatmap::anno_mark()}}. Accepts the same module ids and labels as \code{features}.}
\item{mark_label_by}{Label mode for marked rows. Defaults to \code{label_by}.}
\item{mark_labels_gp}{Font settings for marked row labels. If \code{NULL}, a compact default is chosen from \code{mark_label_by}.}
\item{mark_annotation_width}{Width of the right mark annotation. Numeric values are interpreted as millimeters.}
\item{mark_link_width}{Link width used by \code{\link[ComplexHeatmap:anno_mark]{ComplexHeatmap::anno_mark()}}.}
\item{name}{Heatmap legend title.}
\item{width, height}{Suggested export size. Stored as attributes on the returned heatmap object.}
\item{verbose}{Whether to print messages.}
}
\value{
A \pkg{ComplexHeatmap} heatmap object.
}
\description{
Aggregates scFEA module flux by Seurat metadata groups and draws a heatmap with optional M168 supermodule annotation.
}
