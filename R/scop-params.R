# Inherited via @inheritParams.

#' @title Shared parameters
#'
#' @md
#' @param srt A `Seurat` object.
#' @param assay Assay to use. `NULL` uses the default assay.
#' @param layer Assay layer to use.
#' @param verbose Print progress messages.
#' @param seed Random seed.
#' @param cores Number of CPU cores.
#' @param group.by Metadata column(s) to group cells by.
#' @param split.by Metadata column to split the analysis or plot by.
#' @param reduction Dimensionality reduction to use. `NULL` uses [DefaultReduction].
#' @param features Features to use (assay names or numeric metadata columns).
#' @param cells Cell names to include.
#' @param palette Color palette name. See [thisplot::show_palettes].
#' @param palcolor Custom colors used to build the palette.
#' @param species Species name, e.g. `"Homo_sapiens"` or `"Mus_musculus"`.
#' @param prefix Prefix for stored result names.
#' @param image Spatial image name. Required when multiple images are present;
#' a single image is selected automatically when `NULL`.
#' @param coord.cols Metadata coordinate columns used when no image is available.
#' @param coordinate_space Coordinate space for spatial distances. `"raw"` is
#' the default; `"legacy_display"` remains available for compatibility.
#' @param backend Implementation backend.
#' @param tool_name Name of the `srt@tools` slot used to store results.
#' @param overwrite Overwrite existing results.
#' @param ... Additional arguments passed to the underlying method.
#'
#' @keywords internal
#' @param byrow Whether to fill matrices/plots by row.
#' @param cell_annotation Cell annotations. Palette length should match `cell_annotation`.
#' @param cell_annotation_palcolor Custom colors for cell-type annotations.
#' @param cell_annotation_palette Color palette for cell-type annotations.
#' @param cells.highlight Cells to highlight and their appearance. `TRUE` highlights all cells.
#' @param cluster_columns Whether to cluster heatmap rows/columns. Defaults are both `FALSE`.
#' @param cluster_row_slices Whether to cluster row slices.
#' @param cluster_rows Whether to cluster heatmap rows/columns. Defaults are both `FALSE`.
#' @param cols.highlight Cells to highlight and their appearance. `TRUE` highlights all cells.
#' @param column_title The title for the column names in the heatmap. Default is to use the reference grouping variable.
#' @param combine Combine plots with [patchwork]. `combine = FALSE` returns a list of ggplots.
#' @param complex.mean Whether to use a modified mean test statistic in complex heatmaps.
#' @param coordinate.unit Unit of the raw input coordinates. Automatic unit selection uses technology-specific rules.
#' @param cutoff Maximum selected-cell fraction used by Scissor's alpha search.
#' @param decreasing Order groups decreasingly.
#' @param envname Conda-compatible Python environment. If `NULL`, the environment name will be set to `"scop_env"`.
#' @param feature_annotation_palcolor Custom colors for feature annotations.
#' @param feature_annotation_palette Color palette for feature annotations.
#' @param feature_split Feature splitting. `split_method` is `"kmeans"`, `"hclust"`, or `"mfuzz"`.
#' @param feature_split_palcolor Custom colors for feature-split annotation labels.
#' @param feature_split_palette Split colors.
#' @param global.fdr Global FDR threshold.
#' @param global.threshold Global threshold.
#' @param graph Neighbor-graph edges.
#' @param group_palcolor Custom colors for cell types (groups) in Manhattan plot.
#' @param group_palette Palette for cell types (groups) in Manhattan plot.
#' @param heatmap_palcolor Palette used for CNV heatmap values.
#' @param heatmap_palette Palette used for CNV heatmap values.
#' @param hex Whether to use hexagonal binning for 2D assignments.
#' @param hex.bins Number of hexagonal bins.
#' @param hex.binwidth Width of hexagonal bins.
#' @param hex.linewidth Line width of six-point hexagons.
#' @param label Group labels. `label_insitu = FALSE` uses numbers instead of group names.
#' @param label.bg Group labels. `label_insitu = FALSE` uses numbers instead of group names.
#' @param label.bg.r Group labels. `label_insitu = FALSE` uses numbers instead of group names.
#' @param label.fg Group labels. `label_insitu = FALSE` uses numbers instead of group names.
#' @param label.size Group labels. `label_insitu = FALSE` uses numbers instead of group names.
#' @param legend.direction Legend direction: `"horizontal"` or `"vertical"`.
#' @param legend.position Legend placement (`"none"`, `"left"`, `"right"`, `"bottom"`, `"top"`), direction, and title. `legend.title = NULL` uses the group name.
#' @param lineages Pseudotime columns to use. If `NULL`, lineage-like pseudotime columns such as `Lineage1`, `Lineage2`, `prefix_Lineage1`, or `pseudotime` are detected and merged into one global pseudotime for a single panel. Use `"all"` to plot each detected lineage in separate panels.
#' @param lineages_palette Color palette used for lineage groups.
#' @param linear_reduction Linear reduction (`"pca"`, `"svd"`, `"ica"`, `"nmf"`, `"mds"`, `"glmpca"`). `linear_reduction_dims_use = NULL` uses estimated dimensions, else the first 50.
#' @param local.fdr Local FDR threshold.
#' @param local.threshold Local threshold.
#' @param min_cell Minimum number of cells required.
#' @param n_neighbor_layers Number of neighbor layers to use.
#' @param n_neighbors Number of neighbors to use for constructing the KNN graph.
#' @param n_perm Number of permutations.
#' @param n_split Number of splits for feature dependencies.
#' @param ncol Number of columns of the combined plot.
#' @param nlabel The maximum number of labels to show on each side of the heatmap. If set to 0, no labels will be shown. This can be useful for reducing clutter in large heatmaps.
#' @param nonlinear_reduction Nonlinear reduction (`"umap"`, `"umap-naive"`, `"tsne"`, `"dm"`, `"phate"`, `"pacmap"`, `"trimap"`, `"largevis"`, `"fr"`).
#' @param nproc Number of parallel processes.
#' @param nrow Combine plots with [patchwork]. `combine = FALSE` returns a list of ggplots.
#' @param overlay_image Draw the spatial image beneath spots.
#' @param overwrite Whether to overwrite existing metadata.
#' @param palcolor Custom colors used to create a color palette.
#' @param palette Color palette name.
#' @param pt.alpha Point size and transparency. `pt.size = NULL` scales with `sqrt(n)` (minimum `0.3`). Rasterized points keep at least a two-pixel radius at `raster.dpi = c(512, 512)` and scale with `raster.dpi`.
#' @param pt.size Point size and transparency. `pt.size = NULL` scales with `sqrt(n)` (minimum `0.3`). Rasterized points keep at least a two-pixel radius at `raster.dpi = c(512, 512)` and scale with `raster.dpi`.
#' @param result.name Stored result name.
#' @param row_title The title for the row names in the heatmap. If not provided, the default is to use the query grouping variable.
#' @param run.local Whether to run locally.
#' @param show_column_names Whether to draw row/column names for the heatmap body.
#' @param show_row_names Whether to draw row/column names for the heatmap body.
#' @param single_cell Whether the input is a single-cell (rather than spatial) dataset.
#' @param sizes.highlight Cells to highlight and their appearance. `TRUE` highlights all cells.
#' @param split_method Method for splitting features into groups.
#' @param split_order Order of feature split groups (e.g. `c("A", "B")`).
#' @param subtitle Plot subtitle.
#' @param theme_args Theme name or function, plus extra theme arguments.
#' @param title Plot title. `NULL` hides the title for merged/single panels. When multiple lineages are plotted and `title` is `NULL`, each panel is titled with its lineage column.
#' @param verbose Whether to print messages.
#' @param xlab Plot labels.
#' @param ylab Plot labels.

#' @name scop-params
NULL
