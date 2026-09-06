# Shared parameters

Shared parameters

## Arguments

- srt:

  A `Seurat` object.

- assay:

  Assay to use. `NULL` uses the default assay.

- layer:

  Assay layer to use.

- seed:

  Random seed.

- cores:

  Number of CPU cores.

- group.by:

  Metadata column(s) to group cells by.

- split.by:

  Metadata column to split the analysis or plot by.

- reduction:

  Dimensionality reduction to use. `NULL` uses
  [DefaultReduction](https://mengxu98.github.io/scop/reference/DefaultReduction.md).

- features:

  Features to use (assay names or numeric metadata columns).

- cells:

  Cell names to include.

- species:

  Species name, e.g. `"Homo_sapiens"` or `"Mus_musculus"`.

- prefix:

  Prefix for stored result names.

- image:

  Spatial image name. Required when multiple images are present; a
  single image is selected automatically when `NULL`.

- coord.cols:

  Metadata coordinate columns used when no image is available.

- coordinate_space:

  Coordinate space for spatial distances. `"raw"` is the default;
  `"legacy_display"` remains available for compatibility.

- backend:

  Implementation backend.

- tool_name:

  Name of the `srt@tools` slot used to store results.

- ...:

  Additional arguments passed to the underlying method.

- byrow:

  Whether to fill matrices/plots by row.

- cell_annotation:

  Cell annotations. Palette length should match `cell_annotation`.

- cell_annotation_palcolor:

  Custom colors for cell-type annotations.

- cell_annotation_palette:

  Color palette for cell-type annotations.

- cells.highlight:

  Cells to highlight and their appearance. `TRUE` highlights all cells.

- cluster_columns:

  Whether to cluster heatmap rows/columns. Defaults are both `FALSE`.

- cluster_row_slices:

  Whether to cluster row slices.

- cluster_rows:

  Whether to cluster heatmap rows/columns. Defaults are both `FALSE`.

- cols.highlight:

  Cells to highlight and their appearance. `TRUE` highlights all cells.

- column_title:

  The title for the column names in the heatmap. Default is to use the
  reference grouping variable.

- combine:

  Combine plots with
  [patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html).
  `combine = FALSE` returns a list of ggplots.

- complex.mean:

  Whether to use a modified mean test statistic in complex heatmaps.

- coordinate.unit:

  Unit of the raw input coordinates. Automatic unit selection uses
  technology-specific rules.

- cutoff:

  Maximum selected-cell fraction used by Scissor's alpha search.

- decreasing:

  Order groups decreasingly.

- envname:

  Conda-compatible Python environment. If `NULL`, the environment name
  will be set to `"scop_env"`.

- feature_annotation_palcolor:

  Custom colors for feature annotations.

- feature_annotation_palette:

  Color palette for feature annotations.

- feature_split:

  Feature splitting. `split_method` is `"kmeans"`, `"hclust"`, or
  `"mfuzz"`.

- feature_split_palcolor:

  Custom colors for feature-split annotation labels.

- feature_split_palette:

  Split colors.

- global.fdr:

  Global FDR threshold.

- global.threshold:

  Global threshold.

- graph:

  Neighbor-graph edges.

- group_palcolor:

  Custom colors for cell types (groups) in Manhattan plot.

- group_palette:

  Palette for cell types (groups) in Manhattan plot.

- heatmap_palcolor:

  Palette used for CNV heatmap values.

- heatmap_palette:

  Palette used for CNV heatmap values.

- hex:

  Whether to use hexagonal binning for 2D assignments.

- hex.bins:

  Number of hexagonal bins.

- hex.binwidth:

  Width of hexagonal bins.

- hex.linewidth:

  Line width of six-point hexagons.

- label:

  Group labels. `label_insitu = FALSE` uses numbers instead of group
  names.

- label.bg:

  Group labels. `label_insitu = FALSE` uses numbers instead of group
  names.

- label.bg.r:

  Group labels. `label_insitu = FALSE` uses numbers instead of group
  names.

- label.fg:

  Group labels. `label_insitu = FALSE` uses numbers instead of group
  names.

- label.size:

  Group labels. `label_insitu = FALSE` uses numbers instead of group
  names.

- legend.direction:

  Legend direction: `"horizontal"` or `"vertical"`.

- legend.position:

  Legend placement (`"none"`, `"left"`, `"right"`, `"bottom"`, `"top"`),
  direction, and title. `legend.title = NULL` uses the group name.

- lineages:

  Pseudotime columns to use. If `NULL`, lineage-like pseudotime columns
  such as `Lineage1`, `Lineage2`, `prefix_Lineage1`, or `pseudotime` are
  detected and merged into one global pseudotime for a single panel. Use
  `"all"` to plot each detected lineage in separate panels.

- lineages_palette:

  Color palette used for lineage groups.

- linear_reduction:

  Linear reduction (`"pca"`, `"svd"`, `"ica"`, `"nmf"`, `"mds"`,
  `"glmpca"`). `linear_reduction_dims_use = NULL` uses estimated
  dimensions, else the first 50.

- local.fdr:

  Local FDR threshold.

- local.threshold:

  Local threshold.

- min_cell:

  Minimum number of cells required.

- n_neighbor_layers:

  Number of neighbor layers to use.

- n_neighbors:

  Number of neighbors to use for constructing the KNN graph.

- n_perm:

  Number of permutations.

- n_split:

  Number of splits for feature dependencies.

- ncol:

  Number of columns of the combined plot.

- nlabel:

  The maximum number of labels to show on each side of the heatmap. If
  set to 0, no labels will be shown. This can be useful for reducing
  clutter in large heatmaps.

- nonlinear_reduction:

  Nonlinear reduction (`"umap"`, `"umap-naive"`, `"tsne"`, `"dm"`,
  `"phate"`, `"pacmap"`, `"trimap"`, `"largevis"`, `"fr"`).

- nproc:

  Number of parallel processes.

- nrow:

  Combine plots with
  [patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html).
  `combine = FALSE` returns a list of ggplots.

- overlay_image:

  Draw the spatial image beneath spots.

- overwrite:

  Whether to overwrite existing metadata.

- palcolor:

  Custom colors used to create a color palette.

- palette:

  Color palette name.

- pt.alpha:

  Point size and transparency. `pt.size = NULL` scales with `sqrt(n)`
  (minimum `0.3`). Rasterized points keep at least a two-pixel radius at
  `raster.dpi = c(512, 512)` and scale with `raster.dpi`.

- pt.size:

  Point size and transparency. `pt.size = NULL` scales with `sqrt(n)`
  (minimum `0.3`). Rasterized points keep at least a two-pixel radius at
  `raster.dpi = c(512, 512)` and scale with `raster.dpi`.

- result.name:

  Stored result name.

- row_title:

  The title for the row names in the heatmap. If not provided, the
  default is to use the query grouping variable.

- run.local:

  Whether to run locally.

- show_column_names:

  Whether to draw row/column names for the heatmap body.

- show_row_names:

  Whether to draw row/column names for the heatmap body.

- single_cell:

  Whether the input is a single-cell (rather than spatial) dataset.

- sizes.highlight:

  Cells to highlight and their appearance. `TRUE` highlights all cells.

- split_method:

  Method for splitting features into groups.

- split_order:

  Order of feature split groups (e.g. `c("A", "B")`).

- subtitle:

  Plot subtitle.

- theme_args:

  Theme name or function, plus extra theme arguments.

- title:

  Plot title. `NULL` hides the title for merged/single panels. When
  multiple lineages are plotted and `title` is `NULL`, each panel is
  titled with its lineage column.

- verbose:

  Whether to print messages.

- xlab:

  Plot labels.

- ylab:

  Plot labels.
