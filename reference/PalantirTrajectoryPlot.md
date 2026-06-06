# Plot Palantir trajectories

Plot branch-aware Palantir trajectories on a two-dimensional embedding.

## Usage

``` r
PalantirTrajectoryPlot(
  srt,
  reduction = NULL,
  dims = c(1, 2),
  cells = NULL,
  pseudotime_key = "palantir_pseudotime",
  branch_cols = NULL,
  diff_potential_key = "palantir_diff_potential",
  pseudotime_interval = c(0, 1),
  branch_min_prob = 0.05,
  n_bins = 60,
  min_cells_per_bin = 3,
  smooth = TRUE,
  trajectory_method = c("loess", "bin"),
  smoothness = 1,
  span = 0.75,
  n_path_points = 200,
  cell_color = "pseudotime",
  pt.size = 0.5,
  pt.alpha = 0.8,
  palette = "Dark2",
  palcolor = NULL,
  trajectory_palette = "Dark2",
  trajectory_palcolor = NULL,
  trajectory_linewidth = 1.2,
  trajectory_bg = "black",
  trajectory_bg_stroke = 0.7,
  trajectory_arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed"),
  aspect.ratio = 1,
  title = "Palantir",
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  return_layer = FALSE,
  seed = 11,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- reduction:

  Which dimensionality reduction to use. If not specified, will use the
  reduction returned by
  [DefaultReduction](https://mengxu98.github.io/scop/reference/DefaultReduction.md).

- dims:

  Dimensions to plot, must be a two-length numeric vector specifying x-
  and y-dimensions

- cells:

  A character vector of cell names to use.

- pseudotime_key:

  Name of the metadata column containing Palantir pseudotime.

- branch_cols:

  Metadata columns containing Palantir branch probabilities. If `NULL`,
  columns ending with `"_diff_potential"` are used, excluding
  `pseudotime_key` and `diff_potential_key`.

- diff_potential_key:

  Name of the Palantir entropy/differentiation potential column to
  exclude from branch auto-detection.

- pseudotime_interval:

  Numeric vector of length 2 specifying the pseudotime range to plot.

- branch_min_prob:

  Minimum branch probability used to select cells for a branch
  trajectory.

- n_bins:

  Number of pseudotime bins used to summarize each trajectory.

- min_cells_per_bin:

  Minimum number of cells required in a bin.

- smooth:

  Whether to smooth the trajectory with
  [stats::loess](https://rdrr.io/r/stats/loess.html).

- trajectory_method:

  Method used to fit trajectory coordinates along pseudotime. `"loess"`
  uses a fully R-native smoother over a Palantir-style pseudotime grid;
  `"bin"` uses binned median coordinates.

- smoothness:

  Smoothing multiplier for the R-native loess span. Higher values yield
  smoother curves.

- span:

  Base span used for loess smoothing.

- n_path_points:

  Number of pseudotime points used to draw each smoothed trajectory.

- cell_color:

  Cell coloring mode. Use `"pseudotime"` for Palantir pseudotime,
  `"branch_selection"` for the branch with highest probability, `"none"`
  to hide cell coloring, or any metadata column.

- pt.size:

  Point size for cells.

- pt.alpha:

  Point alpha for cells.

- palette:

  Color palette name. Available palettes can be found in
  [thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html).
  Default is `"Chinese"`.

- palcolor:

  Custom colors used to create a color palette. Default is `NULL`.

- trajectory_palette:

  Color palette for trajectories.

- trajectory_palcolor:

  Custom colors for trajectories.

- trajectory_linewidth:

  Line width of trajectories.

- trajectory_bg:

  Color for the trajectory background stroke.

- trajectory_bg_stroke:

  Width added to the trajectory background stroke.

- trajectory_arrow:

  Arrow used for trajectories. See
  [grid::arrow](https://rdrr.io/r/grid/arrow.html).

- aspect.ratio:

  Aspect ratio of the panel. Default is `1`.

- title:

  The text for the title. Default is `NULL`.

- subtitle:

  The text for the subtitle for the plot which will be displayed below
  the title. Default is `NULL`.

- xlab:

  The x-axis label of the plot. Default is `NULL`.

- ylab:

  The y-axis label of the plot. Default is `NULL`.

- legend.position:

  The position of legends, one of `"none"`, `"left"`, `"right"`,
  `"bottom"`, `"top"`. Default is `"right"`.

- legend.direction:

  The direction of the legend in the plot. Can be one of `"vertical"` or
  `"horizontal"`.

- theme_use:

  Theme used. Can be a character string or a theme function. Default is
  `"theme_scop"`.

- theme_args:

  Other arguments passed to the `theme_use`. Default is
  [`list()`](https://rdrr.io/r/base/list.html).

- return_layer:

  Logical. If `TRUE`, returns ggplot2 layers instead of a complete plot.

- seed:

  Random seed for reproducibility. Default is `11`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## See also

[RunPalantir](https://mengxu98.github.io/scop/reference/RunPalantir.md),
[FeatureDimPlot](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md),
[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
pancreas_sub <- RunPalantir(
  pancreas_sub,
  group.by = "SubCellType",
  linear_reduction = "PCA",
  nonlinear_reduction = "UMAP",
  early_group = "Ductal",
  terminal_groups = c("Alpha", "Beta", "Delta", "Epsilon")
)
PalantirTrajectoryPlot(
  pancreas_sub,
  reduction = "UMAP",
  pseudotime_interval = c(0, 0.9)
)
PalantirTrajectoryPlot(
  pancreas_sub,
  reduction = "UMAP",
  cell_color = "branch_selection",
  pseudotime_interval = c(0, 0.9)
)
} # }
```
