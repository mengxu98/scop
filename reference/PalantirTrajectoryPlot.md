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

  A `Seurat` object.

- reduction:

  Reduction to plot. `NULL` uses
  [DefaultReduction](https://mengxu98.github.io/scop/reference/DefaultReduction.md).

- dims:

  Length-2 vector of dimensions to plot.

- cells:

  Cell names to include.

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

- palette, palcolor:

  Palette name
  ([thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html))
  or custom colors.

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

  Panel aspect ratio.

- title, subtitle, xlab, ylab:

  Plot labels.

- legend.position:

  Legend placement (`"none"`, `"left"`, `"right"`, `"bottom"`, `"top"`),
  direction, and title. `legend.title = NULL` uses the group name.

- legend.direction:

  Legend direction: `"horizontal"` or `"vertical"`.

- theme_use, theme_args:

  Theme name or function, plus extra theme arguments.

- return_layer:

  Logical. If `TRUE`, returns ggplot2 layers instead of a complete plot.

- seed:

  Random seed.

- verbose:

  Whether to print messages.

## See also

[RunPalantir](https://mengxu98.github.io/scop/reference/RunPalantir.md),
[FeatureDimPlot](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md),
[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
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
