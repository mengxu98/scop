# Plot transcription factor activity

Visualize TF activity stored by
[`RunDorothea()`](https://mengxu98.github.io/scop/reference/RunDorothea.md).
Comparison plots (`"bar"`, `"lollipop"`, `"volcano"`) test two groups
and show the mean activity difference `group1 - group2`. `"dim"`
compares TF activity with TF expression on embeddings, and `"targets"`
shows a regulon-target volcano for one TF. Use
[GroupHeatmap](https://mengxu98.github.io/scop/reference/GroupHeatmap.md),
[FeatureStatPlot](https://mengxu98.github.io/scop/reference/FeatureStatPlot.md),
and
[FeatureDimPlot](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md)
directly for ordinary activity visualizations.

## Usage

``` r
DorotheaPlot(
  object,
  group.by,
  group1 = NULL,
  group2 = NULL,
  plot_type = c("bar", "lollipop", "volcano", "dim", "targets"),
  tool_name = "Dorothea",
  assay_name = NULL,
  features = NULL,
  top_n = 30,
  test.use = c("wilcox.test", "t.test"),
  p.adjust.method = "BH",
  color.by = c("p_val", "p_val_adj"),
  rank.by = c("abs_logFC", "p_val", "p_val_adj", "logFC"),
  sort.by = c("logFC", "abs_logFC", "p_val", "p_val_adj"),
  p_floor = .Machine$double.xmin,
  padjustCutoff = 0.05,
  palette = "RdBu",
  palcolor = NULL,
  dim_args = list(),
  compare_expression = TRUE,
  expression_assay = NULL,
  expression_layer = "data",
  nlabel = 10,
  reduction = NULL,
  bar_width = 0.72,
  pt.size = 3.2,
  point_size = NULL,
  aspect.ratio = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  fill.title = NULL,
  flip = TRUE,
  cols = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  return_data = FALSE,
  verbose = TRUE,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object containing results from
  [`RunDorothea()`](https://mengxu98.github.io/scop/reference/RunDorothea.md).

- group.by:

  Metadata column used to define groups.

- group1, group2:

  Two group labels to compare. Required for `"bar"`, `"lollipop"`,
  `"volcano"`, and `"targets"`. Positive logFC means higher TF activity
  (or target expression) in `group1`. Ignored for `"dim"`.

- plot_type:

  Plot type. `"bar"` and `"lollipop"` show signed activity differences,
  `"volcano"` shows TF-level significance versus logFC, `"dim"` compares
  activity and expression on embeddings, and `"targets"` shows a
  regulon-target volcano.

- tool_name:

  Name of the `srt@tools` entry created by
  [`RunDorothea()`](https://mengxu98.github.io/scop/reference/RunDorothea.md).

- assay_name:

  Activity assay used for `"dim"`. If `NULL`, the assay stored by
  [`RunDorothea()`](https://mengxu98.github.io/scop/reference/RunDorothea.md)
  is used, or `"dorothea"`.

- features:

  TFs to plot. If `NULL`, `"bar"`/`"lollipop"` use the top `top_n` TFs,
  `"volcano"` shows all tested TFs, `"dim"` uses TFs present in both
  activity and expression assays, and `"targets"` uses the TF with the
  largest absolute activity difference.

- top_n:

  Number of TFs to show when `features = NULL`. Ignored for `"volcano"`,
  which always plots every tested TF. Set `NULL` to show all tested TFs
  in other comparison plots.

- test.use:

  Statistical test used for each TF or target gene in comparison plots.

- p.adjust.method:

  Method passed to
  [stats::p.adjust](https://rdrr.io/r/stats/p.adjust.html).

- color.by:

  P-value column used for comparison-plot color scales.

- rank.by:

  Metric used to select top TFs when `features = NULL` in comparison
  plots.

- sort.by:

  Metric used to order TFs in `"bar"` and `"lollipop"` plots.

- p_floor:

  Lower bound used before `-log10()` transformation.

- padjustCutoff:

  Adjusted p-value cutoff used to mark significant TFs in `"volcano"`
  plots and non-supporting targets in `"targets"` plots.

- palette, palcolor:

  Palette passed to `palette_colors()` for comparison plots.

- dim_args:

  Additional arguments passed to
  [`FeatureDimPlot()`](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md)
  when `plot_type = "dim"`.

- compare_expression:

  Whether `"dim"` also plots TF gene expression beside TF activity.

- expression_assay, expression_layer:

  Assay and layer used for TF expression in `"dim"` and for target genes
  in `"targets"`.

- nlabel:

  Number of significant TFs labeled in `"volcano"` plots, or significant
  target genes labeled in `"targets"` plots.

- reduction:

  Reduction used by `"dim"` plots. If `NULL`, the default reduction of
  `object` is used.

- bar_width:

  Width of bars in `"bar"` plots.

- pt.size:

  Point size in `"lollipop"`, `"volcano"`, and `"targets"` plots.

- point_size:

  Deprecated alias(es) for `pt.size`; supply exactly one of the two. It
  will be removed in scop 1.0.0.

- aspect.ratio:

  Panel aspect ratio.

- legend.position:

  Legend placement (`"none"`, `"left"`, `"right"`, `"bottom"`, `"top"`),
  direction, and title. `legend.title = NULL` uses the group name.

- legend.direction:

  Legend direction: `"horizontal"` or `"vertical"`.

- title, xlab, ylab, fill.title:

  Axis, plot, and legend titles.

- flip:

  Whether to draw comparison bar/lollipop plots horizontally.

- cols:

  Optional three-color vector used instead of `palette` for diverging
  comparison scales. Kept for backward compatibility.

- theme_use, theme_args:

  Theme name or function, plus extra theme arguments.

- return_data:

  Whether to return a list with the plot and statistics.

- verbose:

  Whether to print the message. Default is `TRUE`.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

A `ggplot` or patchwork object, or a list with `plot` and `data` when
`return_data = TRUE`.

## See also

[RunDorothea](https://mengxu98.github.io/scop/reference/RunDorothea.md),
[GroupHeatmap](https://mengxu98.github.io/scop/reference/GroupHeatmap.md),
[FeatureDimPlot](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md),
[FeatureStatPlot](https://mengxu98.github.io/scop/reference/FeatureStatPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub, verbose = FALSE)
#> ℹ [2026-11-13 10:09:17] Skip `log1p()` because `layer = data` is not "counts"
pancreas_sub <- RunDorothea(
  pancreas_sub,
  layer = "counts",
  species = "Mus_musculus",
  method = "ulm",
  minsize = 5
)
#> ℹ [2026-11-13 10:09:24] Run "DoRothEA"/decoupleR with 12895 regulon edges
#> ℹ [2026-11-13 10:09:42] "DoRothEA" TF activity scores stored in assay "dorothea"
#> ℹ [2026-11-13 10:09:42] "DoRothEA" TF activity scores stored in <Seurat> metadata

DorotheaPlot(
  pancreas_sub,
  group.by = "CellType",
  group1 = "Endocrine",
  group2 = "Ductal",
  plot_type = "bar",
  top_n = 20
)
#> ℹ [2026-11-13 10:09:42] Compare "DoRothEA" TF activity: "Endocrine" vs "Ductal"


DorotheaPlot(
  pancreas_sub,
  group.by = "CellType",
  group1 = "Endocrine",
  group2 = "Ductal",
  plot_type = "lollipop",
  top_n = 20
)
#> ℹ [2026-11-13 10:09:42] Compare "DoRothEA" TF activity: "Endocrine" vs "Ductal"


DorotheaPlot(
  pancreas_sub,
  group.by = "CellType",
  group1 = "Endocrine",
  group2 = "Ductal",
  plot_type = "volcano"
)
#> ℹ [2026-11-13 10:09:43] Compare "DoRothEA" TF activity: "Endocrine" vs "Ductal"


tfs <- head(rownames(pancreas_sub[["dorothea"]]), 20)
GroupHeatmap(
  pancreas_sub,
  group.by = "CellType",
  features = tfs,
  assay = "dorothea",
  layer = "data",
  exp_method = "raw",
  lib_normalize = FALSE
)$plot


DorotheaPlot(
  pancreas_sub,
  group.by = "CellType",
  features = "Sox9",
  plot_type = "dim"
)
#> ℹ [2026-11-13 10:09:43] Draw "DoRothEA" embedding plots for 1 TFs


FeatureStatPlot(
  pancreas_sub,
  group.by = "CellType",
  stat.by = c("Sox9", "Neurod1", "Pdx1"),
  assay = "dorothea",
  layer = "data",
  plot_type = "violin"
)

FeatureDimPlot(
  pancreas_sub,
  features = c("Sox9", "Neurod1", "Pdx1"),
  assay = "dorothea",
  layer = "data"
)


DorotheaPlot(
  pancreas_sub,
  group.by = "CellType",
  group1 = "Endocrine",
  group2 = "Ductal",
  features = "Sox9",
  plot_type = "targets"
)
#> ! [2026-11-13 10:09:45] Dropping 3 "Sox9" targets missing from assay "RNA"
#> ℹ [2026-11-13 10:09:45] Draw "DoRothEA" regulon-target volcano for "Sox9" (10 targets)
```
