# CCC heatmap and dot matrix plot

CCC heatmap and dot matrix plot

## Usage

``` r
CCCHeatmap(
  srt,
  method = NULL,
  condition = NULL,
  dataset = 1,
  comparison = c(1, 2),
  plot_type = c("heatmap", "focused_heatmap", "dot", "matrix_dot", "tile",
    "source_target_dot", "source_target_tile", "sample_dot", "bubble", "bubble_lr",
    "pathway_bubble", "ligand_target", "role_heatmap", "role_network",
    "role_network_marsilea", "diff_heatmap"),
  display_by = c("aggregation", "interaction"),
  sender.use = NULL,
  receiver.use = NULL,
  ligand.use = NULL,
  receptor.use = NULL,
  interaction.use = NULL,
  signaling = NULL,
  pairLR.use = NULL,
  slot.name = "net",
  thresh = 0.05,
  measure = c("count", "weight"),
  pattern = c("outgoing", "incoming", "all"),
  value = "sum",
  top_n = 500,
  top_anno = "bar",
  right_anno = "cell",
  left_anno = "bar",
  bottom_anno = "cell",
  bar_value = "sum",
  add_text = NULL,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  color.by = c("score", "pvalue"),
  x_text_angle = 90,
  facet_by = NULL,
  show_row_names = TRUE,
  show_column_names = TRUE,
  edge_value = c("sum", "mean", "max", "count"),
  border = TRUE,
  value_palette = "RdBu",
  value_palcolor = NULL,
  cell_palette = "Chinese",
  cell_palcolor = NULL,
  palette = "Chinese",
  palcolor = NULL,
  width = NULL,
  height = NULL,
  units = "inch",
  title = NULL,
  subtitle = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list(),
  verbose = TRUE,
  combine_methods = c("separate", "support", "rank", "legacy"),
  resource = NULL,
  sample = NULL,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object.

- method:

  Communication result type to use.

- condition:

  Result name or comparison name.

- dataset:

  Dataset index or name.

- comparison:

  Comparison indices or names.

- plot_type:

  Plot type. One of `"heatmap"` or `"dot"`. `"bubble"` is a
  CellChat-specific interaction bubble matrix. `"ligand_target"` is a
  special heatmap path available only with `NicheNet`/`MultiNicheNet`
  results. `"role_heatmap"` and `"diff_heatmap"` are CellChat-specific
  pathway role views.

- display_by:

  Whether to summarize by `"aggregation"` or `"interaction"`.

- sender.use:

  Sender cell types to keep.

- receiver.use:

  Receiver cell types to keep.

- ligand.use:

  Ligands to keep.

- receptor.use:

  Receptors to keep.

- interaction.use:

  Interaction names to keep.

- signaling:

  Signaling pathway to focus on.

- pairLR.use:

  Specific ligand-receptor pair(s) to keep.

- slot.name:

  CellChat slot name.

- thresh:

  Significance threshold used when extracting communication results.

- measure:

  Summary measure for CellChat objects.

- pattern:

  Pattern used for pathway role plots.

- value:

  Value column or summary statistic to use.

- top_n:

  Number of top records to retain.

- top_anno, bottom_anno:

  Column-side annotations for sender groups. Each side accepts `NULL`,
  `"bar"`, `"box"`, `"point"`, `"line"`, `"histogram"`, `"density"`,
  `"violin"`, `"cell"`, or a character vector containing multiple
  values. Defaults are `top_anno = "bar"` and `bottom_anno = "cell"`.

- left_anno, right_anno:

  Row-side annotations for receiver groups. Each side accepts `NULL`,
  `"bar"`, `"box"`, `"point"`, `"line"`, `"histogram"`, `"density"`,
  `"violin"`, `"cell"`, or a character vector containing multiple
  values. Defaults are `left_anno = "bar"` and `right_anno = "cell"`.

- bar_value:

  Aggregation metric shown in the bar annotations. One or more of
  `"count"`, `"sum"`, `"mean"`, or `"max"`. Multiple values add multiple
  annotation tracks. Default `"sum"`.

- add_text:

  Logical. Show numeric value labels inside each cell (heatmap mode
  only). Default `TRUE` for aggregation mode, `FALSE` for interaction
  mode.

- cluster_rows, cluster_columns:

  Whether to cluster heatmap rows/columns. Defaults are both `FALSE`.

- color.by:

  For interaction heatmaps, value used for tile coloring. Usually
  `"score"` or `"pvalue"`.

- x_text_angle:

  Rotation angle for x-axis labels.

- facet_by:

  Faceting variable for interaction-level plots.

- show_row_names, show_column_names:

  Whether to draw row/column names for the heatmap body.

- edge_value:

  Aggregation statistic for network edges.

- border:

  Logical. Whether to draw borders for the heatmap body and all
  annotation tracks. Default `TRUE`.

- value_palette:

  Palette used for heatmap value fills.

- value_palcolor:

  Optional custom colors for `value_palette`.

- cell_palette:

  Cell annotation palette name.

- cell_palcolor:

  Custom cell annotation colors.

- palette:

  Main palette name.

- palcolor:

  Main custom palette colors.

- width, height:

  Optional heatmap body width and height. When only one is supplied, the
  other is inferred from the matrix dimensions to keep cells square.
  When both are `NULL`, a square-cell size is computed automatically.

- units:

  Units for `width` and `height`. Default `"inch"`.

- title:

  Plot title.

- subtitle:

  Plot subtitle.

- legend.position:

  Legend placement (`"none"`, `"left"`, `"right"`, `"bottom"`, `"top"`),
  direction, and title. `legend.title = NULL` uses the group name.

- legend.direction:

  Legend direction: `"horizontal"` or `"vertical"`.

- font.size:

  Base font size.

- theme_use, theme_args:

  Theme name or function, plus extra theme arguments.

- verbose:

  Whether to print messages.

- combine_methods:

  Behavior when `method = "CCC"`. `"separate"` returns one panel per
  backend, `"support"` counts supporting backends, `"rank"` combines
  within-method percentile ranks for visualization, and `"legacy"`
  retains the deprecated raw-score aggregation.

- resource, sample:

  Optional resource and sample/context filters for unified CCC results.

- ...:

  Additional plot-specific options.

## Value

A ggplot / patchwork object wrapping the ComplexHeatmap grob.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 20:52:08] Start standard processing workflow...
#> ℹ [2026-09-06 20:52:09] Checking a list of <Seurat>...
#> ! [2026-09-06 20:52:09] Data 1/1 of the `srt_list` is "unknown"
#> Warning: Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 20:52:09] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 20:52:09] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 20:52:09] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 20:52:09] Number of available HVF: 2000
#> ℹ [2026-09-06 20:52:09] Finished check
#> ℹ [2026-09-06 20:52:09] Perform `ScaleData()`
#> ℹ [2026-09-06 20:52:10] Perform pca linear dimension reduction
#> ℹ [2026-09-06 20:52:10] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 20:52:10] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 20:52:10] Reorder clusters...
#> ℹ [2026-09-06 20:52:10] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 20:52:10] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 20:52:14] Standard processing workflow completed

pc1 <- Seurat::Embeddings(pancreas_sub, "Standardpca")[, 1]
ct <- as.character(pancreas_sub$CellType)
ct_medians <- tapply(pc1, ct, median)
pancreas_sub$Condition <- ifelse(
  pc1 > ct_medians[ct],
  "ConditionA",
  "ConditionB"
)

pancreas_sub <- RunCellChat(
  pancreas_sub,
  group.by = "CellType",
  group_column = "Condition",
  group_cmp = list(c("ConditionA", "ConditionB")),
  species = "Mus_musculus"
)
#> ℹ [2026-09-06 20:52:14] Start CellChat analysis
#> ℹ [2026-09-06 20:57:07] Processing condition: "ConditionA"
#> [1] "Create a CellChat object from a data matrix"
#> Set cell identities for the new CellChat object 
#> The cell groups used for CellChat analysis are  Ductal, Endocrine, Ngn3-high-EP, Ngn3-low-EP, Pre-endocrine 
#> ! [2026-09-06 20:57:07] Function "CellChatDB.mouse" not found in CellChat namespace
#> Warning: Function "CellChatDB.mouse" not found in CellChat namespace
#> The number of highly variable ligand-receptor pairs used for signaling inference is 542 
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2026-09-06 20:57:09.36946]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2026-09-06 20:57:29.779295]"
#> ℹ [2026-09-06 20:57:29] Processing condition: "ConditionB"
#> [1] "Create a CellChat object from a data matrix"
#> Set cell identities for the new CellChat object 
#> The cell groups used for CellChat analysis are  Ductal, Endocrine, Ngn3-high-EP, Ngn3-low-EP, Pre-endocrine 
#> ! [2026-09-06 20:57:30] Function "CellChatDB.mouse" not found in CellChat namespace
#> Warning: Function "CellChatDB.mouse" not found in CellChat namespace
#> The number of highly variable ligand-receptor pairs used for signaling inference is 597 
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2026-09-06 20:57:31.808889]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2026-09-06 20:57:54.49997]"
#> ℹ [2026-09-06 20:57:54] Merging CellChat objects for comparison "ConditionA_vs_ConditionB"
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
#> ✔ [2026-09-06 20:57:55] CellChat analysis completed

CCCHeatmap(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "heatmap",
  display_by = "aggregation",
  top_n = 20
)


CCCHeatmap(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "dot",
  display_by = "interaction",
  sender.use = "Ductal",
  top_n = 20
)


CCCHeatmap(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "dot",
  display_by = "interaction",
  receiver.use = "Ngn3-low-EP",
  top_n = 20
)


CCCHeatmap(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "heatmap",
  display_by = "aggregation",
  top_n = 20
)


CCCHeatmap(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "heatmap",
  display_by = "interaction",
  facet_by = "sender",
  top_n = 10
)


CCCHeatmap(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "heatmap",
  display_by = "interaction",
  facet_by = "sender",
  color.by = "pvalue",
  top_n = 10
)


CCCHeatmap(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "bubble",
  top_n = 5
)


CCCHeatmap(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA_vs_ConditionB",
  plot_type = "bubble",
  top_n = 5
)


CCCHeatmap(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "role_heatmap",
  #' pattern = "outgoing",
  width = 0.6,
  height = 2.5
)


CCCHeatmap(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "role_heatmap",
  pattern = "outgoing",
  width = 0.6,
  height = 2.5
)


CCCHeatmap(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA_vs_ConditionB",
  plot_type = "role_heatmap",
  pattern = "incoming",
  palette = "Paired",
  width = 0.6,
  height = 3.5
)


CCCHeatmap(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA_vs_ConditionB",
  plot_type = "diff_heatmap",
  pattern = "all",
  palette = "Paired",
  top_n = 20,
  width = 0.6,
  height = 3.5
)
```
