# CCC statistical distribution and summary plots

CCC statistical distribution and summary plots

## Usage

``` r
CCCStatPlot(
  srt,
  method = NULL,
  condition = NULL,
  dataset = 1,
  comparison = c(1, 2),
  plot_type = c("bar", "sankey", "box", "violin", "role_scatter", "role_network",
    "role_network_marsilea", "pathway_summary", "comparison", "lr_contribution", "gene",
    "ranknet", "scatter", "role_change"),
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
  pattern = c("all", "outgoing", "incoming"),
  compare_by = c("overall", "celltype"),
  value = "score",
  stat_type = c("score", "count"),
  top_n = 20,
  x_text_angle = 90,
  min_receiver_flow = 0,
  link_alpha = 0.6,
  facet_by = NULL,
  edge_value = c("sum", "mean", "max", "count"),
  edge_threshold = 0,
  palette = "Chinese",
  palcolor = NULL,
  cell_palette = NULL,
  cell_palcolor = NULL,
  link_palette = NULL,
  link_palcolor = NULL,
  title = NULL,
  subtitle = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list(),
  grid_major = TRUE,
  grid_major_colour = "grey80",
  grid_major_linetype = 2,
  grid_major_linewidth = 0.3,
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
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

  Plot type. One of:

  - `"bar"` — horizontal bar chart of top pairs or interactions
    according to `display_by`.

  - `"sankey"` — alluvial/sankey flow diagram.

  - `"box"` / `"violin"` — distribution of interaction scores across
    sender-receiver pairs.

  - `"comparison"` — comparison bars at overall or celltype level.

  - `"lr_contribution"` — ligand-receptor contribution bar plot.

  - `"gene"` — pathway-related ligand/receptor gene expression panel.

  - `"ranknet"` — pathway ranking comparison plot.

  - `"scatter"` — outgoing vs. incoming signaling strength scatter.

  - `"role_change"` — signaling change scatter for one cell identity.

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

- compare_by:

  Comparison mode for CellChat summary plots.

- value:

  Value column or summary statistic to use.

- stat_type:

  For `"bar"`: what to summarize per interaction. One of `"score"`
  (total aggregated score) or `"count"` (number of significant
  interactions).

- top_n:

  Number of top records to retain.

- x_text_angle:

  Rotation angle for x-axis labels.

- min_receiver_flow:

  For `"sankey"`: minimum total receiver-side flow retained after top-N
  ranking. Useful when many small receiver nodes make the right side
  unreadable.

- link_alpha:

  Alpha used for network edges.

- facet_by:

  Faceting variable for interaction-level plots.

- edge_value:

  Aggregation statistic for network edges.

- edge_threshold:

  Minimum edge value to keep.

- palette:

  Main palette name.

- palcolor:

  Main custom palette colors.

- cell_palette:

  Cell annotation palette name.

- cell_palcolor:

  Custom cell annotation colors.

- link_palette:

  Link palette name.

- link_palcolor:

  Custom link palette colors.

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

- grid_major:

  Whether to show major panel grid lines for applicable statistical
  panels.

- grid_major_colour:

  Color of major panel grid lines.

- grid_major_linetype:

  Linetype of major panel grid lines.

- grid_major_linewidth:

  Line width of major panel grid lines.

- combine:

  Combine plots with
  [patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html).
  `combine = FALSE` returns a list of ggplots.

- nrow:

  Combine plots with
  [patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html).
  `combine = FALSE` returns a list of ggplots.

- ncol:

  Number of columns of the combined plot.

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

A ggplot or recorded base plot object.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 20:59:05] Start standard processing workflow...
#> ℹ [2026-09-06 20:59:06] Checking a list of <Seurat>...
#> ! [2026-09-06 20:59:06] Data 1/1 of the `srt_list` is "unknown"
#> Warning: Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 20:59:06] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 20:59:06] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 20:59:07] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 20:59:07] Number of available HVF: 2000
#> ℹ [2026-09-06 20:59:07] Finished check
#> ℹ [2026-09-06 20:59:07] Perform `ScaleData()`
#> ℹ [2026-09-06 20:59:07] Perform pca linear dimension reduction
#> ℹ [2026-09-06 20:59:07] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 20:59:07] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 20:59:07] Reorder clusters...
#> ℹ [2026-09-06 20:59:07] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 20:59:07] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 20:59:11] Standard processing workflow completed

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
#> ℹ [2026-09-06 20:59:11] Start CellChat analysis
#> ℹ [2026-09-06 20:59:11] Processing condition: "ConditionA"
#> [1] "Create a CellChat object from a data matrix"
#> Set cell identities for the new CellChat object 
#> The cell groups used for CellChat analysis are  Ductal, Endocrine, Ngn3-high-EP, Ngn3-low-EP, Pre-endocrine 
#> ! [2026-09-06 20:59:12] Function "CellChatDB.mouse" not found in CellChat namespace
#> Warning: Function "CellChatDB.mouse" not found in CellChat namespace
#> The number of highly variable ligand-receptor pairs used for signaling inference is 542 
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2026-09-06 20:59:13.725163]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2026-09-06 20:59:33.19493]"
#> ℹ [2026-09-06 20:59:33] Processing condition: "ConditionB"
#> [1] "Create a CellChat object from a data matrix"
#> Set cell identities for the new CellChat object 
#> The cell groups used for CellChat analysis are  Ductal, Endocrine, Ngn3-high-EP, Ngn3-low-EP, Pre-endocrine 
#> ! [2026-09-06 20:59:33] Function "CellChatDB.mouse" not found in CellChat namespace
#> Warning: Function "CellChatDB.mouse" not found in CellChat namespace
#> The number of highly variable ligand-receptor pairs used for signaling inference is 597 
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2026-09-06 20:59:35.201129]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2026-09-06 20:59:57.554552]"
#> ℹ [2026-09-06 20:59:57] Merging CellChat objects for comparison "ConditionA_vs_ConditionB"
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
#> ✔ [2026-09-06 20:59:58] CellChat analysis completed

CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "sankey",
  display_by = "aggregation",
  top_n = 20
)
#> ! [2026-09-06 20:59:58] `thisplot::StatPlot()` sankey is count-based. For `CCCStatPlot()` with `plot_type = 'sankey'`, `edge_value` is used to rank/filter pairs, but flow width is shown by interaction count.
#> Warning: `thisplot::StatPlot()` sankey is count-based. For `CCCStatPlot()` with `plot_type = 'sankey'`, `edge_value` is used to rank/filter pairs, but flow width is shown by interaction count.


CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "sankey",
  display_by = "interaction",
  top_n = 20
)


CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "box",
  facet_by = "sender",
  top_n = 200
)


CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "violin",
  facet_by = "receiver",
  top_n = 200
)
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.


CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "bar",
  palette = "Paired",
  top_n = 100
)


CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "scatter"
)


CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "lr_contribution",
  signaling = "MK"
)


CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "gene",
  signaling = "MK"
)
#> ℹ [2026-09-06 21:00:02] Setting `group.by` to "Features" as `plot.by` is set to "feature"


CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA_vs_ConditionB",
  plot_type = "comparison",
  measure = "count",
  compare_by = "overall"
)


CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA_vs_ConditionB",
  plot_type = "comparison",
  measure = "weight",
  compare_by = "celltype",
  pattern = "all"
)


CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA_vs_ConditionB",
  plot_type = "ranknet"
)


CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA_vs_ConditionB",
  idents.use = "Ductal",
  plot_type = "role_change"
)
```
