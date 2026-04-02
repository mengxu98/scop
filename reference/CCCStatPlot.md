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
  plot_type = c("bar", "sankey", "box", "violin", "comparison", "lr_contribution",
    "gene", "ranknet", "scatter", "role_change"),
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
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  verbose = TRUE,
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

  - `"bar"` — horizontal bar chart of top interactions by aggregated
    score.

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

  Legend position.

- legend.direction:

  Legend direction.

- font.size:

  Base font size.

- theme_use:

  Theme function used for styling.

- theme_args:

  Arguments passed to the theme function.

- combine:

  Whether to combine multiple panels.

- nrow:

  Number of rows in combined layout.

- ncol:

  Number of columns in combined layout.

- verbose:

  Whether to print messages.

- ...:

  Additional plot-specific options.

## Value

A ggplot or recorded base plot object.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-04-02 15:24:49] Start standard processing workflow...
#> ℹ [2026-04-02 15:24:50] Checking a list of <Seurat>...
#> ! [2026-04-02 15:24:50] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-02 15:24:50] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 15:24:51] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 15:24:52] Use the separate HVF from `srt_list`
#> ℹ [2026-04-02 15:24:52] Number of available HVF: 2000
#> ℹ [2026-04-02 15:24:52] Finished check
#> ℹ [2026-04-02 15:24:52] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-02 15:24:53] Perform pca linear dimension reduction
#> ℹ [2026-04-02 15:24:57] Use stored estimated dimensions 1:50 for Standardpca
#> Warning: Caught FutureLaunchError. Canceling all iterations ...
#> ! [2026-04-02 15:24:57] <FutureLaunchError: Caught an unexpected error of class FutureLaunchError when trying to launch future (‘future_lapply-1’) on backend of class SequentialFutureBackend. The reason was: future::evalFuture() failed on runnervmrg6be (pid 85355) at 2026-04-02T15:24:57. Using package 'future' v1.70.0. Possible other reasons: Failed to attach one or more future-backend packages: there is no package called ‘future’ [future <unnamed>; on 4a75d434f7a9a2903adedbeee3372830@runnervmrg6be<85355>] [future ‘future_lapply-1’ (4a75d434f7a9a2903adedbeee3372830-3); on 4a75d434f7a9a2903adedbeee3372830@runnervmrg6be<85355>]>
#> !                       
#> !                       Occurred on: 4a75d434f7a9a2903adedbeee3372830 [runnervmrg6be; pid 85355]
#> !                       Future: 4a75d434f7a9a2903adedbeee3372830-3 (‘future_lapply-1’)
#> !                       
#> !                       DEBUG: BEGIN TROUBLESHOOTING HELP
#> !                       SequentialFuture:
#> !                       Label: ‘future_lapply-1’
#> !                       Expression:
#> Error in glue(str, .envir = .envir, .transformer = transformer, .cli = TRUE,     .trim = .trim): Expecting '}'

pc1 <- Seurat::Embeddings(pancreas_sub, "Standardpca")[, 1]
#> Error in object[[reduction]]: ‘Standardpca’ not found in this Seurat object
#>  
ct <- as.character(pancreas_sub$CellType)
ct_medians <- tapply(pc1, ct, median)
#> Error: object 'pc1' not found
pancreas_sub$Condition <- ifelse(
  pc1 > ct_medians[ct],
  "ConditionA",
  "ConditionB"
)
#> Error: object 'pc1' not found

pancreas_sub <- RunCellChat(
  pancreas_sub,
  group.by = "CellType",
  group_column = "Condition",
  group_cmp = list(c("ConditionA", "ConditionB")),
  species = "Mus_musculus"
)
#> ℹ [2026-04-02 15:24:57] Start CellChat analysis
#> Error in RunCellChat(pancreas_sub, group.by = "CellType", group_column = "Condition",     group_cmp = list(c("ConditionA", "ConditionB")), species = "Mus_musculus"): "Condition" does not exist in <Seurat>

CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "sankey",
  display_by = "aggregation",
  top_n = 20
)
#> Error in get_dataset_object(srt, condition = condition, dataset = dataset): Unable to determine which CellChat object to plot. Please specify
#> `condition`

CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "sankey",
  display_by = "interaction",
  top_n = 20
)
#> Error in get_dataset_object(srt, condition = condition, dataset = dataset): Unable to determine which CellChat object to plot. Please specify
#> `condition`

CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "box",
  facet_by = "sender",
  top_n = 200
)
#> Error in get_dataset_object(srt, condition = condition, dataset = dataset): Unable to determine which CellChat object to plot. Please specify
#> `condition`

CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "violin",
  facet_by = "receiver",
  top_n = 200
)
#> Error in get_dataset_object(srt, condition = condition, dataset = dataset): Unable to determine which CellChat object to plot. Please specify
#> `condition`

CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "bar",
  palette = "Paired",
  top_n = 100
)
#> Error in get_dataset_object(srt, condition = condition, dataset = dataset): Unable to determine which CellChat object to plot. Please specify
#> `condition`

CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "scatter"
)
#> Error in get_dataset_object(srt, condition = condition, dataset = dataset): Unable to determine which CellChat object to plot. Please specify
#> `condition`

CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "lr_contribution",
  signaling = "MK"
)
#> Error in get_dataset_object(srt, condition = condition, dataset = dataset): Unable to determine which CellChat object to plot. Please specify
#> `condition`

CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "gene",
  signaling = "MK"
)
#> Error in get_dataset_object(srt, condition = condition, dataset = dataset): Unable to determine which CellChat object to plot. Please specify
#> `condition`

CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA_vs_ConditionB",
  plot_type = "comparison",
  measure = "count",
  compare_by = "overall"
)
#> Error in get_dataset_object(srt, condition = condition, dataset = dataset): Unable to determine which CellChat object to plot. Please specify
#> `condition`

CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA_vs_ConditionB",
  plot_type = "comparison",
  measure = "weight",
  compare_by = "celltype",
  pattern = "all"
)
#> Error in get_dataset_object(srt, condition = condition, dataset = dataset): Unable to determine which CellChat object to plot. Please specify
#> `condition`

CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA_vs_ConditionB",
  plot_type = "ranknet"
)
#> Error in .cc_get_cmp(srt = srt, condition = condition): Comparison "ConditionA_vs_ConditionB" not found in CellChat results

CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA_vs_ConditionB",
  idents.use = "Ductal",
  plot_type = "role_change"
)
#> Error in .cc_get_cmp(srt = srt, condition = condition): Comparison "ConditionA_vs_ConditionB" not found in CellChat results
```
