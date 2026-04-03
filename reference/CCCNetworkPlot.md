# CCC network and flow plots

CCC network and flow plots

## Usage

``` r
CCCNetworkPlot(
  srt,
  method = NULL,
  condition = NULL,
  dataset = 1,
  comparison = c(1, 2),
  plot_type = c("circle", "chord", "pathway", "individual_lr", "arrow", "sigmoid",
    "bipartite", "embedding_network", "diff_network"),
  display_by = c("aggregation", "interaction"),
  sender.use = NULL,
  receiver.use = NULL,
  ligand.use = NULL,
  receptor.use = NULL,
  interaction.use = NULL,
  group.by = NULL,
  reduction = NULL,
  dims = c(1, 2),
  signaling = NULL,
  pairLR.use = NULL,
  slot.name = "net",
  thresh = 0.05,
  measure = c("count", "weight"),
  value = "score",
  top_n = 20,
  ligand = NULL,
  receptor = NULL,
  reg.by = NULL,
  reg_palette = "Set1",
  reg_palcolor = NULL,
  expr.by = NULL,
  layout = c("circle", "hierarchy", "chord", "kk", "fr", "nicely", "lgl", "mds",
    "graphopt"),
  link_curvature = 0.2,
  link_alpha = 0.6,
  edge_value = c("sum", "mean", "max", "count"),
  edge_threshold = 0,
  edge_size = c(0.5, 1.8),
  edge_color = NULL,
  edge_alpha = 0.6,
  edge_line = c("curved", "straight"),
  edge_curvature = 0.2,
  directed = FALSE,
  arrow_type = "closed",
  arrow_angle = 20,
  arrow_length = grid::unit(0.02, "npc"),
  node_size = 5,
  node_alpha = 0.9,
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
  legend.title = NULL,
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list(),
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

  Plot type. One of `"circle"`, `"chord"`, `"pathway"`,
  `"individual_lr"`, `"arrow"`, `"sigmoid"`, `"bipartite"`,
  `"embedding_network"`, or `"diff_network"`.

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

- group.by:

  For `plot_type = "embedding_network"`: metadata column used to define
  cell groups. If `NULL`, the grouping stored in the CCC result is used
  when available.

- reduction:

  For `plot_type = "embedding_network"`: dimensional reduction to use.
  If `NULL`, the default reduction is used.

- dims:

  For `plot_type = "embedding_network"`: dimensions to plot.

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

- value:

  Value column or summary statistic to use.

- top_n:

  Number of top records to retain.

- ligand:

  For `plot_type = "bipartite"`: the ligand name to focus on. If `NULL`,
  the ligand with the highest total score is used.

- receptor:

  For `plot_type = "bipartite"`: optional receptor names to restrict to.
  If `NULL`, all receptors paired with `ligand` are shown.

- reg.by:

  For `plot_type = "bipartite"`: optional metadata column in `srt` used
  to color edges by regulation status (e.g. up/down). If `NULL`, edges
  are colored by sender cell type.

- reg_palette:

  For `plot_type = "bipartite"`: named character vector or palette name
  for regulation categories.

- reg_palcolor:

  For `plot_type = "bipartite"`: custom colors for regulation palette.

- expr.by:

  For `plot_type = "bipartite"`: optional metadata or score column used
  to scale edge line width. If `NULL`, all edges have equal width.

- layout:

  Layout used for graph-based network views. `"chord"` can also be
  requested via `plot_type = "circle"` for backward compatibility.

- link_curvature:

  Curvature used for circle-like differential links and flow edges.

- link_alpha:

  Alpha used for network edges.

- edge_value:

  Aggregation statistic for network edges.

- edge_threshold:

  Minimum edge value to keep.

- edge_size:

  Range used for scaling edge widths.

- edge_color:

  Optional edge color override. For differential networks, this may also
  be a length-2 vector for negative/positive changes.

- edge_alpha:

  Alpha used for embedding-network edges.

- edge_line:

  Edge geometry for `plot_type = "arrow"`, `"sigmoid"`, and
  `"embedding_network"`.

- edge_curvature:

  Curvature used for curved flow/embedding edges.

- directed:

  Whether to draw arrows for directed networks.

- arrow_type:

  Arrow head type passed to
  [`grid::arrow()`](https://rdrr.io/r/grid/arrow.html).

- arrow_angle:

  Arrow head angle passed to
  [`grid::arrow()`](https://rdrr.io/r/grid/arrow.html).

- arrow_length:

  Arrow length passed to
  [`grid::arrow()`](https://rdrr.io/r/grid/arrow.html).

- node_size:

  Base node size.

- node_alpha:

  Node alpha.

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

- legend.title:

  Legend title.

- font.size:

  Base font size.

- theme_use:

  Theme function used for styling.

- theme_args:

  Arguments passed to the theme function.

- verbose:

  Whether to print messages.

- ...:

  Additional plot-specific options. For chord plots, `reduce`,
  `max.groups`, `small.gap`, `big.gap`, and `lab.cex` can be used to
  adjust the CellChat-like chord layout.

## Value

A ggplot, patchwork, or recorded plot object.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-04-03 08:42:56] Start standard processing workflow...
#> ℹ [2026-04-03 08:42:57] Checking a list of <Seurat>...
#> ! [2026-04-03 08:42:57] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-03 08:42:57] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-03 08:42:58] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-03 08:42:59] Use the separate HVF from `srt_list`
#> ℹ [2026-04-03 08:42:59] Number of available HVF: 2000
#> ℹ [2026-04-03 08:42:59] Finished check
#> ℹ [2026-04-03 08:42:59] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-03 08:43:00] Perform pca linear dimension reduction
#> ℹ [2026-04-03 08:43:01] Use stored estimated dimensions 1:12 for Standardpca
#> ℹ [2026-04-03 08:43:01] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-03 08:43:01] Reorder clusters...
#> ℹ [2026-04-03 08:43:01] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-03 08:43:01] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-03 08:43:01] Perform umap nonlinear dimension reduction using Standardpca (1:12)
#> ℹ [2026-04-03 08:43:04] Perform umap nonlinear dimension reduction using Standardpca (1:12)
#> ✔ [2026-04-03 08:43:06] Standard processing workflow completed

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
#> ℹ [2026-04-03 08:43:06] Start CellChat analysis
#> ℹ [2026-04-03 08:44:34] Processing condition: "ConditionA"
#> Warning: The following arguments are not used: drop
#> [1] "Create a CellChat object from a data matrix"
#> Set cell identities for the new CellChat object 
#> The cell groups used for CellChat analysis are  Ductal, Ngn3-high-EP, Endocrine, Ngn3-low-EP, Pre-endocrine 
#> The number of highly variable ligand-receptor pairs used for signaling inference is 542 
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2026-04-03 08:44:35.886653]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2026-04-03 08:44:54.581426]"
#> ℹ [2026-04-03 08:44:54] Processing condition: "ConditionB"
#> Warning: The following arguments are not used: drop
#> [1] "Create a CellChat object from a data matrix"
#> Set cell identities for the new CellChat object 
#> The cell groups used for CellChat analysis are  Endocrine, Ngn3-high-EP, Ductal, Ngn3-low-EP, Pre-endocrine 
#> The number of highly variable ligand-receptor pairs used for signaling inference is 601 
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2026-04-03 08:44:55.749527]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2026-04-03 08:45:16.247063]"
#> ℹ [2026-04-03 08:45:16] Merging CellChat objects for comparison "ConditionA_vs_ConditionB"
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
#> ✔ [2026-04-03 08:45:16] CellChat analysis completed

CCCNetworkPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "circle",
  display_by = "aggregation",
  value = "count",
  top_n = 20
)



CCCNetworkPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "circle",
  display_by = "aggregation",
  value = "weight",
  top_n = 20
)



CCCNetworkPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "chord",
  display_by = "aggregation",
  top_n = 12
)



CCCNetworkPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "arrow",
  display_by = "interaction",
  sender.use = "Ductal",
  receiver.use = "Ngn3-low-EP",
  edge_line = "straight",
  directed = TRUE,
  top_n = 3
)


CCCNetworkPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "arrow",
  display_by = "interaction",
  top_n = 20
)


CCCNetworkPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "sigmoid",
  display_by = "interaction",
  top_n = 20
)


CCCNetworkPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "bipartite",
  display_by = "aggregation",
  top_n = 20
)


CCCNetworkPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "embedding_network",
  group.by = "CellType",
  reduction = "UMAP",
  top_n = 20,
  label = TRUE
)


CCCNetworkPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "pathway",
  signaling = "MK"
)



CCCNetworkPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA",
  plot_type = "individual_lr",
  signaling = "MK",
  pairLR.use = "MDK_SDC1"
)



CCCNetworkPlot(
  pancreas_sub,
  method = "CellChat",
  condition = "ConditionA_vs_ConditionB",
  plot_type = "diff_network",
  measure = "count"
)

```
