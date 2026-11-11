# Plot SCENIC+ eRegulon activity

[`SCENICPlot()`](https://mengxu98.github.io/scop/reference/SCENICPlot.md)
with `tool_name = "SCENICPlus"` and `assay = "scenicplus"`. Default
`"heatmap_dotplot"`: color = TF expression, size = RSS.

## Usage

``` r
SCENICPlusPlot(
  object,
  group.by,
  tool_name = "SCENICPlus",
  assay = "scenicplus",
  plot_type = c("heatmap_dotplot", "rss_rank", "rss_heatmap", "rss_dotplot",
    "eregulon_dim", "activity_cor_dumbbell", "regulon_size", "network_graph", "network",
    "egrn", "overlap", "target_bar", "coverage"),
  srt = NULL,
  ...
)
```

## Arguments

- object:

  A Seurat object from
  [`RunSCENICPlus()`](https://mengxu98.github.io/scop/reference/RunSCENICPlus.md).

- group.by:

  Metadata column for cell groups.

- tool_name:

  Tools slot. Default `"SCENICPlus"`.

- assay:

  ERegulon AUC assay. Default `"scenicplus"`.

- plot_type:

  Plot type. Default `"heatmap_dotplot"`.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

- ...:

  Arguments passed to other functions.

## Value

Same as
[`SCENICPlot()`](https://mengxu98.github.io/scop/reference/SCENICPlot.md).

## See also

[RunSCENICPlus](https://mengxu98.github.io/scop/reference/RunSCENICPlus.md),
[SCENICPlot](https://mengxu98.github.io/scop/reference/SCENICPlot.md),
[GroupHeatmap](https://mengxu98.github.io/scop/reference/GroupHeatmap.md),
[FeatureStatPlot](https://mengxu98.github.io/scop/reference/FeatureStatPlot.md),
[FeatureDimPlot](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
pancreas_sub <- RunSCENICPlus(
  pancreas_sub,
  backend = "python",
  python_result_dir = "test/scenicplus"
)
scenicplus_dot <- SCENICPlusPlot(pancreas_sub, group.by = "CellType")
example_tfs <- unique(scenicplus_dot$top_table$TF)[1:3]
regulons <- head(rownames(pancreas_sub[["scenicplus"]]), 6)
GroupHeatmap(
  pancreas_sub, group.by = "CellType", features = regulons,
  assay = "scenicplus", layer = "data", exp_method = "raw",
  lib_normalize = FALSE
)$plot
FeatureStatPlot(
  pancreas_sub, group.by = "CellType", stat.by = regulons,
  assay = "scenicplus", layer = "data", plot_type = "violin"
)
FeatureDimPlot(
  pancreas_sub, features = regulons, assay = "scenicplus",
  layer = "data", reduction = "UMAP"
)
SCENICPlusPlot(pancreas_sub, group.by = "CellType", plot_type = "egrn", features = example_tfs)
SCENICPlusPlot(
  pancreas_sub,
  group.by = "CellType", plot_type = "network_graph", features = example_tfs
)
} # }
```
