# Plot SCENIC+ eRegulon activity

[`SCENICPlot()`](https://mengxu98.github.io/scop/reference/SCENICPlot.md)
with `tool_name = "SCENICPlus"` and `assay = "scenicplus"`. Default
`"heatmap_dotplot"`: color = TF expression, size = RSS.

## Usage

``` r
SCENICPlusPlot(
  srt,
  group.by,
  tool_name = "SCENICPlus",
  assay = "scenicplus",
  plot_type = c("heatmap_dotplot", "rss_rank", "rss_heatmap", "rss_dotplot",
    "activity_heatmap", "activity_violin", "activity_dim", "eregulon_dim",
    "activity_cor_dumbbell", "regulon_size", "network_graph", "network", "egrn",
    "overlap", "target_bar", "coverage"),
  ...
)
```

## Arguments

- srt:

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

- ...:

  Passed to
  [`SCENICPlot()`](https://mengxu98.github.io/scop/reference/SCENICPlot.md).

## Value

Same as
[`SCENICPlot()`](https://mengxu98.github.io/scop/reference/SCENICPlot.md).

## See also

[RunSCENICPlus](https://mengxu98.github.io/scop/reference/RunSCENICPlus.md),
[SCENICPlot](https://mengxu98.github.io/scop/reference/SCENICPlot.md)

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
SCENICPlusPlot(pancreas_sub, group.by = "CellType", plot_type = "egrn", features = example_tfs)
SCENICPlusPlot(pancreas_sub, group.by = "CellType", plot_type = "network_graph", features = example_tfs)
} # }
```
