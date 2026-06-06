# Rename clusters for the Seurat object

Rename clusters for the Seurat object

## Usage

``` r
RenameClusters(
  srt,
  group.by,
  nameslist = list(),
  name = "newclusters",
  keep_levels = FALSE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- group.by:

  The old group used to rename cells.

- nameslist:

  A named list of new cluster value.

- name:

  The name of the new cluster stored in the Seurat object.

- keep_levels:

  If the old group is a factor, keep the order of the levels.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Examples

``` r
data(pancreas_sub)

# Rename all clusters
pancreas_sub <- standard_scop(pancreas_sub)
levels(pancreas_sub@meta.data[["SubCellType"]]) <- unique(
  pancreas_sub@meta.data[["SubCellType"]]
)
pancreas_sub <- RenameClusters(
  pancreas_sub,
  group.by = "SubCellType",
  nameslist = letters[1:8]
)
CellDimPlot(pancreas_sub, "newclusters")

# Rename specified clusters
pancreas_sub <- RenameClusters(pancreas_sub,
  group.by = "SubCellType",
  nameslist = list("a" = "Alpha", "b" = "Beta")
)
CellDimPlot(pancreas_sub, "newclusters")

# Merge and rename clusters
pancreas_sub <- RenameClusters(
  pancreas_sub,
  group.by = "SubCellType",
  nameslist = list(
    "EndocrineClusters" = c("Alpha", "Beta", "Epsilon", "Delta")
  ),
  name = "Merged",
  keep_levels = TRUE
)
CellDimPlot(pancreas_sub, "Merged")
```
