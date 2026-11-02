# Plot ATAC reduction-depth correlations

Compute and plot correlations between peak-count totals or explicit
depth metadata and reduced-dimension components, aligning cells
explicitly by name.

## Usage

``` r
ATACDepthCorPlot(object, assay = NULL, reduction = "lsi", n = 10L,
  depth.column = NULL, method = "pearson", cells = NULL,
  return.data = FALSE)
```

## Arguments

- object:

  A Seurat object containing a ChromatinAssay.

- assay:

  Name of the ChromatinAssay; defaults to the active assay.

- reduction:

  Name of a dimensional reduction, normally lsi.

- n:

  Maximum number of leading components; NULL uses all components.

- depth.column:

  Optional numeric non-negative metadata column. NULL recomputes totals
  from the selected assay counts layer, avoiding stale nCount metadata.
  Peak-matrix totals are not necessarily total fragments or sequencing
  reads; select a verified fragment count column explicitly for that
  depth measure.

- method:

  Correlation method: pearson, spearman, or kendall.

- cells:

  Cell names in the selected assay; NULL uses all assay cells.

- return.data:

  Return a correlation table instead of a plot.

## Value

A ggplot, or a data.frame with component, correlation and n.cells when
return.data is TRUE. Constant components have NA correlations. No data
are written into the object and no component is automatically excluded.

## Details

All selected assay cells must occur in the reduction and have finite
depth and embedding values; missing cells are not silently dropped.
Depth must vary. Component numbers are signed correlations; SVD signs
are arbitrary. This diagnostic supports, but does not by itself
determine, LSI selection.

## See also

[`RunATACQC`](https://mengxu98.github.io/scop/reference/RunATACQC.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ATACDepthCorPlot(object, assay = "peaks", reduction = "lsi")
ATACDepthCorPlot(object, assay = "peaks", return.data = TRUE)
} # }
```
