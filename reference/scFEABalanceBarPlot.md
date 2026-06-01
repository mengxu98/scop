# Plot scFEA metabolite balance changes

Computes scFEA balance mean differences between two Seurat metadata
groups and draws a metabolite bar plot.

## Usage

``` r
scFEABalanceBarPlot(
  srt,
  group.by,
  ident.1 = NULL,
  ident.2 = NULL,
  assay = "scFEAbalance",
  layer = "data",
  top_n = NULL,
  p_adj_cutoff = 0.01,
  title = NULL
)
```

## Arguments

- srt:

  A Seurat object returned by
  [`RunscFEA()`](https://mengxu98.github.io/scop/reference/RunscFEA.md).

- group.by:

  Metadata column defining groups. If `ident.1` and `ident.2` are both
  `NULL`, each group is compared against all remaining cells.

- ident.1, ident.2:

  Group names to compare. Difference is `mean(ident.1) - mean(ident.2)`.
  If both are `NULL`, one-vs-rest contrasts are generated automatically.

- assay:

  Balance assay name.

- layer:

  Balance assay layer.

- top_n:

  If `NULL`, plot all metabolites. If positive, plot the top `top_n`
  increased and top `top_n` decreased metabolites after the
  adjusted-p-value filter.

- p_adj_cutoff:

  Adjusted p-value cutoff used when `top_n` is not `NULL`.

- title:

  Optional plot title.

## Value

For a single contrast, a list with `plot` and plotted `data`. For
automatic one-vs-rest mode, a named list of single-contrast results is
returned and combined statistics are stored in the outer `"data"`
attribute.
