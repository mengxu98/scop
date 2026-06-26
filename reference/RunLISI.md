# Compute LISI scores on a Seurat object

Compute per-cell Local Inverse Simpson's Index (LISI) scores from a
dimensional reduction and store them in the `meta.data` and `tools`
slots of a `Seurat` object.

## Usage

``` r
RunLISI(
  srt,
  reductions = NULL,
  reduction = NULL,
  dims = NULL,
  label_colnames = NULL,
  prefix = NULL,
  tool_name = NULL,
  perplexity = 30,
  tol = 1e-05,
  max_iter = 50,
  overwrite = TRUE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A `Seurat` object.

- reductions:

  Character vector of dimensional reductions used to compute LISI. If
  `NULL`,
  [`DefaultReduction()`](https://mengxu98.github.io/scop/reference/DefaultReduction.md)
  is used.

- reduction:

  Deprecated alias of `reductions`.

- dims:

  Dimensions to use from the reduction. Default is `NULL`, which uses
  all available dimensions.

- label_colnames:

  Character vector of metadata columns used for LISI. If `NULL`,
  `RunLISI()` will try to use `srt@misc[["integration_batch"]]`.

- prefix:

  Prefix used for the stored LISI metadata columns. If `NULL`, the
  reduction names are used.

- tool_name:

  Name used to store detailed results in `srt@tools`. Default is
  `"LISI"` when multiple reductions are provided, otherwise
  `paste0(prefix, "_LISI")`.

- perplexity:

  Effective neighborhood size. Default is `30`.

- tol:

  Tolerance used in the binary search for the target perplexity. Default
  is `1e-5`.

- max_iter:

  Maximum number of binary-search iterations. Default is `50`.

- overwrite:

  Whether to overwrite existing metadata columns. Default is `TRUE`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A modified `Seurat` object.

## See also

[thisutils::compute_lisi](https://mengxu98.github.io/thisutils/reference/compute_lisi.html),
[LISIPlot](https://mengxu98.github.io/scop/reference/LISIPlot.md)

## Examples

``` r
data(panc8_sub)
panc8_sub <- integration_scop(
  panc8_sub,
  batch = "tech",
  integration_method = "Harmony5"
)
#> ◌ [2026-06-26 12:03:41] Run integration workflow...
#> Warning: No layers found matching search pattern provided
#> ℹ [2026-06-26 12:03:42] Perform `Seurat::NormalizeData()` on split layers for Seurat v5 integration
#> Error: NormalizeData.Seurat requires one counts layer named 'counts'.
names(panc8_sub@reductions)
#> NULL

panc8_sub <- RunLISI(
  panc8_sub,
  reductions = c("pcaUMAP2D", "Harmony5UMAP2D")
)
#> Error in RunLISI(panc8_sub, reductions = c("pcaUMAP2D", "Harmony5UMAP2D")): Reductions not found in <Seurat>: "pcaUMAP2D" and "Harmony5UMAP2D"
LISIPlot(
  panc8_sub,
  combine = TRUE
)
#> Error in benchmark_feature_plot(srt = srt, features = features, tool_name = tool_name,     reduction = reduction, plot_type = plot_type, plot_boxplot = plot_boxplot,     boxplot_jitter = boxplot_jitter, combine = combine, nrow = nrow,     ncol = ncol, byrow = byrow, pt.size = pt.size, pt.alpha = pt.alpha,     palette = palette, palcolor = palcolor, theme_use = theme_use,     theme_args = theme_args, verbose = verbose, ...): No per-cell benchmark columns found. Please provide `features` or a
#> valid `tool_name`.
```
