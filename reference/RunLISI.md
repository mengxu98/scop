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
  nn_method = c("auto", "exact"),
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

- nn_method:

  Nearest-neighbor backend. One of `"auto"` or `"exact"`. Default is
  `"auto"`, which lets `thisutils` choose the fastest exact backend
  available. Requires the accelerated
  [`thisutils::compute_lisi()`](https://mengxu98.github.io/thisutils/reference/compute_lisi.html)
  interface that exposes `nn_method = c("auto", "exact")`.

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
#> ◌ [2026-05-23 08:55:12] Run integration workflow...
#> Warning: No layers found matching search pattern provided
#> ℹ [2026-05-23 08:55:13] Perform `Seurat::NormalizeData()` on split layers for Seurat v5 integration
#> ℹ [2026-05-23 08:55:16] Perform `Seurat::FindVariableFeatures()` per batch (`HVF_source = 'separate'`)
#> ℹ [2026-05-23 08:55:17] Number of available HVF: 2000
#> Warning: Layer ‘scale.data’ is empty
#> ℹ [2026-05-23 08:55:18] Perform `Seurat::ScaleData()` on split layers for Seurat v5 integration
#> ℹ [2026-05-23 08:55:19] Perform PCA on split layers before `Seurat::IntegrateLayers()`
#> ℹ [2026-05-23 08:55:19] Perform Seurat v5 integration with `HarmonyIntegration()`
#> The `features` argument is ignored by `HarmonyIntegration`.
#> This message is displayed once per session.
#> ! [2026-05-23 08:55:20] No valid estimated dimensions found for Harmony5. Use fallback dimensions 1:50
#> ℹ [2026-05-23 08:55:20] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-05-23 08:55:21] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-05-23 08:55:21] Reorder clusters...
#> ℹ [2026-05-23 08:55:21] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-23 08:55:21] Perform umap nonlinear dimension reduction using Harmony5 (1:50)
#> ℹ [2026-05-23 08:55:28] Perform umap nonlinear dimension reduction using Harmony5 (1:50)
#> ℹ [2026-05-23 08:55:35] Perform umap nonlinear dimension reduction using pca (1:20)
#> ✔ [2026-05-23 08:55:42] Harmony5 integration completed
names(panc8_sub@reductions)
#> [1] "pca"            "Harmony5"       "Harmony5UMAP2D" "Harmony5UMAP3D"
#> [5] "pcaUMAP2D"     

panc8_sub <- RunLISI(
  panc8_sub,
  reductions = c("pcaUMAP2D", "Harmony5UMAP2D")
)
#> ℹ [2026-05-23 08:55:42] Compute LISI scores from reduction "pcaUMAP2D"
#> ℹ [2026-05-23 08:55:42] Compute LISI scores from reduction "Harmony5UMAP2D"
#> ✔ [2026-05-23 08:55:42] Stored LISI scores in metadata: "pcaUMAP2D_tech_LISI" and "Harmony5UMAP2D_tech_LISI"
LISIPlot(
  panc8_sub,
  combine = TRUE
)
```
