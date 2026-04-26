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
  nn_eps = 0,
  use_rann = TRUE,
  nn_method = c("auto", "rann", "fnn", "exact"),
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

  Effective neighborhood size. Defaults to `30`.

- nn_eps:

  Approximation factor passed to
  [RANN::nn2](https://jefferislab.github.io/RANN/reference/nn2.html)
  when `RANN` is available and `use_rann = TRUE`. Defaults to `0`.

- use_rann:

  Whether to prefer
  [RANN::nn2](https://jefferislab.github.io/RANN/reference/nn2.html)
  over [`FNN::get.knn`](https://rdrr.io/pkg/FNN/man/get.knn.html) when
  `nn_method = "auto"` decides not to use the package's built-in exact
  C++ backend. Defaults to `TRUE`.

- nn_method:

  Nearest-neighbor backend. Defaults to `"auto"`, which uses a simple
  heuristic: low-dimensional inputs use the package's exact C++ search,
  while larger/higher-dimensional inputs fall back to `RANN`, then
  `FNN`, then the built-in exact C++ backend. Set to `"hnsw"` to use
  `RcppHNSW` for a faster, approximate search.

- tol:

  Tolerance used in the binary search for the target perplexity.
  Defaults to `1e-5`.

- max_iter:

  Maximum number of binary-search iterations. Defaults to `50`.

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
#> ◌ [2026-04-26 02:22:33] Run integration workflow...
#> Warning: No layers found matching search pattern provided
#> ℹ [2026-04-26 02:22:34] Perform `Seurat::NormalizeData()` on split layers for Seurat v5 integration
#> ℹ [2026-04-26 02:22:36] Perform `Seurat::FindVariableFeatures()` per batch (`HVF_source = 'separate'`)
#> ℹ [2026-04-26 02:22:38] Number of available HVF: 2000
#> Warning: Layer ‘scale.data’ is empty
#> ℹ [2026-04-26 02:22:39] Perform `Seurat::ScaleData()` on split layers for Seurat v5 integration
#> ℹ [2026-04-26 02:22:39] Perform PCA on split layers before `Seurat::IntegrateLayers()`
#> ℹ [2026-04-26 02:22:41] Perform Seurat v5 integration with `HarmonyIntegration()`
#> The `features` argument is ignored by `HarmonyIntegration`.
#> This message is displayed once per session.
#> ℹ [2026-04-26 02:22:42] Estimated dimensions 1:20 for Harmony5
#> ℹ [2026-04-26 02:22:42] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-04-26 02:22:42] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-04-26 02:22:43] Reorder clusters...
#> ℹ [2026-04-26 02:22:43] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-26 02:22:43] Perform umap nonlinear dimension reduction using Harmony5 (1:20)
#> ℹ [2026-04-26 02:22:49] Perform umap nonlinear dimension reduction using Harmony5 (1:20)
#> ℹ [2026-04-26 02:22:55] Perform umap nonlinear dimension reduction using pca (1:20)
#> ✔ [2026-04-26 02:23:02] Harmony5 integration completed
names(panc8_sub@reductions)
#> [1] "pca"            "Harmony5"       "Harmony5UMAP2D" "Harmony5UMAP3D"
#> [5] "pcaUMAP2D"     

panc8_sub <- RunLISI(
  panc8_sub,
  reductions = c("pcaUMAP2D", "Harmony5UMAP2D")
)
#> ℹ [2026-04-26 02:23:02] Compute LISI scores from reduction "pcaUMAP2D"
#> ◌ [2026-04-26 02:23:02] Using "exact" nearest-neighbor backend for compute_lisi
#> ℹ [2026-04-26 02:23:02] Compute LISI scores from reduction "Harmony5UMAP2D"
#> ◌ [2026-04-26 02:23:02] Using "exact" nearest-neighbor backend for compute_lisi
#> ✔ [2026-04-26 02:23:02] Stored LISI scores in metadata: "pcaUMAP2D_tech_LISI" and "Harmony5UMAP2D_tech_LISI"
LISIPlot(
  panc8_sub,
  combine = TRUE
)
```
