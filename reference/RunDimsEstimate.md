# Estimate useful dimensions from a reduction

Estimate useful dimensions from a reduction

## Usage

``` r
RunDimsEstimate(
  srt,
  reduction = NULL,
  reduction_method = NULL,
  k = 30L,
  method = c("scree", "intrinsic", "ensemble"),
  min_dims = 5L,
  fallback_max_dims = 50L,
  variance_threshold = 0.8,
  marginal_gain_threshold = 0.5,
  skip_first = FALSE,
  use_stored = TRUE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A `Seurat` object.

- reduction:

  Name of the dimensional reduction to inspect. Default is `NULL`, which
  automatically selects a PCA-like reduction via
  [`DefaultReduction()`](https://mengxu98.github.io/scop/reference/DefaultReduction.md)
  with `pattern = "pca"`.

- reduction_method:

  Optional reduction method name. When set to `"nmf"` or `"glmpca"`, all
  available dimensions will be retained.

- k:

  Number of neighbors used by
  [intrinsicDimension::maxLikGlobalDimEst](https://rdrr.io/pkg/intrinsicDimension/man/maxLik.html).
  Default is `30`.

- method:

  Dimension-selection method. `"scree"` uses PCA standard deviations
  with broken-stick, elbow, cumulative-variance, and marginal-gain
  criteria. `"intrinsic"` uses
  [intrinsicDimension::maxLikGlobalDimEst](https://rdrr.io/pkg/intrinsicDimension/man/maxLik.html).
  `"ensemble"` keeps the larger recommendation from both methods when
  both are available. Default is `"scree"`.

- min_dims:

  Minimum number of dimensions kept when intrinsic-dimension estimation
  succeeds. Default is `5`.

- fallback_max_dims:

  Maximum number of dimensions kept when no valid estimate is available.
  Default is `50`.

- variance_threshold:

  Cumulative variance threshold used by `method = "scree"`. Default is
  `0.8`.

- marginal_gain_threshold:

  Stop point for marginal variance gain (percentage points) used by
  `method = "scree"`. Default is `0.5`.

- skip_first:

  Whether to drop the first dimension from the returned result. Useful
  for `TFIDF/LSI` workflows. Default is `FALSE`.

- use_stored:

  Whether to use `misc$dims_estimate` already stored in the reduction
  when available. Default is `TRUE`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

An integer vector of dimensions to use.

## See also

[DimsEstimatePlot](https://mengxu98.github.io/scop/reference/DimsEstimatePlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
RunDimsEstimate(pancreas_sub)

DimsEstimatePlot(pancreas_sub)
```
