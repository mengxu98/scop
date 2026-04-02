# Estimate useful dimensions from a reduction

Estimate useful dimensions from a reduction

## Usage

``` r
RunDimsEstimate(
  srt,
  reduction,
  reduction_method = NULL,
  k = 20L,
  min_dims = 10L,
  fallback_max_dims = 50L,
  skip_first = FALSE,
  use_stored = TRUE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A `Seurat` object.

- reduction:

  Name of the dimensional reduction to inspect.

- reduction_method:

  Optional reduction method name. When set to `"nmf"` or `"glmpca"`, all
  available dimensions will be retained.

- k:

  Number of neighbors used by intrinsicDimension::maxLikGlobalDimEst.
  Default is `20`.

- min_dims:

  Minimum number of dimensions kept when intrinsic-dimension estimation
  succeeds. Default is `10`.

- fallback_max_dims:

  Maximum number of dimensions kept when no valid estimate is available.
  Default is `50`.

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
