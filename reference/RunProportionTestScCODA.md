# scCODA differential abundance wrapper

Method-specific scCODA wrapper used by
[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md).
Calls Python helper via reticulate and falls back to sample-level
summaries when unavailable.

## Usage

``` r
RunProportionTestScCODA(
  srt,
  group.by,
  split.by,
  sample.by,
  comparison = NULL,
  reference_cell_type = NULL,
  credible_effect_threshold = 0.95,
  n_mcmc_samples = 20000L,
  n_bootstrap = 500,
  seed = 11,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- group.by:

  Cell-type metadata column.

- split.by:

  Condition metadata column.

- sample.by:

  Sample metadata column.

- comparison:

  Optional comparisons to perform.

- reference_cell_type:

  Optional scCODA reference cell type.

- credible_effect_threshold:

  Inclusion probability threshold for credible effects.

- n_mcmc_samples:

  Requested MCMC sample count.

- n_bootstrap:

  Fallback bootstrap iterations.

- seed:

  Random seed.

- verbose:

  Whether to print messages.

## See also

[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md)
