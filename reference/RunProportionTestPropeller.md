# Propeller differential abundance wrapper

Method-specific propeller-style wrapper used by
[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md).

## Usage

``` r
RunProportionTestPropeller(
  srt,
  group.by,
  split.by,
  sample.by,
  comparison = NULL,
  n_bootstrap = 1000,
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

- n_bootstrap:

  Bootstrap iterations for confidence intervals.

- seed:

  Random seed.

- verbose:

  Whether to print messages.

## See also

[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md)
