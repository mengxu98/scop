# Propeller differential abundance wrapper

Method-specific implementation used by
[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md)
when `proportion_method = "propeller"`. This implementation works on
sample-level proportions using a propeller-style transformed test and
stores standardized outputs for plotting.

## Usage

``` r
RunPropeller(
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

  A `Seurat` object.

- group.by:

  Metadata column(s) used to color cells.

- split.by:

  Metadata column that identifies the condition groups to compare. For
  sample-level methods, if `split.by` is omitted and `sample.by` is
  provided, `sample.by` is treated as the condition column and virtual
  samples are created within each condition.

- sample.by:

  Metadata column that identifies biological samples. For `"milo"`,
  `"sccoda"`, and `"propeller"`, when `sample.by` is omitted or
  identical to `split.by`, virtual samples are created within each
  `split.by` group for convenience.

- comparison:

  Optional pairs of condition labels. Each pair is `c(group1, group2)`.
  `obs_log2FD` is `log2(group1 / group2)`, matching
  [RunDEtest](https://mengxu98.github.io/scop/reference/RunDEtest.md) /
  Seurat `ident.1` vs `ident.2`: positive values mean the cell type is
  more abundant in `group1`. If `NULL`, all pairwise comparisons are
  generated in both directions.

- n_bootstrap:

  Number of bootstrap iterations for confidence intervals.

- seed:

  Random seed.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A method result bundle used internally by
[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md).
