# Permutation-based proportion test

Method-specific implementation used by
[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md)
when `proportion_method = "permutation"`. This method is a
permutation-based statistical test for compositional shift rather than a
dedicated biological model. Bootstrap confidence intervals are reported
as uncertainty estimates for observed log2 fold-differences.

## Usage

``` r
RunPermutation(
  srt,
  group.by,
  split.by,
  comparison = NULL,
  n_permutations = 1000,
  include_all_cells = FALSE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- group.by:

  Name of one or more meta.data columns to group (color) cells by.

- split.by:

  Metadata column that identifies the condition groups to compare. For
  sample-level methods, if `split.by` is omitted and `sample.by` is
  provided, `sample.by` is treated as the condition column and virtual
  samples are created within each condition.

- comparison:

  Optional: specify comparisons to perform.

- n_permutations:

  Number of permutations for permutation-based test.

- include_all_cells:

  Whether to include all cell types in the complete grid for permutation
  mode.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A method result bundle used internally by
[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md).
