# Permutation-based proportion test

Method-specific implementation used by
[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md)
when `proportion_method = "permutation"`. This method is a
permutation-based statistical test for compositional shift rather than a
dedicated biological model. Bootstrap confidence intervals are reported
as uncertainty estimates for observed log2 fold-differences
(`log2(group1 / group2)`).

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

  A `Seurat` object.

- group.by:

  Metadata column(s) used to color cells.

- split.by:

  Metadata column that identifies the condition groups to compare. For
  sample-level methods, if `split.by` is omitted and `sample.by` is
  provided, `sample.by` is treated as the condition column and virtual
  samples are created within each condition.

- comparison:

  Optional pairs of condition labels. Each pair is `c(group1, group2)`.
  `obs_log2FD` is `log2(group1 / group2)`, matching
  [RunDEtest](https://mengxu98.github.io/scop/reference/RunDEtest.md) /
  Seurat `ident.1` vs `ident.2`: positive values mean the cell type is
  more abundant in `group1`. If `NULL`, all pairwise comparisons are
  generated in both directions.

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
