# Permutation-based proportion test

Method-specific permutation implementation used by
[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md).
This is a permutation-based statistical test for compositional shift
(not a separate biological model). Bootstrap confidence intervals are
included as uncertainty estimates for observed log2 fold-differences.

## Usage

``` r
RunProportionTestPermutation(
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

  Cell-type metadata column.

- split.by:

  Condition metadata column.

- comparison:

  Optional comparisons to perform.

- n_permutations:

  Number of permutations.

- include_all_cells:

  Whether to include all cell types in grid completion.

- verbose:

  Whether to print messages.

## See also

[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md)
