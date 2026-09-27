# Propeller differential abundance wrapper

Method-specific implementation used by
[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md)
when `proportion_method = "propeller"`. This implementation uses an
internal logit-transformed sample-proportion t-test (paired when sample
IDs occur in both conditions), not the speckle implementation of
propeller. It stores standardized outputs for plotting.

## Usage

``` r
RunPropeller(
  object,
  group.by,
  split.by,
  sample.by,
  comparison = NULL,
  n_bootstrap = 1000,
  seed = 11,
  verbose = TRUE,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- group.by:

  Metadata column(s) used to color cells.

- split.by:

  Metadata column that identifies the condition groups to compare. For
  sample-level methods, if `split.by` is omitted and `sample.by` is
  provided, `sample.by` is treated as the condition column; a separate
  biological sample column is still required unless descriptive virtual
  samples are explicitly enabled.

- sample.by:

  Metadata column identifying biological samples. Required when calling
  `RunPropeller()` directly. Fully paired sample IDs across two
  conditions use a paired test; partially paired comparisons are
  rejected.

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

  Random seed, including for permutation testing.

- verbose:

  Whether to print the message. Default is `TRUE`.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

A method result bundle used internally by
[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md).
