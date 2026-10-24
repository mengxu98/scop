# Propeller differential abundance wrapper

Method-specific implementation used by
[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md)
when `proportion_method = "propeller"`. Uses the optional speckle
backend: `getTransformedProps()` computes the transformed proportions
and `propeller.ttest()` performs limma empirical Bayes moderated tests.
Each requested comparison is a two-group test; multiple conditions
produce pairwise results, not an omnibus ANOVA. P-values, FDR, and
`Tstatistic` come from speckle. `obs_log2FD` and bootstrap intervals
describe the ratio of mean untransformed sample proportions, with a
`1e-5` pseudocount; these intervals are not moderated-test intervals.

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
  srt = NULL,
  transform = c("logit", "asin"),
  robust = TRUE,
  trend = FALSE
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
  conditions add a sample blocking effect to the speckle design matrix;
  partially paired comparisons are rejected. At least two biological
  samples per condition are required.

- comparison:

  Optional pairs of condition labels. Each pair is `c(group1, group2)`.
  `obs_log2FD` is `log2(group1 / group2)`, matching
  [RunDEtest](https://mengxu98.github.io/scop/reference/RunDEtest.md) /
  Seurat `ident.1` vs `ident.2`: positive values mean the cell type is
  more abundant in `group1`. If `NULL`, all pairwise comparisons are
  generated in both directions.

- n_bootstrap:

  Number of sample bootstrap iterations for descriptive confidence
  intervals. Paired samples are resampled together. Use zero to omit
  bootstrap intervals.

- seed:

  Random seed, including for permutation testing.

- verbose:

  Whether to print the message. Default is `TRUE`.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

- transform:

  Proportion transformation, `"logit"` (default) or `"asin"`. Logit uses
  speckle's count pseudocount of 0.5 before computing proportions.

- robust:

  Whether speckle uses robust empirical Bayes variance shrinkage.

- trend:

  Whether speckle fits a mean-variance trend.

## Value

A method result bundle used internally by
[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md).
