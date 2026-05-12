# Proportion Test

RunProportionTest performs differential abundance testing for cell
proportions. The function acts as a dispatcher and routes to one of the
method-specific implementations: permutation, milo, sccoda, or
propeller.

## Usage

``` r
RunProportionTest(
  srt,
  group.by,
  split.by = NULL,
  comparison = NULL,
  proportion_method,
  sample.by = NULL,
  pseudo_sample_n = 3L,
  n_permutations = 1000,
  FDR_threshold = 0.05,
  log2FD_threshold = log2(1.5),
  include_all_cells = FALSE,
  seed = 11,
  verbose = TRUE,
  ...
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

- proportion_method:

  Differential abundance method. One of `"permutation"`, `"milo"`,
  `"sccoda"`, or `"propeller"`. This argument is required. Alias values
  such as `"permutation_test"` and `"perm"` are accepted and normalized
  to `"permutation"`.

- sample.by:

  Metadata column that identifies biological samples. For `"milo"`,
  `"sccoda"`, and `"propeller"`, when `sample.by` is omitted or
  identical to `split.by`, virtual samples are created within each
  `split.by` group for convenience.

- pseudo_sample_n:

  Number of virtual samples per `split.by` group when a sample-level
  method has no usable `sample.by`.

- n_permutations:

  Number of permutations for permutation-based test.

- FDR_threshold:

  FDR value cutoff for significance.

- log2FD_threshold:

  Absolute value of log2FD cutoff for significance.

- include_all_cells:

  Whether to include all cell types in the complete grid for permutation
  mode.

- seed:

  Random seed.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments passed to the selected method function.

## References

[Miller et al. paper](https://doi.org/10.1158/0008-5472.can-20-3562),
[scProportionTest](https://github.com/rpolicastro/scProportionTest),
[miloR](https://bioconductor.org/packages/miloR),
[scCODA](https://github.com/theislab/scCODA),
[propeller/speckle](https://bioconductor.org/packages/speckle)

## See also

[RunPermutation](https://mengxu98.github.io/scop/reference/RunPermutation.md),
[RunMilo](https://mengxu98.github.io/scop/reference/RunMilo.md),
[RunscCODA](https://mengxu98.github.io/scop/reference/RunscCODA.md),
[RunPropeller](https://mengxu98.github.io/scop/reference/RunPropeller.md),
[ProportionTestPlot](https://mengxu98.github.io/scop/reference/ProportionTestPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunProportionTest(
  pancreas_sub,
  group.by = "CellType",
  split.by = "Phase",
  proportion_method = "permutation",
  comparison = list(c("G2M", "G1"))
)
#> ℹ [2026-05-12 15:54:08] Start proportion test ("permutation")
#> ℹ [2026-05-12 15:54:08] Running comparison: "G1" vs "G2M"
#> ℹ [2026-05-12 15:54:08] Running comparison: "G2M" vs "G1"
#> ✔ [2026-05-12 15:54:08] Proportion test completed ("permutation")

ProportionTestPlot(
  pancreas_sub
)
```
