# Proportion Test

RunProportionTest performs differential abundance testing for cell
proportions through a unified dispatcher.

## Usage

``` r
RunProportionTest(
  srt,
  group.by,
  split.by,
  comparison = NULL,
  proportion_method,
  sample.by = NULL,
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

  Name of a metadata column for cell type/group labels.

- split.by:

  Name of a metadata column defining conditions to compare.

- comparison:

  Optional comparisons to perform. Supports `list(c("A", "B"))` or
  character values like `"A_vs_B"`.

- proportion_method:

  Differential abundance method. Canonical values are `"permutation"`,
  `"milo"`, `"sccoda"`, and `"propeller"`. Alias values (for example
  `"permutation_test"` or `"perm"`) are accepted and normalized to
  `"permutation"`.

- sample.by:

  Metadata column for biological sample IDs. Required for `"milo"`,
  `"sccoda"`, and `"propeller"`.

- n_permutations:

  Number of permutations for permutation mode.

- FDR_threshold:

  FDR value cutoff for significance.

- log2FD_threshold:

  Absolute value of log2FD cutoff for significance.

- include_all_cells:

  Whether to include all cell types in permutation complete grid.

- seed:

  Random seed.

- verbose:

  Whether to print messages.

- ...:

  Additional arguments passed to the selected method function.

## See also

[RunProportionTestPermutation](https://mengxu98.github.io/scop/reference/RunProportionTestPermutation.md),
[RunProportionTestMilo](https://mengxu98.github.io/scop/reference/RunProportionTestMilo.md),
[RunProportionTestScCODA](https://mengxu98.github.io/scop/reference/RunProportionTestScCODA.md),
[RunProportionTestPropeller](https://mengxu98.github.io/scop/reference/RunProportionTestPropeller.md),
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
#> ℹ [2026-04-26 02:25:03] Start proportion test ("permutation")
#> ℹ [2026-04-26 02:25:03] Running comparison: "G1" vs "G2M"
#> ℹ [2026-04-26 02:25:10] Running comparison: "G2M" vs "G1"
#> ✔ [2026-04-26 02:25:17] Proportion test completed ("permutation")

ProportionTestPlot(pancreas_sub)
```
