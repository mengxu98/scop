# Proportion Test

RunProportionTest performs a Monte-carlo permutation test to quantify
the cell proportion differences between each condition.

## Usage

``` r
RunProportionTest(
  srt,
  group.by,
  split.by,
  comparison = NULL,
  n_permutations = 1000,
  FDR_threshold = 0.05,
  log2FD_threshold = log2(1.5),
  include_all_cells = FALSE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object containing count data and metadata.

- group.by:

  The name of the metadata column in `srt` containing cluster labels or
  cell type names.

- split.by:

  The name of the metadata column in `srt` that contains sample
  identifiers.

- comparison:

  Optional: specify comparisons to perform.

- n_permutations:

  Number of permutations for the test.

- FDR_threshold:

  FDR value cutoff for significance.

- log2FD_threshold:

  Absolute value of log2FD cutoff for significance.

- include_all_cells:

  Whether to include all cell types in the complete grid (default:
  FALSE).

- verbose:

  Whether to print the message. Default is `TRUE`.

## References

[Miller et al. paper](https://doi.org/10.1158/0008-5472.can-20-3562),
[scProportionTest](https://github.com/rpolicastro/scProportionTest)

## See also

[ProportionTestPlot](https://mengxu98.github.io/scop/reference/ProportionTestPlot.md)

## Examples

``` r
data(pancreas_sub)
# Default behavior: only include cell types present in comparison groups
pancreas_sub <- RunProportionTest(
  pancreas_sub,
  group.by = "CellType",
  split.by = "Phase",
  comparison = list(c("G2M", "G1"))
)
#> ℹ [2025-11-19 15:00:00] Start proportion test
#> ℹ [2025-11-19 15:00:00] Running comparison: "G1" vs "G2M"
#> ℹ [2025-11-19 15:00:07] Running comparison: "G2M" vs "G1"
#> ✔ [2025-11-19 15:00:14] Proportion test completed

# Include all cell types from the dataset
pancreas_sub <- RunProportionTest(
  pancreas_sub,
  group.by = "CellType",
  split.by = "Phase",
  comparison = list(c("G2M", "G1")),
  include_all_cells = TRUE
)
#> ℹ [2025-11-19 15:00:14] Start proportion test
#> ℹ [2025-11-19 15:00:14] Running comparison: "G1" vs "G2M"
#> ℹ [2025-11-19 15:00:21] Running comparison: "G2M" vs "G1"
#> ✔ [2025-11-19 15:00:30] Proportion test completed

ProportionTestPlot(
  pancreas_sub
)
```
