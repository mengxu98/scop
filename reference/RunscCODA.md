# scCODA differential abundance

scCODA differential abundance

## Usage

``` r
RunscCODA(
  srt,
  group.by,
  split.by,
  sample.by,
  comparison = NULL,
  reference_cell_type = NULL,
  credible_effect_threshold = 0.95,
  n_mcmc_samples = 20000L,
  reuse_reverse_comparisons = TRUE,
  envname = "sccoda_env",
  conda = "auto",
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

- reference_cell_type:

  Optional reference cell type for scCODA.

- credible_effect_threshold:

  Inclusion probability threshold for credible effects.

- n_mcmc_samples:

  Number of MCMC samples requested in scCODA.

- reuse_reverse_comparisons:

  Reuse each fitted pair for the reverse direction by negating effect
  estimates and swapping interval bounds. This avoids fitting the same
  two-group model twice.

- envname:

  Name of the conda-compatible environment used by scCODA. Defaults to
  `"sccoda_env"` to keep the TensorFlow/scCODA stack isolated from the
  default Python environment.

- conda:

  The path or command name of a conda-compatible executable.

- seed:

  Random seed.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A method result bundle used internally by
[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md).
