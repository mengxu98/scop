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
  envname = "scop_sccoda_env",
  conda = "auto",
  seed = 11,
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

- sample.by:

  Metadata column that identifies biological samples. For `"milo"`,
  `"sccoda"`, and `"propeller"`, when `sample.by` is omitted or
  identical to `split.by`, virtual samples are created within each
  `split.by` group for convenience.

- comparison:

  Optional: specify comparisons to perform.

- reference_cell_type:

  Optional reference cell type for scCODA.

- credible_effect_threshold:

  Inclusion probability threshold for credible effects.

- n_mcmc_samples:

  Number of MCMC samples requested in scCODA.

- envname:

  Name of the conda-compatible environment used by scCODA. Defaults to
  `"scop_sccoda_env"` to keep the TensorFlow/scCODA stack isolated from
  the default Python environment.

- conda:

  The path or command name of a conda-compatible executable.

- seed:

  Random seed.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A method result bundle used internally by
[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md).
