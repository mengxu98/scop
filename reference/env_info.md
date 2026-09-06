# Print environment information

Print environment information

## Usage

``` r
env_info(conda, envname, verbose = TRUE)
```

## Arguments

- conda:

  The path or command name of a conda-compatible executable (`conda`,
  `mamba`, or `micromamba`). Use `"auto"` to allow automatically finding
  an appropriate environment manager. If `"micromamba"` is requested and
  micromamba is not available on `PATH`, a package-managed micromamba is
  downloaded automatically.

- envname:

  Conda-compatible Python environment. If `NULL`, the environment name
  will be set to `"scop_env"`.

- verbose:

  Whether to print environment information.
