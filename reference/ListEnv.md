# List conda-compatible Python environments

List conda-compatible Python environments

## Usage

``` r
ListEnv(conda = "auto")
```

## Arguments

- conda:

  The path or command name of a conda-compatible executable (`conda`,
  `mamba`, or `micromamba`). Use `"auto"` to allow automatically finding
  an appropriate environment manager. If `"micromamba"` is requested and
  micromamba is not available on `PATH`, a package-managed micromamba is
  downloaded automatically.

## Value

A data frame of conda-compatible Python environments.
