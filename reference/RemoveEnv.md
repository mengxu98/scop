# Remove a conda-compatible Python environment

Remove a conda-compatible Python environment

## Usage

``` r
RemoveEnv(envname = NULL, conda = "auto", force = FALSE, verbose = TRUE)
```

## Arguments

- envname:

  The name of the conda-compatible Python environment. If `NULL`, the
  environment name will be set to `"scop_env"`. Default is `NULL`.

- conda:

  The path or command name of a conda-compatible executable (`conda`,
  `mamba`, or `micromamba`). Use `"auto"` to allow automatically finding
  an appropriate environment manager. If `"micromamba"` is requested and
  micromamba is not available on `PATH`, a package-managed micromamba is
  downloaded automatically.

- force:

  Whether to force removal without confirmation. Default is `FALSE`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

Invisibly returns `TRUE` if successful, `FALSE` otherwise.

## Examples

``` r
if (FALSE) { # \dontrun{
# Remove default environment
RemoveEnv()

# Remove a specific environment
RemoveEnv("my_old_env")

# Removal without confirmation
RemoveEnv("my_old_env", force = TRUE)
} # }
```
