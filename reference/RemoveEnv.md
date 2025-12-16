# Remove a conda environment

Remove a conda environment

## Usage

``` r
RemoveEnv(envname = NULL, conda = "auto", force = FALSE)
```

## Arguments

- envname:

  The name of the conda environment. If `NULL`, the environment name
  will be set to `"scop_env"`. Default is `NULL`.

- conda:

  The path to a conda executable. Use `"auto"` to allow automatically
  finding an appropriate conda binary.

- force:

  Whether to force removal without confirmation. Default is `FALSE`.

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
