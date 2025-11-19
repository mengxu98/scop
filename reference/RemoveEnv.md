# Remove a conda environment

Remove a conda environment

## Usage

``` r
RemoveEnv(envname = NULL, conda = "auto", force = FALSE)
```

## Arguments

- envname:

  The name of the conda environment to remove. If `NULL`, uses the
  default scop environment name.

- conda:

  The path to a `conda` executable. Use `"auto"` to allow `reticulate`
  to automatically find an appropriate `conda` binary. See **Finding
  Conda** and
  [`conda_binary()`](https://rstudio.github.io/reticulate/reference/conda-tools.html)
  for more details.

- force:

  Whether to force removal without confirmation. Default is `FALSE`.

## Value

Invisibly returns TRUE if successful, FALSE otherwise.

## Examples

``` r
if (FALSE) { # \dontrun{
# Remove the default scop environment
RemoveEnv()

# Remove a specific environment
RemoveEnv("my_old_env")

# Force removal without confirmation
RemoveEnv("my_old_env", force = TRUE)
} # }
```
