# Remove Python packages from conda environment

Remove Python packages from conda environment

## Usage

``` r
remove_python(
  packages,
  envname = NULL,
  conda = "auto",
  pip = FALSE,
  force = FALSE,
  verbose = TRUE
)
```

## Arguments

- packages:

  A character vector of package names to remove.

- envname:

  The name of the conda environment. If `NULL`, the environment name
  will be set to `"scop_env"`. Default is `NULL`.

- conda:

  The path to a conda executable. Use `"auto"` to allow automatically
  finding an appropriate conda binary.

- pip:

  Whether to use pip for package removal. Default is `FALSE` (use
  conda).

- force:

  Whether to force removal without confirmation. Default is `FALSE`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

Invisibly returns.

## Examples

``` r
if (FALSE) { # \dontrun{
# Remove a single package using conda
remove_python("numpy")

# Remove multiple packages using conda
remove_python(c("numpy", "pandas"))

# Remove packages using pip
remove_python("numpy", pip = TRUE)

# Force removal without confirmation
remove_python("numpy", force = TRUE)

# Remove packages from a specific environment
remove_python("numpy", envname = "env")
} # }
```
