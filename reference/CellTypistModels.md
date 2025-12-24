# Get available CellTypist models

Get a list of all available CellTypist models, either downloaded locally
or available from the server.

## Usage

``` r
CellTypistModels(on_the_fly = FALSE, verbose = TRUE)
```

## Arguments

- on_the_fly:

  If `TRUE`, only show downloaded models. If `FALSE`, show all available
  models (fetch list from server). Default is `FALSE`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A data frame containing model names and descriptions.

## See also

[RunCellTypist](https://mengxu98.github.io/scop/reference/RunCellTypist.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Get available models
models <- CellTypistModels()
print(models)
} # }
```
