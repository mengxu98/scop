# List SCOP external datasets

Read a dataset manifest from the external `mengxu98/datasets` repository
or a local mirror. This keeps example data assets out of the SCOP
package while still making them discoverable and reproducible.

## Usage

``` r
ListScopDatasets(
  collection = "Xenium",
  datasets_base_url = "https://raw.githubusercontent.com/mengxu98/datasets/main"
)
```

## Arguments

- collection:

  Dataset collection directory, for example `"Xenium"`.

- datasets_base_url:

  Base URL or local directory containing SCOP dataset collections.

## Value

A data frame parsed from `manifest.tsv`.

## Examples

``` r
if (FALSE) { # \dontrun{
ListScopDatasets("Xenium")
} # }
```
