# Prepare CytoTRACE2 model resources

Download or load the versioned model assets used by \[RunCytoTRACE()\].

## Usage

``` r
PrepareCytoTRACE2(db_update = FALSE, verbose = TRUE)
```

## Arguments

- db_update:

  Force a refresh. `FALSE` loads the cache when available.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A list with a \`CytoTRACE2\` entry containing \`data_dir\`, \`files\`
and \`version\`.

## Details

Model assets are species-independent and cached in the user data
directory. This entry retains the structure used by \[RunCytoTRACE()\].

## See also

\[PrepareDB\], \[RunCytoTRACE\]

## Examples

``` r
if (FALSE) { # \dontrun{
model <- PrepareCytoTRACE2()
model[["CytoTRACE2"]]$version
} # }
```
