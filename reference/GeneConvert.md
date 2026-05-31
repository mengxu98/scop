# Gene ID conversion function using biomart

This function can convert different gene ID types within one species or
between two species using the biomart service.

## Usage

``` r
GeneConvert(
  geneID,
  geneID_from_IDtype = "symbol",
  geneID_to_IDtype = "entrez_id",
  species_from = "Homo_sapiens",
  species_to = NULL,
  Ensembl_version = NULL,
  biomart = NULL,
  mirror = NULL,
  max_tries = 5,
  verbose = TRUE
)
```

## Arguments

- geneID:

  A vector of the geneID character.

- geneID_from_IDtype:

  Gene ID type of the input `geneID`. e.g. `"symbol"`, `"ensembl_id"`,
  `"entrez_id"`

- geneID_to_IDtype:

  Gene ID type(s) to convert to. e.g. `"symbol"`, `"ensembl_id"`,
  `"entrez_id"`.

- species_from:

  Latin names for animals of the input geneID. e.g. `"Homo_sapiens"`,
  `"Mus_musculus"`.

- species_to:

  Latin names for animals of the output geneID. e.g. `"Homo_sapiens"`,
  `"Mus_musculus"`.

- Ensembl_version:

  An integer specifying the Ensembl version. Default is `NULL`. If
  `NULL`, the latest version will be used.

- biomart:

  The name of the BioMart database that you want to connect to. Possible
  options include `"ensembl"`, `"protists_mart"`, `"fungi_mart"`, and
  `"plants_mart"`.

- mirror:

  Specify an Ensembl mirror to connect to. The valid options here are
  `"www"`, `"uswest"`, `"useast"`, `"asia"`.

- max_tries:

  The maximum number of attempts to connect with the BioMart service.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A list with the following elements:

- `geneID_res`: A data.frame contains the all gene IDs mapped in the
  database with columns: `"from_IDtype"`, `"from_geneID"`,
  `"to_IDtype"`, `"to_geneID"`.

- `geneID_collapse`: The data.frame contains all the successfully
  converted gene IDs, and the output gene IDs are collapsed into a list.
  As a result, the `"from_geneID"` column (which is set as the row
  names) of the data.frame is unique.

- `geneID_expand`: The data.frame contains all the successfully
  converted gene IDs, and the output gene IDs are expanded.

- `Ensembl_version`: Ensembl database version.

- `Datasets`: Datasets available in the selected BioMart database.

- `Attributes`: Attributes available in the selected BioMart database.

- `geneID_unmapped`: A character vector of gene IDs that are unmapped in
  the database.

## See also

[AnnotateFeatures](https://mengxu98.github.io/scop/reference/AnnotateFeatures.md),
[ConvertHomologs](https://mengxu98.github.io/scop/reference/ConvertHomologs.md)

## Examples

``` r
res <- GeneConvert(
  geneID = c("CDK1", "MKI67", "TOP2A", "AURKA", "CTCF"),
  species_from = "Homo_sapiens",
  species_to = "Mus_musculus"
)
#> ℹ [2026-05-31 06:08:43] Connect to the Ensembl archives...
#> ℹ [2026-05-31 06:08:43] Using the 115 version of ensembl database...
#> ℹ [2026-05-31 06:08:43] Downloading the ensembl database from https://sep2025.archive.ensembl.org...
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying www mirror
#> Warning: Invalid mirror. Select a mirror from [www, useast, asia].
#> Default when no mirror is specified is to use www.ensembl.org which may be automatically redirected.
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying useast mirror
#> ! [2026-05-31 06:08:48] <simpleError: Your query has been redirected to https://status.ensembl.org indicating this Ensembl service is currently unavailable.
#> !                       Look at ?useEnsembl for details on how to try a mirror site.>
#> ! [2026-05-31 06:08:48] Get errors when connecting with ensembl database...
#> ! [2026-05-31 06:08:49] Retrying...
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying www mirror
#> Warning: Invalid mirror. Select a mirror from [www, useast, asia].
#> Default when no mirror is specified is to use www.ensembl.org which may be automatically redirected.
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying useast mirror
#> ! [2026-05-31 06:08:53] <simpleError: Your query has been redirected to https://status.ensembl.org indicating this Ensembl service is currently unavailable.
#> !                       Look at ?useEnsembl for details on how to try a mirror site.>
#> ! [2026-05-31 06:08:53] Get errors when connecting with ensembl database...
#> ! [2026-05-31 06:08:54] Retrying...
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying asia mirror
#> Warning: Invalid mirror. Select a mirror from [www, useast, asia].
#> Default when no mirror is specified is to use www.ensembl.org which may be automatically redirected.
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying asia mirror
#> ! [2026-05-31 06:08:57] <simpleError: Your query has been redirected to https://status.ensembl.org indicating this Ensembl service is currently unavailable.
#> !                       Look at ?useEnsembl for details on how to try a mirror site.>
#> ! [2026-05-31 06:08:57] Get errors when connecting with ensembl database...
#> ! [2026-05-31 06:08:58] Retrying...
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying asia mirror
#> Warning: Invalid mirror. Select a mirror from [www, useast, asia].
#> Default when no mirror is specified is to use www.ensembl.org which may be automatically redirected.
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying useast mirror
#> ! [2026-05-31 06:09:01] <simpleError: Your query has been redirected to https://status.ensembl.org indicating this Ensembl service is currently unavailable.
#> !                       Look at ?useEnsembl for details on how to try a mirror site.>
#> ! [2026-05-31 06:09:01] Get errors when connecting with ensembl database...
#> ! [2026-05-31 06:09:02] Retrying...
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying asia mirror
#> Warning: Invalid mirror. Select a mirror from [www, useast, asia].
#> Default when no mirror is specified is to use www.ensembl.org which may be automatically redirected.
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying asia mirror
#> ! [2026-05-31 06:09:05] <simpleError: Your query has been redirected to https://status.ensembl.org indicating this Ensembl service is currently unavailable.
#> !                       Look at ?useEnsembl for details on how to try a mirror site.>
#> ! [2026-05-31 06:09:05] Get errors when connecting with ensembl database...
#> Error in try_get(expr = {    if (!is.null(mirror)) {        biomaRt::useEnsembl(biomart = "ensembl", mirror = mirror)    }    else {        mart_try <- tryCatch(biomaRt::useMart(biomart = "ensembl",             host = url), error = function(e) e)        if (inherits(mart_try, "error")) {            mirror_candidates <- c("useast", "uswest", "asia",                 "www")            for (mirror_i in mirror_candidates) {                mart_mirror <- tryCatch(biomaRt::useEnsembl(biomart = "ensembl",                   mirror = mirror_i), error = function(e) e)                if (!inherits(mart_mirror, "error")) {                  mart_try <- mart_mirror                  break                }            }        }        if (inherits(mart_try, "error")) {            stop(mart_try)        }        mart_try    }}, max_tries = max_tries, error_message = "Get errors when connecting with ensembl database..."): <simpleError: Your query has been redirected to
#> https://status.ensembl.org indicating this Ensembl service is currently
#> unavailable. Look at ?useEnsembl for details on how to try a mirror site.>
str(res)
#> Error: object 'res' not found
```
