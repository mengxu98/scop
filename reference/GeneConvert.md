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

[AnnotateFeatures](https://mengxu98.github.io/scop/reference/AnnotateFeatures.md)

## Examples

``` r
if (FALSE) { # \dontrun{
res <- GeneConvert(
  geneID = c("CDK1", "MKI67", "TOP2A", "AURKA", "CTCF"),
  species_from = "Homo_sapiens",
  species_to = "Mus_musculus"
)
str(res)

# Convert the human genes to mouse homologs,
# and replace the raw counts in a Seurat object.
data(pancreas_sub)
counts <- GetAssayData5(
  pancreas_sub,
  assay = "RNA",
  layer = "counts"
)
res <- GeneConvert(
  geneID = rownames(counts),
  geneID_from_IDtype = "symbol",
  geneID_to_IDtype = "symbol",
  species_from = "Mus_musculus",
  species_to = "Homo_sapiens"
)

homologs_counts <- stats::aggregate(
  x = counts[res$geneID_expand[, "from_geneID"], ],
  by = list(res$geneID_expand[, "symbol"]), FUN = sum
)
rownames(homologs_counts) <- homologs_counts[, 1]
homologs_counts <- methods::as(
  thisutils::as_matrix(homologs_counts[, -1]),
  "dgCMatrix"
)
homologs_counts[1:5, 1:5]
} # }
```
