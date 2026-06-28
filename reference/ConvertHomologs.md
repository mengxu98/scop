# Convert homologous gene symbols in expression objects

Convert feature names between species with
[GeneConvert](https://mengxu98.github.io/scop/reference/GeneConvert.md)
and collapse duplicated target homologs by summing expression values.
The Seurat method rebuilds the selected assay from the converted counts
matrix and keeps cell metadata and spatial images when present.

## Usage

``` r
ConvertHomologs(
  object,
  species_from,
  species_to,
  geneID_from_IDtype = "symbol",
  geneID_to_IDtype = "symbol",
  assay = NULL,
  layer = "counts",
  multi_mapping = c("first"),
  keep_unmapped = FALSE,
  collapse_fun = c("sum"),
  Ensembl_version = NULL,
  biomart = NULL,
  mirror = NULL,
  max_tries = 5,
  verbose = TRUE
)

# S3 method for class 'Seurat'
ConvertHomologs(
  object,
  species_from,
  species_to,
  geneID_from_IDtype = "symbol",
  geneID_to_IDtype = "symbol",
  assay = NULL,
  layer = "counts",
  multi_mapping = c("first"),
  keep_unmapped = FALSE,
  collapse_fun = c("sum"),
  Ensembl_version = NULL,
  biomart = NULL,
  mirror = NULL,
  max_tries = 5,
  verbose = TRUE
)

# S3 method for class 'matrix'
ConvertHomologs(
  object,
  species_from,
  species_to,
  geneID_from_IDtype = "symbol",
  geneID_to_IDtype = "symbol",
  assay = NULL,
  layer = "counts",
  multi_mapping = c("first"),
  keep_unmapped = FALSE,
  collapse_fun = c("sum"),
  Ensembl_version = NULL,
  biomart = NULL,
  mirror = NULL,
  max_tries = 5,
  verbose = TRUE
)

# S3 method for class 'Matrix'
ConvertHomologs(
  object,
  species_from,
  species_to,
  geneID_from_IDtype = "symbol",
  geneID_to_IDtype = "symbol",
  assay = NULL,
  layer = "counts",
  multi_mapping = c("first"),
  keep_unmapped = FALSE,
  collapse_fun = c("sum"),
  Ensembl_version = NULL,
  biomart = NULL,
  mirror = NULL,
  max_tries = 5,
  verbose = TRUE
)

# Default S3 method
ConvertHomologs(
  object,
  species_from,
  species_to,
  geneID_from_IDtype = "symbol",
  geneID_to_IDtype = "symbol",
  assay = NULL,
  layer = "counts",
  multi_mapping = c("first"),
  keep_unmapped = FALSE,
  collapse_fun = c("sum"),
  Ensembl_version = NULL,
  biomart = NULL,
  mirror = NULL,
  max_tries = 5,
  verbose = TRUE
)
```

## Arguments

- object:

  A `Seurat` object or a gene-by-cell matrix.

- species_from:

  Latin names for animals of the input geneID. e.g. `"Homo_sapiens"`,
  `"Mus_musculus"`.

- species_to:

  Latin names for animals of the output geneID. e.g. `"Homo_sapiens"`,
  `"Mus_musculus"`.

- geneID_from_IDtype:

  Gene ID type of the input `geneID`. e.g. `"symbol"`, `"ensembl_id"`,
  `"entrez_id"`

- geneID_to_IDtype:

  Gene ID type(s) to convert to. e.g. `"symbol"`, `"ensembl_id"`,
  `"entrez_id"`.

- assay:

  Assay to convert when `object` is a `Seurat` object. If `NULL`, the
  default assay is used.

- layer:

  Assay layer used for conversion. Default `"counts"`.

- multi_mapping:

  How to handle source genes mapped to multiple target homologs.
  `"first"` keeps the first target homolog for each source gene.

- keep_unmapped:

  Whether to keep unmapped source genes with their original names.

- collapse_fun:

  Function used to collapse duplicated target homologs. Currently only
  `"sum"` is supported.

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

A converted object of the same high-level type as `object`. The mapping
table is stored in `@tools$ConvertHomologs` for Seurat objects and in
the `"ConvertHomologs"` attribute for matrix inputs.

## See also

[AnnotateFeatures](https://mengxu98.github.io/scop/reference/AnnotateFeatures.md),
ConvertHomologs\]

## Examples

``` r
data(pancreas_sub)
pancreas_human <- ConvertHomologs(
  pancreas_sub,
  species_from = "Mus_musculus",
  species_to = "Homo_sapiens"
)
#> ℹ [2026-06-28 19:53:10] Connect to the Ensembl archives...
#> ℹ [2026-06-28 19:53:12] Using the 116 version of ensembl database...
#> ℹ [2026-06-28 19:53:12] Downloading the ensembl database from https://jun2026.archive.ensembl.org...
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying asia mirror
#> ℹ [2026-06-28 19:53:24] Searching the dataset mmusculus ...
#> ℹ [2026-06-28 19:53:25] Connecting to the dataset mmusculus_gene_ensembl ...
#> ℹ [2026-06-28 19:53:27] Converting the geneIDs...
#> ℹ [2026-06-28 19:53:40] 14803 genes mapped with "ensembl_symbol"
#> ℹ [2026-06-28 19:53:42] 13 genes mapped with "entrez_symbol"
#> ℹ [2026-06-28 19:53:44] 14 genes mapped with "uniprot_symbol"
#> ℹ [2026-06-28 19:53:45] ==============================
#> ℹ                       14830 genes mapped
#> ℹ                       1168 genes unmapped
#> ℹ                       ==============================
#> ℹ [2026-06-28 19:53:52] Converted 13013 source genes to 12873 target homologs
rownames(pancreas_human)[1:5]
#> [1] "C2orf68"  "C11orf58" "C3orf80"  "C8orf33"  "C9orf85" 
```
