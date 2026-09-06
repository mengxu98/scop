# Annotate Features

Annotate features in a Seurat object with additional metadata from
databases or a GTF file.

## Usage

``` r
AnnotateFeatures(
  srt,
  species = "Homo_sapiens",
  IDtype = c("symbol", "ensembl_id", "entrez_id"),
  db = NULL,
  db_update = FALSE,
  db_version = "latest",
  convert_species = TRUE,
  Ensembl_version = NULL,
  mirror = NULL,
  gtf = NULL,
  merge_gtf_by = "gene_name",
  columns = c("seqname", "feature", "start", "end", "strand", "gene_id", "gene_name",
    "gene_type"),
  assays = "RNA",
  overwrite = FALSE,
  ...,
  verbose = TRUE
)
```

## Arguments

- srt:

  Seurat object to be annotated.

- species:

  `"Homo_sapiens"` or `"Mus_musculus"`.

- IDtype:

  Type of identifier to use for annotation. Options are `"symbol"`,
  `"ensembl_id"`, or `"entrez_id"`.

- db:

  Annotation sources. One or more of `"GO"`, `"GO_BP"`, `"GO_CC"`,
  `"GO_MF"`, `"KEGG"`, `"WikiPathway"`, `"Reactome"`, `"CORUM"`, `"MP"`,
  `"DO"`, `"HPO"`, `"PFAM"`, `"CSPA"`, `"Surfaceome"`, `"SPRomeDB"`,
  `"VerSeDa"`, `"TFLink"`, `"hTFtarget"`, `"TRRUST"`, `"JASPAR"`,
  `"ENCODE"`, `"MSigDB"`, `"CellTalk"`, `"CellChat"`, `"Chromosome"`,
  `"GeneType"`, `"Enzyme"`, `"TF"`, `"CytoTRACE2"`. MSigDB
  subcollections use `"MSigDB_<collection>"` (e.g. `"MSigDB_H"`).
  `"CytoTRACE2"` is species-independent and is required by
  [RunCytoTRACE](https://mengxu98.github.io/scop/reference/RunCytoTRACE.md).

- db_update:

  Force a refresh. `FALSE` loads the cache when available.

- db_version:

  Database version to retrieve.

- convert_species:

  Use a species-converted database when the annotation is missing for
  `species`.

- Ensembl_version:

  Ensembl version. `NULL` uses the latest.

- mirror:

  URL of the mirror to use for Ensembl database.

- gtf:

  Path to the GTF file to be used for annotation.

- merge_gtf_by:

  Column name to merge the GTF file by.

- columns:

  Vector of column names to be used from the GTF file. Default is
  `"seqname"`, `"feature"`, `"start"`, `"end"`, `"strand"`, `"gene_id"`,
  `"gene_name"`, `"gene_type"`.

- assays:

  Character vector of assay names to be annotated.

- overwrite:

  Whether to overwrite existing metadata.

- ...:

  Passed to helper functions.

- verbose:

  Whether to print the message. Default is `TRUE`.

## See also

[PrepareDB](https://mengxu98.github.io/scop/reference/PrepareDB.md),
[ListDB](https://mengxu98.github.io/scop/reference/ListDB.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- AnnotateFeatures(
  pancreas_sub,
  species = "Mus_musculus",
  db = "TF"
)
#> ℹ [2026-09-06 20:52:04] Species: "Mus_musculus"
#> ℹ [2026-09-06 20:52:04] Preparing database: TF
head(
  GetFeaturesData(
    pancreas_sub
  )
)
#>               highly_variable_genes   TF
#> Xkr4                          False <NA>
#> Mrpl15                        False <NA>
#> Npbwr1                         <NA> <NA>
#> 4732440D04Rik                 False <NA>
#> Gm26901                       False <NA>
#> Sntg1                          True <NA>

# Annotate features using a GTF file
gtf_file <- "/refdata-gex-mm10-2020-A/genes/genes.gtf"
if (file.exists(gtf_file)) {
  pancreas_sub <- AnnotateFeatures(
    pancreas_sub,
    gtf = gtf_file
  )
  head(
    GetFeaturesData(
      pancreas_sub
    )
  )
}
```
