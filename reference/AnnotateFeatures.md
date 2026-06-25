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

  A character vector specifying the species for which the gene
  annotation databases should be prepared. Can be `"Homo_sapiens"` or
  `"Mus_musculus"`.

- IDtype:

  Type of identifier to use for annotation. Options are `"symbol"`,
  `"ensembl_id"`, or `"entrez_id"`. Default is `"symbol"`.

- db:

  A character vector specifying the annotation sources to be included in
  the gene annotation databases. Can be one or more of
  `"GO", "GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway", "Reactome", "CORUM", "MP", "DO", "HPO", "PFAM", "CSPA", "Surfaceome", "SPRomeDB", "VerSeDa", "TFLink", "hTFtarget", "TRRUST", "JASPAR", "ENCODE", "MSigDB", "CellTalk", "CellChat", "Chromosome", "GeneType", "Enzyme", "TF", "CytoTRACE2"`.
  MSigDB subcollections can be requested as `"MSigDB_<collection>"`,
  such as `"MSigDB_H"` for human Hallmark and `"MSigDB_MH"` for mouse
  Hallmark. Note: `"CytoTRACE2"` is species-independent and downloads
  pre-trained model data required by
  [RunCytoTRACE](https://mengxu98.github.io/scop/reference/RunCytoTRACE.md).

- db_update:

  Whether the gene annotation databases should be forcefully updated. If
  set to FALSE, the function will attempt to load the cached databases
  instead. Default is `FALSE`.

- db_version:

  A character vector specifying the version of the gene annotation
  databases to be retrieved. Default is `"latest"`.

- convert_species:

  Whether to use a species-converted database when the annotation is
  missing for the specified species. Default is `TRUE`.

- Ensembl_version:

  An integer specifying the Ensembl version. Default is `NULL`. If
  `NULL`, the latest version will be used.

- mirror:

  URL of the mirror to use for Ensembl database. Default is `NULL`.

- gtf:

  Path to the GTF file to be used for annotation. Default is `NULL`.

- merge_gtf_by:

  Column name to merge the GTF file by. Default is `"gene_name"`.

- columns:

  Vector of column names to be used from the GTF file. Default is
  `"seqname"`, `"feature"`, `"start"`, `"end"`, `"strand"`, `"gene_id"`,
  `"gene_name"`, `"gene_type"`.

- assays:

  Character vector of assay names to be annotated. Default is `"RNA"`.

- overwrite:

  Whether to overwrite existing metadata. Default is `FALSE`.

- ...:

  Passed to other functions.

- verbose:

  Whether to print the message. Default is `TRUE`.

## See also

[PrepareDB](https://mengxu98.github.io/scop/reference/PrepareDB.md),
[ListDB](https://mengxu98.github.io/scop/reference/ListDB.md)

## Examples

``` r
if (requireNamespace("R.cache", quietly = TRUE)) {
  data(pancreas_sub)
  pancreas_sub <- AnnotateFeatures(
    pancreas_sub,
    species = "Mus_musculus",
    db = "TF"
  )
  head(
    GetFeaturesData(
      pancreas_sub
    )
  )
}
#> ℹ [2026-06-25 06:29:26] Species: "Mus_musculus"
#> ℹ [2026-06-25 06:29:26] Preparing database: TF
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
