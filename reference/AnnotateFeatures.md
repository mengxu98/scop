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
  overwrite = FALSE
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
  `"GO", "GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway", "Reactome", "CORUM", "MP", "DO", "HPO", "PFAM", "CSPA", "Surfaceome", "SPRomeDB", "VerSeDa", "TFLink", "hTFtarget", "TRRUST", "JASPAR", "ENCODE", "MSigDB", "CellTalk", "CellChat", "Chromosome", "GeneType", "Enzyme", "TF"`.

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

## See also

[PrepareDB](https://mengxu98.github.io/scop/reference/PrepareDB.md),
[ListDB](https://mengxu98.github.io/scop/reference/ListDB.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- AnnotateFeatures(
  pancreas_sub,
  species = "Mus_musculus",
  db = c(
    "Chromosome",
    "GeneType",
    "Enzyme",
    "TF",
    "CSPA",
    "VerSeDa"
  )
)
head(
  GetFeaturesData(
    pancreas_sub
  )
)

# Annotate features using a GTF file
pancreas_sub <- AnnotateFeatures(
  pancreas_sub,
  gtf = "/refdata-gex-mm10-2020-A/genes/genes.gtf"
)
head(
  GetFeaturesData(
    pancreas_sub
  )
)
} # }
```
