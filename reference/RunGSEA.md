# Perform the enrichment analysis (GSEA) on the genes

Perform the enrichment analysis (GSEA) on the genes

## Usage

``` r
RunGSEA(
  srt = NULL,
  group.by = NULL,
  test.use = "wilcox",
  DE_threshold = "p_val_adj < 0.05",
  scoreType = "std",
  geneID = NULL,
  geneScore = NULL,
  geneID_groups = NULL,
  geneID_exclude = NULL,
  IDtype = "symbol",
  result_IDtype = "symbol",
  species = "Homo_sapiens",
  db = "GO_BP",
  db_update = FALSE,
  db_version = "latest",
  db_combine = FALSE,
  convert_species = TRUE,
  Ensembl_version = NULL,
  mirror = NULL,
  features = NULL,
  TERM2GENE = NULL,
  TERM2NAME = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  unlimited_db = c("Chromosome", "GeneType", "TF", "Enzyme", "CSPA"),
  GO_simplify = FALSE,
  GO_simplify_cutoff = "p.adjust < 0.05",
  simplify_method = "Wang",
  simplify_similarityCutoff = 0.7,
  cores = 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object or `SummarizedExperiment` object containing the
  results of differential expression analysis
  ([`RunDEtest()`](https://mengxu98.github.io/scop/reference/RunDEtest.md)).
  If specified, the genes and groups will be extracted from the object
  automatically. If not specified, the `geneID` and `geneID_groups`
  arguments must be provided.

- group.by:

  Metadata column(s) used to color cells.

- test.use:

  Test to be used in differential expression analysis. This argument is
  only used if `srt` is specified.

- DE_threshold:

  Filter condition for differential expression analysis. This argument
  is only used if `srt` is specified.

- scoreType:

  This parameter defines the GSEA score type. Possible options are
  "std", "pos", "neg". By default ("std") the enrichment score is
  computed as in the original GSEA. The "pos" and "neg" score types are
  intended to be used for one-tailed tests (i.e. when one is interested
  only in positive ("pos") or negateive ("neg") enrichment).

- geneID:

  Gene IDs.

- geneScore:

  A numeric vector that specifies the gene scores, for example, the
  log2(fold change) values of gene expression.

- geneID_groups:

  A factor vector specifying the group labels for each gene.

- geneID_exclude:

  Gene IDs to be excluded from the analysis.

- IDtype:

  Type of gene IDs in the `srt` object or `geneID` argument. This
  argument is used to convert the gene IDs to a different type if
  `IDtype` is different from `result_IDtype`.

- result_IDtype:

  Desired type of gene ID to be used in the output. This argument is
  used to convert the gene IDs from `IDtype` to `result_IDtype`.

- species:

  `"Homo_sapiens"` or `"Mus_musculus"`.

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

- db_combine:

  Whether to combine multiple databases into one. If `TRUE`, all
  database specified by `db` will be combined as one named "Combined".

- convert_species:

  Use a species-converted database when the annotation is missing for
  `species`.

- Ensembl_version:

  Ensembl version. `NULL` uses the latest.

- mirror:

  Specify an Ensembl mirror to connect to. The valid options here are
  `"www"`, `"uswest"`, `"useast"`, `"asia"`.

- features:

  A named list of feature lists for custom enrichment gene sets. If
  provided, it takes precedence over `TERM2GENE` and `db`.

- TERM2GENE:

  A data frame specifying the gene-term mapping for a custom database.
  The first column should contain the term IDs, and the second column
  should contain the gene IDs.

- TERM2NAME:

  A data frame specifying the term-name mapping for a custom database.
  The first column should contain the term IDs, and the second column
  should contain the corresponding term names.

- minGSSize:

  The minimum size of a gene set to be considered in the enrichment
  analysis.

- maxGSSize:

  The maximum size of a gene set to be considered in the enrichment
  analysis.

- unlimited_db:

  Names of databases that do not have size restrictions.

- GO_simplify:

  Whether to simplify the GO terms. If `TRUE`, additional results with
  simplified GO terms will be returned.

- GO_simplify_cutoff:

  Filter condition for simplification of GO terms. This argument is only
  used if `GO_simplify` is `TRUE`.

- simplify_method:

  Method to be used for simplification of GO terms. This argument is
  only used if `GO_simplify` is `TRUE`.

- simplify_similarityCutoff:

  The similarity cutoff for simplification of GO terms. This argument is
  only used if `GO_simplify` is `TRUE`.

- cores:

  The number of worker processes to use for parallelization. Default is
  `1`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Passed to helper functions.

## Value

If input is a Seurat object, returns the modified Seurat object with the
enrichment result stored in the tools slot. If input is a geneID vector
with or without geneID_groups, return the enrichment result directly.
Enrichment result is a list with the following component:

- `enrichment`: A data.frame containing all enrichment results.

- `results`: A list of `gseaResult` objects from the DOSE package.

- `geneMap`: A data.frame containing the ID mapping table for input gene
  IDs.

- `input`: A data.frame containing the input gene IDs and gene ID
  groups.

- `DE_threshold`: A specific threshold for differential expression
  analysis (only returned if input is a Seurat object).

## See also

[PrepareDB](https://mengxu98.github.io/scop/reference/PrepareDB.md),
[ListDB](https://mengxu98.github.io/scop/reference/ListDB.md),
[GSEAPlot](https://mengxu98.github.io/scop/reference/GSEAPlot.md),
[RunEnrichment](https://mengxu98.github.io/scop/reference/RunEnrichment.md),
[EnrichmentPlot](https://mengxu98.github.io/scop/reference/EnrichmentPlot.md)

## Examples

``` r
term2gene <- data.frame(
  Term = c(
    rep("Endocrine markers", 5),
    rep("Exocrine markers", 5),
    rep("Ductal markers", 5)
  ),
  symbol = c(
    "INS", "GCG", "SST", "IAPP", "PCSK1",
    "PRSS1", "CPA1", "CELA3A", "REG1A", "CTRB1",
    "KRT19", "SOX9", "MUC1", "CFTR", "KRT7"
  )
)
genes <- unique(term2gene$symbol)
scores <- seq_along(genes)
names(scores) <- genes
gsea_out <- RunGSEA(
  geneID = genes,
  geneScore = scores,
  geneID_groups = rep("Cluster1", length(genes)),
  TERM2GENE = term2gene,
  minGSSize = 2
)
#> ℹ [2026-09-06 22:09:02] Start GSEA analysis
#> ! [2026-09-06 22:09:02] All values in the `geneScore` are greater than zero. Set scoreType = 'pos'
#> ℹ [2026-09-06 22:09:02] Preparing custom GSEA database from `TERM2GENE`/`TERM2NAME` ...
#> ℹ [2026-09-06 22:09:02] Prepared 15 ranked gene rows after ID mapping.
#> ℹ [2026-09-06 22:09:02] Running GSEA for 1 group/database combination(s) using 1 core(s) ...
#> ℹ [2026-09-06 22:09:02] Using 1 core
#> ℹ [2026-09-06 22:09:02] Running GSEA: group "Cluster1", database "custom", genes 15 ...
#> ℹ [2026-09-06 22:09:02] Finished GSEA: group "Cluster1", database "custom".
#> ℹ [2026-09-06 22:09:02] Building results
#> ✔ [2026-09-06 22:09:02] GSEA analysis done
GSEAPlot(
  res = gsea_out,
  db = "custom",
  plot_type = "comparison"
)
```
