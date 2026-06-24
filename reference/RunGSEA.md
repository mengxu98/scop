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

  Name of one or more meta.data columns to group (color) cells by.

- test.use:

  A character vector specifying the test to be used in differential
  expression analysis. This argument is only used if `srt` is specified.

- DE_threshold:

  A character vector specifying the filter condition for differential
  expression analysis. This argument is only used if `srt` is specified.

- scoreType:

  This parameter defines the GSEA score type. Possible options are
  "std", "pos", "neg". By default ("std") the enrichment score is
  computed as in the original GSEA. The "pos" and "neg" score types are
  intended to be used for one-tailed tests (i.e. when one is interested
  only in positive ("pos") or negateive ("neg") enrichment).

- geneID:

  A character vector specifying the gene IDs.

- geneScore:

  A numeric vector that specifies the gene scores, for example, the
  log2(fold change) values of gene expression.

- geneID_groups:

  A factor vector specifying the group labels for each gene.

- geneID_exclude:

  A character vector specifying the gene IDs to be excluded from the
  analysis.

- IDtype:

  A character vector specifying the type of gene IDs in the `srt` object
  or `geneID` argument. This argument is used to convert the gene IDs to
  a different type if `IDtype` is different from `result_IDtype`.

- result_IDtype:

  A character vector specifying the desired type of gene ID to be used
  in the output. This argument is used to convert the gene IDs from
  `IDtype` to `result_IDtype`.

- species:

  A character vector specifying the species for which the gene
  annotation databases should be prepared. Can be `"Homo_sapiens"` or
  `"Mus_musculus"`.

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

- db_combine:

  Whether to combine multiple databases into one. If `TRUE`, all
  database specified by `db` will be combined as one named "Combined".

- convert_species:

  Whether to use a species-converted database when the annotation is
  missing for the specified species. Default is `TRUE`.

- Ensembl_version:

  An integer specifying the Ensembl version. Default is `NULL`. If
  `NULL`, the latest version will be used.

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

  A character vector specifying the names of databases that do not have
  size restrictions.

- GO_simplify:

  Whether to simplify the GO terms. If `TRUE`, additional results with
  simplified GO terms will be returned.

- GO_simplify_cutoff:

  A character vector specifying the filter condition for simplification
  of GO terms. This argument is only used if `GO_simplify` is `TRUE`.

- simplify_method:

  A character vector specifying the method to be used for simplification
  of GO terms. This argument is only used if `GO_simplify` is `TRUE`.

- simplify_similarityCutoff:

  The similarity cutoff for simplification of GO terms. This argument is
  only used if `GO_simplify` is `TRUE`.

- cores:

  The number of cores to use for parallelization with
  [foreach::foreach](https://rdrr.io/pkg/foreach/man/foreach.html).
  Default is `1`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Passed to other functions.

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
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-06-24 04:21:27] Start standard processing workflow...
#> ℹ [2026-06-24 04:21:28] Checking a list of <Seurat>...
#> ! [2026-06-24 04:21:28] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-06-24 04:21:28] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-06-24 04:21:28] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-06-24 04:21:28] Use the separate HVF from `srt_list`
#> ℹ [2026-06-24 04:21:28] Number of available HVF: 2000
#> ℹ [2026-06-24 04:21:28] Finished check
#> ℹ [2026-06-24 04:21:28] Perform `ScaleData()`
#> ℹ [2026-06-24 04:21:28] Perform pca linear dimension reduction
#> ℹ [2026-06-24 04:21:29] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-06-24 04:21:29] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-06-24 04:21:29] Reorder clusters...
#> ℹ [2026-06-24 04:21:29] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-06-24 04:21:29] Perform umap nonlinear dimension reduction
#> ✔ [2026-06-24 04:21:37] Standard processing workflow completed
pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "CellType"
)
#> ℹ [2026-06-24 04:21:37] Data type is log-normalized
#> ℹ [2026-06-24 04:21:37] Start differential expression test
#> ℹ [2026-06-24 04:21:37] Find all markers(wilcox) among [1] 5 groups...
#> ℹ [2026-06-24 04:21:37] Using 1 core
#> ⠙ [2026-06-24 04:21:37] Running for Ductal [1/5] ■■          20% | ETA:  1s
#> ✔ [2026-06-24 04:21:37] Completed 5 tasks in 818ms
#> 
#> ℹ [2026-06-24 04:21:37] Building results
#> ✔ [2026-06-24 04:21:38] Differential expression test completed
pancreas_sub <- RunGSEA(
  pancreas_sub,
  group.by = "CellType",
  DE_threshold = "p_val_adj < 0.05",
  scoreType = "std",
  db = "GO_BP",
  species = "Mus_musculus"
)
#> ℹ [2026-06-24 04:21:38] Start GSEA analysis
#> ℹ [2026-06-24 04:21:38] Resolving `RunDEtest()` results for `group.by` = "CellType" and `test.use` = "wilcox" ...
#> ℹ [2026-06-24 04:21:38] Filtering DE results with `DE_threshold`: "p_val_adj < 0.05" ...
#> ℹ [2026-06-24 04:21:38] Using 3754 ranked genes across 5 group(s) for GSEA.
#> ! [2026-06-24 04:21:38] All values in the `geneScore` are greater than zero. Set scoreType = 'pos'
#> ℹ [2026-06-24 04:21:38] Preparing GSEA database(s): "GO_BP" for species "Mus_musculus" ...
#> ℹ [2026-06-24 04:21:38] Species: "Mus_musculus"
#> ℹ [2026-06-24 04:21:38] Loading cached: GO_BP version: 3.23.0 nterm:14957 created: 2026-06-24 03:23:47
#> ℹ [2026-06-24 04:21:39] Prepared 3754 ranked gene rows after ID mapping.
#> ℹ [2026-06-24 04:21:39] Running GSEA for 5 group/database combination(s) using 1 core(s) ...
#> ℹ [2026-06-24 04:21:39] Using 1 core
#> ℹ [2026-06-24 04:21:39] Running GSEA: group "Ductal", database "GO_BP", genes 1557 ...
#> ℹ [2026-06-24 04:21:48] Finished GSEA: group "Ductal", database "GO_BP".
#> ⠙ [2026-06-24 04:21:39] Running for 1 [1/5] ■■          20% | ETA: 33s
#> ℹ [2026-06-24 04:21:48] Running GSEA: group "Ngn3-high-EP", database "GO_BP", genes 482 ...
#> ⠙ [2026-06-24 04:21:39] Running for 1 [1/5] ■■          20% | ETA: 33s
#> ℹ [2026-06-24 04:21:51] Finished GSEA: group "Ngn3-high-EP", database "GO_BP".
#> ⠙ [2026-06-24 04:21:39] Running for 1 [1/5] ■■          20% | ETA: 33s
#> ⠹ [2026-06-24 04:21:39] Running for 2 [2/5] ■■■■        40% | ETA: 17s
#> ℹ [2026-06-24 04:21:51] Running GSEA: group "Endocrine", database "GO_BP", genes 912 ...
#> ⠹ [2026-06-24 04:21:39] Running for 2 [2/5] ■■■■        40% | ETA: 17s
#> ℹ [2026-06-24 04:21:55] Finished GSEA: group "Endocrine", database "GO_BP".
#> ⠹ [2026-06-24 04:21:39] Running for 2 [2/5] ■■■■        40% | ETA: 17s
#> ⠸ [2026-06-24 04:21:39] Running for 3 [3/5] ■■■■■■      60% | ETA: 11s
#> ℹ [2026-06-24 04:21:55] Running GSEA: group "Ngn3-low-EP", database "GO_BP", genes 423 ...
#> ⠸ [2026-06-24 04:21:39] Running for 3 [3/5] ■■■■■■      60% | ETA: 11s
#> ℹ [2026-06-24 04:21:59] Finished GSEA: group "Ngn3-low-EP", database "GO_BP".
#> ⠸ [2026-06-24 04:21:39] Running for 3 [3/5] ■■■■■■      60% | ETA: 11s
#> ⠼ [2026-06-24 04:21:39] Running for 4 [4/5] ■■■■■■■■    80% | ETA:  5s
#> ℹ [2026-06-24 04:21:59] Running GSEA: group "Pre-endocrine", database "GO_BP", genes 380 ...
#> ⠼ [2026-06-24 04:21:39] Running for 4 [4/5] ■■■■■■■■    80% | ETA:  5s
#> ℹ [2026-06-24 04:22:02] Finished GSEA: group "Pre-endocrine", database "GO_BP".
#> ⠼ [2026-06-24 04:21:39] Running for 4 [4/5] ■■■■■■■■    80% | ETA:  5s
#> ✔ [2026-06-24 04:21:39] Completed 5 tasks in 22.5s
#> 
#> ℹ [2026-06-24 04:21:39] Building results
#> ✔ [2026-06-24 04:22:02] GSEA analysis done
GSEAPlot(
  pancreas_sub,
  db = "GO_BP",
  group.by = "CellType",
  plot_type = "comparison"
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's alpha values.
#> Warning: Removed 9786 rows containing missing values or values outside the scale range
#> (`geom_point()`).

GSEAPlot(
  pancreas_sub,
  db = "GO_BP",
  group.by = "CellType",
  group_use = "Ductal",
  id_use = "GO:0006412"
)
#> Error in `.rowNamesDF<-`(x, value = value): missing values in 'row.names' are not allowed
GSEAPlot(
  pancreas_sub,
  db = "GO_BP",
  group.by = "CellType",
  group_use = "Ductal",
  id_use = c(
    "GO:0046903", "GO:0015031", "GO:0007600"
  )
)
#> Warning: non-unique values when setting 'row.names': 
#> Error in `.rowNamesDF<-`(x, value = value): duplicate 'row.names' are not allowed

# Remove redundant GO terms
pancreas_sub <- RunGSEA(
  pancreas_sub,
  group.by = "CellType",
  db = "GO_BP",
  GO_simplify = TRUE,
  species = "Mus_musculus"
)
#> ℹ [2026-06-24 04:22:03] Start GSEA analysis
#> ℹ [2026-06-24 04:22:03] Resolving `RunDEtest()` results for `group.by` = "CellType" and `test.use` = "wilcox" ...
#> ℹ [2026-06-24 04:22:03] Filtering DE results with `DE_threshold`: "p_val_adj < 0.05" ...
#> ℹ [2026-06-24 04:22:03] Using 3754 ranked genes across 5 group(s) for GSEA.
#> ! [2026-06-24 04:22:03] All values in the `geneScore` are greater than zero. Set scoreType = 'pos'
#> ℹ [2026-06-24 04:22:03] Preparing GSEA database(s): "GO_BP" for species "Mus_musculus" ...
#> ℹ [2026-06-24 04:22:03] Species: "Mus_musculus"
#> ℹ [2026-06-24 04:22:03] Loading cached: GO_BP version: 3.23.0 nterm:14957 created: 2026-06-24 03:23:47
#> ℹ [2026-06-24 04:22:04] Prepared 3754 ranked gene rows after ID mapping.
#> ℹ [2026-06-24 04:22:04] Running GSEA for 5 group/database combination(s) using 1 core(s) ...
#> ℹ [2026-06-24 04:22:04] Using 1 core
#> ℹ [2026-06-24 04:22:04] Running GSEA: group "Ductal", database "GO_BP", genes 1557 ...
#> ⠙ [2026-06-24 04:22:04] Running for 1 [1/5] ■■          20% | ETA: 33s
#> ℹ [2026-06-24 04:22:12] Running GSEA: group "Ngn3-high-EP", database "GO_BP", genes 482 ...
#> ⠙ [2026-06-24 04:22:04] Running for 1 [1/5] ■■          20% | ETA: 33s
#> ⠹ [2026-06-24 04:22:04] Running for 2 [2/5] ■■■■        40% | ETA: 17s
#> ℹ [2026-06-24 04:22:15] Running GSEA: group "Endocrine", database "GO_BP", genes 912 ...
#> ⠹ [2026-06-24 04:22:04] Running for 2 [2/5] ■■■■        40% | ETA: 17s
#> ⠸ [2026-06-24 04:22:04] Running for 3 [3/5] ■■■■■■      60% | ETA: 11s
#> ℹ [2026-06-24 04:22:20] Running GSEA: group "Ngn3-low-EP", database "GO_BP", genes 423 ...
#> ⠸ [2026-06-24 04:22:04] Running for 3 [3/5] ■■■■■■      60% | ETA: 11s
#> ⠼ [2026-06-24 04:22:04] Running for 4 [4/5] ■■■■■■■■    80% | ETA:  5s
#> ℹ [2026-06-24 04:22:23] Running GSEA: group "Pre-endocrine", database "GO_BP", genes 380 ...
#> ⠼ [2026-06-24 04:22:04] Running for 4 [4/5] ■■■■■■■■    80% | ETA:  5s
#> ✔ [2026-06-24 04:22:04] Completed 5 tasks in 22.4s
#> 
#> ℹ [2026-06-24 04:22:04] Building results
#> ! [2026-06-24 04:22:04] Found 5 failed results
#> ℹ [2026-06-24 04:22:26] ✖ Error details:
#> ℹ                       ✖ missing value where TRUE/FALSE needed (5): "1", "2", "3" and 2 more
#> Error in x@result: no applicable method for `@` applied to an object of class "parallelize_error"
GSEAPlot(
  pancreas_sub,
  db = "GO_BP_sim",
  group.by = "CellType",
  plot_type = "comparison"
)
#> Error in resolve_enrichment_plot_db(db = db, enrichment = enrichment): GO_BP_sim is not in the enrichment result

# Or use "geneID", "geneScore" and
# "geneID_groups" as input to run GSEA
de_df <- dplyr::filter(
  pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox,
  p_val_adj < 0.05
)
gsea_out <- RunGSEA(
  geneID = de_df[["gene"]],
  geneScore = de_df[["avg_log2FC"]],
  geneID_groups = de_df[["group1"]],
  db = "GO_BP",
  species = "Mus_musculus"
)
#> ℹ [2026-06-24 04:22:26] Start GSEA analysis
#> ! [2026-06-24 04:22:26] All values in the `geneScore` are greater than zero. Set scoreType = 'pos'
#> ℹ [2026-06-24 04:22:26] Preparing GSEA database(s): "GO_BP" for species "Mus_musculus" ...
#> ℹ [2026-06-24 04:22:26] Species: "Mus_musculus"
#> ℹ [2026-06-24 04:22:26] Loading cached: GO_BP version: 3.23.0 nterm:14957 created: 2026-06-24 03:23:47
#> ℹ [2026-06-24 04:22:28] Prepared 3754 ranked gene rows after ID mapping.
#> ℹ [2026-06-24 04:22:28] Running GSEA for 5 group/database combination(s) using 1 core(s) ...
#> ℹ [2026-06-24 04:22:28] Using 1 core
#> ℹ [2026-06-24 04:22:28] Running GSEA: group "Ductal", database "GO_BP", genes 1557 ...
#> ℹ [2026-06-24 04:22:36] Finished GSEA: group "Ductal", database "GO_BP".
#> ⠙ [2026-06-24 04:22:28] Running for 1 [1/5] ■■          20% | ETA: 33s
#> ℹ [2026-06-24 04:22:36] Running GSEA: group "Ngn3-high-EP", database "GO_BP", genes 482 ...
#> ⠙ [2026-06-24 04:22:28] Running for 1 [1/5] ■■          20% | ETA: 33s
#> ℹ [2026-06-24 04:22:39] Finished GSEA: group "Ngn3-high-EP", database "GO_BP".
#> ⠙ [2026-06-24 04:22:28] Running for 1 [1/5] ■■          20% | ETA: 33s
#> ⠹ [2026-06-24 04:22:28] Running for 2 [2/5] ■■■■        40% | ETA: 17s
#> ℹ [2026-06-24 04:22:39] Running GSEA: group "Endocrine", database "GO_BP", genes 912 ...
#> ⠹ [2026-06-24 04:22:28] Running for 2 [2/5] ■■■■        40% | ETA: 17s
#> ℹ [2026-06-24 04:22:44] Finished GSEA: group "Endocrine", database "GO_BP".
#> ⠹ [2026-06-24 04:22:28] Running for 2 [2/5] ■■■■        40% | ETA: 17s
#> ⠸ [2026-06-24 04:22:28] Running for 3 [3/5] ■■■■■■      60% | ETA: 11s
#> ℹ [2026-06-24 04:22:44] Running GSEA: group "Ngn3-low-EP", database "GO_BP", genes 423 ...
#> ⠸ [2026-06-24 04:22:28] Running for 3 [3/5] ■■■■■■      60% | ETA: 11s
#> ℹ [2026-06-24 04:22:47] Finished GSEA: group "Ngn3-low-EP", database "GO_BP".
#> ⠸ [2026-06-24 04:22:28] Running for 3 [3/5] ■■■■■■      60% | ETA: 11s
#> ⠼ [2026-06-24 04:22:28] Running for 4 [4/5] ■■■■■■■■    80% | ETA:  5s
#> ℹ [2026-06-24 04:22:47] Running GSEA: group "Pre-endocrine", database "GO_BP", genes 380 ...
#> ⠼ [2026-06-24 04:22:28] Running for 4 [4/5] ■■■■■■■■    80% | ETA:  5s
#> ℹ [2026-06-24 04:22:50] Finished GSEA: group "Pre-endocrine", database "GO_BP".
#> ⠼ [2026-06-24 04:22:28] Running for 4 [4/5] ■■■■■■■■    80% | ETA:  5s
#> ✔ [2026-06-24 04:22:28] Completed 5 tasks in 22.5s
#> 
#> ℹ [2026-06-24 04:22:28] Building results
#> ✔ [2026-06-24 04:22:50] GSEA analysis done
GSEAPlot(
  res = gsea_out,
  db = "GO_BP",
  plot_type = "comparison"
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's alpha values.
#> Warning: Removed 9786 rows containing missing values or values outside the scale range
#> (`geom_point()`).


# Use a combined database
pancreas_sub <- RunGSEA(
  pancreas_sub,
  group.by = "CellType",
  db = c(
    "KEGG", "WikiPathway", "Reactome", "PFAM", "MP"
  ),
  db_combine = TRUE,
  species = "Mus_musculus"
)
#> ℹ [2026-06-24 04:22:51] Start GSEA analysis
#> ℹ [2026-06-24 04:22:51] Resolving `RunDEtest()` results for `group.by` = "CellType" and `test.use` = "wilcox" ...
#> ℹ [2026-06-24 04:22:51] Filtering DE results with `DE_threshold`: "p_val_adj < 0.05" ...
#> ℹ [2026-06-24 04:22:51] Using 3754 ranked genes across 5 group(s) for GSEA.
#> ! [2026-06-24 04:22:51] All values in the `geneScore` are greater than zero. Set scoreType = 'pos'
#> ℹ [2026-06-24 04:22:51] Preparing GSEA database(s): "KEGG, WikiPathway, Reactome, PFAM, MP" for species "Mus_musculus" ...
#> ℹ [2026-06-24 04:22:51] Species: "Mus_musculus"
#> ℹ [2026-06-24 04:22:51] Preparing KEGG database
#> ℹ [2026-06-24 04:22:52] <simpleError in utils::download.file(url = url, destfile = destfile, method = method,     quiet = quiet, ...): cannot open URL 'https://rest.kegg.jp/list/organism'>
#> ! [2026-06-24 04:22:52] Failed to download using auto, from <https://rest.kegg.jp/list/organism>
#> ℹ [2026-06-24 04:22:55] <simpleError in utils::download.file(url = url, destfile = destfile, method = method,     quiet = quiet, ...): 'wget' call had nonzero exit status>
#> ! [2026-06-24 04:22:55] Failed to download using wget, from <https://rest.kegg.jp/list/organism>
#> ℹ [2026-06-24 04:22:58] <simpleError in utils::download.file(url = url, destfile = destfile, method = method,     quiet = quiet, ...): cannot open URL 'https://rest.kegg.jp/list/organism'>
#> ! [2026-06-24 04:22:58] Failed to download using libcurl, from <https://rest.kegg.jp/list/organism>
#> Error in `[.data.frame`(orgs, , 3): undefined columns selected
GSEAPlot(
  pancreas_sub,
  db = "Combined",
  group.by = "CellType",
  plot_type = "comparison"
)
#> Error in resolve_enrichment_plot_db(db = db, enrichment = enrichment): Combined is not in the enrichment result
```
