# Perform the enrichment analysis (over-representation) on the genes

Perform the enrichment analysis (over-representation) on the genes

## Usage

``` r
RunEnrichment(
  srt = NULL,
  group.by = NULL,
  test.use = "wilcox",
  DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
  geneID = NULL,
  geneID_groups = NULL,
  geneID_exclude = NULL,
  IDtype = "symbol",
  result_IDtype = "symbol",
  backend = c("cpp", "r"),
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
  verbose = TRUE
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

- geneID:

  A character vector specifying the gene IDs.

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

- backend:

  Enrichment backend. `"cpp"` is the default and uses a fast native
  hypergeometric ORA implementation and returns the enrichment table
  without `enrichResult` objects. `"r"` uses
  [`clusterProfiler::enricher()`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html)
  and returns `enrichResult` objects in `results`. `GO_simplify = TRUE`
  currently uses the R backend.

- species:

  A character vector specifying the species for which the gene
  annotation databases should be prepared. Can be `"Homo_sapiens"` or
  `"Mus_musculus"`.

- db:

  A character vector specifying the annotation sources to be included in
  the gene annotation databases. Can be one or more of
  `"GO", "GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway", "Reactome", "CORUM", "MP", "DO", "HPO", "PFAM", "CSPA", "Surfaceome", "SPRomeDB", "VerSeDa", "TFLink", "hTFtarget", "TRRUST", "JASPAR", "ENCODE", "MSigDB", "CellTalk", "CellChat", "Chromosome", "GeneType", "Enzyme", "TF", "CytoTRACE2"`.
  Note: `"CytoTRACE2"` is species-independent and downloads pre-trained
  model data required by
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

## Value

If input is a Seurat object, returns the modified Seurat object with the
enrichment result stored in the tools slot.

If input is a geneID vector with or without geneID_groups, return the
enrichment result directly.

Enrichment result is a list with the following component:

- `enrichment`: A data.frame containing all enrichment results.

- `results`: A list of `enrichResult` objects from the DOSE package.

- `geneMap`: A data.frame containing the ID mapping table for input gene
  IDs.

- `input`: A data.frame containing the input gene IDs and gene ID
  groups.

- `DE_threshold`: A specific threshold for differential expression
  analysis (only returned if input is a Seurat object).

## See also

[PrepareDB](https://mengxu98.github.io/scop/reference/PrepareDB.md),
[ListDB](https://mengxu98.github.io/scop/reference/ListDB.md),
[EnrichmentPlot](https://mengxu98.github.io/scop/reference/EnrichmentPlot.md),
[RunGSEA](https://mengxu98.github.io/scop/reference/RunGSEA.md),
[GSEAPlot](https://mengxu98.github.io/scop/reference/GSEAPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-05-23 14:30:36] Start standard processing workflow...
#> ℹ [2026-05-23 14:30:36] Checking a list of <Seurat>...
#> ! [2026-05-23 14:30:37] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-05-23 14:30:37] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-05-23 14:30:38] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-23 14:30:39] Use the separate HVF from `srt_list`
#> ℹ [2026-05-23 14:30:39] Number of available HVF: 2000
#> ℹ [2026-05-23 14:30:39] Finished check
#> ℹ [2026-05-23 14:30:39] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-23 14:30:39] Perform pca linear dimension reduction
#> ℹ [2026-05-23 14:30:40] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-05-23 14:30:40] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-23 14:30:40] Reorder clusters...
#> ℹ [2026-05-23 14:30:40] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-23 14:30:40] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-23 14:30:41] Perform umap nonlinear dimension reduction using Standardpca (1:23)
#> ℹ [2026-05-23 14:30:46] Perform umap nonlinear dimension reduction using Standardpca (1:23)
#> ✔ [2026-05-23 14:30:51] Standard processing workflow completed
pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "CellType"
)
#> ℹ [2026-05-23 14:30:51] Data type is log-normalized
#> ℹ [2026-05-23 14:30:51] Start differential expression test
#> ℹ [2026-05-23 14:30:51] Find all markers(wilcox) among [1] 5 groups...
#> ℹ [2026-05-23 14:30:51] Using 1 core
#> ⠙ [2026-05-23 14:30:51] Running for Ductal [1/5] ■■          20% | ETA:  1s
#> ✔ [2026-05-23 14:30:51] Completed 5 tasks in 858ms
#> 
#> ℹ [2026-05-23 14:30:51] Building results
#> ✔ [2026-05-23 14:30:52] Differential expression test completed
pancreas_sub <- RunEnrichment(
  pancreas_sub,
  group.by = "CellType",
  DE_threshold = "p_val_adj < 0.05",
  db = "GO_BP",
  species = "Mus_musculus"
)
#> ℹ [2026-05-23 14:30:52] Start Enrichment analysis
#> ℹ [2026-05-23 14:30:52] Species: "Mus_musculus"
#> ℹ [2026-05-23 14:30:52] Loading cached: GO_BP version: 3.23.0 nterm:14957 created: 2026-05-23 13:33:03
#> ℹ [2026-05-23 14:30:53] Permform enrichment...
#> ℹ [2026-05-23 14:30:53] Using 1 core
#> ⠙ [2026-05-23 14:30:53] Running for 1 [1/5] ■■          20% | ETA: 22s
#> ⠹ [2026-05-23 14:30:53] Running for 2 [2/5] ■■■■        40% | ETA: 13s
#> ⠸ [2026-05-23 14:30:53] Running for 3 [3/5] ■■■■■■      60% | ETA:  8s
#> ⠼ [2026-05-23 14:30:53] Running for 4 [4/5] ■■■■■■■■    80% | ETA:  4s
#> ✔ [2026-05-23 14:30:53] Completed 5 tasks in 18.1s
#> 
#> ℹ [2026-05-23 14:30:53] Building results
#> ✔ [2026-05-23 14:31:11] Enrichment analysis done
EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP",
  group.by = "CellType",
  plot_type = "comparison"
)


pancreas_sub <- RunEnrichment(
  pancreas_sub,
  group.by = "CellType",
  DE_threshold = "p_val_adj < 0.05",
  db = c("MSigDB", "MSigDB_MH"),
  species = "Mus_musculus"
)
#> ℹ [2026-05-23 14:31:12] Start Enrichment analysis
#> ℹ [2026-05-23 14:31:12] Species: "Mus_musculus"
#> ℹ [2026-05-23 14:31:12] Preparing MSigDB database
#> ℹ [2026-05-23 14:31:25] Permform enrichment...
#> ℹ [2026-05-23 14:31:25] Using 1 core
#> ⠙ [2026-05-23 14:31:25] Running for 1 [1/10] ■           10% | ETA:  2m
#> ⠹ [2026-05-23 14:31:25] Running for 2 [2/10] ■■          20% | ETA:  1m
#> ⠸ [2026-05-23 14:31:25] Running for 3 [3/10] ■■■         30% | ETA:  1m
#> ⠼ [2026-05-23 14:31:25] Running for 4 [4/10] ■■■■        40% | ETA: 43s
#> ⠴ [2026-05-23 14:31:25] Running for 5 [5/10] ■■■■■       50% | ETA: 33s
#> ✔ [2026-05-23 14:31:25] Completed 10 tasks in 33.2s
#> 
#> ℹ [2026-05-23 14:31:25] Building results
#> ✔ [2026-05-23 14:31:58] Enrichment analysis done
EnrichmentPlot(
  pancreas_sub,
  db = "MSigDB",
  group.by = "CellType",
  plot_type = "comparison"
)

EnrichmentPlot(
  pancreas_sub,
  db = "MSigDB_MH",
  group.by = "CellType",
  plot_type = "comparison"
)


# Remove redundant GO terms
pancreas_sub <- RunEnrichment(
  pancreas_sub,
  group.by = "CellType",
  db = "GO_BP",
  GO_simplify = TRUE,
  species = "Mus_musculus"
)
#> ℹ [2026-05-23 14:31:59] Start Enrichment analysis
#> ! [2026-05-23 14:31:59] `GO_simplify = TRUE` requires clusterProfiler result objects; using `backend = 'r'` for this run.
#> ℹ [2026-05-23 14:31:59] Species: "Mus_musculus"
#> ℹ [2026-05-23 14:31:59] Loading cached: GO_BP version: 3.23.0 nterm:14957 created: 2026-05-23 13:33:03
#> ℹ [2026-05-23 14:32:00] Permform enrichment...
#> ℹ [2026-05-23 14:32:00] Using 1 core
#> ⠙ [2026-05-23 14:32:00] Running for 1 [1/5] ■■          20% | ETA: 14m
#> ⠹ [2026-05-23 14:32:00] Running for 2 [2/5] ■■■■        40% | ETA:  6m
#> ⠸ [2026-05-23 14:32:00] Running for 3 [3/5] ■■■■■■      60% | ETA:  3m
#> ⠼ [2026-05-23 14:32:00] Running for 4 [4/5] ■■■■■■■■    80% | ETA:  1m
#> ✔ [2026-05-23 14:32:00] Completed 5 tasks in 5m 53.1s
#> 
#> ℹ [2026-05-23 14:32:00] Building results
#> ✔ [2026-05-23 14:37:54] Enrichment analysis done
EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP_sim",
  group.by = "CellType",
  plot_type = "comparison"
)


# Or use "geneID" and "geneID_groups" as input to run enrichment
de_df <- dplyr::filter(
  pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox,
  p_val_adj < 0.05
)
enrich_out <- RunEnrichment(
  geneID = de_df[["gene"]],
  geneID_groups = de_df[["group1"]],
  db = "GO_BP",
  species = "Mus_musculus"
)
#> ℹ [2026-05-23 14:37:54] Start Enrichment analysis
#> ℹ [2026-05-23 14:37:54] Species: "Mus_musculus"
#> ℹ [2026-05-23 14:37:54] Loading cached: GO_BP version: 3.23.0 nterm:14957 created: 2026-05-23 13:33:03
#> ℹ [2026-05-23 14:37:56] Permform enrichment...
#> ℹ [2026-05-23 14:37:56] Using 1 core
#> ⠙ [2026-05-23 14:37:56] Running for 1 [1/5] ■■          20% | ETA: 23s
#> ⠹ [2026-05-23 14:37:56] Running for 2 [2/5] ■■■■        40% | ETA: 13s
#> ⠸ [2026-05-23 14:37:56] Running for 3 [3/5] ■■■■■■      60% | ETA:  8s
#> ⠼ [2026-05-23 14:37:56] Running for 4 [4/5] ■■■■■■■■    80% | ETA:  4s
#> ✔ [2026-05-23 14:37:56] Completed 5 tasks in 18.2s
#> 
#> ℹ [2026-05-23 14:37:56] Building results
#> ✔ [2026-05-23 14:38:14] Enrichment analysis done
EnrichmentPlot(
  res = enrich_out,
  db = "GO_BP",
  plot_type = "comparison"
)


# Use a combined database
pancreas_sub <- RunEnrichment(
  pancreas_sub,
  group.by = "CellType",
  db = c(
    "KEGG", "WikiPathway", "Reactome", "PFAM", "MP"
  ),
  db_combine = TRUE,
  species = "Mus_musculus"
)
#> ℹ [2026-05-23 14:38:14] Start Enrichment analysis
#> ℹ [2026-05-23 14:38:14] Species: "Mus_musculus"
#> ℹ [2026-05-23 14:38:33] Preparing KEGG database
#> ℹ [2026-05-23 14:38:48] Preparing WikiPathway database
#> ℹ [2026-05-23 14:38:51] Preparing Reactome database
#> ℹ [2026-05-23 14:38:56] Preparing MP database
#> ℹ [2026-05-23 14:39:06] Preparing PFAM database
#> 
#> ℹ [2026-05-23 14:39:07] Convert ID types for the KEGG database
#> ℹ [2026-05-23 14:39:07] Converted ID types using local annotation package org.Mm.eg.db
#> ℹ [2026-05-23 14:39:07] Convert ID types for the WikiPathway database
#> ℹ [2026-05-23 14:39:08] Converted ID types using local annotation package org.Mm.eg.db
#> ℹ [2026-05-23 14:39:08] Convert ID types for the Reactome database
#> ℹ [2026-05-23 14:39:08] Converted ID types using local annotation package org.Mm.eg.db
#> ℹ [2026-05-23 14:39:08] Convert ID types for the PFAM database
#> ℹ [2026-05-23 14:39:08] Converted ID types using local annotation package org.Mm.eg.db
#> ℹ [2026-05-23 14:39:08] Create 'Combined' database ...
#> ℹ [2026-05-23 14:39:09] Permform enrichment...
#> ℹ [2026-05-23 14:39:09] Using 1 core
#> ⠙ [2026-05-23 14:39:09] Running for 1 [1/5] ■■          20% | ETA: 13s
#> ⠹ [2026-05-23 14:39:09] Running for 3 [3/5] ■■■■■■      60% | ETA:  5s
#> ⠸ [2026-05-23 14:39:09] Running for 4 [4/5] ■■■■■■■■    80% | ETA:  2s
#> ✔ [2026-05-23 14:39:09] Completed 5 tasks in 9.9s
#> 
#> ℹ [2026-05-23 14:39:09] Building results
#> ✔ [2026-05-23 14:39:19] Enrichment analysis done
EnrichmentPlot(
  pancreas_sub,
  db = "Combined",
  group.by = "CellType",
  plot_type = "comparison"
)
```
