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

  A Seurat object containing the results of differential expression
  analysis (RunDEtest). If specified, the genes and groups will be
  extracted from the Seurat object automatically. If not specified, the
  `geneID` and `geneID_groups` arguments must be provided.

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
#> ℹ [2026-05-11 15:57:57] Start standard processing workflow...
#> ℹ [2026-05-11 15:57:58] Checking a list of <Seurat>...
#> ! [2026-05-11 15:57:58] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-05-11 15:57:58] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-05-11 15:57:59] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-11 15:58:00] Use the separate HVF from `srt_list`
#> ℹ [2026-05-11 15:58:00] Number of available HVF: 2000
#> ℹ [2026-05-11 15:58:00] Finished check
#> ℹ [2026-05-11 15:58:00] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-11 15:58:00] Perform pca linear dimension reduction
#> ℹ [2026-05-11 15:58:01] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-05-11 15:58:01] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-11 15:58:01] Reorder clusters...
#> ℹ [2026-05-11 15:58:02] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-11 15:58:02] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-11 15:58:02] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-05-11 15:58:06] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-05-11 15:58:11] Standard processing workflow completed
pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "CellType"
)
#> ℹ [2026-05-11 15:58:11] Data type is log-normalized
#> ℹ [2026-05-11 15:58:11] Start differential expression test
#> ℹ [2026-05-11 15:58:11] Find all markers(wilcox) among [1] 5 groups...
#> ℹ [2026-05-11 15:58:11] Using 1 core
#> ⠙ [2026-05-11 15:58:11] Running for Ductal [1/5] ■■          20% | ETA:  1s
#> ✔ [2026-05-11 15:58:11] Completed 5 tasks in 896ms
#> 
#> ℹ [2026-05-11 15:58:11] Building results
#> ✔ [2026-05-11 15:58:12] Differential expression test completed
pancreas_sub <- RunEnrichment(
  pancreas_sub,
  group.by = "CellType",
  DE_threshold = "p_val_adj < 0.05",
  db = "GO_BP",
  species = "Mus_musculus"
)
#> ℹ [2026-05-11 15:58:12] Start Enrichment analysis
#> ℹ [2026-05-11 15:58:12] Species: "Mus_musculus"
#> ℹ [2026-05-11 15:58:12] Loading cached: GO_BP version: 3.23.0 nterm:14957 created: 2026-05-11 14:55:43
#> ℹ [2026-05-11 15:58:14] Permform enrichment...
#> ℹ [2026-05-11 15:58:14] Using 1 core
#> ⠙ [2026-05-11 15:58:14] Running for 1 [1/5] ■■          20% | ETA: 19s
#> ⠹ [2026-05-11 15:58:14] Running for 2 [2/5] ■■■■        40% | ETA: 11s
#> ⠸ [2026-05-11 15:58:14] Running for 3 [3/5] ■■■■■■      60% | ETA:  7s
#> ⠼ [2026-05-11 15:58:14] Running for 4 [4/5] ■■■■■■■■    80% | ETA:  3s
#> ✔ [2026-05-11 15:58:14] Completed 5 tasks in 16.8s
#> 
#> ℹ [2026-05-11 15:58:14] Building results
#> ✔ [2026-05-11 15:58:30] Enrichment analysis done
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
#> ℹ [2026-05-11 15:58:31] Start Enrichment analysis
#> ℹ [2026-05-11 15:58:31] Species: "Mus_musculus"
#> ℹ [2026-05-11 15:58:31] Preparing MSigDB database
#> ℹ [2026-05-11 15:58:43] Convert ID types for the MSigDB database
#> Error in switch(tolower(idtype), symbol = "SYMBOL", ensembl_id = "ENSEMBL",     entrez_id = org_key, tair_locus = "TAIR", sgd_gene = "SGD",     NA_character_): EXPR must be a length 1 vector
EnrichmentPlot(
  pancreas_sub,
  db = "MSigDB",
  group.by = "CellType",
  plot_type = "comparison"
)
#> Error in EnrichmentPlot(pancreas_sub, db = "MSigDB", group.by = "CellType",     plot_type = "comparison"): MSigDB is not in the enrichment result
EnrichmentPlot(
  pancreas_sub,
  db = "MSigDB_MH",
  group.by = "CellType",
  plot_type = "comparison"
)
#> Error in EnrichmentPlot(pancreas_sub, db = "MSigDB_MH", group.by = "CellType",     plot_type = "comparison"): MSigDB_MH is not in the enrichment result

# Remove redundant GO terms
pancreas_sub <- RunEnrichment(
  pancreas_sub,
  group.by = "CellType",
  db = "GO_BP",
  GO_simplify = TRUE,
  species = "Mus_musculus"
)
#> ℹ [2026-05-11 15:58:43] Start Enrichment analysis
#> ! [2026-05-11 15:58:43] `GO_simplify = TRUE` requires clusterProfiler result objects; using `backend = 'r'` for this run.
#> ℹ [2026-05-11 15:58:43] Species: "Mus_musculus"
#> ℹ [2026-05-11 15:58:43] Loading cached: GO_BP version: 3.23.0 nterm:14957 created: 2026-05-11 14:55:43
#> ℹ [2026-05-11 15:58:45] Permform enrichment...
#> ℹ [2026-05-11 15:58:45] Using 1 core
#> ⠙ [2026-05-11 15:58:45] Running for 1 [1/5] ■■          20% | ETA: 13m
#> ⠹ [2026-05-11 15:58:45] Running for 2 [2/5] ■■■■        40% | ETA:  6m
#> ⠸ [2026-05-11 15:58:45] Running for 3 [3/5] ■■■■■■      60% | ETA:  3m
#> ⠼ [2026-05-11 15:58:45] Running for 4 [4/5] ■■■■■■■■    80% | ETA:  1m
#> ✔ [2026-05-11 15:58:45] Completed 5 tasks in 5m 54.8s
#> 
#> ℹ [2026-05-11 15:58:45] Building results
#> ✔ [2026-05-11 16:04:40] Enrichment analysis done
EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP_sim",
  group.by = "CellType",
  plot_type = "comparison"
)
#> Error in `[.data.frame`(enrichment, "Groups"): undefined columns selected

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
#> ℹ [2026-05-11 16:04:40] Start Enrichment analysis
#> ℹ [2026-05-11 16:04:40] Species: "Mus_musculus"
#> ℹ [2026-05-11 16:04:40] Loading cached: GO_BP version: 3.23.0 nterm:14957 created: 2026-05-11 14:55:43
#> ℹ [2026-05-11 16:04:41] Permform enrichment...
#> ℹ [2026-05-11 16:04:41] Using 1 core
#> ⠙ [2026-05-11 16:04:41] Running for 1 [1/5] ■■          20% | ETA: 20s
#> ⠹ [2026-05-11 16:04:41] Running for 2 [2/5] ■■■■        40% | ETA: 12s
#> ⠸ [2026-05-11 16:04:41] Running for 3 [3/5] ■■■■■■      60% | ETA:  8s
#> ⠼ [2026-05-11 16:04:41] Running for 4 [4/5] ■■■■■■■■    80% | ETA:  4s
#> ✔ [2026-05-11 16:04:41] Completed 5 tasks in 19.6s
#> 
#> ℹ [2026-05-11 16:04:41] Building results
#> ✔ [2026-05-11 16:05:01] Enrichment analysis done
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
#> ℹ [2026-05-11 16:05:01] Start Enrichment analysis
#> ℹ [2026-05-11 16:05:01] Species: "Mus_musculus"
#> ✔ [2026-05-11 16:05:01] org.Mm.eg.db installed successfully
#> ℹ [2026-05-11 16:05:20] Preparing KEGG database
#> ℹ [2026-05-11 16:05:32] Preparing WikiPathway database
#> ℹ [2026-05-11 16:05:35] Preparing Reactome database
#> ℹ [2026-05-11 16:05:40] Preparing MP database
#> ℹ [2026-05-11 16:05:52] Preparing PFAM database
#> 
#> ℹ [2026-05-11 16:05:53] Convert ID types for the KEGG database
#> ℹ [2026-05-11 16:05:53] Converted ID types using local annotation package org.Mm.eg.db
#> ℹ [2026-05-11 16:05:53] Convert ID types for the WikiPathway database
#> ℹ [2026-05-11 16:05:53] Converted ID types using local annotation package org.Mm.eg.db
#> ℹ [2026-05-11 16:05:53] Convert ID types for the Reactome database
#> ℹ [2026-05-11 16:05:54] Converted ID types using local annotation package org.Mm.eg.db
#> ℹ [2026-05-11 16:05:54] Convert ID types for the PFAM database
#> ℹ [2026-05-11 16:05:54] Converted ID types using local annotation package org.Mm.eg.db
#> ℹ [2026-05-11 16:05:54] Create 'Combined' database ...
#> ℹ [2026-05-11 16:05:54] Permform enrichment...
#> ℹ [2026-05-11 16:05:54] Using 1 core
#> ⠙ [2026-05-11 16:05:54] Running for 1 [1/5] ■■          20% | ETA:  9s
#> ⠹ [2026-05-11 16:05:54] Running for 4 [4/5] ■■■■■■■■    80% | ETA:  1s
#> ✔ [2026-05-11 16:05:54] Completed 5 tasks in 5.1s
#> 
#> ℹ [2026-05-11 16:05:54] Building results
#> ✔ [2026-05-11 16:06:00] Enrichment analysis done
EnrichmentPlot(
  pancreas_sub,
  db = "Combined",
  group.by = "CellType",
  plot_type = "comparison"
)
```
