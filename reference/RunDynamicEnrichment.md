# RunDynamicEnrichment

This function calculates gene-set scores from the specified database
(`db`) for each lineage using the specified scoring method
(`score_method`). It then treats these scores as expression values and
uses them as input to the RunDynamicFeatures function to identify
dynamically enriched terms along the lineage.

## Usage

``` r
RunDynamicEnrichment(
  srt,
  lineages,
  score_method = "AUCell",
  layer = "data",
  assay = NULL,
  min_expcells = 20,
  r.sq = 0.2,
  dev.expl = 0.2,
  padjust = 0.05,
  IDtype = "symbol",
  species = "Homo_sapiens",
  db = "GO_BP",
  db_update = FALSE,
  db_version = "latest",
  convert_species = TRUE,
  Ensembl_version = NULL,
  mirror = NULL,
  TERM2GENE = NULL,
  TERM2NAME = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  cores = 1,
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- srt:

  A Seurat object containing the results of differential expression
  analysis (RunDEtest). If specified, the genes and groups will be
  extracted from the Seurat object automatically. If not specified, the
  `geneID` and `geneID_groups` arguments must be provided.

- lineages:

  A character vector specifying the lineages to plot.

- score_method:

  The method to use for scoring. Can be `"Seurat"`, `"AUCell"`, or
  `"UCell"`. Default is `"Seurat"`.

- layer:

  Which layer to use. Default is `"counts"`.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used.

- min_expcells:

  The minimum number of expected cells. Default is `20`.

- r.sq:

  The R-squared threshold. Default is `0.2`.

- dev.expl:

  The deviance explained threshold. Default is `0.2`.

- padjust:

  The p-value adjustment threshold. Default is `0.05`.

- IDtype:

  A character vector specifying the type of gene IDs in the `srt` object
  or `geneID` argument. This argument is used to convert the gene IDs to
  a different type if `IDtype` is different from `result_IDtype`.

- species:

  A character vector specifying the species for which the gene
  annotation databases should be prepared. Can be `"Homo_sapiens"` or
  `"Mus_musculus"`.

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

- cores:

  The number of cores to use for parallelization with
  [foreach::foreach](https://rdrr.io/pkg/foreach/man/foreach.html).
  Default is `1`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Random seed for reproducibility. Default is `11`.

## See also

[RunDynamicFeatures](https://mengxu98.github.io/scop/reference/RunDynamicFeatures.md),
[DynamicHeatmap](https://mengxu98.github.io/scop/reference/DynamicHeatmap.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-01-30 17:09:21] Start standard scop workflow...
#> ℹ [2026-01-30 17:09:22] Checking a list of <Seurat>...
#> ! [2026-01-30 17:09:22] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-01-30 17:09:22] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-30 17:09:24] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-30 17:09:24] Use the separate HVF from srt_list
#> ℹ [2026-01-30 17:09:24] Number of available HVF: 2000
#> ℹ [2026-01-30 17:09:25] Finished check
#> ℹ [2026-01-30 17:09:25] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-30 17:09:25] Perform pca linear dimension reduction
#> ℹ [2026-01-30 17:09:26] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-01-30 17:09:26] Reorder clusters...
#> ℹ [2026-01-30 17:09:26] Perform umap nonlinear dimension reduction
#> ℹ [2026-01-30 17:09:26] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-01-30 17:09:31] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-01-30 17:09:36] Run scop standard workflow completed
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP"
)

pancreas_sub <- RunDynamicFeatures(
  pancreas_sub,
  lineages = "Lineage1",
  n_candidates = 200
)
#> ℹ [2026-01-30 17:09:37] Start find dynamic features
#> ℹ [2026-01-30 17:09:38] Data type is raw counts
#> ℹ [2026-01-30 17:09:38] Number of candidate features (union): 200
#> ℹ [2026-01-30 17:09:39] Data type is raw counts
#> ℹ [2026-01-30 17:09:39] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-01-30 17:09:39] Using 1 core
#> ⠙ [2026-01-30 17:09:39] Running for Gcg [1/200] ■                              …
#> ⠹ [2026-01-30 17:09:39] Running for Neurog3 [20/200] ■■■■                      …
#> ⠸ [2026-01-30 17:09:39] Running for Cdca8 [72/200] ■■■■■■■■■■■■                …
#> ⠼ [2026-01-30 17:09:39] Running for Clps [124/200] ■■■■■■■■■■■■■■■■■■■■        …
#> ⠴ [2026-01-30 17:09:39] Running for Neurod2 [173/200] ■■■■■■■■■■■■■■■■■■■■■■■■■…
#> ✔ [2026-01-30 17:09:39] Completed 200 tasks in 11.9s
#> 
#> ℹ [2026-01-30 17:09:39] Building results
#> ✔ [2026-01-30 17:09:51] Find dynamic features done
ht1 <- DynamicHeatmap(
  pancreas_sub,
  lineages = "Lineage1",
  cell_annotation = "CellType",
  n_split = 3
)
#> ℹ [2026-01-30 17:09:51] [1] 146 features from Lineage1 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Ins1,Ins2,Nnat,Iapp,Lrpprc,Chgb,Slc38a5,2810417H13Rik,Rbp4...
#> ℹ [2026-01-30 17:09:52] 
#> ℹ                       The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ                       The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ                       If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.


pancreas_sub <- RunDynamicEnrichment(
  pancreas_sub,
  lineages = "Lineage1",
  score_method = "AUCell",
  db = "GO_BP",
  species = "Mus_musculus"
)
#> ℹ [2026-01-30 17:09:53] Species: "Mus_musculus"
#> ℹ [2026-01-30 17:09:53] Loading cached: GO_BP version: 3.22.0 nterm:15169 created: 2026-01-30 16:42:33
#> ℹ [2026-01-30 17:09:56] Start cell scoring
#> ℹ [2026-01-30 17:09:57] Data type is log-normalized
#> ℹ [2026-01-30 17:09:58] Number of feature lists to be scored: 2761
#> ✔ [2026-01-30 17:12:53] Cell scoring completed
#> ℹ [2026-01-30 17:12:53] Start find dynamic features
#> ℹ [2026-01-30 17:12:54] Data type is log-normalized
#> ℹ [2026-01-30 17:12:54] Number of candidate features (union): 2761
#> ℹ [2026-01-30 17:12:54] Data type is log-normalized
#> ℹ [2026-01-30 17:12:54] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-01-30 17:12:54] Using 1 core
#> ⠙ [2026-01-30 17:12:54] Running for GO-BP-2..deoxyribonucleotide.biosynthetic.p…
#> ⠹ [2026-01-30 17:12:54] Running for GO-BP-L.amino.acid.transport [50/2761] ■■  …
#> ⠸ [2026-01-30 17:12:54] Running for GO-BP-cardiac.septum.morphogenesis [263/276…
#> ⠼ [2026-01-30 17:12:54] Running for GO-BP-detection.of.abiotic.stimulus [475/27…
#> ⠴ [2026-01-30 17:12:54] Running for GO-BP-hair.cycle.process [688/2761] ■■■■■■■…
#> ⠦ [2026-01-30 17:12:54] Running for GO-BP-mechanoreceptor.differentiation [900/…
#> ⠧ [2026-01-30 17:12:54] Running for GO-BP-negative.regulation.of.cell.division …
#> ⠇ [2026-01-30 17:12:54] Running for GO-BP-neurotransmitter.receptor.internaliza…
#> ⠏ [2026-01-30 17:12:54] Running for GO-BP-positive.regulation.of.cell.division …
#> ⠋ [2026-01-30 17:12:54] Running for GO-BP-positive.regulation.of.tumor.necrosis…
#> ⠙ [2026-01-30 17:12:54] Running for GO-BP-regulation.of.acute.inflammatory.resp…
#> ⠹ [2026-01-30 17:12:54] Running for GO-BP-regulation.of.ketone.metabolic.proces…
#> ⠸ [2026-01-30 17:12:54] Running for GO-BP-regulation.of.smooth.muscle.cell.diff…
#> ⠼ [2026-01-30 17:12:54] Running for GO-BP-retinol.metabolic.process [2554/2761]…
#> ⠴ [2026-01-30 17:12:54] Running for GO-BP-vitamin.transmembrane.transport [2752…
#> ✔ [2026-01-30 17:12:54] Completed 2761 tasks in 39.9s
#> 
#> ℹ [2026-01-30 17:12:54] Building results
#> ✔ [2026-01-30 17:13:34] Find dynamic features done
#> ✔ [2026-01-30 17:13:34] Dynamic enrichment analysis completed
ht2 <- DynamicHeatmap(
  pancreas_sub,
  assay = "GO_BP",
  lineages = "Lineage1_GO_BP",
  cell_annotation = "CellType",
  n_split = 3,
  split_method = "kmeans-peaktime"
)
#> ℹ [2026-01-30 17:13:34] [1] 1897 features from Lineage1_GO_BP passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       GO-BP-2..deoxyribonucleotide.biosynthetic.process,GO-BP-2..deoxyribonucleotide.metabolic.process,GO-BP-ADP.catabolic.process,GO-BP-ADP.metabolic.process,GO-BP-ATP.biosynthetic.process,GO-BP-ATP.metabolic.process,GO-BP-ATP.synthesis.coupled.electron.transport,GO-BP-B.cell.activation,GO-BP-B.cell.apoptotic.process,GO-BP-B.cell.proliferation...
#> ! [2026-01-30 17:13:35] The values in the 'counts' layer are non-integer. Set the library size to 1.
#> ℹ [2026-01-30 17:13:36] 
#> ℹ                       The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ                       The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ                       If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
```
