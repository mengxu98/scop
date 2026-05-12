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
  backend = c("cpp", "r"),
  cpp_strategy = c("sparse", "topk", "full"),
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

  A character vector specifying the lineage names for which dynamic
  features should be calculated.

- score_method:

  The method to use for scoring. Can be `"Seurat"`, `"AUCell"`, or
  `"UCell"`. Default is `"Seurat"`.

- layer:

  Which layer to use. Default is `"counts"`.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

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

- backend:

  Enrichment backend. `"cpp"` is the default and uses a fast native
  hypergeometric ORA implementation and returns the enrichment table
  without `enrichResult` objects. `"r"` uses
  [`clusterProfiler::enricher()`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html)
  and returns `enrichResult` objects in `results`. `GO_simplify = TRUE`
  currently uses the R backend.

- cpp_strategy:

  C++ AUCell ranking strategy. `"sparse"` ranks non-zero genes and
  approximates zero ties, `"topk"` ranks only genes that can contribute
  to AUCell AUC, and `"full"` ranks all genes.

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
#> ℹ [2026-05-12 04:49:17] Start standard processing workflow...
#> ℹ [2026-05-12 04:49:17] Checking a list of <Seurat>...
#> ! [2026-05-12 04:49:17] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-05-12 04:49:17] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-05-12 04:49:19] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-12 04:49:20] Use the separate HVF from `srt_list`
#> ℹ [2026-05-12 04:49:20] Number of available HVF: 2000
#> ℹ [2026-05-12 04:49:20] Finished check
#> ℹ [2026-05-12 04:49:20] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-12 04:49:20] Perform pca linear dimension reduction
#> ℹ [2026-05-12 04:49:21] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-05-12 04:49:21] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-12 04:49:21] Reorder clusters...
#> ℹ [2026-05-12 04:49:21] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-12 04:49:22] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-12 04:49:22] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-05-12 04:49:26] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-05-12 04:49:31] Standard processing workflow completed
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP"
)
#> Warning: Removed 7 rows containing missing values or values outside the scale range
#> (`geom_path()`).
#> Warning: Removed 7 rows containing missing values or values outside the scale range
#> (`geom_path()`).

pancreas_sub <- RunDynamicFeatures(
  pancreas_sub,
  lineages = "Lineage1",
  fit_method = "pretsa",
  n_candidates = 200
)
#> ℹ [2026-05-12 04:49:32] Start find dynamic features
#> ℹ [2026-05-12 04:49:33] Data type is raw counts
#> ℹ [2026-05-12 04:49:33] Number of candidate features (union): 200
#> ℹ [2026-05-12 04:49:34] Data type is raw counts
#> ℹ [2026-05-12 04:49:34] Calculating dynamic features for "Lineage1"...
#> ✔ [2026-05-12 04:49:34] Find dynamic features done
ht1 <- DynamicHeatmap(
  pancreas_sub,
  lineages = "Lineage1",
  cell_annotation = "CellType",
  n_split = 3
)
#> ℹ [2026-05-12 04:49:34] [1] 130 features from Lineage1 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Ins1,Ins2,Nnat,Iapp,Lrpprc,Npy,Chgb,Slc38a5,2810417H13Rik,Rbp4...
#> ℹ [2026-05-12 04:49:35] 
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
#> ℹ [2026-05-12 04:49:36] Species: "Mus_musculus"
#> ℹ [2026-05-12 04:49:36] Loading cached: GO_BP version: 3.23.0 nterm:14957 created: 2026-05-12 03:50:03
#> ℹ [2026-05-12 04:49:39] Start cell scoring
#> ℹ [2026-05-12 04:49:40] Data type is log-normalized
#> ℹ [2026-05-12 04:49:41] Number of feature lists to be scored: 2717
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#> ✔ [2026-05-12 04:49:44] Cell scoring completed
#> ℹ [2026-05-12 04:49:44] Start find dynamic features
#> ℹ [2026-05-12 04:49:45] Data type is log-normalized
#> ℹ [2026-05-12 04:49:45] Number of candidate features (union): 2717
#> ℹ [2026-05-12 04:49:46] Data type is log-normalized
#> ℹ [2026-05-12 04:49:46] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-05-12 04:49:46] Using 1 core
#> ⠙ [2026-05-12 04:49:46] Running for GO-BP-2..deoxyribonucleotide.biosynthetic.p…
#> ⠹ [2026-05-12 04:49:46] Running for GO-BP-DNA.strand.elongation [25/2717]      …
#> ⠸ [2026-05-12 04:49:46] Running for GO-BP-branching.morphogenesis.of.an.epithel…
#> ⠼ [2026-05-12 04:49:46] Running for GO-BP-cholangiocyte.proliferation [395/2717…
#> ⠴ [2026-05-12 04:49:46] Running for GO-BP-establishment.of.mitotic.spindle.loca…
#> ⠦ [2026-05-12 04:49:46] Running for GO-BP-interferon.beta.production [753/2717]…
#> ⠧ [2026-05-12 04:49:46] Running for GO-BP-metanephric.nephron.tubule.developmen…
#> ⠇ [2026-05-12 04:49:46] Running for GO-BP-negative.regulation.of.cell.cycle.G1.…
#> ⠏ [2026-05-12 04:49:46] Running for GO-BP-nephron.epithelium.morphogenesis [128…
#> ⠋ [2026-05-12 04:49:46] Running for GO-BP-plasminogen.activation [1457/2717] ■■…
#> ⠙ [2026-05-12 04:49:46] Running for GO-BP-positive.regulation.of.nervous.system…
#> ⠹ [2026-05-12 04:49:46] Running for GO-BP-protein.tetramerization [1820/2717] ■…
#> ⠸ [2026-05-12 04:49:46] Running for GO-BP-regulation.of.cellular.response.to.tr…
#> ⠼ [2026-05-12 04:49:46] Running for GO-BP-regulation.of.mitotic.metaphase.anaph…
#> ⠴ [2026-05-12 04:49:46] Running for GO-BP-regulation.of.transferase.activity [2…
#> ⠦ [2026-05-12 04:49:46] Running for GO-BP-serotonin.transport [2527/2717] ■■■■■…
#> ⠧ [2026-05-12 04:49:46] Running for GO-BP-ventricular.cardiac.muscle.cell.diffe…
#> ✔ [2026-05-12 04:49:46] Completed 2717 tasks in 46s
#> 
#> ℹ [2026-05-12 04:49:46] Building results
#> ✔ [2026-05-12 04:50:32] Find dynamic features done
#> ✔ [2026-05-12 04:50:32] Dynamic enrichment analysis completed
ht2 <- DynamicHeatmap(
  pancreas_sub,
  assay = "GO_BP",
  lineages = "Lineage1_GO_BP",
  cell_annotation = "CellType",
  n_split = 3,
  split_method = "kmeans-peaktime"
)
#> ℹ [2026-05-12 04:50:32] [1] 1841 features from Lineage1_GO_BP passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       GO-BP-2..deoxyribonucleotide.biosynthetic.process,GO-BP-2..deoxyribonucleotide.metabolic.process,GO-BP-ADP.catabolic.process,GO-BP-ADP.metabolic.process,GO-BP-ATP.metabolic.process,GO-BP-ATP.synthesis.coupled.electron.transport,GO-BP-B.cell.activation,GO-BP-B.cell.proliferation,GO-BP-CENP.A.containing.chromatin.assembly,GO-BP-D.glucose.import.across.plasma.membrane...
#> ! [2026-05-12 04:50:32] The values in the 'counts' layer are non-integer. Set the library size to 1.
#> ℹ [2026-05-12 04:50:33] 
#> ℹ                       The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ                       The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ                       If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
```
