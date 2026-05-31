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
  features = NULL,
  TERM2GENE = NULL,
  TERM2NAME = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  backend = c("cpp", "r"),
  cpp_strategy = c("sparse", "topk", "full"),
  cores = 1,
  verbose = TRUE,
  seed = 11,
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

- lineages:

  A character vector specifying the lineage names for which dynamic
  features should be calculated.

- score_method:

  The method to use for scoring. Can be `"Seurat"`, `"AUCell"`,
  `"UCell"`, `"GSVA"`, `"ssGSEA"`, `"zscore"`, `"PLAGE"`, or `"VISION"`.
  Multiple methods can be supplied at once; each method will be written
  to a method-suffixed assay before dynamic-feature fitting. Default is
  `"AUCell"`.

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

- ...:

  Passed to other functions.

## See also

[RunDynamicFeatures](https://mengxu98.github.io/scop/reference/RunDynamicFeatures.md),
[DynamicHeatmap](https://mengxu98.github.io/scop/reference/DynamicHeatmap.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-05-31 06:48:28] Start standard processing workflow...
#> ℹ [2026-05-31 06:48:28] Checking a list of <Seurat>...
#> ! [2026-05-31 06:48:28] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-05-31 06:48:28] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-05-31 06:48:30] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-31 06:48:31] Use the separate HVF from `srt_list`
#> ℹ [2026-05-31 06:48:31] Number of available HVF: 2000
#> ℹ [2026-05-31 06:48:31] Finished check
#> ℹ [2026-05-31 06:48:31] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-31 06:48:31] Perform pca linear dimension reduction
#> ℹ [2026-05-31 06:48:32] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-05-31 06:48:32] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-31 06:48:32] Reorder clusters...
#> ℹ [2026-05-31 06:48:33] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-31 06:48:33] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-31 06:48:33] Perform umap nonlinear dimension reduction using Standardpca (1:23)
#> ℹ [2026-05-31 06:48:38] Perform umap nonlinear dimension reduction using Standardpca (1:23)
#> ✔ [2026-05-31 06:48:42] Standard processing workflow completed
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP"
)
#> Warning: Removed 9 rows containing missing values or values outside the scale range
#> (`geom_path()`).
#> Warning: Removed 9 rows containing missing values or values outside the scale range
#> (`geom_path()`).

pancreas_sub <- RunDynamicFeatures(
  pancreas_sub,
  lineages = "Lineage1",
  fit_method = "pretsa",
  n_candidates = 200
)
#> ℹ [2026-05-31 06:48:43] Start find dynamic features
#> ℹ [2026-05-31 06:48:44] Data type is raw counts
#> ℹ [2026-05-31 06:48:45] Number of candidate features (union): 200
#> ℹ [2026-05-31 06:48:45] Data type is raw counts
#> ℹ [2026-05-31 06:48:45] Calculating dynamic features for "Lineage1"...
#> ✔ [2026-05-31 06:48:45] Find dynamic features done
ht1 <- DynamicHeatmap(
  pancreas_sub,
  lineages = "Lineage1",
  cell_annotation = "CellType",
  n_split = 3
)
#> ℹ [2026-05-31 06:48:45] [1] 132 features from Lineage1 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Ins1,Ins2,Nnat,Iapp,Lrpprc,Npy,Chgb,Slc38a5,2810417H13Rik,Rbp4...
#> Error in heatmap_enrichment(geneID = feature_metadata[["features"]], geneID_groups = feature_metadata[["feature_split"]],     feature_split_palette = feature_split_palette, feature_split_palcolor = feature_split_palcolor,     ha_right = ha_right, flip = flip, anno_terms = anno_terms,     anno_keys = anno_keys, anno_features = anno_features, terms_width = terms_width,     terms_fontsize = terms_fontsize, keys_width = keys_width,     keys_fontsize = keys_fontsize, features_width = features_width,     features_fontsize = features_fontsize, IDtype = IDtype, species = species,     db_update = db_update, db_version = db_version, db_combine = db_combine,     convert_species = convert_species, Ensembl_version = Ensembl_version,     mirror = mirror, db = db, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,     minGSSize = minGSSize, maxGSSize = maxGSSize, GO_simplify = GO_simplify,     GO_simplify_cutoff = GO_simplify_cutoff, simplify_method = simplify_method,     simplify_similarityCutoff = simplify_similarityCutoff, pvalueCutoff = pvalueCutoff,     padjustCutoff = padjustCutoff, topTerm = topTerm, show_termid = show_termid,     topWord = topWord, words_excluded = words_excluded, cores = cores,     ...): '...' used in an incorrect context

pancreas_sub <- RunDynamicEnrichment(
  pancreas_sub,
  lineages = "Lineage1",
  score_method = "AUCell",
  db = "GO_BP",
  species = "Mus_musculus"
)
#> ℹ [2026-05-31 06:48:46] Species: "Mus_musculus"
#> ℹ [2026-05-31 06:48:46] Loading cached: GO_BP version: 3.23.0 nterm:14957 created: 2026-05-31 05:58:25
#> ℹ [2026-05-31 06:48:49] Start cell scoring
#> ℹ [2026-05-31 06:48:49] Data type is log-normalized
#> ℹ [2026-05-31 06:48:51] Number of feature lists to be scored: 2719
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#> ✔ [2026-05-31 06:48:57] Cell scoring completed
#> ℹ [2026-05-31 06:48:57] Start find dynamic features
#> ℹ [2026-05-31 06:48:58] Data type is log-normalized
#> ℹ [2026-05-31 06:48:58] Number of candidate features (union): 2719
#> ℹ [2026-05-31 06:48:59] Data type is log-normalized
#> ℹ [2026-05-31 06:48:59] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-05-31 06:48:59] Using 1 core
#> ⠙ [2026-05-31 06:48:59] Running for GO-BP-2..deoxyribonucleotide.biosynthetic.p…
#> ⠹ [2026-05-31 06:48:59] Running for GO-BP-B.cell.differentiation [9/2719]      …
#> ⠸ [2026-05-31 06:48:59] Running for GO-BP-branching.morphogenesis.of.an.epithel…
#> ⠼ [2026-05-31 06:48:59] Running for GO-BP-chromosome.localization [402/2719] ■ …
#> ⠴ [2026-05-31 06:48:59] Running for GO-BP-fatty.acid.catabolic.process [600/271…
#> ⠦ [2026-05-31 06:48:59] Running for GO-BP-left.right.pattern.formation [798/271…
#> ⠧ [2026-05-31 06:48:59] Running for GO-BP-monoatomic.ion.homeostasis [992/2719]…
#> ⠇ [2026-05-31 06:48:59] Running for GO-BP-negative.regulation.of.lymphocyte.apo…
#> ⠏ [2026-05-31 06:48:59] Running for GO-BP-olfactory.bulb.interneuron.differenti…
#> ⠋ [2026-05-31 06:48:59] Running for GO-BP-positive.regulation.of.glutamate.secr…
#> ⠙ [2026-05-31 06:48:59] Running for GO-BP-protein.complex.oligomerization [1774…
#> ⠹ [2026-05-31 06:48:59] Running for GO-BP-regulation.of.cell.cycle.checkpoint […
#> ⠸ [2026-05-31 06:48:59] Running for GO-BP-regulation.of.mitotic.cytokinesis [21…
#> ⠼ [2026-05-31 06:48:59] Running for GO-BP-regulation.of.vasoconstriction [2371/…
#> ⠴ [2026-05-31 06:48:59] Running for GO-BP-sperm.motility [2568/2719] ■■■■■■■■■ …
#> ✔ [2026-05-31 06:48:59] Completed 2719 tasks in 41.4s
#> 
#> ℹ [2026-05-31 06:48:59] Building results
#> ✔ [2026-05-31 06:49:41] Find dynamic features done
#> ✔ [2026-05-31 06:49:41] Dynamic enrichment analysis completed
ht2 <- DynamicHeatmap(
  pancreas_sub,
  assay = "GO_BP",
  lineages = "Lineage1_GO_BP",
  cell_annotation = "CellType",
  n_split = 3,
  split_method = "kmeans-peaktime"
)
#> ℹ [2026-05-31 06:49:41] [1] 1881 features from Lineage1_GO_BP passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       GO-BP-2..deoxyribonucleotide.biosynthetic.process,GO-BP-2..deoxyribonucleotide.metabolic.process,GO-BP-ADP.catabolic.process,GO-BP-ADP.metabolic.process,GO-BP-ATP.metabolic.process,GO-BP-ATP.synthesis.coupled.electron.transport,GO-BP-B.cell.activation,GO-BP-B.cell.proliferation,GO-BP-CENP.A.containing.chromatin.assembly,GO-BP-D.glucose.import.across.plasma.membrane...
#> ! [2026-05-31 06:49:41] The values in the 'counts' layer are non-integer. Set the library size to 1.
#> Error in heatmap_enrichment(geneID = feature_metadata[["features"]], geneID_groups = feature_metadata[["feature_split"]],     feature_split_palette = feature_split_palette, feature_split_palcolor = feature_split_palcolor,     ha_right = ha_right, flip = flip, anno_terms = anno_terms,     anno_keys = anno_keys, anno_features = anno_features, terms_width = terms_width,     terms_fontsize = terms_fontsize, keys_width = keys_width,     keys_fontsize = keys_fontsize, features_width = features_width,     features_fontsize = features_fontsize, IDtype = IDtype, species = species,     db_update = db_update, db_version = db_version, db_combine = db_combine,     convert_species = convert_species, Ensembl_version = Ensembl_version,     mirror = mirror, db = db, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,     minGSSize = minGSSize, maxGSSize = maxGSSize, GO_simplify = GO_simplify,     GO_simplify_cutoff = GO_simplify_cutoff, simplify_method = simplify_method,     simplify_similarityCutoff = simplify_similarityCutoff, pvalueCutoff = pvalueCutoff,     padjustCutoff = padjustCutoff, topTerm = topTerm, show_termid = show_termid,     topWord = topWord, words_excluded = words_excluded, cores = cores,     ...): '...' used in an incorrect context
```
