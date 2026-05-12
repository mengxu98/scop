# Cell scoring

This function performs cell scoring on a Seurat object. It calculates
scores for a given set of features and adds the scores as metadata to
the Seurat object.

## Usage

``` r
CellScoring(
  srt,
  features = NULL,
  layer = "data",
  assay = NULL,
  split.by = NULL,
  IDtype = "symbol",
  species = "Homo_sapiens",
  db = "GO_BP",
  termnames = NULL,
  db_update = FALSE,
  db_version = "latest",
  convert_species = TRUE,
  Ensembl_version = NULL,
  mirror = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  method = "Seurat",
  backend = c("cpp", "r"),
  cpp_strategy = c("sparse", "topk", "full"),
  classification = TRUE,
  name = "",
  new_assay = FALSE,
  seed = 11,
  cores = 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- features:

  A named list of feature lists for scoring. If `NULL`, `db` will be
  used to create features sets.

- layer:

  Which layer to use. Default is `data`.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- split.by:

  Name of a column in meta.data column to split plot by. Default is
  `NULL`.

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

- termnames:

  A vector of term names to be used from the database. Default is
  `NULL`, in which case all features from the database are used.

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

- minGSSize:

  The minimum size of a gene set to be considered in the enrichment
  analysis.

- maxGSSize:

  The maximum size of a gene set to be considered in the enrichment
  analysis.

- method:

  The method to use for scoring. Can be "Seurat", "AUCell", or "UCell".
  Default is `"Seurat"`.

- backend:

  Scoring backend. `"cpp"` is the default for supported methods. `"r"`
  uses the original package implementation. `"cpp"` currently supports
  `method = "Seurat"` and `method = "AUCell"`. `method = "UCell"` falls
  back to `"r"` when `backend` is not explicitly set.

- cpp_strategy:

  C++ AUCell ranking strategy. `"sparse"` ranks non-zero genes and
  approximates zero ties, `"topk"` ranks only genes that can contribute
  to AUCell AUC, and `"full"` ranks all genes.

- classification:

  Whether to perform classification based on the scores. Default is
  `TRUE`.

- name:

  The name of the assay to store the scores in. Only used if new_assay
  is TRUE. Default is `""`.

- new_assay:

  Whether to create a new assay for storing the scores. Default is
  `FALSE`.

- seed:

  Random seed for reproducibility. Default is `11`.

- cores:

  The number of cores to use for parallelization with
  [foreach::foreach](https://rdrr.io/pkg/foreach/man/foreach.html).
  Default is `1`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments to be passed to the scoring methods.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-05-12 03:29:16] Start standard processing workflow...
#> ℹ [2026-05-12 03:29:17] Checking a list of <Seurat>...
#> ! [2026-05-12 03:29:17] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-05-12 03:29:17] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-05-12 03:29:18] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-12 03:29:18] Use the separate HVF from `srt_list`
#> ℹ [2026-05-12 03:29:18] Number of available HVF: 2000
#> ℹ [2026-05-12 03:29:19] Finished check
#> ℹ [2026-05-12 03:29:19] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-12 03:29:19] Perform pca linear dimension reduction
#> ℹ [2026-05-12 03:29:19] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-05-12 03:29:20] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-12 03:29:20] Reorder clusters...
#> ℹ [2026-05-12 03:29:20] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-12 03:29:20] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-12 03:29:20] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-05-12 03:29:23] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-05-12 03:29:26] Standard processing workflow completed
features_all <- rownames(pancreas_sub)
pancreas_sub <- CellScoring(
  pancreas_sub,
  features = list(
    A = features_all[1:100],
    B = features_all[101:200]
  ),
  method = "AUCell",
  name = "test"
)
#> ℹ [2026-05-12 03:29:26] Start cell scoring
#> ℹ [2026-05-12 03:29:26] Data type is log-normalized
#> ℹ [2026-05-12 03:29:26] Number of feature lists to be scored: 2
#> ✔ [2026-05-12 03:29:27] Cell scoring completed
CellDimPlot(pancreas_sub, "test_classification")


FeatureDimPlot(pancreas_sub, "test_A")


data(panc8_sub)
  panc8_sub <- integration_scop(
    panc8_sub,
    batch = "tech",
    integration_method = "Harmony"
  )
#> ◌ [2026-05-12 03:29:27] Run integration workflow...
#> ℹ [2026-05-12 03:29:28] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-05-12 03:29:28] Checking a list of <Seurat>...
#> ! [2026-05-12 03:29:28] Data 1/5 of the `srt_list` is "unknown"
#> ℹ [2026-05-12 03:29:28] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/5 of `srt_list`...
#> ℹ [2026-05-12 03:29:29] Perform `Seurat::FindVariableFeatures()` on 1/5 of `srt_list`...
#> ! [2026-05-12 03:29:30] Data 2/5 of the `srt_list` is "unknown"
#> ℹ [2026-05-12 03:29:30] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 2/5 of `srt_list`...
#> ℹ [2026-05-12 03:29:30] Perform `Seurat::FindVariableFeatures()` on 2/5 of `srt_list`...
#> ! [2026-05-12 03:29:31] Data 3/5 of the `srt_list` is "unknown"
#> ℹ [2026-05-12 03:29:31] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 3/5 of `srt_list`...
#> ℹ [2026-05-12 03:29:32] Perform `Seurat::FindVariableFeatures()` on 3/5 of `srt_list`...
#> ! [2026-05-12 03:29:32] Data 4/5 of the `srt_list` is "unknown"
#> ℹ [2026-05-12 03:29:32] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 4/5 of `srt_list`...
#> ℹ [2026-05-12 03:29:33] Perform `Seurat::FindVariableFeatures()` on 4/5 of `srt_list`...
#> ! [2026-05-12 03:29:33] Data 5/5 of the `srt_list` is "unknown"
#> ℹ [2026-05-12 03:29:33] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 5/5 of `srt_list`...
#> ℹ [2026-05-12 03:29:34] Perform `Seurat::FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-05-12 03:29:35] Use the separate HVF from `srt_list`
#> ℹ [2026-05-12 03:29:35] Number of available HVF: 2000
#> ℹ [2026-05-12 03:29:35] Finished check
#> Warning: Layer ‘scale.data’ is empty
#> ℹ [2026-05-12 03:29:38] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-12 03:29:38] Perform linear dimension reduction("pca")
#> ℹ [2026-05-12 03:29:39] Perform Harmony integration
#> ℹ [2026-05-12 03:29:39] Using "Harmonypca" (1:20) as input
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 't': ‘Z_corr’ is not a valid field or method name for reference class “Rcpp_harmony”

  panc8_sub <- CellScoring(
    panc8_sub,
    layer = "data",
    assay = "RNA",
    db = "GO_BP",
    species = "Homo_sapiens",
    minGSSize = 10,
    maxGSSize = 100,
    method = "AUCell",
    name = "GO",
    new_assay = TRUE
  )
#> ℹ [2026-05-12 03:29:47] Start cell scoring
#> Warning: Layer ‘data’ is empty
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> ! [2026-05-12 03:29:47] Infinite values detected
#> ℹ [2026-05-12 03:29:47] Species: "Homo_sapiens"
#> 
#> ✔ [2026-05-12 03:35:06] org.Hs.eg.db installed successfully
#> ℹ [2026-05-12 03:35:34] Preparing database: GO_BP
#> ℹ [2026-05-12 03:36:03] Convert ID types for the GO_BP database
#> ℹ [2026-05-12 03:36:03] Converted ID types using local annotation package org.Hs.eg.db
#> Warning: Layer ‘data’ is empty
#> ! [2026-05-12 03:36:07] The following features were filtered because not found in the srt assay: "'de novo' NAD+ biosynthetic process from L-tryptophan", "'de novo' protein folding", "'de novo' pyrimidine nucleobase biosynthetic process", "1-phosphatidyl-1D-myo-inositol 4,5-bisphosphate metabolic process", "2'-deoxyribonucleotide biosynthetic process", "2'-deoxyribonucleotide metabolic process", "2-oxoglutarate metabolic process", "3'-UTR-mediated mRNA destabilization", "3'-UTR-mediated mRNA stabilization", "3'-phosphoadenosine 5'-phosphosulfate metabolic process", "7-methylguanosine cap hypermethylation", "ADP catabolic process", "ADP metabolic process", "ADP transport", "AMP biosynthetic process", "AMP metabolic process", "ARF protein signal transduction", "ATF6-mediated unfolded protein response", …, "zymosterol biosynthetic process", and "zymosterol metabolic process"
#> ℹ [2026-05-12 03:36:07] Number of feature lists to be scored: 0
#> Warning: Layer ‘data’ is empty
#> Error in run_aucell_scores(expr_counts = expr_sp, gene_sets = features,     strategy = cpp_strategy): No gene sets retain genes after intersecting with the expression matrix

  panc8_sub <- integration_scop(
    panc8_sub,
    assay = "GO",
    batch = "tech",
    integration_method = "Harmony"
  )
#> ◌ [2026-05-12 03:36:07] Run integration workflow...
#> Error in GetAssay.Seurat(assay_source, assay = assay_use): GO is not an assay present in the given object. Available assays are: RNA
  CellDimPlot(
    panc8_sub,
    group.by = c("tech", "celltype")
  )
#> Error in DefaultReduction(srt): Unable to find any reductions

  pancreas_sub <- CellScoring(
    pancreas_sub,
    layer = "data",
    assay = "RNA",
    db = "GO_BP",
    species = "Mus_musculus",
    termnames = panc8_sub[["GO"]]@meta.features[, "termnames"],
    method = "AUCell",
    name = "GO",
    new_assay = TRUE
  )
#> ℹ [2026-05-12 03:36:07] Start cell scoring
#> ℹ [2026-05-12 03:36:08] Data type is log-normalized
#> ℹ [2026-05-12 03:36:08] Species: "Mus_musculus"
#> 
#> ✔ [2026-05-12 03:39:44] org.Mm.eg.db installed successfully
#> ℹ [2026-05-12 03:40:15] Preparing database: GO_BP
#> ℹ [2026-05-12 03:40:27] Convert ID types for the GO_BP database
#> ℹ [2026-05-12 03:40:27] Converted ID types using local annotation package org.Mm.eg.db
#> Error in panc8_sub[["GO"]]: ‘GO’ not found in this Seurat object
#>  
  pancreas_sub <- standard_scop(
    pancreas_sub,
    assay = "GO"
  )
#> ℹ [2026-05-12 03:40:31] Start standard processing workflow...
#> Error in standard_scop_resolve_assays(srt = srt, assay = assay): `assay` must be present in <Seurat>: "GO"

  pancreas_sub[["tech"]] <- "Mouse"
  panc_merge <- integration_scop(
    srt_list = list(panc8_sub, pancreas_sub),
    assay = "GO",
    batch = "tech",
    integration_method = "Harmony"
  )
#> ◌ [2026-05-12 03:40:31] Run integration workflow...
#> Error in GetAssay.Seurat(assay_source, assay = assay_use): GO is not an assay present in the given object. Available assays are: RNA
  CellDimPlot(
    srt = panc_merge,
    group.by = c("tech", "celltype", "SubCellType", "Phase")
  )
#> Error: object 'panc_merge' not found

genenames <- make.unique(
  thisutils::capitalize(
    rownames(panc8_sub[["RNA"]]),
    force_tolower = TRUE
  )
)
names(genenames) <- rownames(panc8_sub)
panc8_sub <- RenameFeatures(
  panc8_sub,
  newnames = genenames,
  assay = "RNA"
)
#> ℹ [2026-05-12 03:40:31] Rename features for the assay: RNA
panc_merge <- integration_scop(
  srt_list = list(panc8_sub, pancreas_sub),
  assay = "RNA",
  batch = "tech",
  integration_method = "Harmony"
)
#> ◌ [2026-05-12 03:40:31] Run integration workflow...
#> ℹ [2026-05-12 03:40:31] Checking a list of <Seurat>...
#> ! [2026-05-12 03:40:31] `srt_list` have different feature names! Will subset the common features (12928) for downstream analysis
#> Warning: Different features in new layer data than already exists for counts
#> Warning: Different cells and/or features from existing assay RNA
#> Warning: Different features in new layer data than already exists for counts
#> Warning: Different features in new layer data than already exists for data
#> Warning: Different features in new layer data than already exists for scale.data
#> Warning: Different cells and/or features from existing assay RNA
#> ℹ [2026-05-12 03:40:33] Data 1/6 of the `srt_list` has been log-normalized
#> ℹ [2026-05-12 03:40:33] Perform `Seurat::FindVariableFeatures()` on 1/6 of `srt_list`...
#> ! [2026-05-12 03:40:33] Data 2/6 of the `srt_list` is "unknown"
#> ℹ [2026-05-12 03:40:33] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 2/6 of `srt_list`...
#> ℹ [2026-05-12 03:40:35] Perform `Seurat::FindVariableFeatures()` on 2/6 of `srt_list`...
#> ! [2026-05-12 03:40:35] Data 3/6 of the `srt_list` is "unknown"
#> ℹ [2026-05-12 03:40:35] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 3/6 of `srt_list`...
#> ℹ [2026-05-12 03:40:36] Perform `Seurat::FindVariableFeatures()` on 3/6 of `srt_list`...
#> ! [2026-05-12 03:40:37] Data 4/6 of the `srt_list` is "unknown"
#> ℹ [2026-05-12 03:40:37] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 4/6 of `srt_list`...
#> ℹ [2026-05-12 03:40:38] Perform `Seurat::FindVariableFeatures()` on 4/6 of `srt_list`...
#> ! [2026-05-12 03:40:38] Data 5/6 of the `srt_list` is "unknown"
#> ℹ [2026-05-12 03:40:38] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 5/6 of `srt_list`...
#> ℹ [2026-05-12 03:40:39] Perform `Seurat::FindVariableFeatures()` on 5/6 of `srt_list`...
#> ! [2026-05-12 03:40:40] Data 6/6 of the `srt_list` is "unknown"
#> ℹ [2026-05-12 03:40:40] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 6/6 of `srt_list`...
#> ℹ [2026-05-12 03:40:41] Perform `Seurat::FindVariableFeatures()` on 6/6 of `srt_list`...
#> ℹ [2026-05-12 03:40:41] Use the separate HVF from `srt_list`
#> ℹ [2026-05-12 03:40:41] Number of available HVF: 2000
#> ℹ [2026-05-12 03:40:42] Finished check
#> ℹ [2026-05-12 03:40:47] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-12 03:40:48] Perform linear dimension reduction("pca")
#> ℹ [2026-05-12 03:40:49] Perform Harmony integration
#> ℹ [2026-05-12 03:40:49] Using "Harmonypca" (1:20) as input
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 't': ‘Z_corr’ is not a valid field or method name for reference class “Rcpp_harmony”
CellDimPlot(
  srt = panc_merge,
  group.by = c("tech", "celltype", "SubCellType", "Phase")
)
#> Error: object 'panc_merge' not found
```
