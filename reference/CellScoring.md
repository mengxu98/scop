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
  will be used.

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
  `"GO", "GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway", "Reactome", "CORUM", "MP", "DO", "HPO", "PFAM", "CSPA", "Surfaceome", "SPRomeDB", "VerSeDa", "TFLink", "hTFtarget", "TRRUST", "JASPAR", "ENCODE", "MSigDB", "CellTalk", "CellChat", "Chromosome", "GeneType", "Enzyme", "TF"`.

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
#> ℹ [2026-04-02 15:27:31] Start standard processing workflow...
#> ℹ [2026-04-02 15:27:32] Checking a list of <Seurat>...
#> ! [2026-04-02 15:27:32] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-02 15:27:32] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 15:27:33] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 15:27:34] Use the separate HVF from `srt_list`
#> ℹ [2026-04-02 15:27:34] Number of available HVF: 2000
#> ℹ [2026-04-02 15:27:34] Finished check
#> ℹ [2026-04-02 15:27:34] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-02 15:27:35] Perform pca linear dimension reduction
#> ℹ [2026-04-02 15:27:39] Use stored estimated dimensions 1:50 for Standardpca
#> Warning: Caught FutureLaunchError. Canceling all iterations ...
#> ! [2026-04-02 15:27:39] <FutureLaunchError: Caught an unexpected error of class FutureLaunchError when trying to launch future (‘future_lapply-1’) on backend of class SequentialFutureBackend. The reason was: future::evalFuture() failed on runnervmrg6be (pid 85355) at 2026-04-02T15:27:39. Using package 'future' v1.70.0. Possible other reasons: Failed to attach one or more future-backend packages: there is no package called ‘future’ [future <unnamed>; on 4a75d434f7a9a2903adedbeee3372830@runnervmrg6be<85355>] [future ‘future_lapply-1’ (4a75d434f7a9a2903adedbeee3372830-8); on 4a75d434f7a9a2903adedbeee3372830@runnervmrg6be<85355>]>
#> !                       
#> !                       Occurred on: 4a75d434f7a9a2903adedbeee3372830 [runnervmrg6be; pid 85355]
#> !                       Future: 4a75d434f7a9a2903adedbeee3372830-8 (‘future_lapply-1’)
#> !                       
#> !                       DEBUG: BEGIN TROUBLESHOOTING HELP
#> !                       SequentialFuture:
#> !                       Label: ‘future_lapply-1’
#> !                       Expression:
#> Error in glue(str, .envir = .envir, .transformer = transformer, .cli = TRUE,     .trim = .trim): Expecting '}'
features_all <- rownames(pancreas_sub)
pancreas_sub <- CellScoring(
  pancreas_sub,
  features = list(
    A = features_all[1:100],
    B = features_all[101:200]
  ),
  method = "Seurat",
  name = "test"
)
#> ℹ [2026-04-02 15:27:39] Start cell scoring
#> Warning: Layer ‘data’ is empty
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> ! [2026-04-02 15:27:40] Infinite values detected
#> Warning: Layer ‘data’ is empty
#> ! [2026-04-02 15:27:40] The following features were filtered because not found in the srt assay: "A" and "B"
#> ℹ [2026-04-02 15:27:40] Number of feature lists to be scored: 0
#> Warning: Layer ‘data’ is empty
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'rowMeans': subscript out of bounds
CellDimPlot(pancreas_sub, "test_classification")
#> Error in CellDimPlot(pancreas_sub, "test_classification"): "test_classification" is not in the meta.data of srt object

FeatureDimPlot(pancreas_sub, "test_A")
#> Error in DefaultReduction(srt): Unable to find any reductions

if (FALSE) { # \dontrun{
data(panc8_sub)
panc8_sub <- integration_scop(
  panc8_sub,
  batch = "tech",
  integration_method = "Seurat"
)
CellDimPlot(
  panc8_sub,
  group.by = c("tech", "celltype")
)

panc8_sub <- CellScoring(
  panc8_sub,
  layer = "data",
  assay = "RNA",
  db = "GO_BP",
  species = "Homo_sapiens",
  minGSSize = 10,
  maxGSSize = 100,
  method = "Seurat",
  name = "GO",
  new_assay = TRUE
)

panc8_sub <- integration_scop(
  panc8_sub,
  assay = "GO",
  batch = "tech",
  integration_method = "Seurat"
)
CellDimPlot(
  panc8_sub,
  group.by = c("tech", "celltype")
)

pancreas_sub <- CellScoring(
  pancreas_sub,
  layer = "data",
  assay = "RNA",
  db = "GO_BP",
  species = "Mus_musculus",
  termnames = panc8_sub[["GO"]]@meta.features[, "termnames"],
  method = "Seurat",
  name = "GO",
  new_assay = TRUE
)
pancreas_sub <- standard_scop(
  pancreas_sub,
  assay = "GO"
)
CellDimPlot(pancreas_sub, "SubCellType")

pancreas_sub[["tech"]] <- "Mouse"
panc_merge <- integration_scop(
  srt_list = list(panc8_sub, pancreas_sub),
  assay = "GO",
  batch = "tech", integration_method = "Seurat"
)
CellDimPlot(
  srt = panc_merge,
  group.by = c("tech", "celltype", "SubCellType", "Phase")
)

genenames <- make.unique(
  thisutils::capitalize(
    rownames(panc8_sub[["RNA"]])
  ),
  force_tolower = TRUE
)
names(genenames) <- rownames(panc8_sub)
panc8_sub <- RenameFeatures(
  panc8_sub,
  newnames = genenames,
  assay = "RNA"
)
head(rownames(panc8_sub))
panc_merge <- integration_scop(
  srt_list = list(panc8_sub, pancreas_sub),
  assay = "RNA",
  batch = "tech", integration_method = "Seurat"
)
CellDimPlot(
  srt = panc_merge,
  group.by = c("tech", "celltype", "SubCellType", "Phase")
)
} # }
```
