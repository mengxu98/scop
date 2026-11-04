# WNN integration function

Combine paired RNA and ATAC representations with weighted nearest
neighbors, retaining both modality reductions and the WNN graph.

## Usage

``` r
WNN_integrate(
  srt_merge = NULL,
  batch = NULL,
  append = TRUE,
  srt_list = NULL,
  assay = NULL,
  do_normalization = NULL,
  normalization_method = "LogNormalize",
  do_HVF_finding = TRUE,
  HVF_source = "separate",
  HVF_method = "vst",
  nHVF = 2000,
  HVF_min_intersection = 1,
  HVF = NULL,
  do_scaling = TRUE,
  vars_to_regress = NULL,
  regression_model = "linear",
  scale_within_batch = FALSE,
  linear_reduction = "pca",
  linear_reduction_dims = 50,
  linear_reduction_dims_use = NULL,
  linear_reduction_params = list(),
  force_linear_reduction = FALSE,
  nonlinear_reduction = "umap",
  nonlinear_reduction_dims = c(2, 3),
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- srt_merge:

  A merged \`Seurat\` object that includes the batch information.

- batch:

  Accepted for integration-interface compatibility. WNN does not perform
  batch correction; supplying this argument emits a warning.

- append:

  Append integrated results to `srt_merge`.

- srt_list:

  A list of `Seurat` objects to be checked and preprocessed.

- assay:

  One RNA assay and one chromatin assay, or one assay whose partner can
  be inferred. Automatic selection prefers an active RNA/SCT assay, then
  RNA, then SCT; custom RNA assay names must be explicit. Known SCOP
  gene-activity and chromVAR assays cannot substitute for RNA.

- do_normalization:

  Whether data normalization should be performed.

- normalization_method:

  The normalization method to be used. Possible values are
  `"LogNormalize"`, `"SCT"`, `"TFIDF"`, and `"scran"`.

- do_HVF_finding:

  Whether to select variable features during preprocessing.

- HVF_source:

  The source of highly variable features. Possible values are `"global"`
  and `"separate"`.

- HVF_method:

  RNA feature-selection method passed to
  [`RunStandardWorkflow()`](https://mengxu98.github.io/scop/reference/RunStandardWorkflow.md).

- nHVF:

  Maximum number of selected features, subject to modality defaults and
  the number of available features.

- HVF_min_intersection:

  The feature needs to be present in batches for a minimum number of
  times in order to be considered as highly variable.

- HVF:

  Optional RNA feature names. ATAC features are selected independently.

- do_scaling:

  Force scaling via ScaleData.

- vars_to_regress:

  A vector of variable names to include as additional regression
  variables.

- regression_model:

  `"linear"`, `"poisson"`, or `"negativebinomial"`.

- scale_within_batch:

  Scale within each batch. Only used by `"Uncorrected"`, `"Seurat"`,
  `"MNN"`, `"Harmony"`, `"BBKNN"`, `"CSS"`, `"ComBat"`.

- linear_reduction:

  RNA reduction, which must be `"pca"`. ATAC uses SVD.

- linear_reduction_dims:

  Number of components to compute for each modality.

- linear_reduction_dims_use:

  Explicit component indices, or NULL for modality-specific selection.
  Inspect `WNN_dims_list` for the selected indices.

- linear_reduction_params:

  Named arguments for the linear reduction backend.

- force_linear_reduction:

  Whether to recompute existing linear reductions.

- nonlinear_reduction:

  One or more of `"umap"`, `"umap-naive"` and `"fr"`, or
  [`character()`](https://rdrr.io/r/base/character.html) for graphs and
  clusters without an embedding.

- nonlinear_reduction_dims:

  Number(s) of nonlinear dimensions to compute.

- nonlinear_reduction_params:

  Shared parameters, optionally overridden by lists named after
  individual nonlinear reduction methods.

- force_nonlinear_reduction:

  Recompute nonlinear reductions even when stored results exist.

- neighbor_metric, neighbor_k:

  Neighbor graph (`"euclidean"`, `"cosine"`, `"manhattan"`,
  `"hamming"`).

- cluster_algorithm, cluster_resolution:

  Clustering (`"louvain"`, `"slm"`, `"leiden"`). Larger
  `cluster_resolution` yields fewer clusters.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Random seed.

## Value

A Seurat object with RNA PCA and ATAC LSI reductions, `WNN` neighbors,
`WNNKNN` and `WNNSNN` graphs, `WNNclusters` metadata and requested WNN
embeddings. `Misc(object, "WNN_reduction_list")` names the retained
source reductions; `WNN_dims_list` records their selected dimensions.
Both return modes retain complete processed modality assays (including
scaled RNA), ATAC reductions, metadata and representation receipts, with
the processed RNA assay (SCT when selected) active. The caller's
random-number state is restored.

## Details

Cells must represent paired measurements. RNA and ATAC preprocessing is
followed by native multimodal neighbor construction. This function does
not use batch labels to correct either modality. Specify assays
explicitly after adding derived assays, for example
`assay = c("RNA", "peaks")`. Both assays must contain the same cells as
the object. Subset to verified paired cells explicitly before running
WNN. Nonlinear reductions supported here are `"umap"`, `"umap-naive"`
and `"fr"`; use [`character()`](https://rdrr.io/r/base/character.html)
to retain graphs and clusters without computing an embedding. RNA uses
the requested normalization and PCA. ATAC always uses TF-IDF and SVD,
including when RNA uses SCT. `HVF` supplies RNA features; ATAC features
are selected independently from accessibility counts.

## See also

[`ReadMultiomeData()`](https://mengxu98.github.io/scop/reference/ReadMultiomeData.md),
[`RunStandardWorkflow()`](https://mengxu98.github.io/scop/reference/RunStandardWorkflow.md),
[`RunLinkPeaks()`](https://mengxu98.github.io/scop/reference/RunLinkPeaks.md)

## Examples

``` r
# Small paired subset; no external fragments or reference downloads.
data("pbmcmultiome_sub", package = "scop")
paired <- subset(pbmcmultiome_sub, cells = colnames(pbmcmultiome_sub)[1:80])
paired <- WNN_integrate(
  srt_merge = paired, assay = c("RNA", "peaks"),
  linear_reduction_dims = 10, linear_reduction_dims_use = 2:8,
  nHVF = 100, neighbor_k = 10, nonlinear_reduction_dims = 2,
  verbose = FALSE
)
#> ℹ [2026-11-04 19:15:16] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-11-04 19:15:25] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-11-04 19:15:35] Reorder clusters...
#> ℹ [2026-11-04 19:15:35] Skip `log1p()` because `layer = data` is not "counts"
SeuratObject::Reductions(paired)
#> [1] "RNApca"        "RNApcaUMAP2D"  "RNAUMAP2D"     "ATACsvd"      
#> [5] "ATACsvdUMAP2D" "ATACUMAP2D"    "ATAClsi"       "WNNUMAP2D"    
stopifnot(all(paired@misc$WNN_reduction_list %in% SeuratObject::Reductions(paired)))
```
