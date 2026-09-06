# Run spot-level quality control

Calculate common spot-level QC metrics for spatial transcriptomics data
and label failed spots without running single-cell-specific checks such
as doublet calling or ambient RNA decontamination.

## Usage

``` r
RunSpotQC(
  srt = NULL,
  assay = NULL,
  return_filtered = FALSE,
  qc_metrics = c("outlier", "umi", "gene", "mito"),
  outlier_threshold = c("log10_nCount:lower:3", "log10_nFeature:lower:3",
    "spot_featurecount_dist:lower:3"),
  outlier_n = 1,
  UMI_threshold = 500,
  gene_threshold = 200,
  mito_threshold = 20,
  mito_pattern = c("MT-", "Mt-", "mt-"),
  mito_gene = NULL,
  verbose = TRUE,
  seed = 11,
  object = NULL
)
```

## Arguments

- srt:

  A `Seurat` object. The same object may be supplied as `object =` for
  consistency with spatial plotting APIs.

- object:

  Optional alias for `srt`. Supply exactly one of `srt` or `object`.

- assay:

  Assay to use. `NULL` uses the default assay.

- return_filtered:

  Whether to return a spot-filtered Seurat object.

- qc_metrics:

  QC metrics to apply. Available metrics are `"outlier"`, `"umi"`,
  `"gene"`, and `"mito"`.

- outlier_threshold:

  Character vector specifying outlier thresholds as
  `"metric:direction:nmads"`. Available default metrics are
  `"log10_nCount"`, `"log10_nFeature"`, and `"spot_featurecount_dist"`.

- outlier_n:

  Positive integer giving the minimum number of outlier rules required
  to fail a spot. It cannot exceed the number of `outlier_threshold`
  rules.

- UMI_threshold:

  Minimum UMI count required to pass `"umi"` QC.

- gene_threshold:

  Minimum detected gene count required to pass `"gene"` QC.

- mito_threshold:

  Maximum mitochondrial percentage allowed by `"mito"` QC.

- mito_pattern:

  Regex patterns used to identify mitochondrial genes.

- mito_gene:

  Optional explicit mitochondrial gene vector. When provided,
  `mito_pattern` is ignored.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Random seed for reproducibility.

## Value

A `Seurat` object with spot QC metadata columns. Cells that are not
present in the selected assay receive `NA` QC labels. A selected-assay
cell also receives `NA` when at least one requested QC is unresolved and
no requested QC has a known failure. Only cells labelled `Pass` are
retained when `return_filtered = TRUE`.

## Examples

``` r
data(visium_human_pancreas_sub)
spatial <- RunSpotQC(
  visium_human_pancreas_sub,
  assay = "Spatial"
)
#> ◌ [2026-09-06 22:35:26] Running spot-level quality control
#> ✔ [2026-09-06 22:35:26] Spot QC completed: 1986 evaluated, 1907 Pass, 79 Fail
#> ℹ   Scope assay "Spatial", layer "counts"
#> ℹ   Saved metadata column `SpotQC`
#> ℹ   Plot returned object `SpatialSpotPlot(<returned_object>, group.by = "SpotQC")`
SpatialSpotPlot(spatial, group.by = "SpotQC")
```
