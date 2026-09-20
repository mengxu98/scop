# Spatial imports and standard workflow methods

This page is the platform and method reference for the [main spatial
workflow](https://mengxu98.github.io/scop/articles/spatial-main-workflow.md).
It concentrates input contracts: assay, image, observation unit,
coordinate space, and the methods that can consume them. The example
paths below are templates for an existing local vendor output directory.
SCOP does not download large data or install optional backends during
import.

## Platform contracts at a glance

| Input | Observation | Typical assay/image | Coordinate meaning | Natural SCOP entry |
|----|----|----|----|----|
| Visium | spot | `Spatial` / `slice1` | full-resolution image pixels; display scale is separate | [`SpatialSpotPlot()`](https://mengxu98.github.io/scop/reference/SpatialSpotPlot.md), spot QC, Moran/Geary, domain clustering, spot deconvolution |
| Visium HD | bin | `Spatial.008um` or `Spatial.016um` / matching image | full-resolution pixels for the selected bin resolution | [`SpatialSpotPlot()`](https://mengxu98.github.io/scop/reference/SpatialSpotPlot.md), resolution-aware BayesSpace, spatial-variable features |
| Xenium | segmented cell | vendor assay / `fov` | vendor segmentation coordinates in microns | [`SpatialCellPlot()`](https://mengxu98.github.io/scop/reference/SpatialCellPlot.md), cell-level context or communication |

The exact assay and image names are read from the loaded object; do not
copy a name from a different resolution or sample. Check the selected
context before calling a method:

``` r

SpatialCoordinates(object, image = "slice1", space = "raw")
SeuratObject::DefaultAssay(object)
SeuratObject::Images(object)
```

SCOP keeps spot IDs/cell IDs aligned with expression columns. If several
images cover one sample, analysis and plotting require an explicit
selection or an explicit named sample-to-image map. A display
`lowres`/`hires` scale never changes raw distances and does not turn
pixels into microns.

## Import existing vendor output

These paths are deliberately marked as templates. They are not evaluated
by the site build.

``` r

library(scop)

visium <- ReadSpatialData(
  data.dir = "path/to/visium/outs",  # existing local output directory
  technology = "visium",
  sample_id = "S1"
)

visium_hd <- ReadSpatialData(
  data.dir = "path/to/visium-hd/outs",  # existing local output directory
  technology = "visium_hd",
  bin.size = c(8, 16),
  sample_id = "HD1"
)

xenium <- ReadSpatialData(
  data.dir = "path/to/xenium/outs",  # existing local output directory
  technology = "xenium",
  sample_id = "X1"
)
```

[`ReadSpatialData()`](https://mengxu98.github.io/scop/reference/ReadSpatialData.md)
delegates file parsing to Seurat and records import provenance under
`object@misc$scop_spatial_input`. It validates that each returned image,
assay, observation set, resolution, and coordinate-unit label agree. It
does not normalize counts, segment cells, or infer cell types.

Public data sources for testing the input formats are available from the
[10x Genomics Visium datasets](https://www.10xgenomics.com/datasets),
the [10x Visium HD human-pancreas
example](https://www.10xgenomics.com/datasets/visium-hd-cytassist-11mm-human-pancreas),
the [10x Visium HD developer example
data](https://www.10xgenomics.com/support/software/space-ranger/latest/resources/visium-hd-example-data),
and [10x Xenium Explorer demo
datasets](https://www.10xgenomics.com/support/software/xenium-explorer/latest/resources/xenium-explorer-demos).
Download and unpack them outside the package; do not hard-code a machine
path in a reproducible article.

## Visium and Visium HD: keep counts and normalized data separate

Raw counts are the input for count-dependent QC and methods that model
sampling depth. Normalized/log-transformed data are appropriate for
expression-level plots and methods whose contract asks for a normalized
layer. `SpaNorm` is a separate normalization stage; it does not replace
the original count source.

``` r

normalized <- RunStandardWorkflow(
  visium,
  workflow = "spatial",
  assay = "Spatial",
  image = "slice1",
  normalization_method = "SpaNorm",
  spanorm_params = list(new_assay = "SpaNormAnalysis"),
  do_spatial_cluster = TRUE,
  spatial_cluster_method = "BANKSY",
  do_deconvolution = FALSE
)

# Count-dependent method: keep the original assay/layer explicit.
RunBayesSpace(
  normalized,
  assay = "Spatial",
  image = "slice1",
  layer = "counts",
  store_sce = FALSE
)
```

The workflow attaches matching original counts when needed for QC and
uses the separate normalized assay for downstream expression
preprocessing. It rejects an explicit `do_normalization = FALSE`
conflict with `SpaNorm`. Do not relabel normalized values as counts.

## Standard methods by platform

For a Visium spot or HD bin, a compact method selection is:

``` r

baseline <- RunStandardWorkflow(
  visium,
  workflow = "spatial",
  assay = "Spatial",
  image = "slice1",
  do_spot_qc = TRUE,
  do_spatial_variable_features = TRUE,
  spatial_variable_features_params = list(method = "moran", nfeatures = 50),
  do_spatial_cluster = FALSE,
  do_deconvolution = FALSE
)

banksy <- RunStandardWorkflow(
  visium,
  workflow = "spatial",
  assay = "Spatial",
  image = "slice1",
  do_spot_qc = FALSE,
  do_spatial_variable_features = FALSE,
  do_spatial_cluster = TRUE,
  spatial_cluster_method = "BANKSY",
  spatial_cluster_params = list(lambda = 0.2, resolution = 0.6),
  do_deconvolution = FALSE
)

SpatialSpotPlot(banksy, group.by = "BANKSY_cluster", image = "slice1")
```

`BayesSpace` is the spot/bin array branch; `spatial_q` or BayesSpace’s
`q` controls the number of domains. `BANKSY` uses its own resolution
parameters. `SmoothClust` accepts a requested `n_clusters` and its
smoothing parameters. Segmented cell observations are not silently sent
through spot deconvolution or BayesSpace; use cell-level methods and
plots instead.

## Multi-sample PRECAST integration

Integration needs a sample label and a covering image per sample. With a
merged object, the image map names are the merged object’s image keys.
With list input, the map names refer to the original list element names
before Seurat renames duplicate image keys.

``` r

integrated <- RunSpatialIntegration(
  object = list(S1 = visium_S1, S2 = visium_S2),
  method = "PRECAST",
  assay = "Spatial",
  layer = "counts",
  image = c(S1 = "slice1", S2 = "slice1"),
  store_object = FALSE
)

SpatialIntegrationPlot(
  integrated,
  plot_type = "spatial"
)
```

The plotting call resolves the images in the merged object. Inspect
`SeuratObject::Images(integrated)` before supplying an explicit plotting
map: duplicate input names such as `slice1` may have been renamed during
merging. Do not reuse the input-list map as a merged-object map.

The PRECAST package is an execution-time optional dependency. The
template is not a successful run unless PRECAST is installed and the
input objects exist. `store_object = FALSE` omits only the complete
native PRECAST object; domains, embedding, aligned coordinates,
parameters, source information, and standard plots remain stored. Keep
the default `TRUE` when a later native PRECAST analysis needs the full
object.

## Deconvolution and cell2location storage choices

Spot deconvolution is reference-dependent and estimates composition per
spot; it does not turn a spot into a segmented cell.
[`RunRCTD()`](https://mengxu98.github.io/scop/reference/RunRCTD.md),
[`RunSPOTlight()`](https://mengxu98.github.io/scop/reference/RunSPOTlight.md),
and
[`RunCell2location()`](https://mengxu98.github.io/scop/reference/RunCell2location.md)
preserve their own result semantics. Use
[`SpatialDeconvolutionPlot()`](https://mengxu98.github.io/scop/reference/SpatialDeconvolutionPlot.md)
or
[`Cell2locationPlot()`](https://mengxu98.github.io/scop/reference/Cell2locationPlot.md)
after checking the stored result. Cell2location q05 abundance and
normalized proportions are distinct; the abundance view uses its native
range rather than a forced 0–1 scale. Cell2location plotting matrices
retain the original spot IDs. Spots excluded from modeling, such as
zero-count observations, have `NA` values rather than zero abundance.
The saved `cells` field identifies modeled spots, and
`input_summary$dropped_spots` records filtered observations. Historical
compact results are expanded only when those records explain all missing
rows; missing predictions or unknown IDs are not silently accepted. An
explicit `lower_cutoff`, `upper_cutoff`, or quantile control overrides
the default matrix-wide display range. With a one-sided override, each
panel keeps the remaining bound computed by
[`SpatialSpotPlot()`](https://mengxu98.github.io/scop/reference/SpatialSpotPlot.md);
specify both cutoffs when you want identical explicit limits across
panels. `combine` changes only the return form, including for the
dominant point map.

``` r

deconvolved <- RunRCTD(
  visium,
  reference = reference,
  reference_label = "celltype",
  assay = "Spatial",
  reference_assay = "RNA",
  layer = "counts",
  reference_layer = "counts",
  image = "slice1",
  rctd_mode = "full",
  max_cores = 1
)
SpatialDeconvolutionPlot(deconvolved, tool_name = "RCTD", plot_type = "dominant")

Cell2locationPlot(
  cell2location_result,
  plot_type = "abundance",
  combine = FALSE
)
```

## Validation and analysis receipts

Acceptance runs use public slide data outside the package: developer
fixtures check format and computation rather than biological performance
on production slides, and public slide data are not redistributed with
scop.

For a new dataset, preserve the following in the analysis receipt:

- platform and observation unit (`spot`, `bin`, or `cell`);
- assay and layer, selected image, coordinate columns, coordinate space,
  and physical units when known;
- raw-count versus normalized-expression input decisions;
- sample-to-image map for every sample;
- stage status and backend versions for methods that actually ran.
