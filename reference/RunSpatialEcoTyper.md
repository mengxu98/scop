# Run SpatialEcoTyper spatial ecotype analysis

Run SpatialEcoTyper workflows through the optional `SpatialEcoTyper`
package and write spatial ecotype labels or abundances back to a
`Seurat` object when possible.

## Usage

``` r
RunSpatialEcoTyper(
  srt,
  mode = c("single", "multi", "recover", "deconvolute"),
  assay = NULL,
  layer = "data",
  celltype.by = NULL,
  sample.by = NULL,
  x.by = "X",
  y.by = "Y",
  dat = NULL,
  celltypes = NULL,
  features = NULL,
  outprefix = NULL,
  outdir = NULL,
  radius = 50,
  resolution = 0.5,
  nfeatures = 300,
  min.cts.per.region = 2,
  npcs = 20,
  min.cells = 5,
  min.features = 10,
  iterations = 10,
  minibatch = 5000,
  ncores = 4,
  grid.size = round(radius * 1.4),
  filter.region.by.celltypes = NULL,
  k = 20,
  k.sn = 50,
  dropcell = FALSE,
  normalization.method = "None",
  nmf_ranks = 10,
  nrun.per.rank = 30,
  min.coph = 0.95,
  Region = NULL,
  downsample.by.region = TRUE,
  subresolution = 30,
  seed = 1,
  scale = TRUE,
  Ws = NULL,
  ncell.per.run = 500,
  min.score = 0.6,
  W = NULL,
  nsample.per.run = 500,
  sum2one = TRUE,
  prefix = "SpatialEcoTyper",
  tool_name = "SpatialEcoTyper",
  store_results = TRUE,
  allow_partial = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object. For `mode = "deconvolute"`, a numeric expression
  matrix can also be supplied.

- mode:

  SpatialEcoTyper workflow. `"single"` runs single-sample de novo
  discovery, `"multi"` runs conserved ecotype discovery across samples,
  `"recover"` recovers pretrained SE labels, and `"deconvolute"` infers
  SE abundances from bulk or spot-level expression.

- assay:

  Assay used for expression extraction. If `NULL`, the default assay is
  used.

- layer:

  Assay layer used for expression extraction.

- celltype.by:

  Metadata column containing cell type annotations. Required for
  `"single"`, `"multi"`, and `"recover"` unless `celltypes` is supplied
  for `"recover"`.

- sample.by:

  Metadata column identifying samples for `mode = "multi"`.

- x.by, y.by:

  Metadata columns containing single-cell spatial coordinates.

- dat:

  Optional expression matrix used by `"recover"` or `"deconvolute"`. If
  `NULL`, expression is extracted from `srt`.

- celltypes:

  Optional named vector of cell types passed to
  `SpatialEcoTyper::RecoverSE()`.

- features:

  Optional feature vector used to subset the expression matrix.

- outprefix:

  Output prefix passed to `SpatialEcoTyper`. Use `NULL` to avoid writing
  single-sample result files to the working directory.

- outdir:

  Output directory passed to multi-sample SpatialEcoTyper. `NULL`
  creates a temporary directory.

- radius:

  Spatial neighborhood radius, in the same units as `x.by` and `y.by`.

- resolution:

  Louvain clustering resolution used by `SpatialEcoTyper`.

- nfeatures:

  Number of variable features used by `SpatialEcoTyper`.

- min.cts.per.region:

  Minimum number of cell types required in a spatial neighborhood.

- npcs:

  Number of principal components used for similarity networks.

- min.cells:

  Minimum number of cells or spatial meta-cells expressing a feature.

- min.features:

  Minimum number of features detected in a cell or spatial meta-cell.

- iterations:

  Number of similarity network fusion iterations.

- minibatch:

  Number of columns processed per mini-batch in SNF.

- ncores:

  Number of CPU cores used by `SpatialEcoTyper`.

- grid.size:

  Spatial grid size used to discretize coordinates.

- filter.region.by.celltypes:

  Optional cell types used to restrict spatial neighborhoods.

- k:

  Number of spatial nearest neighbors used to construct spatial
  meta-cells.

- k.sn:

  Number of nearest neighbors used to construct similarity networks.

- dropcell:

  Whether cells without spatial ecotype assignments are removed from the
  returned `SpatialEcoTyper` metadata.

- normalization.method, nmf_ranks, nrun.per.rank, min.coph, Region,
  downsample.by.region, subresolution, seed:

  Parameters passed to `SpatialEcoTyper::MultiSpatialEcoTyper()`.

- scale:

  Whether to scale expression for `"recover"` and `"deconvolute"`.

- Ws:

  Pretrained basis matrices passed to `SpatialEcoTyper::RecoverSE()`.

- ncell.per.run:

  Number of cells processed per run by `SpatialEcoTyper::RecoverSE()`.

- min.score:

  Minimum prediction score passed to `SpatialEcoTyper::RecoverSE()`.

- W:

  Pretrained basis matrix passed to `SpatialEcoTyper::DeconvoluteSE()`.

- nsample.per.run:

  Number of samples processed per run by
  `SpatialEcoTyper::DeconvoluteSE()`.

- sum2one:

  Whether inferred SE abundances are normalized to sum to one.

- prefix:

  Prefix used for output metadata columns.

- tool_name:

  Name used to store detailed results in `srt@tools`.

- store_results:

  Whether to store raw results in `srt@tools`.

- allow_partial:

  Whether to allow missing SE labels for cells absent from returned
  `SpatialEcoTyper` metadata. Default is `FALSE` to avoid silent partial
  annotations.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments passed to the selected SpatialEcoTyper function.

## Value

A `Seurat` object with SpatialEcoTyper results in metadata and raw
results stored in `srt@tools[[tool_name]]` when `store_results = TRUE`.
For matrix input with `mode = "deconvolute"`, the abundance matrix is
returned.

## Examples

``` r
if (FALSE) { # \dontrun{
srt <- RunSpatialEcoTyper(
  srt,
  celltype.by = "CellType",
  x.by = "X",
  y.by = "Y"
)

srt <- RunSpatialEcoTyper(
  srt,
  mode = "multi",
  celltype.by = "CellType",
  sample.by = "sample",
  x.by = "X",
  y.by = "Y"
)
} # }
```
