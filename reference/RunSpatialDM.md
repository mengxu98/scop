# Run official SpatialDM spatial ligand-receptor association analysis

Run official SpatialDM spatial ligand-receptor association analysis

## Usage

``` r
RunSpatialDM(
  srt,
  species = c("human", "mouse"),
  lr.database = NULL,
  assay = NULL,
  layer = "data",
  counts.layer = "counts",
  image = NULL,
  coord.cols = c("col", "row"),
  coordinate.unit = NULL,
  method = c("z-score", "permutation"),
  run.local = TRUE,
  l = NULL,
  eff_dist = NULL,
  cutoff = 0.1,
  n_neighbors = NULL,
  n_neighbor_layers = 6,
  single_cell = FALSE,
  complex.mean = c("algebra", "geometric"),
  min_cell = 3,
  n_perm = 1000,
  nproc = 1,
  seed = 1,
  global.fdr = TRUE,
  global.threshold = 0.1,
  local.fdr = FALSE,
  local.threshold = 0.1,
  result.name = "default",
  envname = "spatialdm_env",
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A `Seurat` spatial object.

- species:

  Built-in SpatialDM species (`"human"` or `"mouse"`).

- lr.database:

  Optional custom LR table with `ligand`, `receptor`, and `annotation`
  columns.

- assay, layer, counts.layer:

  Assay and normalized/count expression layers.

- image, coord.cols:

  Spatial image or metadata coordinate columns.

- coordinate.unit:

  Unit of the raw input coordinates. Automatic unit selection uses
  technology-specific rules.

- method:

  Official global/local test mode.

- run.local:

  Whether to run locally.

- l, eff_dist:

  SpatialDM RBF scale; one must be supplied.

- cutoff:

  Maximum selected-cell fraction used by Scissor's alpha search.

- n_neighbors:

  Number of neighbors to use for constructing the KNN graph.

- n_neighbor_layers:

  Number of neighbor layers to use.

- single_cell:

  Whether the input is a single-cell (rather than spatial) dataset.

- complex.mean:

  Whether to use a modified mean test statistic in complex heatmaps.

- min_cell:

  Minimum number of cells required.

- n_perm:

  Number of permutations.

- nproc:

  Number of parallel processes.

- seed:

  Random seed.

- global.fdr:

  Global FDR threshold.

- global.threshold:

  Global threshold.

- local.fdr:

  Local FDR threshold.

- local.threshold:

  Local threshold.

- result.name:

  Stored result name.

- envname:

  Conda-compatible Python environment. If `NULL`, the environment name
  will be set to `"scop_env"`.

- overwrite:

  Whether to overwrite existing metadata.

- verbose:

  Whether to print messages.

## Value

A Seurat object with a schema-v1 `SpatialDM` result bundle.
