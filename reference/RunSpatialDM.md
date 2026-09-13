# Run official SpatialDM spatial ligand-receptor association analysis

Run official SpatialDM spatial ligand-receptor association analysis

## Usage

``` r
RunSpatialDM(
  object,
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
  verbose = TRUE,
  srt = NULL
)
```

## Arguments

- object:

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

  RBF weight cutoff used by SpatialDM to retain local spatial weights.
  This is a fraction between zero and one; it is not a Scissor
  alpha-search parameter.

- n_neighbors:

  Number of neighbors to use for constructing the KNN graph.

- n_neighbor_layers:

  Number of neighbor layers to use.

- single_cell:

  Whether to use the single-cell spatial mode supported by SpatialDM.
  This is not a generic "single-cell rather than spatial" switch.

- complex.mean:

  Mean used for complex ligand-receptor expression: `"algebra"` or
  `"geometric"`.

- min_cell:

  Minimum number of cells required.

- n_perm:

  Number of permutations.

- nproc:

  Number of parallel processes.

- seed:

  Random seed.

- global.fdr, local.fdr:

  Whether the corresponding p-value table uses FDR values for selection.
  These are logical switches, not numeric thresholds.

- global.threshold, local.threshold:

  Selection thresholds for global and local results.

- result.name:

  Stored result name.

- envname:

  Conda-compatible Python environment. If `NULL`, the environment name
  will be set to `"scop_env"`.

- overwrite:

  Whether to overwrite existing metadata.

- verbose:

  Whether to print messages.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

A Seurat object with a schema-v1 `SpatialDM` result bundle.
