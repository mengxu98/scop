# Run official COMMOT spatial communication analysis

Run official Python COMMOT in an isolated subprocess. Raw expression,
metadata, and raw coordinates are exchanged through files; AnnData is
not placed inside the Seurat object.

## Usage

``` r
RunCOMMOT(
  srt,
  group.by,
  species = c("human", "mouse"),
  database = "CellChat",
  assay = NULL,
  layer = "counts",
  image = NULL,
  coord.cols = c("col", "row"),
  distance.threshold = NULL,
  cluster = FALSE,
  direction = FALSE,
  communication.args = list(),
  cluster.args = list(),
  direction.args = list(),
  result.name = "default",
  envname = "commot_env",
  output.dir = NULL,
  store.h5ad = FALSE,
  overwrite = FALSE,
  backend = c("cpp", "r"),
  verbose = TRUE
)
```

## Arguments

- srt:

  A \`Seurat\` spatial object.

- group.by:

  Metadata column used to aggregate communication by group.

- species:

  COMMOT ligand-receptor database species.

- database:

  Official COMMOT ligand-receptor database name.

- assay, layer:

  Assay and raw-count layer.

- image:

  Explicit spatial image. Required for multi-image objects.

- coord.cols:

  Metadata coordinate columns used when no image is present.

- distance.threshold:

  Maximum signaling distance in raw-coordinate units.

- cluster:

  Whether to run official cluster communication permutations.

- direction:

  Whether to run official communication direction analysis.

- communication.args:

  Named arguments for official spatial communication. Runner-only fields
  are \`normalize\`, \`target_sum\`, \`signaling_type\`, and nested
  \`filter_args\`.

- cluster.args:

  Named arguments for official cluster communication.

- direction.args:

  Named arguments for official communication direction.

- result.name:

  Stored result name.

- envname:

  Isolated Python environment name.

- output.dir:

  External directory for a requested H5AD artifact.

- store.h5ad:

  Whether to retain the complete AnnData as an external H5AD.

- overwrite:

  Whether to replace an existing named result or H5AD artifact.

- backend:

  Backend used only for SCOP CCC result aggregation.

- verbose:

  Whether to print progress messages.

## Value

The input \`Seurat\` object with \`srt@tools\$COMMOT\` and unified CCC
results updated.

## References

\<https://github.com/zcang/COMMOT\>
