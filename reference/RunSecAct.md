# Run SecAct secreted protein activity inference

Run the optional `SecAct` R package from `scop` to infer secreted
protein signaling activity from bulk, single-cell, or spatial
transcriptomics profiles. `SecAct` is not bundled with `scop`; it is
checked and installed from `data2intelligence/SecAct` when this function
is called.

## Usage

``` r
RunSecAct(
  srt = NULL,
  inputProfile = NULL,
  inputProfile_control = NULL,
  mode = c("auto", "matrix", "scRNAseq", "ST"),
  assay = NULL,
  layer = "counts",
  cellType_meta = NULL,
  is.singleCellLevel = FALSE,
  is.differential = FALSE,
  is.paired = FALSE,
  is.singleSampleLevel = FALSE,
  scale.factor = NULL,
  sigMatrix = "SecAct",
  is.filter.sig = FALSE,
  is.group.sig = TRUE,
  is.group.cor = 0.9,
  lambda = 5e+05,
  nrand = 1000,
  ncores = 1L,
  backend = "auto",
  rng_method = "mt19937",
  batch_size = NULL,
  output_h5 = NULL,
  activity = c("zscore", "beta", "se", "pvalue"),
  assay_out = "SecAct",
  store_assay = !is.null(srt),
  store_results = TRUE,
  tool_name = "SecAct",
  return_seurat = !is.null(srt),
  verbose = TRUE
)
```

## Arguments

- srt:

  Optional Seurat object. When `mode = "scRNAseq"`, this is passed to
  `SecAct.activity.inference.scRNAseq`. When `mode = "matrix"`,
  expression is extracted from `srt` if `inputProfile` is not supplied.

- inputProfile:

  Expression matrix, Seurat object, or SpaCET object. Matrix input must
  be genes x samples, cells, or spatial spots.

- inputProfile_control:

  Optional control expression matrix or SpaCET object passed to SecAct.

- mode:

  SecAct workflow. `auto` dispatches from the input class.

- assay:

  Assay used when extracting a matrix from `srt`; also sets the default
  assay for Seurat input passed to SecAct.

- layer:

  Assay layer used when `mode = "matrix"` and `inputProfile` is not
  supplied.

- cellType_meta:

  Metadata column containing cell type or state labels for
  `mode = "scRNAseq"` when `is.singleCellLevel = FALSE`.

- is.singleCellLevel:

  Whether SecAct should return single-cell activity rather than
  cell-state activity for Seurat input.

- is.differential, is.paired, is.singleSampleLevel:

  Parameters passed to `SecAct.activity.inference` for matrix or bulk
  input.

- scale.factor:

  Spot-level scale factor passed to `SecAct.activity.inference.ST`.

- sigMatrix:

  SecAct signature matrix name.

- is.filter.sig, is.group.sig, is.group.cor, lambda, nrand, ncores,
  backend, rng_method:

  Parameters passed to SecAct activity inference.

- batch_size, output_h5:

  Optional large-matrix controls passed only to
  `SecAct.activity.inference`.

- activity:

  Activity matrix to store as a Seurat assay when possible.

- assay_out:

  Name of the Seurat assay used to store activity values.

- store_assay:

  Whether to store a SecAct activity matrix as a Seurat assay when its
  columns match Seurat cells.

- store_results:

  Whether to store raw SecAct results in `srt@tools`.

- tool_name:

  Name used in `srt@tools`.

- return_seurat:

  Whether to return a Seurat object when `srt` is supplied.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A Seurat object, SpaCET object, or SecAct result list depending on input
and `return_seurat`.

## Examples

``` r
if (FALSE) { # \dontrun{
srt <- RunSecAct(
  srt,
  mode = "scRNAseq",
  cellType_meta = "celltype",
  is.singleCellLevel = TRUE
)

res <- RunSecAct(
  inputProfile = expr_mat,
  mode = "matrix",
  is.differential = TRUE
)
} # }
```
