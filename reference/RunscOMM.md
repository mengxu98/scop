# Run scOMM label prediction

Run `scOMM` on shared features between a reference object and a query
object, write predicted labels and class scores into query metadata, and
optionally evaluate predictions against a truth label.

## Usage

``` r
RunscOMM(
  srt,
  reference,
  reference_assay = NULL,
  query_assay = NULL,
  reference_label = NULL,
  features = NULL,
  prediction_prefix = "predicted_",
  evaluate = FALSE,
  truth_col = NULL,
  tool_name = "scOMM",
  rare_threshold = 0.05,
  scomm_python = NULL,
  scomm_hidden_nodes = c(128, 64),
  scomm_epochs = 10,
  scomm_batch_size = 32,
  scomm_threshold = 0.5,
  scomm_seed = 11,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- reference:

  Reference `Seurat` object used for supervision.

- reference_assay:

  Assay used in the reference object.

- query_assay:

  Assay used in the query object.

- reference_label:

  Metadata column in the reference used as supervision labels.

- features:

  Shared features passed to `scOMM`. If `NULL`, reference variable
  features are used.

- prediction_prefix:

  Prefix added to prediction metadata columns.

- evaluate:

  Whether to compute prediction metrics against a truth label.

- truth_col:

  Metadata column in `srt` used as the truth label when
  `evaluate = TRUE`.

- tool_name:

  Name used to store detailed results in `srt@tools`.

- rare_threshold:

  Maximum class proportion used to define rare classes when calculating
  `rare_recall`.

- scomm_python:

  Optional Python binary used by the `scOMM` backend. If `NULL`,
  `SCOP_SCOMM_PYTHON` is consulted and reticulate defaults are used
  otherwise.

- scomm_hidden_nodes, scomm_epochs, scomm_batch_size, scomm_threshold,
  scomm_seed:

  Parameters passed to the `scOMM` backend.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A `Seurat` object with `scOMM` predictions stored in metadata and
`tools`.
