# Infer gene regulatory networks with GENIE3

Run GENIE3 regulatory network inference and return a standardized
adjacency table with columns \`TF\`, \`target\`, and \`importance\`.

## Usage

``` r
RunGENIE3(object, ...)

# S3 method for class 'Seurat'
RunGENIE3(
  object,
  assay = NULL,
  layer = "counts",
  regulators = NULL,
  targets = NULL,
  max_edges_per_target = Inf,
  output_file = NULL,
  cores = 1,
  force = FALSE,
  verbose = TRUE,
  ...
)

# S3 method for class 'matrix'
RunGENIE3(object, ...)

# Default S3 method
RunGENIE3(
  object,
  regulators = NULL,
  targets = NULL,
  genes_in = c("rows", "columns"),
  max_edges_per_target = Inf,
  output_file = NULL,
  cores = 1,
  force = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A Seurat object or expression matrix.

- ...:

  Additional backend-specific arguments.

- assay:

  Assay used when \`object\` is a Seurat object.

- layer:

  Assay layer used when \`object\` is a Seurat object.

- regulators:

  Candidate transcription factor genes.

- targets:

  Optional target genes. If \`NULL\`, all genes are considered.

- max_edges_per_target:

  Maximum incoming regulator edges retained per target. The default
  \`Inf\` keeps all positive-importance links.

- output_file:

  Optional path where the adjacency table is written.

- cores:

  Number of workers used by GENIE3.

- force:

  Whether to rebuild existing \`output_file\`.

- verbose:

  Whether to print progress messages.

- genes_in:

  Matrix orientation for matrix inputs. \`"rows"\` means genes x cells;
  \`"columns"\` means cells x genes.

## Value

A data frame with columns \`TF\`, \`target\`, and \`importance\`.
